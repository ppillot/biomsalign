/**
 * @file noalign
 * @author Paul Pillot <paul.pillot@libmol.org>
 * @license MIT
 * @copyright 2021
 * @description This module contains code to make a fast alignment computation
 * without doing the full Dynamic Programming routine as it is overexpensive
 * for long sequences, running in O(n2) time and space. For example aligning
 * a pair of mitochondrial DNA sequences (around 16kpb) requires 128MB for the
 * storage of the trace-back (TB) matrix.
 * It is noticeable that most of the TB is unused, and most of the DP
 * computations are useless when done in regions of low homology between
 * sequences.
 * To speed-up the process for long sequences, we'll try to restrict the DP
 * to the smallest regions possible. This entails to detect the regions where
 * the sequences are identical between sequences to anchor alignments from there
 * without doing any DP recursion.
 * Here is the proposed algorithm for aligning a pair of long genomic sequences
 * SeqA and SeqB:
 * -1. Extract minimizers from SeqA and SeqB and store them in hashes with their
 *     position and the actual sequence they minimize.
 *     Minimizers (Roberts 2004) are efficient constructs based on K-mers that
 *     store in a compressed way a longer sequence for allowing fast comparisons
 *     A K-mer is a substring of fixed length (for example 8 nucleotides).
 *     Storing and comparing all K-mers in a sequence is doable but takes time
 *     and space. A Minimizer is the K-mer with the minimal value in a larger
 *     window (for example a window of 30nt). The minimal value can be the
 *     lowest lexicographic order for example.
 *     If two sequences share an identical window of 30nt, they also share
 *     the same minimizer. Minimizers can be an order of magnitude less
 *     numerous than the exhaustive K-mers enumeration, and they allow to
 *     compare larger segments in one iteration.
 * -2. Find common minimizers and collect the ones minimizing common substrings
 *     At the end of this process a collection of ranges each associated with a
 *     diagonal has been established.
 * -3. Merge/filter ranges. Contiguous ranges on the same diagonal can be merged
 *     A range that is overlapped on its both sides by other ranges can be
 *     removed from the collection. REF to add here
 * -4. Find an optimal order of the remaining diagonals (similarly to a
 *     Longest Increasing Sequence problem).
 * -5. Extend diagonals with exact and approximate matches
 * -6. Solve alignments between diagonals using the DP method (should be very
 *     fast then as the remaining DP matrices should be a fraction of the
 *     initial one)
 *
 *
 * (a) Roberts M. et al. BIOINFORMATICS Vol. 20 no. 18 2004, pages 3363–3369
 * doi:10.1093/bioinformatics/bth408 - Reducing storage requirements for
 * biological sequence comparison
 *
 */

import { TSequence } from "./sequence";
import { DEQueue } from './queue';
import { DEBUG } from "./params";
import Log from "./logger";
import { hammingWeight } from "./bitarray";
import { pairwiseAlignment } from "./align";

export type TMinimizer = {
    kmer: number,
    kmerPos: number,
    winPos: number,
    winPosEnd: number
};

type TKmer = Pick<TMinimizer, 'kmer'|'kmerPos'>;

type TRange = {
    diagId: number,
    begin: number,
    end: number
};

/**
 *
 *
 * @export
 * @param {TSequence} seq Sequence object
 *
 * @see https://doi.org/10.1093/bioinformatics/btaa435 Chirag Jain, Arang Rhie,
 * Haowen Zhang, Claudia Chu, Brian P Walenz, Sergey Koren, Adam M Phillippy,
 * Weighted minimizer sampling improves long read mapping, Bioinformatics,
 * Volume 36, Issue Supplement_1, July 2020, Pages i111–i118
 */
export function extractMinimizers (seq: TSequence, ksize: number, wsize: number): [Map<number, TMinimizer[]>, TMinimizer[], Uint16Array] {

    const lMinzMap: Map<number, TMinimizer[]> = new Map();
    const lMinzArr: TMinimizer[] = [];
    const lQueue = new DEQueue<number>(wsize);
    const lKSIZE = ksize | 0;

        // Create an array of 8-kmers from sequence
        // 8-kmers, with 2 bits per letter in a 4 letter alphabet require 16
        // bits. Let's compute it using a sliding window (in O(n) )

    const lKarr = new Uint16Array(seq.encodedSeq.length - lKSIZE + 1);
    let kval = 0|0;
    let lIdx = 0;
    for (let i = 0; i < lKSIZE; i++) {
        kval |= seq.encodedSeq[lKSIZE - i - 1] << (i * 2);
    }
    lKarr[lIdx++] = kval;

    for (let i = lKSIZE, imax = seq.encodedSeq.length; i < imax; i++) {
        kval = (kval << 2) + seq.encodedSeq[i];
        lKarr[lIdx++] = kval;
    }

        // Compute minimizers

    const lStartStore = wsize - lKSIZE;

    let lPrevMinz: TMinimizer = {kmer: -1, kmerPos: -1, winPosEnd: -1, winPos: -1};

    for (let i = 0; i < lKarr.length; i++) {
        const lKmer = lKarr[i];

            // Remove kmers with lower order than current kmer
            // Note: for sake of symplicity, order is just kmer value.

        while (!lQueue.isEmpty && lKarr[lQueue.getTail()] >= lKmer) {
            lQueue.popTail();
        }
        lQueue.pushTail(i);

        if (i<lStartStore) continue;

            // Remove kmers that are not in this window anymore

        while (lQueue.getHead() <= i - 1 - lStartStore) {
            lQueue.popHead();
        }

            // Store head kmer in minimizer hash

        let lHead = lQueue.getHead();

            // Is it the same minimizer as last time?
            // If yes, only extend the minimized subarray

        if (lHead === lPrevMinz.kmerPos) {
            lPrevMinz.winPosEnd = i + lKSIZE;
        } else {
            let lHeadKmer = lKarr[lHead];
            lPrevMinz = {
                kmer: lHeadKmer,
                kmerPos: lHead,
                winPos: i - lStartStore,
                winPosEnd: i + lKSIZE
            }

                // Add the kmer to the hash, either as a new list or as a new
                // item in a previous list (when the same kmer also minimizes
                // another window)

            if (!lMinzMap.has(lHeadKmer)) {
                lMinzMap.set(lHeadKmer, [lPrevMinz]);
            } else {
                let lList = lMinzMap.get(lHeadKmer) as TMinimizer[];
                lList.push(lPrevMinz);
            }

            lMinzArr.push(lPrevMinz);
        }

    }

    return [lMinzMap, lMinzArr, lKarr];
}


export function noalignPair(seqA: TSequence, seqB: TSequence) {

    const KSIZE = 8;    // fits in 16bits
    const WSIZE = 16;
    const EXTENSION_THRESHOLD = KSIZE/4; // don't extend when more than 25% difference
    let lDebugStats: {[k: string]: Partial<{a: any, b: any, all: any}>} = {};

    const [lMinzA, lMinzAArr, lKmerA] = extractMinimizers(seqA, KSIZE, WSIZE);
    const [lMinzB, lMinzBArr, lKmerB] = extractMinimizers(seqB, KSIZE, WSIZE);
    if (DEBUG) {
        Log.add('Extract Minimizers');
        lDebugStats['Nb Minimizers'] = {a: lMinzAArr.length, b: lMinzBArr.length};
        lDebugStats['Dupl. Minimizers'] = {a: lMinzAArr.length - lMinzA.size, b: lMinzBArr.length - lMinzB.size};
    }

    const lRangesColl = [];
    const lDiagMap = new Map<number, TRange[]>();

    // TODO: break the tie
    for (let i = 0; i < lMinzAArr.length; i++) {
        let kmer = lMinzAArr[i].kmer;

        if (!lMinzB.has(kmer)) continue;
        let listB = lMinzB.get(kmer) as TMinimizer[];

            // Compare minimized string in A with those in listB to retain only
            // the ones that share common minimized strings

        let lMinzSubA = seqA.rawSeq.substring(lMinzAArr[i].winPos, lMinzAArr[i].winPosEnd);

        for (let j = 0; j < listB.length; j++) {
            let lMinzSubB =  seqB.rawSeq.substring(listB[j].winPos, listB[j].winPosEnd);
            let lLen = Math.min(lMinzSubA.length, lMinzSubB.length);
            if (lLen == lMinzSubA.length) {
                if (lMinzSubB.indexOf(lMinzSubA) !== 0) continue;
            } else if (lMinzSubA.indexOf(lMinzSubB) !== 0) continue;

            // common range to store
            let lDiagId = lMinzAArr[i].winPos - listB[j].winPos;
            let lRange = {
                diagId: lDiagId,
                begin : lMinzAArr[i].winPos,
                end: lMinzAArr[i].winPos + lLen
            };
            lRangesColl.push(lRange);
            let lDiagList = lDiagMap.get(lDiagId);
            if (lDiagList === undefined) {
                lDiagMap.set(lDiagId, [lRange])
            } else {
                lDiagList.push(lRange);
            }
        }
    }
    if (DEBUG) {
        Log.add('Filter Minimizers');
        lDebugStats['Nb common'] = {all: lRangesColl.length};
    }

    let lDiagList: TRange[] = [];
    // Merge consecutive or overlapping subarrays on the same diagonals
    // Note: the extensions are made on range references
    lDiagMap.forEach((diags) => {
        let lCurrentSegment = diags[0];
        lDiagList.push(lCurrentSegment);
        diags.forEach(range => {
            if (range.begin <= lCurrentSegment.end) {
                lCurrentSegment.end = range.end;
            } else {
                lCurrentSegment = range;
                lDiagList.push(lCurrentSegment);
            }
        });
    });
    if (DEBUG) {
        Log.add('Merge Minimizers');
        lDebugStats['Nb diagonals'] = {all: lDiagList.length};
    }

    // Sort ranges in diag list from left to right and from top to bottom
    lDiagList.sort((a, b) => {
        let lDelta = a.begin - b.begin;
        return (lDelta === 0) ? a.diagId - b.diagId: lDelta;
    });
    if (DEBUG) Log.add('Sort Minimizers');

    // TODO.Find optimal order of diagonals using LIS like algorithm

    // Extend consecutive segments on same diagonal. We apply a tolerance for
    // extension at 25% difference between Kmers.
    // This is measured by bitwise operations between 16bits numbers.

    let lPrevSegment = lDiagList[0];
    let lExtDiagList = [lPrevSegment];
    let lMissingSegments: any[] = [];

    for (let i = 1; i < lDiagList.length; i++) {
        let lCurrentSegment = lDiagList[i];

        // TODO: when difference between end and begin is >> KSIZE (e.g. 400),
        // there are likely unmapped gaps. It might be useful to speed up
        // matches using shorter Kmers/windows. Test the overhead + noise
        // generated.

        // first case: this segment is on the same diagonal as the previous one
        // Try to extend begin position towards previous segment.

        if (lCurrentSegment.diagId === lPrevSegment.diagId) {
            // Compare kmers in forward direction until score drops
            let m = lPrevSegment.end;   // end position on seq A
            let n = m - lCurrentSegment.diagId;

            while (m < lCurrentSegment.begin
                && dnaHammingDistance(lKmerA[m], lKmerB[n]) <= EXTENSION_THRESHOLD
            ) {
                m += KSIZE;
                n += KSIZE;
            }

            // To avoid boundaries conditions, remove the tip of the extension
            if (m > lPrevSegment.end) {
                lPrevSegment.end = m - EXTENSION_THRESHOLD;
            }

            // Compare in backward direction until score drops or merge has been
            // achieved

            m = lCurrentSegment.begin - KSIZE;
            n = m - lCurrentSegment.diagId;
// WIP Make stats on full window, not just 1 Kmer
            while (m > lPrevSegment.end
                && dnaHammingDistance(lKmerA[m], lKmerB[n]) <= EXTENSION_THRESHOLD
            ) {
                m -= KSIZE;
                n -= KSIZE;
            }

            if (m < lCurrentSegment.begin - KSIZE) {
                lCurrentSegment.begin = m + KSIZE + EXTENSION_THRESHOLD;
            }

            // merge?
            if (lCurrentSegment.begin - lPrevSegment.end <= 3 * KSIZE) {
                lPrevSegment.end = lCurrentSegment.end;
                continue;
            } else {
                lMissingSegments.push({
                    beginDiagId: lPrevSegment.diagId,
                    endDiagId: lCurrentSegment.diagId,
                    begin: lPrevSegment.end,
                    end: lCurrentSegment.begin
                });
                lExtDiagList.push(lCurrentSegment);
                lPrevSegment = lCurrentSegment;
                continue;
            }
        } else {
            // Compare nucleotides in forward direction until score drops
            let m = lPrevSegment.end;   // end position on seq A
            let n = m - lPrevSegment.diagId;

            while (seqA.encodedSeq[m] === seqB.encodedSeq[n]
                && m < lCurrentSegment.begin) {
                m ++;
                n ++;
            }

            lPrevSegment.end = m;

            // Compare in backward direction until sequences differ

            m = lCurrentSegment.begin - 1;
            n = m - lCurrentSegment.diagId;

            while (seqA.encodedSeq[m] === seqB.encodedSeq[n]
                && m > lPrevSegment.end) {
                m --;
                n --;
            }

            lCurrentSegment.begin = m;

            lMissingSegments.push({
                endDiagId: lCurrentSegment.diagId,
                beginDiagId: lPrevSegment.diagId,
                begin: lPrevSegment.end,
                end: lCurrentSegment.begin
            });
            lExtDiagList.push(lCurrentSegment);
            lPrevSegment = lCurrentSegment;
        }
    }

    if (DEBUG) {
        Log.add('Diagonals extension');
        let lCoverage = 0;
        lExtDiagList.forEach(d => {
            lCoverage += d.end - d.begin;
        })
        lCoverage /= seqA.encodedSeq.length;
        lDebugStats['Coverage'] = {all: lCoverage};
        lDebugStats['Extended Diagonals'] = {all: lExtDiagList.length};
        console.table(lDebugStats);

        // console.log(lExtDiagList);
        // console.log(lMissingSegments);
    }

    // Make the alignment from the diagonals and fill the gaps
    let lSeqA = '';
    let lSeqB = '';
    let lDiag: TRange;
    for (let i = 0; i < lExtDiagList.length; i++) {
        lDiag = lExtDiagList[i];
        lSeqA += seqA.rawSeq.substring(lDiag.begin, lDiag.end);
        lSeqB += seqB.rawSeq.substring(lDiag.begin - lDiag.diagId, lDiag.end - lDiag.diagId);

        let lMis = lMissingSegments[i];
        if (lMis) {
            // special case: begin === end
            if (lMis.begin === lMis.end) {
                lSeqA += '-'.repeat(Math.abs(lMis.endDiagId - lMis.beginDiagId));
                lSeqB += seqB.rawSeq.substring(lMis.begin - lMis.beginDiagId, lMis.end - lMis.endDiagId)
                continue;
            }
            if (lMis.begin - lMis.beginDiagId === lMis.end - lMis.endDiagId) {
                lSeqB += '-'.repeat(Math.abs(lMis.endDiagId - lMis.beginDiagId));
                lSeqA += seqA.rawSeq.substring(lMis.begin, lMis.end);
                continue;
            }


            let lResult = pairwiseAlignment({
                rawSeq: seqA.rawSeq.substring(lMis.begin, lMis.end),
                type: seqA.type,
                compressedSeq: new Uint8Array(0),
                encodedSeq: seqA.encodedSeq.slice(lMis.begin, lMis.end)
            }, {
                rawSeq: seqB.rawSeq.substring(lMis.begin - lMis.beginDiagId, lMis.end - lMis.endDiagId),
                type: seqB.type,
                compressedSeq: new Uint8Array(0),
                encodedSeq: seqB.encodedSeq.slice(lMis.begin - lMis.beginDiagId, lMis.end - lMis.endDiagId)
            },[0],[1]);

            lSeqA += lResult.alignment[0];
            lSeqB += lResult.alignment[1];
        }

    }
    if (DEBUG) {
        Log.add('Fill between diagonals');
    }

    // TODO: Finalize alignment here
    return [lSeqA, lSeqB];
}

/**
 * Returns the number of different letters in DNA alphabet between 2 kmers
 * encoded as numbers, using bitwise operations
 *
 * @param {number} a
 * @param {number} b
 */
export function dnaHammingDistance (a: number, b: number) {

    // DNA alphabet encodes 4 symbols using 2 bits. In the following, we'll
    // count the pairs of bits that are different between a and b.
    let c = a ^ b;  // set bit for pair of bits that are different
    c = (c & 0x55555555) | ((c >> 1) & 0x55555555); // intersect with 010101... to retain only the lowest bit in each pair
    return hammingWeight(c);

}
