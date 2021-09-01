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

import { TSequence } from "../sequence/sequence";
import { DEQueue } from '../utils/queue';
import { DEBUG, TAlignmentParam } from "./params";
import Log from "../utils/logger";
import { hammingWeight } from "../utils/bitarray";
import { ALIGNOPT, pairwiseAlignment } from "./align";
import { epath2estring, EPATH_2_STRING, estringCat, estringCountPositive,
    estringDifference, estringLength, estringMerge, estringProduct } from "../utils/estring";

type TRange = {
    diagId: number,
    begin: number,
    end: number
};

type TMinzStore = {
    kmer: Uint16Array,
    kmerPos: Uint16Array,
    winPos: Uint16Array,
    winPosEnd: Uint16Array,
    count: number
}

const KSIZE = 8;    // fits in 16bits
const WSIZE = 16;
const EXTENSION_THRESHOLD = KSIZE/4; // don't extend when more than 25% difference


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
export function extractMinimizers (seq: TSequence, ksize: number, wsize: number): TMinzComp {

    const lMinzMap: Map<number, number[]> = new Map();
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

        // random kmer density is 2/(wsize + 1). This doubles the expected size
        // for the store.

    const lLength = (seq.rawSeq.length * (2 / (wsize + 1)) * 2)|0 ;
    const lStoreMinz: TMinzStore = {
        kmer     : new Uint16Array(lLength),
        kmerPos  : new Uint16Array(lLength),
        winPos   : new Uint16Array(lLength),
        winPosEnd: new Uint16Array(lLength),
        count: 0
    };
    let lCurrMinzPos = NaN;

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

        if (lHead === lCurrMinzPos) {
            lStoreMinz.winPosEnd[lStoreMinz.count - 1] = i + lKSIZE;
            //lPrevMinz.winPosEnd = i + lKSIZE;
        } else {
            const lHeadKmer = lKarr[lHead];
            const lIdx = lStoreMinz.count;
            lStoreMinz.count ++;

            lStoreMinz.kmer     [lIdx] = lHeadKmer;
            lStoreMinz.kmerPos  [lIdx] = lHead;
            lStoreMinz.winPos   [lIdx] = i - lStartStore;
            lStoreMinz.winPosEnd[lIdx] = i + lKSIZE;
            lCurrMinzPos = lHead;

                // Add the kmer to the hash, either as a new list or as a new
                // item in a previous list (when the same kmer also minimizes
                // another window)

            if (!lMinzMap.has(lHeadKmer)) {
                lMinzMap.set(lHeadKmer, [lIdx]);
            } else {
                let lList = lMinzMap.get(lHeadKmer) as number[];
                lList.push(lIdx);
            }

        }

    }

    return [lMinzMap, lStoreMinz, lKarr];
}

/** Map Kmer value to a table of all its indices in the sequence */
type TMinzMap = Map<number, number[]>;
type TMinzComp = [TMinzMap, TMinzStore, Uint16Array];

export function noalignPair(
    seqA: TSequence,
    seqB: TSequence,
    pAlignParam: TAlignmentParam,
    pMinzA?: TMinzComp,
    pMinzB?: TMinzComp
) {

    let lDebugStats: {[k: string]: Partial<{a: any, b: any, all: any}>} = {};

    const [lMinzA, lMinzAStore, lKmerAArr] = pMinzA ?? extractMinimizers(seqA, KSIZE, WSIZE);
    const [lMinzB, lMinzBStore, lKmerBArr] = pMinzB ?? extractMinimizers(seqB, KSIZE, WSIZE);
    if (DEBUG) {
        Log.add('Extract Minimizers');
        lDebugStats['Nb Minimizers'] = {a: lMinzAStore.count, b: lMinzBStore.count};
        lDebugStats['Dupl. Minimizers'] = {a: lMinzAStore.count - lMinzA.size, b: lMinzBStore.count - lMinzB.size};
    }

    let lNbCommonWindows = 0;
    const lDiagMap = new Map<number, TRange[]>();

    // Loop through all kmers fron sequence A to find which ones are matching
    // in sequence B. Doing so, favour the kmers that minimize a common window.
    // Those are better seed candidates for diagonals. They also cover a longer
    // range and thus avoid extra computations.

    for (let i = 0; i < lMinzAStore.count; i++) {
        let kmer = lMinzAStore.kmer[i];

        if (!lMinzB.has(kmer)) continue;
        let listB = lMinzB.get(kmer) as number[];

            // Compare minimized string in A with those in listB to take advantage of
            // the ones that share common minimized strings

        let lMinzSubA = seqA.rawSeq.substring(lMinzAStore.winPos[i], lMinzAStore.winPosEnd[i]);

        for (let j = 0; j < listB.length; j++) {
            let lBidx = listB[j];
            let lMinzSubB =  seqB.rawSeq.substring(lMinzBStore.winPos[lBidx], lMinzBStore.winPosEnd[lBidx]);
            let lLen = Math.min(lMinzSubA.length, lMinzSubB.length);

            let lHasCommonWindow = false;
            if (lLen == lMinzSubA.length) {
                   lHasCommonWindow = (lMinzSubB.indexOf(lMinzSubA) === 0);
            } else lHasCommonWindow = (lMinzSubA.indexOf(lMinzSubB) === 0);

            // common range to store
            let lDiagId = 0;
            let lRange: TRange;
            if (lHasCommonWindow) {
                lDiagId = lMinzAStore.winPos[i] - lMinzBStore.winPos[lBidx];
                lRange = {
                    diagId: lDiagId,
                    begin : lMinzAStore.winPos[i],
                    end   : lMinzAStore.winPos[i] + lLen
                };
                lNbCommonWindows ++;
            } else {
                lDiagId = lMinzAStore.kmerPos[i] - lMinzBStore.kmerPos[lBidx];
                lRange = {
                    diagId: lDiagId,
                    begin : lMinzAStore.kmerPos[i],
                    end   : lMinzAStore.kmerPos[i] + KSIZE - 1
                };
            }


            if (!lDiagMap.has(lDiagId)) {
                lDiagMap.set(lDiagId, [lRange]);
            } else {
                let lDiagList = lDiagMap.get(lDiagId);
                lDiagList!.push(lRange);
            }
        }
    }
    if (DEBUG) {
        Log.add('Filter Minimizers');
        lDebugStats['Nb common windows'] = { all: lNbCommonWindows };
    }

    let lDiagList: TRange[] = [];
    // Merge consecutive or overlapping subarrays on the same diagonals
    // Note: the extensions are made on range references
    lDiagMap.forEach((diags) => {
        if (!diags.length) return;  // safeguard

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

    // Special case: no diagonal found: fallback to NW alignment
    if (lDiagList.length === 0) {
        if (DEBUG) console.table(lDebugStats);

        return pairwiseAlignment(seqA, seqB, pAlignParam).estrings;
    }

    // Sort ranges in diag list from left to right and from top to bottom
    lDiagList.sort((a, b) => {
        let lDelta = a.begin - b.begin;
        return (lDelta === 0) ? a.diagId - b.diagId: lDelta;
    });
    if (DEBUG) Log.add('Sort Minimizers');

    // TODO.Find optimal order of diagonals using LIS like algorithm
    let lDiagIncreasingSequences: TRange[][] = [[lDiagList[0]]];
    let lIncrSeqSizes: number[] = [lDiagList[0].end - lDiagList[0].begin];
    for (let i = 0; i < lDiagList.length; i++) {
        let lDiag = lDiagList[i];
        let lDiagLen = lDiag.end - lDiag.begin;

        // Loop over the list of increasing sequences and see to which one this
        // diagonal could be appended to.
        for (let j = 0; j < lDiagIncreasingSequences.length; j++) {

                // Can lDiag be the start of a new list of ordered diagonals?


        }
    }

    // Extend consecutive segments on same diagonal. We apply a tolerance for
    // extension at 25% difference between Kmers.
    // This is measured by bitwise operations between 16bits numbers.

    let lPrevSegment = lDiagList[0];
    let lExtDiagList = lPrevSegment ? [lPrevSegment] : [];
    let lMissingSegments: any[] = [];

    for (let i = 1; i < lDiagList.length; i++) {
        let lCurrentSegment = lDiagList[i];

        // Temporary fix here: due to crude sorting of diagonals, overlaps are
        // not taken into account. The following prevents intractable situations
        // but it does not guarantee the best outcome.

        if (lCurrentSegment.begin < lPrevSegment.end) continue;

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
            let lWinScore = 0;

            while (m < lCurrentSegment.begin
                && lWinScore < 2
            ) {
                m += KSIZE;
                n += KSIZE;
                lWinScore = dnaHammingDistance(lKmerAArr[m], lKmerBArr[n]) <= EXTENSION_THRESHOLD ? 0 : lWinScore + 1;
            }

            // To avoid boundaries conditions, remove the tip of the extension
            if (m > lPrevSegment.end) {
                lPrevSegment.end = m - EXTENSION_THRESHOLD - KSIZE;
            }

            // Compare in backward direction until score drops or merge has been
            // achieved

            m = lCurrentSegment.begin - KSIZE;
            n = m - lCurrentSegment.diagId;
            lWinScore = 0;
// WIP Make stats on full window, not just 1 Kmer
            while (m > lPrevSegment.end
                && lWinScore < 2
            ) {
                m -= KSIZE;
                n -= KSIZE;
                lWinScore = dnaHammingDistance(lKmerAArr[m], lKmerBArr[n]) <= EXTENSION_THRESHOLD ? 0 : lWinScore + 1;
            }

            if (m < lCurrentSegment.begin - 2 * KSIZE) {
                lCurrentSegment.begin = m + 3 * KSIZE + EXTENSION_THRESHOLD;
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
                    end: lCurrentSegment.begin,
                    size: lCurrentSegment.begin-lPrevSegment.end
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

            lPrevSegment.end = m;   // last before discrepancy. This position is
                                    // exclusive in the prevSegment, and inclusive
                                    // in the missing segment if any.

            // Compare in backward direction until sequences differ

            m = lCurrentSegment.begin;
            n = m - lCurrentSegment.diagId;

            while (seqA.encodedSeq[m] === seqB.encodedSeq[n]
                && m > lPrevSegment.end) {
                m --;
                n --;
            }

            lCurrentSegment.begin = m + 1;  // last before discrepancy

            lMissingSegments.push({
                endDiagId: lCurrentSegment.diagId,
                beginDiagId: lPrevSegment.diagId,
                begin: lPrevSegment.end,
                end: lCurrentSegment.begin,
                size: lCurrentSegment.begin-lPrevSegment.end
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
        });
        lCoverage /= seqA.encodedSeq.length;
        lDebugStats['Coverage'] = {all: lCoverage};
        lDebugStats['Extended Diagonals'] = {all: lExtDiagList.length};
        console.table(lDebugStats);

        //console.table(lMissingSegments);
    }

    // Make the alignment from the diagonals and fill the gaps
    let lDiag: TRange = lExtDiagList[0];
    const lEpathA: number[] = [];
    const lEpathB: number[] = [];
    for (let i = 0; i < lExtDiagList.length; i++) {
        lDiag = lExtDiagList[i];
        lEpathA.push(lDiag.end - lDiag.begin);
        lEpathB.push(lDiag.end - lDiag.begin);

        let lMis = lMissingSegments[i];
        if (lMis) {
            // special case: begin === end
            if (lMis.begin === lMis.end) {
                lEpathA.push(-Math.abs(lMis.endDiagId - lMis.beginDiagId));
                lEpathB.push(0); // lMis.end - lMis.end
                continue;
            }
//TODO: this looks incorrect!? lEpathB twice?
            if (lMis.begin - lMis.beginDiagId === lMis.end - lMis.endDiagId) {
                lEpathB.push(- Math.abs(lMis.endDiagId - lMis.beginDiagId));
                lEpathB.push(lMis.end - lMis.begin);
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
            }, pAlignParam,
            ALIGNOPT.DISABLE_FAVOR_END_GAP | ALIGNOPT.DISABLE_FAVOR_START_GAP);

            lEpathA.push(...lResult.estrings[0]);
            lEpathB.push(...lResult.estrings[1]);
        }

    }

        // Append tail if any
        // At this stage both aligned sequences have the same length.
        // For now, just append the missing part. Ultimately, assess, depending
        // on the tail length if a SW alignmnent is deemed necessary.

    sanitizeEpathEnd(lEpathA, seqA.encodedSeq.length);
    sanitizeEpathEnd(lEpathB, seqB.encodedSeq.length);

    let lEstringA = epath2estring(lEpathA, EPATH_2_STRING.NO_REVERSE);
    let lEstringB = epath2estring(lEpathB, EPATH_2_STRING.NO_REVERSE);

    const lEszA = estringLength(lEstringA);
    const lEszB = estringLength(lEstringB);
    if (lEszA > lEszB) lEstringB = estringCat(lEstringB, [lEszB - lEszA]);
    else if (lEszA < lEszB) lEstringA = estringCat(lEstringA, [lEszA - lEszB]);

    if (DEBUG) {
        Log.add('Fill between diagonals');
    }

    // TODO: Finalize alignment here
    return [lEstringA, lEstringB];
}

/**
 * Given an epath, it ensures that the sum of its positive values equals the
 * target length.
 *
 * @param {number[]} epath
 * @param {number} targetLength
 */
function sanitizeEpathEnd(epath: number[], targetLength: number) {
    let lCount = estringCountPositive(epath);
    let lDelta = targetLength - lCount;
    if (lDelta > 0) epath.push(lDelta);
    else if (lDelta < 0) epath[epath.length - 1] += lDelta;
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

/**
 * Align multiple (long) sequences using the center-star method.
 *
 * @export
 * @param {TSequence[]} pSeq
 * @param {TAlignmentParam} pAlignParam
 * @returns
 */
export function centerStarNoAlign(pSeq: TSequence[], pAlignParam: TAlignmentParam) {

    let lMinz = pSeq.map(seq => { return extractMinimizers(seq, KSIZE, WSIZE) });

    // TODO: find center sequence
    let lCenter = 0;
    let lEStrings: number[][] = [];
    let lCenterEString: number[] = [];
    let lCenterEStringRaw: number[][] = [];

    for (let i = 0; i < lMinz.length; i++) {
        if (i === lCenter) continue;

        let lEsAln = noalignPair(pSeq[lCenter], pSeq[i], pAlignParam, lMinz[lCenter], lMinz[i]);

        // Center Edit strings combines the Edit strings from the alignment to the
        // current sequence with its previous edit strings.
        lCenterEString = lCenterEString.length === 0
            ? lEsAln[0]
            : estringMerge(lCenterEString, lEsAln[0]);

        // Keep track of what was this alignment edit string for center sequence
        // This is used to make a diff, at next iteration.
        lCenterEStringRaw[i] = lEsAln[0];

        // Edit string to align this sequence with the center. It will be
        // modified to account for other sequences introduced editions
        lEStrings[i] = lEsAln[1];
    }

    for (let i = 0; i < lMinz.length; i++) {
        if (i === lCenter) {
            lEStrings[i] = lCenterEString;
            continue;
        }
        let lDiff = estringDifference(lCenterEString, lCenterEStringRaw[i]);
        if (!lDiff) {
            console.error('An error occured while computing the center diff. for sequence #' + i, lCenterEString, lCenterEStringRaw[i]);
            console.dir(lEStrings)
            console.dir(lCenterEStringRaw)
            throw new RangeError();
        }
        lEStrings[i] = estringProduct(lDiff, lEStrings[i]);
    }

    return lEStrings;

}
