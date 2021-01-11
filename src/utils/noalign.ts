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
 *     At the end of this process a collection of ranges each associated to a
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

export type TMinimizer = {
    kmer: number,
    kmerPos: number,
    winPos: number,
    minimizedSubarray: number[]
};

type TKmer = Pick<TMinimizer, 'kmer'|'kmerPos'>;

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
export function extractMinimizers (seq: TSequence, ksize: number, wsize: number) {

    const lMinzMap: Map<number, TMinimizer[]> = new Map();
    const lQueue = new DEQueue<TKmer>(wsize);

        // Create an array of 8-kmers from sequence
        // 8-kmers, with 2 bits per letter in a 4 letter alphabet require 16
        // bits. Let's compute it using a sliding window (in O(n) )

    const lKarr = new Uint16Array(seq.encodedSeq.length - ksize + 1);
    let kval = 0;
    let lIdx = 0;
    for (let i = 0; i < ksize; i++) {
        kval |= seq.encodedSeq[ksize - i - 1] << (i * 2);
    }
    lKarr[lIdx++] = kval;

    for (let i = ksize, imax = seq.encodedSeq.length; i < imax; i++) {
        kval = (kval << 2) + seq.encodedSeq[i];
        lKarr[lIdx++] = kval;
    }

        // Compute minimizers

    const lStartStore = wsize - ksize;

    let lPrevMinz: TMinimizer = {kmer: -1, kmerPos: -1, minimizedSubarray: [], winPos: -1};

    for (let i = 0; i < lKarr.length; i++) {
        const lKmer = lKarr[i];
        let lMinz: TKmer = {
            kmer: lKmer,
            kmerPos: i
        }

            // Remove kmers with lower order than current kmer
            // Note: for sake of symplicity, order is just kmer value.

        while (!lQueue.isEmpty && lQueue.getTail().kmer > lKmer) {
            lQueue.popTail();
        }
        lQueue.pushTail(lMinz);

        if (i<lStartStore) continue;

            // Remove kmers that are not in this window anymore

        while (lQueue.getHead().kmerPos <= i - 1 - lStartStore) {
            lQueue.popHead();
        }

            // Store head kmer in minimizer hash

        let lHead = lQueue.getHead();

            // Is it the same minimizer as last time?
            // If yes, only extend the minimized subarray

        if (lHead.kmerPos === lPrevMinz.kmerPos) {
            lPrevMinz.minimizedSubarray.push(seq.encodedSeq[i + ksize - 1]);
        } else {
            lPrevMinz = {
                kmer: lHead.kmer,
                kmerPos: lHead.kmerPos,
                winPos: i - lStartStore,
                minimizedSubarray: seq.encodedSeq.slice(i - lStartStore, i + ksize)
            }

                // Add the kmer to the hash, either as a new list or as a new
                // item in a previous list (when the same kmer also minimizes
                // another window)

            if (!lMinzMap.has(lHead.kmer)) {
                lMinzMap.set(lHead.kmer, [lPrevMinz]);
            } else {
                let lList = lMinzMap.get(lHead.kmer) as TMinimizer[];
                lList.push(lPrevMinz);
            }

        }

    }

    return lMinzMap;
}
