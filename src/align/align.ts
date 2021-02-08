/**
 * @file align.ts
 * @author Paul Pillot <paul.pillot@libmol.org>
 * @license MIT
 * @copyright 2021
 * @description Routines for NW global alignment of biological sequences.
 * Most of this code is inspired by R.C. Edgar MUSCLE library (a).
 *
 * (a) Edgar, R.C. MUSCLE: a multiple sequence alignment method with reduced
 * time and space complexity. BMC Bioinformatics 5, 113 (2004).
 * https://doi.org/10.1186/1471-2105-5-113
 */

import { getAlignmentParameters } from "./params";
import { profileFromMSA } from "../sequence/profile";
import { TSequence } from "../sequence/sequence";
import { InternalNode, LeafNode } from "../sequence/tree";

/**
 * Trace back matrix transitions values
 * Bits are used to create a bit mask that encodes every transition between the
 * 3 classical matrices
 * @enum {number}
 */
const enum TRACE_BACK {
    MATCH     = 0,          // Match extension
    DEL       = 1,          // Delete extension
    INS       = 1 << 1,     // Insert extension
    MATCH2DEL = 1 << 2,     // Delete opening
    MATCH2INS = 1 << 3      // Insert opening
};

export const enum ALIGNOPT {
    DISABLE_FAVOR_START_GAP = 1,
    DISABLE_FAVOR_END_GAP   = 1 << 1
}


/**
 * Pairwise sequence alignment
 *
 * @export
 * @param {TSequence} seqA
 * @param {TSequence} seqB
 * @param {number} [opt=0]
 * @returns
 */
export function pairwiseAlignment (
    seqA: TSequence,
    seqB: TSequence,
    opt = 0
) {
    const params = getAlignmentParameters();

    const lSeqALen = seqA.rawSeq.length;
    const lSeqBLen = seqB.rawSeq.length;
    let lMatch = 0,      // match score

            // Textbook approach considers 3 n x m matrices for holding
            // the match / delete / insert values computed. These could be
            // reduced to only 2 vectors and 2 scalars

        lMatchArr = [],     // Match column
        lDelArr = [],       // Delete column
        lPrevMatch = 0, // Scalar value containing match value computed at
                        // the previous iteration in the same column (j-1).
                        // Used as a rotation variable.
        lLastInsert = 0,

        tbIdx = 0,
        isOdd = 0;

    // Traceback matrix size: for convenience, an extra column has been added
    // but will remain half empty (TODO: fix iteration boundaries).
    // The traceback values are stored in a typed array of (n+1) * m / 2 bytes.
    // Using bit values for M, I and D transitions, shifted depending on the matrix
    // they encode, it's possible to encode all states using 4 bits.
    // As a cell in the array contains 8 bits, each ultimately contains the
    // results for 2 successive values in the DP matrix.

    const tbM = new Uint8Array(Math.ceil((lSeqALen + 1) * lSeqBLen / 2)); // Trace back matrix
    const sA = seqA.encodedSeq;
    const sB = seqB.encodedSeq;

    let gapOpenA = 0,
        gapOpenB = 0,
        gapExtentA = 0,
        gapExtentB = 0,
        tb = 0,         // trace back value. It's a bitmask of match + del + ins

        matrix = params.scoringMatrix[0];    // memoization of the scoring matrix row for the
                        // i-th amino acid.

    // Fill in first columns in the matrix
    // INITIALISATION
    // half gap open penalty is a tweak to favor gaps at the opening
    // (also at the end) of the sequence, the rational being that
    // it is common to have ragged alignments between distant sequences.
    // (This optimization comes from MAFFT and is applied in MUSCLE).
    // Note that some applications may require on the contrary that the
    // opening gap penalty is similar between internal and external gaps.
    // e.g.: partial alignment between diagonals of a larger alignment.


    const GAP_OPEN = params.gapOP;

    const GAP_START_CORRECTION = opt & ALIGNOPT.DISABLE_FAVOR_START_GAP
        ? 0
        : - GAP_OPEN / 2;

    const GAP_END_CORRECTION = opt & ALIGNOPT.DISABLE_FAVOR_END_GAP
        ? 0
        : - GAP_OPEN / 2;


    lMatchArr[0] = 0;
    lDelArr[0] = -Infinity;
    for (let j = 1; j <= lSeqBLen; j++) {
        lMatchArr[j] = GAP_START_CORRECTION; // No gap extension penalty
        lDelArr[j] = -Infinity;
    }

    // Main DP routine

    for (let i = 1; i <= lSeqALen; i++) {

        lPrevMatch = GAP_START_CORRECTION;
        lLastInsert = -Infinity;
        matrix = params.scoringMatrix[sA[i - 1]];   // memoize.

        for (let j = 1; j <= lSeqBLen; j++) {
            tb = TRACE_BACK.MATCH;


            // Delete [i,j] score computation
            //
            //  lMatchArr was not cleared -> contains the computations
            //      from the previous column (i-1).
            //  lDelArr too.
            //
            //  - Gap extension is not computed due to the constant
            //      added to the scoring matrix (see addCentre).


            gapOpenA = lMatchArr[j] + GAP_OPEN;
            gapExtentA = lDelArr[j];
            if (j === lSeqBLen) {
                gapOpenA += GAP_END_CORRECTION;
            }

            if (gapOpenA >= gapExtentA) {
                lDelArr[j] = gapOpenA;
            } else {
                // No change on gap extent: lDelArr[j] += 0;
                tb += TRACE_BACK.DEL;
            }


            // Insert [i,j] score computation
            //  - lPrevMatch is the value that has just been computed
            //      for the match at this column. We can't use lMatchArr
            //      to store it yet, because for perf reason we are
            //      saving space on this vector.
            //  - ExtentB is still the score for the previous insert
            //      value. Due to the addCentre() trick, no need to add
            //      gap extension penalty here.

            gapOpenB = lPrevMatch + GAP_OPEN;
            gapExtentB = lLastInsert;

            if (i === lSeqALen) {
                gapOpenB += GAP_END_CORRECTION; // Terminal gap correction
            }

            if (gapOpenB >= gapExtentB) {
                lLastInsert = gapOpenB;
            } else {
                // No score change on gap extend: lLastInsert += 0;
                tb += TRACE_BACK.INS;
            }


            // Match [i,j] score computation
            lMatch = lMatchArr[j - 1] + matrix[sB[j - 1]];

            lMatchArr[j - 1] = lPrevMatch;  // var rotation

            if (lMatch >= lLastInsert) {
                if (lMatch >= lDelArr[j]) { // match is optimal
                    lPrevMatch = lMatch;
                } else {                    // delete is optimal
                    lPrevMatch = lDelArr[j];
                    tb += TRACE_BACK.MATCH2DEL;
                }

            } else {
                if (lLastInsert >= lDelArr[j]) { // insert is optimal
                    lPrevMatch = lLastInsert;
                    tb += TRACE_BACK.MATCH2INS;
                } else {                         // delete is optimal
                    lPrevMatch = lDelArr[j];
                    tb += TRACE_BACK.MATCH2DEL;
                }
            }

            // Store trace-back bits
            tbIdx = i * lSeqBLen + j;
            isOdd =  tbIdx % 2;
            tbIdx = tbIdx >>> 1;
            if (isOdd) {
                tbM[tbIdx] += tb
            } else {
                tbM[tbIdx] += tb << 4;
            }

        }

        // fix match vector for its last value
        lMatchArr[lSeqBLen] = lPrevMatch;

    }

    // if (DEBUG) Log.add('End DP computation');

    var score = Math.max(lMatch, lLastInsert, lDelArr[lDelArr.length - 1]);

    // Traceback
    // Move backwards starting from the last optimal scoring position
    // in the TB matrix.

    let i = lSeqALen;
    let j = lSeqBLen;
    var lIdx = lSeqALen * lSeqBLen + lSeqBLen;
    const lSeqA: string[] = [];
    const lSeqB: string[] = [];

    // current matrix is either M (0), D (1) or I(2). Let's have a look
    // at the last value to see if the optimum is coming from DEL
    // (bit value 4) or INS (bit value 8) or Match (no bit value)
    var lCurrentMatrix = (tb & 12) >> 2,
        val = 0;

    while ((i > 0) && (j > 0)) {
        lIdx = i * lSeqBLen + j;
        tbIdx = lIdx >>> 1;
        isOdd = lIdx % 2;
        val = tbM[tbIdx];
        val = isOdd ? val & 0b1111 : val >>> 4;
        if (lCurrentMatrix === 0) {
            val = val >> 2;
            if (val === 0) { //-->Match
                i--;
                j--;
                lSeqA.push(seqA.rawSeq[i]);
                lSeqB.push(seqB.rawSeq[j]);
            }
            // other cases --> run the loop once more to enter the
            // following block.
        } else {
            if (lCurrentMatrix === 2) { //-->Ins
                j--;
                lSeqA.push('-');
                lSeqB.push(seqB.rawSeq[j]);
                val = val & 2;
            } else { //1 --> Del
                i--;
                lSeqA.push(seqA.rawSeq[i]);
                lSeqB.push('-');
                val = val & 1;
            }
        }
        lCurrentMatrix = val;
    }

    let lAlignment = [
        lSeqA.reverse().join(''),
        lSeqB.reverse().join('')
    ];

    // if (DEBUG) Log.add('End traceback');

    // Finish sequences edit by appending the remaining symbols
    if (i > 0) {
        lAlignment[1] = '-'.repeat(i) + lAlignment[1];
        lAlignment[0] = seqA.rawSeq.substring(0, i) + lAlignment[0];
    } else if (j > 0) {
        lAlignment[0] = '-'.repeat(j) + lAlignment[0];
        lAlignment[1] = seqB.rawSeq.substring(0, j) + lAlignment[1];
    }

    return {
        alignment: lAlignment,
        score: score
    };
};


/**
 * Profile to sequence Alignment
 *
 * @param {LeafNode} nodeB a node object containing one sequence and its weight
 * @param {InternalNode} nodeA a node object containing a profile and its weight
 * @return {array) array containing aligned sequences
 */
export function MSASeqAlignment(
    nodeA: InternalNode,
    nodeB: LeafNode,
    opt = 0
) {



    const seqB = nodeB.seq,
        //wA = noeudA.weight, //probably unnecessary
        msaA = nodeA.msa,
        lSeqBLen = seqB.rawSeq.length,
        lProfALen = msaA[0].length;

    let lMatch = 0.0,       // match score
        lMatchArr: number[] = [],
        lDelArr  : number[] = [];
    const tbM = new Uint8Array((lProfALen + 1) * (lSeqBLen + 1)),
        lAlignment: string[] = [];

    let lGapOpenA = 0.0,
        lGapOpenB = 0.0,
        lGapExtendA = 0.0,
        lGapExtendB = 0.0,
        tb = 0,
        lLastInsert = 0.0,
        lPrevLastInsert = 0.0,
        lPrevMatch = 0.0,
        lDeletej_1 = 0.0,
        lInserti_1 = 0.0,
        /** Position specific gap open penaly of profile A at pos i-1 */
        lProfAGapOP = 0.0,
        /** Position specific gap close penalty of profile A at pos i-1 */
        lProfAGapCP = 0.0,
        /** Position specific substitution scores of profile A at pos i-1 */
        lProfASScores = [];

    const params = getAlignmentParameters();

    const GAP_OPEN    = params.gapOP;
    const GAP_OPEN_B  = GAP_OPEN / 2;
    const GAP_CLOSE_B = GAP_OPEN / 2;

    const GAP_START_CORRECTION_B = opt & ALIGNOPT.DISABLE_FAVOR_START_GAP
        ? 0
        : - GAP_OPEN / 4;

    const GAP_END_CORRECTION_B = opt & ALIGNOPT.DISABLE_FAVOR_END_GAP
        ? 0
        : - GAP_OPEN / 4;

    // sequence B as a vector of integers
    const sB = seqB.encodedSeq;

    // Convert MSA in node B to profile
    const profA = profileFromMSA(msaA, GAP_OPEN, nodeA.tabWeight, opt);

    // Note that for performance reason, the outer loop iterates on profile A.
    // This allows caching the values for substitution scores used in the inner
    // loop when iterating over seqB values.

    // INITIALIZATION of DP vectors
    lMatchArr[0] = 0;
    lDelArr[0] = -Infinity;
    for (let j = 1; j <= lSeqBLen; j++) {
        lMatchArr[j] = GAP_OPEN_B + GAP_CLOSE_B + GAP_START_CORRECTION_B; //gap open at j=0 + gap close at j-1
        lDelArr[j] = -Infinity;
    }

    //Main DP routine
    // Note that profile's gap penalties are halved at their extremities

    for (let i = 1; i <= lProfALen; i++) {

        lProfAGapOP = profA[i - 1].m_ScoreGapOpen;
        lProfAGapCP = profA[i - 1].m_ScoreGapClose;
        lProfASScores = profA[i - 1].m_AAScores;

            // lPrevMatch is a scalar, used to compute on the same column (iterations
            // over j) the gap open penalty
        lPrevMatch = GAP_START_CORRECTION_B;
        lLastInsert = -Infinity;
        lDelArr[0] = profA[0].m_ScoreGapOpen;
        lPrevLastInsert = GAP_END_CORRECTION_B;

        for (let j = 1; j <= lSeqBLen; j++) {
            tb = TRACE_BACK.MATCH;

            //Delete i,j score computation
            lGapOpenA = lMatchArr[j] + lProfAGapOP; //
            lGapExtendA = lDelArr[j];
            if (j === lSeqBLen) { //terminal penalties are halved
                lGapOpenA -= lProfAGapOP / 2;
            }

            if (lGapOpenA >= lGapExtendA) {
                lDelArr[j] = lGapOpenA;
            } else {
                //Delete[j] = gapExtendA;
                tb += TRACE_BACK.DEL;
            }

            //Insert i,j score computation
            lGapOpenB = lPrevMatch + GAP_OPEN_B;
            lGapExtendB = lPrevLastInsert = lLastInsert;
            if (i === lProfALen) {
                lGapOpenB += GAP_END_CORRECTION_B;
            }

            if (lGapOpenB >= lGapExtendB) {
                lLastInsert = lGapOpenB;
            } else {
                //gap extend, not score change;
                tb += TRACE_BACK.INS;
            }

            //Match i,j score computation
            lDeletej_1 = lDelArr[j] + lProfAGapCP; //it should be prev
            lInserti_1 = lPrevLastInsert + GAP_CLOSE_B;
            lMatch = lMatchArr[j - 1] + lProfASScores[sB[j - 1]];

            if (j === lSeqBLen) {
                lInserti_1 += GAP_END_CORRECTION_B;
            }

            lMatchArr[j - 1] = lPrevMatch;

            if (lMatch >= lInserti_1) {
                if (lMatch >= lDeletej_1) { //match is optimal
                    lPrevMatch = lMatch;
                } else { //delete is optimal
                    lPrevMatch = lDeletej_1;
                    tb += TRACE_BACK.MATCH2DEL;
                }

            } else {
                if (lInserti_1 >= lDeletej_1) { //insert is optimal
                    lPrevMatch = lInserti_1;
                    tb += TRACE_BACK.MATCH2INS;
                } else { //delete is optimal
                    lPrevMatch = lDeletej_1;
                    tb += TRACE_BACK.MATCH2DEL;
                }
            }
            tbM[i * lSeqBLen + j] = tb;
        }
        lMatchArr[lSeqBLen] = lPrevMatch;
    }
    const score = Math.max(lMatch, lLastInsert, lDelArr[lDelArr.length - 2]);

    //traceback
    let i = lProfALen;
    let j = lSeqBLen;
    let idx = lProfALen * lSeqBLen + lSeqBLen,
        k = 0;

    lAlignment.push( seqB.rawSeq , ...msaA);


    let currentMatrix = (tb & 12) >> 2,
        value = 0;
    while ((i > 0) && (j > 0)) {
        idx = i * lSeqBLen + j;

        if (currentMatrix === 0) {
            value = tbM[idx] >> 2;
            if (value === 0) {
                i--;
                j--;
            }
        } else {
            if (currentMatrix === 2) {
                for (k = 1; k < lAlignment.length; k++) {
                    lAlignment[k] = lAlignment[k].substring(0, i) + '-' + lAlignment[k].substring(i);
                }
                value = tbM[idx] & 2; //@@PP same as in _pairwise ?
                j--;

            } else {
                lAlignment[0] = lAlignment[0].substring(0, j) + '-' + lAlignment[0].substring(j);
                value = tbM[idx] & 1;
                i--;
            }
        }
        currentMatrix = value;
    }

    //Finalize sequences
    if (i === 0) {
        var p = '-'.repeat(j);
        for (k = 1; k < lAlignment.length; k++) {
            lAlignment[k] = p + lAlignment[k];
        }
    } else { //j==0
        lAlignment[0] = '-'.repeat(i) + lAlignment[0];
    }
    return {
        alignment: lAlignment,
        score: score
    };
}

/**
 * profile to profile alignment
 * @param {object} noeudA a node object containing a multiple sequence alignment and its weight
 * @param {object} noeudB a node object containing a multiple sequence alignment and its weight
 */

export function MSAMSAAlignment(
    nodeA: InternalNode,
    nodeB: InternalNode,
    opt = 0
) {

    let lMsaA: string[],
        lMsaB: string[],
        lWtB: number[],
        lWtA: number[];
    let lInverted = false;

    //permutation: make B the smallest, to reduce iterations
    if (nodeA.msa.length < nodeB.msa.length) {
        lMsaB = nodeA.msa;
        lMsaA = nodeB.msa;
        lWtB = nodeA.tabWeight;
        lWtA = nodeB.tabWeight;
        lInverted = true;
    } else {
        lMsaA = nodeA.msa;
        lMsaB = nodeB.msa;
        lWtB = nodeB.tabWeight;
        lWtA = nodeA.tabWeight;
    }

    const n = lMsaA[0].length,
        m = lMsaB[0].length,
        tbM = new Uint8Array((n + 1) * (m + 1));

    let lMatchArr = [],
        lDelArr = [],
        lAlignment: string[] = [],
        lMatch      = 0.0,
        lGapOpenA   = 0.0,
        lGapOpenB   = 0.0,
        lGapExtendA = 0.0,
        lGapExtendB = 0.0,
        tb = 0|0,
        lLastInsert = 0.0,
        lPrevLastInsert = 0,
        lPrevMatch  = 0.0,
        lPrevDelete = 0.0,
        lDeletej_1  = 0.0,
        lInserti_1  = 0.0,
        lProfAGapOP = 0.0,
        lProfAGapCP = 0.0,
        lProfAAAScores = [],
        kmax = 0,
        k = 0;

    const params = getAlignmentParameters();

    const GAP_OPEN = params.gapOP;
    const GAP_START_FACTOR = opt & ALIGNOPT.DISABLE_FAVOR_START_GAP ? 1 : 2;
    const GAP_END_FACTOR   = opt & ALIGNOPT.DISABLE_FAVOR_END_GAP ? 1 : 2;

    //convert to profile
    const profB = profileFromMSA(lMsaB, GAP_OPEN, lWtB, opt);
    const profA = profileFromMSA(lMsaA, GAP_OPEN, lWtA, opt);

    const profBGapOPTab = [],
        profBGapCPTab = [];

    for (var i = 0; i < profB.length; i++) {
        profBGapOPTab.push(profB[i].m_ScoreGapOpen);
        profBGapCPTab.push(profB[i].m_ScoreGapClose);
    }


    //INITIALIZATION of DP vectors
    lMatchArr[0] = 0;
    lDelArr[0] = -Infinity;
    for (var j = 1; j <= m; j++) {
        lMatchArr[j] = profBGapOPTab[0] + profBGapCPTab[j - 1] / GAP_START_FACTOR;
        lDelArr[j] = -Infinity;
    }

    const M0 = profA[0].m_ScoreGapOpen;

    //Main DP routine

    for (i = 1; i <= n; i++) {
        lProfAAAScores = profA[i - 1].m_AAScores;
        lProfAGapOP = profA[i - 1].m_ScoreGapOpen;
        lProfAGapCP = profA[i - 1].m_ScoreGapClose;

        //M0 += profAGapEP / 2;

        lPrevMatch = M0 + lProfAGapCP / GAP_START_FACTOR;
        lLastInsert = -Infinity;
        lDelArr[0] = profA[0].m_ScoreGapOpen;

        for (j = 1; j <= m; j++) {
            tb = TRACE_BACK.MATCH;

            //Delete i,j score computation
            lGapOpenA = lMatchArr[j] + lProfAGapOP; //
            lGapExtendA = lPrevDelete = lDelArr[j];
            if (j === m && !(opt & ALIGNOPT.DISABLE_FAVOR_END_GAP)) {
                lGapOpenA -= lProfAGapOP / 2;
            }

            if (lGapOpenA >= lGapExtendA) {
                lDelArr[j] = lGapOpenA;
            } else {
                //Delete[j] = gapExtendA;
                tb += TRACE_BACK.DEL;
            }

            //Insert i,j score computation
            lGapOpenB = lPrevMatch + profBGapOPTab[j - 1];
            lGapExtendB = lPrevLastInsert = lLastInsert;
            if (i === n && !(opt & ALIGNOPT.DISABLE_FAVOR_END_GAP)) {
                lGapOpenB -= profBGapOPTab[j - 1] / 2;
            }

            if (lGapOpenB >= lGapExtendB) {
                lLastInsert = lGapOpenB;
            } else {
                //lastInsert = gapExtendB;
                tb += TRACE_BACK.INS;
            }

            //Match i,j score computation
            // Note: could be a separate function, but runs faster when inlined
            lMatch = lMatchArr[j - 1];
            kmax = profB[j - 1].m_uResidueGroup;
            k = 0;
            while (k < kmax) {
                var resProfNo = profB[j - 1].m_uSortOrder[k];
                lMatch += profB[j - 1].m_wCounts[resProfNo] * lProfAAAScores[resProfNo];
                k++;
            }

            //match = Match[j - 1] + _utils.sumOfPairsScorePP3(profAAAScores, profB[j - 1]);
            lDeletej_1 = lDelArr[j] + lProfAGapCP;
            lInserti_1 = lPrevLastInsert + profBGapCPTab[j - 1];

            if (j === 1 && !(opt & ALIGNOPT.DISABLE_FAVOR_END_GAP)) { //terminal penalties are halved
                lDeletej_1 -= lProfAGapCP / 2;
            }
            if ((j === m) && !(opt & ALIGNOPT.DISABLE_FAVOR_END_GAP)) {
                lDeletej_1 -= lProfAGapCP / 2;
                lInserti_1 -= profBGapCPTab[m - 1] / 2;
            }

            lMatchArr[j - 1] = lPrevMatch;

            if (lMatch >= lInserti_1) {
                if (lMatch >= lDeletej_1) { //match is optimal
                    lPrevMatch = lMatch;
                } else { //delete is optimal
                    lPrevMatch = lDeletej_1;
                    tb += TRACE_BACK.MATCH2DEL;
                }

            } else {
                if (lInserti_1 >= lDeletej_1) { //insert is optimal
                    lPrevMatch = lInserti_1;
                    tb += TRACE_BACK.MATCH2INS;
                } else { //delete is optimal
                    lPrevMatch = lDeletej_1;
                    tb += TRACE_BACK.MATCH2DEL;
                }
            }
            tbM[i * m + j] = tb;
        }
        lMatchArr[m] = lPrevMatch;
    }
    const score = Math.max(lMatch, lLastInsert, lDelArr[j - 1]);

    //traceback
    i = n;
    j = m;
    k = 0;
    let idx = n * m;

    lAlignment = lMsaA.concat(lMsaB);

    let currentMatrix = (tb & 12) >> 2,
        value = 0;
    while ((i > 0) && (j > 0)) {
        idx = i * m + j;
        if (currentMatrix === 0) {
            value = tbM[idx] >> 2;
            if (value === 0) {
                i--;
                j--;
            }
        } else {
            if (currentMatrix === 2) {
                for (k = 0; k < lMsaA.length; k++) {
                    lAlignment[k] = lAlignment[k].substring(0, i) + '-' + lAlignment[k].substring(i);
                }
                value = tbM[idx] & 2;
                j--;
            } else {
                for (k = lMsaA.length; k < lAlignment.length; k++) {
                    lAlignment[k] = lAlignment[k].substring(0, j) + '-' + lAlignment[k].substring(j);
                }
                value = tbM[idx] & 1;
                i--;
            }
        }
        currentMatrix = value;
    }

    //finalize alignments
    let lPadding: string;
    if (j === 0) {
        lPadding = '-'.repeat(i);
        for (k = lMsaA.length; k < lAlignment.length; k++) {
            lAlignment[k] = lPadding + lAlignment[k];
        }
    } else { //i==0
        lPadding = '-'.repeat(j);
        for (k = 0; k < lMsaA.length; k++) {
            lAlignment[k] = lPadding + lAlignment[k];
        }
    }

    if (lInverted) lAlignment = [
        ...lAlignment.slice(lMsaA.length),
        ...lAlignment.slice(0, lMsaA.length)
    ];

    return {
        alignment: lAlignment,
        score: score
    };
}
