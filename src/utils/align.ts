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
import { profileFromMSA } from "./profile";
import { TSequence } from "./sequence";
import { InternalNode, LeafNode } from "./tree";

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
 * @param {number[]} idA
 * @param {number[]} idB
 * @param {number} [opt=0]
 * @returns
 */
export function pairwiseAlignment (
    seqA: TSequence,
    seqB: TSequence,
    idA: number[],
    idB: number[],
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
        for (var j = 1; j <= lSeqBLen; j++) {
            lMatchArr[j] = GAP_OPEN + GAP_START_CORRECTION; // No gap extension penalty
            lDelArr[j] = -Infinity;
        }

        // Main DP routine

        for (var i = 1; i <= lSeqALen; i++) {

            lPrevMatch = GAP_OPEN + GAP_START_CORRECTION;
            lLastInsert = -Infinity;
            matrix = params.scoringMatrix[sA[i - 1]];   // memoize.

            for (j = 1; j <= lSeqBLen; j++) {
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

        var score = Math.max(lMatch, lLastInsert, lDelArr[j - 1]);

        // Traceback
        // Move backwards starting from the last optimal scoring position
        // in the TB matrix.

        i = lSeqALen;
        j = lSeqBLen;
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
            tSeqNames: idA.concat(idB),
            score: score
        };
};


/**
 * Profile to sequence Alignment
 *
 * @param {LeafNode} nodeA a node object containing one sequence and its weight
 * @param {InternalNode} nodeB a node object containing a multiple sequence alignment and its weight
 * @return {array) array containing aligned sequences
 */
export function MSASeqAlignment(
    nodeA: LeafNode,
    nodeB: InternalNode,
    tSeqANames: number[],
    tSeqBNames: number[]
) {



    var seqA = nodeA.seq,
        //wA = noeudA.weight, //probably unnecessary
        msaB = nodeB.msa,
        m = seqA.rawSeq.length,
        n = msaB[0].length,
        lMatch = 0.0,       // match score
        lMatchArr = [],
        lDelArr = [],
        tbM = new Uint8Array(n * (m + 1)),
        lAlignment: string[] = [],
        gapOpenA = 0.0,
        gapOpenB = 0.0,
        gapExtendA = 0.0,
        gapExtendB = 0.0,
        tb = 0,
        lastInsert = 0.0,
        prevLastInsert = 0.0,
        prevMatch = 0.0,
        prevDelete = 0.0,
        deletej_1 = 0.0,
        inserti_1 = 0.0,
        //M0 = 0.0,
        //profBGapEP = 0.0,
        profBGapOP = 0.0, //gap open penalty for profB at pos i-1
        profBGapCP = 0.0, //gap close penalty for profB at pos i-1
        profBAAScores = []; //substitution scores for profB at pos i-1

    const params = getAlignmentParameters();

    var gapOP = params.gapOP;
    var seqAGapOP = gapOP / 2;
    var seqAGapCP = gapOP / 2;

    // sequence A as number vector
    var sA = seqA.encodedSeq;

    // MSA B as profile
    var profB = profileFromMSA(msaB, gapOP, nodeB.tabWeight);

    //in this alignment for performance reasons, profB is explored from i=0 to n
    //and seqA is explored from j=0 to m

    //INITIALISATION
    lMatchArr[0] = 0; //_utils.sumOfPairsScoreSP(sA[0], profB[0]);
    lDelArr[0] = -Infinity;
    for (var j = 1; j <= m; j++) {
        lMatchArr[j] = (seqAGapOP + seqAGapCP) / 2; //gap open at j=0 + gap close at j-1
        lDelArr[j] = -Infinity;
    }

    //RECURSION
    //$log.debug('seqA: ', tSeqANames, ' penalités gapEP: ', gapEP, ' gapOP: ', gapOP);
    //M0 = Match[0] + profB[0].m_ScoreGapOpen; //est deux fois plus faible qu'une pénalité d'ouverture normale (début de séquence)

    for (var i = 1; i <= n; i++) {

        profBGapOP = profB[i - 1].m_ScoreGapOpen;
        profBGapCP = profB[i - 1].m_ScoreGapClose;
        profBAAScores = profB[i - 1].m_AAScores;

        //M0 += profBGapEP / 2;

        prevMatch = 0; //(i===1) ? 0 : profB[0].m_ScoreGapOpen + profBGapCP/2 ;
        lastInsert = -Infinity;
        lDelArr[0] = profB[0].m_ScoreGapOpen;

        for (j = 1; j <= m; j++) {
            tb = 0;

            //Delete i,j score computation
            gapOpenA = lMatchArr[j] + profBGapOP; //
            gapExtendA = prevDelete = lDelArr[j];
            if (j === m) { //terminal penalties are halved
                gapOpenA -= profBGapOP / 2;
            }

            if (gapOpenA >= gapExtendA) {
                lDelArr[j] = gapOpenA;
            } else {
                //Delete[j] = gapExtendA;
                tb += 1;
            }

            //Insert i,j score computation
            gapOpenB = prevMatch + seqAGapOP;
            gapExtendB = prevLastInsert = lastInsert;
            if (i === n) {
                gapOpenB -= seqAGapOP / 2;
            }

            if (gapOpenB >= gapExtendB) {
                lastInsert = gapOpenB;
            } else {
                //lastInsert = gapExtendB;
                tb += 2;
            }

            //Match i,j score computation
            deletej_1 = lDelArr[j] + profBGapCP; //it should be prev
            inserti_1 = prevLastInsert + seqAGapCP;
            lMatch = lMatchArr[j - 1] + profBAAScores[sA[j - 1]];

            if (j === 1) { //terminal penalties are halved
                deletej_1 -= profBGapCP / 2;
            }
            if (i === 1) {
                inserti_1 -= seqAGapCP / 2;
            }
            if ((j === m) && (i === n)) {
                deletej_1 -= profBGapCP / 2;
                inserti_1 -= seqAGapCP / 2;
            }

            lMatchArr[j - 1] = prevMatch;

            if (lMatch >= inserti_1) {
                if (lMatch >= deletej_1) { //match is optimal
                    prevMatch = lMatch;
                } else { //delete is optimal
                    prevMatch = deletej_1;
                    tb += 4;
                }

            } else {
                if (inserti_1 >= deletej_1) { //insert is optimal
                    prevMatch = inserti_1;
                    tb += 8;
                } else { //delete is optimal
                    prevMatch = deletej_1;
                    tb += 4;
                }
            }
            tbM[i * m + j] = tb;
        }
        lMatchArr[m] = prevMatch;
    }
    var score = Math.max(lMatch, lastInsert, lDelArr[j - 1]);

    //traceback
    i = n;
    j = m;
    var indice = n * m + m,
        k = 0,
        tSeqNames = [];

    lAlignment.push( seqA.rawSeq , ...msaB);
    tSeqNames = [...tSeqANames, ...tSeqBNames];

    var currentMatrix = (tb & 12) >> 2,
        value = 0;
    while ((i > 0) && (j > 0)) {
        indice = i * m + j;

        if (currentMatrix === 0) {
            value = tbM[indice] >> 2;
            if (value === 0) {
                i--;
                j--;
            }
        } else {
            if (currentMatrix === 2) {
                for (k = 1; k < lAlignment.length; k++) {
                    lAlignment[k] = lAlignment[k].substring(0, i) + '-' + lAlignment[k].substring(i);
                }
                value = tbM[indice] & 2; //@@PP same as in _pairwise ?
                j--;

            } else {
                lAlignment[0] = lAlignment[0].substring(0, j) + '-' + lAlignment[0].substring(j);
                value = tbM[indice] & 1;
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
        tSeqNames: tSeqNames,
        score: score
    };
}

/**
 * profile to profile alignment
 * @param {object} noeudA a node object containing a multiple sequence alignment and its weight
 * @param {object} noeudB a node object containing a multiple sequence alignment and its weight
 * @param {integer} gE optional gap extension penalty parameter
 * @param {integer} gO optional gap open penalty parameter
 * @return {array) array containing aligned sequences
 */

export function MSAMSAAlignment(
    nodeA: InternalNode,
    nodeB: InternalNode,
    tSeqANames: number[],
    tSeqBNames: number[]
) {

    var msaA, msaB, wB, wA, tSeqNames;

    if (nodeA.msa.length < nodeB.msa.length) {
        //permutation: make B the smallest, to reduce iterations
        msaB = nodeA.msa;
        msaA = nodeB.msa;
        wB = nodeA.tabWeight;
        wA = nodeB.tabWeight;
        tSeqNames = tSeqBNames.concat(tSeqANames);
    } else {
        msaA = nodeA.msa;
        msaB = nodeB.msa;
        wB = nodeB.tabWeight;
        wA = nodeA.tabWeight;
        tSeqNames = tSeqANames.concat(tSeqBNames);
    }

    var n = msaA[0].length,
        m = msaB[0].length,
        Match = [],
        Delete = [],
        tbM = new Uint8Array(n * (m + 1)),
        alignement = [],
        match = 0,
        gapOpenA = 0,
        gapOpenB = 0,
        gapExtendA = 0,
        gapExtendB = 0,
        tb = 0,
        lastInsert = 0,
        prevLastInsert = 0,
        M0 = 0,
        prevMatch = 0,
        prevDelete = 0,
        deletej_1 = 0,
        inserti_1 = 0,
        profAGapOP = 0,
        profAGapCP = 0,
        profAAAScores = [];

    const _params = getAlignmentParameters();

    var gapOP = _params.gapOP;

    //convert to profile
    var profB = profileFromMSA(msaB, gapOP, wB);
    var profA = profileFromMSA(msaA, gapOP, wA);

    var profBGapOPTab = [],
        profBGapCPTab = [];

    for (var i = 0; i < profB.length; i++) {
        profBGapOPTab.push(profB[i].m_ScoreGapOpen);
        profBGapCPTab.push(profB[i].m_ScoreGapClose);
    }
    /*var profBWCountsTab = [];
    for (var i = 0; i < profB.length; i++) {
        profBWCountsTab.push(profB[i].m_wCounts);
    }*/
    //$log.debug(profBWCountsTab);
    //remplissage de la première ligne et de la première colonne de la matrice
    //INITIALISATION
    Match[0] = 0; //_utils.sumOfPairsScorePP(profA[0], profB[0]);
    Delete[0] = -Infinity;
    for (var j = 1; j <= m; j++) {
        Match[j] = profBGapOPTab[0] + profBGapCPTab[j - 1] / 2;
        Delete[j] = -Infinity;
    }

    M0 = profA[0].m_ScoreGapOpen;
    //remplissage des matrices
    //RECURSION
    for (i = 1; i <= n; i++) {
        profAAAScores = profA[i - 1].m_AAScores;
        profAGapOP = profA[i - 1].m_ScoreGapOpen;
        profAGapCP = profA[i - 1].m_ScoreGapClose;

        //M0 += profAGapEP / 2;

        prevMatch = M0 + profAGapCP / 2;
        lastInsert = -Infinity;
        Delete[0] = profA[0].m_ScoreGapOpen;

        for (j = 1; j <= m; j++) {
            tb = 0;

            //Delete i,j score computation
            gapOpenA = Match[j] + profAGapOP; //
            gapExtendA = prevDelete = Delete[j];
            if (j === m) {
                gapOpenA -= profAGapOP / 2;
            }

            if (gapOpenA >= gapExtendA) {
                Delete[j] = gapOpenA;
            } else {
                //Delete[j] = gapExtendA;
                tb += 1;
            }

            //Insert i,j score computation
            gapOpenB = prevMatch + profBGapOPTab[j - 1];
            gapExtendB = prevLastInsert = lastInsert;
            if (i === n) {
                gapOpenB -= profBGapOPTab[j - 1] / 2;
            }

            if (gapOpenB >= gapExtendB) {
                lastInsert = gapOpenB;
            } else {
                //lastInsert = gapExtendB;
                tb += 2;
            }

            //Match i,j score computation
            match = Match[j - 1];
            var kmax = profB[j - 1].m_uResidueGroup,
                k = 0;
            while (k < kmax) {
                var resProfNo = profB[j - 1].m_uSortOrder[k];
                match += profB[j - 1].m_wCounts[resProfNo] * profAAAScores[resProfNo];
                k++;
            }

            //match = Match[j - 1] + _utils.sumOfPairsScorePP3(profAAAScores, profB[j - 1]);
            deletej_1 = Delete[j] + profAGapCP;
            inserti_1 = prevLastInsert + profBGapCPTab[j - 1];

            if (j === 1) { //terminal penalties are halved
                deletej_1 -= profAGapCP / 2;
            }
            if (i === 1) {
                inserti_1 -= profBGapCPTab[0] / 2;
            }
            if ((j === m) && (i === n)) {
                deletej_1 -= profAGapCP / 2;
                inserti_1 -= profBGapCPTab[m - 1] / 2;
            }

            Match[j - 1] = prevMatch;

            if (match >= inserti_1) {
                if (match >= deletej_1) { //match is optimal
                    prevMatch = match;
                } else { //delete is optimal
                    prevMatch = deletej_1;
                    tb += 4;
                }

            } else {
                if (inserti_1 >= deletej_1) { //insert is optimal
                    prevMatch = inserti_1;
                    tb += 8;
                } else { //delete is optimal
                    prevMatch = deletej_1;
                    tb += 4;
                }
            }
            tbM[i * m + j] = tb;
        }
        Match[m] = prevMatch;
    }
    var score = Math.max(match, lastInsert, Delete[j - 1]);

    //traceback
    i = n;
    j = m;
    k = 0;
    var indice = n * m,
        k = 0,
        p = 0;
    alignement = msaA.concat(msaB);

    var matriceActive = (tb & 12) >> 2,
        value = 0;
    while ((i > 0) && (j > 0)) {
        indice = i * m + j;
        if (matriceActive === 0) {
            value = tbM[indice] >> 2;
            if (value === 0) {
                i--;
                j--;
            }
        } else {
            if (matriceActive === 2) {
                for (k = 0; k < msaA.length; k++) {
                    alignement[k] = alignement[k].substring(0, i) + '-' + alignement[k].substring(i);
                }
                value = tbM[indice] & 2;
                j--;
            } else {
                for (k = msaA.length; k < alignement.length; k++) {
                    alignement[k] = alignement[k].substring(0, j) + '-' + alignement[k].substring(j);
                }
                value = tbM[indice] & 1;
                i--;
            }
        }
        matriceActive = value;
    }
    //fin des profils par ajout direct de la fin de séquence
    let lPadding: string;
    if (j === 0) {
        lPadding = '-'.repeat(i);
        for (k = msaA.length; k < alignement.length; k++) {
            alignement[k] = lPadding + alignement[k];
        }
    } else { //i==0
        lPadding = '-'.repeat(j);
        for (k = 0; k < msaA.length; k++) {
            alignement[k] = lPadding + alignement[k];
        }
    }

    return {
        alignment: alignement,
        tSeqNames: tSeqNames,
        score: score
    };
}
