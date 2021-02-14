/**
 * @file profile
 * @author Paul Pillot <paul.pillot@libmol.org>
 * @license MIT
 * @copyright 2020
 */

import { ALIGNOPT } from '../align/align';
import { SEQUENCE_TYPE, TAlignmentParam } from '../align/params';
import { aaToNum, nucToNum } from './sequence';

export function sumOfPairsScorePP(profA: Profile, profB: Profile) {
    var score = 0,
        resProfNo = 0;

    for (var i = 0; i < profB.m_uResidueGroup; i++) {
        resProfNo = profB.m_uSortOrder[i];
        score += profB.m_fcCounts[resProfNo] * profA.m_AAScores[resProfNo];
    }
    score /= profB.nbSeq;
    //if (score<0) console.log(resNo,prof.m_uSortOrder, score)
    return score;
}

export function sumOfPairsScoreSP(resA: number, prof: Profile) {
    var score = 0,
        resProfNo = 0;

    score = prof.m_AAScores[resA];
    return score;
}

export type Profile = {
    m_bAllGaps: boolean;
    /** Amino Acids by order of frequency */
    m_uSortOrder: number[];
    /** Am. Ac. counts at this position */
    m_fcCounts: number[];
    m_LL: number;
    m_LG: number;
    m_GL: number;
    m_GG: number;
    /** Scores for each Am. Ac at this position */
    m_AAScores: number[];
    /** number of different residues */
    m_uResidueGroup: number;
    /** Frequence of occupency (non indels) */
    m_fOcc: number;
    /** Frequence of gap start */
    m_fcStartOcc: number;
    /** Frequence of gap close */
    m_fcEndOcc: number;
    /** Score for starting gap */
    m_scoreGapOpen: number;
    /** Score for closing gap */
    m_scoreGapClose: number;
    /** Number of sequences */
    nbSeq: number;
};
/*
 * Définition de l'objet ProfPos()
 */
class ProfPos {
    m_bAllGaps = false;
    /** Residues sorted in decreasing frequency order */
    m_uSortOrder: number[] = [];
    /** Residues counts */
    m_fcCounts = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
    /** Weighted contribution per residue */
    m_wCounts = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
    // m_LL = 1;
    // m_LG = 0;
    // m_GL = 0;
    // m_GG = 0;
    /** Substitution score computed at this position */
    m_AAScores = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
    /** Residue diversity (count per residue) */
    m_uResidueGroup = 0;
    /** Occupancy (non-indels) */
    m_fOcc = 1;
    /** Frequency of opening gaps */
    m_fcStartOcc = 0;
    /** Frequency of closing gaps */
    m_fcEndOcc = 0;
    /** Gap Open penalty computed for this position */
    m_ScoreGapOpen = 0;
    /** Gap Close penalty computed for this position */
    m_ScoreGapClose = 0;
    /** Sequence count in this profile */
    nbSeq = 0;
}

/**
 * Compute Profile from MSA
 * Returns an array of object with at each position of the profile the
 * corresponding am. ac. counts, nb of indels for gap closing/opening/extension
 * @param {string[]} pMSA: aligned sequences as an array of strings with `-` as
 *                          indels.
 * @param {number}  pGapO: gap opening penalty
 * @param {number[]} pWeights : weight of each sequence
 */
export function profileFromMSA (pMSA: string[], pGapO: number, pWeights: number[], params: TAlignmentParam,
    opt: number = 0) {
    /** alignment length (columns count) */
    const l = pMSA[0].length;

    let aa = '', // amino acid name (1 letter code)
        na = 0,  // amino acid # (in lexical order)
        nbResMax = 0,
        tabCodeResNum: Uint8Array;

    var prof: ProfPos[] = new Array(l); // profile
    if (params.type === SEQUENCE_TYPE.NUCLEIC) {
        nbResMax = 4;
        tabCodeResNum = nucToNum;
    } else {
        nbResMax = 20;
        tabCodeResNum = aaToNum;
    }

    //uHydrophobicRunLength = 0, // hydrophobic window (to avoid introducing gaps within)
    let w = 0, //sequence weight
        totalW = pWeights.reduce((a, b) =>  a + b, 0);

    const lGapO = pGapO || params.gapOP;

    for (let col = 0; col < l; col++) {
        prof[col] = new ProfPos();
        let fGap = 0;

        for (let numSeq = 0, max = pMSA.length; numSeq < max; numSeq++) {

            prof[col].nbSeq = max;
            w = pWeights[numSeq] / totalW; //sequences individual weights are normalized to the profile total weights
            aa = pMSA[numSeq][col];
            if (aa === '-') {
                fGap += w;

                if (col === 0) {
                    prof[col].m_fcStartOcc += w;
                } else if (pMSA[numSeq][col - 1] !== '-') {
                    prof[col].m_fcStartOcc += w;
                }

                if (col === l - 1) {
                    prof[col].m_fcEndOcc += w;
                } else if (pMSA[numSeq][col + 1] !== '-') {
                    prof[col].m_fcEndOcc += w;
                }

            } else {
                na = tabCodeResNum[aa.charCodeAt(0)];
                prof[col].m_fcCounts[na]++;
                prof[col].m_wCounts[na] += w;
                if (prof[col].m_fcCounts[na] === 1) {
                    prof[col].m_uSortOrder.push(na);
                    prof[col].m_uResidueGroup++;
                }
            }
        }
        prof[col].m_fOcc = totalW - fGap;
        prof[col].m_ScoreGapOpen = lGapO / 2 * (1 - prof[col].m_fcStartOcc);
        prof[col].m_ScoreGapClose = lGapO / 2 * (1 - prof[col].m_fcEndOcc);

        for (var aaNum = 0; aaNum < nbResMax; aaNum++) {
            for (var i = 0; i < prof[col].m_uResidueGroup; i++) {
                var resProfNo = prof[col].m_uSortOrder[i];
                prof[col].m_AAScores[aaNum] += prof[col].m_wCounts[resProfNo] * params.scoringMatrix[aaNum][resProfNo];
            }
        }
    }


    //  Corrections for end/begin gap penalties
    if (opt ^ ALIGNOPT.DISABLE_FAVOR_START_GAP) {
        prof[0]  .m_ScoreGapOpen  /= 2;
        prof[l-1].m_ScoreGapOpen  /= 2;
    }
    if (opt ^ ALIGNOPT.DISABLE_FAVOR_END_GAP) {
        prof[0]  .m_ScoreGapClose /= 2;
        prof[l-1].m_ScoreGapClose /= 2;
    }


    return prof;
}
