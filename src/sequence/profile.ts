/**
 * @file profile
 * @author Paul Pillot <paul.pillot@libmol.org>
 * @license MIT
 * @copyright 2020
 */

import { ALIGNOPT } from '../align/align';
import { SEQUENCE_TYPE, TAlignmentParam } from '../align/params';
import BitArray from '../utils/bitarray';
import { aaToNum, nucToNum } from './sequence';


/*
 * ProfPos definition ()
 */
class ProfPos {
    m_bAllGaps: BitArray;
    /** Residues present at this position unordered */
    m_uSortOrder: Uint8Array;
    /** Residues counts */
    m_fcCounts: Uint8Array;
    /** Weighted contribution per residue */
    m_wCounts: Float32Array;
    // m_LL = 1;
    // m_LG = 0;
    // m_GL = 0;
    // m_GG = 0;
    /** Substitution score computed at this position */
    m_AAScores: Float32Array;
    /** Residue diversity (count per residue) */
    m_uResidueGroup: Uint8Array;
    /** Occupancy (non-indels) */
    m_fOcc: Float32Array;
    /** Frequency of opening gaps */
    m_fcStartOcc: Float32Array;
    /** Frequency of closing gaps */
    m_fcEndOcc: Float32Array;
    /** Gap Open penalty computed for this position */
    m_ScoreGapOpen: Float32Array;
    /** Gap Close penalty computed for this position */
    m_ScoreGapClose: Float32Array;
    /** Sequence count in this profile */
    nbSeq = 0;
    alphaSize: number;
    length = 0;

    constructor (pSize: number, pAlphaSize: number, pNbSeq: number) {
        this.length = pSize;
        this.m_bAllGaps = new BitArray(pSize);
        this.m_uSortOrder = new Uint8Array(pAlphaSize * pSize);
        this.m_fcCounts = new Uint8Array(pAlphaSize * pSize);
        this.m_wCounts = new Float32Array(pAlphaSize * pSize);
        this.m_AAScores = new Float32Array(pAlphaSize * pSize);
        this.m_uResidueGroup = new Uint8Array(pSize);
        this.m_fOcc = new Float32Array(pSize);
        this.m_fcStartOcc = new Float32Array(pSize);
        this.m_fcEndOcc = new Float32Array(pSize);
        this.m_ScoreGapOpen = new Float32Array(pSize);
        this.m_ScoreGapClose = new Float32Array(pSize);
        this.nbSeq = pNbSeq;
        this.alphaSize = pAlphaSize;
    }

    getProxy (id = 0) {
        return new ProfProxy(this, id);
    }
}

class ProfProxy {

    prof: ProfPos;
    idx: number = 0;
    alphaSize = 4;
    offset = 0;

    constructor(pProf: ProfPos, pIdx = 0) {
        this.prof = pProf;
        this.idx = pIdx;
        this.alphaSize = pProf.alphaSize;
        this.offset = this.idx * this.alphaSize;
    }

    setProxy (pIdx: number) {
        this.idx = pIdx;
        this.offset = pIdx * this.alphaSize;
    }

    get m_bAllGaps () {
        return this.prof.m_bAllGaps.get(this.idx);
    }
    set m_bAllGaps (val: boolean) {
        if (val) this.prof.m_bAllGaps.set(this.idx)
        else this.prof.m_bAllGaps.clear(this.idx);
    }

    get m_uSortOrder () {
        return this.prof.m_uSortOrder.subarray(this.offset, this.offset + this.alphaSize);
    }
    set m_uSortOrder (val: Uint8Array) {
        this.prof.m_uSortOrder.set(val, this.offset);
    }

    get m_uResidueGroup () {
        return this.prof.m_uResidueGroup[this.idx];
    }
    set m_uResidueGroup (val: number) {
        this.prof.m_uResidueGroup[this.idx] = val;
    }


    get m_fcCounts () {
        return this.prof.m_fcCounts.subarray(this.offset, this.offset + this.alphaSize);
    }
    set m_fcCounts (val: Uint8Array) {
        this.prof.m_fcCounts.set(val, this.offset);
    }


    get m_wCounts () {
        return this.prof.m_wCounts.subarray(this.offset, this.offset + this.alphaSize);
    }
    set m_wCounts (val: Float32Array) {
        this.prof.m_wCounts.set(val, this.offset);
    }

    get m_AAScores () {
        return this.prof.m_AAScores.subarray(this.offset, this.offset + this.alphaSize);
    }
    set m_AAScores (val: Float32Array) {
        this.prof.m_AAScores.set(val, this.offset);
    }

    get m_fOcc () {
        return this.prof.m_fOcc[this.idx];
    }
    set m_fOcc (val: number) {
        this.prof.m_fOcc[this.idx] = val;
    }

    get m_fcStartOcc () {
        return this.prof.m_fcStartOcc[this.idx];
    }
    set m_fcStartOcc (val: number) {
        this.prof.m_fcStartOcc[this.idx] = val;
    }

    get m_fcEndOcc () {
        return this.prof.m_fcEndOcc[this.idx];
    }
    set m_fcEndOcc (val: number) {
        this.prof.m_fcEndOcc[this.idx] = val;
    }

    get m_ScoreGapOpen () {
        return this.prof.m_ScoreGapOpen[this.idx];
    }
    set m_ScoreGapOpen (val: number) {
        this.prof.m_ScoreGapOpen[this.idx] = val;
    }

    get m_ScoreGapClose () {
        return this.prof.m_ScoreGapClose[this.idx];
    }
    set m_ScoreGapClose (val: number) {
        this.prof.m_ScoreGapClose[this.idx] = val;
    }
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
export function profileFromMSA (pMSA: string[], pWeights: number[], params: TAlignmentParam,
    opt: number = 0) {
    /** alignment length (columns count) */
    const l = pMSA[0].length;

    let aa = '', // amino acid name (1 letter code)
        na = 0,  // amino acid # (in lexical order)
        tabCodeResNum: Uint8Array;
    if (params.type === SEQUENCE_TYPE.NUCLEIC) {
        tabCodeResNum = nucToNum;
    } else {
        tabCodeResNum = aaToNum;
    }

    const lProf = new ProfPos(l, params.abSize, pMSA.length);

    //uHydrophobicRunLength = 0, // hydrophobic window (to avoid introducing gaps within)
    let w = 0, //sequence weight
        totalW = pWeights.reduce((a, b) =>  a + b, 0);

    const lGapO = params.gapOP;
    let lProx = lProf.getProxy();

    for (let col = 0; col < l; col++) {
        lProx.setProxy(col);
        let fGap = 0;

        for (let numSeq = 0, max = pMSA.length; numSeq < max; numSeq++) {

            w = pWeights[numSeq] / totalW; //sequences individual weights are normalized to the profile total weights
            aa = pMSA[numSeq][col];
            if (aa === '-') {
                fGap += w;

                if (col === 0 || pMSA[numSeq][col - 1] !== '-') {
                    lProx.m_fcStartOcc += w;
                }

                if (col === l - 1 || pMSA[numSeq][col + 1] !== '-') {
                    lProx.m_fcEndOcc += w;
                }

            } else {
                na = tabCodeResNum[aa.charCodeAt(0)];
                lProx.m_fcCounts[na]++;
                lProx.m_wCounts[na] += w;
                if (lProx.m_fcCounts[na] === 1) {
                    lProx.m_uSortOrder[lProx.m_uResidueGroup] = na;
                    lProx.m_uResidueGroup++;
                }
            }
        }
        lProx.m_fOcc = totalW - fGap;
        lProx.m_ScoreGapOpen = lGapO / 2 * (1 - lProx.m_fcStartOcc);
        lProx.m_ScoreGapClose = lGapO / 2 * (1 - lProx.m_fcEndOcc);

            // Cache substitution scores at this position

        let imax = lProx.m_uResidueGroup;
        let lSortOrder = lProx.m_uSortOrder;
        let lScore = 0, lResProfNb = 0;
        for (let aaNum = 0; aaNum < params.abSize; aaNum++) {
            lScore = 0;
            for (let i = 0; i < imax; i++) {
                lResProfNb = lSortOrder[i];
                lScore += lProx.m_wCounts[lResProfNb] * params.scoringMatrix[aaNum][lResProfNb];
            }
            lProx.m_AAScores[aaNum] = lScore;
        }
    }


    //  Corrections for end/begin gap penalties
    if (opt ^ ALIGNOPT.DISABLE_FAVOR_START_GAP) {
        lProf.m_ScoreGapOpen[0]  /= 2;
        lProf.m_ScoreGapOpen[l-1]  /= 2;
    }
    if (opt ^ ALIGNOPT.DISABLE_FAVOR_END_GAP) {
        lProf.m_ScoreGapClose[0] /= 2;
        lProf.m_ScoreGapClose[l-1] /= 2;
    }


    return lProf;
}
