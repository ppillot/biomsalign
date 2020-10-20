/**
 * @file params
 * @author Paul Pillot <paul.pillot@libmol.org>
 * @license MIT
 * @copyright 2020
 */

import { SEQUENCE_TYPE } from './sequence';

let typeSeq: SEQUENCE_TYPE = SEQUENCE_TYPE.PROTEIN;

let scoringMatrix: number[][] = [];

let gapOP = 0;

let gapEP = 0;

const BLOSUM62 = {
    matrix: [
        //    A   C   D   E   F   G   H   I   K   L   M   N   P   Q   R   S   T   V   W   Y
        [4, 0, -2, -1, -2, 0, -2, -1, -1, -1, -1, -2, -1, -1, -1, 1, 0, 0, -3, -2, 0], // A
        [0, 9, -3, -4, -2, -3, -3, -1, -3, -1, -1, -3, -3, -3, -3, -1, -1, -1, -2, -2, 0], // C
        [-2, -3, 6, 2, -3, -1, -1, -3, -1, -4, -3, 1, -1, 0, -2, 0, -1, -3, -4, -3, 0], // D
        [-1, -4, 2, 5, -3, -2, 0, -3, 1, -3, -2, 0, -1, 2, 0, 0, -1, -2, -3, -2, 0], // E
        [-2, -2, -3, -3, 6, -3, -1, 0, -3, 0, 0, -3, -4, -3, -3, -2, -2, -1, 1, 3, 0], // F
        [0, -3, -1, -2, -3, 6, -2, -4, -2, -4, -3, 0, -2, -2, -2, 0, -2, -3, -2, -3, 0], // G
        [-2, -3, -1, 0, -1, -2, 8, -3, -1, -3, -2, 1, -2, 0, 0, -1, -2, -3, -2, 2, 0], // H
        [-1, -1, -3, -3, 0, -4, -3, 4, -3, 2, 1, -3, -3, -3, -3, -2, -1, 3, -3, -1, 0], // I
        [-1, -3, -1, 1, -3, -2, -1, -3, 5, -2, -1, 0, -1, 1, 2, 0, -1, -2, -3, -2, 0], // K
        [-1, -1, -4, -3, 0, -4, -3, 2, -2, 4, 2, -3, -3, -2, -2, -2, -1, 1, -2, -1, 0], // L
        [-1, -1, -3, -2, 0, -3, -2, 1, -1, 2, 5, -2, -2, 0, -1, -1, -1, 1, -1, -1, 0], // M
        [-2, -3, 1, 0, -3, 0, 1, -3, 0, -3, -2, 6, -2, 0, 0, 1, 0, -3, -4, -2, 0], // N
        [-1, -3, -1, -1, -4, -2, -2, -3, -1, -3, -2, -2, 7, -1, -2, -1, -1, -2, -4, -3, 0], // P
        [-1, -3, 0, 2, -3, -2, 0, -3, 1, -2, 0, 0, -1, 5, 1, 0, -1, -2, -2, -1, 0], // Q
        [-1, -3, -2, 0, -3, -2, 0, -3, 2, -2, -1, 0, -2, 1, 5, -1, -1, -3, -3, -2, 0], // R
        [1, -1, 0, 0, -2, 0, -1, -2, 0, -2, -1, 1, -1, 0, -1, 4, 1, -2, -3, -2, 0], // S
        [0, -1, -1, -1, -2, -2, -2, -1, -1, -1, -1, 0, -1, -1, -1, 1, 5, 0, -2, -2, 0], // T
        [0, -1, -3, -2, -1, -3, -3, 3, -2, 1, 1, -3, -2, -2, -3, -2, 0, 4, -3, -1, 0], // V
        [-3, -2, -4, -3, 1, -2, -2, -3, -3, -2, -1, -4, -4, -2, -3, -3, -2, -3, 11, 2, 0], // W
        [-2, -2, -3, -2, 3, -3, 2, -1, -2, -1, -1, -2, -3, -1, -2, -2, -2, -1, 2, 7, 0], // Y
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], // *
    ],
    gapOP: -7,
    gapEP: -1,
    gapOPMSA: -15,
    gapEPMSA: -1,
};

const BLASTZ = {
    matrix: [
        //	A		C		G		T
        [91, -114, -31, -123],
        [-114, 100, -125, -31],
        [-31, -125, 100, -114],
        [-123, -31, -114, 91],
    ],
    gapOP: -400,
    gapEP: -30,
    gapOPMSA: -400,
    gapEPMSA: -1,
};

export function setAlignmentParameters(isMSA = false) {
    var m;

    switch (typeSeq) {
        case SEQUENCE_TYPE.PROTEIN:
            m = BLOSUM62;
            scoringMatrix = m.matrix;
            gapOP = isMSA ? m.gapOPMSA : m.gapOP;
            gapEP = isMSA ? m.gapEPMSA : m.gapEP;
            break;
        default:
            m = BLASTZ;
            scoringMatrix = m.matrix;
            gapOP = isMSA ? m.gapOPMSA : m.gapOP;
            gapEP = isMSA ? m.gapEPMSA : m.gapEP;
    }
}

const dotplotWindow = 15;

export function setSeqType(type: SEQUENCE_TYPE) {
    typeSeq = type;
}

export function getAlignmentParameters() {
    return {
        scoringMatrix,
        gapOP,
        gapEP,
    };
}
