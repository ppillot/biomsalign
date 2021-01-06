/**
 * @file params
 * @author Paul Pillot <paul.pillot@libmol.org>
 * @license MIT
 * @copyright 2020
 */

import { SEQUENCE_TYPE } from './sequence';

export const DEBUG = true;

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
    gapOP: -11,
    gapEP: -1,
};

const VTML240 = {
    matrix: [
        //  A    C    D    E    F    G    H    I    K    L    M    N    P    Q    R    S    T    V    W    Y  X
        [58, 23, -12, -7, -44, 10, -23, -14, -14, -27, -17, -8, 1, -9, -22, 23, 15, 5, -74, -45, 0], // A
        [23, 224, -67, -63, -50, -30, -29, 1, -56, -41, -6, -33, -44, -53, -43, 15, 2, 18, -93, -6, 0], // C
        [-12, -67, 111, 59, -104, -4, 4, -84, 6, -88, -65, 48, -13, 18, -29, 5, -7, -63, -105, -73, 0], // D
        [-7, -63, 59, 85, -83, -17, -1, -63, 25, -60, -47, 15, -12, 40, -8, 1, -7, -47, -108, -51, 0], // E
        [-44, -50, -104, -83, 144, -93, 4, 12, -74, 36, 30, -64, -67, -56, -65, -43, -41, -3, 63, 104, 0], // F
        [10, -30, -4, -17, -93, 140, -32, -95, -27, -91, -75, 4, -36, -29, -32, 5, -26, -68, -80, -79, 0], // G
        [-23, -29, 4, -1, 4, -32, 137, -50, 6, -37, -42, 21, -23, 27, 19, -4, -12, -44, -13, 48, 0], // H
        [-14, 1, -84, -63, 12, -95, -50, 86, -53, 53, 47, -62, -60, -47, -55, -43, -8, 69, -27, -24, 0], // I
        [-14, -56, 6, 25, -74, -27, 6, -53, 75, -48, -30, 13, -12, 34, 68, -3, -4, -44, -71, -49, 0], // K
        [-27, -41, -88, -60, 36, -91, -37, 53, -48, 88, 62, -63, -48, -36, -48, -47, -25, 36, -11, -4, 0], // L
        [-17, -6, -65, -47, 30, -75, -42, 47, -30, 62, 103, -45, -54, -21, -31, -35, -9, 31, -46, -20, 0], // M
        [-8, -33, 48, 15, -64, 4, 21, -62, 13, -63, -45, 89, -25, 12, 2, 22, 10, -51, -79, -29, 0], // N
        [1, -44, -13, -12, -67, -36, -23, -60, -12, -48, -54, -25, 160, -6, -20, 5, -12, -42, -76, -83, 0], // P
        [-9, -53, 18, 40, -56, -29, 27, -47, 34, -36, -21, 12, -6, 75, 34, 1, -4, -37, -92, -48, 0], // Q
        [-22, -43, -29, -8, -65, -32, 19, -55, 68, -48, -31, 2, -20, 34, 113, -10, -14, -49, -58, -39, 0], // R
        [23, 15, 5, 1, -43, 5, -4, -43, -3, -47, -35, 22, 5, 1, -10, 53, 32, -28, -62, -31, 0], // S
        [15, 2, -7, -7, -41, -26, -12, -8, -4, -25, -9, 10, -12, -4, -14, 32, 68, 0, -87, -40, 0], // T
        [5, 18, -63, -47, -3, -68, -44, 69, -44, 36, 31, -51, -42, -37, -49, -28, 0, 74, -61, -32, 0], // V
        [-74, -93, -105, -108, 63, -80, -13, -27, -71, -11, -46, -79, -76, -92, -58, -62, -87, -61, 289, 81, 0], // W
        [-45, -6, -73, -51, 104, -79, 48, -24, -49, -4, -20, -29, -83, -48, -39, -31, -40, -32, 81, 162, 0], // Y
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], // X
    ],
    center: 22,
    gapOP: -300,
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
    gapEP: -60,
    center: 120,
};

/**
 * Recursively adds a constant center to all the matrix values
 * @param {Array} matrix Scoring matrix
 * @param {number} center half gap penalty (see Edgar's MUSCLE for details on
 * 			this optimization). In a few words, adding a constant half penalty
 * 			to all substitution scores -i.e. matches-, makes it so that we
 * 			can skip to substract penalties when extending a gap as it will
 * 			be comparatively lower scored).
 */
function addCenter(matrix: number[][], center: number) {
    return matrix.map((row) => row.map((val) => val + center));
}

type AlignParam = {
    gapOpen: number,
    gapExtend: number
}

export function setAlignmentParameters(p?: Partial<AlignParam>) {
    var m;

    switch (typeSeq) {
        case SEQUENCE_TYPE.PROTEIN:
            m = VTML240;
            break;
        default:
            m = BLASTZ;
            break;
    }

    scoringMatrix = addCenter(
        m.matrix,
        p?.gapExtend ? -p.gapExtend * 2 : m.center
    );
    gapOP = m.gapOP;
}

export function setSeqType(type: SEQUENCE_TYPE) {
    typeSeq = type;
}

export function getAlignmentParameters() {
    return {
        type: typeSeq,
        scoringMatrix,
        gapOP,
    };
}