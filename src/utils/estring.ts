/**
 * @file estring.ts
 * @author Paul Pillot <paul.pillot@libmol.org>
 * @license MIT
 * @copyright 2021
 * @description E-strings (Edit strings) are operators that act on a string to
 * add insertion characters.
 * The follofing example is taken from MUSCLE (a) presentation:
 * Sequence: MQTIF
 * Estring: <3, -1, 2>
 * Transformation: MQT-IF
 * Positive integers are the counts of letters conserved in a fragment, negative
 * numbers are the counts of insertions.
 * Note that the sum of the positive numbers is the size of the initial string,
 * and the sum of the absolute values are the final size of the string.
 * The alignment of two sequences can be described by the operators that
 * transform the initial strings in 2 strings of same size.
 *
 * (a) Edgar, R.C. MUSCLE: a multiple sequence alignment method with reduced
 * time and space complexity. BMC Bioinformatics 5, 113 (2004).
 * https://doi.org/10.1186/1471-2105-5-113
 */

/**
 * Edits a string such as a sequence, by inserting gap symbols based on the
 * content of and Edit string.
 *
 * @export
 * @param {string} text
 * @param {number[]} estring
 * @returns {string}    edited string
 */
export function estringTransform (text: string, estring: number[]) {
    const lPieces: string[] = [];
    let lPos = 0;   // cursor in text
    for (let i = 0; i < estring.length; i++) {
        let n = estring[i];
        if (n < 0) {
            lPieces.push('-'.repeat(-n));
        } else {
            lPieces.push(text.substr(lPos, n));
            lPos += n;
        }
    }
    return lPieces.join('');
}

/**
 * returns an edit string equivalent to the composition of the two edit strings
 * passed as an argument
 *
 * @export
 * @param {number[]} estringB
 * @param {number[]} estringA
 * @returns {number[]}
 */
export function estringProduct (estringB: number[], estringA: number[]): number[] {
    const lSteps: number[] = [0];
    let lPosA = 0;
    const lEstringA = estringA.slice(); // make copy, will be modified

    // add new value in lSteps. Collapses tail value if they have same sign
    function pushOperator (val: number) {
        if ((val ^ lSteps[lSteps.length - 1]) >= 0) {    // same sign
            lSteps[lSteps.length - 1] += val;
        } else {    // opposite sign: new item in array
            lSteps.push(val);
        }
    }

    for (let i = 0; i < estringB.length; i++) {
        let n = estringB[i];
        if (n < 0) {
            pushOperator(n);
        } else {
            while (n > 0 && lPosA < lEstringA.length) {
                let lValA = lEstringA[lPosA];
                if (lValA < 0) {
                    if (n < -lValA) {
                        lEstringA[lPosA] += n;
                        pushOperator(-n);
                        n = 0;
                    } else {
                        pushOperator(lValA);
                        n += lValA;
                        lPosA ++;
                    }
                } else {
                    if (n < lValA) {
                        lEstringA[lPosA] -= n;
                        pushOperator(n);
                        n = 0;
                    } else {
                        pushOperator(lValA);
                        n -= lValA;
                        lPosA ++;
                    }
                }
            }
        }
    }

    return lSteps;
}

/**
 * Edit strings catenation.
 *
 * @export
 * @param {number[]} estringA
 * @param {number[]} estringB
 * @returns {number[]}
 */
export function estringCat (estringA: number[], estringB: number[]) {
    if (Math.sign(estringA[estringA.length - 1]) === Math.sign(estringB[0])) {
        let lRet = [...estringA];
        lRet[lRet.length - 1] += estringB[0];
        lRet.push(...estringB.slice(1));    // avoid changing reference
        return lRet;
    }
    return [...estringA, ...estringB];
}

export const enum EPATH_2_STRING {
    NO_REVERSE = 1
};

export function epath2estring (epath: number[], opt = 0): number[] {
    const lEstring: number[] = [];

    let lAcc = 0;
    let lVal = 0;
    for (let i = 0; i < epath.length; i++) {
        lVal = epath[i];
        if ((lVal ^ lAcc) >= 0) lAcc += lVal;
        else {
            lEstring.push(lAcc);
            lAcc = lVal;
        }
    }
    lEstring.push(lAcc);
    if (lEstring[0] === 0) lEstring.shift();

        // Usually, the path comes from a traceback matrix and must be reversed

    if (!(opt & EPATH_2_STRING.NO_REVERSE)) lEstring.reverse();

    return lEstring;
}

export function estringLength (estring: number[]): number {
    let lSize = 0;
    let lVal = 0;
    for (let i = 0; i < estring.length; i++) {
        lVal = estring[i];
        lSize += lVal < 0 ? - lVal : lVal;
    }
    return lSize;
}