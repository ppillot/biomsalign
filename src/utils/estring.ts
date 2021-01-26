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