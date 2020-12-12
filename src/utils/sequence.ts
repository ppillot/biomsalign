/**
 * @file utils
 * @author Paul Pillot <paul.pillot@libmol.org>
 * @license MIT
 * @copyright 2020
 */
import BitArray from './bitarray';

export enum SEQUENCE_TYPE {
    PROTEIN,
    NUCLEIC,
}

export type TSequence = {
    /** Sequence as a string */
    rawSeq: string;
    /** Sequence encoded as a number array */
    encodedSeq: number[];
    /** Sequence encoded as a number array in a compressed alphabet (Dayhoff) */
    compressedSeq: number[];
    /** Sequence type enum (PROTEIN or NUCLEIC) */
    type: SEQUENCE_TYPE;
};

const regAmAc = /[rndbeqyhilkmfpswyv]/gi;
const regNotRNA = /[^acgu]/gi;
const regNotDNA = /[^acgt]/gi;

/**
 * Guess the type of sequence from its raw content
 * @param {string} seq sequence to guess the type from
 */
export function getSequenceType(seq: string) {
    if (regAmAc.test(seq)) return SEQUENCE_TYPE.PROTEIN;

    if (!regNotRNA.test(seq) || !regNotDNA.test(seq)) return SEQUENCE_TYPE.NUCLEIC;

    throw new Error('Unrecognized sequence type: ' + seq);
}

/**
 * Index aa codes based on ASCII code
 */
export const aaToNum = new Uint8Array(96);
aaToNum[65] = 0; // A
aaToNum[67] = 1; // C
aaToNum[68] = 2; // D
aaToNum[69] = 3; // E
aaToNum[70] = 4; // F
aaToNum[71] = 5; // G
aaToNum[72] = 6; // H
aaToNum[73] = 7; // I
aaToNum[75] = 8; // K
aaToNum[76] = 9; // L
aaToNum[77] = 10; // M
aaToNum[78] = 11; // N
aaToNum[80] = 12; // P
aaToNum[81] = 13; // Q
aaToNum[82] = 14; // R
aaToNum[83] = 15; // S
aaToNum[84] = 16; // T
aaToNum[86] = 17; // V
aaToNum[87] = 18; // W
aaToNum[89] = 19; // Y
aaToNum[95] = 20; // _

/**
 * Index nucleotides codes based on ASCII
 */
export const nucToNum = new Uint8Array(96);
nucToNum[65] = 0; // A
nucToNum[67] = 1; // C
nucToNum[71] = 2; // G
nucToNum[84] = 3; // T
nucToNum[85] = 3; // U, uracile aligns with thymine

/**
 * Encode sequence string as a number array. This is paramount to extract
 * substitution scores in O(1) time instead of O(ln) in the innermost loop.
 *
 * @export
 * @param {string} seq
 * @param {SEQUENCE_TYPE} type
 * @returns {number[]}
 */
export function encodeSeqToNum(seq: string, type: SEQUENCE_TYPE) {
    const encodedTab: number[] = [];
    const convTable = type === SEQUENCE_TYPE.PROTEIN ? aaToNum : nucToNum;

    for (let i = 0, imax = seq.length; i < imax; i++) {
        encodedTab.push(convTable[seq.charCodeAt(i)]);
    }

    return encodedTab;
}

/**
 * Converts a proteic sequence to Dayhoff compressed alphabet. This increases
 * the probability of matches between divergent sequences by reducing the
 * number of combinations (20 amino acids to 6 symbols)
 *
 * @export
 * @param {number[]} encodedSeq Sequence encoded as a number array
 * @returns {number[]} compressed sequence encoded as a number array
 */
export function compressToDayhoff(encodedSeq: number[]) {
    const compressedSeq = encodedSeq.map((val) => toDayhoff[val]);
    return compressedSeq;
}

const toDayhoff = [
    0, // A -> A
    1, // C -> B
    2, // D -> C
    2, // E -> C
    3, // F -> D
    0, // G -> A
    4, // H -> E
    5, // I -> F
    4, // K -> E
    5, // L -> F
    5, // M -> F
    2, // N -> C
    0, // P -> A
    2, // Q -> C
    4, // R -> E
    0, // S -> A
    0, // T -> A
    5, // V -> F
    3, // W -> D
    3, // Y -> D
];

export function distanceMatrix(tabSeq: TSequence[]) {
    // A 6-tuple in the Dayhoff alphabet can result in 6^6 combinations.
    // A 6-tuple in the nucleic alphabet can result in 4^6 combinations.
    // To compute a distance between sequences, we evaluate how many 6-tuples
    // they have in common. To speed-up this process, we consider only a binary
    // comparison, i.e. does the tuple exists in the sequence or not (we ignore
    // how many time it actually occurs, which is a limitation of this method
    // for sequences with low complexity).

    // Make the 6/4-kmer bitsets

    const isProtein = tabSeq[0].type === SEQUENCE_TYPE.PROTEIN;
    const alphabetSize = isProtein ? 6 : 4;
    const bitsetLength = Math.pow(6, 6);
    const lKmer = tabSeq.map((seq) => {
        const bitset = new BitArray(bitsetLength);
        const seqAsNum = isProtein ? seq.compressedSeq : seq.encodedSeq;

        // Value for a tuple is 1st letter + 6*2nd letter + 36*3rd letter...

        let tupleval = seqAsNum.slice(0, 5).reduce((acc, val, i) => {
            return acc + val * Math.pow(alphabetSize, i);
        }, 0);
        bitset.set(tupleval);

        // Compute next tuple from previous one

        const topLetterFactor = Math.pow(alphabetSize, 5);
        for (let i = 6, imax = seqAsNum.length; i < imax; i++) {
            // remove lowest letter from tuple

            tupleval -= seqAsNum[i - 6];

            // shift value by 6

            tupleval /= 6;

            // add highest letter to tuple

            tupleval += seqAsNum[i] * topLetterFactor;

            bitset.set(tupleval);
        }

        return bitset;
    });

    // compute distance matrix values

    const l = tabSeq.length;
    const distTab: number[][] = tabSeq.map(() => []);
    let lKmerI: BitArray;
    let lSeqILen: number;
    let lDistance: number;

    for (let i = 0; i < l; i++) {
        distTab[i][i] = 0;
        lKmerI = lKmer[i];
        lSeqILen = tabSeq[i].compressedSeq.length;

        for (let j = i + 1; j < l; j++) {
            lDistance = 1 - lKmerI.getIntersectionSize(lKmer[j]) / lSeqILen;
            distTab[j][i] = distTab[i][j] = lDistance;
        }

    }
    return distTab;
}
