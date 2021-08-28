/**
 * @file utils
 * @author Paul Pillot <paul.pillot@libmol.org>
 * @license MIT
 * @copyright 2020
 */
import BitArray from '../utils/bitarray';
import { DEBUG, SEQUENCE_TYPE } from '../align/params';
import Log from '../utils/logger';

export type TSequence = {
    /** Sequence as a string */
    rawSeq: string;
    /** Sequence encoded as a number array */
    encodedSeq: Uint8Array;
    /** Sequence encoded as a number array in a compressed alphabet (Dayhoff) */
    compressedSeq: Uint8Array;
    /** Sequence type enum (PROTEIN or NUCLEIC) */
    type: SEQUENCE_TYPE;
};

const regAmino    = /^[ACGTRNDEQHILKMFPSWYV]+$/;
const regExtNuc   = /^[ABCDGHKMNRSTUVWY]+$/i;
const regExtAmino = /^[ABCGTRNDEQHIJLKMFPSWYUVXZ]+$/i;
const regNotRNA   = /[^ACGU]/i;
const regNotDNA   = /[^ACGT]/i;
const CODE_NOT_FOUND = 255;

/**
 * Guess the type of sequence from its raw content
 * @param {string} seq sequence to guess the type from
 */
export function getSequenceType(seq: string) {

    if (!regNotRNA.test(seq) || !regNotDNA.test(seq)) {
        return SEQUENCE_TYPE.NUCLEIC
    };

    const isExtendedNuc = regExtNuc.test(seq);
    const isAmino = regAmino.test(seq);

    if (isAmino && !isExtendedNuc) return SEQUENCE_TYPE.PROTEIN;

    if (isAmino && isExtendedNuc) return guessSequenceType(seq);

    if (regExtAmino.test(seq)) return SEQUENCE_TYPE.PROTEIN;

    throw new Error('Unrecognized sequence type: ' + seq);
}

function guessSequenceType (seq: string) {
    const SAMPLE_SIZE = 100;
    let lNucScore = 0;
    const LEN = Math.min(seq.length, SAMPLE_SIZE);

    for (let i = 0; i < LEN; i++) {
        switch (seq[i]) {
            case 'A':
            case 'T':
            case 'U':
            case 'G':
            case 'C':
            case 'N':
                lNucScore ++;
        }
    }

    if (lNucScore / LEN > Math.SQRT1_2) return SEQUENCE_TYPE.NUCLEIC;

    return SEQUENCE_TYPE.PROTEIN;
}

/**
 * Index aa codes based on ASCII code
 */
export const aaToNum = new Uint8Array(96);
aaToNum.fill(CODE_NOT_FOUND);
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
nucToNum.fill(CODE_NOT_FOUND);
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
    const encodedTab = new Uint8Array(seq.length);
    const convTable = type === SEQUENCE_TYPE.PROTEIN ? aaToNum : nucToNum;
    let lEncoding = 0;

    for (let i = 0, imax = seq.length; i < imax; i++) {
        lEncoding = convTable[seq.charCodeAt(i)];
        if (lEncoding === CODE_NOT_FOUND) lEncoding = setRandomCode(seq[i], type);

        encodedTab[i] =  lEncoding;
    }

    return encodedTab;
}

/**
 * Interpolates ambiguous letters in sequences by picking randomly a standard
 * residue amongst the possible choices.
 * @param letter
 * @param type
 * @returns
 */
function setRandomCode(letter: string, type: SEQUENCE_TYPE) {
    const lRand = Math.floor(Math.random() * 100);
    switch (type) {
        case SEQUENCE_TYPE.NUCLEIC:
            switch (letter) {
                case 'M': return [0, 1][lRand % 2]; // A or C (amino)
                case 'R': return [0, 2][lRand % 2]; // A or G (purine)
                case 'W': return [0, 3][lRand % 2]; // A or T (weak)
                case 'S': return [1, 2][lRand % 2]; // C or G (strong)
                case 'Y': return [1, 3][lRand % 2]; // C or T (pyrimidine)
                case 'K': return [2, 3][lRand % 2]; // G or T (keto)
                case 'V': return [0, 1, 2][lRand % 3]; // A or C or G
                case 'H': return [0, 1, 3][lRand % 3]; // A or C or T
                case 'D': return [0, 2, 3][lRand % 3]; // A or G or T
                case 'B': return [1, 2, 3][lRand % 3]; // C or G or T
                case 'N': // passthrough intentional
                default: return lRand % 4;
            }

        case SEQUENCE_TYPE.PROTEIN:
        default:
            switch (letter) {
                case 'B': return [2, 11][lRand % 2]; // Asx
                case 'Z': return [3, 13][lRand % 2]; // Glx
                case 'J': return [7, 9 ][lRand % 2]; // Xle (Leu or Ile)
                case 'U': // passtrhough intentional    Selenocysteine
                case 'X': // passthrough intentional
                default: return lRand % 20;
            }
    }
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
export function compressToDayhoff(encodedSeq: Uint8Array) {
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
    const bitsetLength = Math.pow(alphabetSize, 6);
    const lKmer = tabSeq.map((seq) => {
        const bitset = new BitArray(bitsetLength);
        const seqAsNum = isProtein ? seq.compressedSeq : seq.encodedSeq;

        // Value for a tuple is 1st letter + 6*2nd letter + 36*3rd letter...

        let tupleval = 0;
        for (let i = 0; i <= 5; i++ ) {
            tupleval += seqAsNum[i] * Math.pow(alphabetSize, i);
        }
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

    if (DEBUG) Log.add('K-mer bitset computation');

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

    if (DEBUG) Log.add('K-mer distance computation');
    return distTab;
}

export function sortMSA<K>(msa: K[], order: number[]): K[] {
    let lSorted = order.slice();
    order.forEach((v, k) => lSorted[v] = k)
    return lSorted.map(v => msa[v]);
}

const dayhoffPams = [
    195, /* 75.0% observed d; 195 PAMs estimated = 195% estimated d */
    196, /* 75.1% observed d; 196 PAMs estimated */
    197, 198, 199, 200, 200, 201, 202, 203,
    204, 205, 206, 207, 208, 209, 209, 210, 211, 212,
    213, 214, 215, 216, 217, 218, 219, 220, 221, 222,
    223, 224, 226, 227, 228, 229, 230, 231, 232, 233,
    234, 236, 237, 238, 239, 240, 241, 243, 244, 245,
    246, 248, 249, 250, /* 250 PAMs = 80.3% observed d */
    252, 253, 254, 255, 257, 258,
    260, 261, 262, 264, 265, 267, 268, 270, 271, 273,
    274, 276, 277, 279, 281, 282, 284, 285, 287, 289,
    291, 292, 294, 296, 298, 299, 301, 303, 305, 307,
    309, 311, 313, 315, 317, 319, 321, 323, 325, 328,
    330, 332, 335, 337, 339, 342, 344, 347, 349, 352,
    354, 357, 360, 362, 365, 368, 371, 374, 377, 380,
    383, 386, 389, 393, 396, 399, 403, 407, 410, 414,
    418, 422, 426, 430, 434, 438, 442, 447, 451, 456,
    461, 466, 471, 476, 482, 487, 493, 498, 504, 511,
    517, 524, 531, 538, 545, 553, 560, 569, 577, 586,
    595, 605, 615, 626, 637, 649, 661, 675, 688, 703,
    719, 736, 754, 775, 796, 819, 845, 874, 907, 945, /* 92.9% observed; 945 PAMs */
    988 /* 93.0% observed; 988 PAMs */
];

/**
 * Computes the ratio between the number of matches in both sequences and the shortest sequence length
 * @param   {string} [seqA='',] sequence A, may contain indels
 * @param   {string} seqB       sequence B, may contain indels
 * @returns {number} fractional identity between both sequences
 */
function fractionalIdentity(seqA: string, seqB: string) {
    seqA = seqA || '';
    seqB = seqB || '';

    var match = 0,
        tokenA = '',
        tokenB = '',
        lA = 0,
        lB = 0,
        maxI = seqA.length;

    for (var i = 0; i < maxI; i++) {
        tokenA = seqA.charAt(i);
        tokenB = seqB.charAt(i);

        if ((tokenA === '-') || (tokenB === '-')) {
            lA += (tokenA !== '-') ? 1 : 0;
            lB += (tokenB !== '-') ? 1 : 0;
            continue;
        } else {
            lA++;
            lB++;
            if (tokenA === tokenB) {
                match++;
            }
        }
    }
    return match / Math.min(lA, lB);
}


/**
 * Compute distance matrix between aligned sequences from a MSA.
 * Corrects distances using Kimura method.
 */
export function distanceKimura (msa: string[]) {

    const l = msa.length; // nb of sequences
    const distMatrix = new Array(l); // matrix l x l
    let dPctId = 0; // Percentage of Identity between sequences
    let p = 0;      // Percentage of difference between sequences
    let pIndex = 0; // When p > 0.75 Kimura correction is replaced by a lookup table at this index
    let dk = 0;     // Computed distance

    //parcours de la matrice
    for (var i = 0; i < l; i++) {

        if ( distMatrix[i] === undefined ) {
            distMatrix[i] = new Array(l);
        }

        distMatrix[i][i] = 0;

        for (var j = i + 1; j < l; j++) {

            if (distMatrix[j] === undefined) {
                distMatrix[j] = new Array(l);
            } else if (distMatrix[i][j] !== undefined) {
                continue;
            }

            //  "As sequences diverge, there is an increasing probability of
            // multiple mutations at a single site. To correct for this, we use
            // the following distance estimate" Edgar. BMC Bioinformatics 2004

            dPctId = fractionalIdentity(msa[i], msa[j]);
            p = 0, 1 - dPctId;

            // Kimura formula used for differences < 75%
            if (p < 0.75) {

                dk = - Math.log(1 - p - (p * p) / 5);
                dk = Math.max(0, dk);   // clamp to 0 as the above can return -0!!!

            } else if (p > 0.93) {
                dk = 10;
            } else {

                pIndex = Math.floor((p - 0.75) * 1000 + 0.5);

                if ((pIndex < 0) || (pIndex >= dayhoffPams.length)) {
                    console.trace('Dayhoff parameter not found');
                } else {
                    dk = dayhoffPams[pIndex] / 100;
                }

            }
            distMatrix[i][j] = dk;
            distMatrix[j][i] = dk;
        }
    }

    return distMatrix;
};

export function makeSequence(seq: string, type?: SEQUENCE_TYPE): TSequence {
    const lType = type ?? getSequenceType(seq);
    const encodedSeq = encodeSeqToNum(seq, lType);
    return {
        rawSeq: seq,
        encodedSeq,
        compressedSeq: lType === SEQUENCE_TYPE.PROTEIN ? compressToDayhoff(encodedSeq) : encodedSeq,
        type: lType,
    };
}