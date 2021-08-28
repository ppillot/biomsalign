/**
 * @file index
 * @author Paul Pillot <paul.pillot@libmol.org>
 * @license MIT
 * @copyright 2020-2021
 */

import {
    makeSequence,
    getSequenceType,
    TSequence,
} from './sequence/sequence';
import { DEBUG, SEQUENCE_TYPE, getAlignmentParameters } from './align/params';
import { pairwiseAlignment } from './align/align';
import { progressiveAlignment } from "./align/progressive.alignment";
import Log from './utils/logger';
import { centerStarNoAlign, noalignPair } from './align/noalign';
import { estringTransform } from './utils/estring';

type TAlignOpt = {
    /**
     * Gap open penalty (generally a negative number)
     * @default -300|-400 set based on the sequence type
     */
    gapopen  : number,

    /**
     * Gap extend penalty (generally a negative number)
     * @default -11|-60 set based on the sequence type
     */
    gapextend: number,

    /**
     * Set true to report debug information to the JS console
     * @default false
     */
    debug    : boolean,

    /**
     * Alignment method. `auto` is the default value. For sequences length up
     * to 1600 residues, a `complete` Needleman-Wunsch global alignment is
     * performed. Beyond, a diagonal finding algorithm  (`diag`) is run first
     * and global alignment is restricted to regions between diagonals.
     * `complete` and `diag` parameters force the use of one of those
     * strategies.
     * @default 'auto'
     */
    method   : 'diag'|'complete'|'auto',

    /**
     * Sequence type
     * `auto` performs a detection based on sequences content
     * @default 'auto'
     */
    type     : 'amino'|'nucleic'|'auto',

    /**
     * Substitution scoring matrix. Residues are sorted lexically by their 1
     * letter code (e.g. "A, C, D, E..." for proteins, not "Ala, Arg, Asn,...")
     * which is a discrepancy with BLOSUM matrices found on NCBI.
     * @default set based on sequence type
     */
    matrix   : number[][],

    /**
     * Char that is used as the symbol for gaps in sequences
     * @default '-'
     */
    gapchar  : string
}


class BioMSAClass {
    private sequences: TSequence[] = [];
    private typeSeq: SEQUENCE_TYPE = SEQUENCE_TYPE.NUCLEIC;

    /**
     * Add sequence to the internal list of sequences of the aligner.
     * All sequences must be of same type.
     * @param {(string|string[])} seq sequence or array of sequences to add
     * @returns void
     * @throws Will throw an exception if sequence added is not of same type as
     *      the other internal sequences.
     */
    public addSequence(seq: string | string[]) {

            // Handle Array parameter overload (recursive)

        if (Array.isArray(seq)) {
            seq.forEach((s) => {
                this.addSequence(s);
            });
            return;
        }

            // Sanity checks

        if (typeof seq !== 'string') {
            throw new Error('String type expected for sequences to add.');
        }

        const type = getSequenceType(seq);
        if (this.sequences.length === 0) {
            this.typeSeq = type;
        } else if (this.typeSeq !== type) {
            throw new Error('All sequences must be of same type.');
        }

            // Add sequence to store

        this.sequences.push(makeSequence(seq));
    }

    /**
     * Reset the internal sequence store content
     */
    public reset() {
        this.sequences = [];
        this.typeSeq = SEQUENCE_TYPE.NUCLEIC;
    }


    /**
     * Main alignment function. Takes an array of sequences (strings) as an
     * argument or aligns sequences that were previously added to internal
     * store.
     * This function returns a promise.
     * @param seqArr
     * @returns { Promise<string[]> }
     */
    public align(seqArr?: string[]) {
        if (DEBUG) Log.start();

        const msa = new Promise<any[]>((resolve, reject) => {

                // Sequences prerequisites

            if (seqArr !== undefined) {
                if (Array.isArray(seqArr)) {
                    this.reset();
                    this.addSequence(seqArr);
                } else {
                    return reject('Array of sequences expected');
                }
            } else if (this.sequences.length < 2) {
                return reject('At least 2 sequences are required.');
            }

            if (DEBUG) Log.add('Prepared sequences');

            const lParam = getAlignmentParameters(this.sequences[0].type);

            if (DEBUG) Log.add('Get sequences type');

            const doNoAlign = this.sequences.some(seq => seq.rawSeq.length > 1600);

            if (this.sequences.length == 2) {
                let lEStrings: number[][];
                if (doNoAlign) {
                    lEStrings = noalignPair(
                        this.sequences[0],
                        this.sequences[1],
                        lParam
                    );
                } else {
                    let lResult = pairwiseAlignment(
                        this.sequences[0],
                        this.sequences[1],
                        lParam
                    );
                    lEStrings = lResult.estrings;
                }

                if (DEBUG) Log.summary();

                let lAlignment = [
                    estringTransform(this.sequences[0].rawSeq, lEStrings[0]),
                    estringTransform(this.sequences[1].rawSeq, lEStrings[1])
                ];

                return resolve(lAlignment);
            } else {

                let lResult = doNoAlign ?
                    centerStarNoAlign(this.sequences, lParam)
                    : progressiveAlignment(this.sequences, lParam);

                if (DEBUG) Log.summary();

                if (!lResult) return reject('An error occured');

                return resolve(lResult);
            }
        });

        return msa;
    }

}
const BioMSA = new BioMSAClass();
export default BioMSA;
