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

type TAlignMethod = 'diag'|'complete'|'auto';
type TAlignSequenceType = 'amino'|'nucleic'|'auto';
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
    method   : TAlignMethod,

    /**
     * Sequence type
     * `auto` performs a detection based on sequences content
     * @default 'auto'
     */
    type     : TAlignSequenceType,

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

type TBioMSAConfig = {
    sequenceType: TAlignSequenceType,
    gapchar: string,
    alignmentMethod: TAlignMethod,
    debug: boolean
}

const DEFAULT_CONFIG: TBioMSAConfig = {
    sequenceType: 'auto',
    gapchar: '-',
    alignmentMethod: 'auto',
    debug: false,
}

class BioMSAClass {
    private sequences: TSequence[] = [];
    private typeSeq = SEQUENCE_TYPE.UNSET;
    private config: TBioMSAConfig = {
        ...DEFAULT_CONFIG
    }

    /**
     * Add sequence to the internal list of sequences of the aligner.
     * All sequences must be of same type.
     * @param {(string|string[])} seq sequence or array of sequences to add
     * @returns void
     * @throws Will throw an exception if sequence added is not of same type as
     *      the other internal sequences.
     */
    private addSequence(seq: string | string[]) {

            // Handle Array parameter overload (recursive)

        if (Array.isArray(seq)) {
            seq.forEach((s) => {
                this.addSequence(s);
            });
            return;
        }

            // Sanity checks

        if (typeof seq !== 'string') {
            throw new TypeError('String type expected for sequences to add.');
        }
        seq = seq.toUpperCase();

            // In sequence type automatic detection mode, some checks are added
            // to infer sequence type and sequence type consistency
            // In user defined mode, we let the user do their own checks and we
            // create the sequences eagerly (unrecognized residues are interpolated
            // randomly)

        if (this.config.sequenceType === 'auto') {
            const type = getSequenceType(seq);

                // First one

            if (this.typeSeq === SEQUENCE_TYPE.UNSET) {
                this.typeSeq = type;

                // Once set, all the others must follow

            } else if (this.typeSeq !== type) {
                throw new Error('All sequences must be of same type.');
            }
        }

            // Add sequence to store

        this.sequences.push(makeSequence(seq));
    }

    /**
     * Reset the internal sequence store content
     */
    private reset() {
        this.sequences = [];
        this.typeSeq = SEQUENCE_TYPE.UNSET;
        this.setDefaultConfiguration();
    }

    private setDefaultConfiguration () {
        this.config = {
            ...DEFAULT_CONFIG
        }
    }

    private setUserConfiguration (opt: Partial<TAlignOpt>) {
        function isValid<K>(value: any, allowed: K[]) {
            return (allowed.some(v => v == value));
        }

        if (
            opt.method
            && isValid<TAlignOpt['method']>(opt.method, ['complete', 'diag'])
        ) this.config.alignmentMethod = opt.method;

        if (
            opt.type
            && isValid<TAlignOpt['type']>(opt.type, ['amino', 'nucleic'])
        ) {
            this.config.sequenceType = opt.type;
            this.typeSeq = (opt.type === 'amino') ? SEQUENCE_TYPE.PROTEIN : SEQUENCE_TYPE.NUCLEIC;
        }

        if (opt.gapchar?.length === 1) this.config.gapchar = opt.gapchar;

        if (opt.debug !== undefined) this.config.debug = !!opt.debug;
    }

    /**
     * Main alignment function. Takes an array of sequences (strings) as an
     * argument or aligns sequences that were previously added to internal
     * store.
     * This function returns a promise.
     * @param seqArr
     * @returns { Promise<string[]> }
     */
    public align(seqArr: string[], opt?: Partial<TAlignOpt>) {
        if (DEBUG) Log.start();

        const msa = new Promise<any[]>((resolve, reject) => {

                // Sequences prerequisites

            if (!Array.isArray(seqArr)
                || seqArr.some(v => typeof(v) !== 'string')
            ) {
                return reject('Array of sequences expected');
            }

            if (seqArr.length < 2) {
                return reject('At least 2 sequences are required.');
            }

            this.reset();
            if (opt) this.setUserConfiguration(opt);

            this.addSequence(seqArr);

            if (DEBUG) Log.add('Prepared sequences');

            const lParam = getAlignmentParameters(this.typeSeq);

            if (DEBUG) Log.add('Get sequences type');

            let doNoAlign = false;
            if (this.config.alignmentMethod === 'auto')
                doNoAlign = this.sequences.some(seq => seq.rawSeq.length > 1600);
            else doNoAlign = this.config.alignmentMethod === 'diag';

                // Start alignment

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

                const lAlignment = [
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
