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

class BioMSAClass {
    private sequences: TSequence[] = [];
    private typeSeq: SEQUENCE_TYPE = SEQUENCE_TYPE.NUCLEIC;

    public addSequence(seq: string | string[]) {
        if (Array.isArray(seq)) {
            seq.forEach((s) => {
                this.addSequence(s);
            });
            return;
        }

        if (typeof seq !== 'string') {
            throw new Error('String type expected for sequences to add.');
        }

        const type = getSequenceType(seq);
        if (this.sequences.length === 0) {
            this.typeSeq = type;
        } else if (this.typeSeq !== type) {
            throw new Error('All sequences must be of same type.');
        }

        this.sequences.push(makeSequence(seq));
    }

    public reset() {
        this.sequences = [];
        this.typeSeq = SEQUENCE_TYPE.NUCLEIC;
    }

    public align(seqArr?: string[]) {
        if (DEBUG) Log.start();

        const msa = new Promise<any[]>((resolve, reject) => {
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
