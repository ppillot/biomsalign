/**
 * @file index
 * @author Paul Pillot <paul.pillot@libmol.org>
 * @license MIT
 * @copyright 2020
 */

import {
    compressToDayhoff,
    encodeSeqToNum,
    getSequenceType,
    TSequence,
} from './utils/sequence';
import { setSeqType, DEBUG, SEQUENCE_TYPE, setAlignmentParameters } from './utils/params';
import { pairwiseAlignment, progressiveAlignment } from './utils/align';
import Log from './utils/logger';
import { noalignPair } from './utils/noalign';

export class BioMSA {
    private sequences: TSequence[];
    private typeSeq: SEQUENCE_TYPE = SEQUENCE_TYPE.PROTEIN;

    constructor() {
        this.sequences = [];
    }

    public addSequences(seq: string | string[]) {
        if (Array.isArray(seq)) {
            seq.forEach((s) => {
                this.addSequences(s);
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

        this.sequences.push(this.makeSequenceObj(seq));
    }

    private makeSequenceObj(seq: string): TSequence {
        const encodedSeq = encodeSeqToNum(seq, this.typeSeq);
        return {
            rawSeq: seq,
            encodedSeq,
            compressedSeq: this.typeSeq === SEQUENCE_TYPE.PROTEIN ? compressToDayhoff(encodedSeq) : encodedSeq,
            type: this.typeSeq,
        };
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
                    this.addSequences(seqArr);
                } else {
                    return reject('Array of sequences expected');
                }
            } else if (this.sequences.length < 2) {
                return reject('At least 2 sequences are required.');
            }

            if (DEBUG) Log.add('Prepared sequences');

            setSeqType(this.sequences[0].type);
            setAlignmentParameters();

            if (DEBUG) Log.add('Get sequences type');

            if (this.sequences.length == 2) {
                if (Math.max(this.sequences[0].rawSeq.length, this.sequences[1].rawSeq.length) > 1600) {
                    let lAlignment = noalignPair(this.sequences[0], this.sequences[1]);
                    if (DEBUG) Log.summary();
                    return resolve(lAlignment);
                }
                let lResult = pairwiseAlignment(this.sequences[0], this.sequences[1], [0], [1]);

                if (DEBUG) Log.summary();



                return resolve(lResult.alignment);
            } else {
                let lResult = progressiveAlignment(this.sequences);

                if (DEBUG) Log.summary();

                return resolve(lResult);
            }
        });

        return msa;
    }

}
