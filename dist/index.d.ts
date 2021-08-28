/**
 * @file index
 * @author Paul Pillot <paul.pillot@libmol.org>
 * @license MIT
 * @copyright 2020-2021
 */
declare type TAlignMethod = 'diag' | 'complete' | 'auto';
declare type TAlignSequenceType = 'amino' | 'nucleic' | 'auto';
export declare type TAlignOpt = {
    /**
     * Gap open penalty (generally a negative number)
     * @default -300|-400 set based on the sequence type
     */
    gapopen: number;
    /**
     * Gap extend penalty (generally a negative number)
     * @default -11|-60 set based on the sequence type
     */
    gapextend: number;
    /**
     * Set true to report debug information to the JS console
     * @default false
     */
    debug: boolean;
    /**
     * Alignment method. `auto` is the default value. For sequences length up
     * to 1600 residues, a `complete` Needleman-Wunsch global alignment is
     * performed. Beyond, a diagonal finding algorithm  (`diag`) is run first
     * and global alignment is restricted to regions between diagonals.
     * `complete` and `diag` parameters force the use of one of those
     * strategies.
     * @default 'auto'
     */
    method: TAlignMethod;
    /**
     * Sequence type
     * `auto` performs a detection based on sequences content
     * @default 'auto'
     */
    type: TAlignSequenceType;
    /**
     * Substitution scoring matrix. Residues are sorted lexically by their 1
     * letter code (e.g. "A, C, D, E..." for proteins, not "Ala, Arg, Asn,...")
     * which is a discrepancy with BLOSUM matrices found on NCBI.
     * @default set based on sequence type
     */
    matrix: number[][];
    /**
     * Char that is used as the symbol for gaps in sequences
     * @default '-'
     */
    gapchar: string;
};
declare class BioMSAClass {
    private sequences;
    private typeSeq;
    private config;
    /**
     * Add sequence to the internal list of sequences of the aligner.
     * All sequences must be of same type.
     * @param {(string|string[])} seq sequence or array of sequences to add
     * @returns void
     * @throws Will throw an exception if sequence added is not of same type as
     *      the other internal sequences.
     */
    private addSequence;
    /**
     * Reset the internal sequence store content
     */
    private reset;
    private setDefaultConfiguration;
    private setUserConfiguration;
    /**
     * Main alignment function. Takes an array of sequences (strings) as an
     * argument or aligns sequences that were previously added to internal
     * store.
     * This function returns a promise.
     * @param seqArr
     * @returns { Promise<string[]> }
     */
    align(seqArr: string[], opt?: Partial<TAlignOpt>): Promise<string[]>;
}
declare const BioMSA: BioMSAClass;
export default BioMSA;
