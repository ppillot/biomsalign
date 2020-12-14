import { getAlignmentParameters, setAlignmentParameters, setSeqType } from './utils/params';
import { profileFromMSA, sumOfPairsScorePP, sumOfPairsScoreSP } from './utils/profile';
/**
 * @file biomsa
 * @author Paul Pillot <paul.pillot@libmol.org>
 * @license MIT
 * @copyright 2020
 */

import {
    compressToDayhoff,
    distanceMatrix,
    encodeSeqToNum,
    getSequenceType,
    SEQUENCE_TYPE,
    TSequence,
} from './utils/sequence';
import { isLeafNode, makeTree, LeafNode, clustalWeights, InternalNode, NODE_TYPE } from './utils/tree';

enum TRACE_BACK {
    MATCH     = 0,
    DEL       = 1,
    INS       = 1 << 2,
    MATCH2INS = 1 << 3,
    MATCH2DEL = 1 << 4
};

class BioMSA {
    private sequences: TSequence[];
    private typeSeq: SEQUENCE_TYPE = SEQUENCE_TYPE.PROTEIN;

    constructor() {
        this.sequences = [];
    }

    addSequences(seq: string | string[]) {
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

            if (this.sequences.length == 2) {
                let lResult = this._pairwiseAlignment(this.sequences[0], this.sequences[1], [0], [1]);
                return resolve(lResult.alignment);
            } else {
                return resolve(this._progressiveAlignment(this.sequences));
            }
        });

        return msa;
    }



    private _pairwiseAlignment (
        seqA: TSequence,
        seqB: TSequence,
        idA: number[],
        idB: number[]
    ) {
        setAlignmentParameters();
        const params = getAlignmentParameters();

        var lSeqALen = seqA.rawSeq.length,
            lSeqBLen = seqB.rawSeq.length,
            lMatch = 0,      // match score
                //
                // Textbook approach considers 3 n x m matrices for holding
                // the match / delete / insert values computed. These could be
                // reduces to only 2 vectors and 2 scalars
                //
            lMatchArr = [],     // Match column
            lDelArr = [],       // Delete column
            lPrevMatch = 0, // Scalar value containing match value computed at
                            // the previous iteration in the same column (j-1).
                            // Used as a rotation variable.
            lLastInsert = 0,

            tbM = new Uint8Array(lSeqALen * (lSeqBLen + 1)), // Trace back matrix
            lAlignment: string[] = new Array(2),
            sA = seqA.encodedSeq,
            sB = seqB.encodedSeq,

            gapOpenA = 0,
            gapOpenB = 0,
            gapExtentA = 0,
            gapExtentB = 0,
            tb = 0,         // trace back value. It's a bitmask of match + del + ins

            matrix = [],    // memoization of the scoring matrix row for the
                            // i-th amino acid.

            GAP_OPEN = params.gapOP;

            // Fill in first columns in the matrix
            // INITIALISATION
            // half gap open penalty is a tweak to favor gaps at the opening
            // (also at the end) of the sequence, the rational being that
            // it is common to have ragged alignments between distant sequences.
            // (This optimization comes from MAFFT and is applied in MUSCLE).

            lMatchArr[0] = 0;
            lDelArr[0] = -Infinity;
            for (var j = 1; j <= lSeqBLen; j++) {
                lMatchArr[j] = GAP_OPEN / 2;
                lDelArr[j] = -Infinity;
            }

            // Main DP routine

            for (var i = 1; i <= lSeqALen; i++) {

                lPrevMatch = GAP_OPEN / 2;
                lLastInsert = -Infinity;
                matrix = params.scoringMatrix[sA[i - 1]];   // memoize.

                for (j = 1; j <= lSeqBLen; j++) {
                    tb = TRACE_BACK.MATCH;


                    // Delete [i,j] score computation
                    //
                    //  lMatchArr was not cleared -> contains the computations
                    //      from the previous column (i-1).
                    //  lDelArr too.
                    //
                    //  - Gap extension is not computed due to the constant
                    //      added to the scoring matrix (see addCentre).


                    gapOpenA = lMatchArr[j] + GAP_OPEN;
                    gapExtentA = lDelArr[j];
                    if (j === lSeqBLen) {
                        gapOpenA -= GAP_OPEN / 2; // Terminal gap correction!
                    }

                    if (gapOpenA >= gapExtentA) {
                        lDelArr[j] = gapOpenA;
                    } else {
                        // No change on gap extent: lDelArr[j] += 0;
                        tb += TRACE_BACK.DEL;
                    }


                    // Insert [i,j] score computation
                    //  - lPrevMatch is the value that has just been computed
                    //      for the match at this column. We can't use lMatchArr
                    //      to store it yet, because for perf reason we are
                    //      saving space on this vector.
                    //  - ExtentB is still the score for the previous insert
                    //      value. Due to the addCentre() trick, no need to add
                    //      gap extension penalty here.

                    gapOpenB = lPrevMatch + GAP_OPEN;
                    gapExtentB = lLastInsert;

                    if (i === lSeqALen) {
                        gapOpenB -= GAP_OPEN / 2; // Terminal gap correction!
                    }

                    if (gapOpenB >= gapExtentB) {
                        lLastInsert = gapOpenB;
                    } else {
                        // No change on gap extend: lLastInsert += 0;
                        tb += TRACE_BACK.INS;
                    }


                    // Match [i,j] score computation
                    lMatch = lMatchArr[j - 1] + matrix[sB[j - 1]];

                    lMatchArr[j - 1] = lPrevMatch;  // var rotation

                    if (lMatch >= lLastInsert) {
                        if (lMatch >= lDelArr[j]) { // match is optimal
                            lPrevMatch = lMatch;
                        } else {                    // delete is optimal
                            lPrevMatch = lDelArr[j];
                            tb += TRACE_BACK.MATCH2DEL;
                        }

                    } else {
                        if (lLastInsert >= lDelArr[j]) { // insert is optimal
                            lPrevMatch = lLastInsert;
                            tb += TRACE_BACK.MATCH2INS;
                        } else {                         // delete is optimal
                            lPrevMatch = lDelArr[j];
                            tb += TRACE_BACK.MATCH2DEL;
                        }
                    }

                    // Store trace-back bits
                    tbM[i * lSeqBLen + j] = tb;

                }

                // fix match vector for its last value
                lMatchArr[lSeqBLen] = lPrevMatch;

            }

            var score = Math.max(lMatch, lLastInsert, lDelArr[j - 1]);

            // Traceback
            // Move backwards starting from the last optimal scoring position
            // in the TB matrix.

            i = lSeqALen;
            j = lSeqBLen;
            var lIdx = lSeqALen * lSeqBLen + lSeqBLen;
            lAlignment[0] = seqA.rawSeq;
            lAlignment[1] = seqB.rawSeq;

            // current matrix is either M (0), D (1) or I(2). Let's have a look
            // at the last value to see if the optimum is coming from DEL
            // (bit value 4) or INS (bit value 8) or Match (no bit value)
            var lCurrentMatrix = (tb & 12) >> 2,
                val = 0;

            while ((i > 0) && (j > 0)) {
                lIdx = i * lSeqBLen + j;
                if (lCurrentMatrix === 0) {
                    val = tbM[lIdx] >> 2;
                    if (val === 0) { //-->Match
                        i--;
                        j--;
                    }
                    // other cases --> run the loop once more to enter the
                    // following block.
                } else {
                    if (lCurrentMatrix === 2) { //-->Ins
                        lAlignment[0] = lAlignment[0].substring(0, i) + '-' + lAlignment[0].substring(i);
                        val = tbM[lIdx] & 2;
                        j--;
                    } else { //1 --> Del
                        lAlignment[1] = lAlignment[1].substring(0, j) + '-' + lAlignment[1].substring(j);
                        val = tbM[lIdx] & 1;
                        i--;
                    }
                }
                lCurrentMatrix = val;
            }

            // Finish sequences edit by appending the remaining symbols
            if (j === 0) {
                lAlignment[1] = '-'.repeat(i) + lAlignment[1];
            } else { //i==0
                lAlignment[0] = '-'.repeat(j) + lAlignment[0];
            }


            return {
                alignment: lAlignment,
                tSeqNames: idA.concat(idB),
                score: score
            };
    };

    private _progressiveAlignment(seq: TSequence[]) {
        //calcul des ressemblances par paires (cf. k-mers binary Muscle)
        //console.time("distances");
        var matriceDistances = distanceMatrix(seq);
        //console.timeEnd("distances");
        //console.log("matrice distances : ",matriceDistances);

        //construction d'un arbre
        //console.time("arbre guide");
        var tree = makeTree(matriceDistances.slice(), seq);
        //console.timeEnd("arbre guide")
        //console.log(tree);

        //Alignement selon l'arbre
        //console.time("alignement");

        //parcours des noeuds de l'arbre
        var nbNoeuds = tree.length;

        for (var i = 0; i < nbNoeuds; i++) {
            if ('seq' in tree[i]) continue; //on se déplace jusqu'au premier noeud

            var noeudA = tree[tree[i]['childA']];
            var noeudB = tree[tree[i]['childB']];

            if (isLeafNode(noeudA)) {
                //A est une séquence
                if (isLeafNode(noeudB)) {
                    //B est une séquence

                    tree[i].msa = this._pairwiseAlignment(noeudA.seq, noeudB.seq);
                } else {
                    //B est un alignement

                    tree[i].msa = this._MSASeqAlignment(noeudA.seq, noeudB.msa);
                }
            } else {
                //A et B sont des alignements
                tree[i].msa = this._MSAMSAAlignment(noeudA.msa, noeudB.msa);
            }
        }

        //console.timeEnd("alignement");
        //console.log(msa);

        return tree[tree.length - 1].msa;
    }

    /**
     * Profile to sequence Alignment (see MUSCLE for reference)
     * @param {LeafNode} nodeA a node object containing one sequence and its weight
     * @param {InternalNode} nodeB a node object containing a multiple sequence alignment and its weight
     * @return {array) array containing aligned sequences
     */
    private _MSASeqAlignment(
        nodeA: LeafNode,
        nodeB: InternalNode,
        tSeqANames: number[],
        tSeqBNames: number[]
    ) {



        var seqA = nodeA.seq,
            //wA = noeudA.weight, //probably unnecessary
            msaB = nodeB.msa,
            m = seqA.rawSeq.length,
            n = msaB[0].length,
            lMatch = 0.0,       // match score
            lMatchArr = [],
            lDelArr = [],
            tbM = new Uint8Array(n * (m + 1)),
            lAlignment: string[] = [],
            gapOpenA = 0.0,
            gapOpenB = 0.0,
            gapExtendA = 0.0,
            gapExtendB = 0.0,
            tb = 0,
            lastInsert = 0.0,
            prevLastInsert = 0.0,
            prevMatch = 0.0,
            prevDelete = 0.0,
            deletej_1 = 0.0,
            inserti_1 = 0.0,
            //M0 = 0.0,
            //profBGapEP = 0.0,
            profBGapOP = 0.0, //gap open penalty for profB at pos i-1
            profBGapCP = 0.0, //gap close penalty for profB at pos i-1
            profBAAScores = []; //substitution scores for profB at pos i-1

        const params = getAlignmentParameters();

        var gapOP = params.gapOP;
        var seqAGapOP = gapOP / 2;
        var seqAGapCP = gapOP / 2;

        // sequence A as number vector
        var sA = seqA.encodedSeq;

        // MSA B as profile
        var profB = profileFromMSA(msaB, gapOP, nodeB.tabWeight);

        //in this alignment for performance reasons, profB is explored from i=0 to n
        //and seqA is explored from j=0 to m

        //remplissage de la première ligne et de la première colonne de la matrice
        //INITIALISATION
        lMatchArr[0] = 0; //_utils.sumOfPairsScoreSP(sA[0], profB[0]);
        lDelArr[0] = -Infinity;
        for (var j = 1; j <= m; j++) {
            lMatchArr[j] = (seqAGapOP + seqAGapCP) / 2; //gap open at j=0 + gap close at j-1
            lDelArr[j] = -Infinity;
        }

        //remplissage des matrices
        //RECURSION
        //$log.debug('seqA: ', tSeqANames, ' penalités gapEP: ', gapEP, ' gapOP: ', gapOP);
        //M0 = Match[0] + profB[0].m_ScoreGapOpen; //est deux fois plus faible qu'une pénalité d'ouverture normale (début de séquence)

        for (var i = 1; i <= n; i++) {

            profBGapOP = profB[i - 1].m_ScoreGapOpen;
            profBGapCP = profB[i - 1].m_ScoreGapClose;
            profBAAScores = profB[i - 1].m_AAScores;

            //M0 += profBGapEP / 2;

            prevMatch = 0; //(i===1) ? 0 : profB[0].m_ScoreGapOpen + profBGapCP/2 ;
            lastInsert = -Infinity;
            lDelArr[0] = profB[0].m_ScoreGapOpen;

            for (j = 1; j <= m; j++) {
                tb = 0;

                //Delete i,j score computation
                gapOpenA = lMatchArr[j] + profBGapOP; //
                gapExtendA = prevDelete = lDelArr[j];
                if (j === m) { //terminal penalties are halved
                    gapOpenA -= profBGapOP / 2;
                }

                if (gapOpenA >= gapExtendA) {
                    lDelArr[j] = gapOpenA;
                } else {
                    //Delete[j] = gapExtendA;
                    tb += 1;
                }

                //Insert i,j score computation
                gapOpenB = prevMatch + seqAGapOP;
                gapExtendB = prevLastInsert = lastInsert;
                if (i === n) {
                    gapOpenB -= seqAGapOP / 2;
                }

                if (gapOpenB >= gapExtendB) {
                    lastInsert = gapOpenB;
                } else {
                    //lastInsert = gapExtendB;
                    tb += 2;
                }

                //Match i,j score computation
                deletej_1 = lDelArr[j] + profBGapCP; //it should be prev
                inserti_1 = prevLastInsert + seqAGapCP;
                lMatch = lMatchArr[j - 1] + profBAAScores[sA[j - 1]];

                if (j === 1) { //terminal penalties are halved
                    deletej_1 -= profBGapCP / 2;
                }
                if (i === 1) {
                    inserti_1 -= seqAGapCP / 2;
                }
                if ((j === m) && (i === n)) {
                    deletej_1 -= profBGapCP / 2;
                    inserti_1 -= seqAGapCP / 2;
                }

                lMatchArr[j - 1] = prevMatch;

                if (lMatch >= inserti_1) {
                    if (lMatch >= deletej_1) { //match is optimal
                        prevMatch = lMatch;
                    } else { //delete is optimal
                        prevMatch = deletej_1;
                        tb += 4;
                    }

                } else {
                    if (inserti_1 >= deletej_1) { //insert is optimal
                        prevMatch = inserti_1;
                        tb += 8;
                    } else { //delete is optimal
                        prevMatch = deletej_1;
                        tb += 4;
                    }
                }
                tbM[i * m + j] = tb;

                /* if (j < 4 && i < 5) {
                                        $log.debug('profB:'+ i, ' seqA:'+ j + '/'+sA[j - 1],
                                            'gapOpenB: '+ gapOpenB,
                                            'gapExtentB: '+ gapExtendB,
                                            'gapOpenA: '+ gapOpenA,
                                            'gapExtentA: '+ gapExtendA,
                                            'deletej_1: '+ deletej_1,
                                            'match: '+ match,
                                            'lastInsert: '+ lastInsert,
                                            'delete: '+ Delete[j],
                                            'tb: '+ tb
                                        );
                                    }
                        /* */
            }
            lMatchArr[m] = prevMatch;
        }
        var score = Math.max(lMatch, lastInsert, lDelArr[j - 1]);

        //traceback
        i = n;
        j = m;
        var indice = n * m + m,
            k = 0,
            tSeqNames = [];

        lAlignment.push( seqA.rawSeq , ...msaB);
        tSeqNames = [...tSeqANames, ...tSeqBNames];

        var matriceActive = (tb & 12) >> 2,
            value = 0;
        while ((i > 0) && (j > 0)) {
            indice = i * m + j;
            //$log.debug('i:'+i,'j:'+j,'matrice active:'+value,'tb:'+tbM[indice]);

            if (matriceActive === 0) {
                value = tbM[indice] >> 2;
                if (value === 0) {
                    i--;
                    j--;
                }
            } else {
                if (matriceActive === 2) {
                    for (k = 1; k < lAlignment.length; k++) {
                        lAlignment[k] = lAlignment[k].substring(0, i) + '-' + lAlignment[k].substring(i);
                    }
                    value = tbM[indice] & 2; //@@PP same as in _pairwise ?
                    j--;

                } else {
                    lAlignment[0] = lAlignment[0].substring(0, j) + '-' + lAlignment[0].substring(j);
                    value = tbM[indice] & 1;
                    i--;
                }
            }
            matriceActive = value;
        }
        //fin des profils par ajout direct de la fin de séquence
        if (i === 0) {
            var p = '-'.repeat(j);
            for (k = 1; k < lAlignment.length; k++) {
                lAlignment[k] = p + lAlignment[k];
            }
        } else { //j==0
            lAlignment[0] = '-'.repeat(i) + lAlignment[0];
        }
        return {
            alignment: lAlignment,
            tSeqNames: tSeqNames,
            score: score
        };
    }

    /*
     * fonction alignant 2 profils (alignement multiple x alignement multiple)
     */

    private _MSAMSAAlignment(msaA: string[], msaB: string[], gE?: number, gO?: number) {
        if (msaA.length < msaB.length) {
            //permutation des alignements : B reçoit le plus petit, réduit le nb d'itérations ensuite
            var msa = msaB.slice();
            msaB = msaA.slice();
            msaA = msa;
        }

        var n = msaA[0].length;
        var m = msaB[0].length;
        var Match: number[] = [],
            Delete: number[] = [],
            MatchPrev: number[] = [],
            tbM: number[][] = [],
            alignement = [];
        var gapOpenA, gapOpenB, gapExtendA, gapExtendB, tb, mP, dP, lastInsert, acc, M0;

        setAlignmentParameters(true);
        const params = getAlignmentParameters();

        var gapEP = gE || params.gapEP,
            gapOP = gO || params.gapOP;

        //conversion du MSA en profil
        var profB = profileFromMSA(msaB, gapOP, gapEP);
        var profA = profileFromMSA(msaA, gapOP, gapEP);

        //remplissage de la première ligne et de la première colonne de la matrice
        //INITIALISATION
        MatchPrev[0] = sumOfPairsScorePP(profA[0], profB[0]);
        Delete[0] = -Infinity;
        for (var j = 1; j <= m; j++) {
            MatchPrev[j] = MatchPrev[0] + j * gapEP + profB[0].m_scoreGapOpen;
            Delete[j] = -Infinity;
        }

        M0 = profA[0].m_scoreGapOpen;
        //remplissage des matrices
        //RECURSION
        for (var i = 1; i <= n; i++) {
            M0 += profA[i - 1].m_scoreGapExtend;
            Match[0] = M0;
            tbM[i] = [];
            lastInsert = -Infinity;
            Delete[0] = MatchPrev[0] + profA[i - 1].m_scoreGapOpen;

            for (var j = 1; j <= m; j++) {
                tb = 0;

                //Delete i,j score computation
                gapOpenA = MatchPrev[j] + profA[i - 1].m_scoreGapOpen; //
                gapExtendA = Delete[j] + profA[i - 1].m_scoreGapExtend;
                if (gapOpenA >= gapExtendA) {
                    Delete[j] = gapOpenA;
                } else {
                    Delete[j] = gapExtendA;
                    tb += 1;
                }
                //Insert i,j score computation
                gapOpenB = Match[j - 1] + profB[j - 1].m_scoreGapOpen;
                gapExtendB = lastInsert + profB[j - 1].m_scoreGapExtend;

                if (gapOpenB >= gapExtendB) {
                    lastInsert = gapOpenB;
                } else {
                    lastInsert = gapExtendB;
                    tb += 2;
                }

                //Match i,j score computation
                const match = MatchPrev[j - 1] + sumOfPairsScorePP(profA[i - 1], profB[j - 1]);
                if (match >= lastInsert) {
                    if (match >= Delete[j]) {
                        //match is optimal
                        Match[j] = match;
                    } else {
                        //delete is optimal
                        Match[j] = Delete[j];
                        tb += 4;
                    }
                } else {
                    if (lastInsert >= Delete[j]) {
                        //insert is optimal
                        Match[j] = lastInsert;
                        tb += 8;
                    } else {
                        //delete is optimal
                        Match[j] = Delete[j];
                        tb += 4;
                    }
                }
                tbM[i][j] = tb;
            }
            MatchPrev = Match.slice();
        }

        //traceback
        var i = n,
            j = m;
        alignement = msaA.concat(msaB);

        var matriceActive = 0,
            value = 0;
        while (i > 0 && j > 0) {
            if (matriceActive == 0) {
                value = tbM[i][j] >> 2;
                if (value == 0) {
                    i--;
                    j--;
                }
            } else {
                if (matriceActive == 2) {
                    for (var k = 0; k < msaA.length; k++) {
                        alignement[k] = alignement[k].substring(0, i) + '_' + alignement[k].substring(i);
                    }
                    value = tbM[i][j] & 2;
                    j--;
                } else {
                    for (var k = msaA.length; k < alignement.length; k++) {
                        alignement[k] = alignement[k].substring(0, j) + '_' + alignement[k].substring(j);
                    }
                    value = tbM[i][j] & 1;
                    i--;
                }
            }
            matriceActive = value;
        }
        //fin des profils par ajout direct de la fin de séquence
        if (j == 0) {
            var p = '_'.repeat(i);
            for (var k = msaA.length; k < alignement.length; k++) {
                alignement[k] = p + alignement[k];
            }
        } else {
            //i==0
            var p = '_'.repeat(j);
            for (var k = 0; k < msaA.length; k++) {
                alignement[k] = p + alignement[k];
            }
        }
        return alignement;
    }
}
