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
import { isLeafNode, makeTree, TreeNode } from './utils/tree';

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
        const msa = new Promise<string[]>((resolve, reject) => {
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
                return resolve(this._pairwiseAlignment(this.sequences[0], this.sequences[1]));
            } else {
                return resolve(this._progressiveAlignment(this.sequences));
            }
        });

        return msa;
    }

    private _pairwiseAlignment = function (seqA: TSequence, seqB: TSequence) {
        setAlignmentParameters(false);
        const params = getAlignmentParameters();

        var n = seqA.rawSeq.length;
        var m = seqB.rawSeq.length;
        var Match: number[] = [],
            Delete: number[] = [],
            MatchPrev: number[] = [],
            tbM: number[][] = [],
            alignement: string[] = [],
            sA = [],
            sB = [];
        var gapOpenA, gapOpenB, gapExtentA, gapExtentB, tb, mP, dP, lastInsert, acc;
        var gapEP = params.gapEP,
            gapOP = params.gapOP;

        //conversion des séquences en tableaux de valeurs 1-20
        sA = seqA.encodedSeq;
        sB = seqB.encodedSeq;

        //remplissage de la première ligne et de la première colonne de la matrice
        //INITIALISATION
        MatchPrev[0] = 0;
        Delete[0] = -Infinity;
        for (var j = 1; j <= m; j++) {
            MatchPrev[j] = j * gapEP + gapOP;
            Delete[j] = -Infinity;
        }

        //remplissage des matrices
        //RECURSION
        for (var i = 1; i <= n; i++) {
            Match[0] = i * gapEP + gapOP;
            tbM[i] = [];
            lastInsert = -Infinity;

            for (var j = 1; j <= m; j++) {
                tb = 0;

                //Delete i,j score computation
                gapOpenA = Match[j] + gapOP; //these values have not yet been updated, they reflect the previous column's values
                gapExtentA = Delete[j] + gapEP;
                if (gapOpenA >= gapExtentA) {
                    Delete[j] = gapOpenA;
                } else {
                    Delete[j] = gapExtentA;
                    tb += 1;
                }
                //Insert i,j score computation
                gapOpenB = Match[j - 1] + gapOP;
                gapExtentB = lastInsert + gapEP;

                if (gapOpenB >= gapExtentB) {
                    lastInsert = gapOpenB;
                } else {
                    lastInsert = gapExtentB;
                    tb += 2;
                }

                //Match i,j score computation
                let match = MatchPrev[j - 1] + params.scoringMatrix[sA[i - 1]][sB[j - 1]];
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
        alignement[0] = seqA.rawSeq;
        alignement[1] = seqB.rawSeq;
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
                    alignement[0] = alignement[0].substring(0, i) + '_' + alignement[0].substring(i);
                    value = tbM[i][j] & 2;
                    j--;
                } else {
                    alignement[1] = alignement[1].substring(0, j) + '_' + alignement[1].substring(j);
                    value = tbM[i][j] & 1;
                    i--;
                }
            }
            matriceActive = value;
        }
        //fin des profils par ajout direct de la fin de séquence
        if (j == 0) {
            alignement[1] = '_'.repeat(i) + alignement[1];
        } else {
            //i==0
            alignement[0] = '_'.repeat(j) + alignement[0];
        }
        return alignement;
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

    private _MSASeqAlignment(seqA: TSequence, msaB: string[], gE?: number, gO?: number) {
        var n = seqA.rawSeq.length;
        var m = msaB[0].length;
        var Match: number[] = [],
            Delete: number[] = [],
            MatchPrev: number[] = [],
            tbM: number[][] = [],
            alignement = [];
        var gapOpenA, gapOpenB, gapExtendA, gapExtendB, tb, mP, dP, lastInsert, acc;
        setAlignmentParameters(true);
        const params = getAlignmentParameters();

        var gapEP = gE || params.gapEP,
            gapOP = gO || params.gapOP;

        //conversion des séquences en tableaux de valeurs 1-20
        var sA = seqA.encodedSeq;

        //conversion du MSA en profil
        var profB = profileFromMSA(msaB, gapOP, gapEP);

        //remplissage de la première ligne et de la première colonne de la matrice
        //INITIALISATION
        MatchPrev[0] = sumOfPairsScoreSP(sA[0], profB[0]);
        Delete[0] = -Infinity;
        for (var j = 1; j <= m; j++) {
            MatchPrev[j] = MatchPrev[0] + j * gapEP + profB[0].m_scoreGapOpen;
            Delete[j] = -Infinity;
        }

        //remplissage des matrices
        //RECURSION
        for (var i = 1; i <= n; i++) {
            Match[0] = i * gapEP + gapOP;
            tbM[i] = [];
            lastInsert = -Infinity;
            Delete[0] = MatchPrev[0] + gapOP;

            for (var j = 1; j <= m; j++) {
                tb = 0;

                //Delete i,j score computation
                gapOpenA = MatchPrev[j] + gapOP; //
                gapExtendA = Delete[j] + gapEP;
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
                const match = MatchPrev[j - 1] + sumOfPairsScoreSP(sA[i - 1], profB[j - 1]);
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
        alignement[0] = seqA.rawSeq;
        alignement = alignement.concat(msaB);

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
                    alignement[0] = alignement[0].substring(0, i) + '_' + alignement[0].substring(i);
                    value = tbM[i][j] & 2;
                    j--;
                } else {
                    for (var k = 1; k < alignement.length; k++) {
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
            for (var k = 1; k < alignement.length; k++) {
                alignement[k] = p + alignement[k];
            }
        } else {
            //i==0
            alignement[0] = '_'.repeat(j) + alignement[0];
        }
        return alignement;
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
