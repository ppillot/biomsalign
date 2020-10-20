/**
 * @file profile
 * @author Paul Pillot <paul.pillot@libmol.org>
 * @license MIT
 * @copyright 2020
 */

import { getAlignmentParameters } from './params';
import { aaToNum } from './sequence';

export function sumOfPairsScorePP(profA: Profile, profB: Profile) {
    var score = 0,
        resProfNo = 0;

    for (var i = 0; i < profB.m_uResidueGroup; i++) {
        resProfNo = profB.m_uSortOrder[i];
        score += profB.m_fcCounts[resProfNo] * profA.m_AAScores[resProfNo];
    }
    score /= profB.nbSeq;
    //if (score<0) console.log(resNo,prof.m_uSortOrder, score)
    return score;
}

export function sumOfPairsScoreSP(resA: number, prof: Profile) {
    var score = 0,
        resProfNo = 0;

    score = prof.m_AAScores[resA];
    return score;
}

export type Profile = {
    m_bAllGaps: boolean;
    /** Amino Acids by order of frequency */
    m_uSortOrder: number[];
    /** Am. Ac. counts at this position */
    m_fcCounts: number[];
    m_LL: number;
    m_LG: number;
    m_GL: number;
    m_GG: number;
    /** Scores for each Am. Ac at this position */
    m_AAScores: number[];
    /** number of different residues */
    m_uResidueGroup: number;
    /** Frequence of occupency (non indels) */
    m_fOcc: number;
    /** Frequence of gap start */
    m_fcStartOcc: number;
    /** Frequence of gap close */
    m_fcEndOcc: number;
    /** Score for starting gap */
    m_scoreGapOpen: number;
    /** Score for closing gap */
    m_scoreGapClose: number;
    /** Number of sequences */
    nbSeq: number;
};
/*
 * Définition de l'objet ProfPos()
 */
class ProfPos {
    m_bAllGaps = false;
    m_uSortOrder = []; //acides aminés triés par ordre d'abondance
    m_fcCounts = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]; //effectif de chaque aa dans le profil à cette position
    m_LL = 1;
    m_LG = 0;
    m_GL = 0;
    m_GG = 0;
    m_AAScores = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]; //scores pour chaque acide aminé
    m_uResidueGroup = 0; //nb de résidus différents
    m_fOcc = 1; //fréquence d'occupation de la position (non indels)
    m_fcStartOcc = 0; //fréquence d'occupation par le début d'un gap
    m_fcEndOcc = 0;
    m_scoreGapOpen = 0; //score pour l'ouverture d'un gap à cette position
    m_scoreGapClose = 0; //score pour la fermeture d'un gap à cette position
    nbSeq = 0; //nb de séquences dans le profil
}

/**calcule le profil à partir d'un alignement
 * @msa : alignement reçu sous la forme d'un tableau contenant chaque séquence
 * @weight : tableau indiquant le "poids" de chaque séquence
 */
export function profileFromMSA(msa: string[], gapO: number, gapE: number) {
    //renvoie un tableau à 2D contenant pour chaque position du profil, le nombre d'acides aminés,
    //le nombre de indels dans un gap en ouverture, en fermeture ou en extension
    var longueur = msa[0].length, //taille totale de l'alignement (nb cols)
        aa = '', // acide aminé en code à 1 lettre
        na = 0; // acide aminé en n° (correspondant aux tables)
    var tab = Array(longueur); //tableau contenant le profil
    //var weight = (weight || (function(l){var t=[], q=1/l; for (var i=0;i<l;i++) {t.push(q)} return t})(msa.length)  ); //faire un tableau
    //tableau contenant la pondération de chaque séquence
    var uHydrophobicRunLength = 0, //calcul de la taille d'une fenêtre hydrophobe (si elle existe) pour éviter d'insérer des gaps au sein de celle-ci
        w = 1 / msa.length; //pondération d'une séquence
    //console.log( weight );
    const params = getAlignmentParameters();

    gapO = gapO || params.gapOP;
    gapE = gapE || params.gapEP;

    for (var col = 0; col < longueur; col++) {
        tab[col] = new ProfPos();
        var fGap = 0;

        for (var numSeq = 0, max = msa.length; numSeq < max; numSeq++) {
            tab[col].nbSeq = max;
            //w = weight[numSeq];
            aa = msa[numSeq][col];
            if (aa == '_') {
                fGap += w;
                if (col == 0) {
                    tab[col].m_fcStartOcc += w;
                } else if (msa[numSeq][col - 1] != '_') {
                    tab[col].m_fcStartOcc += w;
                }
                if (col == longueur - 1) {
                    tab[col].m_fcEndOcc += w;
                } else if (msa[numSeq][col + 1] != '_') {
                    tab[col].m_fcEndOcc += w;
                }
            } else {
                na = aaToNum[aa.charCodeAt(0)];
                tab[col].m_fcCounts[na]++;
                if (tab[col].m_fcCounts[na] == 1) {
                    tab[col].m_uSortOrder.push(na);
                    tab[col].m_uResidueGroup++;
                }
            }
        }
        tab[col].m_fOcc = 1 - fGap;
        tab[col].m_scoreGapOpen = gapO * (1 - tab[col].m_fcStartOcc);
        tab[col].m_scoreGapExtend = gapE * tab[col].m_fOcc;
        tab[col].m_scoreGapClose = gapO * (1 - tab[col].m_fcEndOcc);

        for (var aaNum = 0; aaNum < 20; aaNum++) {
            for (var i = 0; i < tab[col].m_uResidueGroup; i++) {
                var resProfNo = tab[col].m_uSortOrder[i];
                tab[col].m_AAScores[aaNum] +=
                    (tab[col].m_fcCounts[resProfNo] * params.scoringMatrix[aaNum][resProfNo]) / msa.length;
            }
        }
    }
    tab[0].m_scoreGapOpen = 0;
    tab[longueur - 1].m_scoreGapClose = 0;

    return tab;
}
