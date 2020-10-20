// /*
//  *	fonction calculant une matrice utilisable pour réaliser un dotplot
//  * prend comme paramètre un tableau de 1 ou deux séquences, ou bien un tableau contenant les n° de séquences déjà ajoutées
//  */
// m.calculeDotPlot = function (tabSeq) {
//     if (tabSeq.length == 1) return this.calculeDotPlot([tabSeq[0], tabSeq[0]]);

//     if (typeof tabSeq[0] == 'string') {
//         this.reset();
//         this.ajouteSequences(tabSeq);
//         this.params.init('', false);
//         tabSeq = [0, 1];
//     }
//     if (typeof tabSeq[0] == 'number') {
//         var seqA = this.sequences[tabSeq[0]],
//             seqB = this.sequences[tabSeq[1]];
//     } else return [];

//     var matrixDotplot = [],
//         vecteurDotPlot = [],
//         matrixScore = [],
//         vecteurScore = [],
//         l = seqA.seq.length,
//         m = seqB.seq.length,
//         diagL = Math.min(l, m),
//         score = 0,
//         lFen = this.params.largeurFenetre,
//         lDemiFen = (lFen / 2) | 0,
//         sNA = [],
//         sNB = [];

//     var tabDummy = (function (l) {
//         var t = [];
//         for (i = 0; i < l; i++) {
//             t.push(20);
//         }
//         return t;
//     })(lDemiFen);

//     sNA = tabDummy.concat(seqA.seqNum).concat(tabDummy);
//     sNB = tabDummy.concat(seqB.seqNum).concat(tabDummy);

//     for (var i = 0, imax = sNA.length; i < imax; i++) {
//         matrixDotplot[i] = [];
//         vecteurScore = [];
//         for (var j = 0, jmax = sNB.length; j < jmax; j++) {
//             vecteurScore.push(this.utils.substitutionNum(sNA[i], sNB[j]));
//         }
//         matrixScore.push(vecteurScore);
//     }
//     console.log(matrixScore);

//     //parcours commençant par la diagonale haut-gauche, vers la droite
//     for (i = lDemiFen; i < diagL + lDemiFen; i++) {
//         //les n° des diagonales partant du haut
//         score = 0;
//         for (j = lDemiFen; j < diagL + lDemiFen - i; j++) {
//             score =
//                 score + matrixScore[i + j + lDemiFen][j + lDemiFen] - matrixScore[i + j - lDemiFen][j - lDemiFen];
//             matrixDotplot[i - lDemiFen][j - lDemiFen] = score;
//         }
//     }

//     //parcours des diagonales de haut en bas
//     for (j = lDemiFen; j < diagL + lDemiFen; j++) {
//         //les n° des diagonales partant du haut
//         score = 0;
//         for (i = lDemiFen; i < diagL + lDemiFen - j; i++) {
//             score =
//                 score + matrixScore[i + j + lDemiFen][j + lDemiFen] - matrixScore[i + j - lDemiFen][j - lDemiFen];
//             matrixDotplot[i - lDemiFen][j - lDemiFen] = score;
//         }
//     }

//     return matrixDotplot;
// };
