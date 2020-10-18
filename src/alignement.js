/*
 * Bibliothèque JS pour réaliser un alignement multiple de séquences par la méthode progressive
 */

var MultAlign = (function() {
	var m = {};
	
	m.sequences = [];
	m.utils = {};
	m.params = {};
	
	m.Sequence = function(seq) {
		this.seq = seq;
		this.seqNum = [];
		this.seqCompresse = "";
	};
	
    
	
	return m;
		
})(document);

/*
 * fonctions principales
 */

(function(m) {
	/*
	 * fonction renvoyant un alignement multiple à partir d'un tableau de séquences
	 * seq : tableau à 1 dimension contenant les séquences à aligner. facultatif si m.sequences est rempli
	 */
	m.aligne = function(seq) {
		
		var msa = []; //alignement multiple : tableau de séquences avec les indels
		
		if (typeof(seq)!="undefined") {
		    this.reset();
			this.ajouteSequences(seq);
		}
		
		/*
		 * Selon le nombre de séquences on adopte des stratégies différentes :
		 * < 2 : pas d'alignement possible
		 * = 2 : alignement par paire
		 * > 2 : alignement progressif
		 */
		
		if (this.sequences.length<2) {
			alert("Un alignement nécessite au moins deux séquences");
			return seq;
		} else {
			if (this.sequences.length == 2) {
				msa = this._pairwiseAlignment (this.sequences[0],this.sequences[1]);
			} else {
				msa = this._alignementProgressif (this.sequences);
			}
		}
		
		return msa;

	};
	
	/*
	 *	fonction calculant une matrice utilisable pour réaliser un dotplot
	 * prend comme paramètre un tableau de 1 ou deux séquences, ou bien un tableau contenant les n° de séquences déjà ajoutées
	 */
	m.calculeDotPlot = function(tabSeq) {
		
		if (tabSeq.length==1) return this.calculeDotPlot([tabSeq[0],tabSeq[0]]);
		
		if (typeof(tabSeq[0])=="string") {
			this.reset();
			this.ajouteSequences(tabSeq);
			this.params.init('',false);
			tabSeq = [0,1]
		} 
		if (typeof(tabSeq[0])=="number") {
			var seqA = this.sequences[tabSeq[0]],
				seqB = this.sequences[tabSeq[1]]
		} else return [];
				
		var matrixDotplot = [],
			vecteurDotPlot = [],
			matrixScore = [],
			vecteurScore = [],
			l = seqA.seq.length,
			m = seqB.seq.length,
			diagL = Math.min(l,m),
			score = 0,
			lFen = this.params.largeurFenetre,
			lDemiFen = (lFen/2) | 0,
			sNA = [],
			sNB = [];
		
		var tabDummy = (function(l){
			var t= [];
			for (i=0; i<l; i++) {
				t.push(20);
			}
			return t;
		})(lDemiFen);
			
		sNA = tabDummy.concat(seqA.seqNum).concat(tabDummy);
		sNB = tabDummy.concat(seqB.seqNum).concat(tabDummy); 
		
		for (var i=0, imax = sNA.length; i<imax; i++ ) {
			matrixDotplot[i]=[];
			vecteurScore = [];
			for (var j=0, jmax = sNB.length; j<jmax ; j++) {
				vecteurScore.push( this.utils.substitutionNum( sNA[i] , sNB[j] ) );				
			}
			matrixScore.push(vecteurScore);
		}
		console.log(matrixScore);
		
		//parcours commençant par la diagonale haut-gauche, vers la droite
		for (i=lDemiFen; i < diagL+lDemiFen; i++) { //les n° des diagonales partant du haut
			score = 0;
			for (j=lDemiFen; j< diagL+lDemiFen-i ; j++) {
				score = score + matrixScore[i+j+lDemiFen][j+lDemiFen] - matrixScore[i+j-lDemiFen][j-lDemiFen];
				matrixDotplot[i-lDemiFen][j-lDemiFen] = score;
			}
		}
		
		//parcours des diagonales de haut en bas
		for (j=lDemiFen; j < diagL+lDemiFen; j++) { //les n° des diagonales partant du haut
			score = 0;
			for (i=lDemiFen; i< diagL +lDemiFen-j ; i++) {
				score = score + matrixScore[i+j+lDemiFen][j+lDemiFen] - matrixScore[i+j-lDemiFen][j-lDemiFen];
				matrixDotplot[i-lDemiFen][j-lDemiFen] = score;
			}
		}
		
		return matrixDotplot;
		
	}
	
	/*
	 * fonction ajoutant des séquences
	 */
	m.ajouteSequences = function(seq) {
		
		if (typeof(seq)=="object") {
			//c'est probablement un tableau de séquences qui a été envoyé
			for (var i=0; i<seq.length; i++) {
				this.ajouteSequences(seq[i]);
			}
		} else {
			if (typeof(seq)=="string") {
				// on vérifie que le type de séquence est le même que le type des séquences précédentes ajoutées (on n'aligne pas ADN avec Protéines)
				var type = this.utils.getSeqType(seq);
				if (type=="error") { return alert("erreur dans la séquence:" + seq); }
				
				if (this.params.typeSeq=="") {
					this.params.setSeqType(type);
					//this.params.typeSeq = type;
				} else {
					if (this.params.typeSeq!=type) {
						return alert("le type de séquence ne correspond pas aux autres séquences déjà ajoutées pour traitement");
					}
				}
				this.sequences.push(this.creationObjetSequence(seq,type));
				
			}
		}
		
	};
	
	/*
	 * fonction créant un objet séquence à partir d'une séquence
	 */
	m.creationObjetSequence = function(seq,type) {
	  
	  var s = new m.Sequence(seq);
	  s.seqNum = m.utils.seqToNum(seq,type);
	  s.seqCompresse = (type=="protein")? m.utils.compression(seq) : seq;
	  
	  return s;
	    
	};
	
	
	/*
	 * fonction réinitialisant l'objet
	 */
	m.reset = function() {
		this.sequences = [];
		this.params.typeSeq = 0;
		this.params.matrix = [];
		this.params.gapEP = 0;
		this.params.gapOP = 0;
	};
	
	/*
	 * fonction pour alignement par paires
	 */
	m._pairwiseAlignment = function (seqA, seqB) {
		
		this.params.init("",false);
		
		var n = seqA.seq.length;
		var m = seqB.seq.length;
		var Match = [], Delete = [], MatchPrev = [], tbM = [], profile = [], alignement = [], sA=[], sB=[];
		var gapOpenA, gapOpenB, gapExtentA, gapExtentB, tb, mP, dP, lastInsert, acc;
		var gapEP = this.params.gapEP,
			gapOP = this.params.gapOP;
		
		//conversion des séquences en tableaux de valeurs 1-20
		sA = seqA.seqNum;
		sB = seqB.seqNum;
		
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
				tb=0;
				
				//Delete i,j score computation
				gapOpenA = Match[j] + gapOP ;//these values have not yet been updated, they reflect the previous column's values
				gapExtentA = Delete[j] + gapEP;
				if (gapOpenA >= gapExtentA) {
					Delete[j] = gapOpenA;
				} else {
					Delete[j] = gapExtentA;
					tb += 1 ;
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
				match = MatchPrev[j - 1] + this.utils.substitutionNum(sA[i - 1], sB[j - 1]);
				if (match >= lastInsert) {
					if (match >= Delete[j]) {//match is optimal
						Match[j] = match;
					} else {//delete is optimal
						Match[j] = Delete[j];
						tb += 4;
					}
	
				} else {
					if (lastInsert >= Delete[j]) {//insert is optimal
						Match[j] = lastInsert;
						tb += 8;
					} else {//delete is optimal
						Match[j] = Delete[j];
						tb += 4;
					}
				}
				tbM[i][j] = tb;
			}
			MatchPrev = Match.slice();
		}
		
		//traceback
		var i = n, j = m;
		alignement[0] = seqA.seq;
		alignement[1] = seqB.seq;
		var matriceActive = 0, value = 0;
		while ((i > 0) && (j > 0)) {
			if (matriceActive ==0) {
				value = tbM[i][j]>>2;
				if (value == 0) {
					i--;
					j--;
				}
			} else {
				if (matriceActive == 2) {
					alignement[0] = alignement[0].substring(0,i) + '_' + alignement[0].substring(i);
					value = tbM[i][j] & 2;
					j--;
				} else {
					alignement[1] = alignement[1].substring(0,j) + '_' + alignement[1].substring(j);
					value = tbM[i][j] & 1;
					i--;
				}
			}
			matriceActive = value;
		}
		//fin des profils par ajout direct de la fin de séquence
		if (j==0) {
			alignement[1] = this.utils.pad('_',i) + alignement[1];
		} else { //i==0
			alignement[0] = this.utils.pad('_',j) + alignement[0];
		}
		return alignement;
	};
    
    /*
     * Alignement d'une séquence avec un profil
     */
	m._MSASeqAlignment = function(seqA, msaB, gE, gO) {
		console.log("MSASeq")
		var n = seqA.seq.length;
		var m = msaB[0].length;
		var Match = [], Delete = [], MatchPrev = [], tbM = [], profile = [], alignement = [];
		var gapOpenA, gapOpenB, gapExtendA, gapExtendB, tb, mP, dP, lastInsert, acc;
		
		this.params.init("",true);
		
		var gapEP = gE || this.params.gapEP, 
			gapOP = gO || this.params.gapOP;
		
		//conversion des séquences en tableaux de valeurs 1-20
		var sA = seqA.seqNum;
		
		//conversion du MSA en profil
		var profB = this.utils.profileFromMSA(msaB , gapOP , gapEP);
		
		//remplissage de la première ligne et de la première colonne de la matrice
		//INITIALISATION
		MatchPrev[0] = this.utils.sumOfPairsScoreSP(sA[0], profB[0]);
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
				tb=0;
				
				//Delete i,j score computation
				gapOpenA = MatchPrev[j] + gapOP; //
				gapExtendA = Delete[j] + gapEP;
				if (gapOpenA >= gapExtendA) {
					Delete[j] = gapOpenA;
				} else {
					Delete[j] = gapExtendA;
					tb += 1 ;
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
				match = MatchPrev[j - 1] + this.utils.sumOfPairsScoreSP(sA[i - 1], profB[j - 1]);
				if (match >= lastInsert) {
					if (match >= Delete[j]) {//match is optimal
						Match[j] = match;
					} else {//delete is optimal
						Match[j] = Delete[j];
						tb += 4;
					}
	
				} else {
					if (lastInsert >= Delete[j]) {//insert is optimal
						Match[j] = lastInsert;
						tb += 8;
					} else {//delete is optimal
						Match[j] = Delete[j];
						tb += 4;
					}
				}
				tbM[i][j] = tb;
			}
			MatchPrev = Match.slice();
		}
		
		//traceback
		var i = n, j = m;
		alignement[0] = seqA.seq;
		alignement = alignement.concat(msaB);
		
		var matriceActive = 0, value = 0;
		while ((i > 0) && (j > 0)) {
			if (matriceActive ==0) {
				value = tbM[i][j]>>2;
				if (value == 0) {
					i--;
					j--;
				}
			} else {
				if (matriceActive == 2) {
					alignement[0] = alignement[0].substring(0,i) + '_' + alignement[0].substring(i);
					value = tbM[i][j] & 2;
					j--;
				} else {
					for (var k=1; k<alignement.length; k++) {
						alignement[k] = alignement[k].substring(0,j) + '_' + alignement[k].substring(j);
					}
					value = tbM[i][j] & 1;
					i--;
				}
			}
			matriceActive = value;
		}
		//fin des profils par ajout direct de la fin de séquence
		if (j==0) {
			var p = this.utils.pad('_',i);
			for (var k=1; k<alignement.length; k++) {
				alignement[k] = p + alignement[k];
			}
		} else { //i==0
			alignement[0] = this.utils.pad('_',j) + alignement[0];
		}
		return alignement;
	};
	
	/*
	 * fonction alignant 2 profils (alignement multiple x alignement multiple)
	 */
	
	m._MSAMSAAlignment = function(msaA, msaB, gE, gO) {
		
		if (msaA.length<msaB.length) {
			//permutation des alignements : B reçoit le plus petit, réduit le nb d'itérations ensuite
			var msa = msaB.slice();
			msaB = msaA.slice();
			msaA = msa;
		}
		
		var n = msaA[0].length;
		var m = msaB[0].length;
		var Match = [], Delete = [], MatchPrev = [], tbM = [], profile = [], alignement = [];
		var gapOpenA, gapOpenB, gapExtendA, gapExtendB, tb, mP, dP, lastInsert, acc, M0;
		
		this.params.init("",true);
		
		var gapEP = gE || this.params.gapEP, 
			gapOP = gO || this.params.gapOP;
		
		//conversion du MSA en profil
		var profB = this.utils.profileFromMSA(msaB , gapOP , gapEP);
		var profA = this.utils.profileFromMSA(msaA , gapOP , gapEP);
		
		//remplissage de la première ligne et de la première colonne de la matrice
		//INITIALISATION
		MatchPrev[0] = this.utils.sumOfPairsScorePP(profA[0], profB[0]);
		Delete[0] = -Infinity;
		for (var j = 1; j <= m; j++) {
			MatchPrev[j] = MatchPrev[0] + j * gapEP + profB[0].m_scoreGapOpen;
			Delete[j] = -Infinity;
		}
	
		M0 = profA[0].m_scoreGapOpen;
		//remplissage des matrices
		//RECURSION
		for (var i = 1; i <= n; i++) {
			M0 += profA[i-1].m_scoreGapExtend;
			Match[0] = M0;
			tbM[i] = [];
			lastInsert = -Infinity;
			Delete[0] = MatchPrev[0] + profA[i-1].m_scoreGapOpen;
			
			for (var j = 1; j <= m; j++) {
				tb=0;
				
				//Delete i,j score computation
				gapOpenA = MatchPrev[j] + profA[i-1].m_scoreGapOpen; //
				gapExtendA = Delete[j] + profA[i-1].m_scoreGapExtend;
				if (gapOpenA >= gapExtendA) {
					Delete[j] = gapOpenA;
				} else {
					Delete[j] = gapExtendA;
					tb += 1 ;
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
				match = MatchPrev[j - 1] + this.utils.sumOfPairsScorePP(profA[i - 1], profB[j - 1]);
				if (match >= lastInsert) {
					if (match >= Delete[j]) {//match is optimal
						Match[j] = match;
					} else {//delete is optimal
						Match[j] = Delete[j];
						tb += 4;
					}
	
				} else {
					if (lastInsert >= Delete[j]) {//insert is optimal
						Match[j] = lastInsert;
						tb += 8;
					} else {//delete is optimal
						Match[j] = Delete[j];
						tb += 4;
					}
				}
				tbM[i][j] = tb;
			}
			MatchPrev = Match.slice();
		}
		
		//traceback
		var i = n, j = m;
		alignement = msaA.concat(msaB);
		
		var matriceActive = 0, value = 0;
		while ((i > 0) && (j > 0)) {
			if (matriceActive ==0) {
				value = tbM[i][j]>>2;
				if (value == 0) {
					i--;
					j--;
				}
			} else {
				if (matriceActive == 2) {
					for (var k=0; k<msaA.length; k++) {
						alignement[k] = alignement[k].substring(0,i) + '_' + alignement[k].substring(i);
					}
					value = tbM[i][j] & 2;
					j--;
				} else {
					for (var k=msaA.length; k<alignement.length; k++) {
						alignement[k] = alignement[k].substring(0,j) + '_' + alignement[k].substring(j);
					}
					value = tbM[i][j] & 1;
					i--;
				}
			}
			matriceActive = value;
		}
		//fin des profils par ajout direct de la fin de séquence
		if (j==0) {
			var p = this.utils.pad('_',i);
			for (var k=msaA.length; k<alignement.length; k++) {
				alignement[k] = p + alignement[k];
			}
		} else { //i==0
			var p = this.utils.pad('_',j);
			for (var k=0; k<msaA.length; k++) {
				alignement[k] = p + alignement[k];
			}
		}
		return alignement;
	};
	

	m._alignementProgressif = function(seq) {
		var ordreSeq=[];
		//conversion des séquences en alphabet compressé (cf. MAFFT)
		//console.time("compression");
		//var seqZ = this.utils.compression(seq.slice());
		//console.timeEnd("compression");
		//console.log(seqZ);
		
		//calcul des ressemblances par paires (cf. k-mers binary Muscle)
		//console.time("distances");
		var matriceDistances = this.utils.distances(seq);
		//console.timeEnd("distances");
		//console.log("matrice distances : ",matriceDistances);
		
		//construction d'un arbre
		//console.time("arbre guide");
		var tree = this.utils.arbre(matriceDistances.slice(),seq);
		//console.timeEnd("arbre guide")
		console.log(tree);
		
		//Alignement selon l'arbre
		//console.time("alignement");
			
		//parcours des noeuds de l'arbre
		var nbNoeuds = tree.length;
		
		for (var i=0; i<nbNoeuds; i++) {
			if (tree[i].seq.seq != "") continue; //on se déplace jusqu'au premier noeud
			
            var noeudA = tree[i]['childA'];
			var noeudB = tree[i]['childB'];
			
			if (tree[noeudA].seq.seq.length>0) { //A est une séquence
				if (tree[noeudB].seq.seq.length>0) { //B est une séquence
				    
                    tree[i].msa = this._pairwiseAlignment (tree[noeudA].seq , tree[noeudB].seq);

				} else { //B est un alignement
				    
					tree[i].msa = this._MSASeqAlignment (tree[noeudA].seq , tree[noeudB].msa);
					
				}
			} else { //A et B sont des alignements
				tree[i].msa = this._MSAMSAAlignment (tree[noeudA].msa , tree[noeudB].msa);
			}
			
		}
		
		//console.timeEnd("alignement");
		//console.log(msa);
		
		return tree[tree.length-1].msa;
	};
	
})(MultAlign);


/*
 * Fonctions utilitaires
 */
(function( u , params ) {
	
	/*
	 * fonction renvoyant un alignement multiple à partir d'un tableau de séquences
	 * seq : tableau à 1 dimension contenant les séquences à aligner. facultatif si m.sequences est rempli
	 */
	u.getSeqType = function(seq) {
		var acidesAmines = /[rndbeqyhilkmfpswyv]/gi;
		var notArn = /[^acgu]/gi;
		var notAdn = /[^acgt]/gi;
		
		if (acidesAmines.test(seq)) return "protein";
		else {
			if (!notArn.test(seq)) return "nucleic";
			else {
				if(!notAdn.test(seq)) return "nucleic";
				else return "error";
			}
		}	
	};
	
	/*
	 * fonction renvoyant une chaîne de nb fois le caractère s
	 */
	u.pad = function(s,nb) {
		var r ="";
		for(var i=0;i<nb;i++) {
			r+=s;
		}
		return r;
	};
	
	/*
	 * tables de conversion d'une chaîne de caractères en suite de nombres
	 * (plus performante)
	 */
	u.codeAANum = {'A':0,'C':1,'D':2,'E':3,'F':4,'G':5,'H':6,'I':7,'K':8,'L':9,'M':10,'N':11,'P':12,'Q':13,'R':14,'S':15,'T':16,'V':17,'W':18,'Y':19,'_':20};
	//u.codeNumAA = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "_"];
	
	u.codeNucNum = {'A':0,'C':1,'G':2,'T':3,'U':3} ;//ARN pouvant être aligné avec ADN, U est synonyme de T
	
	
	u.seqToNum = function(s,type){
		var t =[], tabConversion = [];
		var typeSeq = type || params.typeSeq;
		console.log(typeSeq);
		tabConversion = (type=="protein") ? this.codeAANum : this.codeNucNum;
		 
		for (var i=0,imax=s.length; i < imax; i++){
			t[i] = tabConversion[s[i]];
		}	
		return t;
	};
	
	/*
	 * conversion d'une séquence protéique en alphabet compressé
	 */
	u.compression = function(sequence) {
		var alphabetZ = {'A':'A','G':'A','P':'A','S':'A','S':'A','T':'A',
				'C':'B',
				'D':'C','E':'C','N':'C','Q':'C',
				'F':'D','W':'D','Y':'D',
				'H':'E','K':'E','R':'E',
				'I':'F','L':'F','M':'F','V':'F'					
		}; //Alphabet Dayhoff utilisé par MAFFT
		
		if (typeof(sequence)=="object") {
			for (i in sequence) {
				sequence[i] = this.compression(sequence[i]);
			}
		} else {
			var s='';
			for (var j=0 , k=sequence.length ; j<k ; j++) {
				s += alphabetZ[sequence[j]];
			}
			return s;
		}
		return sequence;
	};
	
	/*
	 * calcul d'une distance entre séquences à partir d'un partage de tuples
	 */
	u.distances = function(tabSeq) {
			
			function nbTuples(sequenceA,sequenceB) {
				//on parcourt la seqA en tranches de 6 caractères
				//on compte combien de ces tranches sont retrouvées en seqB
				var nbBinaire = 0;
				var seqA,seqB;
				
				if ( sequenceA.length < sequenceB.length ) {
					seqA = sequenceA;
					seqB = sequenceB;
				} else {
					seqA = sequenceB;
					seqB = sequenceA;
				}
				
				var nbTotalTuples = seqA.length - 5;
				
				for (var debTranche = 0; debTranche<nbTotalTuples; debTranche++ ){
					var txt = seqA.substr(debTranche,6);
					nbBinaire += (seqB.indexOf(txt)>0)? 1 : 0;
				}
				var score = nbBinaire/(nbTotalTuples);
				return score;
			}
			
			var l = tabSeq.length;
			var tabDistances = [];
			
			//parcours de la matrice
			for (var i=0; i<l ;i++) {
				
				if(!tabDistances[i]) tabDistances[i] = [];
				tabDistances[i][i]=0;
				var min = 1;
				for (var j=i+1 ; j<l ; j++) {
					if(!tabDistances[j]) tabDistances[j] = [];
					var distance = 1 - nbTuples( tabSeq[i].seqCompresse , tabSeq [j].seqCompresse );
					tabDistances[i][j] = distance;
					tabDistances[j][i] = distance;
				}
			}
			return tabDistances;
		};
	
	/*
	 * Construction des clusters à partir de la matrice des distances
	 */
	u.arbre = function(mD,tSeq) {
		
		var clusters = [];
		var nbSeq = mD.length;
		var tabIndexes = []; //(function(nb){ var t=[]; for (var i=0;i<nb;i++) t[i]=i;return t})(nbSeq);
		
		for (var i=0 ; i<nbSeq ; i++) {
			clusters[i] = {
				'seq' : tSeq[i],
				'profile' : [],
				'childA' : i,
				'childB' : i,
				'distance' : 0,
				'numSeq' : [i] 
			};
			tabIndexes[i] = i;
		}
		
		//parcours de la matrice pour retrouver le minimum
		function minMatriceDistances(matrice) {
			var min = 1 , minX = 0 , minY = 0, l = matrice.length;
			
			for (var i=0; i<l ; i++) {
				for (var j=i+1; j<l; j++) { //parcours de la matrice triangulaire uniquement
					//if (matrice[i][j] ==0 ) continue
					if (matrice[i][j] < min) {
						min = matrice[i][j];
						minX = i;
						minY = j;
					}
				}
			}
			//var a = [min , tabIndexes[minX], tabIndexes[minY], minX, minY]
			var a = {
				'seq' : {'seq' : ""},
				'profile' : [],
				'childA' : tabIndexes[minX],
				'childB' : tabIndexes[minY],
				'distance' : min,
				'previousA' : minX,
				'previousB' : minY,
				'msa' : [],
				'numSeq' : []
			};
			return a;
		}
		
		function recalculMatriceDistances( matrice, x, y, tI) {
			var xyMax = Math.max(x,y);
			var xyMin = Math.min(x,y);
			var vecteurMoyen = [];
			var l = matrice.length;
			
			for (var i=0; i<l; i++) { //on parcourt les lignes
				if ((i==x)||(i==y))
					continue;
				var moy = 0.1*(matrice[i][x] + matrice[i][y])/2 + 0.9*Math.min(matrice[i][x] , matrice[i][y]) ; //d'après MAFFT 
				matrice[i].push(moy);
				
				vecteurMoyen.push(moy);
				
				var dummy = matrice[i].splice(xyMax,1);
				dummy = matrice[i].splice(xyMin,1);
			}
			vecteurMoyen.push(0);
			matrice[matrice.length] = vecteurMoyen;
			//console.log(vecteurMoyen)
			dummy = matrice.splice(xyMax,1);
			dummy = matrice.splice(xyMin,1);
			dummy = tI.splice(xyMax,1);
			dummy = tI.splice(xyMin,1);
			
			return matrice;
		};
		
		for (var i=0, nbIterations = nbSeq-1; i<nbIterations; i++) {
			//on ajoute aux clusters les valeurs des noeuds correspondant aux distances les plus faibles entre noeuds et/ou feuilles
			var noeud = minMatriceDistances ( mD );
			tabIndexes.push(nbSeq+i);
			clusters.push( noeud );
			//console.log(i+":",clusters.slice(),"indexes:",tabIndexes);
			//on recalcule la matrice
			mD = recalculMatriceDistances( mD , noeud['previousA'] , noeud['previousB'] , tabIndexes);
			//console.log("matrice distances",mD2.slice())
		}
		
		return clusters;
	};
		
	u.substitutionNum = function(i,j) {
			return params.matrix[i][j];
	};
	
	u.sumOfPairsScorePP = function(profA, profB) {
		var score = 0, resProfNo = 0;
		
		for (var i=0; i< profB.m_uResidueGroup; i++) {
			resProfNo = profB.m_uSortOrder[i];
			score += profB.m_fcCounts[resProfNo] * profA.m_AAScores[resProfNo];
		}
		score /= profB.nbSeq;
		//if (score<0) console.log(resNo,prof.m_uSortOrder, score)
		return score;
	};
	
	u.sumOfPairsScoreSP = function(resA, prof) {
		var score = 0, resProfNo = 0;
		
		score = prof.m_AAScores[resA];
		return score;
	};

	/*
	 * Définition de l'objet ProfPos()
	 */
	u.ProfPos = function() {
		this.m_bAllGaps = false;
		this.m_uSortOrder = []; //acides aminés triés par ordre d'abondance
		this.m_fcCounts = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]; //effectif de chaque aa dans le profil à cette position
		this.m_LL = 1;
		this.m_LG = 0;
		this.m_GL = 0;
		this.m_GG = 0;
		this.m_AAScores = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]; //scores pour chaque acide aminé
		this.m_uResidueGroup = 0; //nb de résidus différents
		this.m_fOcc = 1; //fréquence d'occupation de la position (non indels)
		this.m_fcStartOcc = 0; //fréquence d'occupation par le début d'un gap
		this.m_fcEndOcc = 0;
		this.m_scoreGapOpen =0; //score pour l'ouverture d'un gap à cette position
		this.m_scoreGapClose = 0; //score pour la fermeture d'un gap à cette position
		this.nbSeq = 0; //nb de séquences dans le profil
	};
	
	/**calcule le profil à partir d'un alignement
	 * @msa : alignement reçu sous la forme d'un tableau contenant chaque séquence
	 * @weight : tableau indiquant le "poids" de chaque séquence
	 */
	u.profileFromMSA = function(msa, gapO, gapE) {
		//renvoie un tableau à 2D contenant pour chaque position du profil, le nombre d'acides aminés, 
		//le nombre de indels dans un gap en ouverture, en fermeture ou en extension
		var longueur = msa[0].length, //taille totale de l'alignement (nb cols)
			aa="", 					// acide aminé en code à 1 lettre
			na=0;					// acide aminé en n° (correspondant aux tables)
		var tab = Array(longueur);	//tableau contenant le profil
		//var weight = (weight || (function(l){var t=[], q=1/l; for (var i=0;i<l;i++) {t.push(q)} return t})(msa.length)  ); //faire un tableau
									//tableau contenant la pondération de chaque séquence
		var uHydrophobicRunLength = 0, //calcul de la taille d'une fenêtre hydrophobe (si elle existe) pour éviter d'insérer des gaps au sein de celle-ci
			w=1/msa.length;						//pondération d'une séquence
		//console.log( weight );
		
		var gapO = gapO || params.gapO;
		var gapE = gapE || params.gapE;
		
		for (var col=0; col<longueur; col++) {
			tab[col] = new u.ProfPos();
			var fGap =0;
			
			for (var numSeq=0, max= msa.length; numSeq<max; numSeq++) {
				tab[col].nbSeq = max;
				//w = weight[numSeq];
				aa = msa[numSeq][col];
				if (aa=="_") {
					fGap += w;
					if (col==0) {
						tab[col].m_fcStartOcc += w;
					} else if(msa[numSeq][col-1]!="_") {
						tab[col].m_fcStartOcc += w;
					}
					if (col==longueur-1) {
						tab[col].m_fcEndOcc += w;
					} else if (msa[numSeq][col+1]!="_") {
						tab[col].m_fcEndOcc += w;
					}
				} else {
					na = this.codeAANum[aa];
					tab[col].m_fcCounts[na] ++;
					if (tab[col].m_fcCounts[na] == 1) {
						tab[col].m_uSortOrder.push(na);
						tab[col].m_uResidueGroup ++;
					}
				}	
			}
			tab[col].m_fOcc = 1-fGap;
			tab[col].m_scoreGapOpen = gapO * (1 - tab[col].m_fcStartOcc);
			tab[col].m_scoreGapExtend = gapE * tab[col].m_fOcc;
			tab[col].m_scoreGapClose = gapO * (1 - tab[col].m_fcEndOcc);
			
			for (var aaNum=0; aaNum<20; aaNum++) {
				for (var i=0; i< tab[col].m_uResidueGroup; i++) {
					var resProfNo = tab[col].m_uSortOrder[i];
					tab[col].m_AAScores[aaNum] += tab[col].m_fcCounts[resProfNo] * params.matrix[aaNum][resProfNo] / msa.length;
				}
			}
		}
		tab[0].m_scoreGapOpen = 0;
		tab[longueur-1].m_scoreGapClose = 0;
		
		return tab;
	};
	
})(MultAlign.utils , MultAlign.params);

/*
 * Paramètres de l'alignement
 */
(function(p) {
	
	/*
	 * type de séquence (adn/arn/protein)
	 */
	p.typeSeq = "";
	
	/*
	 * matrice de distance (BLOSUM62, BLASTZ pour adn) et pénalités
	 */
	p.matrix = [];
	p.gapOP = 0;
	p.gapEP = 0;
	
	/*
	 * différentes matrices
	 */
	p.BLOSUM62 = {
		matrix : [
		//    A   C   D   E   F   G   H   I   K   L   M   N   P   Q   R   S   T   V   W   Y 
		    [ 4,  0, -2, -1, -2,  0, -2, -1, -1, -1, -1, -2, -1, -1, -1,  1,  0,  0, -3, -2, 0],  // A
		    [ 0,  9, -3, -4, -2, -3, -3, -1, -3, -1, -1, -3, -3, -3, -3, -1, -1, -1, -2, -2, 0],  // C
		    [-2, -3,  6,  2, -3, -1, -1, -3, -1, -4, -3,  1, -1,  0, -2,  0, -1, -3, -4, -3, 0],  // D
		    [-1, -4,  2,  5, -3, -2,  0, -3,  1, -3, -2,  0, -1,  2,  0,  0, -1, -2, -3, -2, 0],  // E
		    [-2, -2, -3, -3,  6, -3, -1,  0, -3,  0,  0, -3, -4, -3, -3, -2, -2, -1,  1,  3, 0],  // F
		    [ 0, -3, -1, -2, -3,  6, -2, -4, -2, -4, -3,  0, -2, -2, -2,  0, -2, -3, -2, -3, 0],  // G
		    [-2, -3, -1,  0, -1, -2,  8, -3, -1, -3, -2,  1, -2,  0,  0, -1, -2, -3, -2,  2, 0],  // H
		    [-1, -1, -3, -3,  0, -4, -3,  4, -3,  2,  1, -3, -3, -3, -3, -2, -1,  3, -3, -1, 0],  // I
		    [-1, -3, -1,  1, -3, -2, -1, -3,  5, -2, -1,  0, -1,  1,  2,  0, -1, -2, -3, -2, 0],  // K
		    [-1, -1, -4, -3,  0, -4, -3,  2, -2,  4,  2, -3, -3, -2, -2, -2, -1,  1, -2, -1, 0],  // L
		    [-1, -1, -3, -2,  0, -3, -2,  1, -1,  2,  5, -2, -2,  0, -1, -1, -1,  1, -1, -1, 0],  // M
		    [-2, -3,  1,  0, -3,  0,  1, -3,  0, -3, -2,  6, -2,  0,  0,  1,  0, -3, -4, -2, 0],  // N
		    [-1, -3, -1, -1, -4, -2, -2, -3, -1, -3, -2, -2,  7, -1, -2, -1, -1, -2, -4, -3, 0],  // P
		    [-1, -3,  0,  2, -3, -2,  0, -3,  1, -2,  0,  0, -1,  5,  1,  0, -1, -2, -2, -1, 0],  // Q
		    [-1, -3, -2,  0, -3, -2,  0, -3,  2, -2, -1,  0, -2,  1,  5, -1, -1, -3, -3, -2, 0],  // R
		    [ 1, -1,  0,  0, -2,  0, -1, -2,  0, -2, -1,  1, -1,  0, -1,  4,  1, -2, -3, -2, 0],  // S
		    [ 0, -1, -1, -1, -2, -2, -2, -1, -1, -1, -1,  0, -1, -1, -1,  1,  5,  0, -2, -2, 0],  // T
		    [ 0, -1, -3, -2, -1, -3, -3,  3, -2,  1,  1, -3, -2, -2, -3, -2,  0,  4, -3, -1, 0],  // V
		    [-3, -2, -4, -3,  1, -2, -2, -3, -3, -2, -1, -4, -4, -2, -3, -3, -2, -3, 11,  2, 0],  // W
		    [-2, -2, -3, -2,  3, -3,  2, -1, -2, -1, -1, -2, -3, -1, -2, -2, -2, -1,  2,  7, 0],  // Y
		    [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0]   // *
		],
		gapOP : -7,
		gapEP : -1,
		gapOPMSA : -15,
		gapEPMSA : -1
   };
   
   p.BLASTZ = {
   		matrix : [
   			//	A		C		G		T
   			[  91,	 -114,	  -31,	 -123],
   			[-114,	  100,	 -125,	  -31],
   			[ -31,	 -125,	  100,	 -114],
   			[-123,	  -31,	 -114,	   91]
   		],
   		gapOP : -400,
   		gapEP : -30,
   		gapOPMSA : -400,
   		gapEPMSA : -1
   };
   
   /*
    * Fonction d'initialisation des paramètres d'alignement
    */
   p.init = function(matrice,isMSA) {
   		var m,
   			isMSA = isMSA || false;
   		
   		switch (p.typeSeq) {
   			case "protein":
   				m = p.BLOSUM62;
   				p.matrix = m.matrix;
   				p.gapOP = (isMSA)? m.gapOPMSA : m.gapOP;
   				p.gapEP = (isMSA)? m.gapEPMSA : m.gapEP;
   			break;
   			default:
   				m = p.BLASTZ;
   				p.matrix = m.matrix;
   				p.gapOP = (isMSA)? m.gapOPMSA : m.gapOP;
   				p.gapEP = (isMSA)? m.gapEPMSA : m.gapEP;
   		}
   		/*if (matrice=="BLOSUM62") {
   			m = p.BLOSUM62;
   			p.matrix = m.matrix;
   			p.gapOP = (isMSA)? m.gapOPMSA : m.gapOP;
   			p.gapEP = (isMSA)? m.gapEPMSA : m.gapEP;
   		}*/
   };
   
   /*
    * Taille de la fenêtre de calcul du dotplot
    */
   p.largeurFenetre = 15;
   
   p.setSeqType = function(type) {
   		p.typeSeq = type;
   		
   		
   		
   };
})(MultAlign.params);