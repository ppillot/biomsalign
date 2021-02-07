/*global jQuery, MultAlign*/
/**
 * @author Paul Pillot
 */

var SeqView = (function () {
    var s = {};

    s.sequences = {};
    s.pannels = {};
    s.utils = {};
    s.svg = {};

    return s;

})(document);

/**
 * ajout d'un nouveau panneau
 */
(function (s, $) {
    s._Pannel = function (infos) {
        //tailles des objets à afficher
        this.nbSeq = 0; //nb séquences = taille y
        this.longueurMax = 0; //nb max residus à 1 lettre = taille x
        this.longueurMax3 = 0; //nb max residus à 3 lettres = taille x

        this.id = infos.conteneur; //conteneur dans le DOM
        this.selector = $('#' + this.id); //selecteur jQuery correspondant
        this.css = new s.utils.StyleSheet;

        this.sequences = infos.sequences || []; //affectation des séquences transmises au constructeur
        this.sequencesSelected = [];

        this.offsetX = -1; //pour pouvoir initialiser l'affichage
        this.offsetY = -1;

        this.typeRegle3 = (this.sequences[0] == undefined) ? true : (this.sequences[0].type == 'protein') ? true : false; //affichage des graduations par 3 ou 1 en fonction du type de la première séquence

        s.pannels[this.id] = this;
        /*
         * ajoute au conteneur :
         * - un panneau pour les titres (défilement vertical)
         * - un panneau pour les graduations (défilement horizontal)
         * - un panneau pour les séquences (défilement horizontal et vertical)
         */
        var $conteneur = $('<div class="view-seq-conteneur" style=""></div>').appendTo(this.selector),
            $gauche = $('<div style="position: absolute; top:0; left:0; bottom:0; overflow: hidden; width:10em;"></div>').appendTo($conteneur),
            $droite = $('<div style="position: absolute; top:0; left:10em; bottom:0; right:0; overflow: hidden;"></div>').appendTo($conteneur);

        this.$contTitres = $('<div class="view-seq-conteneur-sequences"></div>').appendTo($gauche);
        this.$titresSequences = $('<div class="view-seq-titres-sequences"></div>').appendTo(this.$contTitres)
            .scroll(this, this.scrollSequences);
        this.$regle = $('<div class="view-seq-regle"></div>').appendTo($droite)
            .on('click', this, this.clickRegle);

        this.$contSequences = $('<div class="view-seq-conteneur-sequences"></div>').appendTo($droite);
        this.$sequences = $('<div class="view-seq-sequences"></div>').appendTo(this.$contSequences);

        this.$contDummy = $('<div class="view-seq-conteneur-defilement"></div>').appendTo($droite)
            .scroll(this, this.scrollSequences)
            .on('click', this, this.clickSequences);
        this.$dummy = $('<div style="position: relative; top:0; left:0; "></div>').appendTo(this.$contDummy);

        this.init();
        //gestion du highlight
        this.highlight.init(this);

        //gestion de la coloration
        this.coloration.init(this);

        this.displayOffset(0, 0);



    };

    /*
     * méthode appelée :
     * - à chaque redimensionnement
     * - après chaque ajout ou suppression de séquence
     */
    s._Pannel.prototype.init = function () {
        /**
         * Détermine le nb de séquences et leur longueur max
         */
        this.getSequencesDimensions();

        /**
         * ajoute la règle, au panneau sequences les sequences, et au panneau titre, les titres
         */

        //on prend les mesures des résidus dans le panneau de séquences
        s.utils.getCharSeqDimensions(this);

        //on prend la mesure de la fenêtre de visualisation
        this.nbSeqAffichables = Math.floor(this.$contSequences.height() / s.utils.charH);
        this.nbResAffichables = Math.floor(this.$contSequences.width() / s.utils.charW);
        this.nbRes3Affichables = Math.floor(this.$contSequences.width() / s.utils.charW3);

        //s'il y a trop de résidus à afficher, il y aura une barre de défilement en bas !
        if ((this.nbResAffichables < this.longueurMax) || (this.nbRes3Affichables < this.longueurMax3)) {
            this.nbSeqAffichables = Math.floor((this.$contSequences.height() - s.utils.getScrollBarWidth()) / s.utils.charH);
        }
        //s'il y a trop de séquences à afficher, il y aura une barre de défilement sur le côté !
        if (this.nbSeqAffichables < this.nbSeq) {
            //this.nbResAffichables = Math.floor(($contSequences.width()-s.utils.getScrollBarWidth())/s.utils.charW);
            //this.nbRes3Affichables = Math.floor(($contSequences.width()-s.utils.getScrollBarWidth())/s.utils.charW3);

        }

        var largeurPanneauSequences = Math.max(this.longueurMax * s.utils.charW, this.longueurMax3 * s.utils.charW3) + s.utils.scrollBarSize() * 3;
        var hauteurPanneauSequences = this.nbSeq * s.utils.charH + s.utils.scrollBarSize();

        /**
         * Création de la règle et de ses graduations
         * trois objets svg sont créés :
         * - les tirets pour chaque résidu code 1 lettre
         * - les tirets pour chaque résidu code 3 lettres
         * - les tirets et textes pour les mutiples de 5
         */
        this.hauteurRegle = this.$regle.height();
        var svgRegle = s.svg.createSVG(largeurPanneauSequences, this.hauteurRegle);
        this.$regle.empty()[0].appendChild(svgRegle);

        //calcul des graduations
        this.marquesGraduations1 = s.svg.group();
        svgRegle.appendChild(this.marquesGraduations1);

        for (var i = 0; i < this.nbResAffichables; i++) {
            var x = (i + 0.5) * s.utils.charW;
            var line = s.svg.line(x, this.hauteurRegle, x, Math.round(this.hauteurRegle * 2 / 3), '#444', '2'); //'<line y1="1.9em" x1="'+ (i+0.5)*2.2 +'em" y2="1.5em" x2="'+ (i+0.5)*2.2 +'em" stroke="#444" stroke-width="1px" />'	;
            this.marquesGraduations1.appendChild(line);
        }

        this.marquesGraduations3 = s.svg.group();
        svgRegle.appendChild(this.marquesGraduations3);

        for (i = 0; i < this.nbRes3Affichables; i++) {
            x = (i + 0.5) * s.utils.charW3;
            line = s.svg.line(x, this.hauteurRegle, x, Math.round(this.hauteurRegle * 2 / 3), '#444', '2'); //'<line y1="1.9em" x1="'+ (i+0.5)*2.2 +'em" y2="1.5em" x2="'+ (i+0.5)*2.2 +'em" stroke="#444" stroke-width="1px" />'	;
            this.marquesGraduations3.appendChild(line);
        }

        if (this.typeRegle3) {
            this.marquesGraduations1.style.visibility = 'hidden';
        } else {
            this.marquesGraduations3.style.visibility = 'hidden';
        }

        this.graduations = s.svg.group();
        svgRegle.appendChild(this.graduations);

        this.$dummy.width(largeurPanneauSequences)
            .height(hauteurPanneauSequences);
        this.$contDummy.height(Math.min(this.nbSeqAffichables * s.utils.charH + s.utils.scrollBarSize(), this.nbSeq * s.utils.charH + s.utils.scrollBarSize()))
            .width(Math.max(this.nbResAffichables * s.utils.charW, this.nbRes3Affichables * s.utils.charW3) + s.utils.scrollBarSize());

        /**
         * on vérifie que le toutes les séquences sont affichables ensemble, sinon on prévoit
         * la dimension de la barre de défilement pour le panneau des titres de séquences
         */
        /*if (this.$titresSequences[0].offsetHeight>this.$titresSequences[0].parentNode.offsetHeight) {
        	$contTitres.css('bottom',s.utils.scrollBarSize()+'px');
        }*/
        if (this.nbSeq > this.nbSeqAffichables) {
            this.$contTitres.css('bottom', s.utils.scrollBarSize() + 'px');
        }
    };

    s._Pannel.prototype.displayOffset = function (numResDeb, numSeqDeb) {
        numResDeb = (numResDeb > -1) ? numResDeb : -1;
        numSeqDeb = (numSeqDeb > -1) ? numSeqDeb : -1;

        var noScrollX = (numResDeb == this.offsetX);
        var noScrollY = (numSeqDeb == this.offsetY);

        if (noScrollX && noScrollY) {
            return;
        }
        if (numResDeb == -1) {
            numResDeb = this.offsetX;
            numSeqDeb = this.offsetY;
        }

        //On efface titres et séquences
        this.$sequences.empty();
        if (!noScrollY) this.$titresSequences.empty();
        var seq = '',
            residuCss = '',
            numResFin = 0;

        for (var i = numSeqDeb, imax = Math.min(numSeqDeb + this.nbSeqAffichables, this.nbSeq); i < imax; i++) {
            //console.log(this.nbSeq, imax)
            if (this.sequences[i].type == 'protein') {
                seq = this.sequences[i].seq3;
                residuCss = 'view-seq-residu3';
                var numResDebA = numResDeb;
                numResFin = Math.min(numResDebA + this.nbRes3Affichables, seq.length);
            } else {
                seq = this.sequences[i].seq;
                residuCss = 'view-seq-residu';
                numResDebA = 3 * numResDeb;
                numResFin = Math.min(numResDebA + this.nbResAffichables, seq.length);
            }
            var texteSeq = '';

            for (var j = numResDebA; j < numResFin; j++) {

                texteSeq += '<div class=\'' + residuCss + ' view-seq-residu-' + this.sequences[i].seq[j] + ' ' + residuCss + '-colonne' + j + '\' id=\'s' + i + 'r' + j + '\'>' + seq[j] + '</div>';
            }
            if (texteSeq == '') texteSeq = '<div class=\'view-seq-residu\'>&nbsp;</div>';

            if (!noScrollY) {
                var $titre = $('<div draggable=\'true\'>' + this.sequences[i].titre + '</div>').appendTo(this.$titresSequences);
                $('<input type=\'checkbox\' ' + ((this.sequences[i].checked) ? 'checked' : '') + '/>').prependTo($titre)
                    .on('click', this, this.clickTitreSequence)
                    .attr('id', 'view-seq-' + this.id + '-' + i);
            }
            $('<div class=\'view-seq-ligne-sequences\'>' + texteSeq + '</div>').appendTo(this.$sequences);

        }

        if (!noScrollX) {
            this.afficheGraduations(numResDeb);
        }

        /*if (this.highlight.displayed)
        		this.highlight.show(numResDeb,numResDeb+this.nbRes3Affichables);
        */
        this.offsetX = numResDeb;
        this.offsetY = numSeqDeb;
    };


    s._Pannel.prototype.afficheGraduations = function (numResDeb) {
        var numResDebA, imax, largeur;
        //calcul des graduations des x10 et des x5
        //on vide les graduations
        while (this.graduations.firstChild) {
            this.graduations.removeChild(this.graduations.firstChild);
        }

        //position du début
        if (this.typeRegle3) {
            numResDebA = numResDeb;
            imax = this.nbRes3Affichables;
            largeur = s.utils.charW3;
        } else {
            numResDebA = numResDeb * 3;
            imax = this.nbResAffichables;
            largeur = s.utils.charW;
        }

        //quelle est le decallage de la première graduation multiple de 5 ?
        var offsetGrad5 = 5 - ((numResDebA + 1) % 5);
        if (offsetGrad5 == 5) offsetGrad5 = 0;

        //Affichage des graduations
        for (var i = offsetGrad5; i < imax; i += 5) {
            var y2 = ((numResDebA + i + 1) % 10 == 0) ? Math.round(this.hauteurRegle / 2) : Math.round(this.hauteurRegle / 1.5);
            var x = Math.round((i + 0.5) * largeur * 10) / 10;
            var line = s.svg.line(x, this.hauteurRegle, x, y2, '#444', '3');
            this.graduations.appendChild(line);

            var etiquette = s.svg.text(x, Math.round(this.hauteurRegle / 2) - 2, numResDebA + i + 1, '#444', 'middle');
            this.graduations.appendChild(etiquette);
        }
    };

    s._Pannel.prototype.getSequencesDimensions = function () {
        this.nbSeq = this.sequences.length;
        var max = 0,
            max3 = 0;
        for (var i = 0; i < this.nbSeq; i++) {
            if (this.sequences[i].type == 'protein') {
                max3 = Math.max(max3, this.sequences[i].seq.length);
            } else {
                max = Math.max(max, this.sequences[i].seq.length);
            }
        }
        this.longueurMax = max;
        this.longueurMax3 = max3;
        //console.log(this.nbSeq,this.longueurMax)
    };


    s._Pannel.prototype.scrollSequences = function (event) {
        event.stopPropagation();
        var x = event.target.scrollLeft;
        var y = event.target.scrollTop;
        var that = event.data;

        var offsetX = Math.floor(x / s.utils.charW3);
        var offsetY = Math.floor(y / s.utils.charH);

        that.displayOffset(offsetX, offsetY);

    };

    s._Pannel.prototype.clickSequences = function (event) {
        event.stopPropagation();
        var that = event.data;
        var posX = event.pageX - that.$sequences.offset().left;
        var posY = event.pageY - that.$sequences.offset().top;

        //quel est le n° de séquence ?
        var numSeq = Math.floor(posY / s.utils.charH) + that.offsetY;

        //quel est le n° de résidu ?
        var type = that.sequences[numSeq].type;
        var largeur = (type == 'protein') ? s.utils.charW3 : s.utils.charW;
        var numRes = Math.floor(posX / largeur) + that.offsetX * ((type == 'protein') ? 1 : 3);

        //quel est le résidu ?
        that.highlight.set(numRes, (that.sequences[numSeq].type == 'protein') ? 3 : 1);

        //Mise à jour des graduations
        if (((type == 'protein') && (!that.typeRegle3)) || ((type != 'protein') && (that.typeRegle3))) {
            that.swapRegle();
        }
        //var nomRes = (type == 'protein') ? that.sequences[numSeq].seq3[numRes] : that.sequences[numSeq].seq[numRes];
        //console.log(numSeq, numRes, that.sequences[numSeq].titre, nomRes+(numRes+1));
    };

    var Highlight = function () {
        this.x = [];
        this.x3 = -1;
        this.xCodon = [];

        this.displayed = false;
        this.parent = null;
        this.selector = null;

        this.init = function (parent) {
            this.parent = parent;
        };

        this.set = function (numRes, tailleRes) {
            var x3 = 0,
                x = [],
                xCodon = [];

            if (tailleRes == 3) { //on a cliqué sur un résidu d'une protéine
                x3 = numRes;
                x = [numRes * 3, numRes * 3 + 1, numRes * 3 + 2];
            } else {
                x = [numRes];
                x3 = Math.floor(numRes / 3);
                xCodon = [x3 * 3, x3 * 3 + 1, x3 * 3 + 2].filter(function (value) {
                    return value != x;
                });
            }

            if (this.x.toString() == x.toString()) { //a-t-on cliqué sur une colonne déjà mise en valeur ?
                if (this.displayed) {
                    this.displayed = false;
                    this.hide();
                } else {
                    this.displayed = true;
                    this.show();
                }
            } else {
                this.x = x;
                this.x3 = x3;
                this.xCodon = xCodon;
                this.displayed = true;
                this.show();
            }
        };

        this.show = function () {
            var panelSelector = '#' + this.parent.id;
            var ruleNucleotide = this.x.map(function (value) {
                return panelSelector + ' .view-seq-residu.view-seq-residu-colonne' + value;
            }).join(' , ');

            var rule = ruleNucleotide //"#" + this.parent.id + " .view-seq-residu-colonne" + this.x
                + ', #' + this.parent.id + ' .view-seq-residu3-colonne' + this.x3 + ' {background: #ACF !important}';
            this.parent.css.setRule('highlight', rule);

            if (this.xCodon.length > 0) {
                var ruleNucleotide2 = this.xCodon.map(function (value) {
                    return panelSelector + ' .view-seq-residu.view-seq-residu-colonne' + value;
                }).join(' , ');

                var rule2 = ruleNucleotide2 + ' {background: #CEF !important}';
                this.parent.css.setRule('highlight2', rule2);
            } else {
                this.parent.css.unsetRule('highlight2');
            }


        };

        this.hide = function () {
            this.parent.css.unsetRule('highlight');
            this.parent.css.unsetRule('highlight2');
        };
    };

    s._Pannel.prototype.highlight = new Highlight();

    var Coloration = function () {
        //propriétés
        this.listRules = [];
        this.parent = {};

        //méthodes privées
        this._addRuleName = function (num) {
            var name = 'colorationSequence' + num;
            this.listRules.push(name);
            return name;
        };

        this._rule = function (parentID, res, color) {
            var selector = '#' + parentID + ' .view-seq-residu-';
            selector += (res.toLowerCase() == res) ? res.toUpperCase() + '.view-seq-residu' : res;
            var cssTxt = ' {background: ' + color + '}';
            return selector + cssTxt;
        };

        //méthodes publiques
        this.init = function (parent) {
            this.parent = parent;
        };

        this.set = function (typeColoration) {

            switch (typeColoration) {
            case 'residu':
                var colorScheme = [
                    {
                        res: 'C',
                        color: '#ffed6f'
                    }, {
                        res: 'V',
                        color: '#d9d9d9'
                    }, {
                        res: 'I',
                        color: '#d9d9d9'
					}, {
                        res: 'L',
                        color: '#d9d9d9'
					}, {
                        res: 'A',
                        color: '#d9d9d9'
					}, {
                        res: 'G',
                        color: '#d9d9d9'
					}, {
                        res: 'M',
                        color: '#d9d9d9'
					}, {
                        res: 'T',
                        color: '#b3de69'
					}, {
                        res: 'S',
                        color: '#b3de69'
					}, {
                        res: 'D',
                        color: '#fb8072'
					}, {
                        res: 'E',
                        color: '#fb8072'
					}, {
                        res: 'H',
                        color: '#80b1d3'
					}, {
                        res: 'K',
                        color: '#80b1d3'
					}, {
                        res: 'R',
                        color: '#80b1d3'
					}, {
                        res: 'N',
                        color: '#bc80bd'
					}, {
                        res: 'Q',
                        color: '#bc80bd'
					}, {
                        res: 'F',
                        color: '#fdb462'
					}, {
                        res: 'Y',
                        color: '#fdb462'
					}, {
                        res: 'W',
                        color: '#fdb462'
					}, {
                        res: 'P',
                        color: '#ffffb3'
					}, {
                        res: 'a',
                        color: '#fb8072'
					}, {
                        res: 't',
                        color: '#80b1d3'
					}, {
                        res: 'c',
                        color: '#ffffb3'
					}, {
                        res: 'g',
                        color: '#b3de69'
					}, {
                        res: 'u',
                        color: '#bebada'
					}
				];
                for (var i = 0, imax = colorScheme.length; i < imax; i++) {
                    var txt = this._rule(this.parent.id, colorScheme[i].res, colorScheme[i].color);
                    var ruleName = this._addRuleName(i);
                    this.parent.css.setRule(ruleName, txt);
                    //console.log(txt)
                }
                break;
            case 'propriete':

                break;
            default:
                for (var i = 0, imax = this.listRules.length; i < imax; i++) {
                    this.parent.css.unsetRule(this.listRules[i]);
                }
            }
        };
    };

    s._Pannel.prototype.coloration = new Coloration();


    s._Pannel.prototype.clickRegle = function (event) {
        event.stopPropagation();
        var that = event.data;

        that.swapRegle();
    };

    s._Pannel.prototype.swapRegle = function () {
        this.typeRegle3 = !this.typeRegle3;

        if (this.typeRegle3) {
            this.marquesGraduations1.style.visibility = 'hidden';
            this.marquesGraduations3.style.visibility = 'visible';

        } else {
            this.marquesGraduations3.style.visibility = 'hidden';
            this.marquesGraduations1.style.visibility = 'visible';

        }

        this.afficheGraduations(this.offsetX);
    };

    s._Pannel.prototype.clickTitreSequence = function (event) {
        event.stopPropagation();
        var that = event.data;
        var numSeq = parseInt(event.target.id.substr(event.target.id.lastIndexOf('-') + 1));
        that.sequences[numSeq].checked = !that.sequences[numSeq].checked;
        if (that.sequences[numSeq].checked) {
            that.sequencesSelected.push(numSeq);
        } else {
            that.sequencesSelected = that.sequencesSelected.filter(function (value) {
                return value != numSeq;
            });
        }

        //console.log(that,event.target.id);
    };

    s._Pannel.prototype.ajouteSequence = function (seq) {
        if (typeof (seq) == 'string') seq = s.utils.getSequenceFromText(seq);
        if (seq.length) { //il y a plusieurs séquences à ajouter, dans un tableau
            this.sequences = this.sequences.concat(seq);
            //console.log(seq, this.sequences);

        } else { //il n'y a qu'une séquence à ajouter, à partir d'un objet séquence
            this.sequences.push(seq);
        }
        this.init();
        this.displayOffset();
    };

    s._Pannel.prototype.supprimeSelection = function () {
        var sel = this.sequencesSelected;
        this.sequences = this.sequences.filter(function (value, key) {
            return (sel.indexOf(key) == -1);
        });
        this.sequencesSelected = [];
        this.init();
        this.displayOffset();
    };

    s._Pannel.prototype.supprimeTout = function () {
        this.sequences = [];
        this.sequencesSelected = [];
        this.init();
        this.displayOffset();

        return this;
    };

    s._Pannel.prototype.traduitSelection = function () {
        var sel = this.sequencesSelected,
            //nouvTitre = '',
            nouvSeq = '';

        //vérification du type de séquences de la sélection avant de faire le traitement
        for (var i = 0, imax = sel.length; i < imax; i++) {
            if (this.sequences[sel[i]].type == 'protein') {
                //problème car on ne peut pas traduire des protéines
                return alert('La sélection contient des séquences protéiques qui ne peuvent pas être traduites. Les retirer de la sélection avant de faire le traitement');
            }
        }

        //traduction
        for (var i = 0, imax = sel.length; i < imax; i++) {
            nouvSeq = s.utils.traduit(this.sequences[sel[i]].seq);

            //ajout de la nouvelle séquence à la liste des séqunces existantes
            this.sequences.push({
                titre: 'traduction ' + this.sequences[sel[i]].titre,
                seq: nouvSeq,
                type: 'protein',
                checked: false,
                seq3: s.utils.aa123(nouvSeq)
            });
        }
        this.init();
        this.displayOffset();

    };

    s._Pannel.prototype.aligneSelection = function (target) {
        var sel = this.sequencesSelected,
            typeSeq = '',
            alignement = [];

        //vérification du nb de séquences sélectionnées (>1)
        if (sel.length < 2) return alert('Au moins deux séquences doivent être sélectionnées');

        //vérification de l'unicité du type de séquences de la sélection avant de faire le traitement
        typeSeq = this.sequences[sel[0]].type;
        for (var i = 1, imax = sel.length; i < imax; i++) {
            if (this.sequences[sel[i]].type != typeSeq) {
                //problème car on ne peut pas traduire des protéines
                return alert('Toutes les séquences à aligner ne sont pas du même type. Retirer de la sélection\
				ou convertir les séquences de types différents.');
            }
        }

        //initialisation de l'objet alignement
        MultAlign.reset();

        //création du tableau de séquences à aligner
        for (i = 0; i < imax; i++) {
            MultAlign.ajouteSequences(this.sequences[sel[i]].seq);
        }

        //Alignement
        alignement = MultAlign.aligne();

        //console.log(alignement, this.sequences);
        //création de l'objet séquences
        var sequences = [];
        for (var i = 0, imax = alignement.length; i < imax; i++) {
            sequences[i] = {
                titre: this.sequences[sel[i]].titre,
                seq: alignement[i],
                type: typeSeq,
                checked: false,
                seq3: s.utils.aa123(alignement[i])
            };
        }

        //le panneau existe-t-il ?
        if (typeof (s.pannels[target]) == 'object') {
            //oui ! il existe déjà ! on ajoute les séquences au panneau
            s.pannels[target].ajouteSequence(sequences);
        } else {
            //non ! il n'existe pas ! on crée un nouvel objet SeqView que l'on renvoie
            return s.getPannel({
                conteneur: target,
                sequences: sequences
            });
        }


    };

})(SeqView, jQuery, document);


(function (s, $) {
    s.getPannel = function (infos) {
        var lInfos = {
            conteneur: infos.conteneur,
            sequences: []
        };

        //des vérifications à faire ici
        if (infos.sequences) {
            var l = infos.sequences.length;
            for (var i = 0; i < l; i++) {
                lInfos.sequences[i] = {};

                var sequence = infos.sequences[i];
                if (!sequence.type) {
                    lInfos.sequences[i].type = s.utils.getSeqType(sequence.seq);
                } else {
                    lInfos.sequences[i].type = sequence.type;
                }
                if (lInfos.sequences[i].type == 'protein') {
                    lInfos.sequences[i].seq3 = s.utils.aa123(sequence.seq);
                }
                lInfos.sequences[i].seq = sequence.seq;
                lInfos.sequences[i].checked = false;
                lInfos.sequences[i].titre = ((sequence.titre) ? sequence.titre : 'Séquence ' + i);
            }
        }

        var pannel = new s._Pannel(lInfos);
        return pannel;

    };


})(SeqView, jQuery);

(function (su, $) {
    su.getSeqType = function (seq) {
        var acidesAmines = /[rndbeqyhilkmfpswyv]/gi;
        var notArn = /[^acgu]/gi;
        var notAdn = /[^acgt]/gi;

        if (acidesAmines.test(seq)) return 'protein';
        else {
            if (!notArn.test(seq)) return 'arn';
            else {
                if (!notAdn.test(seq)) return 'adn';
                else return 'error';
            }
        }
    };

    su.getSequenceFromText = function (txt) {
        var seq = [],
            titre = '',
            txtSeq = '',
            typeSeq = '',
            txtSeq3 = [],
            tabTxtSeq = [];
        //tabTxtTitres = [];

        //on recherche le format de la séquence

        if (txt.indexOf('>') == 0) {
            //FASTA ? on ne fait pas de vérification supplémentaire à ce stade !!
            var motif = />(.+)\n/gmi;
            tabTxtSeq = txt.split(motif).slice(1); //tableau où ligne n = titre, ligne n+1 = seq

            for (var i = 0, imax = tabTxtSeq.length; i < imax; i += 2) {

                titre = tabTxtSeq[i];
                if (i + 1 > imax) {
                    alert('erreur dans le texte de la séquence');
                    //console.log(i, seq, tabTxtSeq);
                    break;
                }
                txtSeq = tabTxtSeq[i + 1].split('\n').join(''); //on supprime les retours à la ligne
                typeSeq = this.getSeqType(txtSeq);
                txtSeq3 = (typeSeq == 'protein') ? this.aa123(txtSeq) : [];

                seq.push({
                    titre: titre,
                    seq: txtSeq,
                    type: typeSeq,
                    checked: false,
                    seq3: txtSeq3
                });
            }
        }

        return seq;
    };

    su.codesAa = {
        'A': 'Ala',
        'C': 'Cys',
        'R': 'Arg',
        'N': 'Asn',
        'D': 'Asp',
        'E': 'Glu',
        'Q': 'Gln',
        'G': 'Gly',
        'H': 'His',
        'I': 'Ile',
        'L': 'Leu',
        'K': 'Lys',
        'M': 'Met',
        'F': 'Phe',
        'P': 'Pro',
        'S': 'Ser',
        'T': 'Thr',
        'W': 'Trp',
        'Y': 'Tyr',
        'V': 'Val'
    };

    /**
     * Retourne un tableau contenant la séquence protéique en code à 3 lettres
     * si la lettre n'existe pas en code 3 lettres, elle est retournée telle que (permet de gérer les - !)
     */
    su.aa123 = function (seq) {
        var tabSeq = [],
            l = seq.length;

        for (var i = 0; i < l; i++) {
            tabSeq[i] = (su.codesAa[seq[i]] === undefined) ? seq[i] : su.codesAa[seq[i]];
        }

        return tabSeq;
    };

    su.getScrollBarWidth = function () {
        var inner = document.createElement('p');
        inner.style.width = '100%';
        inner.style.height = '200px';

        var outer = document.createElement('div');
        outer.style.position = 'absolute';
        outer.style.top = '0px';
        outer.style.left = '0px';
        outer.style.visibility = 'hidden';
        outer.style.width = '200px';
        outer.style.height = '150px';
        outer.style.overflow = 'hidden';
        outer.appendChild(inner);

        document.body.appendChild(outer);
        var w1 = inner.offsetWidth;
        outer.style.overflow = 'scroll';
        var w2 = inner.offsetWidth;
        if (w1 == w2) w2 = outer.clientWidth;

        document.body.removeChild(outer);

        return (w1 - w2);
    };

    su.charW = 0;
    su.charW3 = 0;
    su.charH = 0;

    su.getCharSeqDimensions = function (that) {
        var $testMesureSequences = $('<div class="view-seq-residu">a</div><div class="view-seq-residu3">aaa</div>').appendTo(that.$sequences);
        su.charW = $testMesureSequences[0].offsetWidth;
        su.charW3 = $testMesureSequences[1].offsetWidth;

        var $testMesureSequences2 = $('<div class="view-seq-ligne-sequences"><div class="view-seq-residu">a</div></div>').appendTo(that.$sequences);
        su.charH = $testMesureSequences2[0].offsetHeight;

        that.$sequences.empty();
    };

    su.scrollBarW = 0;

    su.scrollBarSize = function () {
        if (su.scrollBarW == 0) {
            su.scrollBarW = su.getScrollBarWidth();
        }
        return su.scrollBarW;
    };

    /**
     * Objet feuille de styles gérant l'ajout ou la suppression des styles
     * Cet objet peut être créé dans un objet panneau de séquences pour gérer
     * la surbrillance, les colorations par résidus, etc...
     */
    su.StyleSheet = function () {
        this.num = -1;
        this.nbRules = 0;
        this.rules = [];
        this.styleSheet = this.create();
    };

    su.StyleSheet.prototype.create = function () {
        //•• retrouver la version pour internet explorer
        var style = document.createElement('style');
        document.getElementsByTagName('head')[0].appendChild(style);

        this.num = document.styleSheets.length - 1;
        return document.styleSheets[this.num];
    };

    su.StyleSheet.prototype.setRule = function (ruleName, ruleText) {
        var ruleIndex = this.rules.indexOf(ruleName);
        if (ruleIndex == -1) { //la règle css n'existe pas encore dans la feuille de style
            ruleIndex = this.nbRules;
            this.nbRules++;
        } else { //la règle existe, il faut la supprimer
            this.styleSheet.deleteRule(ruleIndex);
        }
        var nb = this.styleSheet.insertRule(ruleText, ruleIndex);
        this.rules[nb] = ruleName;
    };

    su.StyleSheet.prototype.unsetRule = function (ruleName) {
        var ruleIndex = (typeof (ruleName) == 'number') ? ruleName : this.rules.indexOf(ruleName);
        if (ruleIndex > -1) {
            this.styleSheet.deleteRule(ruleIndex);
            this.styleSheet.insertRule('#dummy {}', ruleIndex);
        }
    };

    su.codesGenetiques = [
        'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG', //#1 : standard
        'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG', //#2: mitochondries vertébrés
        'FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG', //#3 : mitochondries levures
        'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG', //#4 : mitochondries protozoaires
        'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG', //#5 : mitochondries invertébrés
        'FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG', //#6 : ciliés, dasycladacés
        'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG' //#9 : mitochondries échinodermes et vers plats
    ];

    su.retourneCodons = function (seq) {
        var tabCodons = [];

        var chnSeq = (typeof (seq) == 'object') ? seq.join('') : seq;
        tabCodons = chnSeq.match(/(...)/g);

        return tabCodons;
    };

    su.codeGenetique = su.codesGenetiques[0];

    su.setCodeGenetique = function (nb) {
        if (this.codesGenetiques[nb])
            this.codeGenetique = this.codesGenetiques[nb];
        else return;
    };

    su.codeBase4 = function (seq) {
        var seq4 = '';
        var conv = {
            'T': '0',
            'U': '0',
            'G': '3',
            'C': '1',
            'A': '2',
            't': '0',
            'u': '0',
            'g': '3',
            'c': '1',
            'a': '2'
        };
        for (var i = 0, imax = seq.length; i < imax; i++) {
            seq4 += conv[seq[i]];
        }
        return seq4;
    };

    su.complementaire = function (seq) {
        var comp = {
                'T': 'A',
                'A': 'T',
                'G': 'C',
                'C': 'G'
            },
            seq2 = '';
        var i = seq.length;

        while (i--) {
            seq2 = comp[seq[i]] + seq2;
        }
        return seq2;
    };

    su.transcriptionADNNonCodant = function (seq) {
        var comp = {
                'T': 'A',
                'A': 'U',
                'G': 'C',
                'C': 'G'
            },
            seq2 = '';
        var i = seq.length;

        while (i--) {
            seq2 = comp[seq[i]] + seq2;
        }
        return seq2;
    };

    su.transcriptionADNCodant = function (seq) {
        var reg = /T/g;
        return seq.replace(reg, 'U');
    };

    su.traduit = function (seq) {
        //conversion de la sequence en code base 4
        console.time('traduction');
        /*seq = seq.replace(/[tu]/gi,"0")
        			.replace(/c/gi,"1")
        			.replace(/a/gi,"2")
        			.replace(/g/gi,"3");
        			*/
        seq = this.codeBase4(seq);

        var tabCodons = this.retourneCodons(seq);

        var seqTrad = '',
            aa = '';

        for (var i = 0, imax = tabCodons.length; i < imax; i++) {
            aa = this.codeGenetique[parseInt(tabCodons[i], 4)];
            if (aa == '*') break;
            else seqTrad += aa;
        }

        console.timeEnd('traduction');
        return seqTrad;
    };

    //s.utils.scrollBarSize = s.utils.getScrollBarWidth ();
})(SeqView.utils, jQuery);

/*
 * Fonctions gérant le SVG
 */
(function (svg, $) {

    svg.NS = 'http://www.w3.org/2000/svg';

    svg.createSVG = function (width, height) {
        var svgTag = document.createElementNS(this.NS, 'svg');
        svgTag.setAttribute('width', width);
        svgTag.setAttribute('height', height);
        return svgTag;
    };

    svg.group = function () {
        var svgG = document.createElementNS(this.NS, 'g');
        return svgG;
    };

    svg.line = function (x1, y1, x2, y2, color, epaisseur) {
        var svgL = document.createElementNS(this.NS, 'line');
        svgL.setAttribute('x1', x1);
        svgL.setAttribute('y1', y1);
        svgL.setAttribute('x2', x2);
        svgL.setAttribute('y2', y2);
        svgL.setAttribute('stroke', color);
        svgL.setAttribute('stroke-width', epaisseur);

        return svgL;
    };

    svg.text = function (x, y, text, color, anchor) {
        var svgT = document.createElementNS(this.NS, 'text');
        svgT.setAttribute('x', x);
        svgT.setAttribute('y', y);
        svgT.setAttribute('fill', color);
        svgT.setAttribute('text-anchor', anchor);

        var txt = document.createTextNode(text);
        svgT.appendChild(txt);

        return svgT;
    };

})(SeqView.svg, jQuery);
