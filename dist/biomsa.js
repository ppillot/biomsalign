!function(e,t){"object"==typeof exports&&"undefined"!=typeof module?module.exports=t():"function"==typeof define&&define.amd?define(t):(e="undefined"!=typeof globalThis?globalThis:e||self).biomsa=t()}(this,(function(){"use strict";var e=function(){return(e=Object.assign||function(e){for(var t,r=1,n=arguments.length;r<n;r++)for(var i in t=arguments[r])Object.prototype.hasOwnProperty.call(t,i)&&(e[i]=t[i]);return e}).apply(this,arguments)};function t(e,t){for(var r=0,n=t.length,i=e.length;r<n;r++,i++)e[i]=t[r];return e}function r(e){return 16843009*((e=(858993459&(e-=e>>>1&1431655765))+(e>>>2&858993459))+(e>>>4)&252645135)>>>24}var n,i=function(){function e(e,t){this.length=e,this._words=new Uint32Array(e+32>>>5),!0===t&&this.setAll()}return e.prototype.get=function(e){return 0!=(this._words[e>>>5]&1<<e)},e.prototype.set=function(e){this._words[e>>>5]|=1<<e},e.prototype.clear=function(e){this._words[e>>>5]&=~(1<<e)},e.prototype.flip=function(e){this._words[e>>>5]^=1<<e},e.prototype._assignRange=function(e,t,r){if(!(t<e)){for(var n=this._words,i=!0===r?4294967295:0,o=e>>>5,s=t>>>5,a=o;a<s;++a)n[a]=i;var u=o<<5,h=s<<5;if(!0===r)if(t-e<32)for(var c=e,f=t+1;c<f;++c)n[c>>>5]|=1<<c;else{for(c=e,f=u;c<f;++c)n[c>>>5]|=1<<c;for(c=h,f=t+1;c<f;++c)n[c>>>5]|=1<<c}else if(t-e<32)for(c=e,f=t+1;c<f;++c)n[c>>>5]&=~(1<<c);else{for(c=e,f=u;c<f;++c)n[c>>>5]&=~(1<<c);for(c=h,f=t+1;c<f;++c)n[c>>>5]&=~(1<<c)}return this}},e.prototype.setRange=function(e,t){return this._assignRange(e,t,!0)},e.prototype.clearRange=function(e,t){return this._assignRange(e,t,!1)},e.prototype.setBits=function(){for(var e=[],t=0;t<arguments.length;t++)e[t]=arguments[t];for(var r=this._words,n=e.length,i=0;i<n;++i){var o=e[i];r[o>>>5]|=1<<o}return this},e.prototype.clearBits=function(){for(var e=[],t=0;t<arguments.length;t++)e[t]=arguments[t];for(var r=this._words,n=e.length,i=0;i<n;++i){var o=e[i];r[o>>>5]&=~(1<<o)}return this},e.prototype.setAll=function(){return this._assignRange(0,this.length-1,!0)},e.prototype.clearAll=function(){return this._assignRange(0,this.length-1,!1)},e.prototype.flipAll=function(){for(var e=this._words.length,t=this._words,r=32-this.length%32,n=0;n<e-1;++n)t[n]=~t[n];return t[e-1]=~(t[e-1]<<r)>>>r,this},e.prototype._isRangeValue=function(e,t,r){if(!(t<e)){for(var n=this._words,i=!0===r?4294967295:0,o=e>>>5,s=t>>>5,a=o;a<s;++a)if(n[a]!==i)return!1;if(t-e<32){for(var u=e,h=t+1;u<h;++u)if(!!(n[u>>>5]&1<<u)!==r)return!1}else{var c=s<<5;for(u=e,h=o<<5<<5;u<h;++u)if(!!(n[u>>>5]&1<<u)!==r)return!1;for(u=c,h=t+1;u<h;++u)if(!!(n[u>>>5]&1<<u)!==r)return!1}return!0}},e.prototype.isRangeSet=function(e,t){return this._isRangeValue(e,t,!0)},e.prototype.isRangeClear=function(e,t){return this._isRangeValue(e,t,!1)},e.prototype.isAllSet=function(){return this._isRangeValue(0,this.length-1,!0)},e.prototype.isAllClear=function(){return this._isRangeValue(0,this.length-1,!1)},e.prototype.isSet=function(){for(var e=[],t=0;t<arguments.length;t++)e[t]=arguments[t];for(var r=this._words,n=e.length,i=0;i<n;++i){var o=e[i];if(0==(r[o>>>5]&1<<o))return!1}return!0},e.prototype.isClear=function(){for(var e=[],t=0;t<arguments.length;t++)e[t]=arguments[t];for(var r=this._words,n=e.length,i=0;i<n;++i){var o=e[i];if(0!=(r[o>>>5]&1<<o))return!1}return!0},e.prototype.isEqualTo=function(e){for(var t=this._words,r=e._words,n=Math.min(t.length,r.length),i=0;i<n;++i)if(t[i]!==r[i])return!1;return!0},e.prototype.getSize=function(){for(var e=this._words.length,t=this._words,n=0,i=0;i<e;++i)n+=r(t[i]);return n},e.prototype.difference=function(e){for(var t=this._words,r=e._words,n=Math.min(t.length,r.length),i=0;i<n;++i)t[i]=t[i]&~r[i];for(i=t.length;i<n;++i)t[i]=0;return this},e.prototype.union=function(e){for(var t=this._words,r=e._words,n=Math.min(t.length,r.length),i=0;i<n;++i)t[i]|=r[i];for(i=t.length;i<n;++i)t[i]=0;return this},e.prototype.intersection=function(e){for(var t=this._words,r=e._words,n=Math.min(t.length,r.length),i=0;i<n;++i)t[i]&=r[i];for(i=t.length;i<n;++i)t[i]=0;return this},e.prototype.intersects=function(e){for(var t=this._words,r=e._words,n=Math.min(t.length,r.length),i=0;i<n;++i)if(0!=(t[i]&r[i]))return!0;return!1},e.prototype.getIntersectionSize=function(e){for(var t=this._words,n=e._words,i=Math.min(t.length,n.length),o=0,s=0;s<i;++s)o+=r(t[s]&n[s]);return o},e.prototype.makeIntersection=function(t){var r=this._words,n=t._words,i=Math.min(r.length,n.length),o=new Uint32Array(i),s=Object.create(e.prototype);s._words=o,s.length=Math.min(this.length,t.length);for(var a=0;a<i;++a)o[a]=r[a]&n[a];return s},e.prototype.forEach=function(e){for(var t=this._words.length,n=this._words,i=0,o=0;o<t;++o)for(var s=n[o];0!==s;){var a=s&-s;e((o<<5)+r(a-1),i),s^=a,++i}},e.prototype.toArray=function(){for(var e=this._words,t=new Array(this.getSize()),n=this._words.length,i=0,o=0;o<n;++o)for(var s=e[o];0!==s;){var a=s&-s;t[i++]=(o<<5)+r(a-1),s^=a}return t},e.prototype.toString=function(){return"{"+this.toArray().join(",")+"}"},e.prototype.toSeleString=function(){var e=this.toArray().join(",");return e?"@"+e:"NONE"},e.prototype.clone=function(){var t=Object.create(e.prototype);return t.length=this.length,t._words=new Uint32Array(this._words),t},e}();!function(e){e[e.PROTEIN=0]="PROTEIN",e[e.NUCLEIC=1]="NUCLEIC",e[e.UNSET=2]="UNSET"}(n||(n={}));var o={matrix:[[58,23,-12,-7,-44,10,-23,-14,-14,-27,-17,-8,1,-9,-22,23,15,5,-74,-45,0],[23,224,-67,-63,-50,-30,-29,1,-56,-41,-6,-33,-44,-53,-43,15,2,18,-93,-6,0],[-12,-67,111,59,-104,-4,4,-84,6,-88,-65,48,-13,18,-29,5,-7,-63,-105,-73,0],[-7,-63,59,85,-83,-17,-1,-63,25,-60,-47,15,-12,40,-8,1,-7,-47,-108,-51,0],[-44,-50,-104,-83,144,-93,4,12,-74,36,30,-64,-67,-56,-65,-43,-41,-3,63,104,0],[10,-30,-4,-17,-93,140,-32,-95,-27,-91,-75,4,-36,-29,-32,5,-26,-68,-80,-79,0],[-23,-29,4,-1,4,-32,137,-50,6,-37,-42,21,-23,27,19,-4,-12,-44,-13,48,0],[-14,1,-84,-63,12,-95,-50,86,-53,53,47,-62,-60,-47,-55,-43,-8,69,-27,-24,0],[-14,-56,6,25,-74,-27,6,-53,75,-48,-30,13,-12,34,68,-3,-4,-44,-71,-49,0],[-27,-41,-88,-60,36,-91,-37,53,-48,88,62,-63,-48,-36,-48,-47,-25,36,-11,-4,0],[-17,-6,-65,-47,30,-75,-42,47,-30,62,103,-45,-54,-21,-31,-35,-9,31,-46,-20,0],[-8,-33,48,15,-64,4,21,-62,13,-63,-45,89,-25,12,2,22,10,-51,-79,-29,0],[1,-44,-13,-12,-67,-36,-23,-60,-12,-48,-54,-25,160,-6,-20,5,-12,-42,-76,-83,0],[-9,-53,18,40,-56,-29,27,-47,34,-36,-21,12,-6,75,34,1,-4,-37,-92,-48,0],[-22,-43,-29,-8,-65,-32,19,-55,68,-48,-31,2,-20,34,113,-10,-14,-49,-58,-39,0],[23,15,5,1,-43,5,-4,-43,-3,-47,-35,22,5,1,-10,53,32,-28,-62,-31,0],[15,2,-7,-7,-41,-26,-12,-8,-4,-25,-9,10,-12,-4,-14,32,68,0,-87,-40,0],[5,18,-63,-47,-3,-68,-44,69,-44,36,31,-51,-42,-37,-49,-28,0,74,-61,-32,0],[-74,-93,-105,-108,63,-80,-13,-27,-71,-11,-46,-79,-76,-92,-58,-62,-87,-61,289,81,0],[-45,-6,-73,-51,104,-79,48,-24,-49,-4,-20,-29,-83,-48,-39,-31,-40,-32,81,162,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]],center:22,gapOP:-300},s={matrix:[[91,-114,-31,-123],[-114,100,-125,-31],[-31,-125,100,-114],[-123,-31,-114,91]],gapOP:-400,gapEP:-60,center:120};var a=/^[ACGTRNDEQHILKMFPSWYV]+$/,u=/^[ABCDGHKMNRSTUVWY]+$/i,h=/^[ABCGTRNDEQHIJLKMFPSWYUVXZ]+$/i,c=/[^ACGU]/i,f=/[^ACGT]/i;function p(e){if(!c.test(e)||!f.test(e))return n.NUCLEIC;var t=u.test(e),r=a.test(e);if(r&&!t)return n.PROTEIN;if(r&&t)return function(e){for(var t=100,r=0,i=Math.min(e.length,t),o=0;o<i;o++)switch(e[o]){case"A":case"T":case"U":case"G":case"C":case"N":r++}return r/i>Math.SQRT1_2?n.NUCLEIC:n.PROTEIN}(e);if(h.test(e))return n.PROTEIN;throw new Error("Unrecognized sequence type: "+e)}var d=new Uint8Array(96);d.fill(255),d[65]=0,d[67]=1,d[68]=2,d[69]=3,d[70]=4,d[71]=5,d[72]=6,d[73]=7,d[75]=8,d[76]=9,d[77]=10,d[78]=11,d[80]=12,d[81]=13,d[82]=14,d[83]=15,d[84]=16,d[86]=17,d[87]=18,d[89]=19,d[95]=20;var g=new Uint8Array(96);function l(e,t){var r=Math.floor(100*Math.random());switch(t){case n.NUCLEIC:switch(e){case"M":return[0,1][r%2];case"R":return[0,2][r%2];case"W":return[0,3][r%2];case"S":return[1,2][r%2];case"Y":return[1,3][r%2];case"K":return[2,3][r%2];case"V":return[0,1,2][r%3];case"H":return[0,1,3][r%3];case"D":return[0,2,3][r%3];case"B":return[1,2,3][r%3];case"N":default:return r%4}case n.PROTEIN:default:switch(e){case"B":return[2,11][r%2];case"Z":return[3,13][r%2];case"J":return[7,9][r%2];case"U":case"X":default:return r%20}}}function m(e){return e.map((function(e){return S[e]}))}g.fill(255),g[65]=0,g[67]=1,g[71]=2,g[84]=3,g[85]=3;var S=[0,1,2,2,3,0,4,5,4,5,5,2,0,2,4,0,0,5,3,3];function _(e,t){var r=null!=t?t:p(e),i=function(e,t){for(var r=new Uint8Array(e.length),i=t===n.PROTEIN?d:g,o=0,s=0,a=e.length;s<a;s++)255===(o=i[e.charCodeAt(s)])&&(o=l(e[s],t)),r[s]=o;return r}(e,r);return{rawSeq:e,encodedSeq:i,compressedSeq:r===n.PROTEIN?m(i):i,type:r}}function v(e,t,r){for(var n=[],i=r.gapchar,o=0,s=0;s<t.length;s++){var a=t[s];a<0?n.push(i.repeat(-a)):(n.push(e.substr(o,a)),o+=a)}return n.join("")}function y(e,t){var r=[0],n=0,i=t.slice();function o(e){(e^r[r.length-1])>=0?r[r.length-1]+=e:r.push(e)}for(var s=0;s<e.length;s++){var a=e[s];if(a<0)o(a);else for(;a>0&&n<i.length;){var u=i[n];u<0?a<-u?(i[n]+=a,o(-a),a=0):(o(u),a+=u,n++):a<u?(i[n]-=a,o(a),a=0):(o(u),a-=u,n++)}}return r}function w(e,t){for(var r=[],n=0,i=0,o=e[0],s=t[0];n<e.length||i<t.length;)o!==s?Math.sign(o)!==Math.sign(s)?o<0||void 0===s?(r.push(o),o=e[++n]):(r.push(s),s=t[++i]):o>0?o>s?(r.push(s),o-=s,s=t[++i]):(r.push(o),s-=o,o=e[++n]):(r.push(Math.min(o,s)),o=e[++n],s=t[++i]):(r.push(s),o=e[++n],s=t[++i]);return r}function b(e,t){for(var r=[],n=0,i=0,o=t[0],s=e[0],a=0;n<t.length||i<e.length;)if(a=r[r.length-1],o!==s){if(Math.sign(o)===Math.sign(s)){if(o>0){if(s>o)return;a>0?r[r.length-1]=a+s:r.push(s),o-=s,s=e[++i];continue}if(o>s){a>0?r[r.length-1]=a-o:r.push(Math.abs(o)),s-=o=t[++n];continue}return}if(!(s<0))return;a<0?r[r.length-1]=a+s:r.push(s),s=e[++i]}else a>0?r[r.length-1]=a+Math.abs(o):r.push(Math.abs(o)),o=t[++n],s=e[++i];return r}function O(e,r){if(Math.sign(e[e.length-1])===Math.sign(r[0])){var n=t([],e);return n[n.length-1]+=r[0],n.push.apply(n,r.slice(1)),n}return t(t([],e),r)}function A(e,t){void 0===t&&(t=0);for(var r=[],n=0,i=0,o=0;o<e.length;o++)((i=e[o])^n)>=0?n+=i:(r.push(n),n=i);return r.push(n),0===r[0]&&r.shift(),1&t||r.reverse(),r}function q(e){for(var t=0,r=0,n=0;n<e.length;n++)t+=(r=e[n])<0?-r:r;return t}function C(e){for(var t=[],r=0,n=0,i=0;i<e.length;i++)if((n=e[i])<0)for(var o=0;o<-n;o++)t.push(-1);else for(o=0;o<n;o++)t.push(r),r++;return t}function E(e,t,r,n){void 0===n&&(n=0);var i=e.encodedSeq,o=t.encodedSeq,s=i.length,a=o.length,u=0,h=[],c=[],f=0,p=0,d=0,g=0,l=new Uint8Array(Math.ceil((s+1)*a/2)),m=0,S=0,_=0,v=r.scoringMatrix[0],y=r.gapOP,w=1&n?0:-y/2,b=2&n?0:-y/2;h[0]=0,c[0]=-1/0;for(var O=1;O<=a;O++)h[O]=y+w,c[O]=-1/0;for(var q=1;q<=s;q++){f=w,p=-1/0,v=r.scoringMatrix[i[q-1]];for(O=1;O<=a;O++)_=0,m=h[O]+y,O===a&&(m+=b),m>=c[O]?c[O]=m:_+=1,S=f+y,q===s&&(S+=b),S>=p?p=S:_+=2,u=h[O-1]+v[o[O-1]],h[O-1]=f,u>=p?u>=c[O]?f=u:(f=c[O],_+=4):p>=c[O]?(f=p,_+=8):(f=c[O],_+=4),g=(d=q*a+O)%2,l[d>>>=1]+=g?_:_<<4;h[a]=f}var C=Math.max(u,p,c[c.length-1]),E=P(l,s,a,_),I=E[0],G=E[1];return{estrings:[A(I),A(G)],score:C}}function P(e,t,r,n){for(var i=t,o=r,s=t*r+r,a=1,u=[],h=[],c=(12&n)>>2,f=0;i>0&&o>0;)a=(s=i*r+o)%2,f=e[s>>>1],f=a?15&f:f>>>4,0===c?0===(f>>=2)&&(i--,o--,u.push(1),h.push(1)):2===c?(u.push(-1),h.push(1),f&=2,o--):(u.push(1),h.push(-1),f&=1,i--),c=f;return i>0?(u.push(i),h.push(-i)):o>0&&(u.push(-o),h.push(o)),[u,h]}var I,G=function(){function e(e,t,r,n){void 0===n&&(n=0),this.nbSeq=0,this.length=0,this.weight=0,this.length=e,this.m_bAllGaps=new i(e),this.m_uSortOrder=new Uint8Array(t*e),this.m_fcCounts=new Uint8Array(t*e),this.m_wCounts=new Float32Array(t*e),this.m_AAScores=new Float32Array(t*e),this.m_uResidueGroup=new Uint8Array(e),this.m_fOcc=new Float32Array(e),this.m_fcStartOcc=new Float32Array(e),this.m_fcEndOcc=new Float32Array(e),this.m_ScoreGapOpen=new Float32Array(e),this.m_ScoreGapClose=new Float32Array(e),this.nbSeq=r,this.alphaSize=t,this.weight=n}return e.prototype.getProxy=function(e){return void 0===e&&(e=0),new M(this,e)},e}(),M=function(){function e(e,t){void 0===t&&(t=0),this.idx=0,this.alphaSize=4,this.offset=0,this.prof=e,this.idx=t,this.alphaSize=e.alphaSize,this.offset=this.idx*this.alphaSize}return e.prototype.setProxy=function(e){this.idx=e,this.offset=e*this.alphaSize},Object.defineProperty(e.prototype,"m_bAllGaps",{get:function(){return this.prof.m_bAllGaps.get(this.idx)},set:function(e){e?this.prof.m_bAllGaps.set(this.idx):this.prof.m_bAllGaps.clear(this.idx)},enumerable:!1,configurable:!0}),Object.defineProperty(e.prototype,"m_uSortOrder",{get:function(){return this.prof.m_uSortOrder.subarray(this.offset,this.offset+this.alphaSize)},set:function(e){this.prof.m_uSortOrder.set(e,this.offset)},enumerable:!1,configurable:!0}),Object.defineProperty(e.prototype,"m_uResidueGroup",{get:function(){return this.prof.m_uResidueGroup[this.idx]},set:function(e){this.prof.m_uResidueGroup[this.idx]=e},enumerable:!1,configurable:!0}),Object.defineProperty(e.prototype,"m_fcCounts",{get:function(){return this.prof.m_fcCounts.subarray(this.offset,this.offset+this.alphaSize)},set:function(e){this.prof.m_fcCounts.set(e,this.offset)},enumerable:!1,configurable:!0}),Object.defineProperty(e.prototype,"m_wCounts",{get:function(){return this.prof.m_wCounts.subarray(this.offset,this.offset+this.alphaSize)},set:function(e){this.prof.m_wCounts.set(e,this.offset)},enumerable:!1,configurable:!0}),Object.defineProperty(e.prototype,"m_AAScores",{get:function(){return this.prof.m_AAScores.subarray(this.offset,this.offset+this.alphaSize)},set:function(e){this.prof.m_AAScores.set(e,this.offset)},enumerable:!1,configurable:!0}),Object.defineProperty(e.prototype,"m_fOcc",{get:function(){return this.prof.m_fOcc[this.idx]},set:function(e){this.prof.m_fOcc[this.idx]=e},enumerable:!1,configurable:!0}),Object.defineProperty(e.prototype,"m_fcStartOcc",{get:function(){return this.prof.m_fcStartOcc[this.idx]},set:function(e){this.prof.m_fcStartOcc[this.idx]=e},enumerable:!1,configurable:!0}),Object.defineProperty(e.prototype,"m_fcEndOcc",{get:function(){return this.prof.m_fcEndOcc[this.idx]},set:function(e){this.prof.m_fcEndOcc[this.idx]=e},enumerable:!1,configurable:!0}),Object.defineProperty(e.prototype,"m_ScoreGapOpen",{get:function(){return this.prof.m_ScoreGapOpen[this.idx]},set:function(e){this.prof.m_ScoreGapOpen[this.idx]=e},enumerable:!1,configurable:!0}),Object.defineProperty(e.prototype,"m_ScoreGapClose",{get:function(){return this.prof.m_ScoreGapClose[this.idx]},set:function(e){this.prof.m_ScoreGapClose[this.idx]=e},enumerable:!1,configurable:!0}),e}();function x(e,t,r,n){void 0===n&&(n=0);for(var i=e.encodedSeq,o=i.length,s=new G(o,r.abSize,1,t),a=t,u=r.gapOP,h=s.getProxy(),c=0,f=0;f<o;f++){h.setProxy(f),c=i[f],h.m_fcCounts[c]=1,h.m_wCounts[c]=a,h.m_uSortOrder[h.m_uResidueGroup]=c,h.m_uResidueGroup=1,h.m_fOcc=a,h.m_ScoreGapOpen=u/2*a,h.m_ScoreGapClose=u/2*a;for(var p=0;p<r.abSize;p++)h.m_AAScores[p]=a*r.scoringMatrix[p][c]}return 1&n||(s.m_ScoreGapOpen[0]/=2,s.m_ScoreGapOpen[o-1]/=2),2&n||(s.m_ScoreGapClose[0]/=2,s.m_ScoreGapClose[o-1]/=2),s}function R(e,t,r,n,i,o){void 0===o&&(o=0);var s=q(r),a=e.nbSeq+t.nbSeq,u=e.weight+t.weight,h=new G(s,i.abSize,a,u);console.log(e.weight,t.weight,u);for(var c=i.gapOP,f=h.getProxy(),p=e.getProxy(),d=t.getProxy(),g=C(r),l=C(n),m=0,S=0,_=1,v=1,y=!1,w=!1,b=!1,O=!1,A=0,E=0;E<s;E++){if(f.setProxy(E),m=g[E],p.setProxy(m),-1===m?(y=_>=0,w=E===s-1||g[E+1]>0):(y=!1,w=!1),_=m,S=l[E],d.setProxy(S),-1===S?(b=v>=0,O=E==s-1||l[E+1]>0):(b=!1,O=!1),v=S,A=0,-1!==m&&-1!==S){f.m_fcStartOcc=p.m_fcStartOcc+d.m_fcStartOcc,f.m_fcEndOcc=p.m_fcEndOcc+d.m_fcEndOcc,f.m_fOcc=p.m_fOcc+d.m_fOcc,f.m_fcCounts.set(p.m_fcCounts),f.m_wCounts.set(p.m_wCounts);for(var P=0;P<t.alphaSize;P++){var I=t.m_fcCounts[P];I&&(f.m_fcCounts[P]+=I,f.m_wCounts[P]+=t.m_wCounts[P]),f.m_fcCounts[P]&&(f.m_uSortOrder[A]=P,A++),f.m_AAScores[P]=p.m_AAScores[P]+d.m_AAScores[P]}f.m_uResidueGroup=A}else-1==m?(f.m_fcStartOcc=d.m_fcStartOcc+(y?e.weight:0),f.m_fcEndOcc=d.m_fcEndOcc+(w?e.weight:0),f.m_fOcc=d.m_fOcc,f.m_fcCounts.set(d.m_fcCounts),f.m_wCounts.set(d.m_wCounts),f.m_uSortOrder.set(d.m_uSortOrder),f.m_AAScores.set(d.m_AAScores),f.m_uResidueGroup=d.m_uResidueGroup):-1==S&&(f.m_fcStartOcc=p.m_fcStartOcc+(b?t.weight:0),f.m_fcEndOcc=p.m_fcEndOcc+(O?t.weight:0),f.m_fOcc=p.m_fOcc,f.m_fcCounts.set(p.m_fcCounts),f.m_wCounts.set(p.m_wCounts),f.m_uSortOrder.set(p.m_uSortOrder),f.m_AAScores.set(p.m_AAScores),f.m_uResidueGroup=p.m_uResidueGroup);f.m_ScoreGapOpen=c/2*(1-f.m_fcStartOcc)*f.m_fOcc,f.m_ScoreGapClose=c/2*(1-f.m_fcEndOcc)*f.m_fOcc}return 1&o||(h.m_ScoreGapOpen[0]/=2,h.m_ScoreGapOpen[s-1]/=2),2&o||(h.m_ScoreGapClose[0]/=2,h.m_ScoreGapClose[s-1]/=2),h}function z(e){return"seq"in e}function T(e,t){for(var r=1/0,n=0,i=0,o=e.length,s=0;s<o;s++)for(var a=s+1;a<o;a++)e[s][a]<r&&(r=e[s][a],n=s,i=a);return[{profile:null,childA:t[n],childB:t[i],distance:r,estring:[],msa:[],numSeq:[],type:I.NODE,depth:0,id:"",tabWeight:[],parent:-1,weight:0},n,i]}function U(e,t,r,n){for(var i=t>r?[t,r]:[r,t],o=i[0],s=i[1],a=[],u=e.length,h=0;h<u;h++)if(h!=t&&h!=r){var c=.1*(e[h][t]+e[h][r])/2+.9*Math.min(e[h][t],e[h][r]);e[h].push(c),a.push(c),e[h].splice(o,1),e[h].splice(s,1)}return a.push(0),e.push(a),e.splice(o,1),e.splice(s,1),n.splice(o,1),n.splice(s,1),e}function N(e,r){var o,s,a,u=function(e,n){var i=h[e.childA],o=h[e.childB],s={};if(z(i)||0!==i.msa.length||u(i,n),z(o)||0!==o.msa.length||u(o,n),z(i))if(z(o)){var a=E(i.seq,o.seq,r);s.score=a.score,e.numSeq=[i.numSeq[0],o.numSeq[0]],e.estring=a.estrings,e.profile=R(i.profile,o.profile,a.estrings[0],a.estrings[1],r)}else{o.tabWeight=o.numSeq.map((function(e){return n[e]}));var c=function(e,t,r,n){void 0===n&&(n=0);var i,o,s=t.seq,a=s.rawSeq.length,u=e.profile.length,h=0,c=[],f=[],p=new Uint8Array(Math.ceil((u+1)*(a+1)/2)),d=0,g=0,l=0,m=0,S=0,_=0,v=0,y=0,w=0,b=0,O=new Float32Array(r.abSize),q=r.gapOP,C=q/2,E=q/2,I=1&n?0:-q/4,G=2&n?0:-q/4,M=s.encodedSeq,x=e.profile;c[0]=0,f[0]=-1/0;for(var R=1;R<=a;R++)c[R]=C+E+I,f[R]=-1/0;for(var z=1;z<=u;z++){for(w=x.m_ScoreGapOpen[z-1],b=x.m_ScoreGapClose[z-1],O=x.m_AAScores.subarray((z-1)*r.abSize,z*r.abSize),_=I,m=-1/0,f[0]=x.m_ScoreGapOpen[0],R=1;R<=a;R++)l=0,d=c[R]+w,R===a&&(d-=w/2),d>=f[R]?f[R]=d:l+=1,g=_+C,z===u&&(g+=G),g>=(S=m)?m=g:l+=2,v=f[R]+b,y=S+E,h=c[R-1]+O[M[R-1]],R===a&&(y+=G),c[R-1]=_,h>=y?h>=v?_=h:(_=v,l+=4):y>=v?(_=y,l+=8):(_=v,l+=4),o=(i=z*a+R)%2,p[i>>>=1]+=o?l:l<<4;c[a]=_}var T=Math.max(h,m,f[f.length-2]),U=P(p,u,a,l),N=U[0],j=U[1];return{estrings:[A(N),A(j)],score:T}}(o,i,r);s.score=c.score,e.numSeq=t([i.numSeq[0]],o.numSeq),e.estring=t([c.estrings[1]],o.estring.map((function(e){return y(c.estrings[0],e)}))),e.profile=R(i.profile,o.profile,c.estrings[1],c.estrings[0],r)}else if(!z(o)){o.tabWeight=o.numSeq.map((function(e){return n[e]})),i.tabWeight=i.numSeq.map((function(e){return n[e]}));var f=function(e,t,r,n){void 0===n&&(n=0);var i,o,s=e.profile.length,a=t.profile.length,u=new Uint8Array(Math.ceil((s+1)*(a+1)/2)),h=[],c=[],f=0,p=0,d=0,g=0,l=0,m=0,S=0,_=0,v=0,y=0,w=0,b=0,O=new Float32Array(r.abSize),q=new Uint8Array(r.abSize),C=0,E=0,I=0,G=0,M=0,x=0,R=1&n?1:2,z=t.profile,T=e.profile;for(h[0]=0,c[0]=-1/0,M=1;M<=a;M++)h[M]=z.m_ScoreGapOpen[0]+z.m_ScoreGapClose[M-1]/R,c[M]=-1/0;var U=T.m_ScoreGapOpen[0];for(G=1;G<=s;G++){for(O=T.m_AAScores.subarray(r.abSize*(G-1),r.abSize*G),w=T.m_ScoreGapOpen[G-1],_=U+(b=T.m_ScoreGapClose[G-1])/R,m=-1/0,c[0]=T.m_ScoreGapOpen[0],M=1;M<=a;M++){for(l=0,p=h[M]+w,M!==a||2&n||(p-=w/2),p>=c[M]?c[M]=p:l+=1,d=_+z.m_ScoreGapOpen[M-1],g=S=m,G!==s||2&n||(d-=z.m_ScoreGapOpen[M-1]/2),d>=g?m=d:l+=2,f=h[M-1],I=z.m_uResidueGroup[M-1],x=0,E=(M-1)*r.abSize,q=z.m_uSortOrder.subarray(E,E+I);x<I;)C=q[x],f+=z.m_wCounts[E+C]*O[C],x++;v=c[M]+b,y=S+z.m_ScoreGapClose[M-1],1!==M||2&n||(v-=b/2),M!==a||2&n||(v-=b/2,y-=z.m_ScoreGapClose[a-1]/2),h[M-1]=_,f>=y?f>=v?_=f:(_=v,l+=4):y>=v?(_=y,l+=8):(_=v,l+=4),o=(i=G*a+M)%2,u[i>>>=1]+=o?l:l<<4}h[a]=_}var N=Math.max(f,m,c[M-1]),j=P(u,s,a,l),D=j[0],F=j[1];return{estrings:[A(D),A(F)],score:N}}(i,o,r);s.score=f.score,e.estring=t(t([],i.estring.map((function(e){return y(f.estrings[0],e)}))),o.estring.map((function(e){return y(f.estrings[1],e)}))),e.profile=R(i.profile,o.profile,f.estrings[0],f.estrings[1],r)}return e.numSeq=t(t([],i.numSeq),o.numSeq),s.score},h=function(e,t){for(var r,n,i,o,s=[],a=e.length,u=[],h=0;h<a;h++)s[h]={type:I.LEAF,seq:t[h],profile:null,childA:h,childB:h,distance:0,numSeq:[h],msa:[],id:h.toString(),weight:0,parent:-1,depth:0},u[h]=h;h=0;for(var c=a-1;h<c;h++){n=(r=T(e,u))[0],i=r[1],o=r[2];var f=s[n.childA],p=s[n.childB];n.id=f.id<p.id?"|"+f.id+","+p.id+"|":"|"+p.id+","+f.id+"|",n.depth=Math.max(f.depth,p.depth)+1,p.parent=f.parent=s.length,u.push(a+h),s.push(n),e=U(e,i,o,u)}return s[s.length-1].type=I.ROOT,s}(function(e){for(var t,r,o,s=e[0].type===n.PROTEIN,a=s?6:4,u=Math.pow(a,6),h=e.map((function(e){for(var t=new i(u),r=s?e.compressedSeq:e.encodedSeq,n=0,o=0;o<=5;o++)n+=r[o]*Math.pow(a,o);t.set(n);for(var h=Math.pow(a,5),c=(o=6,r.length);o<c;o++)n-=r[o-6],n/=6,n+=r[o]*h,t.set(n);return t})),c=e.length,f=e.map((function(){return[]})),p=0;p<c;p++){f[p][p]=0,t=h[p],r=e[p].compressedSeq.length;for(var d=p+1;d<c;d++)o=1-t.getIntersectionSize(h[d])/r,f[d][p]=f[p][d]=o}return f}(e),e),c=h[h.length-1],f=function(e){var t=[],r=0;!function t(n){var i=0,o=0;switch(n.type){case I.ROOT:n.weight=0,t(e[n.childA]),t(e[n.childB]);break;case I.NODE:i=e[n.parent].distance-n.distance,o=n.id.split(",").length,n.weight=e[n.parent].weight+i/o,t(e[n.childA]),t(e[n.childB]);break;case I.LEAF:i=e[n.parent].distance,n.weight=e[n.parent].weight+i,r+=n.weight}}(e[e.length-1]);for(var n=0;e[n].type==I.LEAF;)t[n]=e[n].weight/=r,n++;return t}(h);return function(e,t,r){for(var n=0,i=e[n];z(i)&&n<e.length;)i.profile=x(i.seq,i.weight,t,r),i=e[++n]}(h,r),u(c,f),o=c.estring.slice(),s=c.numSeq,a=s.slice(),s.forEach((function(e,t){return a[e]=t})),a.map((function(e){return o[e]}))}!function(e){e[e.LEAF=0]="LEAF",e[e.NODE=1]="NODE",e[e.ROOT=2]="ROOT"}(I||(I={}));var j=function(){function e(e){this.storeSize=e,this.store=new Array(e),this.head=0,this.tail=0}return Object.defineProperty(e.prototype,"size",{get:function(){return(this.tail-this.head+this.storeSize)%this.storeSize},enumerable:!1,configurable:!0}),Object.defineProperty(e.prototype,"isEmpty",{get:function(){return this.head===this.tail},enumerable:!1,configurable:!0}),e.prototype.popHead=function(){if(this.isEmpty)return null;var e=this.store[this.head];return this.head++,this.head>=this.storeSize&&(this.head=0),e},e.prototype.popTail=function(){return this.isEmpty?null:(this.tail--,this.tail<0&&(this.tail=this.storeSize-1),this.store[this.tail])},e.prototype.pushHead=function(e){this.head--,this.head<0&&(this.head=this.storeSize-1),this.store[this.head]=e,this.head===this.tail&&(this.tail=(this.tail-1+this.storeSize)%this.storeSize)},e.prototype.pushTail=function(e){this.store[this.tail]=e,this.tail++,this.tail>=this.storeSize&&(this.tail=0),this.head===this.tail&&(this.head=(this.head-1+this.storeSize)%this.storeSize)},e.prototype.getHead=function(e){void 0===e&&(e=0);var t=0===e?this.head:(this.head+e+this.storeSize)%this.storeSize;return this.store[t]},e.prototype.getTail=function(e){void 0===e&&(e=0);var t=(this.tail-1-e+this.storeSize)%this.storeSize;return this.store[t]},e}();function D(e,t,r){for(var n=new Map,i=new j(r),o=0|t,s=new Uint16Array(e.encodedSeq.length-o+1),a=0,u=0,h=0;h<o;h++)a|=e.encodedSeq[o-h-1]<<2*h;s[u++]=a;h=o;for(var c=e.encodedSeq.length;h<c;h++)a=(a<<2)+e.encodedSeq[h],s[u++]=a;var f=r-o,p=e.rawSeq.length*(2/(r+1))*2|0,d={kmer:new Uint16Array(p),kmerPos:new Uint16Array(p),winPos:new Uint16Array(p),winPosEnd:new Uint16Array(p),count:0},g=NaN,l=0,m=-f-1;for(h=0;h<s.length;h++){var S=s[h];for(m++;!i.isEmpty&&((l=i.getTail())<m||s[l]>=S);)i.popTail();if(i.pushTail(h),!(m<0)){for(;i.getHead()<m;)i.popHead();var _=i.getHead();if(_===g)d.winPosEnd[d.count-1]=h+o;else{var v=s[_],y=d.count;if(d.count++,d.kmer[y]=v,d.kmerPos[y]=_,d.winPos[y]=m,d.winPosEnd[y]=h+o,g=_,n.has(v))n.get(v).push(y);else n.set(v,[y])}}}return[n,d,s]}function F(e,t,r,n,i){var o=null!=n?n:D(e,8,24);o[0];var s=o[1],a=o[2],u=null!=i?i:D(t,8,24),h=u[0],c=u[1],f=u[2];Math.max(e.rawSeq.length,t.rawSeq.length);for(var p=new Map,d=0;d<s.count;d++){var g=s.kmer[d];if(h.has(g)){var l=h.get(g);if(!(l.length>4))for(var m=e.rawSeq.substring(s.winPos[d],s.winPosEnd[d]),S=0;S<l.length;S++){var _=l[S],v=t.rawSeq.substring(c.winPos[_],c.winPosEnd[_]),y=Math.min(m.length,v.length),w=0,b=void 0;if(b=(y==m.length?0===v.indexOf(m):0===m.indexOf(v))?{diagId:w=s.winPos[d]-c.winPos[_],begin:s.winPos[d],end:s.winPos[d]+y}:{diagId:w=s.kmerPos[d]-c.kmerPos[_],begin:s.kmerPos[d],end:s.kmerPos[d]+8-1},p.has(w))p.get(w).push(b);else p.set(w,[b])}}}var C=[];if(p.forEach((function(e){if(e.length){var t=e[0];C.push(t),e.forEach((function(e){e.begin<=t.end?t.end=e.end:(t=e,C.push(t))}))}})),0===C.length)return E(e,t,r).estrings;C.sort((function(e,t){var r=e.begin+e.diagId-t.begin-t.diagId;return 0===r?e.diagId-t.diagId:r})),C=function(e,t,r){var n,i=1/Math.max(t,r),o=t-r,s=[0],a=new Int16Array(e.length);a[0]=-1;for(var u=e[0],h=u.end-u.begin,c=u.end-u.diagId,f=u.begin-u.diagId,p=Math.min(r-c,t-u.end),d=[h],g=[c],l=[p],m=e[0],S=0,_=e[S],v=h,y=0,w=1;w<e.length;w++){if(h=(u=e[w]).end-u.begin,c=u.end-u.diagId,f=u.begin-u.diagId,p=u.diagId<o?r-c:t-u.end,20===y){y=0;for(var b=[],O=s.length-1,A=g[O],q=t-u.begin;O--;)g[O]>A?b.push(O):(A=g[O],d[O]+Math.min(l[O],q)<v&&b.push(O));b.forEach((function(e){s.splice(e,1),d.splice(e,1),g.splice(e,1),l.splice(e,1)}))}if(u.diagId>=m.diagId&&c<g[0]){if(h+p<v)continue;if(a[w]=-1,m=e[w],y++,h<d[0]){s.unshift(w),d.unshift(h),g.unshift(c),l.unshift(p);continue}h>v&&(v=h),s[0]=w,d[0]=h,g[0]=c,l[0]=p}else{for(var C=s.length-1;C>=0;C--)if(_=e[S=s[C]],u.begin>=_.end&&g[C]<f){var E=Math.abs(u.diagId-S)*i,P=d[C]+h-E;if(P+p<v)continue;for(var I=C;I<s.length&&d[I]<P;)I++;if(I===s.length)v=P,l.push(p),d.push(P),g.push(c),s.push(w);else{if(c>=g[g.length-1])continue;l.splice(I,0,p),d.splice(I,0,P),g.splice(I,0,c),s.splice(I,0,w)}a[w]=S,y++;break}h>v&&(a[w]=-1,v=h,l.push(p),d.push(h),g.push(c),s.push(w),y++)}}var G=[],M=null!==(n=s.pop())&&void 0!==n?n:0,x=M;for(;M>-1&&x>-1;)G.push(M),M=a[M],x--;var R=[];for(w=G.length-1;w>=0;w--)R.push(e[G[w]]);return R}(C,e.encodedSeq.length,t.encodedSeq.length);var P={diagId:0,begin:0,end:0},I=[P],G=[];for(d=0;d<C.length;d++){var M=C[d];if(M.diagId!==P.diagId){for(R=(x=P.end)-P.diagId;e.encodedSeq[x]===t.encodedSeq[R]&&x<M.begin&&R<M.begin-M.diagId;)x++,R++;for(P.end=x,R=(x=M.begin)-M.diagId;e.encodedSeq[x]===t.encodedSeq[R]&&x>P.end&&R>P.end-P.diagId;)x--,R--;M.begin=x+1,G.push({endDiagId:M.diagId,beginDiagId:P.diagId,begin:P.end,end:M.begin}),I.push(M),P=M}else{for(var x,R=(x=P.end)-M.diagId,z=0;x<M.begin&&R<M.begin-M.diagId&&z<2;)R+=8,z=L(a[x+=8],f[R])<=2?0:z+1;for(x>P.end&&(P.end=Math.max(x-2-8,P.end)),R=(x=M.begin-8)-M.diagId,z=0;x>P.end&&R>P.end-P.diagId&&z<2;)R-=8,z=L(a[x-=8],f[R])<=2?0:z+1;if(x<M.begin-16&&(M.begin=x+24+2),M.begin-P.end<=24){P.end=M.end;continue}G.push({beginDiagId:P.diagId,endDiagId:M.diagId,begin:P.end,end:M.begin}),I.push(M),P=M}}var T=I[0],U=[],N=[];for(d=0;d<I.length;d++){T=I[d],U.push(T.end-T.begin),N.push(T.end-T.begin);var j=G[d];if(j){if(j.begin===j.end){var F=Math.abs(j.endDiagId-j.beginDiagId);U.push(-F),N.push(F);continue}if(j.begin-j.beginDiagId==j.end-j.endDiagId){F=j.end-j.begin;N.push(-F),U.push(F);continue}var B=E({rawSeq:e.rawSeq.substring(j.begin,j.end),type:e.type,compressedSeq:new Uint8Array(0),encodedSeq:e.encodedSeq.subarray(j.begin,j.end)},{rawSeq:t.rawSeq.substring(j.begin-j.beginDiagId,j.end-j.endDiagId),type:t.type,compressedSeq:new Uint8Array(0),encodedSeq:t.encodedSeq.subarray(j.begin-j.beginDiagId,j.end-j.endDiagId)},r,3);U.push.apply(U,B.estrings[0]),N.push.apply(N,B.estrings[1])}}k(U,e.encodedSeq.length),k(N,t.encodedSeq.length);var H=A(U,1),V=A(N,1),W=q(H),K=q(V);return W>K?V=O(V,[K-W]):W<K&&(H=O(H,[W-K])),[H,V]}function k(e,t){var r=t-function(e){for(var t=0,r=0,n=0;n<e.length;n++)(r=e[n])<0||(t+=r);return t}(e);r>0?e.push(r):r<0&&(e[e.length-1]+=r)}function L(e,t){var n=e^t;return r(n=1431655765&n|n>>1&1431655765)}var B={sequenceType:"auto",gapchar:"-",alignmentMethod:"auto",debug:!1};return new(function(){function t(){this.sequences=[],this.typeSeq=n.UNSET,this.config=e({},B)}return t.prototype.addSequence=function(e){var t=this;if(Array.isArray(e))e.forEach((function(e){t.addSequence(e)}));else{if("string"!=typeof e)throw new TypeError("String type expected for sequences to add.");if(e=e.toUpperCase(),"auto"===this.config.sequenceType){var r=p(e);if(this.typeSeq===n.UNSET)this.typeSeq=r;else if(this.typeSeq!==r)throw new Error("All sequences must be of same type.")}this.sequences.push(_(e))}},t.prototype.reset=function(){this.sequences=[],this.typeSeq=n.UNSET,this.setDefaultConfiguration()},t.prototype.setDefaultConfiguration=function(){this.config=e({},B)},t.prototype.setUserConfiguration=function(e){var t;function r(e,t){return t.some((function(t){return t==e}))}e.method&&r(e.method,["complete","diag"])&&(this.config.alignmentMethod=e.method),e.type&&r(e.type,["amino","nucleic"])&&(this.config.sequenceType=e.type,this.typeSeq="amino"===e.type?n.PROTEIN:n.NUCLEIC),1===(null===(t=e.gapchar)||void 0===t?void 0:t.length)&&(this.config.gapchar=e.gapchar),void 0!==e.debug&&(this.config.debug=!!e.debug)},t.prototype.align=function(e,t){var r=this;return new Promise((function(i,a){if(!Array.isArray(e)||e.some((function(e){return"string"!=typeof e})))return a("Array of sequences expected");if(e.length<2)return a("At least 2 sequences are required.");r.reset(),t&&r.setUserConfiguration(t),r.addSequence(e);var u,h=function(e,t){var r,i,a,u=e===n.PROTEIN?o:s,h=null!==(r=null==t?void 0:t.matrix)&&void 0!==r?r:u.matrix,c=(null==t?void 0:t.gapextend)?2*-t.gapextend:u.center;return{type:e,scoringMatrix:(a=c,h.map((function(e){return e.map((function(e){return e+a}))}))),gapOP:null!==(i=null==t?void 0:t.gapopen)&&void 0!==i?i:u.gapOP,abSize:e===n.PROTEIN?20:4}}(r.typeSeq,t),c=!1;(c="auto"===r.config.alignmentMethod?r.sequences.some((function(e){return e.rawSeq.length>1600})):"diag"===r.config.alignmentMethod,2==r.sequences.length)?u=c?F(r.sequences[0],r.sequences[1],h):E(r.sequences[0],r.sequences[1],h).estrings:u=c?function(e,t){for(var r=e.map((function(e){return D(e,8,24)})),n=[],i=[],o=[],s=0;s<r.length;s++)if(0!==s){var a=F(e[0],e[s],t,r[0],r[s]);i=0===i.length?a[0]:w(i,a[0]),o[s]=a[0],n[s]=a[1]}for(s=0;s<r.length;s++)if(0!==s){var u=b(i,o[s]);if(!u)throw console.error("An error occured while computing the center diff. for sequence #"+s,i,o[s]),console.dir(n),console.dir(o),new RangeError;n[s]=y(u,n[s])}else n[s]=i;return n}(r.sequences,h):N(r.sequences,h);return i(function(e,t,r){for(var n,i=null!==(n=null==r?void 0:r.gapchar)&&void 0!==n?n:"-",o=[],s=0;s<e.length;s++)o.push(v(e[s].rawSeq,t[s],{gapchar:i}));return o}(r.sequences,u,{gapchar:r.config.gapchar}))})).catch((function(e){return Promise.reject(e)}))},t}())}));
//# sourceMappingURL=biomsa.js.map
