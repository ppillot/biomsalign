!function(e,t){"object"==typeof exports&&"undefined"!=typeof module?module.exports=t():"function"==typeof define&&define.amd?define(t):(e="undefined"!=typeof globalThis?globalThis:e||self).biomsa=t()}(this,(function(){"use strict";var e=function(){return(e=Object.assign||function(e){for(var t,r=1,n=arguments.length;r<n;r++)for(var i in t=arguments[r])Object.prototype.hasOwnProperty.call(t,i)&&(e[i]=t[i]);return e}).apply(this,arguments)};function t(e,t){for(var r=0,n=t.length,i=e.length;r<n;r++,i++)e[i]=t[r];return e}function r(e){return 16843009*((e=(858993459&(e-=e>>>1&1431655765))+(e>>>2&858993459))+(e>>>4)&252645135)>>>24}var n,i=function(){function e(e,t){this.length=e,this._words=new Uint32Array(e+32>>>5),!0===t&&this.setAll()}return e.prototype.get=function(e){return 0!=(this._words[e>>>5]&1<<e)},e.prototype.set=function(e){this._words[e>>>5]|=1<<e},e.prototype.clear=function(e){this._words[e>>>5]&=~(1<<e)},e.prototype.flip=function(e){this._words[e>>>5]^=1<<e},e.prototype._assignRange=function(e,t,r){if(!(t<e)){for(var n=this._words,i=!0===r?4294967295:0,o=e>>>5,s=t>>>5,a=o;a<s;++a)n[a]=i;var u=o<<5,c=s<<5;if(!0===r)if(t-e<32)for(var h=e,f=t+1;h<f;++h)n[h>>>5]|=1<<h;else{for(h=e,f=u;h<f;++h)n[h>>>5]|=1<<h;for(h=c,f=t+1;h<f;++h)n[h>>>5]|=1<<h}else if(t-e<32)for(h=e,f=t+1;h<f;++h)n[h>>>5]&=~(1<<h);else{for(h=e,f=u;h<f;++h)n[h>>>5]&=~(1<<h);for(h=c,f=t+1;h<f;++h)n[h>>>5]&=~(1<<h)}return this}},e.prototype.setRange=function(e,t){return this._assignRange(e,t,!0)},e.prototype.clearRange=function(e,t){return this._assignRange(e,t,!1)},e.prototype.setBits=function(){for(var e=[],t=0;t<arguments.length;t++)e[t]=arguments[t];for(var r=this._words,n=e.length,i=0;i<n;++i){var o=e[i];r[o>>>5]|=1<<o}return this},e.prototype.clearBits=function(){for(var e=[],t=0;t<arguments.length;t++)e[t]=arguments[t];for(var r=this._words,n=e.length,i=0;i<n;++i){var o=e[i];r[o>>>5]&=~(1<<o)}return this},e.prototype.setAll=function(){return this._assignRange(0,this.length-1,!0)},e.prototype.clearAll=function(){return this._assignRange(0,this.length-1,!1)},e.prototype.flipAll=function(){for(var e=this._words.length,t=this._words,r=32-this.length%32,n=0;n<e-1;++n)t[n]=~t[n];return t[e-1]=~(t[e-1]<<r)>>>r,this},e.prototype._isRangeValue=function(e,t,r){if(!(t<e)){for(var n=this._words,i=!0===r?4294967295:0,o=e>>>5,s=t>>>5,a=o;a<s;++a)if(n[a]!==i)return!1;if(t-e<32){for(var u=e,c=t+1;u<c;++u)if(!!(n[u>>>5]&1<<u)!==r)return!1}else{var h=s<<5;for(u=e,c=o<<5<<5;u<c;++u)if(!!(n[u>>>5]&1<<u)!==r)return!1;for(u=h,c=t+1;u<c;++u)if(!!(n[u>>>5]&1<<u)!==r)return!1}return!0}},e.prototype.isRangeSet=function(e,t){return this._isRangeValue(e,t,!0)},e.prototype.isRangeClear=function(e,t){return this._isRangeValue(e,t,!1)},e.prototype.isAllSet=function(){return this._isRangeValue(0,this.length-1,!0)},e.prototype.isAllClear=function(){return this._isRangeValue(0,this.length-1,!1)},e.prototype.isSet=function(){for(var e=[],t=0;t<arguments.length;t++)e[t]=arguments[t];for(var r=this._words,n=e.length,i=0;i<n;++i){var o=e[i];if(0==(r[o>>>5]&1<<o))return!1}return!0},e.prototype.isClear=function(){for(var e=[],t=0;t<arguments.length;t++)e[t]=arguments[t];for(var r=this._words,n=e.length,i=0;i<n;++i){var o=e[i];if(0!=(r[o>>>5]&1<<o))return!1}return!0},e.prototype.isEqualTo=function(e){for(var t=this._words,r=e._words,n=Math.min(t.length,r.length),i=0;i<n;++i)if(t[i]!==r[i])return!1;return!0},e.prototype.getSize=function(){for(var e=this._words.length,t=this._words,n=0,i=0;i<e;++i)n+=r(t[i]);return n},e.prototype.difference=function(e){for(var t=this._words,r=e._words,n=Math.min(t.length,r.length),i=0;i<n;++i)t[i]=t[i]&~r[i];for(i=t.length;i<n;++i)t[i]=0;return this},e.prototype.union=function(e){for(var t=this._words,r=e._words,n=Math.min(t.length,r.length),i=0;i<n;++i)t[i]|=r[i];for(i=t.length;i<n;++i)t[i]=0;return this},e.prototype.intersection=function(e){for(var t=this._words,r=e._words,n=Math.min(t.length,r.length),i=0;i<n;++i)t[i]&=r[i];for(i=t.length;i<n;++i)t[i]=0;return this},e.prototype.intersects=function(e){for(var t=this._words,r=e._words,n=Math.min(t.length,r.length),i=0;i<n;++i)if(0!=(t[i]&r[i]))return!0;return!1},e.prototype.getIntersectionSize=function(e){for(var t=this._words,n=e._words,i=Math.min(t.length,n.length),o=0,s=0;s<i;++s)o+=r(t[s]&n[s]);return o},e.prototype.makeIntersection=function(t){var r=this._words,n=t._words,i=Math.min(r.length,n.length),o=new Uint32Array(i),s=Object.create(e.prototype);s._words=o,s.length=Math.min(this.length,t.length);for(var a=0;a<i;++a)o[a]=r[a]&n[a];return s},e.prototype.forEach=function(e){for(var t=this._words.length,n=this._words,i=0,o=0;o<t;++o)for(var s=n[o];0!==s;){var a=s&-s;e((o<<5)+r(a-1),i),s^=a,++i}},e.prototype.toArray=function(){for(var e=this._words,t=new Array(this.getSize()),n=this._words.length,i=0,o=0;o<n;++o)for(var s=e[o];0!==s;){var a=s&-s;t[i++]=(o<<5)+r(a-1),s^=a}return t},e.prototype.toString=function(){return"{"+this.toArray().join(",")+"}"},e.prototype.toSeleString=function(){var e=this.toArray().join(",");return e?"@"+e:"NONE"},e.prototype.clone=function(){var t=Object.create(e.prototype);return t.length=this.length,t._words=new Uint32Array(this._words),t},e}();!function(e){e[e.PROTEIN=0]="PROTEIN",e[e.NUCLEIC=1]="NUCLEIC",e[e.UNSET=2]="UNSET"}(n||(n={}));var o,s={matrix:[[58,23,-12,-7,-44,10,-23,-14,-14,-27,-17,-8,1,-9,-22,23,15,5,-74,-45,0],[23,224,-67,-63,-50,-30,-29,1,-56,-41,-6,-33,-44,-53,-43,15,2,18,-93,-6,0],[-12,-67,111,59,-104,-4,4,-84,6,-88,-65,48,-13,18,-29,5,-7,-63,-105,-73,0],[-7,-63,59,85,-83,-17,-1,-63,25,-60,-47,15,-12,40,-8,1,-7,-47,-108,-51,0],[-44,-50,-104,-83,144,-93,4,12,-74,36,30,-64,-67,-56,-65,-43,-41,-3,63,104,0],[10,-30,-4,-17,-93,140,-32,-95,-27,-91,-75,4,-36,-29,-32,5,-26,-68,-80,-79,0],[-23,-29,4,-1,4,-32,137,-50,6,-37,-42,21,-23,27,19,-4,-12,-44,-13,48,0],[-14,1,-84,-63,12,-95,-50,86,-53,53,47,-62,-60,-47,-55,-43,-8,69,-27,-24,0],[-14,-56,6,25,-74,-27,6,-53,75,-48,-30,13,-12,34,68,-3,-4,-44,-71,-49,0],[-27,-41,-88,-60,36,-91,-37,53,-48,88,62,-63,-48,-36,-48,-47,-25,36,-11,-4,0],[-17,-6,-65,-47,30,-75,-42,47,-30,62,103,-45,-54,-21,-31,-35,-9,31,-46,-20,0],[-8,-33,48,15,-64,4,21,-62,13,-63,-45,89,-25,12,2,22,10,-51,-79,-29,0],[1,-44,-13,-12,-67,-36,-23,-60,-12,-48,-54,-25,160,-6,-20,5,-12,-42,-76,-83,0],[-9,-53,18,40,-56,-29,27,-47,34,-36,-21,12,-6,75,34,1,-4,-37,-92,-48,0],[-22,-43,-29,-8,-65,-32,19,-55,68,-48,-31,2,-20,34,113,-10,-14,-49,-58,-39,0],[23,15,5,1,-43,5,-4,-43,-3,-47,-35,22,5,1,-10,53,32,-28,-62,-31,0],[15,2,-7,-7,-41,-26,-12,-8,-4,-25,-9,10,-12,-4,-14,32,68,0,-87,-40,0],[5,18,-63,-47,-3,-68,-44,69,-44,36,31,-51,-42,-37,-49,-28,0,74,-61,-32,0],[-74,-93,-105,-108,63,-80,-13,-27,-71,-11,-46,-79,-76,-92,-58,-62,-87,-61,289,81,0],[-45,-6,-73,-51,104,-79,48,-24,-49,-4,-20,-29,-83,-48,-39,-31,-40,-32,81,162,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]],center:22,gapOP:-300},a={matrix:[[91,-114,-31,-123],[-114,100,-125,-31],[-31,-125,100,-114],[-123,-31,-114,91]],gapOP:-400,gapEP:-60,center:120};function u(e){o.has(e)?u(e+"_"):o.set(e,performance.now())}var c={start:function(){o=new Map,u("START")},add:u,summary:function(){var e={start:0};u("END");var t=0;o.forEach((function(r,n){"START"!==n?(e[n]=r-t,t=r):t=r})),console.table(e)}},h=/^[ACGTRNDEQHILKMFPSWYV]+$/,f=/^[ABCDGHKMNRSTUVWY]+$/i,p=/^[ABCGTRNDEQHIJLKMFPSWYUVXZ]+$/i,d=/[^ACGU]/i,l=/[^ACGT]/i;function g(e){if(!d.test(e)||!l.test(e))return n.NUCLEIC;var t=f.test(e),r=h.test(e);if(r&&!t)return n.PROTEIN;if(r&&t)return function(e){for(var t=100,r=0,i=Math.min(e.length,t),o=0;o<i;o++)switch(e[o]){case"A":case"T":case"U":case"G":case"C":case"N":r++}return r/i>Math.SQRT1_2?n.NUCLEIC:n.PROTEIN}(e);if(p.test(e))return n.PROTEIN;throw new Error("Unrecognized sequence type: "+e)}var m=new Uint8Array(96);m.fill(255),m[65]=0,m[67]=1,m[68]=2,m[69]=3,m[70]=4,m[71]=5,m[72]=6,m[73]=7,m[75]=8,m[76]=9,m[77]=10,m[78]=11,m[80]=12,m[81]=13,m[82]=14,m[83]=15,m[84]=16,m[86]=17,m[87]=18,m[89]=19,m[95]=20;var S=new Uint8Array(96);function _(e,t){var r=Math.floor(100*Math.random());switch(t){case n.NUCLEIC:switch(e){case"M":return[0,1][r%2];case"R":return[0,2][r%2];case"W":return[0,3][r%2];case"S":return[1,2][r%2];case"Y":return[1,3][r%2];case"K":return[2,3][r%2];case"V":return[0,1,2][r%3];case"H":return[0,1,3][r%3];case"D":return[0,2,3][r%3];case"B":return[1,2,3][r%3];case"N":default:return r%4}case n.PROTEIN:default:switch(e){case"B":return[2,11][r%2];case"Z":return[3,13][r%2];case"J":return[7,9][r%2];case"U":case"X":default:return r%20}}}function y(e){return e.map((function(e){return v[e]}))}S.fill(255),S[65]=0,S[67]=1,S[71]=2,S[84]=3,S[85]=3;var v=[0,1,2,2,3,0,4,5,4,5,5,2,0,2,4,0,0,5,3,3];function w(e,t){var r=null!=t?t:g(e),i=function(e,t){for(var r=new Uint8Array(e.length),i=t===n.PROTEIN?m:S,o=0,s=0,a=e.length;s<a;s++)255===(o=i[e.charCodeAt(s)])&&(o=_(e[s],t)),r[s]=o;return r}(e,r);return{rawSeq:e,encodedSeq:i,compressedSeq:r===n.PROTEIN?y(i):i,type:r}}function b(e,t,r){for(var n=[],i=r.gapchar,o=0,s=0;s<t.length;s++){var a=t[s];a<0?n.push(i.repeat(-a)):(n.push(e.substr(o,a)),o+=a)}return n.join("")}function O(e,t){var r=[0],n=0,i=t.slice();function o(e){(e^r[r.length-1])>=0?r[r.length-1]+=e:r.push(e)}for(var s=0;s<e.length;s++){var a=e[s];if(a<0)o(a);else for(;a>0&&n<i.length;){var u=i[n];u<0?a<-u?(i[n]+=a,o(-a),a=0):(o(u),a+=u,n++):a<u?(i[n]-=a,o(a),a=0):(o(u),a-=u,n++)}}return r}function q(e,t){for(var r=[],n=0,i=0,o=e[0],s=t[0];n<e.length||i<t.length;)o!==s?Math.sign(o)!==Math.sign(s)?o<0||void 0===s?(r.push(o),o=e[++n]):(r.push(s),s=t[++i]):o>0?o>s?(r.push(s),o-=s,s=t[++i]):(r.push(o),s-=o,o=e[++n]):(r.push(Math.min(o,s)),o=e[++n],s=t[++i]):(r.push(s),o=e[++n],s=t[++i]);return r}function A(e,t){for(var r=[],n=0,i=0,o=t[0],s=e[0],a=0;n<t.length||i<e.length;)if(a=r[r.length-1],o!==s){if(Math.sign(o)===Math.sign(s)){if(o>0){if(s>o)return;a>0?r[r.length-1]=a+s:r.push(s),o-=s,s=e[++i];continue}if(o>s){a>0?r[r.length-1]=a-o:r.push(Math.abs(o)),s-=o=t[++n];continue}return}if(!(s<0))return;a<0?r[r.length-1]=a+s:r.push(s),s=e[++i]}else a>0?r[r.length-1]=a+Math.abs(o):r.push(Math.abs(o)),o=t[++n],s=e[++i];return r}function E(e,r){if(Math.sign(e[e.length-1])===Math.sign(r[0])){var n=t([],e);return n[n.length-1]+=r[0],n.push.apply(n,r.slice(1)),n}return t(t([],e),r)}function C(e,t){void 0===t&&(t=0);for(var r=[],n=0,i=0,o=0;o<e.length;o++)((i=e[o])^n)>=0?n+=i:(r.push(n),n=i);return r.push(n),0===r[0]&&r.shift(),1&t||r.reverse(),r}function P(e){for(var t=0,r=0,n=0;n<e.length;n++)t+=(r=e[n])<0?-r:r;return t}function M(e){for(var t=[],r=0,n=0,i=0;i<e.length;i++)if((n=e[i])<0)for(var o=0;o<-n;o++)t.push(-1);else for(o=0;o<n;o++)t.push(r),r++;return t}function G(e,t,r,n){void 0===n&&(n=0);var i=e.rawSeq.length,o=t.rawSeq.length,s=0,a=[],u=[],c=0,h=0,f=0,p=0,d=new Uint8Array(Math.ceil((i+1)*o/2)),l=e.encodedSeq,g=t.encodedSeq,m=0,S=0,_=0,y=r.scoringMatrix[0],v=r.gapOP,w=1&n?0:-v/2,b=2&n?0:-v/2;a[0]=0,u[0]=-1/0;for(var O=1;O<=o;O++)a[O]=w,u[O]=-1/0;for(var q=1;q<=i;q++){c=w,h=-1/0,y=r.scoringMatrix[l[q-1]];for(O=1;O<=o;O++)_=0,m=a[O]+v,O===o&&(m+=b),m>=u[O]?u[O]=m:_+=1,S=c+v,q===i&&(S+=b),S>=h?h=S:_+=2,s=a[O-1]+y[g[O-1]],a[O-1]=c,s>=h?s>=u[O]?c=s:(c=u[O],_+=4):h>=u[O]?(c=h,_+=8):(c=u[O],_+=4),p=(f=q*o+O)%2,d[f>>>=1]+=p?_:_<<4;a[o]=c}var A=Math.max(s,h,u[u.length-1]),E=z(d,i,o,_),P=E[0],M=E[1];return{estrings:[C(P),C(M)],score:A}}function z(e,t,r,n){for(var i=t,o=r,s=t*r+r,a=s%2,u=[],c=[],h=(12&n)>>2,f=0;i>0&&o>0;)a=(s=i*r+o)%2,f=e[s>>>1],f=a?15&f:f>>>4,0===h?0===(f>>=2)&&(i--,o--,u.push(1),c.push(1)):2===h?(u.push(-1),c.push(1),f&=2,o--):(u.push(1),c.push(-1),f&=1,i--),h=f;return i>0?(u.push(i),c.push(-i)):o>0&&(u.push(-o),c.push(o)),[u,c]}var x,I=function(){function e(e,t,r,n){void 0===n&&(n=0),this.nbSeq=0,this.length=0,this.weight=0,this.length=e,this.m_bAllGaps=new i(e),this.m_uSortOrder=new Uint8Array(t*e),this.m_fcCounts=new Uint8Array(t*e),this.m_wCounts=new Float32Array(t*e),this.m_AAScores=new Float32Array(t*e),this.m_uResidueGroup=new Uint8Array(e),this.m_fOcc=new Float32Array(e),this.m_fcStartOcc=new Float32Array(e),this.m_fcEndOcc=new Float32Array(e),this.m_ScoreGapOpen=new Float32Array(e),this.m_ScoreGapClose=new Float32Array(e),this.nbSeq=r,this.alphaSize=t,this.weight=n}return e.prototype.getProxy=function(e){return void 0===e&&(e=0),new R(this,e)},e}(),R=function(){function e(e,t){void 0===t&&(t=0),this.idx=0,this.alphaSize=4,this.offset=0,this.prof=e,this.idx=t,this.alphaSize=e.alphaSize,this.offset=this.idx*this.alphaSize}return e.prototype.setProxy=function(e){this.idx=e,this.offset=e*this.alphaSize},Object.defineProperty(e.prototype,"m_bAllGaps",{get:function(){return this.prof.m_bAllGaps.get(this.idx)},set:function(e){e?this.prof.m_bAllGaps.set(this.idx):this.prof.m_bAllGaps.clear(this.idx)},enumerable:!1,configurable:!0}),Object.defineProperty(e.prototype,"m_uSortOrder",{get:function(){return this.prof.m_uSortOrder.subarray(this.offset,this.offset+this.alphaSize)},set:function(e){this.prof.m_uSortOrder.set(e,this.offset)},enumerable:!1,configurable:!0}),Object.defineProperty(e.prototype,"m_uResidueGroup",{get:function(){return this.prof.m_uResidueGroup[this.idx]},set:function(e){this.prof.m_uResidueGroup[this.idx]=e},enumerable:!1,configurable:!0}),Object.defineProperty(e.prototype,"m_fcCounts",{get:function(){return this.prof.m_fcCounts.subarray(this.offset,this.offset+this.alphaSize)},set:function(e){this.prof.m_fcCounts.set(e,this.offset)},enumerable:!1,configurable:!0}),Object.defineProperty(e.prototype,"m_wCounts",{get:function(){return this.prof.m_wCounts.subarray(this.offset,this.offset+this.alphaSize)},set:function(e){this.prof.m_wCounts.set(e,this.offset)},enumerable:!1,configurable:!0}),Object.defineProperty(e.prototype,"m_AAScores",{get:function(){return this.prof.m_AAScores.subarray(this.offset,this.offset+this.alphaSize)},set:function(e){this.prof.m_AAScores.set(e,this.offset)},enumerable:!1,configurable:!0}),Object.defineProperty(e.prototype,"m_fOcc",{get:function(){return this.prof.m_fOcc[this.idx]},set:function(e){this.prof.m_fOcc[this.idx]=e},enumerable:!1,configurable:!0}),Object.defineProperty(e.prototype,"m_fcStartOcc",{get:function(){return this.prof.m_fcStartOcc[this.idx]},set:function(e){this.prof.m_fcStartOcc[this.idx]=e},enumerable:!1,configurable:!0}),Object.defineProperty(e.prototype,"m_fcEndOcc",{get:function(){return this.prof.m_fcEndOcc[this.idx]},set:function(e){this.prof.m_fcEndOcc[this.idx]=e},enumerable:!1,configurable:!0}),Object.defineProperty(e.prototype,"m_ScoreGapOpen",{get:function(){return this.prof.m_ScoreGapOpen[this.idx]},set:function(e){this.prof.m_ScoreGapOpen[this.idx]=e},enumerable:!1,configurable:!0}),Object.defineProperty(e.prototype,"m_ScoreGapClose",{get:function(){return this.prof.m_ScoreGapClose[this.idx]},set:function(e){this.prof.m_ScoreGapClose[this.idx]=e},enumerable:!1,configurable:!0}),e}();function T(e,t,r,n){void 0===n&&(n=0);for(var i=e.encodedSeq,o=i.length,s=new I(o,r.abSize,1,t),a=t,u=r.gapOP,c=s.getProxy(),h=0,f=0;f<o;f++){c.setProxy(f),h=i[f],c.m_fcCounts[h]=1,c.m_wCounts[h]=a,c.m_uSortOrder[c.m_uResidueGroup]=h,c.m_uResidueGroup=1,c.m_fOcc=a,c.m_ScoreGapOpen=u/2,c.m_ScoreGapClose=u/2;for(var p=0;p<r.abSize;p++)c.m_AAScores[p]=a*r.scoringMatrix[p][h]}return 1^n&&(s.m_ScoreGapOpen[0]/=2,s.m_ScoreGapOpen[o-1]/=2),2^n&&(s.m_ScoreGapClose[0]/=2,s.m_ScoreGapClose[o-1]/=2),s}function N(e,t,r,n,i,o){void 0===o&&(o=0);for(var s=P(r),a=e.nbSeq+t.nbSeq,u=e.weight+t.weight,c=new I(s,i.abSize,a,u),h=i.gapOP,f=c.getProxy(),p=e.getProxy(),d=t.getProxy(),l=M(r),g=M(n),m=0,S=0,_=!1,y=!1,v=!1,w=!1,b=0,O=0;O<s;O++){if(f.setProxy(O),m=l[O],p.setProxy(m),-1===m?(_=!0,y=O===s-1||l[O+1]>0):(_=!1,y=!1),S=g[O],d.setProxy(S),-1===S?(v=!0,w=O==s-1||g[O+1]>0):(v=!1,w=!1),b=0,-1!==m&&-1!==S){f.m_fcStartOcc=p.m_fcStartOcc+d.m_fcStartOcc,f.m_fcEndOcc=p.m_fcEndOcc+d.m_fcEndOcc,f.m_fOcc=p.m_fOcc+d.m_fOcc,f.m_fcCounts.set(p.m_fcCounts),f.m_wCounts.set(p.m_wCounts);for(var q=0;q<t.alphaSize;q++){var A=t.m_fcCounts[q];A&&(f.m_fcCounts[q]+=A,f.m_wCounts[q]+=t.m_wCounts[q]),f.m_fcCounts[q]&&(f.m_uSortOrder[b]=q,b++),f.m_AAScores[q]=p.m_AAScores[q]+d.m_AAScores[q]}f.m_uResidueGroup=b}-1==m&&(f.m_fcStartOcc=d.m_fcStartOcc+(_?e.weight:0),f.m_fcEndOcc=d.m_fcEndOcc+(y?e.weight:0),f.m_fOcc=d.m_fOcc,f.m_fcCounts.set(d.m_fcCounts),f.m_wCounts.set(d.m_wCounts),f.m_uSortOrder.set(d.m_uSortOrder),f.m_AAScores.set(d.m_AAScores),f.m_uResidueGroup=d.m_uResidueGroup),-1==S&&(f.m_fcStartOcc=p.m_fcStartOcc+(v?t.weight:0),f.m_fcEndOcc=p.m_fcEndOcc+(w?t.weight:0),f.m_fOcc=p.m_fOcc,f.m_fcCounts.set(p.m_fcCounts),f.m_wCounts.set(p.m_wCounts),f.m_uSortOrder.set(p.m_uSortOrder),f.m_AAScores.set(p.m_AAScores),f.m_uResidueGroup=p.m_uResidueGroup),f.m_ScoreGapOpen=h/2*(1-f.m_fcStartOcc),f.m_ScoreGapClose=h/2*(1-f.m_fcEndOcc)}return 1^o&&(c.m_ScoreGapOpen[0]/=2,c.m_ScoreGapOpen[s-1]/=2),2^o&&(c.m_ScoreGapClose[0]/=2,c.m_ScoreGapClose[s-1]/=2),c}function U(e){return"seq"in e}function D(e,t){for(var r=1/0,n=0,i=0,o=e.length,s=0;s<o;s++)for(var a=s+1;a<o;a++)e[s][a]<r&&(r=e[s][a],n=s,i=a);return[{profile:null,childA:t[n],childB:t[i],distance:r,estring:[],msa:[],numSeq:[],type:x.NODE,depth:0,id:"",tabWeight:[],parent:-1,weight:0},n,i]}function j(e,t,r,n){for(var i=t>r?[t,r]:[r,t],o=i[0],s=i[1],a=[],u=e.length,c=0;c<u;c++)if(c!=t&&c!=r){var h=.1*(e[c][t]+e[c][r])/2+.9*Math.min(e[c][t],e[c][r]);e[c].push(h),a.push(h),e[c].splice(o,1),e[c].splice(s,1)}return a.push(0),e.push(a),e.splice(o,1),e.splice(s,1),n.splice(o,1),n.splice(s,1),e}function F(e,r){var o=function(e,n){var i=a[e.childA],s=a[e.childB],u={};if(U(i)||0!==i.msa.length||o(i,n),U(s)||0!==s.msa.length||o(s,n),U(i))if(U(s)){var h=G(i.seq,s.seq,r);u.score=h.score,e.numSeq=[i.numSeq[0],s.numSeq[0]],e.estring=h.estrings,e.profile=N(i.profile,s.profile,h.estrings[0],h.estrings[1],r),c.add("seq "+i.numSeq+" - seq "+s.numSeq)}else{s.tabWeight=n.filter((function(e,t){return s.numSeq.includes(t)}));var f=function(e,t,r,n){void 0===n&&(n=0);var i,o,s=t.seq,a=s.rawSeq.length,u=e.profile.length,c=0,h=[],f=[],p=new Uint8Array(Math.ceil((u+1)*(a+1)/2)),d=0,l=0,g=0,m=0,S=0,_=0,y=0,v=0,w=0,b=0,O=new Float32Array(r.abSize),q=r.gapOP,A=q/2,E=q/2,P=1&n?0:-q/4,M=2&n?0:-q/4,G=s.encodedSeq,x=e.profile;h[0]=0,f[0]=-1/0;for(var I=1;I<=a;I++)h[I]=A+E+P,f[I]=-1/0;for(var R=1;R<=u;R++){for(w=x.m_ScoreGapOpen[R-1],b=x.m_ScoreGapClose[R-1],O=x.m_AAScores.subarray((R-1)*r.abSize,R*r.abSize),_=P,m=-1/0,f[0]=x.m_ScoreGapOpen[0],S=M,I=1;I<=a;I++)g=0,d=h[I]+w,I===a&&(d-=w/2),d>=f[I]?f[I]=d:g+=1,l=_+A,R===u&&(l+=M),l>=(S=m)?m=l:g+=2,y=f[I]+b,v=S+E,c=h[I-1]+O[G[I-1]],I===a&&(v+=M),h[I-1]=_,c>=v?c>=y?_=c:(_=y,g+=4):v>=y?(_=v,g+=8):(_=y,g+=4),o=(i=R*a+I)%2,p[i>>>=1]+=o?g:g<<4;h[a]=_}var T=Math.max(c,m,f[f.length-2]),N=z(p,u,a,g),U=N[0],D=N[1];return{estrings:[C(U),C(D)],score:T}}(s,i,r);u.score=f.score,e.numSeq=t([i.numSeq[0]],s.numSeq),e.estring=t([f.estrings[1]],s.estring.map((function(e){return O(f.estrings[0],e)}))),e.profile=N(i.profile,s.profile,f.estrings[1],f.estrings[0],r),c.add("seq "+i.numSeq+" - MSA "+s.numSeq)}else if(!U(s)){s.tabWeight=n.filter((function(e,t){return s.numSeq.includes(t)})),i.tabWeight=n.filter((function(e,t){return i.numSeq.includes(t)}));var p=function(e,t,r,n){void 0===n&&(n=0);var i,o,s=e.profile.length,a=t.profile.length,u=new Uint8Array(Math.ceil((s+1)*(a+1)/2)),c=[],h=[],f=0,p=0,d=0,l=0,g=0,m=0,S=0,_=0,y=0,v=0,w=0,b=0,O=new Float32Array(r.abSize),q=new Uint8Array(r.abSize),A=0,E=0,P=0,M=0,G=0,x=0;r.gapOP;var I=1&n?1:2,R=t.profile,T=e.profile;for(c[0]=0,h[0]=-1/0,G=1;G<=a;G++)c[G]=R.m_ScoreGapOpen[0]+R.m_ScoreGapClose[G-1]/I,h[G]=-1/0;var N=T.m_ScoreGapOpen[0];for(M=1;M<=s;M++){for(O=T.m_AAScores.subarray(r.abSize*(M-1),r.abSize*M),w=T.m_ScoreGapOpen[M-1],_=N+(b=T.m_ScoreGapClose[M-1])/I,m=-1/0,h[0]=T.m_ScoreGapOpen[0],G=1;G<=a;G++){for(g=0,p=c[G]+w,G!==a||2&n||(p-=w/2),p>=h[G]?h[G]=p:g+=1,d=_+R.m_ScoreGapOpen[G-1],l=S=m,M!==s||2&n||(d-=R.m_ScoreGapOpen[G-1]/2),d>=l?m=d:g+=2,f=c[G-1],P=R.m_uResidueGroup[G-1],x=0,E=(G-1)*r.abSize,q=R.m_uSortOrder.subarray(E,E+P);x<P;)A=q[x],f+=R.m_wCounts[E+A]*O[A],x++;y=h[G]+b,v=S+R.m_ScoreGapClose[G-1],1!==G||2&n||(y-=b/2),G!==a||2&n||(y-=b/2,v-=R.m_ScoreGapClose[a-1]/2),c[G-1]=_,f>=v?f>=y?_=f:(_=y,g+=4):v>=y?(_=v,g+=8):(_=y,g+=4),o=(i=M*a+G)%2,u[i>>>=1]+=o?g:g<<4}c[a]=_}var U=Math.max(f,m,h[G-1]),D=z(u,s,a,g),j=D[0],F=D[1];return{estrings:[C(j),C(F)],score:U}}(i,s,r);u.score=p.score,e.estring=t(t([],i.estring.map((function(e){return O(p.estrings[0],e)}))),s.estring.map((function(e){return O(p.estrings[1],e)}))),e.profile=N(i.profile,s.profile,p.estrings[0],p.estrings[1],r),c.add("MSA "+i.numSeq+" - MSA "+s.numSeq)}return e.numSeq=t(t([],i.numSeq),s.numSeq),u.score};c.add("Start Progressive Alignment");var s=function(e){var t=e[0].type===n.PROTEIN,r=t?6:4,o=Math.pow(r,6),s=e.map((function(e){for(var n=new i(o),s=t?e.compressedSeq:e.encodedSeq,a=0,u=0;u<=5;u++)a+=s[u]*Math.pow(r,u);n.set(a);for(var c=Math.pow(r,5),h=(u=6,s.length);u<h;u++)a-=s[u-6],a/=6,a+=s[u]*c,n.set(a);return n}));c.add("K-mer bitset computation");for(var a,u,h,f=e.length,p=e.map((function(){return[]})),d=0;d<f;d++){p[d][d]=0,a=s[d],u=e[d].compressedSeq.length;for(var l=d+1;l<f;l++)h=1-a.getIntersectionSize(s[l])/u,p[l][d]=p[d][l]=h}return c.add("K-mer distance computation"),p}(e);c.add("K-mer distance matrix");var a=function(e,t){for(var r,n,i,o,s=[],a=e.length,u=[],c=0;c<a;c++)s[c]={type:x.LEAF,seq:t[c],profile:null,childA:c,childB:c,distance:0,numSeq:[c],msa:[],id:c.toString(),weight:0,parent:-1,depth:0},u[c]=c;c=0;for(var h=a-1;c<h;c++){n=(r=D(e,u))[0],i=r[1],o=r[2];var f=s[n.childA],p=s[n.childB];n.id=f.id<p.id?"|"+f.id+","+p.id+"|":"|"+p.id+","+f.id+"|",n.depth=Math.max(f.depth,p.depth)+1,p.parent=f.parent=s.length,u.push(a+c),s.push(n),e=j(e,i,o,u)}return s[s.length-1].type=x.ROOT,s}(s,e),u=a[a.length-1];c.add("Build Tree");var h=function(e){var t=[],r=0;!function t(n){var i=0,o=0;switch(n.type){case x.ROOT:n.weight=0,t(e[n.childA]),t(e[n.childB]);break;case x.NODE:i=e[n.parent].distance-n.distance,o=n.id.split(",").length,n.weight=e[n.parent].weight+i/o,t(e[n.childA]),t(e[n.childB]);break;case x.LEAF:i=e[n.parent].distance,n.weight=e[n.parent].weight+i,r+=n.weight}}(e[e.length-1]);for(var n=0;e[n].type==x.LEAF;)t[n]=e[n].weight/=r,n++;return t}(a);c.add("Compute Weights - Start MSA"),function(e,t,r){for(var n=0,i=e[n];U(i)&&n<e.length;)i.profile=T(i.seq,i.weight,t,r),i=e[++n]}(a,r),o(u,h);var f,p,d,l=(f=u.estring.slice(),p=u.numSeq,d=p.slice(),p.forEach((function(e,t){return d[e]=t})),d.map((function(e){return f[e]})));return c.add("End MSA computation"),l}!function(e){e[e.LEAF=0]="LEAF",e[e.NODE=1]="NODE",e[e.ROOT=2]="ROOT"}(x||(x={}));var B=function(){function e(e){this.storeSize=e,this.store=new Array(e),this.head=0,this.tail=0}return Object.defineProperty(e.prototype,"size",{get:function(){return(this.tail-this.head+this.storeSize)%this.storeSize},enumerable:!1,configurable:!0}),Object.defineProperty(e.prototype,"isEmpty",{get:function(){return this.head===this.tail},enumerable:!1,configurable:!0}),e.prototype.popHead=function(){if(this.isEmpty)return null;var e=this.store[this.head];return this.head++,this.head>=this.storeSize&&(this.head=0),e},e.prototype.popTail=function(){return this.isEmpty?null:(this.tail--,this.tail<0&&(this.tail=this.storeSize-1),this.store[this.tail])},e.prototype.pushHead=function(e){this.head--,this.head<0&&(this.head=this.storeSize-1),this.store[this.head]=e,this.head===this.tail&&(this.tail=(this.tail-1+this.storeSize)%this.storeSize)},e.prototype.pushTail=function(e){this.store[this.tail]=e,this.tail++,this.tail>=this.storeSize&&(this.tail=0),this.head===this.tail&&(this.head=(this.head-1+this.storeSize)%this.storeSize)},e.prototype.getHead=function(e){void 0===e&&(e=0);var t=0===e?this.head:(this.head+e+this.storeSize)%this.storeSize;return this.store[t]},e.prototype.getTail=function(e){void 0===e&&(e=0);var t=(this.tail-1-e+this.storeSize)%this.storeSize;return this.store[t]},e}();function L(e,t,r){for(var n=new Map,i=new B(r),o=0|t,s=new Uint16Array(e.encodedSeq.length-o+1),a=0,u=0,c=0;c<o;c++)a|=e.encodedSeq[o-c-1]<<2*c;s[u++]=a;c=o;for(var h=e.encodedSeq.length;c<h;c++)a=(a<<2)+e.encodedSeq[c],s[u++]=a;var f=r-o,p=e.rawSeq.length*(2/(r+1))*2|0,d={kmer:new Uint16Array(p),kmerPos:new Uint16Array(p),winPos:new Uint16Array(p),winPosEnd:new Uint16Array(p),count:0},l=NaN;for(c=0;c<s.length;c++){for(var g=s[c];!i.isEmpty&&s[i.getTail()]>=g;)i.popTail();if(i.pushTail(c),!(c<f)){for(;i.getHead()<=c-1-f;)i.popHead();var m=i.getHead();if(m===l)d.winPosEnd[d.count-1]=c+o;else{var S=s[m],_=d.count;if(d.count++,d.kmer[_]=S,d.kmerPos[_]=m,d.winPos[_]=c-f,d.winPosEnd[_]=c+o,l=m,n.has(S))n.get(S).push(_);else n.set(S,[_])}}}return[n,d,s]}function H(e,t,r,n,i){var o={},s=null!=n?n:L(e,8,16),a=s[0],u=s[1],h=s[2],f=null!=i?i:L(t,8,16),p=f[0],d=f[1],l=f[2];c.add("Extract Minimizers"),o["Nb Minimizers"]={a:u.count,b:d.count},o["Dupl. Minimizers"]={a:u.count-a.size,b:d.count-p.size};for(var g=[],m=new Map,S=0;S<u.count;S++){var _=u.kmer[S];if(p.has(_))for(var y=p.get(_),v=e.rawSeq.substring(u.winPos[S],u.winPosEnd[S]),w=0;w<y.length;w++){var b=y[w],O=t.rawSeq.substring(d.winPos[b],d.winPosEnd[b]),q=Math.min(v.length,O.length);if(q==v.length){if(0!==O.indexOf(v))continue}else if(0!==v.indexOf(O))continue;var A=u.winPos[S]-d.winPos[b],M={diagId:A,begin:u.winPos[S],end:u.winPos[S]+q};if(g.push(M),m.has(A))m.get(A).push(M);else m.set(A,[M])}}c.add("Filter Minimizers"),o["Nb common"]={all:g.length};var z=[];m.forEach((function(e){var t=e[0];z.push(t),e.forEach((function(e){e.begin<=t.end?t.end=e.end:(t=e,z.push(t))}))})),c.add("Merge Minimizers"),o["Nb diagonals"]={all:z.length},z.sort((function(e,t){var r=e.begin-t.begin;return 0===r?e.diagId-t.diagId:r})),c.add("Sort Minimizers");var x=z[0],I=[x],R=[];for(S=1;S<z.length;S++){var T=z[S];if(!(T.begin<x.end))if(T.diagId!==x.diagId){for(U=(N=x.end)-x.diagId;e.encodedSeq[N]===t.encodedSeq[U]&&N<T.begin;)N++,U++;for(x.end=N,U=(N=T.begin)-T.diagId;e.encodedSeq[N]===t.encodedSeq[U]&&N>x.end;)N--,U--;T.begin=N+1,R.push({endDiagId:T.diagId,beginDiagId:x.diagId,begin:x.end,end:T.begin,size:T.begin-x.end}),I.push(T),x=T}else{for(var N,U=(N=x.end)-T.diagId,D=0;N<T.begin&&D<2;)U+=8,D=W(h[N+=8],l[U])<=2?0:D+1;for(N>x.end&&(x.end=N-2-8),U=(N=T.begin-8)-T.diagId,D=0;N>x.end&&D<2;)U-=8,D=W(h[N-=8],l[U])<=2?0:D+1;if(N<T.begin-16&&(T.begin=N+24+2),T.begin-x.end<=24){x.end=T.end;continue}R.push({beginDiagId:x.diagId,endDiagId:T.diagId,begin:x.end,end:T.begin,size:T.begin-x.end}),I.push(T),x=T}}c.add("Diagonals extension");var j=0;I.forEach((function(e){j+=e.end-e.begin})),j/=e.encodedSeq.length,o.Coverage={all:j},o["Extended Diagonals"]={all:I.length},console.table(o);var F=I[0],B=[],H=[];for(S=0;S<I.length;S++){F=I[S],B.push(F.end-F.begin),H.push(F.end-F.begin);var k=R[S];if(k){if(k.begin===k.end){B.push(-Math.abs(k.endDiagId-k.beginDiagId)),H.push(k.end-k.begin);continue}if(k.begin-k.beginDiagId==k.end-k.endDiagId){H.push(-Math.abs(k.endDiagId-k.beginDiagId)),H.push(k.end-k.begin);continue}var K=G({rawSeq:e.rawSeq.substring(k.begin,k.end),type:e.type,compressedSeq:new Uint8Array(0),encodedSeq:e.encodedSeq.slice(k.begin,k.end)},{rawSeq:t.rawSeq.substring(k.begin-k.beginDiagId,k.end-k.endDiagId),type:t.type,compressedSeq:new Uint8Array(0),encodedSeq:t.encodedSeq.slice(k.begin-k.beginDiagId,k.end-k.endDiagId)},r,3);B.push.apply(B,K.estrings[0]),H.push.apply(H,K.estrings[1])}}V(B,e.encodedSeq.length),V(H,t.encodedSeq.length);var Y=C(B,1),Q=C(H,1),$=P(Y),J=P(Q);return $>J?Q=E(Q,[J-$]):$<J&&(Y=E(Y,[$-J])),c.add("Fill between diagonals"),[Y,Q]}function V(e,t){var r=t-function(e){for(var t=0,r=0,n=0;n<e.length;n++)(r=e[n])<0||(t+=r);return t}(e);r>0?e.push(r):r<0&&(e[e.length-1]+=r)}function W(e,t){var n=e^t;return r(n=1431655765&n|n>>1&1431655765)}var k={sequenceType:"auto",gapchar:"-",alignmentMethod:"auto",debug:!1};return new(function(){function t(){this.sequences=[],this.typeSeq=n.UNSET,this.config=e({},k)}return t.prototype.addSequence=function(e){var t=this;if(Array.isArray(e))e.forEach((function(e){t.addSequence(e)}));else{if("string"!=typeof e)throw new TypeError("String type expected for sequences to add.");if(e=e.toUpperCase(),"auto"===this.config.sequenceType){var r=g(e);if(this.typeSeq===n.UNSET)this.typeSeq=r;else if(this.typeSeq!==r)throw new Error("All sequences must be of same type.")}this.sequences.push(w(e))}},t.prototype.reset=function(){this.sequences=[],this.typeSeq=n.UNSET,this.setDefaultConfiguration()},t.prototype.setDefaultConfiguration=function(){this.config=e({},k)},t.prototype.setUserConfiguration=function(e){var t;function r(e,t){return t.some((function(t){return t==e}))}e.method&&r(e.method,["complete","diag"])&&(this.config.alignmentMethod=e.method),e.type&&r(e.type,["amino","nucleic"])&&(this.config.sequenceType=e.type,this.typeSeq="amino"===e.type?n.PROTEIN:n.NUCLEIC),1===(null===(t=e.gapchar)||void 0===t?void 0:t.length)&&(this.config.gapchar=e.gapchar),void 0!==e.debug&&(this.config.debug=!!e.debug)},t.prototype.align=function(e,t){var r=this;return c.start(),new Promise((function(i,o){if(!Array.isArray(e)||e.some((function(e){return"string"!=typeof e})))return o("Array of sequences expected");if(e.length<2)return o("At least 2 sequences are required.");r.reset(),t&&r.setUserConfiguration(t),r.addSequence(e),c.add("Prepared sequences");var u=function(e,t){var r,i,o,u=e===n.PROTEIN?s:a,c=null!==(r=null==t?void 0:t.matrix)&&void 0!==r?r:u.matrix,h=(null==t?void 0:t.gapextend)?2*-t.gapextend:u.center;return{type:e,scoringMatrix:(o=h,c.map((function(e){return e.map((function(e){return e+o}))}))),gapOP:null!==(i=null==t?void 0:t.gapopen)&&void 0!==i?i:u.gapOP,abSize:e===n.PROTEIN?20:4}}(r.typeSeq,t);c.add("Get sequences type");var h,f=!1;(f="auto"===r.config.alignmentMethod?r.sequences.some((function(e){return e.rawSeq.length>1600})):"diag"===r.config.alignmentMethod,2==r.sequences.length)?h=f?H(r.sequences[0],r.sequences[1],u):G(r.sequences[0],r.sequences[1],u).estrings:h=f?function(e,t){for(var r=e.map((function(e){return L(e,8,16)})),n=[],i=[],o=[],s=0;s<r.length;s++)if(0!==s){var a=H(e[0],e[s],t,r[0],r[s]);i=0===i.length?a[0]:q(i,a[0]),o[s]=a[0],n[s]=a[1]}for(s=0;s<r.length;s++)if(0!==s){var u=A(i,o[s]);if(!u)throw console.error("An error occured while computing the center diff. for sequence #"+s,i,o[s]),console.dir(n),console.dir(o),new RangeError;n[s]=O(u,n[s])}else n[s]=i;return n}(r.sequences,u):F(r.sequences,u);return c.summary(),i(function(e,t,r){for(var n,i=null!==(n=null==r?void 0:r.gapchar)&&void 0!==n?n:"-",o=[],s=0;s<e.length;s++)o.push(b(e[s].rawSeq,t[s],{gapchar:i}));return o}(r.sequences,h,{gapchar:r.config.gapchar}))}))},t}())}));
//# sourceMappingURL=biomsa.js.map
