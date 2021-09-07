var e=function(){return(e=Object.assign||function(e){for(var t,r=1,n=arguments.length;r<n;r++)for(var i in t=arguments[r])Object.prototype.hasOwnProperty.call(t,i)&&(e[i]=t[i]);return e}).apply(this,arguments)};function t(e,t){for(var r=0,n=t.length,i=e.length;r<n;r++,i++)e[i]=t[r];return e}function r(e){return 16843009*((e=(858993459&(e-=e>>>1&1431655765))+(e>>>2&858993459))+(e>>>4)&252645135)>>>24}var n,i=function(){function e(e,t){this.length=e,this._words=new Uint32Array(e+32>>>5),!0===t&&this.setAll()}return e.prototype.get=function(e){return 0!=(this._words[e>>>5]&1<<e)},e.prototype.set=function(e){this._words[e>>>5]|=1<<e},e.prototype.clear=function(e){this._words[e>>>5]&=~(1<<e)},e.prototype.flip=function(e){this._words[e>>>5]^=1<<e},e.prototype._assignRange=function(e,t,r){if(!(t<e)){for(var n=this._words,i=!0===r?4294967295:0,o=e>>>5,s=t>>>5,a=o;a<s;++a)n[a]=i;var u=o<<5,c=s<<5;if(!0===r)if(t-e<32)for(var h=e,f=t+1;h<f;++h)n[h>>>5]|=1<<h;else{for(h=e,f=u;h<f;++h)n[h>>>5]|=1<<h;for(h=c,f=t+1;h<f;++h)n[h>>>5]|=1<<h}else if(t-e<32)for(h=e,f=t+1;h<f;++h)n[h>>>5]&=~(1<<h);else{for(h=e,f=u;h<f;++h)n[h>>>5]&=~(1<<h);for(h=c,f=t+1;h<f;++h)n[h>>>5]&=~(1<<h)}return this}},e.prototype.setRange=function(e,t){return this._assignRange(e,t,!0)},e.prototype.clearRange=function(e,t){return this._assignRange(e,t,!1)},e.prototype.setBits=function(){for(var e=[],t=0;t<arguments.length;t++)e[t]=arguments[t];for(var r=this._words,n=e.length,i=0;i<n;++i){var o=e[i];r[o>>>5]|=1<<o}return this},e.prototype.clearBits=function(){for(var e=[],t=0;t<arguments.length;t++)e[t]=arguments[t];for(var r=this._words,n=e.length,i=0;i<n;++i){var o=e[i];r[o>>>5]&=~(1<<o)}return this},e.prototype.setAll=function(){return this._assignRange(0,this.length-1,!0)},e.prototype.clearAll=function(){return this._assignRange(0,this.length-1,!1)},e.prototype.flipAll=function(){for(var e=this._words.length,t=this._words,r=32-this.length%32,n=0;n<e-1;++n)t[n]=~t[n];return t[e-1]=~(t[e-1]<<r)>>>r,this},e.prototype._isRangeValue=function(e,t,r){if(!(t<e)){for(var n=this._words,i=!0===r?4294967295:0,o=e>>>5,s=t>>>5,a=o;a<s;++a)if(n[a]!==i)return!1;if(t-e<32){for(var u=e,c=t+1;u<c;++u)if(!!(n[u>>>5]&1<<u)!==r)return!1}else{var h=s<<5;for(u=e,c=o<<5<<5;u<c;++u)if(!!(n[u>>>5]&1<<u)!==r)return!1;for(u=h,c=t+1;u<c;++u)if(!!(n[u>>>5]&1<<u)!==r)return!1}return!0}},e.prototype.isRangeSet=function(e,t){return this._isRangeValue(e,t,!0)},e.prototype.isRangeClear=function(e,t){return this._isRangeValue(e,t,!1)},e.prototype.isAllSet=function(){return this._isRangeValue(0,this.length-1,!0)},e.prototype.isAllClear=function(){return this._isRangeValue(0,this.length-1,!1)},e.prototype.isSet=function(){for(var e=[],t=0;t<arguments.length;t++)e[t]=arguments[t];for(var r=this._words,n=e.length,i=0;i<n;++i){var o=e[i];if(0==(r[o>>>5]&1<<o))return!1}return!0},e.prototype.isClear=function(){for(var e=[],t=0;t<arguments.length;t++)e[t]=arguments[t];for(var r=this._words,n=e.length,i=0;i<n;++i){var o=e[i];if(0!=(r[o>>>5]&1<<o))return!1}return!0},e.prototype.isEqualTo=function(e){for(var t=this._words,r=e._words,n=Math.min(t.length,r.length),i=0;i<n;++i)if(t[i]!==r[i])return!1;return!0},e.prototype.getSize=function(){for(var e=this._words.length,t=this._words,n=0,i=0;i<e;++i)n+=r(t[i]);return n},e.prototype.difference=function(e){for(var t=this._words,r=e._words,n=Math.min(t.length,r.length),i=0;i<n;++i)t[i]=t[i]&~r[i];for(i=t.length;i<n;++i)t[i]=0;return this},e.prototype.union=function(e){for(var t=this._words,r=e._words,n=Math.min(t.length,r.length),i=0;i<n;++i)t[i]|=r[i];for(i=t.length;i<n;++i)t[i]=0;return this},e.prototype.intersection=function(e){for(var t=this._words,r=e._words,n=Math.min(t.length,r.length),i=0;i<n;++i)t[i]&=r[i];for(i=t.length;i<n;++i)t[i]=0;return this},e.prototype.intersects=function(e){for(var t=this._words,r=e._words,n=Math.min(t.length,r.length),i=0;i<n;++i)if(0!=(t[i]&r[i]))return!0;return!1},e.prototype.getIntersectionSize=function(e){for(var t=this._words,n=e._words,i=Math.min(t.length,n.length),o=0,s=0;s<i;++s)o+=r(t[s]&n[s]);return o},e.prototype.makeIntersection=function(t){var r=this._words,n=t._words,i=Math.min(r.length,n.length),o=new Uint32Array(i),s=Object.create(e.prototype);s._words=o,s.length=Math.min(this.length,t.length);for(var a=0;a<i;++a)o[a]=r[a]&n[a];return s},e.prototype.forEach=function(e){for(var t=this._words.length,n=this._words,i=0,o=0;o<t;++o)for(var s=n[o];0!==s;){var a=s&-s;e((o<<5)+r(a-1),i),s^=a,++i}},e.prototype.toArray=function(){for(var e=this._words,t=new Array(this.getSize()),n=this._words.length,i=0,o=0;o<n;++o)for(var s=e[o];0!==s;){var a=s&-s;t[i++]=(o<<5)+r(a-1),s^=a}return t},e.prototype.toString=function(){return"{"+this.toArray().join(",")+"}"},e.prototype.toSeleString=function(){var e=this.toArray().join(",");return e?"@"+e:"NONE"},e.prototype.clone=function(){var t=Object.create(e.prototype);return t.length=this.length,t._words=new Uint32Array(this._words),t},e}();!function(e){e[e.PROTEIN=0]="PROTEIN",e[e.NUCLEIC=1]="NUCLEIC",e[e.UNSET=2]="UNSET"}(n||(n={}));var o,s={matrix:[[58,23,-12,-7,-44,10,-23,-14,-14,-27,-17,-8,1,-9,-22,23,15,5,-74,-45,0],[23,224,-67,-63,-50,-30,-29,1,-56,-41,-6,-33,-44,-53,-43,15,2,18,-93,-6,0],[-12,-67,111,59,-104,-4,4,-84,6,-88,-65,48,-13,18,-29,5,-7,-63,-105,-73,0],[-7,-63,59,85,-83,-17,-1,-63,25,-60,-47,15,-12,40,-8,1,-7,-47,-108,-51,0],[-44,-50,-104,-83,144,-93,4,12,-74,36,30,-64,-67,-56,-65,-43,-41,-3,63,104,0],[10,-30,-4,-17,-93,140,-32,-95,-27,-91,-75,4,-36,-29,-32,5,-26,-68,-80,-79,0],[-23,-29,4,-1,4,-32,137,-50,6,-37,-42,21,-23,27,19,-4,-12,-44,-13,48,0],[-14,1,-84,-63,12,-95,-50,86,-53,53,47,-62,-60,-47,-55,-43,-8,69,-27,-24,0],[-14,-56,6,25,-74,-27,6,-53,75,-48,-30,13,-12,34,68,-3,-4,-44,-71,-49,0],[-27,-41,-88,-60,36,-91,-37,53,-48,88,62,-63,-48,-36,-48,-47,-25,36,-11,-4,0],[-17,-6,-65,-47,30,-75,-42,47,-30,62,103,-45,-54,-21,-31,-35,-9,31,-46,-20,0],[-8,-33,48,15,-64,4,21,-62,13,-63,-45,89,-25,12,2,22,10,-51,-79,-29,0],[1,-44,-13,-12,-67,-36,-23,-60,-12,-48,-54,-25,160,-6,-20,5,-12,-42,-76,-83,0],[-9,-53,18,40,-56,-29,27,-47,34,-36,-21,12,-6,75,34,1,-4,-37,-92,-48,0],[-22,-43,-29,-8,-65,-32,19,-55,68,-48,-31,2,-20,34,113,-10,-14,-49,-58,-39,0],[23,15,5,1,-43,5,-4,-43,-3,-47,-35,22,5,1,-10,53,32,-28,-62,-31,0],[15,2,-7,-7,-41,-26,-12,-8,-4,-25,-9,10,-12,-4,-14,32,68,0,-87,-40,0],[5,18,-63,-47,-3,-68,-44,69,-44,36,31,-51,-42,-37,-49,-28,0,74,-61,-32,0],[-74,-93,-105,-108,63,-80,-13,-27,-71,-11,-46,-79,-76,-92,-58,-62,-87,-61,289,81,0],[-45,-6,-73,-51,104,-79,48,-24,-49,-4,-20,-29,-83,-48,-39,-31,-40,-32,81,162,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]],center:22,gapOP:-300},a={matrix:[[91,-114,-31,-123],[-114,100,-125,-31],[-31,-125,100,-114],[-123,-31,-114,91]],gapOP:-400,gapEP:-60,center:120};function u(e){o.has(e)?u(e+"_"):o.set(e,performance.now())}var c={start:function(){o=new Map,u("START")},add:u,summary:function(){var e={start:0};u("END");var t=0;o.forEach((function(r,n){"START"!==n?(e[n]=r-t,t=r):t=r})),e.TOTAL=o.get("END")-o.get("START"),console.table(e)}},h=/^[ACGTRNDEQHILKMFPSWYV]+$/,f=/^[ABCDGHKMNRSTUVWY]+$/i,p=/^[ABCGTRNDEQHIJLKMFPSWYUVXZ]+$/i,d=/[^ACGU]/i,l=/[^ACGT]/i;function g(e){if(!d.test(e)||!l.test(e))return n.NUCLEIC;var t=f.test(e),r=h.test(e);if(r&&!t)return n.PROTEIN;if(r&&t)return function(e){for(var t=100,r=0,i=Math.min(e.length,t),o=0;o<i;o++)switch(e[o]){case"A":case"T":case"U":case"G":case"C":case"N":r++}return r/i>Math.SQRT1_2?n.NUCLEIC:n.PROTEIN}(e);if(p.test(e))return n.PROTEIN;throw new Error("Unrecognized sequence type: "+e)}var m=new Uint8Array(96);m.fill(255),m[65]=0,m[67]=1,m[68]=2,m[69]=3,m[70]=4,m[71]=5,m[72]=6,m[73]=7,m[75]=8,m[76]=9,m[77]=10,m[78]=11,m[80]=12,m[81]=13,m[82]=14,m[83]=15,m[84]=16,m[86]=17,m[87]=18,m[89]=19,m[95]=20;var S=new Uint8Array(96);function v(e,t){var r=Math.floor(100*Math.random());switch(t){case n.NUCLEIC:switch(e){case"M":return[0,1][r%2];case"R":return[0,2][r%2];case"W":return[0,3][r%2];case"S":return[1,2][r%2];case"Y":return[1,3][r%2];case"K":return[2,3][r%2];case"V":return[0,1,2][r%3];case"H":return[0,1,3][r%3];case"D":return[0,2,3][r%3];case"B":return[1,2,3][r%3];case"N":default:return r%4}case n.PROTEIN:default:switch(e){case"B":return[2,11][r%2];case"Z":return[3,13][r%2];case"J":return[7,9][r%2];case"U":case"X":default:return r%20}}}function _(e){return e.map((function(e){return y[e]}))}S.fill(255),S[65]=0,S[67]=1,S[71]=2,S[84]=3,S[85]=3;var y=[0,1,2,2,3,0,4,5,4,5,5,2,0,2,4,0,0,5,3,3];function w(e,t){var r=null!=t?t:g(e),i=function(e,t){for(var r=new Uint8Array(e.length),i=t===n.PROTEIN?m:S,o=0,s=0,a=e.length;s<a;s++)255===(o=i[e.charCodeAt(s)])&&(o=v(e[s],t)),r[s]=o;return r}(e,r);return{rawSeq:e,encodedSeq:i,compressedSeq:r===n.PROTEIN?_(i):i,type:r}}function b(e,t,r){for(var n=[],i=r.gapchar,o=0,s=0;s<t.length;s++){var a=t[s];a<0?n.push(i.repeat(-a)):(n.push(e.substr(o,a)),o+=a)}return n.join("")}function O(e,t){var r=[0],n=0,i=t.slice();function o(e){(e^r[r.length-1])>=0?r[r.length-1]+=e:r.push(e)}for(var s=0;s<e.length;s++){var a=e[s];if(a<0)o(a);else for(;a>0&&n<i.length;){var u=i[n];u<0?a<-u?(i[n]+=a,o(-a),a=0):(o(u),a+=u,n++):a<u?(i[n]-=a,o(a),a=0):(o(u),a-=u,n++)}}return r}function A(e,t){for(var r=[],n=0,i=0,o=e[0],s=t[0];n<e.length||i<t.length;)o!==s?Math.sign(o)!==Math.sign(s)?o<0||void 0===s?(r.push(o),o=e[++n]):(r.push(s),s=t[++i]):o>0?o>s?(r.push(s),o-=s,s=t[++i]):(r.push(o),s-=o,o=e[++n]):(r.push(Math.min(o,s)),o=e[++n],s=t[++i]):(r.push(s),o=e[++n],s=t[++i]);return r}function q(e,t){for(var r=[],n=0,i=0,o=t[0],s=e[0],a=0;n<t.length||i<e.length;)if(a=r[r.length-1],o!==s){if(Math.sign(o)===Math.sign(s)){if(o>0){if(s>o)return;a>0?r[r.length-1]=a+s:r.push(s),o-=s,s=e[++i];continue}if(o>s){a>0?r[r.length-1]=a-o:r.push(Math.abs(o)),s-=o=t[++n];continue}return}if(!(s<0))return;a<0?r[r.length-1]=a+s:r.push(s),s=e[++i]}else a>0?r[r.length-1]=a+Math.abs(o):r.push(Math.abs(o)),o=t[++n],s=e[++i];return r}function E(e,r){if(Math.sign(e[e.length-1])===Math.sign(r[0])){var n=t([],e);return n[n.length-1]+=r[0],n.push.apply(n,r.slice(1)),n}return t(t([],e),r)}function C(e,t){void 0===t&&(t=0);for(var r=[],n=0,i=0,o=0;o<e.length;o++)((i=e[o])^n)>=0?n+=i:(r.push(n),n=i);return r.push(n),0===r[0]&&r.shift(),1&t||r.reverse(),r}function P(e){for(var t=0,r=0,n=0;n<e.length;n++)t+=(r=e[n])<0?-r:r;return t}function M(e){for(var t=[],r=0,n=0,i=0;i<e.length;i++)if((n=e[i])<0)for(var o=0;o<-n;o++)t.push(-1);else for(o=0;o<n;o++)t.push(r),r++;return t}function I(e,t,r,n){void 0===n&&(n=0);var i=e.encodedSeq,o=t.encodedSeq,s=i.length,a=o.length,u=0,c=[],h=[],f=0,p=0,d=0,l=0,g=new Uint8Array(Math.ceil((s+1)*a/2)),m=0,S=0,v=0,_=r.scoringMatrix[0],y=r.gapOP,w=1&n?0:-y/2,b=2&n?0:-y/2;c[0]=0,h[0]=-1/0;for(var O=1;O<=a;O++)c[O]=y+w,h[O]=-1/0;for(var A=1;A<=s;A++){f=w,p=-1/0,_=r.scoringMatrix[i[A-1]];for(O=1;O<=a;O++)v=0,m=c[O]+y,O===a&&(m+=b),m>=h[O]?h[O]=m:v+=1,S=f+y,A===s&&(S+=b),S>=p?p=S:v+=2,u=c[O-1]+_[o[O-1]],c[O-1]=f,u>=p?u>=h[O]?f=u:(f=h[O],v+=4):p>=h[O]?(f=p,v+=8):(f=h[O],v+=4),l=(d=A*a+O)%2,g[d>>>=1]+=l?v:v<<4;c[a]=f}var q=Math.max(u,p,h[h.length-1]),E=G(g,s,a,v),P=E[0],M=E[1];return{estrings:[C(P),C(M)],score:q}}function G(e,t,r,n){for(var i=t,o=r,s=t*r+r,a=1,u=[],c=[],h=(12&n)>>2,f=0;i>0&&o>0;)a=(s=i*r+o)%2,f=e[s>>>1],f=a?15&f:f>>>4,0===h?0===(f>>=2)&&(i--,o--,u.push(1),c.push(1)):2===h?(u.push(-1),c.push(1),f&=2,o--):(u.push(1),c.push(-1),f&=1,i--),h=f;return i>0?(u.push(i),c.push(-i)):o>0&&(u.push(-o),c.push(o)),[u,c]}var x,z=function(){function e(e,t,r,n){void 0===n&&(n=0),this.nbSeq=0,this.length=0,this.weight=0,this.length=e,this.m_bAllGaps=new i(e),this.m_uSortOrder=new Uint8Array(t*e),this.m_fcCounts=new Uint8Array(t*e),this.m_wCounts=new Float32Array(t*e),this.m_AAScores=new Float32Array(t*e),this.m_uResidueGroup=new Uint8Array(e),this.m_fOcc=new Float32Array(e),this.m_fcStartOcc=new Float32Array(e),this.m_fcEndOcc=new Float32Array(e),this.m_ScoreGapOpen=new Float32Array(e),this.m_ScoreGapClose=new Float32Array(e),this.nbSeq=r,this.alphaSize=t,this.weight=n}return e.prototype.getProxy=function(e){return void 0===e&&(e=0),new R(this,e)},e}(),R=function(){function e(e,t){void 0===t&&(t=0),this.idx=0,this.alphaSize=4,this.offset=0,this.prof=e,this.idx=t,this.alphaSize=e.alphaSize,this.offset=this.idx*this.alphaSize}return e.prototype.setProxy=function(e){this.idx=e,this.offset=e*this.alphaSize},Object.defineProperty(e.prototype,"m_bAllGaps",{get:function(){return this.prof.m_bAllGaps.get(this.idx)},set:function(e){e?this.prof.m_bAllGaps.set(this.idx):this.prof.m_bAllGaps.clear(this.idx)},enumerable:!1,configurable:!0}),Object.defineProperty(e.prototype,"m_uSortOrder",{get:function(){return this.prof.m_uSortOrder.subarray(this.offset,this.offset+this.alphaSize)},set:function(e){this.prof.m_uSortOrder.set(e,this.offset)},enumerable:!1,configurable:!0}),Object.defineProperty(e.prototype,"m_uResidueGroup",{get:function(){return this.prof.m_uResidueGroup[this.idx]},set:function(e){this.prof.m_uResidueGroup[this.idx]=e},enumerable:!1,configurable:!0}),Object.defineProperty(e.prototype,"m_fcCounts",{get:function(){return this.prof.m_fcCounts.subarray(this.offset,this.offset+this.alphaSize)},set:function(e){this.prof.m_fcCounts.set(e,this.offset)},enumerable:!1,configurable:!0}),Object.defineProperty(e.prototype,"m_wCounts",{get:function(){return this.prof.m_wCounts.subarray(this.offset,this.offset+this.alphaSize)},set:function(e){this.prof.m_wCounts.set(e,this.offset)},enumerable:!1,configurable:!0}),Object.defineProperty(e.prototype,"m_AAScores",{get:function(){return this.prof.m_AAScores.subarray(this.offset,this.offset+this.alphaSize)},set:function(e){this.prof.m_AAScores.set(e,this.offset)},enumerable:!1,configurable:!0}),Object.defineProperty(e.prototype,"m_fOcc",{get:function(){return this.prof.m_fOcc[this.idx]},set:function(e){this.prof.m_fOcc[this.idx]=e},enumerable:!1,configurable:!0}),Object.defineProperty(e.prototype,"m_fcStartOcc",{get:function(){return this.prof.m_fcStartOcc[this.idx]},set:function(e){this.prof.m_fcStartOcc[this.idx]=e},enumerable:!1,configurable:!0}),Object.defineProperty(e.prototype,"m_fcEndOcc",{get:function(){return this.prof.m_fcEndOcc[this.idx]},set:function(e){this.prof.m_fcEndOcc[this.idx]=e},enumerable:!1,configurable:!0}),Object.defineProperty(e.prototype,"m_ScoreGapOpen",{get:function(){return this.prof.m_ScoreGapOpen[this.idx]},set:function(e){this.prof.m_ScoreGapOpen[this.idx]=e},enumerable:!1,configurable:!0}),Object.defineProperty(e.prototype,"m_ScoreGapClose",{get:function(){return this.prof.m_ScoreGapClose[this.idx]},set:function(e){this.prof.m_ScoreGapClose[this.idx]=e},enumerable:!1,configurable:!0}),e}();function T(e,t,r,n){void 0===n&&(n=0);for(var i=e.encodedSeq,o=i.length,s=new z(o,r.abSize,1,t),a=t,u=r.gapOP,c=s.getProxy(),h=0,f=0;f<o;f++){c.setProxy(f),h=i[f],c.m_fcCounts[h]=1,c.m_wCounts[h]=a,c.m_uSortOrder[c.m_uResidueGroup]=h,c.m_uResidueGroup=1,c.m_fOcc=a,c.m_ScoreGapOpen=u/2*a,c.m_ScoreGapClose=u/2*a;for(var p=0;p<r.abSize;p++)c.m_AAScores[p]=a*r.scoringMatrix[p][h]}return 1&n||(s.m_ScoreGapOpen[0]/=2,s.m_ScoreGapOpen[o-1]/=2),2&n||(s.m_ScoreGapClose[0]/=2,s.m_ScoreGapClose[o-1]/=2),s}function N(e,t,r,n,i,o){void 0===o&&(o=0);for(var s=P(r),a=e.nbSeq+t.nbSeq,u=e.weight+t.weight,c=new z(s,i.abSize,a,u),h=i.gapOP,f=c.getProxy(),p=e.getProxy(),d=t.getProxy(),l=M(r),g=M(n),m=0,S=0,v=!1,_=!1,y=!1,w=!1,b=0,O=0;O<s;O++){if(f.setProxy(O),m=l[O],p.setProxy(m),-1===m?(v=!0,_=O===s-1||l[O+1]>0):(v=!1,_=!1),S=g[O],d.setProxy(S),-1===S?(y=!0,w=O==s-1||g[O+1]>0):(y=!1,w=!1),b=0,-1!==m&&-1!==S){f.m_fcStartOcc=p.m_fcStartOcc+d.m_fcStartOcc,f.m_fcEndOcc=p.m_fcEndOcc+d.m_fcEndOcc,f.m_fOcc=p.m_fOcc+d.m_fOcc,f.m_fcCounts.set(p.m_fcCounts),f.m_wCounts.set(p.m_wCounts);for(var A=0;A<t.alphaSize;A++){var q=t.m_fcCounts[A];q&&(f.m_fcCounts[A]+=q,f.m_wCounts[A]+=t.m_wCounts[A]),f.m_fcCounts[A]&&(f.m_uSortOrder[b]=A,b++),f.m_AAScores[A]=p.m_AAScores[A]+d.m_AAScores[A]}f.m_uResidueGroup=b}-1==m&&(f.m_fcStartOcc=d.m_fcStartOcc+(v?e.weight:0),f.m_fcEndOcc=d.m_fcEndOcc+(_?e.weight:0),f.m_fOcc=d.m_fOcc,f.m_fcCounts.set(d.m_fcCounts),f.m_wCounts.set(d.m_wCounts),f.m_uSortOrder.set(d.m_uSortOrder),f.m_AAScores.set(d.m_AAScores),f.m_uResidueGroup=d.m_uResidueGroup),-1==S&&(f.m_fcStartOcc=p.m_fcStartOcc+(y?t.weight:0),f.m_fcEndOcc=p.m_fcEndOcc+(w?t.weight:0),f.m_fOcc=p.m_fOcc,f.m_fcCounts.set(p.m_fcCounts),f.m_wCounts.set(p.m_wCounts),f.m_uSortOrder.set(p.m_uSortOrder),f.m_AAScores.set(p.m_AAScores),f.m_uResidueGroup=p.m_uResidueGroup),f.m_ScoreGapOpen=h/2*(1-f.m_fcStartOcc)*u,f.m_ScoreGapClose=h/2*(1-f.m_fcEndOcc)*u}return 1&o||(c.m_ScoreGapOpen[0]/=2,c.m_ScoreGapOpen[s-1]/=2),2&o||(c.m_ScoreGapClose[0]/=2,c.m_ScoreGapClose[s-1]/=2),c}function U(e){return"seq"in e}function D(e,t){for(var r=1/0,n=0,i=0,o=e.length,s=0;s<o;s++)for(var a=s+1;a<o;a++)e[s][a]<r&&(r=e[s][a],n=s,i=a);return[{profile:null,childA:t[n],childB:t[i],distance:r,estring:[],msa:[],numSeq:[],type:x.NODE,depth:0,id:"",tabWeight:[],parent:-1,weight:0},n,i]}function j(e,t,r,n){for(var i=t>r?[t,r]:[r,t],o=i[0],s=i[1],a=[],u=e.length,c=0;c<u;c++)if(c!=t&&c!=r){var h=.1*(e[c][t]+e[c][r])/2+.9*Math.min(e[c][t],e[c][r]);e[c].push(h),a.push(h),e[c].splice(o,1),e[c].splice(s,1)}return a.push(0),e.push(a),e.splice(o,1),e.splice(s,1),n.splice(o,1),n.splice(s,1),e}function F(e,r){var o=function(e,n){var i=a[e.childA],s=a[e.childB],u={};if(U(i)||0!==i.msa.length||o(i,n),U(s)||0!==s.msa.length||o(s,n),U(i))if(U(s)){var h=I(i.seq,s.seq,r);u.score=h.score,e.numSeq=[i.numSeq[0],s.numSeq[0]],e.estring=h.estrings,e.profile=N(i.profile,s.profile,h.estrings[0],h.estrings[1],r),c.add("seq "+i.numSeq+" - seq "+s.numSeq)}else{s.tabWeight=n.filter((function(e,t){return s.numSeq.includes(t)}));var f=function(e,t,r,n){void 0===n&&(n=0);var i,o,s=t.seq,a=s.rawSeq.length,u=e.profile.length,c=0,h=[],f=[],p=new Uint8Array(Math.ceil((u+1)*(a+1)/2)),d=0,l=0,g=0,m=0,S=0,v=0,_=0,y=0,w=0,b=0,O=new Float32Array(r.abSize),A=r.gapOP,q=A/2,E=A/2,P=1&n?0:-A/4,M=2&n?0:-A/4,I=s.encodedSeq,x=e.profile;h[0]=0,f[0]=-1/0;for(var z=1;z<=a;z++)h[z]=q+E+P,f[z]=-1/0;for(var R=1;R<=u;R++){for(w=x.m_ScoreGapOpen[R-1],b=x.m_ScoreGapClose[R-1],O=x.m_AAScores.subarray((R-1)*r.abSize,R*r.abSize),v=P,m=-1/0,f[0]=x.m_ScoreGapOpen[0],z=1;z<=a;z++)g=0,d=h[z]+w,z===a&&(d-=w/2),d>=f[z]?f[z]=d:g+=1,l=v+q,R===u&&(l+=M),l>=(S=m)?m=l:g+=2,_=f[z]+b,y=S+E,c=h[z-1]+O[I[z-1]],z===a&&(y+=M),h[z-1]=v,c>=y?c>=_?v=c:(v=_,g+=4):y>=_?(v=y,g+=8):(v=_,g+=4),o=(i=R*a+z)%2,p[i>>>=1]+=o?g:g<<4;h[a]=v}var T=Math.max(c,m,f[f.length-2]),N=G(p,u,a,g),U=N[0],D=N[1];return{estrings:[C(U),C(D)],score:T}}(s,i,r);u.score=f.score,e.numSeq=t([i.numSeq[0]],s.numSeq),e.estring=t([f.estrings[1]],s.estring.map((function(e){return O(f.estrings[0],e)}))),e.profile=N(i.profile,s.profile,f.estrings[1],f.estrings[0],r),c.add("seq "+i.numSeq+" - MSA "+s.numSeq)}else if(!U(s)){s.tabWeight=n.filter((function(e,t){return s.numSeq.includes(t)})),i.tabWeight=n.filter((function(e,t){return i.numSeq.includes(t)}));var p=function(e,t,r,n){void 0===n&&(n=0);var i,o,s=e.profile.length,a=t.profile.length,u=new Uint8Array(Math.ceil((s+1)*(a+1)/2)),c=[],h=[],f=0,p=0,d=0,l=0,g=0,m=0,S=0,v=0,_=0,y=0,w=0,b=0,O=new Float32Array(r.abSize),A=new Uint8Array(r.abSize),q=0,E=0,P=0,M=0,I=0,x=0,z=1&n?1:2,R=t.profile,T=e.profile;for(c[0]=0,h[0]=-1/0,I=1;I<=a;I++)c[I]=R.m_ScoreGapOpen[0]+R.m_ScoreGapClose[I-1]/z,h[I]=-1/0;var N=T.m_ScoreGapOpen[0];for(M=1;M<=s;M++){for(O=T.m_AAScores.subarray(r.abSize*(M-1),r.abSize*M),w=T.m_ScoreGapOpen[M-1],v=N+(b=T.m_ScoreGapClose[M-1])/z,m=-1/0,h[0]=T.m_ScoreGapOpen[0],I=1;I<=a;I++){for(g=0,p=c[I]+w,I!==a||2&n||(p-=w/2),p>=h[I]?h[I]=p:g+=1,d=v+R.m_ScoreGapOpen[I-1],l=S=m,M!==s||2&n||(d-=R.m_ScoreGapOpen[I-1]/2),d>=l?m=d:g+=2,f=c[I-1],P=R.m_uResidueGroup[I-1],x=0,E=(I-1)*r.abSize,A=R.m_uSortOrder.subarray(E,E+P);x<P;)q=A[x],f+=R.m_wCounts[E+q]*O[q],x++;_=h[I]+b,y=S+R.m_ScoreGapClose[I-1],1!==I||2&n||(_-=b/2),I!==a||2&n||(_-=b/2,y-=R.m_ScoreGapClose[a-1]/2),c[I-1]=v,f>=y?f>=_?v=f:(v=_,g+=4):y>=_?(v=y,g+=8):(v=_,g+=4),o=(i=M*a+I)%2,u[i>>>=1]+=o?g:g<<4}c[a]=v}var U=Math.max(f,m,h[I-1]),D=G(u,s,a,g),j=D[0],F=D[1];return{estrings:[C(j),C(F)],score:U}}(i,s,r);u.score=p.score,e.estring=t(t([],i.estring.map((function(e){return O(p.estrings[0],e)}))),s.estring.map((function(e){return O(p.estrings[1],e)}))),e.profile=N(i.profile,s.profile,p.estrings[0],p.estrings[1],r),c.add("MSA "+i.numSeq+" - MSA "+s.numSeq)}return e.numSeq=t(t([],i.numSeq),s.numSeq),u.score};c.add("Start Progressive Alignment");var s=function(e){var t=e[0].type===n.PROTEIN,r=t?6:4,o=Math.pow(r,6),s=e.map((function(e){for(var n=new i(o),s=t?e.compressedSeq:e.encodedSeq,a=0,u=0;u<=5;u++)a+=s[u]*Math.pow(r,u);n.set(a);for(var c=Math.pow(r,5),h=(u=6,s.length);u<h;u++)a-=s[u-6],a/=6,a+=s[u]*c,n.set(a);return n}));c.add("K-mer bitset computation");for(var a,u,h,f=e.length,p=e.map((function(){return[]})),d=0;d<f;d++){p[d][d]=0,a=s[d],u=e[d].compressedSeq.length;for(var l=d+1;l<f;l++)h=1-a.getIntersectionSize(s[l])/u,p[l][d]=p[d][l]=h}return c.add("K-mer distance computation"),p}(e);c.add("K-mer distance matrix");var a=function(e,t){for(var r,n,i,o,s=[],a=e.length,u=[],c=0;c<a;c++)s[c]={type:x.LEAF,seq:t[c],profile:null,childA:c,childB:c,distance:0,numSeq:[c],msa:[],id:c.toString(),weight:0,parent:-1,depth:0},u[c]=c;c=0;for(var h=a-1;c<h;c++){n=(r=D(e,u))[0],i=r[1],o=r[2];var f=s[n.childA],p=s[n.childB];n.id=f.id<p.id?"|"+f.id+","+p.id+"|":"|"+p.id+","+f.id+"|",n.depth=Math.max(f.depth,p.depth)+1,p.parent=f.parent=s.length,u.push(a+c),s.push(n),e=j(e,i,o,u)}return s[s.length-1].type=x.ROOT,s}(s,e),u=a[a.length-1];c.add("Build Tree");var h=function(e){var t=[],r=0;!function t(n){var i=0,o=0;switch(n.type){case x.ROOT:n.weight=0,t(e[n.childA]),t(e[n.childB]);break;case x.NODE:i=e[n.parent].distance-n.distance,o=n.id.split(",").length,n.weight=e[n.parent].weight+i/o,t(e[n.childA]),t(e[n.childB]);break;case x.LEAF:i=e[n.parent].distance,n.weight=e[n.parent].weight+i,r+=n.weight}}(e[e.length-1]);for(var n=0;e[n].type==x.LEAF;)t[n]=e[n].weight/=r,n++;return t}(a);c.add("Compute Weights - Start MSA"),function(e,t,r){for(var n=0,i=e[n];U(i)&&n<e.length;)i.profile=T(i.seq,i.weight,t,r),i=e[++n]}(a,r),o(u,h);var f,p,d,l=(f=u.estring.slice(),p=u.numSeq,d=p.slice(),p.forEach((function(e,t){return d[e]=t})),d.map((function(e){return f[e]})));return c.add("End MSA computation"),l}!function(e){e[e.LEAF=0]="LEAF",e[e.NODE=1]="NODE",e[e.ROOT=2]="ROOT"}(x||(x={}));var L=function(){function e(e){this.storeSize=e,this.store=new Array(e),this.head=0,this.tail=0}return Object.defineProperty(e.prototype,"size",{get:function(){return(this.tail-this.head+this.storeSize)%this.storeSize},enumerable:!1,configurable:!0}),Object.defineProperty(e.prototype,"isEmpty",{get:function(){return this.head===this.tail},enumerable:!1,configurable:!0}),e.prototype.popHead=function(){if(this.isEmpty)return null;var e=this.store[this.head];return this.head++,this.head>=this.storeSize&&(this.head=0),e},e.prototype.popTail=function(){return this.isEmpty?null:(this.tail--,this.tail<0&&(this.tail=this.storeSize-1),this.store[this.tail])},e.prototype.pushHead=function(e){this.head--,this.head<0&&(this.head=this.storeSize-1),this.store[this.head]=e,this.head===this.tail&&(this.tail=(this.tail-1+this.storeSize)%this.storeSize)},e.prototype.pushTail=function(e){this.store[this.tail]=e,this.tail++,this.tail>=this.storeSize&&(this.tail=0),this.head===this.tail&&(this.head=(this.head-1+this.storeSize)%this.storeSize)},e.prototype.getHead=function(e){void 0===e&&(e=0);var t=0===e?this.head:(this.head+e+this.storeSize)%this.storeSize;return this.store[t]},e.prototype.getTail=function(e){void 0===e&&(e=0);var t=(this.tail-1-e+this.storeSize)%this.storeSize;return this.store[t]},e}();function k(e,t,r){for(var n=new Map,i=new L(r),o=0|t,s=new Uint16Array(e.encodedSeq.length-o+1),a=0,u=0,c=0;c<o;c++)a|=e.encodedSeq[o-c-1]<<2*c;s[u++]=a;c=o;for(var h=e.encodedSeq.length;c<h;c++)a=(a<<2)+e.encodedSeq[c],s[u++]=a;var f=r-o,p=e.rawSeq.length*(2/(r+1))*2|0,d={kmer:new Uint16Array(p),kmerPos:new Uint16Array(p),winPos:new Uint16Array(p),winPosEnd:new Uint16Array(p),count:0},l=NaN;for(c=0;c<s.length;c++){for(var g=s[c];!i.isEmpty&&s[i.getTail()]>=g;)i.popTail();if(i.pushTail(c),!(c<f)){for(;i.getHead()<=c-1-f;)i.popHead();var m=i.getHead();if(m===l)d.winPosEnd[d.count-1]=c+o;else{var S=s[m],v=d.count;if(d.count++,d.kmer[v]=S,d.kmerPos[v]=m,d.winPos[v]=c-f,d.winPosEnd[v]=c+o,l=m,n.has(S))n.get(S).push(v);else n.set(S,[v])}}}return[n,d,s]}function B(e,t,r,n,i){var o,s,a={},u=null!=n?n:k(e,8,16),h=u[0],f=u[1],p=u[2],d=null!=i?i:k(t,8,16),l=d[0],g=d[1],m=d[2];c.add("Extract Minimizers"),a["Nb Minimizers"]={a:f.count,b:g.count},a["Dupl. Minimizers"]={a:f.count-h.size,b:g.count-l.size};for(var S=0,v=Math.max(e.rawSeq.length,t.rawSeq.length),_=new Map,y=0;y<f.count;y++){var w=f.kmer[y];if(l.has(w)){var b=l.get(w);if(!(b.length>4))for(var O=e.rawSeq.substring(f.winPos[y],f.winPosEnd[y]),A=0;A<b.length;A++){var q=b[A],M=t.rawSeq.substring(g.winPos[q],g.winPosEnd[q]),G=Math.min(O.length,M.length),x=0,z=void 0;if((G==O.length?0===M.indexOf(O):0===O.indexOf(M))?(z={diagId:x=f.winPos[y]-g.winPos[q],begin:f.winPos[y],end:f.winPos[y]+G},S++):z={diagId:x=f.kmerPos[y]-g.kmerPos[q],begin:f.kmerPos[y],end:f.kmerPos[y]+8-1},_.has(x))_.get(x).push(z);else _.set(x,[z])}}}c.add("Filter Common Minimizers"),a["Nb common windows"]={all:S};var R=[];if(_.forEach((function(e){if(e.length){var t=e[0];R.push(t),e.forEach((function(e){e.begin<=t.end?t.end=e.end:(t=e,R.push(t))}))}})),c.add("Merge Minimizers"),a["Nb diagonals"]={all:R.length},0===R.length)return console.table(a),I(e,t,r).estrings;R.sort((function(e,t){var r=e.begin+e.diagId-t.begin-t.diagId;return 0===r?e.diagId-t.diagId:r})),c.add("Sort Minimizers");var T=[[0]],N=R[0],U=N.end-N.begin,D=N.end-N.diagId,j=N.begin-N.diagId,F=[U],L=[D];for(y=1;y<R.length;y++){U=(N=R[y]).end-N.begin,D=N.end-N.diagId,j=N.begin-N.diagId;var B=R[T[0][0]];if(N.diagId>=B.diagId&&D<L[0]){if(U<F[0]){T.unshift([y]),F.unshift(U),L.unshift(D);continue}T[0]=[y],F[0]=U,L[0]=D}else for(A=T.length-1;A>=0;A--){var W=T[A],K=R[W[W.length-1]];if(N.begin>=K.end&&L[A]<j){for(var Y=Math.abs(N.diagId-K.diagId)/v,Q=F[A]+U-Y,$=A;$<T.length&&F[$]<Q;)$++;F.splice($,0,Q),L.splice($,0,D);var J=W.slice();for(J.push(y),T.splice($,0,J);$--;)L[$]>D&&F[$]<Q&&(T.splice($,1),F.splice($,1),L.splice($,1));break}}}R=null!==(s=null===(o=T.pop())||void 0===o?void 0:o.map((function(e){return R[e]})))&&void 0!==s?s:[],c.add("Filter Minimizers - Optimal list");var X=0;R.forEach((function(e){X+=e.end-e.begin})),X/=v,a["Strict Coverage"]={all:X},a["Nb filtered diagonals"]={all:R.length};var Z=R[0],ee=Z?[Z]:[],te=[];for(y=1;y<R.length;y++){var re=R[y];if(re.diagId!==Z.diagId){for(ie=(ne=Z.end)-Z.diagId;e.encodedSeq[ne]===t.encodedSeq[ie]&&ne<re.begin&&ie<re.begin-re.diagId;)ne++,ie++;for(Z.end=ne,ie=(ne=re.begin)-re.diagId;e.encodedSeq[ne]===t.encodedSeq[ie]&&ne>Z.end&&ie>Z.end-Z.diagId;)ne--,ie--;re.begin=ne+1,te.push({endDiagId:re.diagId,beginDiagId:Z.diagId,begin:Z.end,end:re.begin}),ee.push(re),Z=re}else{for(var ne,ie=(ne=Z.end)-re.diagId,oe=0;ne<re.begin&&ie<re.begin-re.diagId&&oe<2;)ie+=8,oe=V(p[ne+=8],m[ie])<=2?0:oe+1;for(ne>Z.end&&(Z.end=Math.max(ne-2-8,Z.end)),ie=(ne=re.begin-8)-re.diagId,oe=0;ne>Z.end&&ie>Z.end-Z.diagId&&oe<2;)ie-=8,oe=V(p[ne-=8],m[ie])<=2?0:oe+1;if(ne<re.begin-16&&(re.begin=ne+24+2),re.begin-Z.end<=24){Z.end=re.end;continue}te.push({beginDiagId:Z.diagId,endDiagId:re.diagId,begin:Z.end,end:re.begin}),ee.push(re),Z=re}}c.add("Diagonals extension");var se=0;ee.forEach((function(e){se+=e.end-e.begin})),se/=v,a.Coverage={all:se},a["Extended Diagonals"]={all:ee.length},console.table(a),N=ee[0];var ae=[],ue=[];for(y=0;y<ee.length;y++){N=ee[y],ae.push(N.end-N.begin),ue.push(N.end-N.begin);var ce=te[y];if(ce){if(ce.begin===ce.end){var he=Math.abs(ce.endDiagId-ce.beginDiagId);ae.push(-he),ue.push(he);continue}if(ce.begin-ce.beginDiagId==ce.end-ce.endDiagId){he=ce.end-ce.begin;ue.push(-he),ae.push(he);continue}var fe=I({rawSeq:e.rawSeq.substring(ce.begin,ce.end),type:e.type,compressedSeq:new Uint8Array(0),encodedSeq:e.encodedSeq.subarray(ce.begin,ce.end)},{rawSeq:t.rawSeq.substring(ce.begin-ce.beginDiagId,ce.end-ce.endDiagId),type:t.type,compressedSeq:new Uint8Array(0),encodedSeq:t.encodedSeq.subarray(ce.begin-ce.beginDiagId,ce.end-ce.endDiagId)},r,3);ae.push.apply(ae,fe.estrings[0]),ue.push.apply(ue,fe.estrings[1])}}H(ae,e.encodedSeq.length),H(ue,t.encodedSeq.length);var pe=C(ae,1),de=C(ue,1),le=P(pe),ge=P(de);return le>ge?de=E(de,[ge-le]):le<ge&&(pe=E(pe,[le-ge])),c.add("Fill between diagonals"),[pe,de]}function H(e,t){var r=t-function(e){for(var t=0,r=0,n=0;n<e.length;n++)(r=e[n])<0||(t+=r);return t}(e);r>0?e.push(r):r<0&&(e[e.length-1]+=r)}function V(e,t){var n=e^t;return r(n=1431655765&n|n>>1&1431655765)}var W={sequenceType:"auto",gapchar:"-",alignmentMethod:"auto",debug:!1},K=new(function(){function t(){this.sequences=[],this.typeSeq=n.UNSET,this.config=e({},W)}return t.prototype.addSequence=function(e){var t=this;if(Array.isArray(e))e.forEach((function(e){t.addSequence(e)}));else{if("string"!=typeof e)throw new TypeError("String type expected for sequences to add.");if(e=e.toUpperCase(),"auto"===this.config.sequenceType){var r=g(e);if(this.typeSeq===n.UNSET)this.typeSeq=r;else if(this.typeSeq!==r)throw new Error("All sequences must be of same type.")}this.sequences.push(w(e))}},t.prototype.reset=function(){this.sequences=[],this.typeSeq=n.UNSET,this.setDefaultConfiguration()},t.prototype.setDefaultConfiguration=function(){this.config=e({},W)},t.prototype.setUserConfiguration=function(e){var t;function r(e,t){return t.some((function(t){return t==e}))}e.method&&r(e.method,["complete","diag"])&&(this.config.alignmentMethod=e.method),e.type&&r(e.type,["amino","nucleic"])&&(this.config.sequenceType=e.type,this.typeSeq="amino"===e.type?n.PROTEIN:n.NUCLEIC),1===(null===(t=e.gapchar)||void 0===t?void 0:t.length)&&(this.config.gapchar=e.gapchar),void 0!==e.debug&&(this.config.debug=!!e.debug)},t.prototype.align=function(e,t){var r=this;return c.start(),new Promise((function(i,o){if(!Array.isArray(e)||e.some((function(e){return"string"!=typeof e})))return o("Array of sequences expected");if(e.length<2)return o("At least 2 sequences are required.");r.reset(),t&&r.setUserConfiguration(t),r.addSequence(e),c.add("Prepared sequences");var u=function(e,t){var r,i,o,u=e===n.PROTEIN?s:a,c=null!==(r=null==t?void 0:t.matrix)&&void 0!==r?r:u.matrix,h=(null==t?void 0:t.gapextend)?2*-t.gapextend:u.center;return{type:e,scoringMatrix:(o=h,c.map((function(e){return e.map((function(e){return e+o}))}))),gapOP:null!==(i=null==t?void 0:t.gapopen)&&void 0!==i?i:u.gapOP,abSize:e===n.PROTEIN?20:4}}(r.typeSeq,t);c.add("Get sequences type");var h,f=!1;(f="auto"===r.config.alignmentMethod?r.sequences.some((function(e){return e.rawSeq.length>1600})):"diag"===r.config.alignmentMethod,2==r.sequences.length)?h=f?B(r.sequences[0],r.sequences[1],u):I(r.sequences[0],r.sequences[1],u).estrings:h=f?function(e,t){for(var r=e.map((function(e){return k(e,8,16)})),n=[],i=[],o=[],s=0;s<r.length;s++)if(0!==s){var a=B(e[0],e[s],t,r[0],r[s]);i=0===i.length?a[0]:A(i,a[0]),o[s]=a[0],n[s]=a[1]}for(s=0;s<r.length;s++)if(0!==s){var u=q(i,o[s]);if(!u)throw console.error("An error occured while computing the center diff. for sequence #"+s,i,o[s]),console.dir(n),console.dir(o),new RangeError;n[s]=O(u,n[s])}else n[s]=i;return n}(r.sequences,u):F(r.sequences,u);return c.summary(),i(function(e,t,r){for(var n,i=null!==(n=null==r?void 0:r.gapchar)&&void 0!==n?n:"-",o=[],s=0;s<e.length;s++)o.push(b(e[s].rawSeq,t[s],{gapchar:i}));return o}(r.sequences,h,{gapchar:r.config.gapchar}))})).catch((function(e){return c.summary(),Promise.reject(e)}))},t}());export{K as default};
//# sourceMappingURL=biomsa.esm.js.map
