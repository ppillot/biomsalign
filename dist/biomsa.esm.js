function t(t){return 16843009*((t=(858993459&(t-=t>>>1&1431655765))+(t>>>2&858993459))+(t>>>4)&252645135)>>>24}class e{constructor(t,e){this.length=t,this.t=new Uint32Array(t+32>>>5),!0===e&&this.i()}get(t){return 0!=(this.t[t>>>5]&1<<t)}set(t){this.t[t>>>5]|=1<<t}clear(t){this.t[t>>>5]&=~(1<<t)}o(t){this.t[t>>>5]^=1<<t}h(t,e,r){if(e<t)return;const n=this.t,i=!0===r?4294967295:0,s=t>>>5,o=e>>>5;for(let t=s;t<o;++t)n[t]=i;const h=s<<5,c=o<<5;if(!0===r)if(e-t<32)for(let r=t,i=e+1;r<i;++r)n[r>>>5]|=1<<r;else{for(let e=t,r=h;e<r;++e)n[e>>>5]|=1<<e;for(let t=c,r=e+1;t<r;++t)n[t>>>5]|=1<<t}else if(e-t<32)for(let r=t,i=e+1;r<i;++r)n[r>>>5]&=~(1<<r);else{for(let e=t,r=h;e<r;++e)n[e>>>5]&=~(1<<e);for(let t=c,r=e+1;t<r;++t)n[t>>>5]&=~(1<<t)}return this}l(t,e){return this.h(t,e,!0)}u(t,e){return this.h(t,e,!1)}A(...t){const e=this.t,r=t.length;for(let n=0;n<r;++n){const r=t[n];e[r>>>5]|=1<<r}return this}p(...t){const e=this.t,r=t.length;for(let n=0;n<r;++n){const r=t[n];e[r>>>5]&=~(1<<r)}return this}i(){return this.h(0,this.length-1,!0)}M(){return this.h(0,this.length-1,!1)}m(){const t=this.t.length,e=this.t,r=32-this.length%32;for(let r=0;r<t-1;++r)e[r]=~e[r];return e[t-1]=~(e[t-1]<<r)>>>r,this}S(t,e,r){if(e<t)return;const n=this.t,i=!0===r?4294967295:0,s=t>>>5,o=e>>>5;for(let t=s;t<o;++t)if(n[t]!==i)return!1;if(e-t<32){for(let i=t,s=e+1;i<s;++i)if(!!(n[i>>>5]&1<<i)!==r)return!1}else{const i=o<<5;for(let e=t,i=s<<5<<5;e<i;++e)if(!!(n[e>>>5]&1<<e)!==r)return!1;for(let t=i,s=e+1;t<s;++t)if(!!(n[t>>>5]&1<<t)!==r)return!1}return!0}U(t,e){return this.S(t,e,!0)}_(t,e){return this.S(t,e,!1)}q(){return this.S(0,this.length-1,!0)}v(){return this.S(0,this.length-1,!1)}O(...t){const e=this.t,r=t.length;for(let n=0;n<r;++n){const r=t[n];if(0==(e[r>>>5]&1<<r))return!1}return!0}C(...t){const e=this.t,r=t.length;for(let n=0;n<r;++n){const r=t[n];if(0!=(e[r>>>5]&1<<r))return!1}return!0}G(t){const e=this.t,r=t.t,n=Math.min(e.length,r.length);for(let t=0;t<n;++t)if(e[t]!==r[t])return!1;return!0}P(){const e=this.t.length,r=this.t;let n=0;for(let i=0;i<e;++i)n+=t(r[i]);return n}R(t){const e=this.t,r=t.t,n=Math.min(e.length,r.length);for(let t=0;t<n;++t)e[t]=e[t]&~r[t];for(let t=e.length;t<n;++t)e[t]=0;return this}I(t){const e=this.t,r=t.t,n=Math.min(e.length,r.length);for(let t=0;t<n;++t)e[t]|=r[t];for(let t=e.length;t<n;++t)e[t]=0;return this}F(t){const e=this.t,r=t.t,n=Math.min(e.length,r.length);for(let t=0;t<n;++t)e[t]&=r[t];for(let t=e.length;t<n;++t)e[t]=0;return this}T(t){const e=this.t,r=t.t,n=Math.min(e.length,r.length);for(let t=0;t<n;++t)if(0!=(e[t]&r[t]))return!0;return!1}D(e){const r=this.t,n=e.t,i=Math.min(r.length,n.length);let s=0;for(let e=0;e<i;++e)s+=t(r[e]&n[e]);return s}B(t){const r=this.t,n=t.t,i=Math.min(r.length,n.length),s=new Uint32Array(i),o=Object.create(e.prototype);o.t=s,o.length=Math.min(this.length,t.length);for(let t=0;t<i;++t)s[t]=r[t]&n[t];return o}forEach(e){const r=this.t.length,n=this.t;let i=0;for(let s=0;s<r;++s){let r=n[s];for(;0!==r;){const n=r&-r;e((s<<5)+t(n-1),i),r^=n,++i}}}toArray(){const e=this.t,r=Array(this.P()),n=this.t.length;let i=0;for(let s=0;s<n;++s){let n=e[s];for(;0!==n;){const e=n&-n;r[i++]=(s<<5)+t(e-1),n^=e}}return r}toString(){return"{"+this.toArray().join(",")+"}"}N(){const t=this.toArray().join(",");return t?"@"+t:"NONE"}clone(){const t=Object.create(e.prototype);return t.length=this.length,t.t=new Uint32Array(this.t),t}}const r={matrix:[[58,23,-12,-7,-44,10,-23,-14,-14,-27,-17,-8,1,-9,-22,23,15,5,-74,-45,0],[23,224,-67,-63,-50,-30,-29,1,-56,-41,-6,-33,-44,-53,-43,15,2,18,-93,-6,0],[-12,-67,111,59,-104,-4,4,-84,6,-88,-65,48,-13,18,-29,5,-7,-63,-105,-73,0],[-7,-63,59,85,-83,-17,-1,-63,25,-60,-47,15,-12,40,-8,1,-7,-47,-108,-51,0],[-44,-50,-104,-83,144,-93,4,12,-74,36,30,-64,-67,-56,-65,-43,-41,-3,63,104,0],[10,-30,-4,-17,-93,140,-32,-95,-27,-91,-75,4,-36,-29,-32,5,-26,-68,-80,-79,0],[-23,-29,4,-1,4,-32,137,-50,6,-37,-42,21,-23,27,19,-4,-12,-44,-13,48,0],[-14,1,-84,-63,12,-95,-50,86,-53,53,47,-62,-60,-47,-55,-43,-8,69,-27,-24,0],[-14,-56,6,25,-74,-27,6,-53,75,-48,-30,13,-12,34,68,-3,-4,-44,-71,-49,0],[-27,-41,-88,-60,36,-91,-37,53,-48,88,62,-63,-48,-36,-48,-47,-25,36,-11,-4,0],[-17,-6,-65,-47,30,-75,-42,47,-30,62,103,-45,-54,-21,-31,-35,-9,31,-46,-20,0],[-8,-33,48,15,-64,4,21,-62,13,-63,-45,89,-25,12,2,22,10,-51,-79,-29,0],[1,-44,-13,-12,-67,-36,-23,-60,-12,-48,-54,-25,160,-6,-20,5,-12,-42,-76,-83,0],[-9,-53,18,40,-56,-29,27,-47,34,-36,-21,12,-6,75,34,1,-4,-37,-92,-48,0],[-22,-43,-29,-8,-65,-32,19,-55,68,-48,-31,2,-20,34,113,-10,-14,-49,-58,-39,0],[23,15,5,1,-43,5,-4,-43,-3,-47,-35,22,5,1,-10,53,32,-28,-62,-31,0],[15,2,-7,-7,-41,-26,-12,-8,-4,-25,-9,10,-12,-4,-14,32,68,0,-87,-40,0],[5,18,-63,-47,-3,-68,-44,69,-44,36,31,-51,-42,-37,-49,-28,0,74,-61,-32,0],[-74,-93,-105,-108,63,-80,-13,-27,-71,-11,-46,-79,-76,-92,-58,-62,-87,-61,289,81,0],[-45,-6,-73,-51,104,-79,48,-24,-49,-4,-20,-29,-83,-48,-39,-31,-40,-32,81,162,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]],H:22,$:-300},n={matrix:[[91,-114,-31,-123],[-114,100,-125,-31],[-31,-125,100,-114],[-123,-31,-114,91]],$:-400,k:-60,H:120},i=/^[ACGTRNDEQHILKMFPSWYV]+$/,s=/^[ABCDGHKMNRSTUVWY]+$/i,o=/^[ABCGTRNDEQHIJLKMFPSWYUVXZ]+$/i,h=/[^ACGU]/i,c=/[^ACGT]/i;function l(t){if(!h.test(t)||!c.test(t))return 1;const e=s.test(t),r=i.test(t);if(r&&!e)return 0;if(r&&e)return u(t);if(o.test(t))return u(t);throw Error("Unrecognized sequence type: "+t)}function u(t){let e=0;const r=Math.min(t.length,100);for(let n=0;n<r;n++)switch(t[n]){case"A":case"T":case"U":case"G":case"C":case"N":e++}return e/r>Math.SQRT1_2?1:0}const f=new Uint8Array(96);f.fill(255),f[65]=0,f[67]=1,f[68]=2,f[69]=3,f[70]=4,f[71]=5,f[72]=6,f[73]=7,f[75]=8,f[76]=9,f[77]=10,f[78]=11,f[80]=12,f[81]=13,f[82]=14,f[83]=15,f[84]=16,f[86]=17,f[87]=18,f[89]=19,f[95]=20;const a=new Uint8Array(96);function d(t,e){const r=Math.floor(100*Math.random());if(1===e)switch(t){case"M":return[0,1][r%2];case"R":return[0,2][r%2];case"W":return[0,3][r%2];case"S":return[1,2][r%2];case"Y":return[1,3][r%2];case"K":return[2,3][r%2];case"V":return[0,1,2][r%3];case"H":return[0,1,3][r%3];case"D":return[0,2,3][r%3];case"B":return[1,2,3][r%3];default:return r%4}switch(t){case"B":return[2,11][r%2];case"Z":return[3,13][r%2];case"J":return[7,9][r%2];default:return r%20}}function g(t){return t.map((t=>w[t]))}a.fill(255),a[65]=0,a[67]=1,a[71]=2,a[84]=3,a[85]=3;const w=[0,1,2,2,3,0,4,5,4,5,5,2,0,2,4,0,0,5,3,3];function A(t,e,r){const n=[],i=r.gapchar;let s=0;for(let r=0;r<e.length;r++){let o=e[r];o<0?n.push(i.repeat(-o)):(n.push(t.substr(s,o)),s+=o)}return n.join("")}function p(t,e){const r=[0];let n=0;const i=e.slice();function s(t){(t^r[r.length-1])>=0?r[r.length-1]+=t:r.push(t)}for(let e=0;e<t.length;e++){let r=t[e];if(r<0)s(r);else for(;r>0&&n<i.length;){let t=i[n];t<0?r<-t?(i[n]+=r,s(-r),r=0):(s(t),r+=t,n++):r<t?(i[n]-=r,s(r),r=0):(s(t),r-=t,n++)}}return r}function y(t,e){const r=[];let n=0,i=0,s=t[0],o=e[0];for(;n<t.length||i<e.length;)s!==o?Math.sign(s)!==Math.sign(o)?s<0||void 0===o?(r.push(s),s=t[++n]):(r.push(o),o=e[++i]):s>0?s>o?(r.push(o),s-=o,o=e[++i]):(r.push(s),o-=s,s=t[++n]):(r.push(Math.min(s,o)),s=t[++n],o=e[++i]):(r.push(o),s=t[++n],o=e[++i]);return r}function M(t,e){const r=[];let n=0,i=0,s=e[0],o=t[0],h=0;for(;n<e.length||i<t.length;)if(h=r[r.length-1],s!==o){if(Math.sign(s)===Math.sign(o)){if(s>0){if(o>s)return;h>0?r[r.length-1]=h+o:r.push(o),s-=o,o=t[++i];continue}if(s>o){h>0?r[r.length-1]=h-s:r.push(Math.abs(s)),o-=s,s=e[++n];continue}return}if(!(o<0))return;h<0?r[r.length-1]=h+o:r.push(o),o=t[++i]}else h>0?r[r.length-1]=h+Math.abs(s):r.push(Math.abs(s)),s=e[++n],o=t[++i];return r}function m(t,e){if(Math.sign(t[t.length-1])===Math.sign(e[0])){let r=[...t];return r[r.length-1]+=e[0],r.push(...e.slice(1)),r}return[...t,...e]}function S(t,e=0){const r=[];let n=0,i=0;for(let e=0;e<t.length;e++)i=t[e],(i^n)>=0?n+=i:(r.push(n),n=i);return r.push(n),0===r[0]&&r.shift(),1&e||r.reverse(),r}function U(t){let e=0,r=0;for(let n=0;n<t.length;n++)r=t[n],e+=r<0?-r:r;return e}function _(t){const e=[];let r=0,n=0;for(let i=0;i<t.length;i++)if(n=t[i],n<0)for(let t=0;t<-n;t++)e.push(-1);else for(let t=0;t<n;t++)e.push(r),r++;return e}function q(t,e,r,n=0){const i=t.V,s=e.V,o=i.length,h=s.length;let c=0,l=[],u=[],f=0,a=0,d=0,g=0;const w=new Uint8Array(Math.ceil((o+1)*h/2));let A=0,p=0,y=0,M=0,m=0,U=r.W[0];const _=r.$,q=1&n?0:-_/2,b=2&n?0:-_/2;l[0]=0,u[0]=-1/0;for(let t=1;t<=h;t++)l[t]=_+q,u[t]=-1/0;for(let t=1;t<=o;t++){f=q,a=-1/0,U=r.W[i[t-1]];for(let e=1;e<=h;e++)m=0,A=l[e]+_,y=u[e],e===h&&(A+=b),A>=y?u[e]=A:m+=1,p=f+_,M=a,t===o&&(p+=b),p>=M?a=p:m+=2,c=l[e-1]+U[s[e-1]],l[e-1]=f,c>=a?c>=u[e]?f=c:(f=u[e],m+=4):a>=u[e]?(f=a,m+=8):(f=u[e],m+=4),d=t*h+e,g=d%2,d>>>=1,w[d]+=g?m:m<<4;l[h]=f}var O=Math.max(c,a,u[u.length-1]);const[C,E]=v(w,o,h,m);return{j:[S(C),S(E)],K:O}}function v(t,e,r,n){let i=e,s=r,o=e*r+r,h=o,c=1;const l=[],u=[];let f=(12&n)>>2,a=0;for(;i>0&&s>0;)o=i*r+s,h=o>>>1,c=o%2,a=t[h],a=c?15&a:a>>>4,0===f?(a>>=2,0===a&&(i--,s--,l.push(1),u.push(1))):2===f?(l.push(-1),u.push(1),a&=2,s--):(l.push(1),u.push(-1),a&=1,i--),f=a;return i>0?(l.push(i),u.push(-i)):s>0&&(l.push(-s),u.push(s)),[l,u]}class b{constructor(t,r,n,i=0){this.Y=0,this.length=0,this.weight=0,this.length=t,this.J=new e(t),this.L=new Uint8Array(r*t),this.Z=new Uint8Array(r*t),this.X=new Float32Array(r*t),this.tt=new Float32Array(r*t),this.et=new Uint8Array(t),this.rt=new Float32Array(t),this.nt=new Float32Array(t),this.it=new Float32Array(t),this.st=new Float32Array(t),this.ot=new Float32Array(t),this.Y=n,this.ht=r,this.weight=i}ct(t=0){return new O(this,t)}}class O{constructor(t,e=0){this.lt=0,this.ht=4,this.offset=0,this.ut=t,this.lt=e,this.ht=t.ht,this.offset=this.lt*this.ht}ft(t){this.lt=t,this.offset=t*this.ht}get J(){return this.ut.J.get(this.lt)}set J(t){t?this.ut.J.set(this.lt):this.ut.J.clear(this.lt)}get L(){return this.ut.L.subarray(this.offset,this.offset+this.ht)}set L(t){this.ut.L.set(t,this.offset)}get et(){return this.ut.et[this.lt]}set et(t){this.ut.et[this.lt]=t}get Z(){return this.ut.Z.subarray(this.offset,this.offset+this.ht)}set Z(t){this.ut.Z.set(t,this.offset)}get X(){return this.ut.X.subarray(this.offset,this.offset+this.ht)}set X(t){this.ut.X.set(t,this.offset)}get tt(){return this.ut.tt.subarray(this.offset,this.offset+this.ht)}set tt(t){this.ut.tt.set(t,this.offset)}get rt(){return this.ut.rt[this.lt]}set rt(t){this.ut.rt[this.lt]=t}get nt(){return this.ut.nt[this.lt]}set nt(t){this.ut.nt[this.lt]=t}get it(){return this.ut.it[this.lt]}set it(t){this.ut.it[this.lt]=t}get st(){return this.ut.st[this.lt]}set st(t){this.ut.st[this.lt]=t}get ot(){return this.ut.ot[this.lt]}set ot(t){this.ut.ot[this.lt]=t}}function C(t,e,r,n=0){const i=t.V,s=i.length,o=new b(s,r.dt,1,e);let h=e;const c=r.$;let l=o.ct(),u=0;for(let t=0;t<s;t++){l.ft(t),u=i[t],l.Z[u]=1,l.X[u]=h,l.L[l.et]=u,l.et=1,l.rt=h,l.st=c/2*h,l.ot=c/2*h;for(let t=0;t<r.dt;t++)l.tt[t]=h*r.W[t][u]}return 1&n||(o.st[0]/=2,o.st[s-1]/=2),2&n||(o.ot[0]/=2,o.ot[s-1]/=2),o}function E(t,e,r,n,i,s=0){const o=U(r),h=t.Y+e.Y,c=t.weight+e.weight,l=new b(o,i.dt,h,c),u=i.$;let f=l.ct(),a=t.ct(),d=e.ct(),g=_(r),w=_(n),A=0,p=0,y=1,M=1,m=!1,S=!1,q=!1,v=!1,O=0;for(let r=0;r<o;r++){if(f.ft(r),A=g[r],a.ft(A),-1===A?(m=y>=0,S=r===o-1||g[r+1]>0):(m=!1,S=!1),y=A,p=w[r],d.ft(p),-1===p?(q=M>=0,v=r==o-1||w[r+1]>0):(q=!1,v=!1),M=p,O=0,-1!==A&&-1!==p){f.nt=a.nt+d.nt,f.it=a.it+d.it,f.rt=a.rt+d.rt,f.Z.set(a.Z),f.X.set(a.X);for(let t=0;t<e.ht;t++){let r=e.Z[t];r&&(f.Z[t]+=r,f.X[t]+=e.X[t]),f.Z[t]&&(f.L[O]=t,O++),f.tt[t]=a.tt[t]+d.tt[t]}f.et=O}else-1==A?(f.nt=d.nt+(m?t.weight:0),f.it=d.it+(S?t.weight:0),f.rt=d.rt,f.Z.set(d.Z),f.X.set(d.X),f.L.set(d.L),f.tt.set(d.tt),f.et=d.et):-1==p&&(f.nt=a.nt+(q?e.weight:0),f.it=a.it+(v?e.weight:0),f.rt=a.rt,f.Z.set(a.Z),f.X.set(a.X),f.L.set(a.L),f.tt.set(a.tt),f.et=a.et);f.st=u/2*(1-f.nt)*f.rt,f.ot=u/2*(1-f.it)*f.rt}return 1&s||(l.st[0]/=2,l.st[o-1]/=2),2&s||(l.ot[0]/=2,l.ot[o-1]/=2),l}function G(t){return"gt"in t}function P(t,e){let r=1/0,n=0,i=0,s=t.length;for(let e=0;e<s;e++)for(let o=e+1;o<s;o++)t[e][o]<r&&(r=t[e][o],n=e,i=o);return[{profile:null,wt:e[n],At:e[i],yt:r,Mt:[],St:[],Ut:[],type:1,_t:0,id:"",qt:[],parent:-1,weight:0},n,i]}function R(t,e,r,n){const[i,s]=e>r?[e,r]:[r,e],o=[],h=t.length;for(var c=0;c<h;c++){if(c==e||c==r)continue;const n=.1*(t[c][e]+t[c][r])/2+.9*Math.min(t[c][e],t[c][r]);t[c].push(n),o.push(n),t[c].splice(i,1),t[c].splice(s,1)}return o.push(0),t.push(o),t.splice(i,1),t.splice(s,1),n.splice(i,1),n.splice(s,1),t}class I{constructor(t){this.vt=t,this.store=Array(t),this.head=0,this.bt=0}get size(){return(this.bt-this.head+this.vt)%this.vt}get Ot(){return this.head===this.bt}Ct(){if(this.Ot)return null;const t=this.store[this.head];return this.head++,this.head>=this.vt&&(this.head=0),t}Et(){return this.Ot?null:(this.bt--,this.bt<0&&(this.bt=this.vt-1),this.store[this.bt])}Gt(t){this.head--,this.head<0&&(this.head=this.vt-1),this.store[this.head]=t,this.head===this.bt&&(this.bt=(this.bt-1+this.vt)%this.vt)}Pt(t){this.store[this.bt]=t,this.bt++,this.bt>=this.vt&&(this.bt=0),this.head===this.bt&&(this.head=(this.head-1+this.vt)%this.vt)}Rt(t=0){const e=0===t?this.head:(this.head+t+this.vt)%this.vt;return this.store[e]}It(t=0){const e=(this.bt-1-t+this.vt)%this.vt;return this.store[e]}}function F(t,e,r){const n=new Map,i=new I(r),s=0|e,o=new Uint16Array(t.V.length-s+1);let h=0,c=0;for(let e=0;e<s;e++)h|=t.V[s-e-1]<<2*e;o[c++]=h;for(let e=s,r=t.V.length;e<r;e++)h=(h<<2)+t.V[e],o[c++]=h;const l=r-s,u=t.Ft.length*(2/(r+1))*2|0,f={Tt:new Uint16Array(u),Dt:new Uint16Array(u),xt:new Uint16Array(u),Bt:new Uint16Array(u),count:0};let a=NaN,d=0,g=-l-1;for(let t=0;t<o.length;t++){const e=o[t];for(g++;!i.Ot&&((d=i.It())<g||o[d]>=e);)i.Et();if(i.Pt(t),g<0)continue;for(;i.Rt()<g;)i.Ct();let r=i.Rt();if(r===a)f.Bt[f.count-1]=t+s;else{const e=o[r],i=f.count;f.count++,f.Tt[i]=e,f.Dt[i]=r,f.xt[i]=g,f.Bt[i]=t+s,a=r,n.has(e)?n.get(e).push(i):n.set(e,[i])}}return[n,f,o]}function T(t,e,r,n,i){const[s,o,h]=null!=n?n:F(t,8,24),[c,l,u]=null!=i?i:F(e,8,24);t.Ft.length,e.Ft.length;const f=new Map;for(let r=0;r<o.count;r++){let n=o.Tt[r];if(!c.has(n))continue;let i=c.get(n);if(i.length>4)continue;let s=t.Ft.substring(o.xt[r],o.Bt[r]);for(let t=0;t<i.length;t++){let n=i[t],h=e.Ft.substring(l.xt[n],l.Bt[n]),c=Math.min(s.length,h.length),u=!1;u=c==s.length?0===h.indexOf(s):0===s.indexOf(h);let a,d=0;u?(d=o.xt[r]-l.xt[n],a={Nt:d,Ht:o.xt[r],end:o.xt[r]+c}):(d=o.Dt[r]-l.Dt[n],a={Nt:d,Ht:o.Dt[r],end:o.Dt[r]+8-1}),f.has(d)?f.get(d).push(a):f.set(d,[a])}}let a=[];if(f.forEach((t=>{if(!t.length)return;let e=t[0];a.push(e),t.forEach((t=>{t.Ht<=e.end?e.end=t.end:(e=t,a.push(e))}))})),0===a.length)return q(t,e,r).j;a.sort(((t,e)=>{let r=t.Ht+t.Nt-e.Ht-e.Nt;return 0===r?t.Nt-e.Nt:r})),a=function(t,e,r){var n;const i=1/Math.max(e,r),s=e-r;let o=[0],h=new Int16Array(t.length);h[0]=-1;let c=t[0],l=c.end-c.Ht,u=c.end-c.Nt,f=c.Ht-c.Nt,a=Math.min(r-u,e-c.end),d=[l],g=[u],w=[a],A=t[0],p=0,y=t[p],M=l,m=0;for(let n=1;n<t.length;n++){if(c=t[n],l=c.end-c.Ht,u=c.end-c.Nt,f=c.Ht-c.Nt,a=c.Nt<s?r-u:e-c.end,20===m){m=0;let t=[],r=o.length-1,n=g[r],i=e-c.Ht;for(;r--;)g[r]>n?t.push(r):(n=g[r],d[r]+Math.min(w[r],i)<M&&t.push(r));t.forEach((t=>{o.splice(t,1),d.splice(t,1),g.splice(t,1),w.splice(t,1)}))}if(c.Nt>=A.Nt&&u<g[0]){if(l+a<M)continue;if(h[n]=-1,A=t[n],m++,l<d[0]){o.unshift(n),d.unshift(l),g.unshift(u),w.unshift(a);continue}l>M&&(M=l),o[0]=n,d[0]=l,g[0]=u,w[0]=a}else{for(let e=o.length-1;e>=0;e--)if(p=o[e],y=t[p],c.Ht>=y.end&&g[e]<f){const t=Math.abs(c.Nt-p)*i,r=d[e]+l-t;if(r+a<M)continue;let s=e;for(;s<o.length&&d[s]<r;)s++;if(s===o.length)M=r,w.push(a),d.push(r),g.push(u),o.push(n);else{if(u>=g[g.length-1])continue;w.splice(s,0,a),d.splice(s,0,r),g.splice(s,0,u),o.splice(s,0,n)}h[n]=p,m++;break}l>M&&(h[n]=-1,M=l,w.push(a),d.push(l),g.push(u),o.push(n),m++)}}let S=[],U=null!==(n=o.pop())&&void 0!==n?n:0,_=U;for(;U>-1&&_>-1;)S.push(U),U=h[U],_--;let q=[];for(let e=S.length-1;e>=0;e--)q.push(t[S[e]]);return q}(a,t.V.length,e.V.length);let d={Nt:0,Ht:0,end:0},g=[d],w=[];for(let r=0;r<a.length;r++){let n=a[r];if(n.Nt!==d.Nt){let r=d.end,i=r-d.Nt;for(;t.V[r]===e.V[i]&&r<n.Ht&&i<n.Ht-n.Nt;)r++,i++;for(d.end=r,r=n.Ht,i=r-n.Nt;t.V[r]===e.V[i]&&r>d.end&&i>d.end-d.Nt;)r--,i--;n.Ht=r+1,w.push({$t:n.Nt,kt:d.Nt,Ht:d.end,end:n.Ht}),g.push(n),d=n}else{let t=d.end,e=t-n.Nt,r=0;for(;t<n.Ht&&e<n.Ht-n.Nt&&r<2;)t+=8,e+=8,r=x(h[t],u[e])<=2?0:r+1;for(t>d.end&&(d.end=Math.max(t-2-8,d.end)),t=n.Ht-8,e=t-n.Nt,r=0;t>d.end&&e>d.end-d.Nt&&r<2;)t-=8,e-=8,r=x(h[t],u[e])<=2?0:r+1;if(t<n.Ht-16&&(n.Ht=t+24+2),n.Ht-d.end<=24){d.end=n.end;continue}w.push({kt:d.Nt,$t:n.Nt,Ht:d.end,end:n.Ht}),g.push(n),d=n}}let A=g[0];const p=[],y=[];for(let n=0;n<g.length;n++){A=g[n],p.push(A.end-A.Ht),y.push(A.end-A.Ht);let i=w[n];if(i){if(i.Ht===i.end){const t=Math.abs(i.$t-i.kt);p.push(-t),y.push(t);continue}if(i.Ht-i.kt==i.end-i.$t){const t=i.end-i.Ht;y.push(-t),p.push(t);continue}let n=q({Ft:t.Ft.substring(i.Ht,i.end),type:t.type,zt:new Uint8Array(0),V:t.V.subarray(i.Ht,i.end)},{Ft:e.Ft.substring(i.Ht-i.kt,i.end-i.$t),type:e.type,zt:new Uint8Array(0),V:e.V.subarray(i.Ht-i.kt,i.end-i.$t)},r,3);p.push(...n.j[0]),y.push(...n.j[1])}}D(p,t.V.length),D(y,e.V.length);let M=S(p,1),_=S(y,1);const v=U(M),b=U(_);return v>b?_=m(_,[b-v]):v<b&&(M=m(M,[v-b])),[M,_]}function D(t,e){let r=e-function(t){let e=0,r=0;for(let n=0;n<t.length;n++)r=t[n],r<0||(e+=r);return e}(t);r>0?t.push(r):r<0&&(t[t.length-1]+=r)}function x(e,r){let n=e^r;return n=1431655765&n|n>>1&1431655765,t(n)}const B={Vt:"auto",gapchar:"-",Wt:"auto",debug:!1},N=new class{constructor(){this.jt=[],this.Kt=2,this.Yt=Object.assign({},B)}Jt(t){if(Array.isArray(t))t.forEach((t=>{this.Jt(t)}));else{if("string"!=typeof t)throw new TypeError("String type expected for sequences to add.");if(t=t.toUpperCase(),"auto"===this.Yt.Vt){const e=l(t);if(2===this.Kt)this.Kt=e;else if(this.Kt!==e)throw Error("All sequences must be of same type.")}this.jt.push(function(t,e){const r=null!=e?e:l(t),n=function(t,e){const r=new Uint8Array(t.length),n=0===e?f:a;let i=0;for(let s=0,o=t.length;s<o;s++)i=n[t.charCodeAt(s)],255===i&&(i=d(t[s],e)),r[s]=i;return r}(t,r);return{Ft:t,V:n,zt:0===r?g(n):n,type:r}}(t,this.Kt))}}reset(){this.jt=[],this.Kt=2,this.Lt()}Lt(){this.Yt=Object.assign({},B)}Qt(t){var e;function r(t,e){return e.some((e=>e==t))}t.method&&r(t.method,["complete","diag"])&&(this.Yt.Wt=t.method),t.type&&r(t.type,["amino","nucleic"])&&(this.Yt.Vt=t.type,this.Kt="amino"===t.type?0:1),1===(null===(e=t.gapchar)||void 0===e?void 0:e.length)&&(this.Yt.gapchar=t.gapchar),void 0!==t.debug&&(this.Yt.debug=!!t.debug)}align(t,i){const s=new Promise(((s,o)=>{if(!Array.isArray(t)||t.some((t=>"string"!=typeof t)))return o("Array of sequences expected");if(t.length<2)return o("At least 2 sequences are required.");this.reset(),i&&this.Qt(i),this.Jt(t);const h=function(t,e){var i,s;const o=0===t?r:n;let h=null!==(i=null==e?void 0:e.matrix)&&void 0!==i?i:o.matrix,c=(null==e?void 0:e.gapextend)?2*-e.gapextend:o.H;var l;return{type:t,W:(l=c,h.map((t=>t.map((t=>t+l))))),$:null!==(s=null==e?void 0:e.gapopen)&&void 0!==s?s:o.$,dt:0===t?20:4}}(this.Kt,i);let c,l=!1;l="auto"===this.Yt.Wt?this.jt.some((t=>t.Ft.length>1600)):"diag"===this.Yt.Wt,c=2==this.jt.length?l?T(this.jt[0],this.jt[1],h):q(this.jt[0],this.jt[1],h).j:l?function(t,e){let r=t.map((t=>F(t,8,24))),n=[],i=[],s=[];for(let o=0;o<r.length;o++){if(0===o)continue;let h=T(t[0],t[o],e,r[0],r[o]);i=0===i.length?h[0]:y(i,h[0]),s[o]=h[0],n[o]=h[1]}for(let t=0;t<r.length;t++){if(0===t){n[t]=i;continue}let e=M(i,s[t]);if(!e)throw console.error("An error occured while computing the center diff. for sequence #"+t,i,s[t]),console.dir(n),console.dir(s),new RangeError;n[t]=p(e,n[t])}return n}(this.jt,h):function(t,r){const n=(t,e)=>{const s=i[t.wt],o=i[t.At];let h={};if(G(s)||0!==s.St.length||n(s,e),G(o)||0!==o.St.length||n(o,e),G(s))if(G(o)){const e=q(s.gt,o.gt,r);h.K=e.K,t.Ut=[s.Ut[0],o.Ut[0]],t.Mt=e.j,t.profile=E(s.profile,o.profile,e.j[0],e.j[1],r)}else{o.qt=o.Ut.map((t=>e[t]));const n=function(t,e,r,n=0){const i=e.gt,s=i.Ft.length,o=t.profile.length;let h=0,c=[],l=[];const u=new Uint8Array(Math.ceil((o+1)*(s+1)/2));let f,a,d=0,g=0,w=0,A=0,p=0,y=0,M=0,m=0,U=0,_=0,q=0,b=0,O=new Float32Array(r.dt);const C=r.$,E=C/2,G=C/2,P=1&n?0:-C/4,R=2&n?0:-C/4,I=i.V,F=t.profile;c[0]=0,l[0]=-1/0;for(let t=1;t<=s;t++)c[t]=E+G+P,l[t]=-1/0;for(let t=1;t<=o;t++){q=F.st[t-1],b=F.ot[t-1],O=F.tt.subarray((t-1)*r.dt,t*r.dt),m=P,y=-1/0,l[0]=F.st[0];for(let e=1;e<=s;e++)p=0,d=c[e]+q,w=l[e],e===s&&(d-=q/2),d>=w?l[e]=d:p+=1,g=m+E,A=M=y,t===o&&(g+=R),g>=A?y=g:p+=2,U=l[e]+b,_=M+G,h=c[e-1]+O[I[e-1]],e===s&&(_+=R),c[e-1]=m,h>=_?h>=U?m=h:(m=U,p+=4):_>=U?(m=_,p+=8):(m=U,p+=4),f=t*s+e,a=f%2,f>>>=1,u[f]+=a?p:p<<4;c[s]=m}const T=Math.max(h,y,l[l.length-2]),[D,x]=v(u,o,s,p);return{j:[S(D),S(x)],K:T}}(o,s,r);h.K=n.K,t.Ut=[s.Ut[0],...o.Ut],t.Mt=[n.j[1],...o.Mt.map((t=>p(n.j[0],t)))],t.profile=E(s.profile,o.profile,n.j[1],n.j[0],r)}else if(!G(o)){o.qt=o.Ut.map((t=>e[t])),s.qt=s.Ut.map((t=>e[t]));const n=function(t,e,r,n=0){const i=t.profile.length,s=e.profile.length,o=new Uint8Array(Math.ceil((i+1)*(s+1)/2));let h,c,l=[],u=[],f=0,a=0,d=0,g=0,w=0,A=0,p=0,y=0,M=0,m=0,U=0,_=0,q=0,b=new Float32Array(r.dt),O=new Uint8Array(r.dt),C=0,E=0,G=0,P=0,R=0,I=0;const F=1&n?1:2,T=e.profile,D=t.profile;for(l[0]=0,u[0]=-1/0,R=1;R<=s;R++)l[R]=T.st[0]+T.ot[R-1]/F,u[R]=-1/0;const x=D.st[0];for(P=1;P<=i;P++){for(b=D.tt.subarray(r.dt*(P-1),r.dt*P),_=D.st[P-1],q=D.ot[P-1],M=x+q/F,p=-1/0,u[0]=D.st[0],R=1;R<=s;R++){for(A=0,a=l[R]+_,g=u[R],R!==s||2&n||(a-=_/2),a>=g?u[R]=a:A+=1,d=M+T.st[R-1],w=y=p,P!==i||2&n||(d-=T.st[R-1]/2),d>=w?p=d:A+=2,f=l[R-1],G=T.et[R-1],I=0,E=(R-1)*r.dt,O=T.L.subarray(E,E+G);I<G;)C=O[I],f+=T.X[E+C]*b[C],I++;m=u[R]+q,U=y+T.ot[R-1],1!==R||2&n||(m-=q/2),R!==s||2&n||(m-=q/2,U-=T.ot[s-1]/2),l[R-1]=M,f>=U?f>=m?M=f:(M=m,A+=4):U>=m?(M=U,A+=8):(M=m,A+=4),h=P*s+R,c=h%2,h>>>=1,o[h]+=c?A:A<<4}l[s]=M}const B=Math.max(f,p,u[R-1]),[N,H]=v(o,i,s,A);return{j:[S(N),S(H)],K:B}}(s,o,r);h.K=n.K,t.Mt=[...s.Mt.map((t=>p(n.j[0],t))),...o.Mt.map((t=>p(n.j[1],t)))],t.profile=E(s.profile,o.profile,n.j[0],n.j[1],r)}return t.Ut=[...s.Ut,...o.Ut],h.K},i=function(t,e){let r,n,i,s=[],o=t.length,h=[];for(var c=0;c<o;c++)s[c]={type:0,gt:e[c],profile:null,wt:c,At:c,yt:0,Ut:[c],St:[],id:c.toString(),weight:0,parent:-1,_t:0},h[c]=c;c=0;for(var l=o-1;c<l;c++){[r,n,i]=P(t,h);const e=s[r.wt],l=s[r.At];r.id=e.id<l.id?`|${e.id},${l.id}|`:`|${l.id},${e.id}|`,r._t=Math.max(e._t,l._t)+1,l.parent=e.parent=s.length,h.push(o+c),s.push(r),t=R(t,n,i,h)}return s[s.length-1].type=2,s}(function(t){const r=0===t[0].type,n=r?6:4,i=Math.pow(n,6),s=t.map((t=>{const s=new e(i),o=r?t.zt:t.V;let h=0;for(let t=0;t<6;t++)h+=o[t]*Math.pow(n,t);s.set(h);const c=Math.pow(n,5);for(let t=6,e=o.length;t<e;t++)h-=o[t-6],h/=n,h+=o[t]*c,s.set(h);return s})),o=t.length,h=t.map((()=>[]));let c,l,u,f,a,d,g;for(let t=0;t<o;t++){h[t][t]=0,c=s[t],u=c.P(),d=1/(-Math.LN2/Math.log(1-u/i)*2);for(let e=t+1;e<o;e++)f=s[e].P(),a=c.D(s[e]),a+=i-(u+f-a),g=Math.ceil(d*f),a-=g,a=Math.max(a,0),l=1-a/i,h[e][t]=h[t][e]=l}return h}(t),t),s=i[i.length-1],o=function(t){let e=0;return function r(n){var i=0,s=0;switch(n.type){case 2:n.weight=0,r(t[n.wt]),r(t[n.At]);break;case 1:i=t[n.parent].yt-n.yt,s=n.id.split(",").length,n.weight=t[n.parent].weight+i/s,r(t[n.wt]),r(t[n.At]);break;case 0:i=t[n.parent].yt,n.weight=t[n.parent].weight+i,e+=n.weight}}(t[t.length-1]),t.map((t=>t.weight))}(i);return function(t,e,r){let n=0,i=t[n];for(;G(i)&&n<t.length;)i.profile=C(i.gt,i.weight,e,void 0),n++,i=t[n]}(i,r),n(s,o),function(t,e){let r=e.slice();return e.forEach(((t,e)=>r[t]=e)),r.map((e=>t[e]))}(s.Mt.slice(),s.Ut)}(this.jt,h);const u=function(t,e,r){var n;let i=null!==(n=null==r?void 0:r.gapchar)&&void 0!==n?n:"-",s=[];for(let r=0;r<t.length;r++)s.push(A(t[r].Ft,e[r],{gapchar:i}));return s}(this.jt,c,{gapchar:this.Yt.gapchar});return s(u)})).catch((t=>Promise.reject(t)));return s}};export{N as default};
//# sourceMappingURL=biomsa.esm.js.map
