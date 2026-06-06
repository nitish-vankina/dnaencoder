const BASE_TABLE=['A','C','G','T'];
const BASE_BITS={A:0b00,C:0b01,G:0b10,T:0b11,a:0b00,c:0b01,g:0b10,t:0b11};
const MAX_HOMOPOLYMER=4;
const GC_MIN=0.40,GC_MAX=0.60;

function textToBytes(text){return new TextEncoder().encode(text)}
function bytesToText(bytes){return new TextDecoder().decode(bytes)}
function bytesToBinaryString(bytes){return Array.from(bytes).map(b=>b.toString(2).padStart(8,'0')).join(' ')}
function explainByte(byte){const binary=byte.toString(2).padStart(8,'0');const diads=[binary.slice(0,2),binary.slice(2,4),binary.slice(4,6),binary.slice(6,8)];const bases=diads.map(d=>BASE_TABLE[parseInt(d,2)]);return{byte,binary,diads,bases}}
function encodeByte(b){return BASE_TABLE[(b>>6)&0b11]+BASE_TABLE[(b>>4)&0b11]+BASE_TABLE[(b>>2)&0b11]+BASE_TABLE[b&0b11]}
function decodeQuad(q){let r=0;for(let i=0;i<4;i++){const bits=BASE_BITS[q[i]];if(bits===undefined)throw new Error(`Invalid nucleotide '${q[i]}' at position ${i}`);r=(r<<2)|bits}return r}
function encodeBytes(bytes){const out=new Array(bytes.length*4);for(let i=0;i<bytes.length;i++){const b=bytes[i];out[i*4]=BASE_TABLE[(b>>6)&0b11];out[i*4+1]=BASE_TABLE[(b>>4)&0b11];out[i*4+2]=BASE_TABLE[(b>>2)&0b11];out[i*4+3]=BASE_TABLE[b&0b11]}return out.join('')}
function decodeStrand(s){if(s.length%4!==0)throw new Error(`Strand length ${s.length} is not a multiple of 4`);const out=new Uint8Array(s.length/4);for(let i=0;i<out.length;i++)out[i]=decodeQuad(s.slice(i*4,i*4+4));return out}
function encode(text){return encodeBytes(textToBytes(text))}
function decode(s){return bytesToText(decodeStrand(s))}
function gcContent(s){if(s.length===0)return 0;let gc=0;for(const ch of s)if(ch==='G'||ch==='C'||ch==='g'||ch==='c')gc++;return gc/s.length}
function findHomopolymerRuns(s,minLength=MAX_HOMOPOLYMER+1){if(s.length===0)return[];const runs=[];let start=0,len=1;const up=s.toUpperCase();for(let i=1;i<up.length;i++){if(up[i]===up[i-1])len++;else{if(len>=minLength)runs.push({base:up[start],position:start,length:len});start=i;len=1}}if(len>=minLength)runs.push({base:up[start],position:start,length:len});return runs}
function reverseComplement(s){const comp={A:'T',T:'A',C:'G',G:'C',a:'t',t:'a',c:'g',g:'c'};return s.split('').reverse().map(ch=>comp[ch]??ch).join('')}
function meltingTemp(s){if(s.length===0)return null;const up=s.toUpperCase();for(const ch of up)if(!'ACGT'.includes(ch))return null;if(up.length<14){const at=up.split('').filter(b=>b==='A'||b==='T').length;const gc=up.split('').filter(b=>b==='G'||b==='C').length;return 2*at+4*gc}const NN={AA:{h:-7.9,s:-22.2},TT:{h:-7.9,s:-22.2},AT:{h:-7.2,s:-20.4},TA:{h:-7.2,s:-21.3},CA:{h:-8.5,s:-22.7},TG:{h:-8.5,s:-22.7},GT:{h:-8.4,s:-22.4},AC:{h:-8.4,s:-22.4},CT:{h:-7.8,s:-21.0},AG:{h:-7.8,s:-21.0},GA:{h:-8.2,s:-22.2},TC:{h:-8.2,s:-22.2},CG:{h:-10.6,s:-27.2},GC:{h:-9.8,s:-24.4},GG:{h:-8.0,s:-19.9},CC:{h:-8.0,s:-19.9}};let dh=0,ds=0;for(let i=0;i<up.length-1;i++){const pair=up[i]+up[i+1];if(NN[pair]){dh+=NN[pair].h;ds+=NN[pair].s}}const firstGC=up[0]==='G'||up[0]==='C';const lastGC=up[up.length-1]==='G'||up[up.length-1]==='C';dh+=firstGC?0.1:2.3;dh+=lastGC?0.1:2.3;ds+=firstGC?-2.8:4.1;ds+=lastGC?-2.8:4.1;const R=1.987,CT=250e-9;return(dh*1000)/(ds+R*Math.log(CT/4))-273.15}
function analyzeStrand(s){const gc=gcContent(s);const counts={A:0,C:0,G:0,T:0};for(const ch of s.toUpperCase())if(counts[ch]!==undefined)counts[ch]++;return{length:s.length,gcFraction:gc,gcOk:gc>=GC_MIN&&gc<=GC_MAX,hpRuns:findHomopolymerRuns(s),tm:meltingTemp(s.slice(0,200)),revComp:reverseComplement(s),baseCounts:counts}}

const HelixEngine={encode,decode,textToBytes,bytesToText,bytesToBinaryString,explainByte,encodeByte,decodeQuad,encodeBytes,decodeStrand,gcContent,findHomopolymerRuns,reverseComplement,meltingTemp,analyzeStrand,BASE_TABLE,BASE_BITS,MAX_HOMOPOLYMER,GC_MIN,GC_MAX};
if(typeof module!=='undefined'&&module.exports)module.exports=HelixEngine;else if(typeof window!=='undefined')window.HelixEngine=HelixEngine;
if(typeof window!=='undefined'){window.textToBytes=textToBytes;window.bytesToText=bytesToText;window.encodeBytes=encodeBytes;window.decodeStrand=decodeStrand;window.gcContent=gcContent;window.analyzeStrand=analyzeStrand;window.encode=encode}