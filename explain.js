function explainChar(char) {
  const BASE_TABLE = ['A', 'C', 'G', 'T'];
  const byte   = char.charCodeAt(0);
  const binary = byte.toString(2).padStart(8, '0');
  const diads  = [binary.slice(0,2), binary.slice(2,4), binary.slice(4,6), binary.slice(6,8)];
  const bases  = diads.map(d => BASE_TABLE[parseInt(d, 2)]);
  return { char, byte, binary, diads, bases };
}