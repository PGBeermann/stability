'use strict';

// Masas monoisotópicas (u), masa nominal, electrones de valencia y valencias
// habituales de los elementos soportados.
const ELEMENTOS = {
  H:  { masa: 1.00782503207,  nominal: 1,   eVal: 1, valencias: [1] },
  B:  { masa: 11.0093054,     nominal: 11,  eVal: 3, valencias: [3] },
  C:  { masa: 12.0,           nominal: 12,  eVal: 4, valencias: [4] },
  N:  { masa: 14.0030740048,  nominal: 14,  eVal: 5, valencias: [3, 5] },
  O:  { masa: 15.99491461956, nominal: 16,  eVal: 6, valencias: [2] },
  F:  { masa: 18.99840322,    nominal: 19,  eVal: 7, valencias: [1] },
  Si: { masa: 27.9769265325,  nominal: 28,  eVal: 4, valencias: [4] },
  P:  { masa: 30.97376163,    nominal: 31,  eVal: 5, valencias: [3, 5] },
  S:  { masa: 31.972071,      nominal: 32,  eVal: 6, valencias: [2, 4, 6] },
  Cl: { masa: 34.96885268,    nominal: 35,  eVal: 7, valencias: [1] },
  Br: { masa: 78.9183371,     nominal: 79,  eVal: 7, valencias: [1] },
  I:  { masa: 126.904473,     nominal: 127, eVal: 7, valencias: [1] },
};

const MASA_ELECTRON = 0.00054857990946;

module.exports = { ELEMENTOS, MASA_ELECTRON };
