'use strict';

const { analizarFragmento, compararFragmentos, indiceDesdeHia, CATEGORIAS } = require('./estabilidad');
const { parsearFormula, formulaHill } = require('./formula');
const { prepararMolecula, escribirSmiles } = require('./smiles');
const { formatearInforme, formatearComparacion } = require('./informe');
const ejemplos = require('./ejemplos');

module.exports = {
  analizarFragmento,
  compararFragmentos,
  parsearFormula,
  formulaHill,
  prepararMolecula,
  escribirSmiles,
  formatearInforme,
  formatearComparacion,
  indiceDesdeHia,
  CATEGORIAS,
  ejemplos,
};
