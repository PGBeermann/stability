'use strict';

const test = require('node:test');
const assert = require('node:assert/strict');
const { analizarFragmento, parsearFormula } = require('../src');

const hia = (formula, smiles) => analizarFragmento({ formula, smiles }).hia;

test('fórmula: m/z, RDB y tipo de ion', () => {
  const f = parsearFormula('C7H7+');
  assert.equal(f.tipo, 'EE');
  assert.equal(f.mzNominal, 91);
  assert.ok(Math.abs(f.mz - 91.0542) < 0.0005);
  assert.equal(f.rdb, 4.5);
  assert.equal(parsearFormula('C6H6+•').tipo, 'OE');
  assert.equal(parsearFormula('(CH3)3C+').hill, 'C4H9');
});

test('orden de carbocationes alquílicos: 3° > 2° > 1° > metilo', () => {
  const metilo = hia('CH3+', '[CH3+]');
  const etilo = hia('C2H5+', 'C[CH2+]');
  const isopropilo = hia('C3H7+', 'C[CH+]C');
  const tbutilo = hia('C4H9+', 'C[C+](C)C');
  assert.ok(tbutilo < isopropilo && isopropilo < etilo && etilo < metilo);
  assert.ok(Math.abs(metilo - 312) < 1);
  assert.ok(Math.abs(tbutilo - 233) < 5);
});

test('aromaticidad: tropilio más estable que bencilo; ciclopentadienilo antiaromático', () => {
  assert.ok(hia('C7H7+', '[cH+]1cccccc1') < hia('C7H7+', '[CH2+]c1ccccc1'));
  const cp = analizarFragmento({ formula: 'C5H5+', smiles: '[CH+]1C=CC=C1' });
  assert.match(cp.clase, /antiaromático/);
});

test('cationes estabilizados por heteroátomos y acilio', () => {
  const iminio = analizarFragmento({ formula: 'CH4N+', smiles: 'C=[NH2+]' });
  assert.match(iminio.clase, /iminio/);
  assert.ok(iminio.hia < hia('CH3O+', 'C=[OH+]'));
  const acetilo = analizarFragmento({ formula: 'C2H3O+', smiles: 'CC#[O+]' });
  assert.equal(acetilo.modo, 'acilio');
  assert.ok(Math.abs(acetilo.hia - 231) < 5);
});

test('fenilo y vinilo son poco estables', () => {
  assert.ok(hia('C6H5+', '[c+]1ccccc1') > 280);
  assert.ok(hia('C2H3+', 'C=[CH+]') > 280);
});

test('n-propilo se reordena a isopropilo por migración de hidruro', () => {
  const r = analizarFragmento({ formula: 'C3H7+', smiles: 'CC[CH2+]' });
  assert.ok(r.reordenamientos.length > 0);
  assert.equal(r.reordenamientos[0].smiles, 'C[CH+]C');
});

test('iones radical: benceno estable, alcohol inestable, McLafferty', () => {
  assert.ok(analizarFragmento({ formula: 'C6H6+•', smiles: 'c1ccccc1' }).indice >= 80);
  assert.ok(analizarFragmento({ formula: 'C2H6O+•', smiles: 'CCO' }).indice < 30);
  const cetona = analizarFragmento({ formula: 'C5H10O+•', smiles: 'CCCC(C)=O' });
  assert.ok(cetona.fragmentaciones.some((t) => /McLafferty posible/.test(t)));
});

test('detecta fórmula y estructura incompatibles', () => {
  const r = analizarFragmento({ formula: 'C7H9+', smiles: '[CH2+]c1ccccc1' });
  assert.equal(r.estructura.coincideConFormula, false);
  assert.throws(() => analizarFragmento({ formula: 'C7H7+', smiles: 'Cc1ccccc1' }), /neutra/);
});
