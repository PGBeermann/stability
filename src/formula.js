'use strict';

const { ELEMENTOS, MASA_ELECTRON } = require('./elementos');

/** Fórmula en orden de Hill: C, H y luego alfabético (sin C: todo alfabético). */
function formulaHill(atomos) {
  const els = Object.keys(atomos).filter((e) => atomos[e] > 0);
  const orden = [];
  if (els.includes('C')) {
    orden.push('C');
    if (els.includes('H')) orden.push('H');
    orden.push(...els.filter((e) => e !== 'C' && e !== 'H').sort());
  } else {
    orden.push(...els.sort());
  }
  return orden.map((e) => e + (atomos[e] > 1 ? atomos[e] : '')).join('');
}

function parsearCuerpo(cuerpo) {
  const pila = [{}];
  let i = 0;
  const leerNumero = () => {
    let num = '';
    while (i < cuerpo.length && /\d/.test(cuerpo[i])) num += cuerpo[i++];
    return num ? parseInt(num, 10) : 1;
  };
  while (i < cuerpo.length) {
    const c = cuerpo[i];
    if (c === '(' || c === '[') {
      pila.push({});
      i++;
    } else if (c === ')' || c === ']') {
      if (pila.length < 2) throw new Error(`Paréntesis desbalanceados en la fórmula "${cuerpo}".`);
      i++;
      const mult = leerNumero();
      const grupo = pila.pop();
      const tope = pila[pila.length - 1];
      for (const [el, n] of Object.entries(grupo)) tope[el] = (tope[el] || 0) + n * mult;
    } else if (/[A-Z]/.test(c)) {
      let el = c;
      i++;
      if (i < cuerpo.length && /[a-z]/.test(cuerpo[i])) {
        el += cuerpo[i];
        i++;
      }
      if (!ELEMENTOS[el]) throw new Error(`Elemento no soportado en la fórmula: "${el}".`);
      const n = leerNumero();
      const tope = pila[pila.length - 1];
      tope[el] = (tope[el] || 0) + n;
    } else {
      throw new Error(`Carácter inesperado "${c}" en la fórmula "${cuerpo}".`);
    }
  }
  if (pila.length !== 1) throw new Error(`Paréntesis desbalanceados en la fórmula "${cuerpo}".`);
  return pila[0];
}

/**
 * Interpreta fórmulas de iones como "C7H7+", "[C7H7]+", "C6H6+•", "C6H6.+",
 * "(CH3)3C+". Si no se indica carga, se asume +1 (modo positivo de EI-MS).
 */
function parsearFormula(texto) {
  if (!texto || !String(texto).trim()) throw new Error('La fórmula está vacía.');
  let s = String(texto).trim().replace(/\s+/g, '').replace(/\^/g, '');
  let fin = s.length;
  while (fin > 0 && '+-•·.*'.includes(s[fin - 1])) fin--;
  const cola = s.slice(fin);
  let cuerpo = s.slice(0, fin);
  const radicalMarcado = /[•·.*]/.test(cola);
  const mas = (cola.match(/\+/g) || []).length;
  const menos = (cola.match(/-/g) || []).length;
  let carga = mas - menos;
  const mCorchete = cuerpo.match(/^\[(.+)\](\d+)$/);
  if (mCorchete && carga !== 0) {
    cuerpo = mCorchete[1];
    carga = Math.sign(carga) * parseInt(mCorchete[2], 10);
  } else if (/^\[[^\]]+\]$/.test(cuerpo)) {
    cuerpo = cuerpo.slice(1, -1);
  }
  let cargaAsumida = false;
  if (mas === 0 && menos === 0) {
    carga = 1;
    cargaAsumida = true;
  }
  if (!cuerpo) throw new Error(`No se encontraron elementos en la fórmula "${texto}".`);

  const atomos = parsearCuerpo(cuerpo);
  const n = (el) => atomos[el] || 0;

  let masa = 0;
  let nominal = 0;
  let eVal = 0;
  for (const [el, k] of Object.entries(atomos)) {
    masa += k * ELEMENTOS[el].masa;
    nominal += k * ELEMENTOS[el].nominal;
    eVal += k * ELEMENTOS[el].eVal;
  }
  masa -= carga * MASA_ELECTRON;
  const electrones = eVal - carga;
  const tipo = Math.abs(electrones) % 2 === 1 ? 'OE' : 'EE';

  // RDB (anillos + dobles enlaces) = 1 + Σ n_i (v_i − 2) / 2
  const rdb = 1 + n('C') + n('Si')
    - (n('H') + n('F') + n('Cl') + n('Br') + n('I')) / 2
    + (n('N') + n('P') + n('B')) / 2;

  // Intensidades isotópicas aproximadas (relativas a M = 100 %)
  const m1 = 1.07 * n('C') + 0.0115 * n('H') + 0.369 * n('N') + 0.038 * n('O') + 0.79 * n('S') + 5.08 * n('Si');
  const m2 = (1.07 * n('C')) ** 2 / 200 + 0.205 * n('O') + 4.47 * n('S') + 3.35 * n('Si')
    + 31.96 * n('Cl') + 97.28 * n('Br');

  const notas = [];
  if (cargaAsumida) notas.push('No se indicó la carga en la fórmula: se asumió +1.');
  if (radicalMarcado && tipo === 'EE') {
    notas.push('La fórmula se marcó como radical (•) pero tiene un número par de electrones: se trata como ion EE⁺.');
  }
  if (!radicalMarcado && tipo === 'OE') {
    notas.push('La fórmula tiene un número impar de electrones: es un ion radical (OE⁺•).');
  }
  if (Math.abs(carga) !== 1) notas.push(`Carga ${carga}: el análisis de estabilidad está pensado para iones monocargados.`);

  const mzNominal = nominal / Math.abs(carga || 1);
  const nN = n('N');
  const reglaN = tipo === 'OE'
    ? `Ion OE⁺• con m/z nominal ${mzNominal % 2 === 0 ? 'par' : 'impar'}: ${nN} átomo(s) de N (en iones OE, masa impar ⇔ número impar de N).`
    : `Ion EE⁺ con m/z nominal ${mzNominal % 2 === 0 ? 'par' : 'impar'}: ${nN} átomo(s) de N (en iones EE la regla se invierte: masa par ⇔ número impar de N).`;

  return {
    texto: String(texto).trim(),
    atomos,
    hill: formulaHill(atomos),
    carga,
    cargaAsumida,
    radicalMarcado,
    tipo,
    electrones,
    masaExacta: masa,
    mz: masa / Math.abs(carga || 1),
    nominal,
    mzNominal,
    rdb,
    rdbEstructural: tipo === 'EE' ? rdb - 0.5 : rdb,
    isotopos: { m1, m2 },
    reglaNitrogeno: reglaN,
    notas,
  };
}

module.exports = { parsearFormula, formulaHill };
