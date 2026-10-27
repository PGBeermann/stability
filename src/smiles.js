'use strict';

/*
 * Lector de SMILES sencillo (subconjunto de OpenSMILES) pensado para
 * fragmentos iónicos pequeños: átomos orgánicos y entre corchetes, cargas,
 * H explícitos, ramas, cierres de anillo y aromaticidad (con kekulización).
 */

const { ELEMENTOS } = require('./elementos');

const ORGANICOS = new Set(['B', 'C', 'N', 'O', 'P', 'S', 'F', 'I', 'Cl', 'Br']);
const AROMATICOS = new Set(['b', 'c', 'n', 'o', 'p', 's']);
const SEGUNDO_PERIODO = new Set(['B', 'C', 'N', 'O', 'F']);
const MAX_TAM_ANILLO = 8;

function parsearCorchete(txt, pos) {
  const m = txt.match(/^(\d+)?(se|as|[bcnops]|[A-Z][a-z]?)(@{1,2}(?:TH|AL|SP|TB|OH)?\d*)?(H\d*)?([+-]+\d*)?(:\d+)?$/);
  if (!m) throw new Error(`Átomo entre corchetes no válido: "[${txt}]" (posición ${pos + 1}).`);
  const simbolo = m[2];
  const aromatico = simbolo[0] === simbolo[0].toLowerCase();
  const el = aromatico ? simbolo[0].toUpperCase() + simbolo.slice(1) : simbolo;
  if (el === 'H') {
    throw new Error('Escriba los hidrógenos dentro del corchete del átomo pesado, p. ej. [CH3+] en lugar de [H].');
  }
  if (!ELEMENTOS[el]) throw new Error(`Elemento no soportado en el SMILES: "${el}".`);
  const h = m[4] ? (m[4].length > 1 ? parseInt(m[4].slice(1), 10) : 1) : 0;
  let carga = 0;
  if (m[5]) {
    const signo = m[5][0] === '+' ? 1 : -1;
    const dig = m[5].match(/\d+$/);
    carga = dig ? signo * parseInt(dig[0], 10) : signo * m[5].length;
  }
  return { el, aromatico, carga, hExplicito: h };
}

function parsearSmiles(smiles) {
  const s = String(smiles || '').trim();
  if (!s) throw new Error('El SMILES está vacío.');
  const mol = { atoms: [], bonds: [] };
  const pila = [];
  const anillos = new Map();
  let prev = null;
  let enlacePend = null;
  let i = 0;

  const agregarEnlace = (a, b, simbolo) => {
    if (a === b) throw new Error('Un átomo no puede enlazarse consigo mismo (revise los cierres de anillo).');
    if (mol.bonds.some((e) => (e.a === a && e.b === b) || (e.a === b && e.b === a))) {
      throw new Error(`Enlace duplicado entre los átomos ${a + 1} y ${b + 1}.`);
    }
    let orden;
    switch (simbolo) {
      case '=': orden = 2; break;
      case '#': orden = 3; break;
      case ':': orden = 1.5; break;
      case '-': case '/': case '\\': orden = 1; break;
      default: orden = mol.atoms[a].aromatico && mol.atoms[b].aromatico ? 1.5 : 1;
    }
    mol.bonds.push({ id: mol.bonds.length, a, b, orden, aromatico: orden === 1.5 });
  };

  while (i < s.length) {
    const c = s[i];
    if (c === '(') {
      if (prev === null) throw new Error('Una rama "(" no puede abrir el SMILES.');
      pila.push(prev);
      i++;
      continue;
    }
    if (c === ')') {
      if (!pila.length) throw new Error('Paréntesis ")" sin abrir en el SMILES.');
      prev = pila.pop();
      i++;
      continue;
    }
    if ('-=#:/\\'.includes(c)) {
      enlacePend = c;
      i++;
      continue;
    }
    if (c === '$') throw new Error('Los enlaces cuádruples ($) no están soportados.');
    if (c === '.') {
      throw new Error('El SMILES contiene fragmentos desconectados ("."): describa un solo fragmento iónico.');
    }
    if (/\d/.test(c) || c === '%') {
      let num;
      if (c === '%') {
        num = s.substr(i + 1, 2);
        if (!/^\d\d$/.test(num)) throw new Error('Cierre de anillo "%" debe ir seguido de dos dígitos.');
        i += 3;
      } else {
        num = c;
        i++;
      }
      if (prev === null) throw new Error('Un cierre de anillo no puede ir antes del primer átomo.');
      if (anillos.has(num)) {
        const abierto = anillos.get(num);
        anillos.delete(num);
        agregarEnlace(abierto.atomo, prev, enlacePend || abierto.enlace);
      } else {
        anillos.set(num, { atomo: prev, enlace: enlacePend });
      }
      enlacePend = null;
      continue;
    }

    let atomo;
    if (c === '[') {
      const j = s.indexOf(']', i);
      if (j < 0) throw new Error('Falta "]" en el SMILES.');
      atomo = parsearCorchete(s.slice(i + 1, j), i);
      i = j + 1;
    } else {
      const dos = s.substr(i, 2);
      if (dos === 'Cl' || dos === 'Br') {
        atomo = { el: dos, aromatico: false };
        i += 2;
      } else if (ORGANICOS.has(c)) {
        atomo = { el: c, aromatico: false };
        i++;
      } else if (AROMATICOS.has(c)) {
        atomo = { el: c.toUpperCase(), aromatico: true };
        i++;
      } else {
        throw new Error(`Carácter no reconocido en el SMILES: "${c}" (posición ${i + 1}).`);
      }
      atomo.carga = 0;
      atomo.hExplicito = null;
    }
    atomo.idx = mol.atoms.length;
    mol.atoms.push(atomo);
    if (prev !== null) agregarEnlace(prev, atomo.idx, enlacePend);
    enlacePend = null;
    prev = atomo.idx;
  }
  if (anillos.size) throw new Error(`Cierre(s) de anillo sin pareja: ${[...anillos.keys()].join(', ')}.`);
  if (pila.length) throw new Error('Falta cerrar un paréntesis en el SMILES.');
  return mol;
}

function construirAdyacencia(mol) {
  mol.adj = mol.atoms.map(() => []);
  for (const b of mol.bonds) {
    mol.adj[b.a].push({ vecino: b.b, enlace: b });
    mol.adj[b.b].push({ vecino: b.a, enlace: b });
  }
}

function contarHImplicitos(mol) {
  for (const a of mol.atoms) {
    if (a.hExplicito !== null) {
      a.h = a.hExplicito;
      continue;
    }
    const val = ELEMENTOS[a.el].valencias;
    if (a.aromatico) {
      if (a.el === 'O' || a.el === 'S') {
        a.h = 0;
        continue;
      }
      let suma = 1;
      for (const x of mol.adj[a.idx]) suma += x.enlace.aromatico ? 1 : x.enlace.orden;
      a.h = Math.max(0, val[0] - suma);
    } else {
      let suma = 0;
      for (const x of mol.adj[a.idx]) suma += x.enlace.aromatico ? 1 : x.enlace.orden;
      const v = val.find((x) => x >= suma);
      a.h = v === undefined ? 0 : v - suma;
    }
  }
}

function valenciaObjetivo(el, q) {
  switch (el) {
    case 'C': case 'Si': return 4 - Math.abs(q);
    case 'N': case 'P': return 3 + q;
    case 'O': case 'S': return 2 + q;
    case 'B': return 3 - q;
    default: return 1;
  }
}

/** Convierte los enlaces aromáticos en una estructura de Kekulé (simples/dobles). */
function kekulizar(mol) {
  const aromaticos = mol.bonds.filter((b) => b.aromatico);
  if (!aromaticos.length) return;
  const necesita = mol.atoms.map((a) => {
    if (!a.aromatico) return false;
    let sigma = a.h;
    for (const x of mol.adj[a.idx]) sigma += x.enlace.aromatico ? 1 : x.enlace.orden;
    return valenciaObjetivo(a.el, a.carga) - sigma >= 1;
  });
  const pareja = mol.atoms.map(() => -1);
  const doble = new Set();
  const pendientes = mol.atoms.filter((a) => necesita[a.idx]).map((a) => a.idx);
  const candidatos = (i) => mol.adj[i].filter((x) => x.enlace.aromatico && necesita[x.vecino] && pareja[x.vecino] === -1);

  const resolver = () => {
    let mejor = -1;
    let mejorC = null;
    for (const i of pendientes) {
      if (pareja[i] !== -1) continue;
      const c = candidatos(i);
      if (mejor === -1 || c.length < mejorC.length) {
        mejor = i;
        mejorC = c;
        if (!c.length) break;
      }
    }
    if (mejor === -1) return true;
    for (const x of mejorC) {
      pareja[mejor] = x.vecino;
      pareja[x.vecino] = mejor;
      doble.add(x.enlace.id);
      if (resolver()) return true;
      pareja[mejor] = -1;
      pareja[x.vecino] = -1;
      doble.delete(x.enlace.id);
    }
    return false;
  };
  if (!resolver()) {
    throw new Error('No se pudo asignar una estructura de Kekulé al sistema aromático. Revise los H explícitos (p. ej. [nH] en el pirrol) y las cargas.');
  }
  for (const b of aromaticos) {
    b.orden = doble.has(b.id) ? 2 : 1;
    b.aromatico = false;
  }
}

function calcularPropiedades(mol) {
  for (const a of mol.atoms) {
    let suma = a.h;
    for (const x of mol.adj[a.idx]) suma += x.enlace.orden;
    const nb = ELEMENTOS[a.el].eVal - a.carga - suma;
    if (nb < 0) {
      throw new Error(`Valencia imposible en el átomo ${a.idx + 1} (${a.el}, carga ${a.carga >= 0 ? '+' : ''}${a.carga}): ${suma} enlaces.`);
    }
    const capa = 2 * suma + nb;
    if (SEGUNDO_PERIODO.has(a.el) && capa > 8) {
      throw new Error(`El átomo ${a.idx + 1} (${a.el}) supera el octeto: revise enlaces, H y carga.`);
    }
    a.sumaEnlaces = suma;
    a.noEnlazantes = nb;
    a.paresLibres = Math.floor(nb / 2);
    a.radical = nb % 2;
    a.electronesCapa = capa;
    a.sexteto = SEGUNDO_PERIODO.has(a.el) && capa === 6 && a.radical === 0;
    a.tieneMultiple = mol.adj[a.idx].some((x) => x.enlace.orden >= 2);
  }
}

function verificarConexion(mol) {
  if (!mol.atoms.length) return;
  const visto = new Set([0]);
  const cola = [0];
  while (cola.length) {
    const u = cola.shift();
    for (const { vecino } of mol.adj[u]) {
      if (!visto.has(vecino)) {
        visto.add(vecino);
        cola.push(vecino);
      }
    }
  }
  if (visto.size !== mol.atoms.length) throw new Error('La estructura no está conectada.');
}

/** Ciclos simples de hasta MAX_TAM_ANILLO átomos. */
function encontrarCiclos(mol) {
  const ciclos = [];
  const n = mol.atoms.length;
  for (let s = 0; s < n; s++) {
    const camino = [s];
    const enCamino = new Set([s]);
    const dfs = (u) => {
      for (const { vecino: v } of mol.adj[u]) {
        if (ciclos.length > 5000) return;
        if (v === s && camino.length >= 3) {
          if (camino[1] < camino[camino.length - 1]) ciclos.push([...camino]);
          continue;
        }
        if (v <= s || enCamino.has(v) || camino.length >= MAX_TAM_ANILLO) continue;
        camino.push(v);
        enCamino.add(v);
        dfs(v);
        camino.pop();
        enCamino.delete(v);
      }
    };
    dfs(s);
  }
  return ciclos;
}

/** Cuenta electrones π de un ciclo (regla de Hückel). */
function evaluarHuckel(mol, ciclo) {
  const enCiclo = new Set(ciclo);
  let e = 0;
  for (const i of ciclo) {
    const a = mol.atoms[i];
    const adj = mol.adj[i];
    const dobleDentro = adj.some((x) => x.enlace.orden === 2 && enCiclo.has(x.vecino));
    const multipleFuera = adj.find((x) => x.enlace.orden >= 2 && !enCiclo.has(x.vecino));
    if (dobleDentro) e += 1;
    else if (multipleFuera) e += mol.atoms[multipleFuera.vecino].el === 'C' ? 1 : 0;
    else if (a.sexteto) e += 0;
    else if (a.radical) e += 1;
    else if (a.paresLibres >= 1) e += 2;
    else return { conjugado: false, electrones: null, aromatico: false, antiaromatico: false };
  }
  return { conjugado: true, electrones: e, aromatico: e % 4 === 2, antiaromatico: e > 0 && e % 4 === 0 };
}

function recalcular(mol) {
  construirAdyacencia(mol);
  calcularPropiedades(mol);
  mol.ciclos = encontrarCiclos(mol);
  mol.huckel = mol.ciclos.map((c) => evaluarHuckel(mol, c));
  return mol;
}

function prepararMolecula(smiles) {
  const mol = parsearSmiles(smiles);
  construirAdyacencia(mol);
  contarHImplicitos(mol);
  kekulizar(mol);
  verificarConexion(mol);
  recalcular(mol);
  mol.smilesOriginal = String(smiles).trim();
  return mol;
}

function clonar(mol) {
  const m = {
    atoms: mol.atoms.map((a) => ({ idx: a.idx, el: a.el, aromatico: a.aromatico, carga: a.carga, h: a.h, hExplicito: a.h })),
    bonds: mol.bonds.map((b) => ({ id: b.id, a: b.a, b: b.b, orden: b.orden, aromatico: false })),
  };
  return recalcular(m);
}

/** ¿El átomo i está en un anillo aromático (Hückel) que no contiene al átomo `excluir`? */
function enAnilloAromatico(mol, i, excluir = -1) {
  return mol.ciclos.some((c, k) => mol.huckel[k].aromatico && c.includes(i) && !c.includes(excluir));
}

function enAnilloDeTam(mol, i, tam, excluir = -1) {
  return mol.ciclos.some((c) => c.length === tam && c.includes(i) && !c.includes(excluir));
}

function enlaceEnAnillo(mol, a, b) {
  return mol.ciclos.some((c) => {
    const ia = c.indexOf(a);
    const ib = c.indexOf(b);
    if (ia < 0 || ib < 0) return false;
    const d = Math.abs(ia - ib);
    return d === 1 || d === c.length - 1;
  });
}

/**
 * Posiciones a las que se deslocaliza la carga de un carbocatión por
 * resonancia alílica: C⁺–C=C → C=C–C⁺ (se recorren caminos alternantes).
 */
function sitiosResonancia(mol, sitio) {
  const sitios = [sitio];
  const visto = new Set([sitio]);
  const cola = [sitio];
  while (cola.length) {
    const p = cola.shift();
    for (const { vecino: a, enlace } of mol.adj[p]) {
      if (enlace.orden !== 1 || mol.atoms[a].el !== 'C') continue;
      for (const { vecino: b, enlace: e2 } of mol.adj[a]) {
        if (b === p || e2.orden !== 2 || mol.atoms[b].el !== 'C' || visto.has(b)) continue;
        visto.add(b);
        sitios.push(b);
        cola.push(b);
      }
    }
  }
  return sitios;
}

function conteoAtomos(mol) {
  const c = {};
  for (const a of mol.atoms) {
    c[a.el] = (c[a.el] || 0) + 1;
    if (a.h) c.H = (c.H || 0) + a.h;
  }
  return c;
}

function cargaTotal(mol) {
  return mol.atoms.reduce((s, a) => s + a.carga, 0);
}

function textoAtomo(mol, a) {
  let suma = 0;
  for (const x of mol.adj[a.idx]) suma += x.enlace.orden;
  const v = ELEMENTOS[a.el].valencias.find((x) => x >= suma);
  const hDefecto = v === undefined ? null : v - suma;
  if (ORGANICOS.has(a.el) && a.carga === 0 && hDefecto === a.h) return a.el;
  let t = '[' + a.el;
  if (a.h) t += 'H' + (a.h > 1 ? a.h : '');
  if (a.carga) t += (a.carga > 0 ? '+' : '-') + (Math.abs(a.carga) > 1 ? Math.abs(a.carga) : '');
  return t + ']';
}

/** Escribe un SMILES (forma de Kekulé) de la molécula. */
function escribirSmiles(mol) {
  const n = mol.atoms.length;
  if (!n) return '';
  const visitado = new Array(n).fill(false);
  const hijos = Array.from({ length: n }, () => []);
  const cierres = Array.from({ length: n }, () => []);
  const usado = new Set();
  let etiqueta = 0;
  const dfs = (u) => {
    visitado[u] = true;
    for (const { vecino: v, enlace: e } of mol.adj[u]) {
      if (usado.has(e.id)) continue;
      usado.add(e.id);
      if (visitado[v]) {
        etiqueta++;
        cierres[v].push({ etiqueta, e });
        cierres[u].push({ etiqueta, e });
      } else {
        hijos[u].push({ v, e });
        dfs(v);
      }
    }
  };
  dfs(0);
  const sim = (e) => (e.orden === 2 ? '=' : e.orden === 3 ? '#' : '');
  const lab = (k) => (k < 10 ? String(k) : '%' + k);
  const escribir = (u) => {
    let t = textoAtomo(mol, mol.atoms[u]);
    for (const c of cierres[u]) t += sim(c.e) + lab(c.etiqueta);
    hijos[u].forEach((h, k) => {
      const sub = sim(h.e) + escribir(h.v);
      t += k < hijos[u].length - 1 ? '(' + sub + ')' : sub;
    });
    return t;
  };
  return escribir(0);
}

module.exports = {
  prepararMolecula,
  recalcular,
  clonar,
  enAnilloAromatico,
  enAnilloDeTam,
  enlaceEnAnillo,
  sitiosResonancia,
  evaluarHuckel,
  conteoAtomos,
  cargaTotal,
  escribirSmiles,
};
