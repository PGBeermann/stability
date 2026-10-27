'use strict';

/*
 * Modelo de estabilidad de cationes fragmento de EI-MS.
 *
 * Cationes de electrones pares (EE⁺): se estima la afinidad por hidruro en
 * fase gaseosa, HIA = ΔH(R⁺ + H⁻ → R–H) en kcal/mol. Cuanto MENOR es la HIA,
 * MÁS estable es el catión. Se parte de un catión de referencia y se suman
 * contribuciones de sustituyentes (hiperconjugación, resonancia, donación de
 * pares libres, efectos inductivos, aromaticidad, tensión de anillo).
 * Los incrementos están calibrados con valores experimentales aproximados:
 *   CH3⁺ 312 · C2H5⁺ 270 · i-C3H7⁺ 250 · t-C4H9⁺ 233 · alilo 256 ·
 *   bencilo ≈ 238 · tropilio ≈ 201 · vinilo 288 · fenilo ≈ 290 ·
 *   CH3CO⁺ ≈ 231 · CH2=NH2⁺ ≈ 218 · CH2=OH⁺ ≈ 250.
 *
 * Iones radical (OE⁺•): no tienen HIA; se estima cualitativamente la
 * abundancia esperada del ion según los grupos funcionales presentes
 * (orden empírico clásico de estabilidad de M⁺• en EI).
 */

const sm = require('./smiles');
const { parsearFormula, formulaHill } = require('./formula');
const { buscarIonConocido } = require('./iones-conocidos');

const HIA_MIN = 190;
const HIA_MAX = 320;
const D_SATURACION = 180;
const ATENUACION = [1, 0.65, 0.55, 0.5];

const BASES = {
  carbenio: { valor: 312, texto: 'catión metilo CH3⁺' },
  heteroatomo: { valor: 330, texto: 'catión con sexteto sobre heteroátomo (nitrenio/oxenio)' },
  vinilico: { valor: 288, texto: 'catión vinilo CH2=CH⁺' },
  arilo: { valor: 290, texto: 'catión fenilo C6H5⁺ (orbital vacío sp² en el plano del anillo)' },
  acilio: { valor: 260, texto: 'ion formilo HC≡O⁺ (el O completa el octeto: C≡O⁺)' },
  nitrilio: { valor: 240, texto: 'ion nitrilio HC≡NH⁺' },
};
const FACTOR_MODO = { carbenio: 1, heteroatomo: 1, vinilico: 0.6, arilo: 0, acilio: 0.6, nitrilio: 0.6 };

const CATEGORIAS = [
  { min: 80, nombre: 'MUY ESTABLE',
    ee: 'Catión fuertemente estabilizado: suele dar picos muy intensos y con frecuencia el pico base del espectro.',
    oe: 'Ion radical muy estable: ion molecular intenso, a menudo pico base.' },
  { min: 60, nombre: 'ESTABLE',
    ee: 'Catión estable: se espera un pico intenso; las rupturas que lo forman están favorecidas.',
    oe: 'Ion radical estable: M⁺• claramente visible.' },
  { min: 40, nombre: 'MODERADAMENTE ESTABLE',
    ee: 'Estabilidad intermedia: pico de intensidad media; compite con otras rutas de fragmentación.',
    oe: 'Ion radical de estabilidad media: M⁺• de intensidad moderada.' },
  { min: 20, nombre: 'POCO ESTABLE',
    ee: 'Catión poco estabilizado: pico débil; tiende a reordenarse o a seguir fragmentándose.',
    oe: 'Ion radical poco estable: M⁺• débil; fragmenta con facilidad.' },
  { min: 0, nombre: 'MUY INESTABLE',
    ee: 'Catión muy inestable: rara vez aparece como pico importante en EI-MS.',
    oe: 'Ion radical muy inestable: M⁺• muy débil o ausente.' },
];

function categoriaDe(indice) {
  return CATEGORIAS.find((c) => indice >= c.min);
}

function indiceDesdeHia(hia) {
  const x = Math.round(((HIA_MAX - hia) / (HIA_MAX - HIA_MIN)) * 100);
  return Math.max(0, Math.min(100, x));
}

const r1 = (x) => Math.round(x * 10) / 10;

function esCarbonoAcilo(mol, i) {
  const a = mol.atoms[i];
  return a.el === 'C' && mol.adj[i].some((x) => x.enlace.orden === 2 && mol.atoms[x.vecino].el === 'O');
}

/** Clasifica un sustituyente unido al centro catiónico. s: estabiliza (kcal/mol), d: desestabiliza. */
function clasificarSustituyente(mol, sitio, v) {
  const a = mol.atoms[v];
  const otros = mol.adj[v].filter((x) => x.vecino !== sitio);
  const A = mol.atoms;
  switch (a.el) {
    case 'C': {
      if (a.sexteto || a.carga > 0) return { tipo: 'cation', etiqueta: 'carbono catiónico adyacente (repulsión de cargas)', d: 40 };
      const triple = otros.find((x) => x.enlace.orden === 3);
      if (triple) {
        return A[triple.vecino].el === 'C'
          ? { tipo: 'alquinilo', etiqueta: 'alquinilo (conjugación con C≡C, propargílico)', s: 41 }
          : { tipo: 'aceptor', etiqueta: 'nitrilo –C≡N (atractor de electrones)', d: 22 };
      }
      const dobleHet = otros.find((x) => x.enlace.orden === 2 && A[x.vecino].el !== 'C');
      if (dobleHet) {
        const elx = A[dobleHet.vecino].el;
        return { tipo: 'aceptor', etiqueta: elx === 'O' ? 'carbonilo C=O (atractor de electrones)' : `grupo C=${elx} (atractor)`, d: elx === 'O' ? 18 : 12 };
      }
      if (otros.some((x) => x.enlace.orden === 2)) {
        if (sm.enAnilloAromatico(mol, v, sitio)) {
          return { tipo: 'arilo', etiqueta: 'arilo (resonancia con el anillo aromático, bencílico)', s: 74 };
        }
        return { tipo: 'vinilo', etiqueta: 'vinilo (resonancia con C=C, alílico)', s: 56 };
      }
      if (sm.enAnilloDeTam(mol, v, 3, sitio)) {
        return { tipo: 'ciclopropilo', etiqueta: 'ciclopropilo (conjugación σ, ciclopropilcarbinilo)', s: 62 };
      }
      const electroneg = otros.filter((x) => ['F', 'Cl', 'Br', 'O', 'N'].includes(A[x.vecino].el)).length;
      return {
        tipo: 'alquilo',
        etiqueta: `alquilo (hiperconjugación: ${a.h} C–H y ${otros.length} enlace(s) C–C/C–X)`,
        s: 42,
        dExtra: electroneg * 5,
        nElectroneg: electroneg,
      };
    }
    case 'N': {
      if (a.carga > 0 || a.paresLibres === 0) return { tipo: 'aceptor', etiqueta: 'N sin par libre (atractor)', d: 25 };
      if (otros.some((x) => x.enlace.orden >= 2)) return { tipo: 'donador-N', etiqueta: 'N de imina (donación débil del par libre)', s: 30 };
      if (otros.some((x) => esCarbonoAcilo(mol, x.vecino))) return { tipo: 'donador-N', etiqueta: 'N de amida (par libre compartido con C=O)', s: 60 };
      if (sm.enAnilloAromatico(mol, v, sitio)) return { tipo: 'donador-N', etiqueta: 'N aromático (par libre parte del sexteto)', s: 50 };
      if (otros.some((x) => sm.enAnilloAromatico(mol, x.vecino, sitio))) return { tipo: 'donador-N', etiqueta: 'N–arilo (par libre parcialmente deslocalizado)', s: 85 };
      const nC = Math.min(otros.length, 2);
      return {
        tipo: 'donador-N',
        etiqueta: `amino ${['–NH2', '–NHR', '–NR2'][nC]} (donación del par libre → forma iminio C=N⁺)`,
        s: [94, 100, 105][nC],
      };
    }
    case 'O': {
      if (a.carga > 0 || a.paresLibres === 0) return { tipo: 'aceptor', etiqueta: 'O sin par libre disponible (atractor)', d: 25 };
      if (otros.some((x) => esCarbonoAcilo(mol, x.vecino))) return { tipo: 'donador-O', etiqueta: 'O de éster/ácido (par libre compartido con C=O)', s: 40 };
      if (sm.enAnilloAromatico(mol, v, sitio)) return { tipo: 'donador-O', etiqueta: 'O aromático (furano)', s: 45 };
      if (otros.some((x) => sm.enAnilloAromatico(mol, x.vecino, sitio))) return { tipo: 'donador-O', etiqueta: 'O–arilo (fenoxilo)', s: 60 };
      if (a.h >= 1) return { tipo: 'donador-O', etiqueta: 'hidroxilo –OH (donación del par libre → forma oxonio C=OH⁺)', s: 62 };
      return { tipo: 'donador-O', etiqueta: 'alcoxilo –OR (donación del par libre → ion oxocarbenio)', s: 70 };
    }
    case 'S': {
      if (a.carga > 0 || a.sumaEnlaces > 2) return { tipo: 'aceptor', etiqueta: 'S oxidado (atractor)', d: 20 };
      return a.h >= 1
        ? { tipo: 'donador-S', etiqueta: 'sulfanilo –SH (donación del par libre)', s: 65 }
        : { tipo: 'donador-S', etiqueta: 'alquiltio –SR (donación del par libre)', s: 72 };
    }
    case 'F': return { tipo: 'halogeno', etiqueta: 'flúor (donación π débil, atracción σ fuerte)', s: 22 };
    case 'Cl': return { tipo: 'halogeno', etiqueta: 'cloro (donación π del par libre)', s: 28 };
    case 'Br': return { tipo: 'halogeno', etiqueta: 'bromo (donación π del par libre)', s: 28 };
    case 'I': return { tipo: 'halogeno', etiqueta: 'yodo (donación π del par libre)', s: 28 };
    case 'Si': return { tipo: 'sililo', etiqueta: 'sililo (efecto α moderado)', s: 30 };
    case 'P': return a.paresLibres > 0
      ? { tipo: 'donador-P', etiqueta: 'fósforo con par libre', s: 40 }
      : { tipo: 'aceptor', etiqueta: 'fósforo sin par libre (atractor)', d: 15 };
    case 'B': return { tipo: 'aceptor', etiqueta: 'boro (orbital vacío, no dona)', d: 10 };
    default: return { tipo: 'otro', etiqueta: `sustituyente ${a.el}` };
  }
}

/** Suma estabilizaciones con atenuación (los efectos no son aditivos: se saturan). */
function combinarEstabilizaciones(lista) {
  const ord = [...lista].sort((x, y) => y.s - x.s);
  let acum = 0;
  ord.forEach((it, k) => {
    const f = k === 0 ? 1 : ATENUACION[Math.min(k, ATENUACION.length - 1)] * (1 - acum / D_SATURACION);
    it.efectivo = it.s * Math.max(f, 0.1);
    acum += it.efectivo;
  });
  return ord;
}

function nombrarClase(modo, info) {
  const { tipos, anilloAromatico, anilloAnti, grado, formaHet } = info;
  if (modo === 'acilio') return 'ion acilio (R–C≡O⁺)';
  if (modo === 'nitrilio') return 'ion nitrilio (R–C≡N⁺–R)';
  if (modo === 'arilo') return 'catión arilo (tipo fenilo)';
  if (modo === 'vinilico') return 'catión vinílico';
  if (modo === 'heteroatomo') return 'catión con sexteto en heteroátomo (nitrenio/oxenio)';
  if (anilloAromatico) return 'catión aromático (Hückel 4n+2, tipo tropilio/ciclopropenilo)';
  if (anilloAnti) return 'catión antiaromático (Hückel 4n)';
  const gTxt = ['metilo', 'primario', 'secundario', 'terciario'][Math.min(grado, 3)];
  if (tipos.includes('donador-N')) return formaHet ? 'ion iminio (C=N⁺)' : 'α-aminocarbocatión ⇄ ion iminio (C=N⁺)';
  if (tipos.includes('donador-O')) return 'ion oxocarbenio / oxonio (C=O⁺)';
  if (tipos.includes('donador-S')) return 'ion tiocarbenio (C=S⁺)';
  if (tipos.includes('arilo')) return `carbocatión bencílico ${gTxt}`;
  if (tipos.includes('vinilo')) return `carbocatión alílico ${gTxt}`;
  if (tipos.includes('ciclopropilo')) return 'catión ciclopropilcarbinilo';
  if (tipos.includes('alquinilo')) return `carbocatión propargílico ${gTxt}`;
  if (tipos.includes('halogeno')) return `α-halocarbocatión ${gTxt}`;
  if (grado === 0) return 'catión metilo';
  return `carbocatión alquílico ${gTxt}`;
}

function analizarOnio(mol, pos) {
  const nombres = { N: 'amonio', P: 'fosfonio', O: 'oxonio', S: 'sulfonio' };
  const nombre = nombres[pos.el] || 'halonio';
  const indice = { N: 72, P: 72, S: 68, O: 55 }[pos.el] ?? 35;
  const cat = categoriaDe(indice);
  return {
    tipoIon: 'EE',
    modo: 'onio',
    clase: `ion ${nombre} saturado (carga en heteroátomo con octeto completo)`,
    sitio: { atomo: pos.idx + 1, elemento: pos.el, h: pos.h },
    contribuciones: [],
    hia: null,
    indice,
    categoria: cat.nombre,
    interpretacion: cat.ee,
    estructurasResonantes: 1,
    sitiosDeslocalizacion: [pos.idx + 1],
    reordenamientos: [],
    notas: [
      'La carga reside en un heteroátomo con octeto completo (sin orbital vacío): no se aplica el modelo de afinidad por hidruro y el índice es cualitativo.',
      'Los iones onio saturados (R–NH3⁺, R–OH2⁺…) son más típicos de la ionización química (CI); en EI suelen perder moléculas neutras (NH3, H2O, ROH).',
    ],
    smilesAnalizado: sm.escribirSmiles(mol),
  };
}

function buscarReordenamientos(trabajo, sitio, hia) {
  const res = [];
  const probar = (m, descripcion) => {
    try {
      sm.recalcular(m);
      const r = analizarPar(m, { sinReordenamientos: true });
      if (r.hia !== null && r.hia < hia - 4) {
        res.push({ descripcion, smiles: sm.escribirSmiles(m), clase: r.clase, hia: r.hia, delta: r1(r.hia - hia), indice: r.indice });
      }
    } catch (_) { /* estructura no válida: se descarta */ }
  };
  for (const { vecino: v } of trabajo.adj[sitio]) {
    const av = trabajo.atoms[v];
    if (av.el !== 'C' || av.tieneMultiple || av.sexteto || av.carga !== 0) continue;
    if (av.h >= 1) {
      const m = sm.clonar(trabajo);
      m.atoms[sitio].h += 1;
      m.atoms[sitio].carga -= 1;
      m.atoms[v].h -= 1;
      m.atoms[v].carga += 1;
      probar(m, `Migración 1,2 de hidruro (H⁻) del átomo ${v + 1} al átomo ${sitio + 1}`);
    }
    for (const { vecino: w, enlace: e } of trabajo.adj[v]) {
      if (w === sitio) continue;
      const aw = trabajo.atoms[w];
      if (aw.el !== 'C' || aw.tieneMultiple || aw.sexteto) continue;
      if (trabajo.adj[sitio].some((x) => x.vecino === w)) continue;
      if (sm.enlaceEnAnillo(trabajo, v, w)) continue;
      const m = sm.clonar(trabajo);
      const b = m.bonds[e.id];
      if (b.a === v) b.a = sitio; else b.b = sitio;
      m.atoms[sitio].carga -= 1;
      m.atoms[v].carga += 1;
      probar(m, `Migración 1,2 de ${aw.h === 3 ? 'metilo' : 'alquilo'} (átomo ${w + 1}) del átomo ${v + 1} al átomo ${sitio + 1}`);
    }
  }
  const vistos = new Set();
  return res
    .sort((x, y) => x.hia - y.hia)
    .filter((r) => (vistos.has(r.smiles) ? false : vistos.add(r.smiles)))
    .slice(0, 3);
}

/** Análisis de un catión de electrones pares (EE⁺). */
function analizarPar(mol, opciones = {}) {
  const notas = [];
  const positivos = mol.atoms.filter((a) => a.carga > 0);
  if (!positivos.length) throw new Error('La estructura no tiene ningún átomo con carga positiva.');
  if (positivos.length > 1 || mol.atoms.some((a) => a.carga < 0)) {
    notas.push('La estructura tiene varias cargas formales; se analiza el centro positivo principal.');
  }
  let pos = positivos.find((a) => a.sexteto) || positivos[0];
  let trabajo = mol;
  let formaHet = null;

  if (!pos.sexteto) {
    // Ion onio con enlace π a carbono (C=N⁺, C=O⁺, C≡O⁺…): se analiza su forma carbenio.
    const pi = mol.adj[pos.idx].find((x) => x.enlace.orden >= 2 && mol.atoms[x.vecino].el === 'C');
    if (!pi || !['N', 'O', 'S', 'P'].includes(pos.el)) return analizarOnio(mol, pos);
    trabajo = sm.clonar(mol);
    trabajo.bonds[pi.enlace.id].orden -= 1;
    trabajo.atoms[pos.idx].carga -= 1;
    trabajo.atoms[pi.vecino].carga += 1;
    sm.recalcular(trabajo);
    formaHet = { el: pos.el, ordenOriginal: pi.enlace.orden, heteroatomo: pos.idx + 1, carbono: pi.vecino + 1, smiles: sm.escribirSmiles(trabajo) };
    pos = trabajo.atoms[pi.vecino];
    if (!pos.sexteto) return analizarOnio(mol, mol.atoms[formaHet.heteroatomo - 1]);
  }

  const A = trabajo.atoms;
  const sitio = pos.idx;
  const vecinos = trabajo.adj[sitio];
  const multiples = vecinos.filter((x) => x.enlace.orden >= 2);
  let modo;
  if (pos.el !== 'C') modo = 'heteroatomo';
  else if (multiples.some((x) => A[x.vecino].el === 'O')) modo = 'acilio';
  else if (multiples.some((x) => A[x.vecino].el === 'N')) modo = 'nitrilio';
  else if (multiples.length) modo = sm.enAnilloAromatico(trabajo, sitio) ? 'arilo' : 'vinilico';
  else modo = 'carbenio';

  const base = BASES[modo];
  const contribuciones = [{ concepto: `Base: ${base.texto}`, valor: base.valor, tipo: 'base' }];
  const estabilizantes = [];
  const desestabilizantes = [];
  const tipos = [];
  const sustituyentes = (modo === 'carbenio' || modo === 'heteroatomo') ? vecinos : vecinos.filter((x) => x.enlace.orden === 1);
  const escalaD = (modo === 'carbenio' || modo === 'heteroatomo') ? 1 : 0.5;

  for (const x of sustituyentes) {
    const c = clasificarSustituyente(trabajo, sitio, x.vecino);
    tipos.push(c.tipo);
    const etiqueta = `${c.etiqueta} [átomo ${x.vecino + 1}]`;
    if (c.s) {
      let s = c.s * FACTOR_MODO[modo];
      if (modo === 'acilio' || modo === 'nitrilio') {
        if (c.tipo === 'alquilo') s = 29;
        if (c.tipo === 'arilo') s = 42;
      }
      if (s > 0) estabilizantes.push({ concepto: etiqueta, s });
    }
    if (c.d) desestabilizantes.push({ concepto: etiqueta, valor: c.d * escalaD });
    if (c.dExtra) {
      desestabilizantes.push({ concepto: `efecto inductivo β: ${c.nElectroneg} átomo(s) electronegativo(s) en el C α [átomo ${x.vecino + 1}]`, valor: c.dExtra * escalaD });
    }
  }
  if (modo === 'arilo') {
    notas.push('En el catión arilo el orbital vacío es sp² y está en el plano del anillo: no se conjuga con el sistema π, por eso el anillo no lo estabiliza.');
  }

  for (const it of combinarEstabilizaciones(estabilizantes)) {
    const aten = it.efectivo < it.s - 0.05 ? ` (atenuado por saturación; aislado −${r1(it.s)})` : '';
    contribuciones.push({ concepto: it.concepto + aten, valor: -it.efectivo, tipo: 'estabiliza' });
  }
  for (const d of desestabilizantes) contribuciones.push({ ...d, tipo: 'desestabiliza' });

  // Deslocalización por resonancia (solo carbenio con orbital p vacío).
  let sitiosRes = [sitio];
  let formasDonador = estabilizantes.filter((e) => /par libre|halóg|flúor|cloro|bromo|yodo/.test(e.concepto)).length;
  const ciclosSitio = trabajo.ciclos.map((c, k) => k).filter((k) => trabajo.ciclos[k].includes(sitio));
  let anilloAromatico = null;
  let anilloAnti = null;
  if (modo === 'carbenio') {
    sitiosRes = sm.sitiosResonancia(trabajo, sitio);
    const ar = ciclosSitio.find((k) => trabajo.huckel[k].aromatico);
    if (ar !== undefined) anilloAromatico = { tam: trabajo.ciclos[ar].length, e: trabajo.huckel[ar].electrones };
    else {
      const an = ciclosSitio.find((k) => trabajo.huckel[k].antiaromatico);
      if (an !== undefined) anilloAnti = { tam: trabajo.ciclos[an].length, e: trabajo.huckel[an].electrones };
    }

    const setSitios = new Set(sitiosRes);
    let remoto = 0;
    const detalle = [];
    for (const r of sitiosRes.slice(1)) {
      for (const { vecino: y, enlace } of trabajo.adj[r]) {
        if (enlace.orden !== 1 || setSitios.has(y)) continue;
        const ay = A[y];
        if (['N', 'O', 'S'].includes(ay.el) && ay.paresLibres > 0 && ay.carga === 0 && !ay.tieneMultiple && !sm.enAnilloAromatico(trabajo, y)) {
          remoto += ay.el === 'N' ? 22 : 15;
          detalle.push(`${ay.el} en átomo ${y + 1}`);
          formasDonador++;
        } else if (ay.el === 'C' && !ay.tieneMultiple && !ay.sexteto && ay.carga === 0) {
          remoto += 6;
          detalle.push(`alquilo en átomo ${y + 1}`);
        }
      }
    }
    if (remoto > 0) {
      contribuciones.push({ concepto: `Sustituyentes donadores en posiciones de resonancia (${detalle.join(', ')})`, valor: -Math.min(remoto, 35), tipo: 'estabiliza' });
    }

    const nVin = tipos.filter((t) => t === 'vinilo').length;
    const nAr = tipos.filter((t) => t === 'arilo').length;
    const extra = sitiosRes.length - 1 - nVin - 3 * nAr;
    if (!anilloAromatico && !anilloAnti && extra > 0) {
      contribuciones.push({ concepto: `Conjugación extendida (${extra} posición(es) adicional(es) para la carga)`, valor: -8 * Math.min(extra, 3), tipo: 'estabiliza' });
    }
    if (anilloAromatico) {
      contribuciones.push({ concepto: `Aromaticidad: anillo de ${anilloAromatico.tam} átomos con ${anilloAromatico.e} electrones π (Hückel 4n+2)`, valor: -30, tipo: 'estabiliza' });
    } else if (anilloAnti) {
      contribuciones.push({ concepto: `Antiaromaticidad: anillo de ${anilloAnti.tam} átomos con ${anilloAnti.e} electrones π (4n)`, valor: 45, tipo: 'desestabiliza' });
    }
  }

  const tamAnillo = ciclosSitio.length ? Math.min(...ciclosSitio.map((k) => trabajo.ciclos[k].length)) : null;
  if (tamAnillo === 3) contribuciones.push({ concepto: 'Tensión de anillo de 3 miembros en el centro catiónico', valor: 20, tipo: 'desestabiliza' });
  if (tamAnillo === 4) contribuciones.push({ concepto: 'Tensión de anillo de 4 miembros en el centro catiónico', valor: 10, tipo: 'desestabiliza' });

  // Polarizabilidad: átomos pesados en posición β.
  const alfa = new Set(vecinos.map((x) => x.vecino));
  const beta = new Set();
  for (const a of alfa) for (const { vecino: b } of trabajo.adj[a]) if (b !== sitio && !alfa.has(b)) beta.add(b);
  if (beta.size) {
    contribuciones.push({ concepto: `Polarizabilidad: ${beta.size} átomo(s) pesado(s) en posición β`, valor: -Math.min(9, 1.5 * beta.size), tipo: 'estabiliza' });
  }

  const hia = r1(contribuciones.reduce((s, c) => s + c.valor, 0));
  contribuciones.forEach((c) => { c.valor = r1(c.valor); });
  const indice = indiceDesdeHia(hia);
  const cat = categoriaDe(indice);
  const grado = vecinos.filter((x) => A[x.vecino].el === 'C').length;
  const clase = nombrarClase(modo, { tipos, anilloAromatico, anilloAnti, grado, formaHet });

  if (formaHet) {
    const tipoTxt = formaHet.ordenOriginal === 3 ? 'triple' : 'doble';
    notas.push(`La carga se dio sobre el ${formaHet.el} (átomo ${formaHet.heteroatomo}) con enlace ${tipoTxt} a C. Se analiza la forma resonante con la carga en el carbono (${formaHet.smiles}); el ${formaHet.el} dona su par libre y la forma con octeto completo es la contribuyente principal.`);
  }
  if (anilloAromatico) notas.push(`La carga está deslocalizada por igual en los ${sitiosRes.length} átomos del anillo aromático.`);
  if (modo === 'carbenio' && tipos.includes('arilo') && pos.h === 2) {
    notas.push('Los iones bencilo (C7H7⁺ y derivados) se isomerizan en gran parte al ion tropilio (anillo de 7 miembros, aromático), aún más estable (HIA ≈ 201 kcal/mol).');
  }
  if (modo === 'carbenio' && grado === 1 && tipos.every((t) => t === 'alquilo' || t === 'otro') && pos.h === 2 && trabajo.atoms.length > 2) {
    notas.push('Los carbocationes alquílicos primarios (salvo metilo y etilo) apenas existen como tales en fase gaseosa: se reordenan a iones secundarios/terciarios o a ciclopropano protonado.');
  }

  const reordenamientos = !opciones.sinReordenamientos && modo === 'carbenio' && !anilloAromatico
    ? buscarReordenamientos(trabajo, sitio, hia)
    : [];

  return {
    tipoIon: 'EE',
    modo,
    clase,
    sitio: { atomo: sitio + 1, elemento: pos.el, h: pos.h },
    formaHeteroatomo: formaHet,
    estructurasResonantes: sitiosRes.length + formasDonador,
    sitiosDeslocalizacion: sitiosRes.map((i) => i + 1),
    contribuciones,
    hia,
    indice,
    categoria: cat.nombre,
    interpretacion: cat.ee,
    reordenamientos,
    notas,
    smilesAnalizado: sm.escribirSmiles(trabajo),
  };
}

// ---------------------------------------------------------------- iones OE⁺•

const GRUPOS_OE = {
  aromatico: { nombre: 'anillo aromático', puntaje: 95, ie: 9.0, sitio: 'sistema π aromático',
    frag: 'M⁺• intenso. Con cadena lateral: ruptura bencílica → ion C7H7⁺ (m/z 91, tropilio).' },
  fenol: { nombre: 'fenol', puntaje: 85, ie: 8.5, sitio: 'sistema π aromático conjugado con el O',
    frag: 'M⁺• intenso; pérdida de CO (M−28) y de CHO• (M−29).' },
  conjugado: { nombre: 'sistema π conjugado', puntaje: 85, ie: 8.8, sitio: 'sistema π conjugado',
    frag: 'M⁺• relativamente intenso; rupturas alílicas.' },
  aliciclico: { nombre: 'anillo alicíclico', puntaje: 75, ie: 9.8, sitio: 'enlaces σ C–C del anillo',
    frag: 'M⁺• apreciable (romper un enlace no separa la molécula); pérdida de alquenos neutros (p. ej. C2H4).' },
  sulfuro: { nombre: 'sulfuro / tiol', puntaje: 72, ie: 8.7, sitio: 'par libre del azufre',
    frag: 'Ruptura α al azufre → ion tiocarbenio (CH2=SR⁺).' },
  alqueno: { nombre: 'alqueno', puntaje: 70, ie: 9.7, sitio: 'enlace π C=C',
    frag: 'Ruptura alílica → catión alilo estabilizado.' },
  alquino: { nombre: 'alquino', puntaje: 65, ie: 10.2, sitio: 'enlace π C≡C',
    frag: 'Ruptura propargílica → catión propargilo (m/z 39).' },
  lineal: { nombre: 'alcano lineal', puntaje: 62, ie: 10.4, sitio: 'enlace σ C–C',
    frag: 'Serie CnH2n+1⁺ (m/z 29, 43, 57, 71…) con máximo en C3–C4; picos separados por 14 u.' },
  amida: { nombre: 'amida', puntaje: 60, ie: 9.3, sitio: 'par libre del O/N de la amida',
    frag: 'Ruptura α y reordenamiento de McLafferty; m/z 44 (H2N–C≡O⁺) en amidas primarias.' },
  cetona: { nombre: 'cetona / aldehído', puntaje: 60, ie: 9.6, sitio: 'par libre (n) del O carbonílico',
    frag: 'Ruptura α → ion acilio R–C≡O⁺; reordenamiento de McLafferty si hay H en γ.' },
  yoduro: { nombre: 'yoduro', puntaje: 60, ie: 9.3, sitio: 'par libre del yodo',
    frag: 'Pérdida de I• (M−127); ion I⁺ a m/z 127.' },
  amina: { nombre: 'amina', puntaje: 55, ie: 8.5, sitio: 'par libre (n) del nitrógeno',
    frag: 'Ruptura α al N → ion iminio (CH2=NH2⁺, m/z 30, en aminas primarias con CH2 en α).' },
  bromuro: { nombre: 'bromuro', puntaje: 50, ie: 10.3, sitio: 'par libre del bromo',
    frag: 'Doblete M/M+2 ≈ 1:1; pérdida de Br•.' },
  ester: { nombre: 'éster', puntaje: 48, ie: 10.2, sitio: 'par libre del O carbonílico',
    frag: 'Ruptura α (pérdida de •OR) → ion acilio; reordenamiento de McLafferty.' },
  eter: { nombre: 'éter', puntaje: 45, ie: 9.8, sitio: 'par libre del oxígeno',
    frag: 'Ruptura α → ion oxonio (CH2=OR⁺); ruptura del enlace C–O.' },
  acido: { nombre: 'ácido carboxílico', puntaje: 40, ie: 10.5, sitio: 'par libre del O carbonílico',
    frag: 'Pérdida de •OH (M−17) y •COOH (M−45); McLafferty → m/z 60.' },
  cloruro: { nombre: 'cloruro', puntaje: 40, ie: 11.0, sitio: 'par libre del cloro',
    frag: 'Doblete M/M+2 ≈ 3:1; pérdida de Cl• o de HCl.' },
  nitrilo: { nombre: 'nitrilo', puntaje: 30, ie: 11.5, sitio: 'enlace π C≡N',
    frag: 'M⁺• débil; M−1 frecuente; McLafferty → m/z 41.' },
  fluoruro: { nombre: 'fluoruro', puntaje: 30, ie: 12.0, sitio: 'par libre del flúor',
    frag: 'Pérdida de HF (M−20).' },
  ramificado: { nombre: 'alcano ramificado', puntaje: 30, ie: 10.0, sitio: 'enlace σ C–C en el punto de ramificación',
    frag: 'Ruptura preferente en los carbonos ramificados → carbocationes secundarios/terciarios.' },
  alcohol: { nombre: 'alcohol', puntaje: 15, ie: 10.4, sitio: 'par libre del oxígeno',
    frag: 'M⁺• débil o ausente; pérdida de H2O (M−18); ruptura α → CH2=OH⁺ (m/z 31) en alcoholes primarios.' },
};

function detectarGrupos(mol) {
  const A = mol.atoms;
  const g = {};
  const add = (k, n = 1) => { g[k] = (g[k] || 0) + n; };
  const aromaticos = mol.ciclos.filter((c, k) => mol.huckel[k].aromatico && c.length >= 5 && c.length <= 7);
  if (aromaticos.length) add('aromatico', aromaticos.length);
  const enAromatico = new Set(aromaticos.flat());
  if (mol.ciclos.some((c, k) => !mol.huckel[k].aromatico && !c.some((x) => enAromatico.has(x)))) add('aliciclico');

  const conjugadoCon = (i, excluir) => mol.adj[i].some((x) => x.vecino !== excluir && x.enlace.orden === 1
    && (A[x.vecino].tieneMultiple || enAromatico.has(x.vecino)));
  for (const b of mol.bonds) {
    const x = A[b.a];
    const y = A[b.b];
    if (b.orden === 2 && x.el === 'C' && y.el === 'C' && !(enAromatico.has(b.a) && enAromatico.has(b.b))) {
      add(conjugadoCon(b.a, b.b) || conjugadoCon(b.b, b.a) ? 'conjugado' : 'alqueno');
    }
    if (b.orden === 3 && x.el === 'C' && y.el === 'C') add('alquino');
    if (b.orden === 3 && [x.el, y.el].sort().join('') === 'CN') add('nitrilo');
  }

  for (const a of A) {
    const adj = mol.adj[a.idx];
    if (a.el === 'C' && esCarbonoAcilo(mol, a.idx)) {
      const simples = adj.filter((x) => x.enlace.orden === 1);
      const oSimple = simples.find((x) => A[x.vecino].el === 'O');
      if (oSimple) add(A[oSimple.vecino].h >= 1 ? 'acido' : 'ester');
      else if (simples.some((x) => A[x.vecino].el === 'N')) add('amida');
      else add('cetona');
    }
    if (a.el === 'O' && !a.tieneMultiple && a.carga === 0 && !enAromatico.has(a.idx)) {
      if (adj.some((x) => esCarbonoAcilo(mol, x.vecino))) continue;
      if (a.h >= 1) add(adj.some((x) => enAromatico.has(x.vecino)) ? 'fenol' : 'alcohol');
      else if (adj.length === 2) add('eter');
    }
    if (a.el === 'N' && !a.tieneMultiple && a.paresLibres > 0 && a.carga === 0 && !enAromatico.has(a.idx)
      && !adj.some((x) => esCarbonoAcilo(mol, x.vecino))) add('amina');
    if (a.el === 'S' && !a.tieneMultiple && a.sumaEnlaces === 2 && !enAromatico.has(a.idx)) add('sulfuro');
    if (a.el === 'I') add('yoduro');
    if (a.el === 'Br') add('bromuro');
    if (a.el === 'Cl') add('cloruro');
    if (a.el === 'F') add('fluoruro');
  }

  let terciarios = 0;
  let cuaternarios = 0;
  for (const a of A) {
    if (a.el !== 'C' || a.tieneMultiple) continue;
    const nC = mol.adj[a.idx].filter((x) => A[x.vecino].el === 'C').length;
    if (nC === 3) terciarios++;
    if (nC === 4) cuaternarios++;
  }
  const soloCH = A.every((a) => a.el === 'C');
  if (soloCH && !Object.keys(g).length) add(terciarios + cuaternarios > 0 ? 'ramificado' : 'lineal');
  return { grupos: g, terciarios, cuaternarios, enAromatico };
}

function tieneMcLafferty(mol) {
  const A = mol.atoms;
  for (const c of A) {
    if (!esCarbonoAcilo(mol, c.idx)) continue;
    for (const { vecino: al, enlace: e1 } of mol.adj[c.idx]) {
      if (e1.orden !== 1 || A[al].el !== 'C') continue;
      for (const { vecino: be, enlace: e2 } of mol.adj[al]) {
        if (be === c.idx || e2.orden !== 1 || A[be].el !== 'C') continue;
        for (const { vecino: ga, enlace: e3 } of mol.adj[be]) {
          if (ga === al || e3.orden !== 1 || !['C', 'O', 'N'].includes(A[ga].el)) continue;
          if (A[ga].h >= 1) return true;
        }
      }
    }
  }
  return false;
}

/** Análisis cualitativo de un ion radical (OE⁺•). */
function analizarRadical(mol, { ionMolecular }) {
  const { grupos, terciarios, cuaternarios, enAromatico } = detectarGrupos(mol);
  const claves = Object.keys(grupos);
  const notas = [];
  if (ionMolecular) notas.push('La estructura dada es neutra: se analiza como ion molecular M⁺• (la molécula ha perdido un electrón).');
  if (!claves.length) {
    claves.push('lineal');
    grupos.lineal = 1;
  }
  const lista = claves.map((k) => ({ clave: k, ...GRUPOS_OE[k], n: grupos[k] }));
  const principal = [...lista].sort((x, y) => y.puntaje - x.puntaje)[0];
  const contribuciones = [{ concepto: `Grupo más estabilizante: ${principal.nombre}`, valor: principal.puntaje, tipo: 'base' }];
  const hayAromatico = claves.includes('aromatico') || claves.includes('fenol');

  if (grupos.aromatico > 1) contribuciones.push({ concepto: `${grupos.aromatico} anillos aromáticos (mayor deslocalización)`, valor: Math.min(10, 3 * (grupos.aromatico - 1)), tipo: 'estabiliza' });
  if (!hayAromatico) {
    const pen = Math.min(25, 4 * terciarios + 8 * cuaternarios);
    if (pen) contribuciones.push({ concepto: `Ramificación (${terciarios} C terciario(s), ${cuaternarios} C cuaternario(s)): favorece la fragmentación`, valor: -pen, tipo: 'desestabiliza' });
    const pesados = mol.atoms.length;
    if (pesados > 12) contribuciones.push({ concepto: `Tamaño (${pesados} átomos pesados): más canales de fragmentación`, valor: -Math.min(15, pesados - 12), tipo: 'desestabiliza' });
  }
  let indice = contribuciones.reduce((s, c) => s + c.valor, 0);
  if (claves.includes('alcohol') && !hayAromatico && indice > 30) {
    contribuciones.push({ concepto: 'Grupo alcohol: pierde H2O con mucha facilidad', valor: 30 - indice, tipo: 'desestabiliza' });
    indice = 30;
  }
  indice = Math.max(0, Math.min(100, Math.round(indice)));

  // Sitio probable de ionización: el grupo con menor energía de ionización.
  const masFacil = [...lista].sort((x, y) => x.ie - y.ie)[0];
  let ie = masFacil.ie;
  let sitio = masFacil.sitio;
  const aminaArilo = mol.atoms.some((a) => a.el === 'N' && !a.tieneMultiple && a.paresLibres > 0
    && mol.adj[a.idx].some((x) => enAromatico.has(x.vecino)));
  if (aminaArilo) {
    ie = 7.7;
    sitio = 'par libre del N conjugado con el anillo aromático (tipo anilina)';
  }
  ie = r1(Math.max(7, ie - Math.min(1, 0.05 * Math.max(0, mol.atoms.length - 6))));

  const fragmentaciones = lista.map((g) => `${g.nombre}: ${g.frag}`);
  if (tieneMcLafferty(mol)) {
    fragmentaciones.push('Reordenamiento de McLafferty posible: hay un H en posición γ respecto al C=O (pérdida de un alqueno neutro y formación de un ion radical enol).');
  }
  const cat = categoriaDe(indice);
  return {
    tipoIon: 'OE',
    modo: 'radical',
    clase: ionMolecular ? 'ion molecular radical (M⁺•)' : 'ion radical (OE⁺•)',
    grupos: lista.map((g) => ({ grupo: g.nombre, cantidad: g.n })),
    sitioIonizacion: sitio,
    energiaIonizacion: ie,
    contribuciones,
    hia: null,
    indice,
    categoria: cat.nombre,
    interpretacion: cat.oe,
    fragmentaciones,
    reordenamientos: [],
    notas,
    smilesAnalizado: sm.escribirSmiles(mol),
  };
}

// ---------------------------------------------------------------- entrada principal

/**
 * Analiza un fragmento: { formula: 'C7H7+', smiles: '[cH+]1cccccc1', nombre?: '...' }.
 */
function analizarFragmento({ formula, smiles, nombre } = {}) {
  if (!smiles || !String(smiles).trim()) throw new Error('Debe indicar la estructura del fragmento en SMILES.');
  const mol = sm.prepararMolecula(smiles);
  const carga = sm.cargaTotal(mol);
  const conteo = sm.conteoAtomos(mol);
  const hillEstructura = formulaHill(conteo);
  const advertencias = [];

  if (carga < 0) throw new Error('La estructura tiene carga neta negativa; este programa analiza cationes (EI-MS en modo positivo).');

  let fi;
  if (formula && String(formula).trim()) {
    fi = parsearFormula(formula);
  } else {
    fi = parsearFormula(`${hillEstructura}${'+'.repeat(Math.max(carga, 1))}`);
    advertencias.push('No se indicó la fórmula: se calculó a partir de la estructura.');
  }
  const coincide = Object.keys({ ...fi.atomos, ...conteo }).every((el) => (fi.atomos[el] || 0) === (conteo[el] || 0));
  if (!coincide) {
    advertencias.push(`La fórmula indicada (${fi.hill}) no coincide con la de la estructura (${hillEstructura}). El análisis se basa en la estructura.`);
  }

  let resultado;
  if (carga === 0) {
    if (fi.tipo === 'EE' && coincide) {
      throw new Error('La estructura es neutra. Para un catión EE⁺ indique el centro catiónico en el SMILES, p. ej. C[CH2+] o CC#[O+].');
    }
    if (fi.tipo === 'EE') {
      throw new Error(`La fórmula ${fi.texto} es de un catión EE⁺ pero la estructura es neutra (${hillEstructura}). Indique el centro catiónico en el SMILES, p. ej. [CH2+]c1ccccc1.`);
    }
    resultado = analizarRadical(mol, { ionMolecular: true });
  } else {
    if (carga > 1) advertencias.push(`La estructura tiene carga +${carga}; se analiza el centro catiónico principal.`);
    if (fi.carga !== carga) advertencias.push(`La carga de la fórmula (${fi.carga > 0 ? '+' : ''}${fi.carga}) no coincide con la de la estructura (+${carga}).`);
    const radicales = mol.atoms.reduce((s, a) => s + a.radical, 0);
    if (radicales % 2 === 1) {
      if (fi.tipo === 'EE' && coincide) advertencias.push('La estructura tiene un electrón desapareado: se analiza como ion radical.');
      resultado = analizarRadical(mol, { ionMolecular: false });
    } else {
      if (fi.tipo === 'OE' && coincide) advertencias.push('La fórmula indica un ion radical, pero la estructura no tiene electrones desapareados.');
      resultado = analizarPar(mol);
    }
  }

  const tipoAnalisis = resultado.tipoIon;
  return {
    nombre: nombre || null,
    entrada: { formula: formula || null, smiles: String(smiles).trim() },
    formula: fi,
    estructura: { formula: hillEstructura, carga, coincideConFormula: coincide, atomosPesados: mol.atoms.length },
    ionConocido: buscarIonConocido(hillEstructura, tipoAnalisis),
    advertencias: [...advertencias, ...fi.notas],
    ...resultado,
  };
}

/** Analiza varios fragmentos y los ordena de más a menos estable. */
function compararFragmentos(lista) {
  const resultados = lista.map((f) => {
    try {
      return analizarFragmento(f);
    } catch (e) {
      return { nombre: f.nombre || null, entrada: f, error: e.message, indice: -1 };
    }
  });
  return resultados.sort((a, b) => b.indice - a.indice);
}

module.exports = {
  analizarFragmento,
  compararFragmentos,
  indiceDesdeHia,
  CATEGORIAS,
};
