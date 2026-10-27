'use strict';

const ANCHO = 74;
const linea = (c = '-') => c.repeat(ANCHO);

function barra(indice, largo = 30) {
  const llenos = Math.round((indice / 100) * largo);
  return '[' + '#'.repeat(llenos) + '.'.repeat(largo - llenos) + ']';
}

function fila(etiqueta, valor) {
  return `  ${etiqueta.padEnd(22)}: ${valor}`;
}

function signo(v) {
  return (v >= 0 ? '+' : '−') + Math.abs(v).toFixed(1);
}

function ajustar(texto, sangria) {
  const palabras = String(texto).split(/\s+/);
  const lineas = [];
  let actual = '';
  for (const p of palabras) {
    if ((actual + ' ' + p).trim().length > ANCHO - sangria) {
      lineas.push(actual.trim());
      actual = p;
    } else {
      actual += ' ' + p;
    }
  }
  if (actual.trim()) lineas.push(actual.trim());
  return lineas.map((l, i) => (i === 0 ? l : ' '.repeat(sangria) + l)).join('\n');
}

function ajustarAncho(texto, ancho) {
  const lineas = [];
  let actual = '';
  for (const p of texto.split(/\s+/)) {
    if (actual && (actual + ' ' + p).length > ancho) {
      lineas.push(actual);
      actual = '  ' + p;
    } else {
      actual = actual ? actual + ' ' + p : p;
    }
  }
  if (actual) lineas.push(actual);
  return lineas;
}

function vinetas(titulo, lista) {
  if (!lista || !lista.length) return [];
  return ['', `  ${titulo}`, ...lista.map((t) => `    - ${ajustar(t, 6)}`)];
}

/** Informe de texto de un resultado de analizarFragmento. */
function formatearInforme(r) {
  if (r.error) {
    return [linea('='), `  ${r.nombre || r.entrada.smiles}: ERROR – ${r.error}`, linea('=')].join('\n');
  }
  const f = r.formula;
  const tipoTxt = r.tipoIon === 'EE' ? 'ion de electrones pares, EE⁺' : 'ion radical, OE⁺•';
  const out = [
    linea('='),
    '  ESTABILIDAD DEL CATIÓN FRAGMENTO (EI-MS)',
    linea('='),
  ];
  if (r.nombre) out.push(fila('Fragmento', r.nombre));
  out.push(
    fila('Fórmula', `${f.texto}  →  ${f.hill} (${tipoTxt})`),
    fila('Estructura (SMILES)', `${r.entrada.smiles}  (fórmula: ${r.estructura.formula} ${r.estructura.coincideConFormula ? '✓ coincide' : '✗ NO coincide'})`),
    fila('m/z (monoisotópica)', `${f.mz.toFixed(4)}  (nominal ${f.mzNominal})`),
    fila('RDB', r.tipoIon === 'EE'
      ? `${f.rdb} → ${f.rdbEstructural} anillos + enlaces π en la estructura del catión`
      : `${f.rdb} (anillos + enlaces π)`),
    fila('Isótopos (aprox.)', `M+1 ≈ ${f.isotopos.m1.toFixed(1)} %   M+2 ≈ ${f.isotopos.m2.toFixed(1)} %`),
    `  ${'Regla del nitrógeno'.padEnd(22)}: ${ajustar(f.reglaNitrogeno, 26)}`,
  );
  if (r.ionConocido) out.push(`  ${'Ion diagnóstico'.padEnd(22)}: ${ajustar(r.ionConocido, 26)}`);

  out.push('', linea(), fila('Clasificación', r.clase));
  if (r.tipoIon === 'EE') {
    out.push(fila('Centro catiónico', `átomo ${r.sitio.atomo} (${r.sitio.elemento}${r.sitio.h ? 'H' + (r.sitio.h > 1 ? r.sitio.h : '') : ''})`));
    if (r.smilesAnalizado !== r.entrada.smiles) out.push(fila('Forma analizada', r.smilesAnalizado));
    out.push(fila('Formas resonantes', `${r.estructurasResonantes} (carga de carbono sobre el/los átomo(s) ${r.sitiosDeslocalizacion.join(', ')})`));
  } else {
    out.push(
      fila('Grupos detectados', r.grupos.map((g) => (g.cantidad > 1 ? `${g.grupo} ×${g.cantidad}` : g.grupo)).join(', ')),
      `  ${'Sitio de ionización'.padEnd(22)}: ${ajustar(r.sitioIonizacion, 26)}`,
      fila('EI estimada', `≈ ${r.energiaIonizacion} eV (energía para arrancar el electrón)`),
    );
  }

  out.push('', r.tipoIon === 'EE' && r.hia !== null
    ? '  Contribuciones a la afinidad por hidruro estimada, HIA (kcal/mol):'
    : '  Contribuciones al índice de estabilidad:');
  for (const c of r.contribuciones) {
    const valor = c.tipo === 'base' ? c.valor.toFixed(1).padStart(7) : signo(c.valor).padStart(7);
    const partes = ajustarAncho(c.concepto, ANCHO - 16);
    partes.forEach((p, k) => {
      out.push(k === partes.length - 1 ? `    ${p.padEnd(ANCHO - 14)}${valor}` : `    ${p}`);
    });
  }
  if (r.hia !== null) {
    out.push(`    ${'-'.repeat(ANCHO - 6)}`, `    ${'HIA estimada (menor = más estable)'.padEnd(ANCHO - 14)}${r.hia.toFixed(1).padStart(7)}`);
  }

  out.push(
    '',
    `  ÍNDICE DE ESTABILIDAD: ${barra(r.indice)} ${r.indice}/100  →  ${r.categoria}`,
    `  ${ajustar(r.interpretacion, 2)}`,
  );

  if (r.reordenamientos && r.reordenamientos.length) {
    out.push(...vinetas('Reordenamientos favorables (el catión puede isomerizarse):',
      r.reordenamientos.map((x) => `${x.descripcion} → ${x.smiles} (${x.clase}); HIA ${x.hia} (${x.delta} kcal/mol), índice ${x.indice}.`)));
  }
  out.push(...vinetas('Fragmentaciones esperadas:', r.fragmentaciones));
  out.push(...vinetas('Notas:', r.notas));
  out.push(...vinetas('Advertencias:', r.advertencias));
  out.push(linea('='));
  return out.join('\n');
}

/** Tabla comparativa para varios fragmentos (ya ordenados). */
function formatearComparacion(resultados) {
  const out = [
    linea('='),
    '  COMPARACIÓN DE ESTABILIDAD (de más a menos estable)',
    linea('='),
    `  ${'#'.padEnd(3)}${'Fragmento'.padEnd(26)}${'m/z'.padStart(9)}  ${'Tipo'.padEnd(5)}${'HIA'.padStart(7)}  ${'Índice'.padStart(6)}  Categoría`,
    `  ${linea().slice(2)}`,
  ];
  resultados.forEach((r, i) => {
    const nombre = (r.nombre || r.entrada.smiles).slice(0, 25);
    if (r.error) {
      out.push(`  ${String(i + 1).padEnd(3)}${nombre.padEnd(26)}  ERROR: ${r.error}`);
      return;
    }
    out.push(`  ${String(i + 1).padEnd(3)}${nombre.padEnd(26)}${r.formula.mz.toFixed(2).padStart(9)}  ${r.tipoIon.padEnd(5)}${(r.hia === null ? '—' : r.hia.toFixed(1)).padStart(7)}  ${String(r.indice).padStart(6)}  ${r.categoria}`);
  });
  out.push(
    linea('='),
    '  HIA = afinidad por hidruro estimada (kcal/mol), solo para iones EE⁺: menor HIA = más estable.',
    '  En iones OE⁺• el índice es cualitativo (abundancia esperada del ion radical).',
    linea('='),
  );
  return out.join('\n');
}

module.exports = { formatearInforme, formatearComparacion };
