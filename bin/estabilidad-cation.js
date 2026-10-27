#!/usr/bin/env node
'use strict';

const fs = require('fs');
const path = require('path');
const readline = require('readline/promises');
const {
  analizarFragmento, compararFragmentos, formatearInforme, formatearComparacion, ejemplos,
} = require('../src');

const AYUDA = `
Estabilidad de cationes fragmento de EI-MS

Uso:
  estabilidad-cation -f <fórmula> -s <SMILES> [-n <nombre>] [--json]
  estabilidad-cation <fórmula> <SMILES>
  estabilidad-cation --comparar <archivo.json> [--json]
  estabilidad-cation --ejemplos [--detalle] [--json]
  estabilidad-cation                      (modo interactivo)

Opciones:
  -f, --formula   Fórmula del ion: C7H7+, C2H3O+, C6H6+• (o C6H6.+ para radical)
  -s, --smiles    Estructura en SMILES con el centro catiónico marcado:
                    [CH2+]c1ccccc1   bencilo        C[C+](C)C   tert-butilo
                    CC#[O+]          acetilo        C=[NH2+]    iminio
                    [cH+]1cccccc1    tropilio       c1ccccc1    benceno M+• (neutro)
  -n, --nombre    Nombre del fragmento (opcional)
  --comparar      Archivo JSON con [{ "nombre", "formula", "smiles" }, ...]
  --ejemplos      Analiza y ordena fragmentos típicos de EI-MS
  --detalle       Con --ejemplos/--comparar: muestra también cada informe
  --json          Salida en JSON
  -h, --ayuda     Muestra esta ayuda
`;

function leerArgs(argv) {
  const op = { posicionales: [] };
  for (let i = 0; i < argv.length; i++) {
    const a = argv[i];
    const sig = () => {
      if (i + 1 >= argv.length) throw new Error(`Falta el valor de ${a}.`);
      return argv[++i];
    };
    switch (a) {
      case '-f': case '--formula': op.formula = sig(); break;
      case '-s': case '--smiles': op.smiles = sig(); break;
      case '-n': case '--nombre': op.nombre = sig(); break;
      case '--comparar': op.comparar = sig(); break;
      case '--ejemplos': op.ejemplos = true; break;
      case '--detalle': op.detalle = true; break;
      case '--json': op.json = true; break;
      case '-h': case '--ayuda': case '--help': op.ayuda = true; break;
      default:
        if (a.startsWith('--')) throw new Error(`Opción desconocida: ${a}`);
        op.posicionales.push(a);
    }
  }
  if (!op.formula && !op.smiles && op.posicionales.length) {
    if (op.posicionales.length === 1) op.smiles = op.posicionales[0];
    else [op.formula, op.smiles] = op.posicionales;
  }
  return op;
}

function mostrarLista(lista, op) {
  const res = compararFragmentos(lista);
  if (op.json) {
    console.log(JSON.stringify(res, null, 2));
    return;
  }
  if (op.detalle) for (const r of res) console.log(formatearInforme(r) + '\n');
  console.log(formatearComparacion(res));
}

async function interactivo() {
  const rl = readline.createInterface({ input: process.stdin, output: process.stdout });
  console.log('Estabilidad de cationes fragmento (EI-MS). Deje la fórmula vacía para salir.');
  console.log('Ejemplo: fórmula C7H7+ y SMILES [cH+]1cccccc1\n');
  try {
    for (;;) {
      const formula = (await rl.question('Fórmula del fragmento: ')).trim();
      if (!formula) break;
      const smiles = (await rl.question('Estructura (SMILES)   : ')).trim();
      try {
        console.log('\n' + formatearInforme(analizarFragmento({ formula, smiles })) + '\n');
      } catch (e) {
        console.error(`Error: ${e.message}\n`);
      }
    }
  } finally {
    rl.close();
  }
}

async function main() {
  const op = leerArgs(process.argv.slice(2));
  if (op.ayuda) {
    console.log(AYUDA);
    return;
  }
  if (op.ejemplos) {
    mostrarLista(ejemplos, op);
    return;
  }
  if (op.comparar) {
    const ruta = path.resolve(op.comparar);
    const lista = JSON.parse(fs.readFileSync(ruta, 'utf8'));
    if (!Array.isArray(lista)) throw new Error('El archivo debe contener un arreglo JSON de fragmentos.');
    mostrarLista(lista, op);
    return;
  }
  if (op.smiles) {
    const r = analizarFragmento({ formula: op.formula, smiles: op.smiles, nombre: op.nombre });
    console.log(op.json ? JSON.stringify(r, null, 2) : formatearInforme(r));
    return;
  }
  if (op.formula) throw new Error('Falta la estructura: use -s "<SMILES>".');
  await interactivo();
}

main().catch((e) => {
  console.error(`Error: ${e.message}`);
  process.exitCode = 1;
});
