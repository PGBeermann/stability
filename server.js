'use strict';

/*
 * Servidor web del analizador de estabilidad de cationes (sin dependencias).
 *   GET  /               página web
 *   POST /api/analizar   { formula, smiles, nombre? } → resultado
 *   POST /api/comparar   [{ formula, smiles, nombre? }, ...] → resultados ordenados
 *   GET  /api/ejemplos   fragmentos de ejemplo ya analizados
 *   GET  /salud          comprobación de estado
 */

const http = require('http');
const fs = require('fs');
const path = require('path');
const { analizarFragmento, compararFragmentos, ejemplos } = require('./src');

const PUERTO = parseInt(process.env.PORT || '3000', 10);
const MAX_CUERPO = 20 * 1024;
const MAX_SMILES = 200;
const MAX_FRAGMENTOS = 30;

const PAGINA = fs.readFileSync(path.join(__dirname, 'public', 'index.html'));

const CABECERAS_SEGURIDAD = {
  'X-Content-Type-Options': 'nosniff',
  'X-Frame-Options': 'DENY',
  'Referrer-Policy': 'no-referrer',
  'Content-Security-Policy': "default-src 'self'; style-src 'self' 'unsafe-inline'; script-src 'self' 'unsafe-inline'; img-src 'self' data:",
};

function enviar(res, estado, cuerpo, tipo = 'application/json; charset=utf-8') {
  const datos = typeof cuerpo === 'string' || Buffer.isBuffer(cuerpo) ? cuerpo : JSON.stringify(cuerpo);
  res.writeHead(estado, { 'Content-Type': tipo, 'Cache-Control': 'no-store', ...CABECERAS_SEGURIDAD });
  res.end(datos);
}

function leerJson(req) {
  return new Promise((resolve, reject) => {
    let tam = 0;
    const partes = [];
    req.on('data', (c) => {
      tam += c.length;
      if (tam > MAX_CUERPO) {
        reject(Object.assign(new Error('La solicitud es demasiado grande.'), { estado: 413 }));
        req.destroy();
        return;
      }
      partes.push(c);
    });
    req.on('end', () => {
      try {
        resolve(JSON.parse(Buffer.concat(partes).toString('utf8') || '{}'));
      } catch (_) {
        reject(Object.assign(new Error('El cuerpo no es JSON válido.'), { estado: 400 }));
      }
    });
    req.on('error', reject);
  });
}

function validarFragmento(f) {
  if (!f || typeof f !== 'object') throw Object.assign(new Error('Fragmento no válido.'), { estado: 400 });
  const smiles = String(f.smiles || '').trim();
  const formula = String(f.formula || '').trim();
  const nombre = f.nombre ? String(f.nombre).slice(0, 80) : undefined;
  if (!smiles) throw Object.assign(new Error('Indique la estructura (SMILES).'), { estado: 400 });
  if (smiles.length > MAX_SMILES || formula.length > 60) {
    throw Object.assign(new Error(`La estructura es demasiado larga (máximo ${MAX_SMILES} caracteres).`), { estado: 400 });
  }
  return { smiles, formula, nombre };
}

const ejemplosAnalizados = compararFragmentos(ejemplos);

const servidor = http.createServer(async (req, res) => {
  const url = new URL(req.url, 'http://localhost');
  try {
    if (req.method === 'GET' && (url.pathname === '/' || url.pathname === '/index.html')) {
      return enviar(res, 200, PAGINA, 'text/html; charset=utf-8');
    }
    if (req.method === 'GET' && url.pathname === '/salud') return enviar(res, 200, { estado: 'ok' });
    if (req.method === 'GET' && url.pathname === '/api/ejemplos') {
      return enviar(res, 200, { ejemplos, resultados: ejemplosAnalizados });
    }
    if (req.method === 'POST' && url.pathname === '/api/analizar') {
      const f = validarFragmento(await leerJson(req));
      try {
        return enviar(res, 200, analizarFragmento(f));
      } catch (e) {
        return enviar(res, 422, { error: e.message });
      }
    }
    if (req.method === 'POST' && url.pathname === '/api/comparar') {
      const lista = await leerJson(req);
      if (!Array.isArray(lista) || !lista.length || lista.length > MAX_FRAGMENTOS) {
        return enviar(res, 400, { error: `Envíe un arreglo de 1 a ${MAX_FRAGMENTOS} fragmentos.` });
      }
      return enviar(res, 200, compararFragmentos(lista.map(validarFragmento)));
    }
    return enviar(res, 404, { error: 'No encontrado.' });
  } catch (e) {
    return enviar(res, e.estado || 500, { error: e.estado ? e.message : 'Error interno del servidor.' });
  }
});

servidor.requestTimeout = 15000;
servidor.listen(PUERTO, '0.0.0.0', () => {
  console.log(`Analizador de estabilidad de cationes escuchando en el puerto ${PUERTO}`);
});

const cerrar = () => servidor.close(() => process.exit(0));
process.on('SIGTERM', cerrar);
process.on('SIGINT', cerrar);
