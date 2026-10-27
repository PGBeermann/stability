# Estabilidad de cationes fragmento en EI-MS

Programa en Node.js (sin dependencias) que estima la **estabilidad de un catión fragmento** producido por ionización electrónica (EI-MS) a partir de su **fórmula** y su **estructura** (SMILES).

Está pensado para docencia: además del resultado numérico, muestra *por qué* el catión es estable o inestable (hiperconjugación, resonancia, donación de pares libres, aromaticidad, efectos inductivos, tensión de anillo) y qué significa eso en el espectro de masas.

## Requisitos

- Node.js 18 o superior.

## Uso rápido

```bash
node bin/estabilidad-cation.js -f "C7H7+" -s "[cH+]1cccccc1" -n "Ion tropilio"
```

```bash
node bin/estabilidad-cation.js --ejemplos
```

```bash
node bin/estabilidad-cation.js --comparar ejemplos/fragmentos-tolueno.json
```

Sin argumentos se abre un **modo interactivo** que pide la fórmula y el SMILES.

| Opción | Descripción |
|---|---|
| `-f`, `--formula` | Fórmula del ion: `C7H7+`, `C2H3O+`, `(CH3)3C+`, `C6H6+•` (o `C6H6.+`) |
| `-s`, `--smiles` | Estructura en SMILES con el centro catiónico marcado |
| `-n`, `--nombre` | Nombre del fragmento (opcional) |
| `--comparar archivo.json` | Ordena varios fragmentos `[{ "nombre", "formula", "smiles" }]` |
| `--ejemplos` | Analiza y ordena fragmentos típicos de EI-MS |
| `--detalle` | Con `--ejemplos`/`--comparar`, muestra el informe de cada fragmento |
| `--json` | Salida en JSON |

> En Windows, si los acentos o símbolos (⁺, •) salen mal en la consola, ejecute `chcp 65001` antes.

## Cómo escribir la estructura (SMILES)

| Fragmento | Fórmula | SMILES |
|---|---|---|
| tert-butilo | `C4H9+` | `C[C+](C)C` |
| bencilo | `C7H7+` | `[CH2+]c1ccccc1` |
| tropilio | `C7H7+` | `[cH+]1cccccc1` |
| acetilo (acilio) | `C2H3O+` | `CC#[O+]` o `C[C+]=O` |
| benzoilo | `C7H5O+` | `[O+]#Cc1ccccc1` |
| iminio | `CH4N+` | `C=[NH2+]` o `[CH2+]N` |
| oxonio | `CH3O+` | `C=[OH+]` |
| fenilo | `C6H5+` | `[c+]1ccccc1` |
| benceno M⁺• | `C6H6+•` | `c1ccccc1` (estructura neutra = ion molecular) |

- El átomo con la carga se escribe entre corchetes con sus H: `[CH2+]`, `[C+]`, `[NH2+]`.
- Para iones radical (OE⁺•), como el ion molecular, se puede dar la estructura neutra de la molécula.
- El programa comprueba que la fórmula y la estructura coincidan y avisa si no.

## Qué calcula

**Datos de la fórmula**
- m/z monoisotópica (restando la masa del electrón) y nominal.
- RDB (anillos + dobles enlaces); valor semientero ⇒ ion EE⁺, entero ⇒ ion OE⁺•.
- Paridad electrónica (EE⁺ / OE⁺•) y regla del nitrógeno.
- Intensidades isotópicas aproximadas M+1 y M+2.
- Identificación de iones diagnósticos frecuentes (m/z 30, 31, 43, 57, 77, 91, 105…).

**Cationes de electrones pares (EE⁺)**: se estima la **afinidad por hidruro en fase gaseosa** (HIA, kcal/mol), `R⁺ + H⁻ → R–H`. **Menor HIA = catión más estable.** Se parte de un catión de referencia (CH3⁺ 312, vinilo 288, fenilo 290, HCO⁺ 260) y se suman contribuciones:

| Efecto | Ejemplo | Contribución aproximada |
|---|---|---|
| Hiperconjugación (alquilo) | etilo, isopropilo, t-butilo | −42, atenuado al acumularse |
| Resonancia alílica / bencílica | alilo, bencilo | −56 / −74 |
| Par libre de N, O, S | CH2=NH2⁺, CH2=OH⁺ | −94 a −105 (N), −62 a −70 (O) |
| Aromaticidad (Hückel 4n+2) | tropilio, ciclopropenilo | −30 |
| Antiaromaticidad (4n) | ciclopentadienilo | +45 |
| Grupos atractores en α | C=O, C≡N | +18 a +25 |
| Tensión de anillo | anillo de 3 / 4 | +20 / +10 |
| Polarizabilidad | átomos pesados en β | −1,5 cada uno |

Los efectos se combinan con saturación (no son aditivos), de modo que el modelo reproduce valores experimentales aproximados: CH3⁺ 312, C2H5⁺ 270, i-C3H7⁺ 250, t-C4H9⁺ 233, alilo 256, bencilo ≈ 238, tropilio ≈ 201, CH3CO⁺ ≈ 231, CH2=NH2⁺ ≈ 218.

La HIA se convierte en un **índice de 0 a 100** (HIA 320 → 0, HIA 190 → 100) con cinco categorías, de *muy inestable* a *muy estable*, y una interpretación para el espectro.

El programa también:
- Convierte iones onio con enlace π (C=N⁺, C=O⁺, C≡O⁺) a su forma resonante carbenio y analiza la donación del heteroátomo.
- Cuenta las formas resonantes y los átomos sobre los que se deslocaliza la carga.
- Busca **reordenamientos favorables** (migraciones 1,2 de hidruro o de alquilo): por ejemplo, n-propilo → isopropilo.

**Iones radical (OE⁺•)**: no tienen HIA. Se detectan los grupos funcionales y se estima de forma cualitativa la abundancia esperada del ion según el orden empírico de estabilidad del ion molecular en EI (aromáticos > alquenos conjugados > alicíclicos > sulfuros > alcanos lineales > cetonas > aminas > ésteres > éteres > ácidos > alcanos ramificados > alcoholes). También se indica el sitio probable de ionización, una energía de ionización aproximada y las fragmentaciones esperadas (ruptura α, bencílica, McLafferty…).

## Uso como biblioteca

```js
const { analizarFragmento, formatearInforme } = require('./src');

const r = analizarFragmento({ formula: 'C4H9+', smiles: 'C[C+](C)C' });
console.log(r.hia, r.indice, r.categoria); // 234 66 'ESTABLE'
console.log(formatearInforme(r));
```

## Pruebas

```bash
npm test
```

## Limitaciones

- Es un **modelo empírico educativo**: la HIA estimada tiene una incertidumbre típica de ±10 kcal/mol y no sustituye a cálculos cuánticos ni a datos termoquímicos.
- La intensidad real de un pico depende también de la energía de activación y de la competencia entre rutas de fragmentación (regla de Stevenson, efectos cinéticos), no solo de la estabilidad del catión.
- Lector de SMILES simplificado: no admite estereoquímica, fragmentos desconectados (`.`) ni elementos fuera de H, B, C, N, O, F, Si, P, S, Cl, Br, I.
