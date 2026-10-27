'use strict';

// Iones diagnósticos frecuentes en espectros EI-MS (clave: fórmula de Hill|tipo).
const IONES = {
  'CH3|EE': 'm/z 15: catión metilo CH3⁺ (pérdida o presencia de grupos metilo).',
  'C2H5|EE': 'm/z 29: catión etilo C2H5⁺ (serie alquílica).',
  'CHO|EE': 'm/z 29: ion formilo HC≡O⁺ (aldehídos).',
  'C2H3O|EE': 'm/z 43: ion acetilo CH3C≡O⁺ (metilcetonas, acetatos); muy frecuente.',
  'C3H7|EE': 'm/z 43: catión propilo/isopropilo C3H7⁺ (serie alquílica).',
  'C3H5|EE': 'm/z 41: catión alilo C3H5⁺.',
  'C3H3|EE': 'm/z 39: catión ciclopropenilo C3H3⁺ (aromático, 2 e π).',
  'C4H9|EE': 'm/z 57: catión butilo/tert-butilo C4H9⁺ (alcanos, compuestos con t-Bu).',
  'CH4N|EE': 'm/z 30: ion iminio CH2=NH2⁺ (pico base típico de aminas primarias R–CH2–NH2).',
  'C2H6N|EE': 'm/z 44: ion iminio CH3CH=NH2⁺ / CH2=NHCH3⁺ (aminas).',
  'CH3O|EE': 'm/z 31: ion oxonio CH2=OH⁺ (alcoholes primarios, éteres metílicos).',
  'CH3S|EE': 'm/z 47: ion tiocarbenio CH2=SH⁺ (tioles, sulfuros).',
  'CHO2|EE': 'm/z 45: ion COOH⁺ (ácidos carboxílicos).',
  'CH2NO|EE': 'm/z 44: ion H2N–C≡O⁺ (amidas primarias).',
  'C7H7|EE': 'm/z 91: ion tropilio/bencilo C7H7⁺ (compuestos con grupo bencilo); suele ser pico base.',
  'C6H5|EE': 'm/z 77: catión fenilo C6H5⁺ (derivados monosustituidos del benceno).',
  'C7H5O|EE': 'm/z 105: ion benzoilo C6H5C≡O⁺ (benzofenonas, ésteres benzoicos).',
  'C5H5|EE': 'm/z 65: C5H5⁺ (pérdida de C2H2 desde el ion tropilio).',
  'C6H6|OE': 'm/z 78: ion molecular del benceno C6H6⁺•.',
  'C3H6O|OE': 'm/z 58: ion enol de McLafferty CH2=C(OH)CH3⁺• (metilcetonas con H en γ).',
  'C2H4O2|OE': 'm/z 60: ion de McLafferty CH2=C(OH)2⁺• (ácidos carboxílicos).',
  'C2H4O|OE': 'm/z 44: ion de McLafferty CH2=CHOH⁺• (aldehídos con H en γ).',
};

function buscarIonConocido(hill, tipo) {
  return IONES[`${hill}|${tipo}`] || null;
}

module.exports = { buscarIonConocido };
