# Model for renal rivaroxaban excretion.

Model `rivaroxaban_kidney`, SBML Level 3 Version 2, read from rivaroxaban\_kidney.xml, written by sbmlode 0.3.0.

Model for renal rivaroxaban excretion.

The renal clearance of rivaroxaban is the filtration of its unbound fraction, GFR \* fu\_riv, times RIVEX\_f for the active secretion by P-gp and BCRP (f\_pgp); the glomerular filtration rate is the renal function times the healthy GFR of 101 ml/min. About 36 % of an oral dose of 10 mg is excreted unchanged in the urine (Weinz2009, Kubitza2013).

The content of this model has been carefully created in a manual research effort.

Terms of use

Copyright © 2026 Matthias König.

This work is licensed under a Creative Commons Attribution-ShareAlike 4.0 International License.

Redistribution and use of any part of this model, with or without modification, are permitted provided that the following conditions are met:

Redistributions of this SBML file must retain the above copyright notice, this list of conditions and the following disclaimer.

Redistributions in a different form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

This model is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

## Units

| Time | Substance | Extent | Volume | Area | Length |
| --- | --- | --- | --- | --- | --- |
| min | mmol | mmol | l | m<sup>2</sup> | m |

## Compartments

| Symbol | Id | Name | Size | Unit | Constant |
| --- | --- | --- | --- | --- | :---: |
| $\mathrm{Vext}$ | `Vext` | plasma | $1.5$ | l | ✓ |
| $\mathrm{Vki}$ | `Vki` | kidney | $0.3$ | l | ✓ |
| $\mathrm{Vurine}$ | `Vurine` | urine | $1$ | l | ✓ |

## Species

| Symbol | Id | Name | Compartment | Value | Unit | Properties |
| --- | --- | --- | --- | --- | --- | --- |
| $\mathrm{riv}_{\mathrm{ext}}$ | `riv_ext` | rivaroxaban (plasma) | $\mathrm{Vext}$ | $0$ | mmol/l | concentration |
| $\mathrm{riv}_{\mathrm{urine}}$ | `riv_urine` | rivaroxaban (urine) | $\mathrm{Vurine}$ |  | mmol | amount |
| $\mathrm{rx}_{\mathrm{ext}}$ | `rx_ext` | rivaroxaban metabolites (plasma) | $\mathrm{Vext}$ | $0$ | mmol/l | concentration |
| $\mathrm{rx}_{\mathrm{urine}}$ | `rx_urine` | rivaroxaban metabolites (urine) | $\mathrm{Vurine}$ |  | mmol | amount |

## Parameters

| Symbol | Id | Name | Value | Unit | Constant |
| --- | --- | --- | --- | --- | :---: |
| $f_{\mathrm{renal\_function}}$ | `f_renal_function` | renal function | $1$ | \- | ✓ |
| $\mathrm{GFR}_{\mathrm{base}}$ | `GFR_base` | glomerular filtration rate (healthy) | $101$ | ml/min | ✓ |
| $\mathrm{GFR}$ | `GFR` | glomerular filtration rate |  | ml/min |  |
| $\mathrm{fu}_{\mathrm{riv}}$ | `fu_riv` | fraction unbound in plasma rivaroxaban | $0.079$ | \- | ✓ |
| $\mathrm{RIVEX}_{\mathrm{f}}$ | `RIVEX_f` | renal clearance of unbound rivaroxaban relative to the GFR | $6.39514735276726$ | \- | ✓ |
| $f_{\mathrm{pgp}}$ | `f_pgp` | P-gp/BCRP activity | $1$ | \- | ✓ |
| $\mathrm{RXEX}_{\mathrm{f}}$ | `RXEX_f` | renal clearance of the rivaroxaban metabolites relative to the GFR | $7$ | \- | ✓ |

## Initial assignments and assignment rules

The initial assignments set the values at $t = 0$:

$$
\begin{aligned}
\mathrm{riv}_{\mathrm{urine}} &= 0 \cdot \mathrm{Vurine} \\[1ex]
\mathrm{rx}_{\mathrm{urine}} &= 0 \cdot \mathrm{Vurine}
\end{aligned}
$$

The assignment rules hold at every time $t$:

$$
\begin{aligned}
\mathrm{GFR} &= f_{\mathrm{renal\_function}} \cdot \mathrm{GFR}_{\mathrm{base}}
\end{aligned}
$$

## Reactions

| Rate | Id | Name | Equation |
| --- | --- | --- | --- |
| $v_{\mathrm{RIVEX}}$ | `RIVEX` | rivaroxaban excretion (RIVEX) | $\mathrm{riv}_{\mathrm{ext}} \longrightarrow \mathrm{riv}_{\mathrm{urine}}$ |
| $v_{\mathrm{RXEX}}$ | `RXEX` | rivaroxaban metabolites excretion (RXEX) | $\mathrm{rx}_{\mathrm{ext}} \longrightarrow \mathrm{rx}_{\mathrm{urine}}$ |

The rates of the reactions are:

$$
\begin{aligned}
v_{\mathrm{RIVEX}} &= \mathrm{GFR} \cdot 0.001 \cdot \mathrm{fu}_{\mathrm{riv}} \cdot \mathopen{}\left(1 + \mathopen{}\left(\mathrm{RIVEX}_{\mathrm{f}} - 1\right) \cdot f_{\mathrm{pgp}}\right) \cdot \mathrm{riv}_{\mathrm{ext}} \\[1ex]
v_{\mathrm{RXEX}} &= \mathrm{GFR} \cdot 0.001 \cdot \mathrm{RXEX}_{\mathrm{f}} \cdot \mathrm{rx}_{\mathrm{ext}}
\end{aligned}
$$

## ODE system

The states change in time with the rates of the reactions:

$$
\begin{aligned}
\frac{\mathrm{d} \,\mathrm{riv}_{\mathrm{ext}}}{\mathrm{d} t} &= \frac{-v_{\mathrm{RIVEX}}}{\mathrm{Vext}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{riv}_{\mathrm{urine}}}{\mathrm{d} t} &= v_{\mathrm{RIVEX}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{rx}_{\mathrm{ext}}}{\mathrm{d} t} &= \frac{-v_{\mathrm{RXEX}}}{\mathrm{Vext}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{rx}_{\mathrm{urine}}}{\mathrm{d} t} &= v_{\mathrm{RXEX}}
\end{aligned}
$$
