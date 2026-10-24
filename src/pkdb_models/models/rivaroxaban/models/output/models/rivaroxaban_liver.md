# Model for hepatic rivaroxaban metabolism.

Model `rivaroxaban_liver`, SBML Level 3 Version 2, read from rivaroxaban\_liver.xml, written by sbmlode 0.3.0.

Model for hepatic rivaroxaban metabolism.

Rivaroxaban is metabolized in the liver into a variety of metabolites, which are lumped (rx): oxidation by CYP3A4 (f\_cyp3a4) and CYP2J2 and CYP independent hydrolysis (Weinz2009, Mueck2014). The metabolites leave the liver into the plasma, from which they are excreted in the urine, and with the bile into the feces.

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
| $\mathrm{Vli}$ | `Vli` | liver | $1.5$ | l | ✓ |
| $\mathrm{Vmem}$ | `Vmem` | membrane | $\mathrm{NaN}$ | m<sup>2</sup> | ✓ |
| $\mathrm{Vbi}$ | `Vbi` | bile | $1$ | l | ✓ |
| $\mathrm{Vlumen}$ | `Vlumen` | gut lumen | $1$ | l | ✓ |

## Species

| Symbol | Id | Name | Compartment | Value | Unit | Properties |
| --- | --- | --- | --- | --- | --- | --- |
| $\mathrm{riv}_{\mathrm{ext}}$ | `riv_ext` | rivaroxaban (plasma) | $\mathrm{Vext}$ | $0$ | mmol/l | concentration |
| $\mathrm{riv}$ | `riv` | rivaroxaban (liver) | $\mathrm{Vli}$ | $0$ | mmol/l | concentration |
| $\mathrm{rx}_{\mathrm{ext}}$ | `rx_ext` | rivaroxaban metabolites (plasma) | $\mathrm{Vext}$ | $0$ | mmol/l | concentration |
| $\mathrm{rx}$ | `rx` | rivaroxaban metabolites (liver) | $\mathrm{Vli}$ | $0$ | mmol/l | concentration |
| $\mathrm{rx}_{\mathrm{bi}}$ | `rx_bi` | rivaroxaban metabolites (bile) | $\mathrm{Vbi}$ |  | mmol | amount |
| $\mathrm{rx}_{\mathrm{lumen}}$ | `rx_lumen` | rivaroxaban metabolites (gut) | $\mathrm{Vlumen}$ |  | mmol | amount |

## Parameters

| Symbol | Id | Name | Value | Unit | Constant |
| --- | --- | --- | --- | --- | :---: |
| $\mathrm{RXEX}_{\mathrm{BI\_k}}$ | `RXEX_BI_k` | rate for rivaroxaban metabolite export in bile | $0.17953060417325$ | 1/min | ✓ |
| $f_{\mathrm{cyp3a4}}$ | `f_cyp3a4` | CYP3A4 activity | $1$ | \- | ✓ |
| $\mathrm{fm}_{\mathrm{cyp3a4}}$ | `fm_cyp3a4` | fraction of the hepatic metabolism by CYP3A4 | $0.4$ | \- | ✓ |
| $\mathrm{RIVIM}_{\mathrm{k}}$ | `RIVIM_k` | rate of rivaroxaban import | $1000$ | 1/min | ✓ |
| $\mathrm{RIV2RX}_{\mathrm{k}}$ | `RIV2RX_k` | rate of rivaroxaban conversion | $0.0482107281329536$ | 1/min | ✓ |
| $\mathrm{RXEX}_{\mathrm{k}}$ | `RXEX_k` | rate of rivaroxaban metabolite export | $1000$ | 1/min | ✓ |

## Initial assignments and assignment rules

The initial assignments set the values at $t = 0$:

$$
\begin{aligned}
\mathrm{rx}_{\mathrm{bi}} &= 0 \cdot \mathrm{Vbi} \\[1ex]
\mathrm{rx}_{\mathrm{lumen}} &= 0 \cdot \mathrm{Vlumen}
\end{aligned}
$$

## Reactions

| Rate | Id | Name | Equation |
| --- | --- | --- | --- |
| $v_{\mathrm{RIVIM}}$ | `RIVIM` | rivaroxaban import (RIVIM) | $\mathrm{riv}_{\mathrm{ext}} \rightleftharpoons \mathrm{riv}$ |
| $v_{\mathrm{RIV2RX}}$ | `RIV2RX` | rivaroxaban conversion (RIV2RX) | $\mathrm{riv} \longrightarrow \mathrm{rx}$ |
| $v_{\mathrm{RXEX}}$ | `RXEX` | rivaroxaban metabolite export (RXEX) | $\mathrm{rx} \rightleftharpoons \mathrm{rx}_{\mathrm{ext}}$ |
| $v_{\mathrm{RXEX\_BI}}$ | `RXEX_BI` | rivaroxaban metabolite bile export | $\mathrm{rx} \longrightarrow \mathrm{rx}_{\mathrm{bi}}$ |
| $v_{\mathrm{RXEXEHC}}$ | `RXEXEHC` | rivaroxaban metabolite enterohepatic circulation | $\mathrm{rx}_{\mathrm{bi}} \longrightarrow \mathrm{rx}_{\mathrm{lumen}}$ |

The rates of the reactions are:

$$
\begin{aligned}
v_{\mathrm{RIVIM}} &= \mathrm{RIVIM}_{\mathrm{k}} \cdot \mathrm{Vli} \cdot \mathopen{}\left(\mathrm{riv}_{\mathrm{ext}} - \mathrm{riv}\right) \\[1ex]
v_{\mathrm{RIV2RX}} &= \mathopen{}\left(1 - \mathrm{fm}_{\mathrm{cyp3a4}} + \mathrm{fm}_{\mathrm{cyp3a4}} \cdot f_{\mathrm{cyp3a4}}\right) \cdot \mathrm{RIV2RX}_{\mathrm{k}} \cdot \mathrm{Vli} \cdot \mathrm{riv} \\[1ex]
v_{\mathrm{RXEX}} &= \mathrm{RXEX}_{\mathrm{k}} \cdot \mathrm{Vli} \cdot \mathopen{}\left(\mathrm{rx} - \mathrm{rx}_{\mathrm{ext}}\right) \\[1ex]
v_{\mathrm{RXEX\_BI}} &= \mathrm{RXEX}_{\mathrm{BI\_k}} \cdot \mathrm{Vli} \cdot \mathrm{rx} \\[1ex]
v_{\mathrm{RXEXEHC}} &= v_{\mathrm{RXEX\_BI}}
\end{aligned}
$$

## ODE system

The states change in time with the rates of the reactions:

$$
\begin{aligned}
\frac{\mathrm{d} \,\mathrm{riv}_{\mathrm{ext}}}{\mathrm{d} t} &= \frac{-v_{\mathrm{RIVIM}}}{\mathrm{Vext}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{riv}}{\mathrm{d} t} &= \frac{v_{\mathrm{RIVIM}} - v_{\mathrm{RIV2RX}}}{\mathrm{Vli}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{rx}_{\mathrm{ext}}}{\mathrm{d} t} &= \frac{v_{\mathrm{RXEX}}}{\mathrm{Vext}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{rx}}{\mathrm{d} t} &= \frac{v_{\mathrm{RIV2RX}} - v_{\mathrm{RXEX}} - v_{\mathrm{RXEX\_BI}}}{\mathrm{Vli}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{rx}_{\mathrm{bi}}}{\mathrm{d} t} &= v_{\mathrm{RXEX\_BI}} - v_{\mathrm{RXEXEHC}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{rx}_{\mathrm{lumen}}}{\mathrm{d} t} &= v_{\mathrm{RXEXEHC}}
\end{aligned}
$$
