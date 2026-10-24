# Model for coagulation factors and readouts.

Model `rivaroxaban_coagulation`, SBML Level 3 Version 2, read from rivaroxaban\_coagulation.xml, written by sbmlode 0.3.0.

Model for coagulation factors and readouts.

The effects depend on the total concentration of rivaroxaban in the venous plasma Cve\_riv, as the assays measure the total drug in a plasma sample. The prothrombin time (PT) is linear in the concentration (Neoplastin, up to the concentrations of 80 mg, Kubitza2005a, Mueck2014); the activated partial thromboplastin time (aPTT) and the inhibition of factor Xa follow an Emax model of the concentration, the PT and the aPTT relative to their baseline (PT\_ref, aPTT\_ref), which an experiment sets to the baseline of its subjects.

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
| $\mathrm{Vplasma}$ | `Vplasma` | plasma | $5$ | l | ✓ |

## Species

| Symbol | Id | Name | Compartment | Value | Unit | Properties |
| --- | --- | --- | --- | --- | --- | --- |
| $\mathrm{Cve}_{\mathrm{riv}}$ | `Cve_riv` | rivaroxaban | $\mathrm{Vplasma}$ | $0$ | mmol/l | concentration |

## Parameters

| Symbol | Id | Name | Value | Unit | Constant |
| --- | --- | --- | --- | --- | :---: |
| $\mathrm{PT}_{\mathrm{ref}}$ | `PT_ref` | prothrombin time reference | $12.5$ | s | ✓ |
| $\mathrm{aPTT}_{\mathrm{ref}}$ | `aPTT_ref` | activated partial thromboplastin time reference | $28.4$ | s | ✓ |
| $\mathrm{PT}$ | `PT` | prothrombin time |  | s |  |
| $\mathrm{aPTT}$ | `aPTT` | activated partial thromboplastin time |  | s |  |
| $\mathrm{Xa}_{\mathrm{inhibition}}$ | `Xa_inhibition` | inhibition of factor Xa |  | \- |  |
| $\mathrm{PT}_{\mathrm{change}}$ | `PT_change` | prothrombin time change |  | s |  |
| $\mathrm{PT}_{\mathrm{ratio}}$ | `PT_ratio` | prothrombin time ratio |  | \- |  |
| $\mathrm{aPTT}_{\mathrm{change}}$ | `aPTT_change` | activated partial thromboplastin time change |  | s |  |
| $\mathrm{aPTT}_{\mathrm{ratio}}$ | `aPTT_ratio` | activated partial thromboplastin time ratio |  | \- |  |
| $\mathrm{PT}_{\mathrm{slope}}$ | `PT_slope` | relative increase of the PT per rivaroxaban concentration | $628.063358133609$ | l/mmol | ✓ |
| $\mathrm{Emax}_{\mathrm{aPTT}}$ | `Emax_aPTT` | maximum effect aPTT | $1.59525021000745$ | \- | ✓ |
| $\mathrm{EC50}_{\mathrm{riv\_aPTT}}$ | `EC50_riv_aPTT` | EC50 rivaroxaban aPTT | $0.000845487215970562$ | mmol/l | ✓ |
| $\mathrm{Emax}_{\mathrm{Xa}}$ | `Emax_Xa` | maximum effect Xa inhibition | $1$ | \- | ✓ |
| $\mathrm{EC50}_{\mathrm{riv\_Xa}}$ | `EC50_riv_Xa` | EC50 rivaroxaban Xa inhibition | $0.000558135544034145$ | mmol/l | ✓ |

## Initial assignments and assignment rules

The assignment rules hold at every time $t$:

$$
\begin{aligned}
\mathrm{PT} &= \mathrm{PT}_{\mathrm{ref}} \cdot \mathopen{}\left(1 + \mathrm{PT}_{\mathrm{slope}} \cdot \mathrm{Cve}_{\mathrm{riv}}\right) \\[1ex]
\mathrm{PT}_{\mathrm{change}} &= \mathrm{PT} - \mathrm{PT}_{\mathrm{ref}} \\[1ex]
\mathrm{PT}_{\mathrm{ratio}} &= \frac{\mathrm{PT}}{\mathrm{PT}_{\mathrm{ref}}} \\[1ex]
\mathrm{aPTT} &= \mathrm{aPTT}_{\mathrm{ref}} \cdot \mathopen{}\left(1 + \frac{\mathrm{Emax}_{\mathrm{aPTT}} \cdot \mathrm{Cve}_{\mathrm{riv}}}{\mathrm{Cve}_{\mathrm{riv}} + \mathrm{EC50}_{\mathrm{riv\_aPTT}}}\right) \\[1ex]
\mathrm{aPTT}_{\mathrm{change}} &= \mathrm{aPTT} - \mathrm{aPTT}_{\mathrm{ref}} \\[1ex]
\mathrm{aPTT}_{\mathrm{ratio}} &= \frac{\mathrm{aPTT}}{\mathrm{aPTT}_{\mathrm{ref}}} \\[1ex]
\mathrm{Xa}_{\mathrm{inhibition}} &= \frac{\mathrm{Emax}_{\mathrm{Xa}} \cdot \mathrm{Cve}_{\mathrm{riv}}}{\mathrm{Cve}_{\mathrm{riv}} + \mathrm{EC50}_{\mathrm{riv\_Xa}}}
\end{aligned}
$$
