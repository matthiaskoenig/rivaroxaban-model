# Model for rivaroxaban absorption in the small intestine

Model `rivaroxaban_intestine`, SBML Level 3 Version 2, read from rivaroxaban\_intestine.xml, written by sbmlode 0.3.0.

Model for rivaroxaban absorption

An oral dose of a tablet (PODOSE\_riv) enters the stomach (POSTOMACH\_riv) and is emptied into the intestine in three stages (emptying&lt;k&gt;\_riv, POINTESTINE\_riv), where it dissolves into the lumen (dissolution\_riv). The dissolution of rivaroxaban, a compound of low solubility (BCS class II), saturates at a small undissolved amount (Km\_dis\_riv), its capacity is Vmax\_dis\_riv; the undissolved drug leaves the absorptive part of the intestine with the transit (ktransit\_riv) and is excreted with the feces (transit\_riv). The fraction of the dose which dissolves therefore falls with the dose: the bioavailability is 80-100 % at 10 mg and about 66 % at 20 mg in the fasted state, the AUC per dose falls from 10 to 80 mg to 40 % and the urinary excretion of unchanged drug from 42 % at 1.25 mg to 11 % of the dose at 80 mg (Kubitza2005a); the Cmax of 10, 20, 40 and 80 mg (141, 173, 234 and 316 µg/l) follows from the capacity of the dissolution.

Food (fed) increases the solubility (FED\_f\_solubility, bile micelles) and slows the gastric emptying (FED\_f\_emptying): with food the AUC of 20 mg increases by 39 % and the Cmax by 76 %, the tmax is delayed by 1.25 hr; at 10 mg the AUC increases by 28 % (Kubitza2006).

An oral solution (POSOLUTION\_riv) enters the lumen of the intestine without the dissolution; the passage of the solution from the stomach is lumped into the absorption (tmax 0.5 hr, Kubitza2005a, Weinz2009).

The dissolved rivaroxaban is absorbed (RIVABS); a fraction 1 - F\_riv\_abs is not absorbed and excreted with the feces (RIVEXC): 7 % of an oral solution of \[14C\]rivaroxaban is excreted unchanged with the feces (Weinz2009). The metabolites which are excreted with the bile are excreted with the feces (RXEXC).

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
| $\mathrm{Vext}$ | `Vext` | plasma | $1$ | l | ✓ |
| $\mathrm{Vlumen}$ | `Vlumen` | lumen | $1$ | l |  |
| $\mathrm{Vfeces}$ | `Vfeces` | feces | $1$ | l | ✓ |
| $\mathrm{Vapical}$ | `Vapical` | apical membrane (intestinal membrane enterocytes) | $\mathrm{NaN}$ | m<sup>2</sup> | ✓ |

## Species

| Symbol | Id | Name | Compartment | Value | Unit | Properties |
| --- | --- | --- | --- | --- | --- | --- |
| $\mathrm{riv}_{\mathrm{ext}}$ | `riv_ext` | rivaroxaban (plasma) | $\mathrm{Vext}$ | $0$ | mmol/l | concentration |
| $\mathrm{riv}_{\mathrm{lumen}}$ | `riv_lumen` | rivaroxaban (intestinal volume) | $\mathrm{Vlumen}$ | $0$ | mmol/l | concentration |
| $\mathrm{riv}_{\mathrm{feces}}$ | `riv_feces` | rivaroxaban (feces) | $\mathrm{Vfeces}$ |  | mmol | amount |
| $\mathrm{rx}_{\mathrm{lumen}}$ | `rx_lumen` | rivaroxaban metabolites (intestinal volume) | $\mathrm{Vlumen}$ | $0$ | mmol/l | concentration |
| $\mathrm{rx}_{\mathrm{feces}}$ | `rx_feces` | rivaroxaban metabolites (feces) | $\mathrm{Vfeces}$ |  | mmol | amount |

## Parameters

| Symbol | Id | Name | Value | Unit | Constant |
| --- | --- | --- | --- | --- | :---: |
| $F_{\mathrm{riv\_abs}}$ | `F_riv_abs` | fraction absorbed of the dissolved rivaroxaban | $0.93$ | \- | ✓ |
| $\mathrm{RIVABS}_{\mathrm{k}}$ | `RIVABS_k` | rate of rivaroxaban absorption | $0.0818778258685152$ | 1/min | ✓ |
| $\mathrm{fed}$ | `fed` | fed state (0 fasted, 1 fed) | $0$ | \- | ✓ |
| $\mathrm{FED}_{\mathrm{f\_solubility}}$ | `FED_f_solubility` | fed effect on the solubility | $1.60753941070449$ | \- | ✓ |
| $\mathrm{FED}_{\mathrm{f\_emptying}}$ | `FED_f_emptying` | fed effect on the gastric emptying | $0.202833867463379$ | \- | ✓ |
| $f_{\mathrm{solubility}}$ | `f_solubility` | factor of the solubility |  | \- |  |
| $f_{\mathrm{emptying}}$ | `f_emptying` | factor of the gastric emptying |  | \- |  |
| $\mathrm{PODOSE}_{\mathrm{riv}}$ | `PODOSE_riv` | oral dose rivaroxaban \[mg\] | $0$ | mg |  |
| $\mathrm{POSOLUTION}_{\mathrm{riv}}$ | `POSOLUTION_riv` | oral dose as solution rivaroxaban \[mg\] | $0$ | mg |  |
| $\mathrm{POSTOMACH}_{\mathrm{riv}}$ | `POSTOMACH_riv` | oral dose in stomach rivaroxaban \[mg\] |  | mg |  |
| $\mathrm{POSTOMACH1}_{\mathrm{riv}}$ | `POSTOMACH1_riv` | oral dose in stomach rivaroxaban, stage 1 \[mg\] | $0$ | mg |  |
| $\mathrm{POSTOMACH2}_{\mathrm{riv}}$ | `POSTOMACH2_riv` | oral dose in stomach rivaroxaban, stage 2 \[mg\] | $0$ | mg |  |
| $\mathrm{POSTOMACH3}_{\mathrm{riv}}$ | `POSTOMACH3_riv` | oral dose in stomach rivaroxaban, stage 3 \[mg\] | $0$ | mg |  |
| $\mathrm{POINTESTINE}_{\mathrm{riv}}$ | `POINTESTINE_riv` | undissolved oral dose in intestine rivaroxaban \[mg\] | $0$ | mg |  |
| $\mathrm{Ka}_{\mathrm{application\_riv}}$ | `Ka_application_riv` | Ka \[1/hr\] application rivaroxaban | $1000$ | 1/hr | ✓ |
| $\mathrm{application}_{\mathrm{riv}}$ | `application_riv` | PO application rivaroxaban |  | mmol/min |  |
| $\mathrm{Ka}_{\mathrm{emptying\_riv}}$ | `Ka_emptying_riv` | rate of a stage of the gastric emptying of the tablet | $11.7784725255136$ | 1/hr | ✓ |
| $\mathrm{Vmax}_{\mathrm{dis\_riv}}$ | `Vmax_dis_riv` | capacity of the dissolution rivaroxaban | $4.03658587067912$ | mg/hr | ✓ |
| $\mathrm{Km}_{\mathrm{dis\_riv}}$ | `Km_dis_riv` | undissolved amount of half-maximal dissolution rivaroxaban | $1$ | mg | ✓ |
| $\mathrm{ktransit}_{\mathrm{riv}}$ | `ktransit_riv` | rate of the transit of undissolved rivaroxaban | $0.5$ | 1/hr | ✓ |
| $\mathrm{Mr}_{\mathrm{riv}}$ | `Mr_riv` | Molecular weight rivaroxaban \[g/mole\] | $435.88$ | g/mol | ✓ |
| $\mathrm{emptying1}_{\mathrm{riv}}$ | `emptying1_riv` | gastric emptying rivaroxaban, stage 1 |  | mmol/min |  |
| $\mathrm{emptying2}_{\mathrm{riv}}$ | `emptying2_riv` | gastric emptying rivaroxaban, stage 2 |  | mmol/min |  |
| $\mathrm{emptying3}_{\mathrm{riv}}$ | `emptying3_riv` | gastric emptying rivaroxaban, stage 3 |  | mmol/min |  |
| $\mathrm{RXEXC}_{\mathrm{k}}$ | `RXEXC_k` | rate of rivaroxaban metabolites excretion | $0.000462962962962963$ | 1/min | ✓ |
| $\mathrm{absorption}_{\mathrm{riv}}$ | `absorption_riv` | absorption rivaroxaban |  | mmol/min |  |

## Initial assignments and assignment rules

The initial assignments set the values at $t = 0$:

$$
\begin{aligned}
\mathrm{riv}_{\mathrm{feces}} &= 0 \cdot \mathrm{Vfeces} \\[1ex]
\mathrm{rx}_{\mathrm{feces}} &= 0 \cdot \mathrm{Vfeces}
\end{aligned}
$$

The assignment rules hold at every time $t$:

$$
\begin{aligned}
\mathrm{absorption}_{\mathrm{riv}} &= \mathrm{RIVABS}_{\mathrm{k}} \cdot \mathrm{Vlumen} \cdot \mathrm{riv}_{\mathrm{lumen}} \\[1ex]
f_{\mathrm{solubility}} &= 1 + \mathrm{fed} \cdot \mathopen{}\left(\mathrm{FED}_{\mathrm{f\_solubility}} - 1\right) \\[1ex]
f_{\mathrm{emptying}} &= 1 + \mathrm{fed} \cdot \mathopen{}\left(\mathrm{FED}_{\mathrm{f\_emptying}} - 1\right) \\[1ex]
\mathrm{application}_{\mathrm{riv}} &= \frac{\frac{\mathrm{Ka}_{\mathrm{application\_riv}}}{60} \cdot \mathrm{PODOSE}_{\mathrm{riv}}}{\mathrm{Mr}_{\mathrm{riv}}} \\[1ex]
\mathrm{emptying1}_{\mathrm{riv}} &= \frac{\frac{\mathrm{Ka}_{\mathrm{emptying\_riv}}}{60} \cdot f_{\mathrm{emptying}} \cdot \mathrm{POSTOMACH1}_{\mathrm{riv}}}{\mathrm{Mr}_{\mathrm{riv}}} \\[1ex]
\mathrm{emptying2}_{\mathrm{riv}} &= \frac{\frac{\mathrm{Ka}_{\mathrm{emptying\_riv}}}{60} \cdot f_{\mathrm{emptying}} \cdot \mathrm{POSTOMACH2}_{\mathrm{riv}}}{\mathrm{Mr}_{\mathrm{riv}}} \\[1ex]
\mathrm{emptying3}_{\mathrm{riv}} &= \frac{\frac{\mathrm{Ka}_{\mathrm{emptying\_riv}}}{60} \cdot f_{\mathrm{emptying}} \cdot \mathrm{POSTOMACH3}_{\mathrm{riv}}}{\mathrm{Mr}_{\mathrm{riv}}} \\[1ex]
\mathrm{POSTOMACH}_{\mathrm{riv}} &= \mathrm{POSTOMACH1}_{\mathrm{riv}} + \mathrm{POSTOMACH2}_{\mathrm{riv}} + \mathrm{POSTOMACH3}_{\mathrm{riv}}
\end{aligned}
$$

## Reactions

| Rate | Id | Name | Equation |
| --- | --- | --- | --- |
| $v_{\mathrm{RIVABS}}$ | `RIVABS` | absorption rivaroxaban | $\mathrm{riv}_{\mathrm{lumen}} \longrightarrow \mathrm{riv}_{\mathrm{ext}}$ |
| $v_{\mathrm{RIVEXC}}$ | `RIVEXC` | excretion rivaroxaban (feces) | $\mathrm{riv}_{\mathrm{lumen}} \longrightarrow \mathrm{riv}_{\mathrm{feces}}$ |
| $v_{\mathrm{RXEXC}}$ | `RXEXC` | excretion rivaroxaban metabolites (feces) | $\mathrm{rx}_{\mathrm{lumen}} \longrightarrow \mathrm{rx}_{\mathrm{feces}}$ |
| $v_{\mathrm{dissolution\_riv}}$ | `dissolution_riv` | dissolution rivaroxaban | $\varnothing \longrightarrow \mathrm{riv}_{\mathrm{lumen}}$ |
| $v_{\mathrm{transit\_riv}}$ | `transit_riv` | transit undissolved rivaroxaban (feces) | $\varnothing \longrightarrow \mathrm{riv}_{\mathrm{feces}}$ |
| $v_{\mathrm{application\_solution\_riv}}$ | `application_solution_riv` | application of the solution rivaroxaban | $\varnothing \longrightarrow \mathrm{riv}_{\mathrm{lumen}}$ |

The rates of the reactions are:

$$
\begin{aligned}
v_{\mathrm{RIVABS}} &= F_{\mathrm{riv\_abs}} \cdot \mathrm{absorption}_{\mathrm{riv}} \\[1ex]
v_{\mathrm{RIVEXC}} &= \mathopen{}\left(1 - F_{\mathrm{riv\_abs}}\right) \cdot \mathrm{absorption}_{\mathrm{riv}} \\[1ex]
v_{\mathrm{RXEXC}} &= \mathrm{RXEXC}_{\mathrm{k}} \cdot \mathrm{Vlumen} \cdot \mathrm{rx}_{\mathrm{lumen}} \\[1ex]
v_{\mathrm{dissolution\_riv}} &= \frac{\frac{\frac{f_{\mathrm{solubility}} \cdot \mathrm{Vmax}_{\mathrm{dis\_riv}}}{60}}{\mathrm{Mr}_{\mathrm{riv}}} \cdot \mathrm{POINTESTINE}_{\mathrm{riv}}}{\mathrm{POINTESTINE}_{\mathrm{riv}} + \mathrm{Km}_{\mathrm{dis\_riv}}} \\[1ex]
v_{\mathrm{transit\_riv}} &= \frac{\frac{\mathrm{ktransit}_{\mathrm{riv}}}{60} \cdot \mathrm{POINTESTINE}_{\mathrm{riv}}}{\mathrm{Mr}_{\mathrm{riv}}} \\[1ex]
v_{\mathrm{application\_solution\_riv}} &= \frac{\frac{\mathrm{Ka}_{\mathrm{application\_riv}}}{60} \cdot \mathrm{POSOLUTION}_{\mathrm{riv}}}{\mathrm{Mr}_{\mathrm{riv}}}
\end{aligned}
$$

## ODE system

The states change in time with the rates of the reactions and the rate rules:

$$
\begin{aligned}
\frac{\mathrm{d} \,\mathrm{riv}_{\mathrm{ext}}}{\mathrm{d} t} &= \frac{v_{\mathrm{RIVABS}}}{\mathrm{Vext}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{riv}_{\mathrm{lumen}}}{\mathrm{d} t} &= \frac{-v_{\mathrm{RIVABS}} - v_{\mathrm{RIVEXC}} + v_{\mathrm{dissolution\_riv}} + v_{\mathrm{application\_solution\_riv}}}{\mathrm{Vlumen}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{riv}_{\mathrm{feces}}}{\mathrm{d} t} &= v_{\mathrm{RIVEXC}} + v_{\mathrm{transit\_riv}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{rx}_{\mathrm{lumen}}}{\mathrm{d} t} &= \frac{-v_{\mathrm{RXEXC}}}{\mathrm{Vlumen}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{rx}_{\mathrm{feces}}}{\mathrm{d} t} &= v_{\mathrm{RXEXC}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{PODOSE}_{\mathrm{riv}}}{\mathrm{d} t} &= -\mathrm{application}_{\mathrm{riv}} \cdot \mathrm{Mr}_{\mathrm{riv}} \qquad \text{(rate rule)} \\[1ex]
\frac{\mathrm{d} \,\mathrm{POSOLUTION}_{\mathrm{riv}}}{\mathrm{d} t} &= -v_{\mathrm{application\_solution\_riv}} \cdot \mathrm{Mr}_{\mathrm{riv}} \qquad \text{(rate rule)} \\[1ex]
\frac{\mathrm{d} \,\mathrm{POSTOMACH1}_{\mathrm{riv}}}{\mathrm{d} t} &= \mathopen{}\left(\mathrm{application}_{\mathrm{riv}} - \mathrm{emptying1}_{\mathrm{riv}}\right) \cdot \mathrm{Mr}_{\mathrm{riv}} \qquad \text{(rate rule)} \\[1ex]
\frac{\mathrm{d} \,\mathrm{POSTOMACH2}_{\mathrm{riv}}}{\mathrm{d} t} &= \mathopen{}\left(\mathrm{emptying1}_{\mathrm{riv}} - \mathrm{emptying2}_{\mathrm{riv}}\right) \cdot \mathrm{Mr}_{\mathrm{riv}} \qquad \text{(rate rule)} \\[1ex]
\frac{\mathrm{d} \,\mathrm{POSTOMACH3}_{\mathrm{riv}}}{\mathrm{d} t} &= \mathopen{}\left(\mathrm{emptying2}_{\mathrm{riv}} - \mathrm{emptying3}_{\mathrm{riv}}\right) \cdot \mathrm{Mr}_{\mathrm{riv}} \qquad \text{(rate rule)} \\[1ex]
\frac{\mathrm{d} \,\mathrm{POINTESTINE}_{\mathrm{riv}}}{\mathrm{d} t} &= \mathopen{}\left(\mathrm{emptying3}_{\mathrm{riv}} - v_{\mathrm{dissolution\_riv}} - v_{\mathrm{transit\_riv}}\right) \cdot \mathrm{Mr}_{\mathrm{riv}} \qquad \text{(rate rule)}
\end{aligned}
$$
