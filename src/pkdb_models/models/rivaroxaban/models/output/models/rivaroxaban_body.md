# Rivaroxaban body model

Model `rivaroxaban_body`, SBML Level 3 Version 2, read from rivaroxaban\_body.xml, written by sbmlode 0.3.0.

Physiologically based pharmacokinetic/ pharmacodynamic (PBPK/PD) model rivaroxaban (body model)

Rivaroxaban, a direct oral anticoagulant (DOAC), is a widely used drug for the prevention and treatment of thromboembolic disorders. While its fixed dosing regimen and predictable anticoagulation profile improve clinical usability, significant inter-individual variability in pharmacokinetics (PK) and pharmacodynamics (PD) exists, particularly in patients with altered physiology or comorbidities. This study established an open database of rivaroxaban PK and PD data from 14 clinical studies. The data were used to develop and validate a physiologically based pharmacokinetic/pharmacodynamic (PBPK/PD) model of rivaroxaban that captures its absorption, distribution, metabolism, excretion, and coagulation effects. The model incorporates the effects of dose, food intake, hepatic and renal function, and allows simulation of activated Factor X (Factor Xa) inhibition, prothrombin time (PT), and activated partial thromboplastin time (aPTT). Simulations accurately reproduced observed plasma concentration-time and coagulation profiles across a range of doses, fed/fasting conditions, as well as with different degrees of cirrhosis and renal impairment. This mechanistic PBPK/PD model provides a framework for exploring variability in rivaroxaban response across diverse clinical scenarios and supports its use in personalized dosing strategies, clinical trial design, and regulatory evaluation.

License

Source Code: MIT

Documentation: CC BY-SA 4.0

Models: CC BY-SA 4.0

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
| $\mathrm{Vre}$ | `Vre` | rest of body |  | l |  |
| $\mathrm{Vgu}$ | `Vgu` | gut |  | l |  |
| $\mathrm{Vki}$ | `Vki` | kidney |  | l |  |
| $\mathrm{Vli}$ | `Vli` | liver |  | l |  |
| $\mathrm{Vlu}$ | `Vlu` | lung |  | l |  |
| $\mathrm{Vve}$ | `Vve` | venous blood |  | l |  |
| $\mathrm{Var}$ | `Var` | arterial blood |  | l |  |
| $\mathrm{Vurine}$ | `Vurine` | urine | $1$ | l | ✓ |
| $\mathrm{Vfeces}$ | `Vfeces` | feces | $1$ | l | ✓ |
| $\mathrm{Vstomach}$ | `Vstomach` | stomach | $1$ | l | ✓ |
| $\mathrm{Vpo}$ | `Vpo` | portal plasma |  | l |  |
| $\mathrm{Vhv}$ | `Vhv` | hepatic venous plasma |  | l |  |
| $\mathrm{Vki}_{\mathrm{tissue}}$ | `Vki_tissue` | ki tissue |  | l |  |
| $\mathrm{Vki}_{\mathrm{plasma}}$ | `Vki_plasma` | ki plasma |  | l |  |
| $\mathrm{Vli}_{\mathrm{tissue}}$ | `Vli_tissue` | li tissue |  | l |  |
| $\mathrm{Vli}_{\mathrm{plasma}}$ | `Vli_plasma` | li plasma |  | l |  |
| $\mathrm{Vlu}_{\mathrm{tissue}}$ | `Vlu_tissue` | lu tissue |  | l |  |
| $\mathrm{Vlu}_{\mathrm{plasma}}$ | `Vlu_plasma` | lu plasma |  | l |  |
| $\mathrm{Vgu}_{\mathrm{tissue}}$ | `Vgu_tissue` | gu tissue |  | l |  |
| $\mathrm{Vgu}_{\mathrm{plasma}}$ | `Vgu_plasma` | gu plasma |  | l |  |
| $\mathrm{Vre}_{\mathrm{tissue}}$ | `Vre_tissue` | re tissue |  | l |  |
| $\mathrm{Vre}_{\mathrm{plasma}}$ | `Vre_plasma` | re plasma |  | l |  |
| $\mathrm{LI}_{\mathrm{\_Vmem}}$ | `LI__Vmem` | membrane | $\mathrm{NaN}$ | m<sup>2</sup> | ✓ |
| $\mathrm{LI}_{\mathrm{\_Vbi}}$ | `LI__Vbi` | bile | $1$ | l | ✓ |
| $\mathrm{GU}_{\mathrm{\_Vapical}}$ | `GU__Vapical` | apical membrane (intestinal membrane enterocytes) | $\mathrm{NaN}$ | m<sup>2</sup> | ✓ |

## Species

| Symbol | Id | Name | Compartment | Value | Unit | Properties |
| --- | --- | --- | --- | --- | --- | --- |
| $\mathrm{Cki}_{\mathrm{plasma\_riv}}$ | `Cki_plasma_riv` | rivaroxaban (kidney plasma) | $\mathrm{Vki}_{\mathrm{plasma}}$ | $0$ | mmol/l | concentration |
| $\mathrm{Cli}_{\mathrm{plasma\_riv}}$ | `Cli_plasma_riv` | rivaroxaban (liver plasma) | $\mathrm{Vli}_{\mathrm{plasma}}$ | $0$ | mmol/l | concentration |
| $\mathrm{Clu}_{\mathrm{plasma\_riv}}$ | `Clu_plasma_riv` | rivaroxaban (lung plasma) | $\mathrm{Vlu}_{\mathrm{plasma}}$ | $0$ | mmol/l | concentration |
| $\mathrm{Cgu}_{\mathrm{plasma\_riv}}$ | `Cgu_plasma_riv` | rivaroxaban (gut plasma) | $\mathrm{Vgu}_{\mathrm{plasma}}$ | $0$ | mmol/l | concentration |
| $\mathrm{Cre}_{\mathrm{plasma\_riv}}$ | `Cre_plasma_riv` | rivaroxaban (rest plasma) | $\mathrm{Vre}_{\mathrm{plasma}}$ | $0$ | mmol/l | concentration |
| $\mathrm{Car}_{\mathrm{riv}}$ | `Car_riv` | rivaroxaban (arterial blood plasma) | $\mathrm{Var}$ | $0$ | mmol/l | concentration |
| $\mathrm{Cve}_{\mathrm{riv}}$ | `Cve_riv` | rivaroxaban (venous blood plasma) | $\mathrm{Vve}$ | $0$ | mmol/l | concentration |
| $\mathrm{Cpo}_{\mathrm{riv}}$ | `Cpo_riv` | rivaroxaban (portal vein plasma) | $\mathrm{Vpo}$ | $0$ | mmol/l | concentration |
| $\mathrm{Chv}_{\mathrm{riv}}$ | `Chv_riv` | rivaroxaban (hepatic vein plasma) | $\mathrm{Vhv}$ | $0$ | mmol/l | concentration |
| $\mathrm{Clu}_{\mathrm{riv}}$ | `Clu_riv` | rivaroxaban (lung) | $\mathrm{Vlu}_{\mathrm{tissue}}$ | $0$ | mmol/l | concentration |
| $\mathrm{Cre}_{\mathrm{riv}}$ | `Cre_riv` | rivaroxaban (rest) | $\mathrm{Vre}_{\mathrm{tissue}}$ | $0$ | mmol/l | concentration |
| $\mathrm{Aurine}_{\mathrm{riv}}$ | `Aurine_riv` | rivaroxaban (urine) | $\mathrm{Vurine}$ |  | mmol | amount |
| $\mathrm{Afeces}_{\mathrm{riv}}$ | `Afeces_riv` | rivaroxaban (feces) | $\mathrm{Vfeces}$ |  | mmol | amount |
| $\mathrm{Cki}_{\mathrm{plasma\_rx}}$ | `Cki_plasma_rx` | rivaroxaban metabolites (kidney plasma) | $\mathrm{Vki}_{\mathrm{plasma}}$ | $0$ | mmol/l | concentration |
| $\mathrm{Cli}_{\mathrm{plasma\_rx}}$ | `Cli_plasma_rx` | rivaroxaban metabolites (liver plasma) | $\mathrm{Vli}_{\mathrm{plasma}}$ | $0$ | mmol/l | concentration |
| $\mathrm{Clu}_{\mathrm{plasma\_rx}}$ | `Clu_plasma_rx` | rivaroxaban metabolites (lung plasma) | $\mathrm{Vlu}_{\mathrm{plasma}}$ | $0$ | mmol/l | concentration |
| $\mathrm{Cgu}_{\mathrm{plasma\_rx}}$ | `Cgu_plasma_rx` | rivaroxaban metabolites (gut plasma) | $\mathrm{Vgu}_{\mathrm{plasma}}$ | $0$ | mmol/l | concentration |
| $\mathrm{Cre}_{\mathrm{plasma\_rx}}$ | `Cre_plasma_rx` | rivaroxaban metabolites (rest plasma) | $\mathrm{Vre}_{\mathrm{plasma}}$ | $0$ | mmol/l | concentration |
| $\mathrm{Car}_{\mathrm{rx}}$ | `Car_rx` | rivaroxaban metabolites (arterial blood plasma) | $\mathrm{Var}$ | $0$ | mmol/l | concentration |
| $\mathrm{Cve}_{\mathrm{rx}}$ | `Cve_rx` | rivaroxaban metabolites (venous blood plasma) | $\mathrm{Vve}$ | $0$ | mmol/l | concentration |
| $\mathrm{Cpo}_{\mathrm{rx}}$ | `Cpo_rx` | rivaroxaban metabolites (portal vein plasma) | $\mathrm{Vpo}$ | $0$ | mmol/l | concentration |
| $\mathrm{Chv}_{\mathrm{rx}}$ | `Chv_rx` | rivaroxaban metabolites (hepatic vein plasma) | $\mathrm{Vhv}$ | $0$ | mmol/l | concentration |
| $\mathrm{Aurine}_{\mathrm{rx}}$ | `Aurine_rx` | rivaroxaban metabolites (urine) | $\mathrm{Vurine}$ |  | mmol | amount |
| $\mathrm{Afeces}_{\mathrm{rx}}$ | `Afeces_rx` | rivaroxaban metabolites (feces) | $\mathrm{Vfeces}$ |  | mmol | amount |
| $\mathrm{Cve}_{\mathrm{riv\_rx}}$ | `Cve_riv_rx` | RIV + RX (plasma) | $\mathrm{Vve}$ |  | mmol/l | concentration |
| $\mathrm{Aurine}_{\mathrm{riv\_rx}}$ | `Aurine_riv_rx` | RIV + RX (urine) | $\mathrm{Vurine}$ |  | mmol | amount |
| $\mathrm{Afeces}_{\mathrm{riv\_rx}}$ | `Afeces_riv_rx` | RIV + RX (feces) | $\mathrm{Vfeces}$ |  | mmol | amount |
| $\mathrm{Cgu}_{\mathrm{rx}}$ | `Cgu_rx` | rivaroxaban metabolites (intestinal volume) | $\mathrm{Vgu}$ | $0$ | mmol/l | concentration |
| $\mathrm{LI}_{\mathrm{\_riv}}$ | `LI__riv` | rivaroxaban (liver) | $\mathrm{Vli}_{\mathrm{tissue}}$ | $0$ | mmol/l | concentration |
| $\mathrm{LI}_{\mathrm{\_rx}}$ | `LI__rx` | rivaroxaban metabolites (liver) | $\mathrm{Vli}_{\mathrm{tissue}}$ | $0$ | mmol/l | concentration |
| $\mathrm{LI}_{\mathrm{\_rx\_bi}}$ | `LI__rx_bi` | rivaroxaban metabolites (bile) | $\mathrm{LI}_{\mathrm{\_Vbi}}$ |  | mmol | amount |
| $\mathrm{GU}_{\mathrm{\_riv\_lumen}}$ | `GU__riv_lumen` | rivaroxaban (intestinal volume) | $\mathrm{Vgu}$ | $0$ | mmol/l | concentration |

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
| $\mathrm{BW}$ | `BW` | body weight | $75$ | kg | ✓ |
| $\mathrm{BW}_{\mathrm{ref}}$ | `BW_ref` | reference body weight \[kg\] | $75$ | kg | ✓ |
| $\mathrm{COBW}$ | `COBW` | cardiac output per bodyweight | $1.548$ | ml/s/kg | ✓ |
| $f_{\mathrm{cardiac\_function}}$ | `f_cardiac_function` | heart function | $1$ | \- |  |
| $\mathrm{CO}$ | `CO` | cardiac output |  | ml/s |  |
| $\mathrm{QC}$ | `QC` | cardiac output |  | l/min |  |
| $\mathrm{Fblood}$ | `Fblood` | blood fraction of organ volume | $0.02$ | \- |  |
| $\mathrm{HCT}$ | `HCT` | hematocrit | $0.51$ | \- | ✓ |
| $\mathrm{FVgu}$ | `FVgu` | gut fractional tissue volume | $0.0171$ | l/kg | ✓ |
| $\mathrm{FVki}$ | `FVki` | kidney fractional tissue volume | $0.0044$ | l/kg | ✓ |
| $\mathrm{FVli}$ | `FVli` | liver fractional tissue volume | $0.021$ | l/kg | ✓ |
| $\mathrm{FVlu}$ | `FVlu` | lung fractional tissue volume | $0.0076$ | l/kg | ✓ |
| $\mathrm{FVre}$ | `FVre` | rest of body fractional tissue volume |  | l/kg |  |
| $\mathrm{FVve}$ | `FVve` | venous fractional tissue volume | $0.0514$ | l/kg | ✓ |
| $\mathrm{FVar}$ | `FVar` | arterial fractional tissue volume | $0.0257$ | l/kg | ✓ |
| $\mathrm{FVpo}$ | `FVpo` | portal fractional tissue volume | $0.001$ | l/kg | ✓ |
| $\mathrm{FVhv}$ | `FVhv` | hepatic venous fractional tissue volume | $0.001$ | l/kg | ✓ |
| $\mathrm{FQgu}$ | `FQgu` | gut fractional tissue blood flow | $0.18$ | \- | ✓ |
| $\mathrm{FQki}$ | `FQki` | kidney fractional tissue blood flow | $0.19$ | \- | ✓ |
| $\mathrm{FQh}$ | `FQh` | hepatic (venous side) fractional tissue blood flow | $0.215$ | \- | ✓ |
| $\mathrm{FQlu}$ | `FQlu` | lung fractional tissue blood flow | $1$ | \- | ✓ |
| $\mathrm{FQre}$ | `FQre` | rest of body fractional tissue blood flow |  | \- |  |
| $f_{\mathrm{cirrhosis}}$ | `f_cirrhosis` | severity of cirrhosis \[0, 0.95\] | $0$ | \- | ✓ |
| $f_{\mathrm{shunts}}$ | `f_shunts` | fraction of portal venous blood shunted by the liver |  | \- |  |
| $f_{\mathrm{tissue\_loss}}$ | `f_tissue_loss` | fraction of lost parenchymal liver volume |  | \- |  |
| $\mathrm{PODOSE}_{\mathrm{riv}}$ | `PODOSE_riv` | oral dose riv \[mg\] | $0$ | mg |  |
| $\mathrm{POSOLUTION}_{\mathrm{riv}}$ | `POSOLUTION_riv` | oral dose as solution riv \[mg\] | $0$ | mg |  |
| $\mathrm{Mr}_{\mathrm{riv}}$ | `Mr_riv` | Molecular weight riv \[g/mole\] | $435.88$ | g/mol | ✓ |
| $\mathrm{ftissue}_{\mathrm{riv}}$ | `ftissue_riv` | tissue distribution riv | $5.80764575579028$ | l/min | ✓ |
| $\mathrm{Kp}_{\mathrm{riv}}$ | `Kp_riv` | tissue/plasma partition coefficient riv | $0.407068700304121$ | \- | ✓ |
| $\mathrm{IVDOSE}_{\mathrm{riv}}$ | `IVDOSE_riv` | IV bolus dose riv \[mg\] | $0$ | mg |  |
| $\mathrm{ti}_{\mathrm{riv}}$ | `ti_riv` | injection time riv \[s\] | $10$ | s | ✓ |
| $\mathrm{Ki}_{\mathrm{riv}}$ | `Ki_riv` | Ki \[1/min\] injection riv |  | 1/min |  |
| $\mathrm{Ri}_{\mathrm{riv}}$ | `Ri_riv` | Ri \[mg/min\] rate of infusion riv | $0$ | mg/min | ✓ |
| $\mathrm{Mr}_{\mathrm{rx}}$ | `Mr_rx` | Molecular weight rx \[g/mole\] | $435.88$ | g/mol | ✓ |
| $\mathrm{Qgu}$ | `Qgu` | gut blood flow |  | l/min |  |
| $\mathrm{Qki}$ | `Qki` | kidney blood flow |  | l/min |  |
| $\mathrm{Qh}$ | `Qh` | hepatic (venous side) blood flow |  | l/min |  |
| $\mathrm{Qha}$ | `Qha` | hepatic artery blood flow |  | l/min |  |
| $\mathrm{Qlu}$ | `Qlu` | lung blood flow |  | l/min |  |
| $\mathrm{Qre}$ | `Qre` | rest of body blood flow |  | l/min |  |
| $\mathrm{Qpo}$ | `Qpo` | portal blood flow |  | l/min |  |
| $\mathrm{Afeces}_{\mathrm{urine\_riv\_rx}}$ | `Afeces_urine_riv_rx` | Sum of RIV and RX feces and urine |  | mmol |  |
| $\mathrm{KI}_{\mathrm{\_f\_renal\_function}}$ | `KI__f_renal_function` | renal function | $1$ | \- | ✓ |
| $\mathrm{KI}_{\mathrm{\_GFR\_base}}$ | `KI__GFR_base` | glomerular filtration rate (healthy) | $101$ | ml/min | ✓ |
| $\mathrm{KI}_{\mathrm{\_GFR}}$ | `KI__GFR` | glomerular filtration rate |  | ml/min |  |
| $\mathrm{KI}_{\mathrm{\_fu\_riv}}$ | `KI__fu_riv` | fraction unbound in plasma rivaroxaban | $0.079$ | \- | ✓ |
| $\mathrm{KI}_{\mathrm{\_RIVEX\_f}}$ | `KI__RIVEX_f` | renal clearance of unbound rivaroxaban relative to the GFR | $6.39514735276726$ | \- | ✓ |
| $\mathrm{KI}_{\mathrm{\_f\_pgp}}$ | `KI__f_pgp` | P-gp/BCRP activity | $1$ | \- | ✓ |
| $\mathrm{KI}_{\mathrm{\_RXEX\_f}}$ | `KI__RXEX_f` | renal clearance of the rivaroxaban metabolites relative to the GFR | $7$ | \- | ✓ |
| $\mathrm{LI}_{\mathrm{\_RXEX\_BI\_k}}$ | `LI__RXEX_BI_k` | rate for rivaroxaban metabolite export in bile | $0.17953060417325$ | 1/min | ✓ |
| $\mathrm{LI}_{\mathrm{\_f\_cyp3a4}}$ | `LI__f_cyp3a4` | CYP3A4 activity | $1$ | \- | ✓ |
| $\mathrm{LI}_{\mathrm{\_fm\_cyp3a4}}$ | `LI__fm_cyp3a4` | fraction of the hepatic metabolism by CYP3A4 | $0.4$ | \- | ✓ |
| $\mathrm{LI}_{\mathrm{\_RIVIM\_k}}$ | `LI__RIVIM_k` | rate of rivaroxaban import | $1000$ | 1/min | ✓ |
| $\mathrm{LI}_{\mathrm{\_RIV2RX\_k}}$ | `LI__RIV2RX_k` | rate of rivaroxaban conversion | $0.0482107281329536$ | 1/min | ✓ |
| $\mathrm{LI}_{\mathrm{\_RXEX\_k}}$ | `LI__RXEX_k` | rate of rivaroxaban metabolite export | $1000$ | 1/min | ✓ |
| $\mathrm{GU}_{\mathrm{\_F\_riv\_abs}}$ | `GU__F_riv_abs` | fraction absorbed of the dissolved rivaroxaban | $0.93$ | \- | ✓ |
| $\mathrm{GU}_{\mathrm{\_RIVABS\_k}}$ | `GU__RIVABS_k` | rate of rivaroxaban absorption | $0.0818778258685152$ | 1/min | ✓ |
| $\mathrm{GU}_{\mathrm{\_fed}}$ | `GU__fed` | fed state (0 fasted, 1 fed) | $0$ | \- | ✓ |
| $\mathrm{GU}_{\mathrm{\_FED\_f\_solubility}}$ | `GU__FED_f_solubility` | fed effect on the solubility | $1.60753941070449$ | \- | ✓ |
| $\mathrm{GU}_{\mathrm{\_FED\_f\_emptying}}$ | `GU__FED_f_emptying` | fed effect on the gastric emptying | $0.202833867463379$ | \- | ✓ |
| $\mathrm{GU}_{\mathrm{\_f\_solubility}}$ | `GU__f_solubility` | factor of the solubility |  | \- |  |
| $\mathrm{GU}_{\mathrm{\_f\_emptying}}$ | `GU__f_emptying` | factor of the gastric emptying |  | \- |  |
| $\mathrm{GU}_{\mathrm{\_POSTOMACH\_riv}}$ | `GU__POSTOMACH_riv` | oral dose in stomach rivaroxaban \[mg\] |  | mg |  |
| $\mathrm{GU}_{\mathrm{\_POSTOMACH1\_riv}}$ | `GU__POSTOMACH1_riv` | oral dose in stomach rivaroxaban, stage 1 \[mg\] | $0$ | mg |  |
| $\mathrm{GU}_{\mathrm{\_POSTOMACH2\_riv}}$ | `GU__POSTOMACH2_riv` | oral dose in stomach rivaroxaban, stage 2 \[mg\] | $0$ | mg |  |
| $\mathrm{GU}_{\mathrm{\_POSTOMACH3\_riv}}$ | `GU__POSTOMACH3_riv` | oral dose in stomach rivaroxaban, stage 3 \[mg\] | $0$ | mg |  |
| $\mathrm{GU}_{\mathrm{\_POINTESTINE\_riv}}$ | `GU__POINTESTINE_riv` | undissolved oral dose in intestine rivaroxaban \[mg\] | $0$ | mg |  |
| $\mathrm{GU}_{\mathrm{\_Ka\_application\_riv}}$ | `GU__Ka_application_riv` | Ka \[1/hr\] application rivaroxaban | $1000$ | 1/hr | ✓ |
| $\mathrm{GU}_{\mathrm{\_application\_riv}}$ | `GU__application_riv` | PO application rivaroxaban |  | mmol/min |  |
| $\mathrm{GU}_{\mathrm{\_Ka\_emptying\_riv}}$ | `GU__Ka_emptying_riv` | rate of a stage of the gastric emptying of the tablet | $11.7784725255136$ | 1/hr | ✓ |
| $\mathrm{GU}_{\mathrm{\_Vmax\_dis\_riv}}$ | `GU__Vmax_dis_riv` | capacity of the dissolution rivaroxaban | $4.03658587067912$ | mg/hr | ✓ |
| $\mathrm{GU}_{\mathrm{\_Km\_dis\_riv}}$ | `GU__Km_dis_riv` | undissolved amount of half-maximal dissolution rivaroxaban | $1$ | mg | ✓ |
| $\mathrm{GU}_{\mathrm{\_ktransit\_riv}}$ | `GU__ktransit_riv` | rate of the transit of undissolved rivaroxaban | $0.5$ | 1/hr | ✓ |
| $\mathrm{GU}_{\mathrm{\_emptying1\_riv}}$ | `GU__emptying1_riv` | gastric emptying rivaroxaban, stage 1 |  | mmol/min |  |
| $\mathrm{GU}_{\mathrm{\_emptying2\_riv}}$ | `GU__emptying2_riv` | gastric emptying rivaroxaban, stage 2 |  | mmol/min |  |
| $\mathrm{GU}_{\mathrm{\_emptying3\_riv}}$ | `GU__emptying3_riv` | gastric emptying rivaroxaban, stage 3 |  | mmol/min |  |
| $\mathrm{GU}_{\mathrm{\_RXEXC\_k}}$ | `GU__RXEXC_k` | rate of rivaroxaban metabolites excretion | $0.000462962962962963$ | 1/min | ✓ |
| $\mathrm{GU}_{\mathrm{\_absorption\_riv}}$ | `GU__absorption_riv` | absorption rivaroxaban |  | mmol/min |  |

## Initial assignments and assignment rules

The initial assignments set the values at $t = 0$:

$$
\begin{aligned}
\mathrm{Aurine}_{\mathrm{riv}} &= 0 \cdot \mathrm{Vurine} \\[1ex]
\mathrm{Afeces}_{\mathrm{riv}} &= 0 \cdot \mathrm{Vfeces} \\[1ex]
\mathrm{Aurine}_{\mathrm{rx}} &= 0 \cdot \mathrm{Vurine} \\[1ex]
\mathrm{Afeces}_{\mathrm{rx}} &= 0 \cdot \mathrm{Vfeces} \\[1ex]
\mathrm{LI}_{\mathrm{\_rx\_bi}} &= 0 \cdot \mathrm{LI}_{\mathrm{\_Vbi}}
\end{aligned}
$$

The assignment rules hold at every time $t$:

$$
\begin{aligned}
\mathrm{PT} &= \mathrm{PT}_{\mathrm{ref}} \cdot \mathopen{}\left(1 + \mathrm{PT}_{\mathrm{slope}} \cdot \mathrm{Cve}_{\mathrm{riv}}\right) \\[1ex]
\mathrm{PT}_{\mathrm{change}} &= \mathrm{PT} - \mathrm{PT}_{\mathrm{ref}} \\[1ex]
\mathrm{PT}_{\mathrm{ratio}} &= \frac{\mathrm{PT}}{\mathrm{PT}_{\mathrm{ref}}} \\[1ex]
\mathrm{aPTT} &= \mathrm{aPTT}_{\mathrm{ref}} \cdot \mathopen{}\left(1 + \frac{\mathrm{Emax}_{\mathrm{aPTT}} \cdot \mathrm{Cve}_{\mathrm{riv}}}{\mathrm{Cve}_{\mathrm{riv}} + \mathrm{EC50}_{\mathrm{riv\_aPTT}}}\right) \\[1ex]
\mathrm{aPTT}_{\mathrm{change}} &= \mathrm{aPTT} - \mathrm{aPTT}_{\mathrm{ref}} \\[1ex]
\mathrm{aPTT}_{\mathrm{ratio}} &= \frac{\mathrm{aPTT}}{\mathrm{aPTT}_{\mathrm{ref}}} \\[1ex]
\mathrm{Xa}_{\mathrm{inhibition}} &= \frac{\mathrm{Emax}_{\mathrm{Xa}} \cdot \mathrm{Cve}_{\mathrm{riv}}}{\mathrm{Cve}_{\mathrm{riv}} + \mathrm{EC50}_{\mathrm{riv\_Xa}}} \\[1ex]
f_{\mathrm{shunts}} &= f_{\mathrm{cirrhosis}} \\[1ex]
f_{\mathrm{tissue\_loss}} &= f_{\mathrm{cirrhosis}} \\[1ex]
\mathrm{FVre} &= 1 - \mathrm{FVgu} - \mathrm{FVki} - \mathrm{FVli} \\
&\quad {} - \mathrm{FVlu} - \mathrm{FVve} - \mathrm{FVar} \\[1ex]
\mathrm{FQre} &= 1 - \mathrm{FQki} - \mathrm{FQh} \\[1ex]
\mathrm{CO} &= f_{\mathrm{cardiac\_function}} \cdot \mathrm{BW} \cdot \mathrm{COBW} \\[1ex]
\mathrm{QC} &= \frac{\mathrm{CO}}{1000} \cdot 60 \\[1ex]
\mathrm{Vgu} &= \mathrm{BW} \cdot \mathrm{FVgu} \\[1ex]
\mathrm{Vki} &= \mathrm{BW} \cdot \mathrm{FVki} \\[1ex]
\mathrm{Vli} &= \mathrm{BW} \cdot \mathrm{FVli} \\[1ex]
\mathrm{Vlu} &= \mathrm{BW} \cdot \mathrm{FVlu} \\[1ex]
\mathrm{Vre} &= \mathrm{BW} \cdot \mathrm{FVre} \\[1ex]
\mathrm{Vve} &= \mathrm{BW} \cdot \mathrm{FVve} - \frac{\mathrm{FVve}}{\mathrm{FVar} + \mathrm{FVve}} \cdot \mathrm{BW} \cdot \mathrm{Fblood} \cdot \mathopen{}\left(1 - \mathrm{FVve} - \mathrm{FVar}\right) \\[1ex]
\mathrm{Var} &= \mathrm{BW} \cdot \mathrm{FVar} - \frac{\mathrm{FVar}}{\mathrm{FVar} + \mathrm{FVve}} \cdot \mathrm{BW} \cdot \mathrm{Fblood} \cdot \mathopen{}\left(1 - \mathrm{FVve} - \mathrm{FVar}\right) \\[1ex]
\mathrm{Vpo} &= \mathopen{}\left(1 - \mathrm{HCT}\right) \cdot \mathopen{}\left(\mathrm{BW} \cdot \mathrm{FVpo} - \frac{\mathrm{FVpo}}{\mathrm{FVar} + \mathrm{FVve} + \mathrm{FVpo} + \mathrm{FVhv}} \cdot \mathrm{BW} \cdot \mathrm{Fblood} \cdot \mathopen{}\left(1 - \mathrm{FVar} - \mathrm{FVve} - \mathrm{FVpo} - \mathrm{FVhv}\right)\right) \\[1ex]
\mathrm{Vhv} &= \mathopen{}\left(1 - \mathrm{HCT}\right) \cdot \mathopen{}\left(\mathrm{BW} \cdot \mathrm{FVhv} - \frac{\mathrm{FVhv}}{\mathrm{FVar} + \mathrm{FVve} + \mathrm{FVpo} + \mathrm{FVhv}} \cdot \mathrm{BW} \cdot \mathrm{Fblood} \cdot \mathopen{}\left(1 - \mathrm{FVar} - \mathrm{FVve} - \mathrm{FVpo} - \mathrm{FVhv}\right)\right) \\[1ex]
\mathrm{Qgu} &= \mathrm{QC} \cdot \mathrm{FQgu} \\[1ex]
\mathrm{Qki} &= \mathrm{QC} \cdot \mathrm{FQki} \\[1ex]
\mathrm{Qh} &= \mathrm{QC} \cdot \mathrm{FQh} \\[1ex]
\mathrm{Qha} &= \mathrm{Qh} - \mathrm{Qgu} \\[1ex]
\mathrm{Qlu} &= \mathrm{QC} \cdot \mathrm{FQlu} \\[1ex]
\mathrm{Qre} &= \mathrm{QC} \cdot \mathrm{FQre} \\[1ex]
\mathrm{Qpo} &= \mathrm{Qgu} \\[1ex]
\mathrm{Vki}_{\mathrm{plasma}} &= \mathrm{Vki} \cdot \mathrm{Fblood} \cdot \mathopen{}\left(1 - \mathrm{HCT}\right) \\[1ex]
\mathrm{Vki}_{\mathrm{tissue}} &= \mathrm{Vki} \cdot \mathopen{}\left(1 - \mathrm{Fblood}\right) \\[1ex]
\mathrm{Vli}_{\mathrm{plasma}} &= \mathrm{Vli} \cdot \mathrm{Fblood} \cdot \mathopen{}\left(1 - \mathrm{HCT}\right) \\[1ex]
\mathrm{Vli}_{\mathrm{tissue}} &= \mathrm{Vli} \cdot \mathopen{}\left(1 - f_{\mathrm{tissue\_loss}}\right) \cdot \mathopen{}\left(1 - \mathrm{Fblood}\right) \\[1ex]
\mathrm{Vlu}_{\mathrm{plasma}} &= \mathrm{Vlu} \cdot \mathrm{Fblood} \cdot \mathopen{}\left(1 - \mathrm{HCT}\right) \\[1ex]
\mathrm{Vlu}_{\mathrm{tissue}} &= \mathrm{Vlu} \cdot \mathopen{}\left(1 - \mathrm{Fblood}\right) \\[1ex]
\mathrm{Vgu}_{\mathrm{plasma}} &= \mathrm{Vgu} \cdot \mathrm{Fblood} \cdot \mathopen{}\left(1 - \mathrm{HCT}\right) \\[1ex]
\mathrm{Vgu}_{\mathrm{tissue}} &= \mathrm{Vgu} \cdot \mathopen{}\left(1 - \mathrm{Fblood}\right) \\[1ex]
\mathrm{Vre}_{\mathrm{plasma}} &= \mathrm{Vre} \cdot \mathrm{Fblood} \cdot \mathopen{}\left(1 - \mathrm{HCT}\right) \\[1ex]
\mathrm{Vre}_{\mathrm{tissue}} &= \mathrm{Vre} \cdot \mathopen{}\left(1 - \mathrm{Fblood}\right) \\[1ex]
\mathrm{Ki}_{\mathrm{riv}} &= \frac{0.693}{\mathrm{ti}_{\mathrm{riv}}} \cdot 60 \\[1ex]
\mathrm{Cve}_{\mathrm{riv\_rx}} &= \mathrm{Cve}_{\mathrm{riv}} + \mathrm{Cve}_{\mathrm{rx}} \\[1ex]
\mathrm{Aurine}_{\mathrm{riv\_rx}} &= \mathrm{Aurine}_{\mathrm{riv}} + \mathrm{Aurine}_{\mathrm{rx}} \\[1ex]
\mathrm{Afeces}_{\mathrm{riv\_rx}} &= \mathrm{Afeces}_{\mathrm{riv}} + \mathrm{Afeces}_{\mathrm{rx}} \\[1ex]
\mathrm{Afeces}_{\mathrm{urine\_riv\_rx}} &= \mathrm{Afeces}_{\mathrm{riv\_rx}} + \mathrm{Aurine}_{\mathrm{riv\_rx}} \\[1ex]
\mathrm{KI}_{\mathrm{\_GFR}} &= \mathrm{KI}_{\mathrm{\_f\_renal\_function}} \cdot \mathrm{KI}_{\mathrm{\_GFR\_base}} \\[1ex]
\mathrm{GU}_{\mathrm{\_absorption\_riv}} &= \mathrm{GU}_{\mathrm{\_RIVABS\_k}} \cdot \mathrm{Vgu} \cdot \mathrm{GU}_{\mathrm{\_riv\_lumen}} \\[1ex]
\mathrm{GU}_{\mathrm{\_f\_solubility}} &= 1 + \mathrm{GU}_{\mathrm{\_fed}} \cdot \mathopen{}\left(\mathrm{GU}_{\mathrm{\_FED\_f\_solubility}} - 1\right) \\[1ex]
\mathrm{GU}_{\mathrm{\_f\_emptying}} &= 1 + \mathrm{GU}_{\mathrm{\_fed}} \cdot \mathopen{}\left(\mathrm{GU}_{\mathrm{\_FED\_f\_emptying}} - 1\right) \\[1ex]
\mathrm{GU}_{\mathrm{\_application\_riv}} &= \frac{\frac{\mathrm{GU}_{\mathrm{\_Ka\_application\_riv}}}{60} \cdot \mathrm{PODOSE}_{\mathrm{riv}}}{\mathrm{Mr}_{\mathrm{riv}}} \\[1ex]
\mathrm{GU}_{\mathrm{\_emptying1\_riv}} &= \frac{\frac{\mathrm{GU}_{\mathrm{\_Ka\_emptying\_riv}}}{60} \cdot \mathrm{GU}_{\mathrm{\_f\_emptying}} \cdot \mathrm{GU}_{\mathrm{\_POSTOMACH1\_riv}}}{\mathrm{Mr}_{\mathrm{riv}}} \\[1ex]
\mathrm{GU}_{\mathrm{\_emptying2\_riv}} &= \frac{\frac{\mathrm{GU}_{\mathrm{\_Ka\_emptying\_riv}}}{60} \cdot \mathrm{GU}_{\mathrm{\_f\_emptying}} \cdot \mathrm{GU}_{\mathrm{\_POSTOMACH2\_riv}}}{\mathrm{Mr}_{\mathrm{riv}}} \\[1ex]
\mathrm{GU}_{\mathrm{\_emptying3\_riv}} &= \frac{\frac{\mathrm{GU}_{\mathrm{\_Ka\_emptying\_riv}}}{60} \cdot \mathrm{GU}_{\mathrm{\_f\_emptying}} \cdot \mathrm{GU}_{\mathrm{\_POSTOMACH3\_riv}}}{\mathrm{Mr}_{\mathrm{riv}}} \\[1ex]
\mathrm{GU}_{\mathrm{\_POSTOMACH\_riv}} &= \mathrm{GU}_{\mathrm{\_POSTOMACH1\_riv}} + \mathrm{GU}_{\mathrm{\_POSTOMACH2\_riv}} + \mathrm{GU}_{\mathrm{\_POSTOMACH3\_riv}}
\end{aligned}
$$

## Reactions

| Rate | Id | Name | Equation |
| --- | --- | --- | --- |
| $v_{\mathrm{transport\_lu\_riv}}$ | `transport_lu_riv` | transport rivaroxaban | $\mathrm{Clu}_{\mathrm{plasma\_riv}} \rightleftharpoons \mathrm{Clu}_{\mathrm{riv}}$ |
| $v_{\mathrm{transport\_re\_riv}}$ | `transport_re_riv` | transport rivaroxaban | $\mathrm{Cre}_{\mathrm{plasma\_riv}} \rightleftharpoons \mathrm{Cre}_{\mathrm{riv}}$ |
| $v_{\mathrm{iv\_riv}}$ | `iv_riv` | iv rivaroxaban | $\varnothing \longrightarrow \mathrm{Cve}_{\mathrm{riv}}$ |
| $v_{\mathrm{Flow\_ar\_ki\_riv}}$ | `Flow_ar_ki_riv` | inflow kidney rivaroxaban | $\mathrm{Car}_{\mathrm{riv}} \longrightarrow \mathrm{Cki}_{\mathrm{plasma\_riv}}$ |
| $v_{\mathrm{Flow\_ki\_ve\_riv}}$ | `Flow_ki_ve_riv` | outflow kidney rivaroxaban | $\mathrm{Cki}_{\mathrm{plasma\_riv}} \longrightarrow \mathrm{Cve}_{\mathrm{riv}}$ |
| $v_{\mathrm{Flow\_arli\_li\_riv}}$ | `Flow_arli_li_riv` | arterial inflow liver rivaroxaban | $\mathrm{Car}_{\mathrm{riv}} \longrightarrow \mathrm{Cli}_{\mathrm{plasma\_riv}}$ |
| $v_{\mathrm{Flow\_arli\_hv\_riv}}$ | `Flow_arli_hv_riv` | flow arterial shunts | $\mathrm{Car}_{\mathrm{riv}} \longrightarrow \mathrm{Chv}_{\mathrm{riv}}$ |
| $v_{\mathrm{Flow\_po\_li\_riv}}$ | `Flow_po_li_riv` | outflow po rivaroxaban | $\mathrm{Cpo}_{\mathrm{riv}} \longrightarrow \mathrm{Cli}_{\mathrm{plasma\_riv}}$ |
| $v_{\mathrm{Flow\_po\_hv\_riv}}$ | `Flow_po_hv_riv` | flow portal shunts | $\mathrm{Cpo}_{\mathrm{riv}} \longrightarrow \mathrm{Chv}_{\mathrm{riv}}$ |
| $v_{\mathrm{Flow\_li\_hv\_riv}}$ | `Flow_li_hv_riv` | outflow liver rivaroxaban | $\mathrm{Cli}_{\mathrm{plasma\_riv}} \longrightarrow \mathrm{Chv}_{\mathrm{riv}}$ |
| $v_{\mathrm{Flow\_hv\_ve\_riv}}$ | `Flow_hv_ve_riv` | outflow hepatic vein rivaroxaban | $\mathrm{Chv}_{\mathrm{riv}} \longrightarrow \mathrm{Cve}_{\mathrm{riv}}$ |
| $v_{\mathrm{Flow\_ve\_lu\_riv}}$ | `Flow_ve_lu_riv` | inflow lung rivaroxaban | $\mathrm{Cve}_{\mathrm{riv}} \longrightarrow \mathrm{Clu}_{\mathrm{plasma\_riv}}$ |
| $v_{\mathrm{Flow\_lu\_ar\_riv}}$ | `Flow_lu_ar_riv` | outflow lung rivaroxaban | $\mathrm{Clu}_{\mathrm{plasma\_riv}} \longrightarrow \mathrm{Car}_{\mathrm{riv}}$ |
| $v_{\mathrm{Flow\_ar\_gu\_riv}}$ | `Flow_ar_gu_riv` | inflow gut rivaroxaban | $\mathrm{Car}_{\mathrm{riv}} \longrightarrow \mathrm{Cgu}_{\mathrm{plasma\_riv}}$ |
| $v_{\mathrm{Flow\_gu\_po\_riv}}$ | `Flow_gu_po_riv` | outflow gut rivaroxaban | $\mathrm{Cgu}_{\mathrm{plasma\_riv}} \longrightarrow \mathrm{Cpo}_{\mathrm{riv}}$ |
| $v_{\mathrm{Flow\_ar\_re\_riv}}$ | `Flow_ar_re_riv` | inflow rest rivaroxaban | $\mathrm{Car}_{\mathrm{riv}} \longrightarrow \mathrm{Cre}_{\mathrm{plasma\_riv}}$ |
| $v_{\mathrm{Flow\_re\_ve\_riv}}$ | `Flow_re_ve_riv` | outflow rest rivaroxaban | $\mathrm{Cre}_{\mathrm{plasma\_riv}} \longrightarrow \mathrm{Cve}_{\mathrm{riv}}$ |
| $v_{\mathrm{Flow\_ar\_ki\_rx}}$ | `Flow_ar_ki_rx` | inflow kidney rivaroxaban metabolites | $\mathrm{Car}_{\mathrm{rx}} \longrightarrow \mathrm{Cki}_{\mathrm{plasma\_rx}}$ |
| $v_{\mathrm{Flow\_ki\_ve\_rx}}$ | `Flow_ki_ve_rx` | outflow kidney rivaroxaban metabolites | $\mathrm{Cki}_{\mathrm{plasma\_rx}} \longrightarrow \mathrm{Cve}_{\mathrm{rx}}$ |
| $v_{\mathrm{Flow\_arli\_li\_rx}}$ | `Flow_arli_li_rx` | arterial inflow liver rivaroxaban metabolites | $\mathrm{Car}_{\mathrm{rx}} \longrightarrow \mathrm{Cli}_{\mathrm{plasma\_rx}}$ |
| $v_{\mathrm{Flow\_arli\_hv\_rx}}$ | `Flow_arli_hv_rx` | flow arterial shunts | $\mathrm{Car}_{\mathrm{rx}} \longrightarrow \mathrm{Chv}_{\mathrm{rx}}$ |
| $v_{\mathrm{Flow\_po\_li\_rx}}$ | `Flow_po_li_rx` | outflow po rivaroxaban metabolites | $\mathrm{Cpo}_{\mathrm{rx}} \longrightarrow \mathrm{Cli}_{\mathrm{plasma\_rx}}$ |
| $v_{\mathrm{Flow\_po\_hv\_rx}}$ | `Flow_po_hv_rx` | flow portal shunts | $\mathrm{Cpo}_{\mathrm{rx}} \longrightarrow \mathrm{Chv}_{\mathrm{rx}}$ |
| $v_{\mathrm{Flow\_li\_hv\_rx}}$ | `Flow_li_hv_rx` | outflow liver rivaroxaban metabolites | $\mathrm{Cli}_{\mathrm{plasma\_rx}} \longrightarrow \mathrm{Chv}_{\mathrm{rx}}$ |
| $v_{\mathrm{Flow\_hv\_ve\_rx}}$ | `Flow_hv_ve_rx` | outflow hepatic vein rivaroxaban metabolites | $\mathrm{Chv}_{\mathrm{rx}} \longrightarrow \mathrm{Cve}_{\mathrm{rx}}$ |
| $v_{\mathrm{Flow\_ve\_lu\_rx}}$ | `Flow_ve_lu_rx` | inflow lung rivaroxaban metabolites | $\mathrm{Cve}_{\mathrm{rx}} \longrightarrow \mathrm{Clu}_{\mathrm{plasma\_rx}}$ |
| $v_{\mathrm{Flow\_lu\_ar\_rx}}$ | `Flow_lu_ar_rx` | outflow lung rivaroxaban metabolites | $\mathrm{Clu}_{\mathrm{plasma\_rx}} \longrightarrow \mathrm{Car}_{\mathrm{rx}}$ |
| $v_{\mathrm{Flow\_ar\_gu\_rx}}$ | `Flow_ar_gu_rx` | inflow gut rivaroxaban metabolites | $\mathrm{Car}_{\mathrm{rx}} \longrightarrow \mathrm{Cgu}_{\mathrm{plasma\_rx}}$ |
| $v_{\mathrm{Flow\_gu\_po\_rx}}$ | `Flow_gu_po_rx` | outflow gut rivaroxaban metabolites | $\mathrm{Cgu}_{\mathrm{plasma\_rx}} \longrightarrow \mathrm{Cpo}_{\mathrm{rx}}$ |
| $v_{\mathrm{Flow\_ar\_re\_rx}}$ | `Flow_ar_re_rx` | inflow rest rivaroxaban metabolites | $\mathrm{Car}_{\mathrm{rx}} \longrightarrow \mathrm{Cre}_{\mathrm{plasma\_rx}}$ |
| $v_{\mathrm{Flow\_re\_ve\_rx}}$ | `Flow_re_ve_rx` | outflow rest rivaroxaban metabolites | $\mathrm{Cre}_{\mathrm{plasma\_rx}} \longrightarrow \mathrm{Cve}_{\mathrm{rx}}$ |
| $v_{\mathrm{KI\_\_RIVEX}}$ | `KI__RIVEX` | rivaroxaban excretion (RIVEX) | $\mathrm{Cki}_{\mathrm{plasma\_riv}} \longrightarrow \mathrm{Aurine}_{\mathrm{riv}}$ |
| $v_{\mathrm{KI\_\_RXEX}}$ | `KI__RXEX` | rivaroxaban metabolites excretion (RXEX) | $\mathrm{Cki}_{\mathrm{plasma\_rx}} \longrightarrow \mathrm{Aurine}_{\mathrm{rx}}$ |
| $v_{\mathrm{LI\_\_RIVIM}}$ | `LI__RIVIM` | rivaroxaban import (RIVIM) | $\mathrm{Cli}_{\mathrm{plasma\_riv}} \rightleftharpoons \mathrm{LI}_{\mathrm{\_riv}}$ |
| $v_{\mathrm{LI\_\_RIV2RX}}$ | `LI__RIV2RX` | rivaroxaban conversion (RIV2RX) | $\mathrm{LI}_{\mathrm{\_riv}} \longrightarrow \mathrm{LI}_{\mathrm{\_rx}}$ |
| $v_{\mathrm{LI\_\_RXEX}}$ | `LI__RXEX` | rivaroxaban metabolite export (RXEX) | $\mathrm{LI}_{\mathrm{\_rx}} \rightleftharpoons \mathrm{Cli}_{\mathrm{plasma\_rx}}$ |
| $v_{\mathrm{LI\_\_RXEX\_BI}}$ | `LI__RXEX_BI` | rivaroxaban metabolite bile export | $\mathrm{LI}_{\mathrm{\_rx}} \longrightarrow \mathrm{LI}_{\mathrm{\_rx\_bi}}$ |
| $v_{\mathrm{LI\_\_RXEXEHC}}$ | `LI__RXEXEHC` | rivaroxaban metabolite enterohepatic circulation | $\mathrm{LI}_{\mathrm{\_rx\_bi}} \longrightarrow \mathrm{Cgu}_{\mathrm{rx}}$ |
| $v_{\mathrm{GU\_\_RIVABS}}$ | `GU__RIVABS` | absorption rivaroxaban | $\mathrm{GU}_{\mathrm{\_riv\_lumen}} \longrightarrow \mathrm{Cgu}_{\mathrm{plasma\_riv}}$ |
| $v_{\mathrm{GU\_\_RIVEXC}}$ | `GU__RIVEXC` | excretion rivaroxaban (feces) | $\mathrm{GU}_{\mathrm{\_riv\_lumen}} \longrightarrow \mathrm{Afeces}_{\mathrm{riv}}$ |
| $v_{\mathrm{GU\_\_RXEXC}}$ | `GU__RXEXC` | excretion rivaroxaban metabolites (feces) | $\mathrm{Cgu}_{\mathrm{rx}} \longrightarrow \mathrm{Afeces}_{\mathrm{rx}}$ |
| $v_{\mathrm{GU\_\_dissolution\_riv}}$ | `GU__dissolution_riv` | dissolution rivaroxaban | $\varnothing \longrightarrow \mathrm{GU}_{\mathrm{\_riv\_lumen}}$ |
| $v_{\mathrm{GU\_\_transit\_riv}}$ | `GU__transit_riv` | transit undissolved rivaroxaban (feces) | $\varnothing \longrightarrow \mathrm{Afeces}_{\mathrm{riv}}$ |
| $v_{\mathrm{GU\_\_application\_solution\_riv}}$ | `GU__application_solution_riv` | application of the solution rivaroxaban | $\varnothing \longrightarrow \mathrm{GU}_{\mathrm{\_riv\_lumen}}$ |

The rates of the reactions are:

$$
\begin{aligned}
v_{\mathrm{transport\_lu\_riv}} &= \frac{\mathrm{ftissue}_{\mathrm{riv}} \cdot \mathrm{BW}}{\mathrm{BW}_{\mathrm{ref}}} \cdot \mathopen{}\left(\mathrm{Clu}_{\mathrm{plasma\_riv}} \cdot \mathrm{Kp}_{\mathrm{riv}} - \mathrm{Clu}_{\mathrm{riv}}\right) \\[1ex]
v_{\mathrm{transport\_re\_riv}} &= \frac{\mathrm{ftissue}_{\mathrm{riv}} \cdot \mathrm{BW}}{\mathrm{BW}_{\mathrm{ref}}} \cdot \mathopen{}\left(\mathrm{Cre}_{\mathrm{plasma\_riv}} \cdot \mathrm{Kp}_{\mathrm{riv}} - \mathrm{Cre}_{\mathrm{riv}}\right) \\[1ex]
v_{\mathrm{iv\_riv}} &= \frac{\mathrm{Ki}_{\mathrm{riv}} \cdot \mathrm{IVDOSE}_{\mathrm{riv}}}{\mathrm{Mr}_{\mathrm{riv}}} \\[1ex]
v_{\mathrm{Flow\_ar\_ki\_riv}} &= \mathrm{Qki} \cdot \mathrm{Car}_{\mathrm{riv}} \\[1ex]
v_{\mathrm{Flow\_ki\_ve\_riv}} &= \mathrm{Qki} \cdot \mathrm{Cki}_{\mathrm{plasma\_riv}} \\[1ex]
v_{\mathrm{Flow\_arli\_li\_riv}} &= \mathopen{}\left(1 - f_{\mathrm{shunts}}\right) \cdot \mathrm{Qha} \cdot \mathrm{Car}_{\mathrm{riv}} \\[1ex]
v_{\mathrm{Flow\_arli\_hv\_riv}} &= f_{\mathrm{shunts}} \cdot \mathrm{Qha} \cdot \mathrm{Car}_{\mathrm{riv}} \\[1ex]
v_{\mathrm{Flow\_po\_li\_riv}} &= \mathopen{}\left(1 - f_{\mathrm{shunts}}\right) \cdot \mathrm{Qpo} \cdot \mathrm{Cpo}_{\mathrm{riv}} \\[1ex]
v_{\mathrm{Flow\_po\_hv\_riv}} &= f_{\mathrm{shunts}} \cdot \mathrm{Qpo} \cdot \mathrm{Cpo}_{\mathrm{riv}} \\[1ex]
v_{\mathrm{Flow\_li\_hv\_riv}} &= \mathopen{}\left(1 - f_{\mathrm{shunts}}\right) \cdot \mathopen{}\left(\mathrm{Qpo} + \mathrm{Qha}\right) \cdot \mathrm{Cli}_{\mathrm{plasma\_riv}} \\[1ex]
v_{\mathrm{Flow\_hv\_ve\_riv}} &= \mathrm{Qh} \cdot \mathrm{Chv}_{\mathrm{riv}} \\[1ex]
v_{\mathrm{Flow\_ve\_lu\_riv}} &= \mathrm{Qlu} \cdot \mathrm{Cve}_{\mathrm{riv}} \\[1ex]
v_{\mathrm{Flow\_lu\_ar\_riv}} &= \mathrm{Qlu} \cdot \mathrm{Clu}_{\mathrm{plasma\_riv}} \\[1ex]
v_{\mathrm{Flow\_ar\_gu\_riv}} &= \mathrm{Qgu} \cdot \mathrm{Car}_{\mathrm{riv}} \\[1ex]
v_{\mathrm{Flow\_gu\_po\_riv}} &= \mathrm{Qgu} \cdot \mathrm{Cgu}_{\mathrm{plasma\_riv}} \\[1ex]
v_{\mathrm{Flow\_ar\_re\_riv}} &= \mathrm{Qre} \cdot \mathrm{Car}_{\mathrm{riv}} \\[1ex]
v_{\mathrm{Flow\_re\_ve\_riv}} &= \mathrm{Qre} \cdot \mathrm{Cre}_{\mathrm{plasma\_riv}} \\[1ex]
v_{\mathrm{Flow\_ar\_ki\_rx}} &= \mathrm{Qki} \cdot \mathrm{Car}_{\mathrm{rx}} \\[1ex]
v_{\mathrm{Flow\_ki\_ve\_rx}} &= \mathrm{Qki} \cdot \mathrm{Cki}_{\mathrm{plasma\_rx}} \\[1ex]
v_{\mathrm{Flow\_arli\_li\_rx}} &= \mathopen{}\left(1 - f_{\mathrm{shunts}}\right) \cdot \mathrm{Qha} \cdot \mathrm{Car}_{\mathrm{rx}} \\[1ex]
v_{\mathrm{Flow\_arli\_hv\_rx}} &= f_{\mathrm{shunts}} \cdot \mathrm{Qha} \cdot \mathrm{Car}_{\mathrm{rx}} \\[1ex]
v_{\mathrm{Flow\_po\_li\_rx}} &= \mathopen{}\left(1 - f_{\mathrm{shunts}}\right) \cdot \mathrm{Qpo} \cdot \mathrm{Cpo}_{\mathrm{rx}} \\[1ex]
v_{\mathrm{Flow\_po\_hv\_rx}} &= f_{\mathrm{shunts}} \cdot \mathrm{Qpo} \cdot \mathrm{Cpo}_{\mathrm{rx}} \\[1ex]
v_{\mathrm{Flow\_li\_hv\_rx}} &= \mathopen{}\left(1 - f_{\mathrm{shunts}}\right) \cdot \mathopen{}\left(\mathrm{Qpo} + \mathrm{Qha}\right) \cdot \mathrm{Cli}_{\mathrm{plasma\_rx}} \\[1ex]
v_{\mathrm{Flow\_hv\_ve\_rx}} &= \mathrm{Qh} \cdot \mathrm{Chv}_{\mathrm{rx}} \\[1ex]
v_{\mathrm{Flow\_ve\_lu\_rx}} &= \mathrm{Qlu} \cdot \mathrm{Cve}_{\mathrm{rx}} \\[1ex]
v_{\mathrm{Flow\_lu\_ar\_rx}} &= \mathrm{Qlu} \cdot \mathrm{Clu}_{\mathrm{plasma\_rx}} \\[1ex]
v_{\mathrm{Flow\_ar\_gu\_rx}} &= \mathrm{Qgu} \cdot \mathrm{Car}_{\mathrm{rx}} \\[1ex]
v_{\mathrm{Flow\_gu\_po\_rx}} &= \mathrm{Qgu} \cdot \mathrm{Cgu}_{\mathrm{plasma\_rx}} \\[1ex]
v_{\mathrm{Flow\_ar\_re\_rx}} &= \mathrm{Qre} \cdot \mathrm{Car}_{\mathrm{rx}} \\[1ex]
v_{\mathrm{Flow\_re\_ve\_rx}} &= \mathrm{Qre} \cdot \mathrm{Cre}_{\mathrm{plasma\_rx}} \\[1ex]
v_{\mathrm{KI\_\_RIVEX}} &= \mathrm{KI}_{\mathrm{\_GFR}} \cdot 0.001 \cdot \mathrm{KI}_{\mathrm{\_fu\_riv}} \cdot \mathopen{}\left(1 + \mathopen{}\left(\mathrm{KI}_{\mathrm{\_RIVEX\_f}} - 1\right) \cdot \mathrm{KI}_{\mathrm{\_f\_pgp}}\right) \cdot \mathrm{Cki}_{\mathrm{plasma\_riv}} \\[1ex]
v_{\mathrm{KI\_\_RXEX}} &= \mathrm{KI}_{\mathrm{\_GFR}} \cdot 0.001 \cdot \mathrm{KI}_{\mathrm{\_RXEX\_f}} \cdot \mathrm{Cki}_{\mathrm{plasma\_rx}} \\[1ex]
v_{\mathrm{LI\_\_RIVIM}} &= \mathrm{LI}_{\mathrm{\_RIVIM\_k}} \cdot \mathrm{Vli}_{\mathrm{tissue}} \cdot \mathopen{}\left(\mathrm{Cli}_{\mathrm{plasma\_riv}} - \mathrm{LI}_{\mathrm{\_riv}}\right) \\[1ex]
v_{\mathrm{LI\_\_RIV2RX}} &= \mathopen{}\left(1 - \mathrm{LI}_{\mathrm{\_fm\_cyp3a4}} + \mathrm{LI}_{\mathrm{\_fm\_cyp3a4}} \cdot \mathrm{LI}_{\mathrm{\_f\_cyp3a4}}\right) \cdot \mathrm{LI}_{\mathrm{\_RIV2RX\_k}} \cdot \mathrm{Vli}_{\mathrm{tissue}} \cdot \mathrm{LI}_{\mathrm{\_riv}} \\[1ex]
v_{\mathrm{LI\_\_RXEX}} &= \mathrm{LI}_{\mathrm{\_RXEX\_k}} \cdot \mathrm{Vli}_{\mathrm{tissue}} \cdot \mathopen{}\left(\mathrm{LI}_{\mathrm{\_rx}} - \mathrm{Cli}_{\mathrm{plasma\_rx}}\right) \\[1ex]
v_{\mathrm{LI\_\_RXEX\_BI}} &= \mathrm{LI}_{\mathrm{\_RXEX\_BI\_k}} \cdot \mathrm{Vli}_{\mathrm{tissue}} \cdot \mathrm{LI}_{\mathrm{\_rx}} \\[1ex]
v_{\mathrm{LI\_\_RXEXEHC}} &= v_{\mathrm{LI\_\_RXEX\_BI}} \\[1ex]
v_{\mathrm{GU\_\_RIVABS}} &= \mathrm{GU}_{\mathrm{\_F\_riv\_abs}} \cdot \mathrm{GU}_{\mathrm{\_absorption\_riv}} \\[1ex]
v_{\mathrm{GU\_\_RIVEXC}} &= \mathopen{}\left(1 - \mathrm{GU}_{\mathrm{\_F\_riv\_abs}}\right) \cdot \mathrm{GU}_{\mathrm{\_absorption\_riv}} \\[1ex]
v_{\mathrm{GU\_\_RXEXC}} &= \mathrm{GU}_{\mathrm{\_RXEXC\_k}} \cdot \mathrm{Vgu} \cdot \mathrm{Cgu}_{\mathrm{rx}} \\[1ex]
v_{\mathrm{GU\_\_dissolution\_riv}} &= \frac{\frac{\frac{\mathrm{GU}_{\mathrm{\_f\_solubility}} \cdot \mathrm{GU}_{\mathrm{\_Vmax\_dis\_riv}}}{60}}{\mathrm{Mr}_{\mathrm{riv}}} \cdot \mathrm{GU}_{\mathrm{\_POINTESTINE\_riv}}}{\mathrm{GU}_{\mathrm{\_POINTESTINE\_riv}} + \mathrm{GU}_{\mathrm{\_Km\_dis\_riv}}} \\[1ex]
v_{\mathrm{GU\_\_transit\_riv}} &= \frac{\frac{\mathrm{GU}_{\mathrm{\_ktransit\_riv}}}{60} \cdot \mathrm{GU}_{\mathrm{\_POINTESTINE\_riv}}}{\mathrm{Mr}_{\mathrm{riv}}} \\[1ex]
v_{\mathrm{GU\_\_application\_solution\_riv}} &= \frac{\frac{\mathrm{GU}_{\mathrm{\_Ka\_application\_riv}}}{60} \cdot \mathrm{POSOLUTION}_{\mathrm{riv}}}{\mathrm{Mr}_{\mathrm{riv}}}
\end{aligned}
$$

## ODE system

The states change in time with the rates of the reactions and the rate rules:

$$
\begin{aligned}
\frac{\mathrm{d} \,\mathrm{Cki}_{\mathrm{plasma\_riv}}}{\mathrm{d} t} &= \frac{v_{\mathrm{Flow\_ar\_ki\_riv}} - v_{\mathrm{Flow\_ki\_ve\_riv}} - v_{\mathrm{KI\_\_RIVEX}}}{\mathrm{Vki}_{\mathrm{plasma}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Cli}_{\mathrm{plasma\_riv}}}{\mathrm{d} t} &= \frac{v_{\mathrm{Flow\_arli\_li\_riv}} + v_{\mathrm{Flow\_po\_li\_riv}} - v_{\mathrm{Flow\_li\_hv\_riv}} - v_{\mathrm{LI\_\_RIVIM}}}{\mathrm{Vli}_{\mathrm{plasma}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Clu}_{\mathrm{plasma\_riv}}}{\mathrm{d} t} &= \frac{-v_{\mathrm{transport\_lu\_riv}} + v_{\mathrm{Flow\_ve\_lu\_riv}} - v_{\mathrm{Flow\_lu\_ar\_riv}}}{\mathrm{Vlu}_{\mathrm{plasma}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Cgu}_{\mathrm{plasma\_riv}}}{\mathrm{d} t} &= \frac{v_{\mathrm{Flow\_ar\_gu\_riv}} - v_{\mathrm{Flow\_gu\_po\_riv}} + v_{\mathrm{GU\_\_RIVABS}}}{\mathrm{Vgu}_{\mathrm{plasma}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Cre}_{\mathrm{plasma\_riv}}}{\mathrm{d} t} &= \frac{-v_{\mathrm{transport\_re\_riv}} + v_{\mathrm{Flow\_ar\_re\_riv}} - v_{\mathrm{Flow\_re\_ve\_riv}}}{\mathrm{Vre}_{\mathrm{plasma}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Car}_{\mathrm{riv}}}{\mathrm{d} t} &= \frac{1}{\mathrm{Var}} \Bigl( -v_{\mathrm{Flow\_ar\_ki\_riv}} - v_{\mathrm{Flow\_arli\_li\_riv}} - v_{\mathrm{Flow\_arli\_hv\_riv}} + v_{\mathrm{Flow\_lu\_ar\_riv}} \\
&\quad {} - v_{\mathrm{Flow\_ar\_gu\_riv}} - v_{\mathrm{Flow\_ar\_re\_riv}} \Bigr) \\[1ex]
\frac{\mathrm{d} \,\mathrm{Cve}_{\mathrm{riv}}}{\mathrm{d} t} &= \frac{1}{\mathrm{Vve}} \Bigl( v_{\mathrm{iv\_riv}} + v_{\mathrm{Flow\_ki\_ve\_riv}} + v_{\mathrm{Flow\_hv\_ve\_riv}} - v_{\mathrm{Flow\_ve\_lu\_riv}} \\
&\quad {} + v_{\mathrm{Flow\_re\_ve\_riv}} \Bigr) \\[1ex]
\frac{\mathrm{d} \,\mathrm{Cpo}_{\mathrm{riv}}}{\mathrm{d} t} &= \frac{-v_{\mathrm{Flow\_po\_li\_riv}} - v_{\mathrm{Flow\_po\_hv\_riv}} + v_{\mathrm{Flow\_gu\_po\_riv}}}{\mathrm{Vpo}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Chv}_{\mathrm{riv}}}{\mathrm{d} t} &= \frac{v_{\mathrm{Flow\_arli\_hv\_riv}} + v_{\mathrm{Flow\_po\_hv\_riv}} + v_{\mathrm{Flow\_li\_hv\_riv}} - v_{\mathrm{Flow\_hv\_ve\_riv}}}{\mathrm{Vhv}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Clu}_{\mathrm{riv}}}{\mathrm{d} t} &= \frac{v_{\mathrm{transport\_lu\_riv}}}{\mathrm{Vlu}_{\mathrm{tissue}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Cre}_{\mathrm{riv}}}{\mathrm{d} t} &= \frac{v_{\mathrm{transport\_re\_riv}}}{\mathrm{Vre}_{\mathrm{tissue}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Aurine}_{\mathrm{riv}}}{\mathrm{d} t} &= v_{\mathrm{KI\_\_RIVEX}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Afeces}_{\mathrm{riv}}}{\mathrm{d} t} &= v_{\mathrm{GU\_\_RIVEXC}} + v_{\mathrm{GU\_\_transit\_riv}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Cki}_{\mathrm{plasma\_rx}}}{\mathrm{d} t} &= \frac{v_{\mathrm{Flow\_ar\_ki\_rx}} - v_{\mathrm{Flow\_ki\_ve\_rx}} - v_{\mathrm{KI\_\_RXEX}}}{\mathrm{Vki}_{\mathrm{plasma}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Cli}_{\mathrm{plasma\_rx}}}{\mathrm{d} t} &= \frac{v_{\mathrm{Flow\_arli\_li\_rx}} + v_{\mathrm{Flow\_po\_li\_rx}} - v_{\mathrm{Flow\_li\_hv\_rx}} + v_{\mathrm{LI\_\_RXEX}}}{\mathrm{Vli}_{\mathrm{plasma}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Clu}_{\mathrm{plasma\_rx}}}{\mathrm{d} t} &= \frac{v_{\mathrm{Flow\_ve\_lu\_rx}} - v_{\mathrm{Flow\_lu\_ar\_rx}}}{\mathrm{Vlu}_{\mathrm{plasma}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Cgu}_{\mathrm{plasma\_rx}}}{\mathrm{d} t} &= \frac{v_{\mathrm{Flow\_ar\_gu\_rx}} - v_{\mathrm{Flow\_gu\_po\_rx}}}{\mathrm{Vgu}_{\mathrm{plasma}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Cre}_{\mathrm{plasma\_rx}}}{\mathrm{d} t} &= \frac{v_{\mathrm{Flow\_ar\_re\_rx}} - v_{\mathrm{Flow\_re\_ve\_rx}}}{\mathrm{Vre}_{\mathrm{plasma}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Car}_{\mathrm{rx}}}{\mathrm{d} t} &= \frac{1}{\mathrm{Var}} \Bigl( -v_{\mathrm{Flow\_ar\_ki\_rx}} - v_{\mathrm{Flow\_arli\_li\_rx}} - v_{\mathrm{Flow\_arli\_hv\_rx}} + v_{\mathrm{Flow\_lu\_ar\_rx}} \\
&\quad {} - v_{\mathrm{Flow\_ar\_gu\_rx}} - v_{\mathrm{Flow\_ar\_re\_rx}} \Bigr) \\[1ex]
\frac{\mathrm{d} \,\mathrm{Cve}_{\mathrm{rx}}}{\mathrm{d} t} &= \frac{v_{\mathrm{Flow\_ki\_ve\_rx}} + v_{\mathrm{Flow\_hv\_ve\_rx}} - v_{\mathrm{Flow\_ve\_lu\_rx}} + v_{\mathrm{Flow\_re\_ve\_rx}}}{\mathrm{Vve}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Cpo}_{\mathrm{rx}}}{\mathrm{d} t} &= \frac{-v_{\mathrm{Flow\_po\_li\_rx}} - v_{\mathrm{Flow\_po\_hv\_rx}} + v_{\mathrm{Flow\_gu\_po\_rx}}}{\mathrm{Vpo}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Chv}_{\mathrm{rx}}}{\mathrm{d} t} &= \frac{v_{\mathrm{Flow\_arli\_hv\_rx}} + v_{\mathrm{Flow\_po\_hv\_rx}} + v_{\mathrm{Flow\_li\_hv\_rx}} - v_{\mathrm{Flow\_hv\_ve\_rx}}}{\mathrm{Vhv}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Aurine}_{\mathrm{rx}}}{\mathrm{d} t} &= v_{\mathrm{KI\_\_RXEX}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Afeces}_{\mathrm{rx}}}{\mathrm{d} t} &= v_{\mathrm{GU\_\_RXEXC}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Cgu}_{\mathrm{rx}}}{\mathrm{d} t} &= \frac{v_{\mathrm{LI\_\_RXEXEHC}} - v_{\mathrm{GU\_\_RXEXC}}}{\mathrm{Vgu}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{LI}_{\mathrm{\_riv}}}{\mathrm{d} t} &= \frac{v_{\mathrm{LI\_\_RIVIM}} - v_{\mathrm{LI\_\_RIV2RX}}}{\mathrm{Vli}_{\mathrm{tissue}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{LI}_{\mathrm{\_rx}}}{\mathrm{d} t} &= \frac{v_{\mathrm{LI\_\_RIV2RX}} - v_{\mathrm{LI\_\_RXEX}} - v_{\mathrm{LI\_\_RXEX\_BI}}}{\mathrm{Vli}_{\mathrm{tissue}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{LI}_{\mathrm{\_rx\_bi}}}{\mathrm{d} t} &= v_{\mathrm{LI\_\_RXEX\_BI}} - v_{\mathrm{LI\_\_RXEXEHC}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{GU}_{\mathrm{\_riv\_lumen}}}{\mathrm{d} t} &= \frac{-v_{\mathrm{GU\_\_RIVABS}} - v_{\mathrm{GU\_\_RIVEXC}} + v_{\mathrm{GU\_\_dissolution\_riv}} + v_{\mathrm{GU\_\_application\_solution\_riv}}}{\mathrm{Vgu}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{PODOSE}_{\mathrm{riv}}}{\mathrm{d} t} &= -\mathrm{GU}_{\mathrm{\_application\_riv}} \cdot \mathrm{Mr}_{\mathrm{riv}} \qquad \text{(rate rule)} \\[1ex]
\frac{\mathrm{d} \,\mathrm{POSOLUTION}_{\mathrm{riv}}}{\mathrm{d} t} &= -v_{\mathrm{GU\_\_application\_solution\_riv}} \cdot \mathrm{Mr}_{\mathrm{riv}} \qquad \text{(rate rule)} \\[1ex]
\frac{\mathrm{d} \,\mathrm{IVDOSE}_{\mathrm{riv}}}{\mathrm{d} t} &= -v_{\mathrm{iv\_riv}} \cdot \mathrm{Mr}_{\mathrm{riv}} + \mathrm{Ri}_{\mathrm{riv}} \qquad \text{(rate rule)} \\[1ex]
\frac{\mathrm{d} \,\mathrm{GU}_{\mathrm{\_POSTOMACH1\_riv}}}{\mathrm{d} t} &= \mathopen{}\left(\mathrm{GU}_{\mathrm{\_application\_riv}} - \mathrm{GU}_{\mathrm{\_emptying1\_riv}}\right) \cdot \mathrm{Mr}_{\mathrm{riv}} \qquad \text{(rate rule)} \\[1ex]
\frac{\mathrm{d} \,\mathrm{GU}_{\mathrm{\_POSTOMACH2\_riv}}}{\mathrm{d} t} &= \mathopen{}\left(\mathrm{GU}_{\mathrm{\_emptying1\_riv}} - \mathrm{GU}_{\mathrm{\_emptying2\_riv}}\right) \cdot \mathrm{Mr}_{\mathrm{riv}} \qquad \text{(rate rule)} \\[1ex]
\frac{\mathrm{d} \,\mathrm{GU}_{\mathrm{\_POSTOMACH3\_riv}}}{\mathrm{d} t} &= \mathopen{}\left(\mathrm{GU}_{\mathrm{\_emptying2\_riv}} - \mathrm{GU}_{\mathrm{\_emptying3\_riv}}\right) \cdot \mathrm{Mr}_{\mathrm{riv}} \qquad \text{(rate rule)} \\[1ex]
\frac{\mathrm{d} \,\mathrm{GU}_{\mathrm{\_POINTESTINE\_riv}}}{\mathrm{d} t} &= \mathopen{}\left(\mathrm{GU}_{\mathrm{\_emptying3\_riv}} - v_{\mathrm{GU\_\_dissolution\_riv}} - v_{\mathrm{GU\_\_transit\_riv}}\right) \cdot \mathrm{Mr}_{\mathrm{riv}} \qquad \text{(rate rule)}
\end{aligned}
$$
