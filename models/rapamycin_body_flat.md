# Rapamycin body model

Model `rapamycin_body`, SBML Level 3 Version 2, read from rapamycin\_body\_flat.xml, written by sbmlode 0.3.0.

Rapamycin body model

Whole-body model of rapamycin (sirolimus) with the blood circulation, the tissue distribution and the submodels of the intestine, the liver and the kidneys.

The circulating species are whole blood concentrations: rapamycin is bound in the erythrocytes (FKBP12, blood to plasma ratio 36-142), nearly all studies measure whole blood. The blood flows transport the blood concentrations, the plasma concentration is the blood concentration over the blood to plasma ratio BP\_rap.

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
| $\mathrm{Vpo}$ | `Vpo` | portal blood |  | l |  |
| $\mathrm{Vhv}$ | `Vhv` | hepatic venous blood |  | l |  |
| $\mathrm{Vki}_{\mathrm{tissue}}$ | `Vki_tissue` | ki tissue |  | l |  |
| $\mathrm{Vki}_{\mathrm{blood}}$ | `Vki_blood` | ki blood |  | l |  |
| $\mathrm{Vli}_{\mathrm{tissue}}$ | `Vli_tissue` | li tissue |  | l |  |
| $\mathrm{Vli}_{\mathrm{blood}}$ | `Vli_blood` | li blood |  | l |  |
| $\mathrm{Vlu}_{\mathrm{tissue}}$ | `Vlu_tissue` | lu tissue |  | l |  |
| $\mathrm{Vlu}_{\mathrm{blood}}$ | `Vlu_blood` | lu blood |  | l |  |
| $\mathrm{Vgu}_{\mathrm{tissue}}$ | `Vgu_tissue` | gu tissue |  | l |  |
| $\mathrm{Vgu}_{\mathrm{blood}}$ | `Vgu_blood` | gu blood |  | l |  |
| $\mathrm{Vre}_{\mathrm{tissue}}$ | `Vre_tissue` | re tissue |  | l |  |
| $\mathrm{Vre}_{\mathrm{blood}}$ | `Vre_blood` | re blood |  | l |  |
| $\mathrm{Vgulumen}$ | `Vgulumen` | gut lumen |  | l |  |
| $\mathrm{LI}_{\mathrm{\_Vmem}}$ | `LI__Vmem` | plasma membrane | $\mathrm{NaN}$ | m<sup>2</sup> | ✓ |
| $\mathrm{LI}_{\mathrm{\_Vapical}}$ | `LI__Vapical` | apical membrane | $\mathrm{NaN}$ | m<sup>2</sup> | ✓ |
| $\mathrm{GU}_{\mathrm{\_Ventero}}$ | `GU__Ventero` | intestinal lining (enterocytes) |  | l |  |
| $\mathrm{GU}_{\mathrm{\_Vapical}}$ | `GU__Vapical` | apical membrane (intestinal membrane enterocytes) | $\mathrm{NaN}$ | m<sup>2</sup> | ✓ |
| $\mathrm{GU}_{\mathrm{\_Vbaso}}$ | `GU__Vbaso` | basolateral membrane (intestinal membrane enterocytes) | $\mathrm{NaN}$ | m<sup>2</sup> | ✓ |
| $\mathrm{GU}_{\mathrm{\_Vint\_0}}$ | `GU__Vint_0` | intestinal lumen 0 (inner part of intestine) |  | l |  |
| $\mathrm{GU}_{\mathrm{\_Vint\_1}}$ | `GU__Vint_1` | intestinal lumen 1 (inner part of intestine) |  | l |  |
| $\mathrm{GU}_{\mathrm{\_Vint\_2}}$ | `GU__Vint_2` | intestinal lumen 2 (inner part of intestine) |  | l |  |
| $\mathrm{GU}_{\mathrm{\_Vint\_3}}$ | `GU__Vint_3` | intestinal lumen 3 (inner part of intestine) |  | l |  |
| $\mathrm{GU}_{\mathrm{\_Vint\_4}}$ | `GU__Vint_4` | intestinal lumen 4 (inner part of intestine) |  | l |  |

## Species

| Symbol | Id | Name | Compartment | Value | Unit | Properties |
| --- | --- | --- | --- | --- | --- | --- |
| $\mathrm{Cki}_{\mathrm{blood\_rap}}$ | `Cki_blood_rap` | rapamycin (kidney blood) | $\mathrm{Vki}_{\mathrm{blood}}$ | $0$ | mmol/l | concentration |
| $\mathrm{Cli}_{\mathrm{blood\_rap}}$ | `Cli_blood_rap` | rapamycin (liver blood) | $\mathrm{Vli}_{\mathrm{blood}}$ | $0$ | mmol/l | concentration |
| $\mathrm{Clu}_{\mathrm{blood\_rap}}$ | `Clu_blood_rap` | rapamycin (lung blood) | $\mathrm{Vlu}_{\mathrm{blood}}$ | $0$ | mmol/l | concentration |
| $\mathrm{Cgu}_{\mathrm{blood\_rap}}$ | `Cgu_blood_rap` | rapamycin (gut blood) | $\mathrm{Vgu}_{\mathrm{blood}}$ | $0$ | mmol/l | concentration |
| $\mathrm{Cre}_{\mathrm{blood\_rap}}$ | `Cre_blood_rap` | rapamycin (rest blood) | $\mathrm{Vre}_{\mathrm{blood}}$ | $0$ | mmol/l | concentration |
| $\mathrm{Car}_{\mathrm{rap}}$ | `Car_rap` | rapamycin (arterial blood blood) | $\mathrm{Var}$ | $0$ | mmol/l | concentration |
| $\mathrm{Cve}_{\mathrm{rap}}$ | `Cve_rap` | rapamycin (venous blood blood) | $\mathrm{Vve}$ | $0$ | mmol/l | concentration |
| $\mathrm{Cpo}_{\mathrm{rap}}$ | `Cpo_rap` | rapamycin (portal vein blood) | $\mathrm{Vpo}$ | $0$ | mmol/l | concentration |
| $\mathrm{Chv}_{\mathrm{rap}}$ | `Chv_rap` | rapamycin (hepatic vein blood) | $\mathrm{Vhv}$ | $0$ | mmol/l | concentration |
| $\mathrm{Clu}_{\mathrm{rap}}$ | `Clu_rap` | rapamycin (lung) | $\mathrm{Vlu}_{\mathrm{tissue}}$ | $0$ | mmol/l | concentration |
| $\mathrm{Cre}_{\mathrm{rap}}$ | `Cre_rap` | rapamycin (rest) | $\mathrm{Vre}_{\mathrm{tissue}}$ | $0$ | mmol/l | concentration |
| $\mathrm{Cre}_{\mathrm{deep\_rap}}$ | `Cre_deep_rap` | rapamycin bound (rest) | $\mathrm{Vre}_{\mathrm{tissue}}$ | $0$ | mmol/l | concentration |
| $\mathrm{Afeces}_{\mathrm{rap}}$ | `Afeces_rap` | rapamycin (feces) | $\mathrm{Vfeces}$ |  | mmol | amount |
| $\mathrm{Clumen}_{\mathrm{rap}}$ | `Clumen_rap` | rapamycin (lumen) | $\mathrm{Vgulumen}$ | $0$ | mmol/l | concentration |
| $\mathrm{Cki}_{\mathrm{blood\_rx}}$ | `Cki_blood_rx` | rapamycin metabolites (kidney blood) | $\mathrm{Vki}_{\mathrm{blood}}$ | $0$ | mmol/l | concentration |
| $\mathrm{Cli}_{\mathrm{blood\_rx}}$ | `Cli_blood_rx` | rapamycin metabolites (liver blood) | $\mathrm{Vli}_{\mathrm{blood}}$ | $0$ | mmol/l | concentration |
| $\mathrm{Clu}_{\mathrm{blood\_rx}}$ | `Clu_blood_rx` | rapamycin metabolites (lung blood) | $\mathrm{Vlu}_{\mathrm{blood}}$ | $0$ | mmol/l | concentration |
| $\mathrm{Cgu}_{\mathrm{blood\_rx}}$ | `Cgu_blood_rx` | rapamycin metabolites (gut blood) | $\mathrm{Vgu}_{\mathrm{blood}}$ | $0$ | mmol/l | concentration |
| $\mathrm{Cre}_{\mathrm{blood\_rx}}$ | `Cre_blood_rx` | rapamycin metabolites (rest blood) | $\mathrm{Vre}_{\mathrm{blood}}$ | $0$ | mmol/l | concentration |
| $\mathrm{Car}_{\mathrm{rx}}$ | `Car_rx` | rapamycin metabolites (arterial blood blood) | $\mathrm{Var}$ | $0$ | mmol/l | concentration |
| $\mathrm{Cve}_{\mathrm{rx}}$ | `Cve_rx` | rapamycin metabolites (venous blood blood) | $\mathrm{Vve}$ | $0$ | mmol/l | concentration |
| $\mathrm{Cpo}_{\mathrm{rx}}$ | `Cpo_rx` | rapamycin metabolites (portal vein blood) | $\mathrm{Vpo}$ | $0$ | mmol/l | concentration |
| $\mathrm{Chv}_{\mathrm{rx}}$ | `Chv_rx` | rapamycin metabolites (hepatic vein blood) | $\mathrm{Vhv}$ | $0$ | mmol/l | concentration |
| $\mathrm{Aurine}_{\mathrm{rx}}$ | `Aurine_rx` | rapamycin metabolites (urine) | $\mathrm{Vurine}$ |  | mmol | amount |
| $\mathrm{Afeces}_{\mathrm{rx}}$ | `Afeces_rx` | rapamycin metabolites (feces) | $\mathrm{Vfeces}$ |  | mmol | amount |
| $\mathrm{Clumen}_{\mathrm{rx}}$ | `Clumen_rx` | rapamycin metabolites (lumen) | $\mathrm{Vgulumen}$ | $0$ | mmol/l | concentration |
| $\mathrm{Cveblood}_{\mathrm{rap}}$ | `Cveblood_rap` | rapamycin (venous blood) | $\mathrm{Vve}$ |  | mmol/l | concentration |
| $\mathrm{Cveplasma}_{\mathrm{rap}}$ | `Cveplasma_rap` | rapamycin (venous plasma) | $\mathrm{Vve}$ |  | mmol/l | concentration |
| $\mathrm{Cveblood}_{\mathrm{raptot}}$ | `Cveblood_raptot` | rapamycin and metabolites (venous blood) | $\mathrm{Vve}$ |  | mmol/l | concentration |
| $\mathrm{Cveplasma}_{\mathrm{raptot}}$ | `Cveplasma_raptot` | rapamycin and metabolites (venous plasma) | $\mathrm{Vve}$ |  | mmol/l | concentration |
| $\mathrm{Afeces}_{\mathrm{raptot}}$ | `Afeces_raptot` | rapamycin and metabolites (feces) | $\mathrm{Vfeces}$ |  | mmol | amount |
| $\mathrm{LI}_{\mathrm{\_rap}}$ | `LI__rap` | rapamycin (liver) | $\mathrm{Vli}_{\mathrm{tissue}}$ | $0$ | mmol/l | concentration |
| $\mathrm{LI}_{\mathrm{\_rx}}$ | `LI__rx` | rapamycin metabolites (liver) | $\mathrm{Vli}_{\mathrm{tissue}}$ | $0$ | mmol/l | concentration |
| $\mathrm{GU}_{\mathrm{\_rap\_stomach}}$ | `GU__rap_stomach` | rapamycin (stomach) | $\mathrm{Vstomach}$ |  | mmol | amount, boundary |
| $\mathrm{GU}_{\mathrm{\_rap\_entero}}$ | `GU__rap_entero` | rapamycin (enterocytes) | $\mathrm{GU}_{\mathrm{\_Ventero}}$ | $0$ | mmol/l | concentration |
| $\mathrm{GU}_{\mathrm{\_rx\_entero}}$ | `GU__rx_entero` | rapamycin metabolites (enterocytes) | $\mathrm{GU}_{\mathrm{\_Ventero}}$ | $0$ | mmol/l | concentration |
| $\mathrm{GU}_{\mathrm{\_rap\_colon}}$ | `GU__rap_colon` | rapamycin (colon) | $\mathrm{Vgulumen}$ | $0$ | mmol/l | concentration |
| $\mathrm{GU}_{\mathrm{\_rx\_int\_0}}$ | `GU__rx_int_0` | rapamycin metabolites (intestine) 0 | $\mathrm{GU}_{\mathrm{\_Vint\_0}}$ | $0$ | mmol/l | concentration |
| $\mathrm{GU}_{\mathrm{\_rx\_int\_1}}$ | `GU__rx_int_1` | rapamycin metabolites (intestine) 1 | $\mathrm{GU}_{\mathrm{\_Vint\_1}}$ | $0$ | mmol/l | concentration |
| $\mathrm{GU}_{\mathrm{\_rx\_int\_2}}$ | `GU__rx_int_2` | rapamycin metabolites (intestine) 2 | $\mathrm{GU}_{\mathrm{\_Vint\_2}}$ | $0$ | mmol/l | concentration |
| $\mathrm{GU}_{\mathrm{\_rx\_int\_3}}$ | `GU__rx_int_3` | rapamycin metabolites (intestine) 3 | $\mathrm{GU}_{\mathrm{\_Vint\_3}}$ | $0$ | mmol/l | concentration |
| $\mathrm{GU}_{\mathrm{\_rx\_int\_4}}$ | `GU__rx_int_4` | rapamycin metabolites (intestine) 4 | $\mathrm{GU}_{\mathrm{\_Vint\_4}}$ | $0$ | mmol/l | concentration |

## Parameters

| Symbol | Id | Name | Value | Unit | Constant |
| --- | --- | --- | --- | --- | :---: |
| $f_{\mathrm{lumen}}$ | `f_lumen` | fraction lumen of intestine | $0.9$ | \- | ✓ |
| $\mathrm{BW}$ | `BW` | body weight \[kg\] | $75$ | kg | ✓ |
| $\mathrm{BW}_{\mathrm{ref}}$ | `BW_ref` | reference body weight \[kg\] | $75$ | kg | ✓ |
| $\mathrm{HEIGHT}$ | `HEIGHT` | height \[cm\] | $170$ | cm | ✓ |
| $\mathrm{HR}$ | `HR` | heart rate \[1/min\] | $70$ | 1/min | ✓ |
| $\mathrm{HRrest}$ | `HRrest` | resting heart rate \[1/min\] | $70$ | 1/min | ✓ |
| $\mathrm{BSA}$ | `BSA` | body surface area \[m^2\] |  | m<sup>2</sup> |  |
| $\mathrm{COBW}$ | `COBW` | cardiac output per bodyweight \[ml/s/kg\] | $1.548$ | ml/s/kg | ✓ |
| $f_{\mathrm{cardiac\_function}}$ | `f_cardiac_function` | heart function | $1$ | \- |  |
| $\mathrm{CO}$ | `CO` | cardiac output \[ml/s\] |  | ml/s |  |
| $\mathrm{QC}$ | `QC` | cardiac output \[l/min\] |  | l/min |  |
| $\mathrm{COHRI}$ | `COHRI` | increase of cardiac output per heartbeat \[ml\] | $150$ | ml | ✓ |
| $\mathrm{Fblood}$ | `Fblood` | blood fraction of organ volume | $0.02$ | \- |  |
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
| $r_{\mathrm{cyp3a5}}$ | `r_cyp3a5` | capacity of CYP3A5 relative to CYP3A4 | $0.457628938029239$ | \- | ✓ |
| $\mathrm{PODOSE}_{\mathrm{rap}}$ | `PODOSE_rap` | oral dose rap \[mg\] | $0$ | mg |  |
| $\mathrm{Mr}_{\mathrm{rap}}$ | `Mr_rap` | Molecular weight rap \[g/mole\] | $914.1719$ | g/mol | ✓ |
| $\mathrm{ftissue}_{\mathrm{rap}}$ | `ftissue_rap` | tissue distribution rap | $3.71570319809393$ | l/min | ✓ |
| $\mathrm{Kp}_{\mathrm{rap}}$ | `Kp_rap` | tissue/blood partition coefficient rap | $0.648135556432825$ | \- | ✓ |
| $\mathrm{Kpdeep}_{\mathrm{rap}}$ | `Kpdeep_rap` | bound/free partition coefficient rap (tissue) | $3.50211544434062$ | \- | ✓ |
| $\mathrm{kdeep}_{\mathrm{rap}}$ | `kdeep_rap` | rate of the slow binding rap (tissue) | $0.000648410866448561$ | 1/min | ✓ |
| $\mathrm{BP}_{\mathrm{rap}}$ | `BP_rap` | blood to plasma ratio rap | $142$ | \- | ✓ |
| $\mathrm{IVDOSE}_{\mathrm{rap}}$ | `IVDOSE_rap` | IV bolus dose rap \[mg\] | $0$ | mg |  |
| $\mathrm{ti}_{\mathrm{rap}}$ | `ti_rap` | injection time rap \[s\] | $10$ | s | ✓ |
| $\mathrm{Ki}_{\mathrm{rap}}$ | `Ki_rap` | Ki \[1/min\] injection rap |  | 1/min |  |
| $\mathrm{Ri}_{\mathrm{rap}}$ | `Ri_rap` | Ri \[mg/min\] rate of infusion rap | $0$ | mg/min | ✓ |
| $\mathrm{Qgu}$ | `Qgu` | gut blood flow |  | l/min |  |
| $\mathrm{Qki}$ | `Qki` | kidney blood flow |  | l/min |  |
| $\mathrm{Qh}$ | `Qh` | hepatic (venous side) blood flow |  | l/min |  |
| $\mathrm{Qha}$ | `Qha` | hepatic artery blood flow |  | l/min |  |
| $\mathrm{Qlu}$ | `Qlu` | lung blood flow |  | l/min |  |
| $\mathrm{Qre}$ | `Qre` | rest of body blood flow |  | l/min |  |
| $\mathrm{Qpo}$ | `Qpo` | portal blood flow |  | l/min |  |
| $\mathrm{KI}_{\mathrm{\_f\_renal\_function}}$ | `KI__f_renal_function` | parameter for renal function | $1$ | \- | ✓ |
| $\mathrm{KI}_{\mathrm{\_RXEX\_k}}$ | `KI__RXEX_k` | rate urinary excretion of rapamycin metabolites | $0.12$ | 1/min | ✓ |
| $\mathrm{LI}_{\mathrm{\_RXBEX\_k}}$ | `LI__RXBEX_k` | rate of the biliary excretion of the rapamycin metabolites | $0.15$ | 1/min | ✓ |
| $\mathrm{LI}_{\mathrm{\_RAPIM\_k}}$ | `LI__RAPIM_k` | rate rapamycin import | $100$ | 1/min | ✓ |
| $\mathrm{LI}_{\mathrm{\_RAP2RX\_k}}$ | `LI__RAP2RX_k` | rate of rapamycin metabolism by CYP3A4 | $0.0261505045623823$ | 1/min | ✓ |
| $\mathrm{LI}_{\mathrm{\_f\_cyp3a4}}$ | `LI__f_cyp3a4` | CYP3A4 activity | $1$ | \- | ✓ |
| $\mathrm{LI}_{\mathrm{\_f\_cyp3a5}}$ | `LI__f_cyp3a5` | CYP3A5 activity | $0$ | \- | ✓ |
| $\mathrm{LI}_{\mathrm{\_RXEX\_k}}$ | `LI__RXEX_k` | rate rapamycin metabolites export | $100$ | 1/min | ✓ |
| $\mathrm{GU}_{\mathrm{\_f\_entero}}$ | `GU__f_entero` | fraction of the enterocytes of the intestine volume | $0.1$ | \- | ✓ |
| $\mathrm{GU}_{\mathrm{\_Vchain}}$ | `GU__Vchain` | volume of chain compartment | $0.1$ | l | ✓ |
| $\mathrm{GU}_{\mathrm{\_F\_rap\_abs}}$ | `GU__F_rap_abs` | fraction absorbed rapamycin (fasted) | $0.5$ | \- | ✓ |
| $\mathrm{GU}_{\mathrm{\_RAPABS\_k}}$ | `GU__RAPABS_k` | rate of rapamycin absorption | $0.05$ | 1/min | ✓ |
| $\mathrm{GU}_{\mathrm{\_fed}}$ | `GU__fed` | fed state (0 fasted, 1 fed) | $0$ | \- | ✓ |
| $\mathrm{GU}_{\mathrm{\_FED\_f\_absorption}}$ | `GU__FED_f_absorption` | fed effect on the fraction absorbed | $1.20035322885454$ | \- | ✓ |
| $\mathrm{GU}_{\mathrm{\_FED\_f\_dissolution}}$ | `GU__FED_f_dissolution` | fed effect on the dissolution and gastric emptying | $0.171348443283468$ | \- | ✓ |
| $\mathrm{GU}_{\mathrm{\_f\_absorption}}$ | `GU__f_absorption` | factor of the fraction absorbed |  | \- |  |
| $\mathrm{GU}_{\mathrm{\_f\_dissolution}}$ | `GU__f_dissolution` | factor of the dissolution and gastric emptying |  | \- |  |
| $\mathrm{GU}_{\mathrm{\_absorption}}$ | `GU__absorption` | rapamycin leaving the lumen |  | mmol/min |  |
| $\mathrm{GU}_{\mathrm{\_RAPBASO\_k}}$ | `GU__RAPBASO_k` | rate of rapamycin transport from the enterocytes into the blood | $0.05$ | 1/min | ✓ |
| $\mathrm{GU}_{\mathrm{\_RAP2RX\_k}}$ | `GU__RAP2RX_k` | rate of rapamycin metabolism by CYP3A4 (gut wall) | $0.1$ | 1/min | ✓ |
| $\mathrm{GU}_{\mathrm{\_f\_cyp3a4}}$ | `GU__f_cyp3a4` | CYP3A4 activity | $1$ | \- | ✓ |
| $\mathrm{GU}_{\mathrm{\_f\_cyp3a5}}$ | `GU__f_cyp3a5` | CYP3A5 activity | $0$ | \- | ✓ |
| $\mathrm{GU}_{\mathrm{\_RAPPG\_k}}$ | `GU__RAPPG_k` | rate of rapamycin efflux by P-glycoprotein | $0.02$ | 1/min | ✓ |
| $\mathrm{GU}_{\mathrm{\_f\_pg}}$ | `GU__f_pg` | P-glycoprotein activity | $1$ | \- | ✓ |
| $\mathrm{GU}_{\mathrm{\_RXEFF\_k}}$ | `GU__RXEFF_k` | rate of the release of the metabolites into the lumen | $1$ | 1/min | ✓ |
| $\mathrm{GU}_{\mathrm{\_RXEXC\_k}}$ | `GU__RXEXC_k` | rate of the transit of the metabolites to the feces | $0.0033$ | 1/min | ✓ |
| $\mathrm{GU}_{\mathrm{\_Ka\_dis\_rap}}$ | `GU__Ka_dis_rap` | Ka\_dis \[1/hr\] dissolution rapamycin | $3.64645771834172$ | 1/hr | ✓ |

## Initial assignments and assignment rules

The initial assignments set the values at $t = 0$:

$$
\begin{aligned}
\mathrm{Afeces}_{\mathrm{rap}} &= 0 \cdot \mathrm{Vfeces} \\[1ex]
\mathrm{Aurine}_{\mathrm{rx}} &= 0 \cdot \mathrm{Vurine} \\[1ex]
\mathrm{Afeces}_{\mathrm{rx}} &= 0 \cdot \mathrm{Vfeces} \\[1ex]
\mathrm{GU}_{\mathrm{\_rap\_stomach}} &= 0 \cdot \mathrm{Vstomach}
\end{aligned}
$$

The assignment rules hold at every time $t$:

$$
\begin{aligned}
f_{\mathrm{shunts}} &= f_{\mathrm{cirrhosis}} \\[1ex]
f_{\mathrm{tissue\_loss}} &= f_{\mathrm{cirrhosis}} \\[1ex]
\mathrm{FVre} &= 1 - \mathrm{FVgu} - \mathrm{FVki} - \mathrm{FVli} \\
&\quad {} - \mathrm{FVlu} - \mathrm{FVve} - \mathrm{FVar} \\[1ex]
\mathrm{FQre} &= 1 - \mathrm{FQki} - \mathrm{FQh} \\[1ex]
\mathrm{BSA} &= 0.024265 \cdot \mathopen{}\left(\frac{\mathrm{BW}}{1}\right)^{0.5378} \cdot \mathopen{}\left(\frac{\mathrm{HEIGHT}}{1}\right)^{0.3964} \\[1ex]
\mathrm{CO} &= f_{\mathrm{cardiac\_function}} \cdot \mathrm{BW} \cdot \mathrm{COBW} + \frac{\mathopen{}\left(\mathrm{HR} - \mathrm{HRrest}\right) \cdot \mathrm{COHRI}}{60} \\[1ex]
\mathrm{QC} &= \frac{\mathrm{CO}}{1000} \cdot 60 \\[1ex]
\mathrm{Vgu} &= \mathrm{BW} \cdot \mathrm{FVgu} \\[1ex]
\mathrm{Vgulumen} &= f_{\mathrm{lumen}} \cdot \mathrm{Vgu} \\[1ex]
\mathrm{Vki} &= \mathrm{BW} \cdot \mathrm{FVki} \\[1ex]
\mathrm{Vli} &= \mathrm{BW} \cdot \mathrm{FVli} \\[1ex]
\mathrm{Vlu} &= \mathrm{BW} \cdot \mathrm{FVlu} \\[1ex]
\mathrm{Vre} &= \mathrm{BW} \cdot \mathrm{FVre} \\[1ex]
\mathrm{Vve} &= \mathrm{BW} \cdot \mathrm{FVve} - \frac{\mathrm{FVve}}{\mathrm{FVar} + \mathrm{FVve}} \cdot \mathrm{BW} \cdot \mathrm{Fblood} \cdot \mathopen{}\left(1 - \mathrm{FVve} - \mathrm{FVar}\right) \\[1ex]
\mathrm{Var} &= \mathrm{BW} \cdot \mathrm{FVar} - \frac{\mathrm{FVar}}{\mathrm{FVar} + \mathrm{FVve}} \cdot \mathrm{BW} \cdot \mathrm{Fblood} \cdot \mathopen{}\left(1 - \mathrm{FVve} - \mathrm{FVar}\right) \\[1ex]
\mathrm{Vpo} &= \mathrm{BW} \cdot \mathrm{FVpo} - \frac{\mathrm{FVpo}}{\mathrm{FVar} + \mathrm{FVve} + \mathrm{FVpo} + \mathrm{FVhv}} \cdot \mathrm{BW} \cdot \mathrm{Fblood} \cdot \mathopen{}\left(1 - \mathrm{FVar} - \mathrm{FVve} - \mathrm{FVpo} - \mathrm{FVhv}\right) \\[1ex]
\mathrm{Vhv} &= \mathrm{BW} \cdot \mathrm{FVhv} - \frac{\mathrm{FVhv}}{\mathrm{FVar} + \mathrm{FVve} + \mathrm{FVpo} + \mathrm{FVhv}} \cdot \mathrm{BW} \cdot \mathrm{Fblood} \cdot \mathopen{}\left(1 - \mathrm{FVar} - \mathrm{FVve} - \mathrm{FVpo} - \mathrm{FVhv}\right) \\[1ex]
\mathrm{Qgu} &= \mathrm{QC} \cdot \mathrm{FQgu} \\[1ex]
\mathrm{Qki} &= \mathrm{QC} \cdot \mathrm{FQki} \\[1ex]
\mathrm{Qh} &= \mathrm{QC} \cdot \mathrm{FQh} \\[1ex]
\mathrm{Qha} &= \mathrm{Qh} - \mathrm{Qgu} \\[1ex]
\mathrm{Qlu} &= \mathrm{QC} \cdot \mathrm{FQlu} \\[1ex]
\mathrm{Qre} &= \mathrm{QC} \cdot \mathrm{FQre} \\[1ex]
\mathrm{Qpo} &= \mathrm{Qgu} \\[1ex]
\mathrm{Vki}_{\mathrm{blood}} &= \mathrm{Vki} \cdot \mathrm{Fblood} \\[1ex]
\mathrm{Vki}_{\mathrm{tissue}} &= \mathrm{Vki} \cdot \mathopen{}\left(1 - \mathrm{Fblood}\right) \\[1ex]
\mathrm{Vli}_{\mathrm{blood}} &= \mathrm{Vli} \cdot \mathrm{Fblood} \\[1ex]
\mathrm{Vli}_{\mathrm{tissue}} &= \mathrm{Vli} \cdot \mathopen{}\left(1 - f_{\mathrm{tissue\_loss}}\right) \cdot \mathopen{}\left(1 - \mathrm{Fblood}\right) \\[1ex]
\mathrm{Vlu}_{\mathrm{blood}} &= \mathrm{Vlu} \cdot \mathrm{Fblood} \\[1ex]
\mathrm{Vlu}_{\mathrm{tissue}} &= \mathrm{Vlu} \cdot \mathopen{}\left(1 - \mathrm{Fblood}\right) \\[1ex]
\mathrm{Vgu}_{\mathrm{blood}} &= \mathrm{Vgu} \cdot \mathrm{Fblood} \\[1ex]
\mathrm{Vgu}_{\mathrm{tissue}} &= \mathrm{Vgu} \cdot \mathopen{}\left(1 - \mathrm{Fblood}\right) \\[1ex]
\mathrm{Vre}_{\mathrm{blood}} &= \mathrm{Vre} \cdot \mathrm{Fblood} \\[1ex]
\mathrm{Vre}_{\mathrm{tissue}} &= \mathrm{Vre} \cdot \mathopen{}\left(1 - \mathrm{Fblood}\right) \\[1ex]
\mathrm{Ki}_{\mathrm{rap}} &= \frac{0.693}{\mathrm{ti}_{\mathrm{rap}}} \cdot 60 \\[1ex]
\mathrm{Cveblood}_{\mathrm{rap}} &= \mathrm{Cve}_{\mathrm{rap}} \\[1ex]
\mathrm{Cveplasma}_{\mathrm{rap}} &= \frac{\mathrm{Cve}_{\mathrm{rap}}}{\mathrm{BP}_{\mathrm{rap}}} \\[1ex]
\mathrm{Cveblood}_{\mathrm{raptot}} &= \mathrm{Cve}_{\mathrm{rap}} + \mathrm{Cve}_{\mathrm{rx}} \\[1ex]
\mathrm{Cveplasma}_{\mathrm{raptot}} &= \frac{\mathrm{Cve}_{\mathrm{rap}}}{\mathrm{BP}_{\mathrm{rap}}} + \mathrm{Cve}_{\mathrm{rx}} \\[1ex]
\mathrm{Afeces}_{\mathrm{raptot}} &= \mathrm{Afeces}_{\mathrm{rap}} + \mathrm{Afeces}_{\mathrm{rx}} \\[1ex]
\mathrm{GU}_{\mathrm{\_Vint\_0}} &= \mathrm{GU}_{\mathrm{\_Vchain}} \\[1ex]
\mathrm{GU}_{\mathrm{\_Vint\_1}} &= \mathrm{GU}_{\mathrm{\_Vchain}} \\[1ex]
\mathrm{GU}_{\mathrm{\_Vint\_2}} &= \mathrm{GU}_{\mathrm{\_Vchain}} \\[1ex]
\mathrm{GU}_{\mathrm{\_Vint\_3}} &= \mathrm{GU}_{\mathrm{\_Vchain}} \\[1ex]
\mathrm{GU}_{\mathrm{\_Vint\_4}} &= \mathrm{GU}_{\mathrm{\_Vchain}} \\[1ex]
\mathrm{GU}_{\mathrm{\_Ventero}} &= \mathrm{GU}_{\mathrm{\_f\_entero}} \cdot \mathrm{Vgu} \\[1ex]
\mathrm{GU}_{\mathrm{\_f\_absorption}} &= 1 + \mathrm{GU}_{\mathrm{\_fed}} \cdot \mathopen{}\left(\mathrm{GU}_{\mathrm{\_FED\_f\_absorption}} - 1\right) \\[1ex]
\mathrm{GU}_{\mathrm{\_f\_dissolution}} &= 1 + \mathrm{GU}_{\mathrm{\_fed}} \cdot \mathopen{}\left(\mathrm{GU}_{\mathrm{\_FED\_f\_dissolution}} - 1\right) \\[1ex]
\mathrm{GU}_{\mathrm{\_absorption}} &= \mathrm{GU}_{\mathrm{\_RAPABS\_k}} \cdot \mathrm{Vgulumen} \cdot \mathrm{Clumen}_{\mathrm{rap}}
\end{aligned}
$$

## Reactions

| Rate | Id | Name | Equation |
| --- | --- | --- | --- |
| $v_{\mathrm{transport\_lu\_rap}}$ | `transport_lu_rap` | transport rapamycin | $\mathrm{Clu}_{\mathrm{blood\_rap}} \rightleftharpoons \mathrm{Clu}_{\mathrm{rap}}$ |
| $v_{\mathrm{transport\_re\_rap}}$ | `transport_re_rap` | transport rapamycin | $\mathrm{Cre}_{\mathrm{blood\_rap}} \rightleftharpoons \mathrm{Cre}_{\mathrm{rap}}$ |
| $v_{\mathrm{deep\_re\_rap}}$ | `deep_re_rap` | slow binding rapamycin (tissue) | $\mathrm{Cre}_{\mathrm{rap}} \rightleftharpoons \mathrm{Cre}_{\mathrm{deep\_rap}}$ |
| $v_{\mathrm{iv\_rap}}$ | `iv_rap` | iv rapamycin | $\varnothing \longrightarrow \mathrm{Cve}_{\mathrm{rap}}$ |
| $v_{\mathrm{Flow\_ar\_ki\_rap}}$ | `Flow_ar_ki_rap` | inflow kidney rapamycin | $\mathrm{Car}_{\mathrm{rap}} \longrightarrow \mathrm{Cki}_{\mathrm{blood\_rap}}$ |
| $v_{\mathrm{Flow\_ki\_ve\_rap}}$ | `Flow_ki_ve_rap` | outflow kidney rapamycin | $\mathrm{Cki}_{\mathrm{blood\_rap}} \longrightarrow \mathrm{Cve}_{\mathrm{rap}}$ |
| $v_{\mathrm{Flow\_arli\_li\_rap}}$ | `Flow_arli_li_rap` | arterial inflow liver rapamycin | $\mathrm{Car}_{\mathrm{rap}} \longrightarrow \mathrm{Cli}_{\mathrm{blood\_rap}}$ |
| $v_{\mathrm{Flow\_arli\_hv\_rap}}$ | `Flow_arli_hv_rap` | flow arterial shunts | $\mathrm{Car}_{\mathrm{rap}} \longrightarrow \mathrm{Chv}_{\mathrm{rap}}$ |
| $v_{\mathrm{Flow\_po\_li\_rap}}$ | `Flow_po_li_rap` | outflow po rapamycin | $\mathrm{Cpo}_{\mathrm{rap}} \longrightarrow \mathrm{Cli}_{\mathrm{blood\_rap}}$ |
| $v_{\mathrm{Flow\_po\_hv\_rap}}$ | `Flow_po_hv_rap` | flow portal shunts | $\mathrm{Cpo}_{\mathrm{rap}} \longrightarrow \mathrm{Chv}_{\mathrm{rap}}$ |
| $v_{\mathrm{Flow\_li\_hv\_rap}}$ | `Flow_li_hv_rap` | outflow liver rapamycin | $\mathrm{Cli}_{\mathrm{blood\_rap}} \longrightarrow \mathrm{Chv}_{\mathrm{rap}}$ |
| $v_{\mathrm{Flow\_hv\_ve\_rap}}$ | `Flow_hv_ve_rap` | outflow hepatic vein rapamycin | $\mathrm{Chv}_{\mathrm{rap}} \longrightarrow \mathrm{Cve}_{\mathrm{rap}}$ |
| $v_{\mathrm{Flow\_ve\_lu\_rap}}$ | `Flow_ve_lu_rap` | inflow lung rapamycin | $\mathrm{Cve}_{\mathrm{rap}} \longrightarrow \mathrm{Clu}_{\mathrm{blood\_rap}}$ |
| $v_{\mathrm{Flow\_lu\_ar\_rap}}$ | `Flow_lu_ar_rap` | outflow lung rapamycin | $\mathrm{Clu}_{\mathrm{blood\_rap}} \longrightarrow \mathrm{Car}_{\mathrm{rap}}$ |
| $v_{\mathrm{Flow\_ar\_gu\_rap}}$ | `Flow_ar_gu_rap` | inflow gut rapamycin | $\mathrm{Car}_{\mathrm{rap}} \longrightarrow \mathrm{Cgu}_{\mathrm{blood\_rap}}$ |
| $v_{\mathrm{Flow\_gu\_po\_rap}}$ | `Flow_gu_po_rap` | outflow gut rapamycin | $\mathrm{Cgu}_{\mathrm{blood\_rap}} \longrightarrow \mathrm{Cpo}_{\mathrm{rap}}$ |
| $v_{\mathrm{Flow\_ar\_re\_rap}}$ | `Flow_ar_re_rap` | inflow rest rapamycin | $\mathrm{Car}_{\mathrm{rap}} \longrightarrow \mathrm{Cre}_{\mathrm{blood\_rap}}$ |
| $v_{\mathrm{Flow\_re\_ve\_rap}}$ | `Flow_re_ve_rap` | outflow rest rapamycin | $\mathrm{Cre}_{\mathrm{blood\_rap}} \longrightarrow \mathrm{Cve}_{\mathrm{rap}}$ |
| $v_{\mathrm{Flow\_ar\_ki\_rx}}$ | `Flow_ar_ki_rx` | inflow kidney rapamycin metabolites | $\mathrm{Car}_{\mathrm{rx}} \longrightarrow \mathrm{Cki}_{\mathrm{blood\_rx}}$ |
| $v_{\mathrm{Flow\_ki\_ve\_rx}}$ | `Flow_ki_ve_rx` | outflow kidney rapamycin metabolites | $\mathrm{Cki}_{\mathrm{blood\_rx}} \longrightarrow \mathrm{Cve}_{\mathrm{rx}}$ |
| $v_{\mathrm{Flow\_arli\_li\_rx}}$ | `Flow_arli_li_rx` | arterial inflow liver rapamycin metabolites | $\mathrm{Car}_{\mathrm{rx}} \longrightarrow \mathrm{Cli}_{\mathrm{blood\_rx}}$ |
| $v_{\mathrm{Flow\_arli\_hv\_rx}}$ | `Flow_arli_hv_rx` | flow arterial shunts | $\mathrm{Car}_{\mathrm{rx}} \longrightarrow \mathrm{Chv}_{\mathrm{rx}}$ |
| $v_{\mathrm{Flow\_po\_li\_rx}}$ | `Flow_po_li_rx` | outflow po rapamycin metabolites | $\mathrm{Cpo}_{\mathrm{rx}} \longrightarrow \mathrm{Cli}_{\mathrm{blood\_rx}}$ |
| $v_{\mathrm{Flow\_po\_hv\_rx}}$ | `Flow_po_hv_rx` | flow portal shunts | $\mathrm{Cpo}_{\mathrm{rx}} \longrightarrow \mathrm{Chv}_{\mathrm{rx}}$ |
| $v_{\mathrm{Flow\_li\_hv\_rx}}$ | `Flow_li_hv_rx` | outflow liver rapamycin metabolites | $\mathrm{Cli}_{\mathrm{blood\_rx}} \longrightarrow \mathrm{Chv}_{\mathrm{rx}}$ |
| $v_{\mathrm{Flow\_hv\_ve\_rx}}$ | `Flow_hv_ve_rx` | outflow hepatic vein rapamycin metabolites | $\mathrm{Chv}_{\mathrm{rx}} \longrightarrow \mathrm{Cve}_{\mathrm{rx}}$ |
| $v_{\mathrm{Flow\_ve\_lu\_rx}}$ | `Flow_ve_lu_rx` | inflow lung rapamycin metabolites | $\mathrm{Cve}_{\mathrm{rx}} \longrightarrow \mathrm{Clu}_{\mathrm{blood\_rx}}$ |
| $v_{\mathrm{Flow\_lu\_ar\_rx}}$ | `Flow_lu_ar_rx` | outflow lung rapamycin metabolites | $\mathrm{Clu}_{\mathrm{blood\_rx}} \longrightarrow \mathrm{Car}_{\mathrm{rx}}$ |
| $v_{\mathrm{Flow\_ar\_gu\_rx}}$ | `Flow_ar_gu_rx` | inflow gut rapamycin metabolites | $\mathrm{Car}_{\mathrm{rx}} \longrightarrow \mathrm{Cgu}_{\mathrm{blood\_rx}}$ |
| $v_{\mathrm{Flow\_gu\_po\_rx}}$ | `Flow_gu_po_rx` | outflow gut rapamycin metabolites | $\mathrm{Cgu}_{\mathrm{blood\_rx}} \longrightarrow \mathrm{Cpo}_{\mathrm{rx}}$ |
| $v_{\mathrm{Flow\_ar\_re\_rx}}$ | `Flow_ar_re_rx` | inflow rest rapamycin metabolites | $\mathrm{Car}_{\mathrm{rx}} \longrightarrow \mathrm{Cre}_{\mathrm{blood\_rx}}$ |
| $v_{\mathrm{Flow\_re\_ve\_rx}}$ | `Flow_re_ve_rx` | outflow rest rapamycin metabolites | $\mathrm{Cre}_{\mathrm{blood\_rx}} \longrightarrow \mathrm{Cve}_{\mathrm{rx}}$ |
| $v_{\mathrm{KI\_\_RXEX}}$ | `KI__RXEX` | rapamycin metabolites excretion (RXEX) | $\mathrm{Cki}_{\mathrm{blood\_rx}} \longrightarrow \mathrm{Aurine}_{\mathrm{rx}}$ |
| $v_{\mathrm{LI\_\_RAPIM}}$ | `LI__RAPIM` | rapamycin import (RAPIM) | $\mathrm{Cli}_{\mathrm{blood\_rap}} \rightleftharpoons \mathrm{LI}_{\mathrm{\_rap}}$ |
| $v_{\mathrm{LI\_\_RAP2RX}}$ | `LI__RAP2RX` | rapamycin metabolism (CYP3A4/5) | $\mathrm{LI}_{\mathrm{\_rap}} \longrightarrow \mathrm{LI}_{\mathrm{\_rx}}$ |
| $v_{\mathrm{LI\_\_RXEX}}$ | `LI__RXEX` | rapamycin metabolites export (RXEX) | $\mathrm{LI}_{\mathrm{\_rx}} \rightleftharpoons \mathrm{Cli}_{\mathrm{blood\_rx}}$ |
| $v_{\mathrm{LI\_\_RXBEX}}$ | `LI__RXBEX` | rapamycin metabolites biliary excretion | $\mathrm{LI}_{\mathrm{\_rx}} \longrightarrow \mathrm{Clumen}_{\mathrm{rx}}$ |
| $v_{\mathrm{GU\_\_RAPABS}}$ | `GU__RAPABS` | absorption rapamycin (enterocytes) | $\mathrm{Clumen}_{\mathrm{rap}} \longrightarrow \mathrm{GU}_{\mathrm{\_rap\_entero}}$ |
| $v_{\mathrm{GU\_\_RAPEXC}}$ | `GU__RAPEXC` | excretion rapamycin (colon) | $\mathrm{Clumen}_{\mathrm{rap}} \longrightarrow \mathrm{GU}_{\mathrm{\_rap\_colon}}$ |
| $v_{\mathrm{GU\_\_RAPFEC}}$ | `GU__RAPFEC` | excretion rapamycin (feces) | $\mathrm{GU}_{\mathrm{\_rap\_colon}} \longrightarrow \mathrm{Afeces}_{\mathrm{rap}}$ |
| $v_{\mathrm{GU\_\_RAPBASO}}$ | `GU__RAPBASO` | transport rapamycin (blood) | $\mathrm{GU}_{\mathrm{\_rap\_entero}} \longrightarrow \mathrm{Cgu}_{\mathrm{blood\_rap}}$ |
| $v_{\mathrm{GU\_\_RAP2RX}}$ | `GU__RAP2RX` | rapamycin metabolism via CYP3A4/5 | $\mathrm{GU}_{\mathrm{\_rap\_entero}} \longrightarrow \mathrm{GU}_{\mathrm{\_rx\_entero}}$ |
| $v_{\mathrm{GU\_\_RAPPG}}$ | `GU__RAPPG` | efflux rapamycin by P-glycoprotein | $\mathrm{GU}_{\mathrm{\_rap\_entero}} \longrightarrow \mathrm{Clumen}_{\mathrm{rap}}$ |
| $v_{\mathrm{GU\_\_RXEFF}}$ | `GU__RXEFF` | release rapamycin metabolites (lumen) | $\mathrm{GU}_{\mathrm{\_rx\_entero}} \longrightarrow \mathrm{Clumen}_{\mathrm{rx}}$ |
| $v_{\mathrm{GU\_\_RXEXC}}$ | `GU__RXEXC` | excretion rapamycin metabolites (feces) | $\mathrm{Clumen}_{\mathrm{rx}} \longrightarrow \mathrm{GU}_{\mathrm{\_rx\_int\_0}}$ |
| $v_{\mathrm{GU\_\_RXEXC\_0}}$ | `GU__RXEXC_0` | excretion rapamycin metabolites 0 | $\mathrm{GU}_{\mathrm{\_rx\_int\_0}} \longrightarrow \mathrm{GU}_{\mathrm{\_rx\_int\_1}}$ |
| $v_{\mathrm{GU\_\_RXEXC\_1}}$ | `GU__RXEXC_1` | excretion rapamycin metabolites 1 | $\mathrm{GU}_{\mathrm{\_rx\_int\_1}} \longrightarrow \mathrm{GU}_{\mathrm{\_rx\_int\_2}}$ |
| $v_{\mathrm{GU\_\_RXEXC\_2}}$ | `GU__RXEXC_2` | excretion rapamycin metabolites 2 | $\mathrm{GU}_{\mathrm{\_rx\_int\_2}} \longrightarrow \mathrm{GU}_{\mathrm{\_rx\_int\_3}}$ |
| $v_{\mathrm{GU\_\_RXEXC\_3}}$ | `GU__RXEXC_3` | excretion rapamycin metabolites 3 | $\mathrm{GU}_{\mathrm{\_rx\_int\_3}} \longrightarrow \mathrm{GU}_{\mathrm{\_rx\_int\_4}}$ |
| $v_{\mathrm{GU\_\_RXEXC\_4}}$ | `GU__RXEXC_4` | excretion rapamycin metabolites 4 | $\mathrm{GU}_{\mathrm{\_rx\_int\_4}} \longrightarrow \mathrm{Afeces}_{\mathrm{rx}}$ |
| $v_{\mathrm{GU\_\_dissolution\_rap}}$ | `GU__dissolution_rap` | dissolution rapamycin | $\mathrm{GU}_{\mathrm{\_rap\_stomach}} \longrightarrow \mathrm{Clumen}_{\mathrm{rap}}$ |

The rates of the reactions are:

$$
\begin{aligned}
v_{\mathrm{transport\_lu\_rap}} &= \frac{\mathrm{ftissue}_{\mathrm{rap}} \cdot \mathrm{BW}}{\mathrm{BW}_{\mathrm{ref}}} \cdot \mathopen{}\left(\mathrm{Clu}_{\mathrm{blood\_rap}} \cdot \mathrm{Kp}_{\mathrm{rap}} - \mathrm{Clu}_{\mathrm{rap}}\right) \\[1ex]
v_{\mathrm{transport\_re\_rap}} &= \frac{\mathrm{ftissue}_{\mathrm{rap}} \cdot \mathrm{BW}}{\mathrm{BW}_{\mathrm{ref}}} \cdot \mathopen{}\left(\mathrm{Cre}_{\mathrm{blood\_rap}} \cdot \mathrm{Kp}_{\mathrm{rap}} - \mathrm{Cre}_{\mathrm{rap}}\right) \\[1ex]
v_{\mathrm{deep\_re\_rap}} &= \mathrm{kdeep}_{\mathrm{rap}} \cdot \mathrm{Vre}_{\mathrm{tissue}} \cdot \mathopen{}\left(\mathrm{Kpdeep}_{\mathrm{rap}} \cdot \mathrm{Cre}_{\mathrm{rap}} - \mathrm{Cre}_{\mathrm{deep\_rap}}\right) \\[1ex]
v_{\mathrm{iv\_rap}} &= \frac{\mathrm{Ki}_{\mathrm{rap}} \cdot \mathrm{IVDOSE}_{\mathrm{rap}}}{\mathrm{Mr}_{\mathrm{rap}}} \\[1ex]
v_{\mathrm{Flow\_ar\_ki\_rap}} &= \mathrm{Qki} \cdot \mathrm{Car}_{\mathrm{rap}} \\[1ex]
v_{\mathrm{Flow\_ki\_ve\_rap}} &= \mathrm{Qki} \cdot \mathrm{Cki}_{\mathrm{blood\_rap}} \\[1ex]
v_{\mathrm{Flow\_arli\_li\_rap}} &= \mathopen{}\left(1 - f_{\mathrm{shunts}}\right) \cdot \mathrm{Qha} \cdot \mathrm{Car}_{\mathrm{rap}} \\[1ex]
v_{\mathrm{Flow\_arli\_hv\_rap}} &= f_{\mathrm{shunts}} \cdot \mathrm{Qha} \cdot \mathrm{Car}_{\mathrm{rap}} \\[1ex]
v_{\mathrm{Flow\_po\_li\_rap}} &= \mathopen{}\left(1 - f_{\mathrm{shunts}}\right) \cdot \mathrm{Qpo} \cdot \mathrm{Cpo}_{\mathrm{rap}} \\[1ex]
v_{\mathrm{Flow\_po\_hv\_rap}} &= f_{\mathrm{shunts}} \cdot \mathrm{Qpo} \cdot \mathrm{Cpo}_{\mathrm{rap}} \\[1ex]
v_{\mathrm{Flow\_li\_hv\_rap}} &= \mathopen{}\left(1 - f_{\mathrm{shunts}}\right) \cdot \mathopen{}\left(\mathrm{Qpo} + \mathrm{Qha}\right) \cdot \mathrm{Cli}_{\mathrm{blood\_rap}} \\[1ex]
v_{\mathrm{Flow\_hv\_ve\_rap}} &= \mathrm{Qh} \cdot \mathrm{Chv}_{\mathrm{rap}} \\[1ex]
v_{\mathrm{Flow\_ve\_lu\_rap}} &= \mathrm{Qlu} \cdot \mathrm{Cve}_{\mathrm{rap}} \\[1ex]
v_{\mathrm{Flow\_lu\_ar\_rap}} &= \mathrm{Qlu} \cdot \mathrm{Clu}_{\mathrm{blood\_rap}} \\[1ex]
v_{\mathrm{Flow\_ar\_gu\_rap}} &= \mathrm{Qgu} \cdot \mathrm{Car}_{\mathrm{rap}} \\[1ex]
v_{\mathrm{Flow\_gu\_po\_rap}} &= \mathrm{Qgu} \cdot \mathrm{Cgu}_{\mathrm{blood\_rap}} \\[1ex]
v_{\mathrm{Flow\_ar\_re\_rap}} &= \mathrm{Qre} \cdot \mathrm{Car}_{\mathrm{rap}} \\[1ex]
v_{\mathrm{Flow\_re\_ve\_rap}} &= \mathrm{Qre} \cdot \mathrm{Cre}_{\mathrm{blood\_rap}} \\[1ex]
v_{\mathrm{Flow\_ar\_ki\_rx}} &= \mathrm{Qki} \cdot \mathrm{Car}_{\mathrm{rx}} \\[1ex]
v_{\mathrm{Flow\_ki\_ve\_rx}} &= \mathrm{Qki} \cdot \mathrm{Cki}_{\mathrm{blood\_rx}} \\[1ex]
v_{\mathrm{Flow\_arli\_li\_rx}} &= \mathopen{}\left(1 - f_{\mathrm{shunts}}\right) \cdot \mathrm{Qha} \cdot \mathrm{Car}_{\mathrm{rx}} \\[1ex]
v_{\mathrm{Flow\_arli\_hv\_rx}} &= f_{\mathrm{shunts}} \cdot \mathrm{Qha} \cdot \mathrm{Car}_{\mathrm{rx}} \\[1ex]
v_{\mathrm{Flow\_po\_li\_rx}} &= \mathopen{}\left(1 - f_{\mathrm{shunts}}\right) \cdot \mathrm{Qpo} \cdot \mathrm{Cpo}_{\mathrm{rx}} \\[1ex]
v_{\mathrm{Flow\_po\_hv\_rx}} &= f_{\mathrm{shunts}} \cdot \mathrm{Qpo} \cdot \mathrm{Cpo}_{\mathrm{rx}} \\[1ex]
v_{\mathrm{Flow\_li\_hv\_rx}} &= \mathopen{}\left(1 - f_{\mathrm{shunts}}\right) \cdot \mathopen{}\left(\mathrm{Qpo} + \mathrm{Qha}\right) \cdot \mathrm{Cli}_{\mathrm{blood\_rx}} \\[1ex]
v_{\mathrm{Flow\_hv\_ve\_rx}} &= \mathrm{Qh} \cdot \mathrm{Chv}_{\mathrm{rx}} \\[1ex]
v_{\mathrm{Flow\_ve\_lu\_rx}} &= \mathrm{Qlu} \cdot \mathrm{Cve}_{\mathrm{rx}} \\[1ex]
v_{\mathrm{Flow\_lu\_ar\_rx}} &= \mathrm{Qlu} \cdot \mathrm{Clu}_{\mathrm{blood\_rx}} \\[1ex]
v_{\mathrm{Flow\_ar\_gu\_rx}} &= \mathrm{Qgu} \cdot \mathrm{Car}_{\mathrm{rx}} \\[1ex]
v_{\mathrm{Flow\_gu\_po\_rx}} &= \mathrm{Qgu} \cdot \mathrm{Cgu}_{\mathrm{blood\_rx}} \\[1ex]
v_{\mathrm{Flow\_ar\_re\_rx}} &= \mathrm{Qre} \cdot \mathrm{Car}_{\mathrm{rx}} \\[1ex]
v_{\mathrm{Flow\_re\_ve\_rx}} &= \mathrm{Qre} \cdot \mathrm{Cre}_{\mathrm{blood\_rx}} \\[1ex]
v_{\mathrm{KI\_\_RXEX}} &= \mathrm{KI}_{\mathrm{\_f\_renal\_function}} \cdot \mathrm{KI}_{\mathrm{\_RXEX\_k}} \cdot \mathrm{Vki}_{\mathrm{tissue}} \cdot \mathrm{Cki}_{\mathrm{blood\_rx}} \\[1ex]
v_{\mathrm{LI\_\_RAPIM}} &= \mathrm{LI}_{\mathrm{\_RAPIM\_k}} \cdot \mathrm{Vli}_{\mathrm{tissue}} \cdot \mathopen{}\left(\mathrm{Cli}_{\mathrm{blood\_rap}} - \mathrm{LI}_{\mathrm{\_rap}}\right) \\[1ex]
v_{\mathrm{LI\_\_RAP2RX}} &= \mathopen{}\left(\mathrm{LI}_{\mathrm{\_f\_cyp3a4}} + r_{\mathrm{cyp3a5}} \cdot \mathrm{LI}_{\mathrm{\_f\_cyp3a5}}\right) \cdot \mathrm{LI}_{\mathrm{\_RAP2RX\_k}} \cdot \mathrm{Vli}_{\mathrm{tissue}} \cdot \mathrm{LI}_{\mathrm{\_rap}} \\[1ex]
v_{\mathrm{LI\_\_RXEX}} &= \mathrm{LI}_{\mathrm{\_RXEX\_k}} \cdot \mathrm{Vli}_{\mathrm{tissue}} \cdot \mathopen{}\left(\mathrm{LI}_{\mathrm{\_rx}} - \mathrm{Cli}_{\mathrm{blood\_rx}}\right) \\[1ex]
v_{\mathrm{LI\_\_RXBEX}} &= \mathrm{LI}_{\mathrm{\_RXBEX\_k}} \cdot \mathrm{Vli}_{\mathrm{tissue}} \cdot \mathrm{LI}_{\mathrm{\_rx}} \\[1ex]
v_{\mathrm{GU\_\_RAPABS}} &= \mathrm{GU}_{\mathrm{\_F\_rap\_abs}} \cdot \mathrm{GU}_{\mathrm{\_f\_absorption}} \cdot \mathrm{GU}_{\mathrm{\_absorption}} \\[1ex]
v_{\mathrm{GU\_\_RAPEXC}} &= \mathopen{}\left(1 - \mathrm{GU}_{\mathrm{\_F\_rap\_abs}} \cdot \mathrm{GU}_{\mathrm{\_f\_absorption}}\right) \cdot \mathrm{GU}_{\mathrm{\_absorption}} \\[1ex]
v_{\mathrm{GU\_\_RAPFEC}} &= \frac{\mathrm{GU}_{\mathrm{\_RXEXC\_k}}}{6} \cdot \mathrm{Vgulumen} \cdot \mathrm{GU}_{\mathrm{\_rap\_colon}} \\[1ex]
v_{\mathrm{GU\_\_RAPBASO}} &= \mathrm{GU}_{\mathrm{\_RAPBASO\_k}} \cdot \mathrm{GU}_{\mathrm{\_Ventero}} \cdot \mathrm{GU}_{\mathrm{\_rap\_entero}} \\[1ex]
v_{\mathrm{GU\_\_RAP2RX}} &= \mathopen{}\left(\mathrm{GU}_{\mathrm{\_f\_cyp3a4}} + r_{\mathrm{cyp3a5}} \cdot \mathrm{GU}_{\mathrm{\_f\_cyp3a5}}\right) \cdot \mathrm{GU}_{\mathrm{\_RAP2RX\_k}} \cdot \mathrm{GU}_{\mathrm{\_Ventero}} \cdot \mathrm{GU}_{\mathrm{\_rap\_entero}} \\[1ex]
v_{\mathrm{GU\_\_RAPPG}} &= \mathrm{GU}_{\mathrm{\_f\_pg}} \cdot \mathrm{GU}_{\mathrm{\_RAPPG\_k}} \cdot \mathrm{GU}_{\mathrm{\_Ventero}} \cdot \mathrm{GU}_{\mathrm{\_rap\_entero}} \\[1ex]
v_{\mathrm{GU\_\_RXEFF}} &= \mathrm{GU}_{\mathrm{\_RXEFF\_k}} \cdot \mathrm{GU}_{\mathrm{\_Ventero}} \cdot \mathrm{GU}_{\mathrm{\_rx\_entero}} \\[1ex]
v_{\mathrm{GU\_\_RXEXC}} &= \mathrm{GU}_{\mathrm{\_RXEXC\_k}} \cdot \mathrm{Vgulumen} \cdot \mathrm{Clumen}_{\mathrm{rx}} \\[1ex]
v_{\mathrm{GU\_\_RXEXC\_0}} &= \mathrm{GU}_{\mathrm{\_RXEXC\_k}} \cdot \mathrm{GU}_{\mathrm{\_Vint\_0}} \cdot \mathrm{GU}_{\mathrm{\_rx\_int\_0}} \\[1ex]
v_{\mathrm{GU\_\_RXEXC\_1}} &= \mathrm{GU}_{\mathrm{\_RXEXC\_k}} \cdot \mathrm{GU}_{\mathrm{\_Vint\_1}} \cdot \mathrm{GU}_{\mathrm{\_rx\_int\_1}} \\[1ex]
v_{\mathrm{GU\_\_RXEXC\_2}} &= \mathrm{GU}_{\mathrm{\_RXEXC\_k}} \cdot \mathrm{GU}_{\mathrm{\_Vint\_2}} \cdot \mathrm{GU}_{\mathrm{\_rx\_int\_2}} \\[1ex]
v_{\mathrm{GU\_\_RXEXC\_3}} &= \mathrm{GU}_{\mathrm{\_RXEXC\_k}} \cdot \mathrm{GU}_{\mathrm{\_Vint\_3}} \cdot \mathrm{GU}_{\mathrm{\_rx\_int\_3}} \\[1ex]
v_{\mathrm{GU\_\_RXEXC\_4}} &= \mathrm{GU}_{\mathrm{\_RXEXC\_k}} \cdot \mathrm{GU}_{\mathrm{\_Vint\_4}} \cdot \mathrm{GU}_{\mathrm{\_rx\_int\_4}} \\[1ex]
v_{\mathrm{GU\_\_dissolution\_rap}} &= \frac{\frac{\mathrm{GU}_{\mathrm{\_f\_dissolution}} \cdot \mathrm{GU}_{\mathrm{\_Ka\_dis\_rap}}}{60} \cdot \mathrm{PODOSE}_{\mathrm{rap}}}{\mathrm{Mr}_{\mathrm{rap}}}
\end{aligned}
$$

## ODE system

The states change in time with the rates of the reactions and the rate rules:

$$
\begin{aligned}
\frac{\mathrm{d} \,\mathrm{Cki}_{\mathrm{blood\_rap}}}{\mathrm{d} t} &= \frac{v_{\mathrm{Flow\_ar\_ki\_rap}} - v_{\mathrm{Flow\_ki\_ve\_rap}}}{\mathrm{Vki}_{\mathrm{blood}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Cli}_{\mathrm{blood\_rap}}}{\mathrm{d} t} &= \frac{v_{\mathrm{Flow\_arli\_li\_rap}} + v_{\mathrm{Flow\_po\_li\_rap}} - v_{\mathrm{Flow\_li\_hv\_rap}} - v_{\mathrm{LI\_\_RAPIM}}}{\mathrm{Vli}_{\mathrm{blood}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Clu}_{\mathrm{blood\_rap}}}{\mathrm{d} t} &= \frac{-v_{\mathrm{transport\_lu\_rap}} + v_{\mathrm{Flow\_ve\_lu\_rap}} - v_{\mathrm{Flow\_lu\_ar\_rap}}}{\mathrm{Vlu}_{\mathrm{blood}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Cgu}_{\mathrm{blood\_rap}}}{\mathrm{d} t} &= \frac{v_{\mathrm{Flow\_ar\_gu\_rap}} - v_{\mathrm{Flow\_gu\_po\_rap}} + v_{\mathrm{GU\_\_RAPBASO}}}{\mathrm{Vgu}_{\mathrm{blood}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Cre}_{\mathrm{blood\_rap}}}{\mathrm{d} t} &= \frac{-v_{\mathrm{transport\_re\_rap}} + v_{\mathrm{Flow\_ar\_re\_rap}} - v_{\mathrm{Flow\_re\_ve\_rap}}}{\mathrm{Vre}_{\mathrm{blood}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Car}_{\mathrm{rap}}}{\mathrm{d} t} &= \frac{1}{\mathrm{Var}} \Bigl( -v_{\mathrm{Flow\_ar\_ki\_rap}} - v_{\mathrm{Flow\_arli\_li\_rap}} - v_{\mathrm{Flow\_arli\_hv\_rap}} + v_{\mathrm{Flow\_lu\_ar\_rap}} \\
&\quad {} - v_{\mathrm{Flow\_ar\_gu\_rap}} - v_{\mathrm{Flow\_ar\_re\_rap}} \Bigr) \\[1ex]
\frac{\mathrm{d} \,\mathrm{Cve}_{\mathrm{rap}}}{\mathrm{d} t} &= \frac{1}{\mathrm{Vve}} \Bigl( v_{\mathrm{iv\_rap}} + v_{\mathrm{Flow\_ki\_ve\_rap}} + v_{\mathrm{Flow\_hv\_ve\_rap}} - v_{\mathrm{Flow\_ve\_lu\_rap}} \\
&\quad {} + v_{\mathrm{Flow\_re\_ve\_rap}} \Bigr) \\[1ex]
\frac{\mathrm{d} \,\mathrm{Cpo}_{\mathrm{rap}}}{\mathrm{d} t} &= \frac{-v_{\mathrm{Flow\_po\_li\_rap}} - v_{\mathrm{Flow\_po\_hv\_rap}} + v_{\mathrm{Flow\_gu\_po\_rap}}}{\mathrm{Vpo}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Chv}_{\mathrm{rap}}}{\mathrm{d} t} &= \frac{v_{\mathrm{Flow\_arli\_hv\_rap}} + v_{\mathrm{Flow\_po\_hv\_rap}} + v_{\mathrm{Flow\_li\_hv\_rap}} - v_{\mathrm{Flow\_hv\_ve\_rap}}}{\mathrm{Vhv}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Clu}_{\mathrm{rap}}}{\mathrm{d} t} &= \frac{v_{\mathrm{transport\_lu\_rap}}}{\mathrm{Vlu}_{\mathrm{tissue}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Cre}_{\mathrm{rap}}}{\mathrm{d} t} &= \frac{v_{\mathrm{transport\_re\_rap}} - v_{\mathrm{deep\_re\_rap}}}{\mathrm{Vre}_{\mathrm{tissue}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Cre}_{\mathrm{deep\_rap}}}{\mathrm{d} t} &= \frac{v_{\mathrm{deep\_re\_rap}}}{\mathrm{Vre}_{\mathrm{tissue}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Afeces}_{\mathrm{rap}}}{\mathrm{d} t} &= v_{\mathrm{GU\_\_RAPFEC}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Clumen}_{\mathrm{rap}}}{\mathrm{d} t} &= \frac{-v_{\mathrm{GU\_\_RAPABS}} - v_{\mathrm{GU\_\_RAPEXC}} + v_{\mathrm{GU\_\_RAPPG}} + v_{\mathrm{GU\_\_dissolution\_rap}}}{\mathrm{Vgulumen}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Cki}_{\mathrm{blood\_rx}}}{\mathrm{d} t} &= \frac{v_{\mathrm{Flow\_ar\_ki\_rx}} - v_{\mathrm{Flow\_ki\_ve\_rx}} - v_{\mathrm{KI\_\_RXEX}}}{\mathrm{Vki}_{\mathrm{blood}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Cli}_{\mathrm{blood\_rx}}}{\mathrm{d} t} &= \frac{v_{\mathrm{Flow\_arli\_li\_rx}} + v_{\mathrm{Flow\_po\_li\_rx}} - v_{\mathrm{Flow\_li\_hv\_rx}} + v_{\mathrm{LI\_\_RXEX}}}{\mathrm{Vli}_{\mathrm{blood}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Clu}_{\mathrm{blood\_rx}}}{\mathrm{d} t} &= \frac{v_{\mathrm{Flow\_ve\_lu\_rx}} - v_{\mathrm{Flow\_lu\_ar\_rx}}}{\mathrm{Vlu}_{\mathrm{blood}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Cgu}_{\mathrm{blood\_rx}}}{\mathrm{d} t} &= \frac{v_{\mathrm{Flow\_ar\_gu\_rx}} - v_{\mathrm{Flow\_gu\_po\_rx}}}{\mathrm{Vgu}_{\mathrm{blood}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Cre}_{\mathrm{blood\_rx}}}{\mathrm{d} t} &= \frac{v_{\mathrm{Flow\_ar\_re\_rx}} - v_{\mathrm{Flow\_re\_ve\_rx}}}{\mathrm{Vre}_{\mathrm{blood}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Car}_{\mathrm{rx}}}{\mathrm{d} t} &= \frac{1}{\mathrm{Var}} \Bigl( -v_{\mathrm{Flow\_ar\_ki\_rx}} - v_{\mathrm{Flow\_arli\_li\_rx}} - v_{\mathrm{Flow\_arli\_hv\_rx}} + v_{\mathrm{Flow\_lu\_ar\_rx}} \\
&\quad {} - v_{\mathrm{Flow\_ar\_gu\_rx}} - v_{\mathrm{Flow\_ar\_re\_rx}} \Bigr) \\[1ex]
\frac{\mathrm{d} \,\mathrm{Cve}_{\mathrm{rx}}}{\mathrm{d} t} &= \frac{v_{\mathrm{Flow\_ki\_ve\_rx}} + v_{\mathrm{Flow\_hv\_ve\_rx}} - v_{\mathrm{Flow\_ve\_lu\_rx}} + v_{\mathrm{Flow\_re\_ve\_rx}}}{\mathrm{Vve}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Cpo}_{\mathrm{rx}}}{\mathrm{d} t} &= \frac{-v_{\mathrm{Flow\_po\_li\_rx}} - v_{\mathrm{Flow\_po\_hv\_rx}} + v_{\mathrm{Flow\_gu\_po\_rx}}}{\mathrm{Vpo}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Chv}_{\mathrm{rx}}}{\mathrm{d} t} &= \frac{v_{\mathrm{Flow\_arli\_hv\_rx}} + v_{\mathrm{Flow\_po\_hv\_rx}} + v_{\mathrm{Flow\_li\_hv\_rx}} - v_{\mathrm{Flow\_hv\_ve\_rx}}}{\mathrm{Vhv}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Aurine}_{\mathrm{rx}}}{\mathrm{d} t} &= v_{\mathrm{KI\_\_RXEX}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Afeces}_{\mathrm{rx}}}{\mathrm{d} t} &= v_{\mathrm{GU\_\_RXEXC\_4}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Clumen}_{\mathrm{rx}}}{\mathrm{d} t} &= \frac{v_{\mathrm{LI\_\_RXBEX}} + v_{\mathrm{GU\_\_RXEFF}} - v_{\mathrm{GU\_\_RXEXC}}}{\mathrm{Vgulumen}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{LI}_{\mathrm{\_rap}}}{\mathrm{d} t} &= \frac{v_{\mathrm{LI\_\_RAPIM}} - v_{\mathrm{LI\_\_RAP2RX}}}{\mathrm{Vli}_{\mathrm{tissue}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{LI}_{\mathrm{\_rx}}}{\mathrm{d} t} &= \frac{v_{\mathrm{LI\_\_RAP2RX}} - v_{\mathrm{LI\_\_RXEX}} - v_{\mathrm{LI\_\_RXBEX}}}{\mathrm{Vli}_{\mathrm{tissue}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{GU}_{\mathrm{\_rap\_entero}}}{\mathrm{d} t} &= \frac{v_{\mathrm{GU\_\_RAPABS}} - v_{\mathrm{GU\_\_RAPBASO}} - v_{\mathrm{GU\_\_RAP2RX}} - v_{\mathrm{GU\_\_RAPPG}}}{\mathrm{GU}_{\mathrm{\_Ventero}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{GU}_{\mathrm{\_rx\_entero}}}{\mathrm{d} t} &= \frac{v_{\mathrm{GU\_\_RAP2RX}} - v_{\mathrm{GU\_\_RXEFF}}}{\mathrm{GU}_{\mathrm{\_Ventero}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{GU}_{\mathrm{\_rap\_colon}}}{\mathrm{d} t} &= \frac{v_{\mathrm{GU\_\_RAPEXC}} - v_{\mathrm{GU\_\_RAPFEC}}}{\mathrm{Vgulumen}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{GU}_{\mathrm{\_rx\_int\_0}}}{\mathrm{d} t} &= \frac{v_{\mathrm{GU\_\_RXEXC}} - v_{\mathrm{GU\_\_RXEXC\_0}}}{\mathrm{GU}_{\mathrm{\_Vint\_0}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{GU}_{\mathrm{\_rx\_int\_1}}}{\mathrm{d} t} &= \frac{v_{\mathrm{GU\_\_RXEXC\_0}} - v_{\mathrm{GU\_\_RXEXC\_1}}}{\mathrm{GU}_{\mathrm{\_Vint\_1}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{GU}_{\mathrm{\_rx\_int\_2}}}{\mathrm{d} t} &= \frac{v_{\mathrm{GU\_\_RXEXC\_1}} - v_{\mathrm{GU\_\_RXEXC\_2}}}{\mathrm{GU}_{\mathrm{\_Vint\_2}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{GU}_{\mathrm{\_rx\_int\_3}}}{\mathrm{d} t} &= \frac{v_{\mathrm{GU\_\_RXEXC\_2}} - v_{\mathrm{GU\_\_RXEXC\_3}}}{\mathrm{GU}_{\mathrm{\_Vint\_3}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{GU}_{\mathrm{\_rx\_int\_4}}}{\mathrm{d} t} &= \frac{v_{\mathrm{GU\_\_RXEXC\_3}} - v_{\mathrm{GU\_\_RXEXC\_4}}}{\mathrm{GU}_{\mathrm{\_Vint\_4}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{PODOSE}_{\mathrm{rap}}}{\mathrm{d} t} &= -v_{\mathrm{GU\_\_dissolution\_rap}} \cdot \mathrm{Mr}_{\mathrm{rap}} \qquad \text{(rate rule)} \\[1ex]
\frac{\mathrm{d} \,\mathrm{IVDOSE}_{\mathrm{rap}}}{\mathrm{d} t} &= -v_{\mathrm{iv\_rap}} \cdot \mathrm{Mr}_{\mathrm{rap}} + \mathrm{Ri}_{\mathrm{rap}} \qquad \text{(rate rule)}
\end{aligned}
$$
