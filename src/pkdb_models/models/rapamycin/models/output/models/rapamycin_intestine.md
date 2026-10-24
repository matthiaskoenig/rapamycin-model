# Model for rapamycin absorption in the small intestine

Model `rapamycin_intestine`, SBML Level 3 Version 2, read from rapamycin\_intestine.xml, written by sbmlode 0.3.0.

Model for rapamycin absorption

The dose dissolves in the stomach and passes into the lumen of the intestine (dissolution\_rap). A fraction F\_rap\_abs of the rapamycin which leaves the lumen is absorbed into the enterocytes, the rest is excreted with the feces (RAPEXC). In the enterocytes rapamycin is metabolized by CYP3A4 and CYP3A5 (RAP2RX), pumped back into the lumen by P-glycoprotein (RAPPG) or transported into the blood (RAPBASO); the metabolites formed in the enterocytes are released into the lumen and excreted with the feces after the transit through the colon.

The fraction absorbed and the metabolism of the gut wall set the oral bioavailability of about 14 % (oral solution, Rapamune label); the radioactivity of an oral dose of 14C-rapamycin is recovered to 91 % in the feces, mainly as metabolites (Leung2006).

Food (fed) slows the dissolution and the gastric emptying (FED\_f\_dissolution) and increases the fraction absorbed (FED\_f\_absorption): a high fat meal increases the AUC of rapamycin by 35 % and delays its tmax from 0.8 to 3.1 hr (Zimmerman1999).

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
| $\mathrm{Vext}$ | `Vext` | blood | $1$ | l | ✓ |
| $\mathrm{Vgu}$ | `Vgu` | intestine | $1.2825$ | l | ✓ |
| $\mathrm{Vlumen}$ | `Vlumen` | intestinal lumen (inner part of intestine) | $1.15425$ | l |  |
| $\mathrm{Vfeces}$ | `Vfeces` | feces | $1$ | l | ✓ |
| $\mathrm{Ventero}$ | `Ventero` | intestinal lining (enterocytes) |  | l |  |
| $\mathrm{Vapical}$ | `Vapical` | apical membrane (intestinal membrane enterocytes) | $\mathrm{NaN}$ | m<sup>2</sup> | ✓ |
| $\mathrm{Vbaso}$ | `Vbaso` | basolateral membrane (intestinal membrane enterocytes) | $\mathrm{NaN}$ | m<sup>2</sup> | ✓ |
| $\mathrm{Vstomach}$ | `Vstomach` | stomach | $1$ | l | ✓ |
| $\mathrm{Vint}_{0}$ | `Vint_0` | intestinal lumen 0 (inner part of intestine) |  | l |  |
| $\mathrm{Vint}_{1}$ | `Vint_1` | intestinal lumen 1 (inner part of intestine) |  | l |  |
| $\mathrm{Vint}_{2}$ | `Vint_2` | intestinal lumen 2 (inner part of intestine) |  | l |  |
| $\mathrm{Vint}_{3}$ | `Vint_3` | intestinal lumen 3 (inner part of intestine) |  | l |  |
| $\mathrm{Vint}_{4}$ | `Vint_4` | intestinal lumen 4 (inner part of intestine) |  | l |  |

## Species

| Symbol | Id | Name | Compartment | Value | Unit | Properties |
| --- | --- | --- | --- | --- | --- | --- |
| $\mathrm{rap}_{\mathrm{stomach}}$ | `rap_stomach` | rapamycin (stomach) | $\mathrm{Vstomach}$ |  | mmol | amount, boundary |
| $\mathrm{rap}_{\mathrm{lumen}}$ | `rap_lumen` | rapamycin (intestinal volume) | $\mathrm{Vlumen}$ | $0$ | mmol/l | concentration |
| $\mathrm{rx}_{\mathrm{lumen}}$ | `rx_lumen` | rapamycin metabolites (intestinal volume) | $\mathrm{Vlumen}$ | $0$ | mmol/l | concentration |
| $\mathrm{rap}_{\mathrm{entero}}$ | `rap_entero` | rapamycin (enterocytes) | $\mathrm{Ventero}$ | $0$ | mmol/l | concentration |
| $\mathrm{rx}_{\mathrm{entero}}$ | `rx_entero` | rapamycin metabolites (enterocytes) | $\mathrm{Ventero}$ | $0$ | mmol/l | concentration |
| $\mathrm{rap}_{\mathrm{ext}}$ | `rap_ext` | rapamycin (blood) | $\mathrm{Vext}$ | $0$ | mmol/l | concentration |
| $\mathrm{rap}_{\mathrm{colon}}$ | `rap_colon` | rapamycin (colon) | $\mathrm{Vlumen}$ | $0$ | mmol/l | concentration |
| $\mathrm{rap}_{\mathrm{feces}}$ | `rap_feces` | rapamycin (feces) | $\mathrm{Vfeces}$ |  | mmol | amount |
| $\mathrm{rx}_{\mathrm{feces}}$ | `rx_feces` | rapamycin metabolites (feces) | $\mathrm{Vfeces}$ |  | mmol | amount |
| $\mathrm{rx}_{\mathrm{int\_0}}$ | `rx_int_0` | rapamycin metabolites (intestine) 0 | $\mathrm{Vint}_{0}$ | $0$ | mmol/l | concentration |
| $\mathrm{rx}_{\mathrm{int\_1}}$ | `rx_int_1` | rapamycin metabolites (intestine) 1 | $\mathrm{Vint}_{1}$ | $0$ | mmol/l | concentration |
| $\mathrm{rx}_{\mathrm{int\_2}}$ | `rx_int_2` | rapamycin metabolites (intestine) 2 | $\mathrm{Vint}_{2}$ | $0$ | mmol/l | concentration |
| $\mathrm{rx}_{\mathrm{int\_3}}$ | `rx_int_3` | rapamycin metabolites (intestine) 3 | $\mathrm{Vint}_{3}$ | $0$ | mmol/l | concentration |
| $\mathrm{rx}_{\mathrm{int\_4}}$ | `rx_int_4` | rapamycin metabolites (intestine) 4 | $\mathrm{Vint}_{4}$ | $0$ | mmol/l | concentration |

## Parameters

| Symbol | Id | Name | Value | Unit | Constant |
| --- | --- | --- | --- | --- | :---: |
| $f_{\mathrm{entero}}$ | `f_entero` | fraction of the enterocytes of the intestine volume | $0.1$ | \- | ✓ |
| $\mathrm{Vchain}$ | `Vchain` | volume of chain compartment | $0.1$ | l | ✓ |
| $F_{\mathrm{rap\_abs}}$ | `F_rap_abs` | fraction absorbed rapamycin (fasted) | $0.5$ | \- | ✓ |
| $\mathrm{RAPABS}_{\mathrm{k}}$ | `RAPABS_k` | rate of rapamycin absorption | $0.05$ | 1/min | ✓ |
| $\mathrm{fed}$ | `fed` | fed state (0 fasted, 1 fed) | $0$ | \- | ✓ |
| $\mathrm{FED}_{\mathrm{f\_absorption}}$ | `FED_f_absorption` | fed effect on the fraction absorbed | $1.20035322885454$ | \- | ✓ |
| $\mathrm{FED}_{\mathrm{f\_dissolution}}$ | `FED_f_dissolution` | fed effect on the dissolution and gastric emptying | $0.171348443283468$ | \- | ✓ |
| $f_{\mathrm{absorption}}$ | `f_absorption` | factor of the fraction absorbed |  | \- |  |
| $f_{\mathrm{dissolution}}$ | `f_dissolution` | factor of the dissolution and gastric emptying |  | \- |  |
| $\mathrm{absorption}$ | `absorption` | rapamycin leaving the lumen |  | mmol/min |  |
| $\mathrm{RAPBASO}_{\mathrm{k}}$ | `RAPBASO_k` | rate of rapamycin transport from the enterocytes into the blood | $0.05$ | 1/min | ✓ |
| $\mathrm{RAP2RX}_{\mathrm{k}}$ | `RAP2RX_k` | rate of rapamycin metabolism by CYP3A4 (gut wall) | $0.1$ | 1/min | ✓ |
| $r_{\mathrm{cyp3a5}}$ | `r_cyp3a5` | capacity of CYP3A5 relative to CYP3A4 | $1$ | \- | ✓ |
| $f_{\mathrm{cyp3a4}}$ | `f_cyp3a4` | CYP3A4 activity | $1$ | \- | ✓ |
| $f_{\mathrm{cyp3a5}}$ | `f_cyp3a5` | CYP3A5 activity | $0$ | \- | ✓ |
| $\mathrm{RAPPG}_{\mathrm{k}}$ | `RAPPG_k` | rate of rapamycin efflux by P-glycoprotein | $0.02$ | 1/min | ✓ |
| $f_{\mathrm{pg}}$ | `f_pg` | P-glycoprotein activity | $1$ | \- | ✓ |
| $\mathrm{RXEFF}_{\mathrm{k}}$ | `RXEFF_k` | rate of the release of the metabolites into the lumen | $1$ | 1/min | ✓ |
| $\mathrm{RXEXC}_{\mathrm{k}}$ | `RXEXC_k` | rate of the transit of the metabolites to the feces | $0.0033$ | 1/min | ✓ |
| $\mathrm{PODOSE}_{\mathrm{rap}}$ | `PODOSE_rap` | oral dose rapamycin \[mg\] | $0$ | mg |  |
| $\mathrm{Ka}_{\mathrm{dis\_rap}}$ | `Ka_dis_rap` | Ka\_dis \[1/hr\] dissolution rapamycin | $3.64645771834172$ | 1/hr | ✓ |
| $\mathrm{Mr}_{\mathrm{rap}}$ | `Mr_rap` | Molecular weight rapamycin \[g/mole\] | $914.1719$ | g/mol | ✓ |

## Initial assignments and assignment rules

The initial assignments set the values at $t = 0$:

$$
\begin{aligned}
\mathrm{rap}_{\mathrm{stomach}} &= 0 \cdot \mathrm{Vstomach} \\[1ex]
\mathrm{rap}_{\mathrm{feces}} &= 0 \cdot \mathrm{Vfeces} \\[1ex]
\mathrm{rx}_{\mathrm{feces}} &= 0 \cdot \mathrm{Vfeces}
\end{aligned}
$$

The assignment rules hold at every time $t$:

$$
\begin{aligned}
\mathrm{Vint}_{0} &= \mathrm{Vchain} \\[1ex]
\mathrm{Vint}_{1} &= \mathrm{Vchain} \\[1ex]
\mathrm{Vint}_{2} &= \mathrm{Vchain} \\[1ex]
\mathrm{Vint}_{3} &= \mathrm{Vchain} \\[1ex]
\mathrm{Vint}_{4} &= \mathrm{Vchain} \\[1ex]
\mathrm{Ventero} &= f_{\mathrm{entero}} \cdot \mathrm{Vgu} \\[1ex]
f_{\mathrm{absorption}} &= 1 + \mathrm{fed} \cdot \mathopen{}\left(\mathrm{FED}_{\mathrm{f\_absorption}} - 1\right) \\[1ex]
f_{\mathrm{dissolution}} &= 1 + \mathrm{fed} \cdot \mathopen{}\left(\mathrm{FED}_{\mathrm{f\_dissolution}} - 1\right) \\[1ex]
\mathrm{absorption} &= \mathrm{RAPABS}_{\mathrm{k}} \cdot \mathrm{Vlumen} \cdot \mathrm{rap}_{\mathrm{lumen}}
\end{aligned}
$$

## Reactions

| Rate | Id | Name | Equation |
| --- | --- | --- | --- |
| $v_{\mathrm{RAPABS}}$ | `RAPABS` | absorption rapamycin (enterocytes) | $\mathrm{rap}_{\mathrm{lumen}} \longrightarrow \mathrm{rap}_{\mathrm{entero}}$ |
| $v_{\mathrm{RAPEXC}}$ | `RAPEXC` | excretion rapamycin (colon) | $\mathrm{rap}_{\mathrm{lumen}} \longrightarrow \mathrm{rap}_{\mathrm{colon}}$ |
| $v_{\mathrm{RAPFEC}}$ | `RAPFEC` | excretion rapamycin (feces) | $\mathrm{rap}_{\mathrm{colon}} \longrightarrow \mathrm{rap}_{\mathrm{feces}}$ |
| $v_{\mathrm{RAPBASO}}$ | `RAPBASO` | transport rapamycin (blood) | $\mathrm{rap}_{\mathrm{entero}} \longrightarrow \mathrm{rap}_{\mathrm{ext}}$ |
| $v_{\mathrm{RAP2RX}}$ | `RAP2RX` | rapamycin metabolism via CYP3A4/5 | $\mathrm{rap}_{\mathrm{entero}} \longrightarrow \mathrm{rx}_{\mathrm{entero}}$ |
| $v_{\mathrm{RAPPG}}$ | `RAPPG` | efflux rapamycin by P-glycoprotein | $\mathrm{rap}_{\mathrm{entero}} \longrightarrow \mathrm{rap}_{\mathrm{lumen}}$ |
| $v_{\mathrm{RXEFF}}$ | `RXEFF` | release rapamycin metabolites (lumen) | $\mathrm{rx}_{\mathrm{entero}} \longrightarrow \mathrm{rx}_{\mathrm{lumen}}$ |
| $v_{\mathrm{RXEXC}}$ | `RXEXC` | excretion rapamycin metabolites (feces) | $\mathrm{rx}_{\mathrm{lumen}} \longrightarrow \mathrm{rx}_{\mathrm{int\_0}}$ |
| $v_{\mathrm{RXEXC\_0}}$ | `RXEXC_0` | excretion rapamycin metabolites 0 | $\mathrm{rx}_{\mathrm{int\_0}} \longrightarrow \mathrm{rx}_{\mathrm{int\_1}}$ |
| $v_{\mathrm{RXEXC\_1}}$ | `RXEXC_1` | excretion rapamycin metabolites 1 | $\mathrm{rx}_{\mathrm{int\_1}} \longrightarrow \mathrm{rx}_{\mathrm{int\_2}}$ |
| $v_{\mathrm{RXEXC\_2}}$ | `RXEXC_2` | excretion rapamycin metabolites 2 | $\mathrm{rx}_{\mathrm{int\_2}} \longrightarrow \mathrm{rx}_{\mathrm{int\_3}}$ |
| $v_{\mathrm{RXEXC\_3}}$ | `RXEXC_3` | excretion rapamycin metabolites 3 | $\mathrm{rx}_{\mathrm{int\_3}} \longrightarrow \mathrm{rx}_{\mathrm{int\_4}}$ |
| $v_{\mathrm{RXEXC\_4}}$ | `RXEXC_4` | excretion rapamycin metabolites 4 | $\mathrm{rx}_{\mathrm{int\_4}} \longrightarrow \mathrm{rx}_{\mathrm{feces}}$ |
| $v_{\mathrm{dissolution\_rap}}$ | `dissolution_rap` | dissolution rapamycin | $\mathrm{rap}_{\mathrm{stomach}} \longrightarrow \mathrm{rap}_{\mathrm{lumen}}$ |

The rates of the reactions are:

$$
\begin{aligned}
v_{\mathrm{RAPABS}} &= F_{\mathrm{rap\_abs}} \cdot f_{\mathrm{absorption}} \cdot \mathrm{absorption} \\[1ex]
v_{\mathrm{RAPEXC}} &= \mathopen{}\left(1 - F_{\mathrm{rap\_abs}} \cdot f_{\mathrm{absorption}}\right) \cdot \mathrm{absorption} \\[1ex]
v_{\mathrm{RAPFEC}} &= \frac{\mathrm{RXEXC}_{\mathrm{k}}}{6} \cdot \mathrm{Vlumen} \cdot \mathrm{rap}_{\mathrm{colon}} \\[1ex]
v_{\mathrm{RAPBASO}} &= \mathrm{RAPBASO}_{\mathrm{k}} \cdot \mathrm{Ventero} \cdot \mathrm{rap}_{\mathrm{entero}} \\[1ex]
v_{\mathrm{RAP2RX}} &= \mathopen{}\left(f_{\mathrm{cyp3a4}} + r_{\mathrm{cyp3a5}} \cdot f_{\mathrm{cyp3a5}}\right) \cdot \mathrm{RAP2RX}_{\mathrm{k}} \cdot \mathrm{Ventero} \cdot \mathrm{rap}_{\mathrm{entero}} \\[1ex]
v_{\mathrm{RAPPG}} &= f_{\mathrm{pg}} \cdot \mathrm{RAPPG}_{\mathrm{k}} \cdot \mathrm{Ventero} \cdot \mathrm{rap}_{\mathrm{entero}} \\[1ex]
v_{\mathrm{RXEFF}} &= \mathrm{RXEFF}_{\mathrm{k}} \cdot \mathrm{Ventero} \cdot \mathrm{rx}_{\mathrm{entero}} \\[1ex]
v_{\mathrm{RXEXC}} &= \mathrm{RXEXC}_{\mathrm{k}} \cdot \mathrm{Vlumen} \cdot \mathrm{rx}_{\mathrm{lumen}} \\[1ex]
v_{\mathrm{RXEXC\_0}} &= \mathrm{RXEXC}_{\mathrm{k}} \cdot \mathrm{Vint}_{0} \cdot \mathrm{rx}_{\mathrm{int\_0}} \\[1ex]
v_{\mathrm{RXEXC\_1}} &= \mathrm{RXEXC}_{\mathrm{k}} \cdot \mathrm{Vint}_{1} \cdot \mathrm{rx}_{\mathrm{int\_1}} \\[1ex]
v_{\mathrm{RXEXC\_2}} &= \mathrm{RXEXC}_{\mathrm{k}} \cdot \mathrm{Vint}_{2} \cdot \mathrm{rx}_{\mathrm{int\_2}} \\[1ex]
v_{\mathrm{RXEXC\_3}} &= \mathrm{RXEXC}_{\mathrm{k}} \cdot \mathrm{Vint}_{3} \cdot \mathrm{rx}_{\mathrm{int\_3}} \\[1ex]
v_{\mathrm{RXEXC\_4}} &= \mathrm{RXEXC}_{\mathrm{k}} \cdot \mathrm{Vint}_{4} \cdot \mathrm{rx}_{\mathrm{int\_4}} \\[1ex]
v_{\mathrm{dissolution\_rap}} &= \frac{\frac{f_{\mathrm{dissolution}} \cdot \mathrm{Ka}_{\mathrm{dis\_rap}}}{60} \cdot \mathrm{PODOSE}_{\mathrm{rap}}}{\mathrm{Mr}_{\mathrm{rap}}}
\end{aligned}
$$

## ODE system

The states change in time with the rates of the reactions and the rate rules:

$$
\begin{aligned}
\frac{\mathrm{d} \,\mathrm{rap}_{\mathrm{lumen}}}{\mathrm{d} t} &= \frac{-v_{\mathrm{RAPABS}} - v_{\mathrm{RAPEXC}} + v_{\mathrm{RAPPG}} + v_{\mathrm{dissolution\_rap}}}{\mathrm{Vlumen}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{rx}_{\mathrm{lumen}}}{\mathrm{d} t} &= \frac{v_{\mathrm{RXEFF}} - v_{\mathrm{RXEXC}}}{\mathrm{Vlumen}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{rap}_{\mathrm{entero}}}{\mathrm{d} t} &= \frac{v_{\mathrm{RAPABS}} - v_{\mathrm{RAPBASO}} - v_{\mathrm{RAP2RX}} - v_{\mathrm{RAPPG}}}{\mathrm{Ventero}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{rx}_{\mathrm{entero}}}{\mathrm{d} t} &= \frac{v_{\mathrm{RAP2RX}} - v_{\mathrm{RXEFF}}}{\mathrm{Ventero}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{rap}_{\mathrm{ext}}}{\mathrm{d} t} &= \frac{v_{\mathrm{RAPBASO}}}{\mathrm{Vext}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{rap}_{\mathrm{colon}}}{\mathrm{d} t} &= \frac{v_{\mathrm{RAPEXC}} - v_{\mathrm{RAPFEC}}}{\mathrm{Vlumen}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{rap}_{\mathrm{feces}}}{\mathrm{d} t} &= v_{\mathrm{RAPFEC}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{rx}_{\mathrm{feces}}}{\mathrm{d} t} &= v_{\mathrm{RXEXC\_4}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{rx}_{\mathrm{int\_0}}}{\mathrm{d} t} &= \frac{v_{\mathrm{RXEXC}} - v_{\mathrm{RXEXC\_0}}}{\mathrm{Vint}_{0}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{rx}_{\mathrm{int\_1}}}{\mathrm{d} t} &= \frac{v_{\mathrm{RXEXC\_0}} - v_{\mathrm{RXEXC\_1}}}{\mathrm{Vint}_{1}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{rx}_{\mathrm{int\_2}}}{\mathrm{d} t} &= \frac{v_{\mathrm{RXEXC\_1}} - v_{\mathrm{RXEXC\_2}}}{\mathrm{Vint}_{2}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{rx}_{\mathrm{int\_3}}}{\mathrm{d} t} &= \frac{v_{\mathrm{RXEXC\_2}} - v_{\mathrm{RXEXC\_3}}}{\mathrm{Vint}_{3}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{rx}_{\mathrm{int\_4}}}{\mathrm{d} t} &= \frac{v_{\mathrm{RXEXC\_3}} - v_{\mathrm{RXEXC\_4}}}{\mathrm{Vint}_{4}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{PODOSE}_{\mathrm{rap}}}{\mathrm{d} t} &= -v_{\mathrm{dissolution\_rap}} \cdot \mathrm{Mr}_{\mathrm{rap}} \qquad \text{(rate rule)}
\end{aligned}
$$
