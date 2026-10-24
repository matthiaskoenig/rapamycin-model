# Model for hepatic rapamycin metabolism.

Model `rapamycin_liver`, SBML Level 3 Version 2, read from rapamycin\_liver.xml, written by sbmlode 0.3.0.

Model for hepatic rapamycin metabolism

Rapamycin equilibrates fast between the blood and the liver and is metabolized by CYP3A4 and CYP3A5 (RAP2RX), linear in the concentrations of the therapeutic doses. The concentrations are those of whole blood: about 95 % of the rapamycin of the blood is bound in the erythrocytes (FKBP12), the intrinsic clearance is a clearance of whole blood. Rapamycin has a low hepatic extraction (blood clearance about 3 l/hr against a hepatic blood flow of about 90 l/hr).

The metabolites are released into the blood (RXEX) and excreted into the bile (RXBEX), from where they reach the lumen of the intestine and the feces.

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
| $\mathrm{Vext}$ | `Vext` | blood | $1.5$ | l | ✓ |
| $\mathrm{Vli}$ | `Vli` | liver | $1.5$ | l | ✓ |
| $\mathrm{Vmem}$ | `Vmem` | plasma membrane | $\mathrm{NaN}$ | m<sup>2</sup> | ✓ |
| $\mathrm{Vapical}$ | `Vapical` | apical membrane | $\mathrm{NaN}$ | m<sup>2</sup> | ✓ |
| $\mathrm{Vlumen}$ | `Vlumen` | intestinal lumen (inner part of intestine) | $1.15425$ | l |  |

## Species

| Symbol | Id | Name | Compartment | Value | Unit | Properties |
| --- | --- | --- | --- | --- | --- | --- |
| $\mathrm{rap}_{\mathrm{ext}}$ | `rap_ext` | rapamycin (blood) | $\mathrm{Vext}$ | $0$ | mmol/l | concentration |
| $\mathrm{rx}_{\mathrm{ext}}$ | `rx_ext` | rapamycin metabolites (blood) | $\mathrm{Vext}$ | $0$ | mmol/l | concentration |
| $\mathrm{rap}$ | `rap` | rapamycin (liver) | $\mathrm{Vli}$ | $0$ | mmol/l | concentration |
| $\mathrm{rx}$ | `rx` | rapamycin metabolites (liver) | $\mathrm{Vli}$ | $0$ | mmol/l | concentration |
| $\mathrm{rx}_{\mathrm{lumen}}$ | `rx_lumen` | rapamycin metabolites (lumen) | $\mathrm{Vlumen}$ | $0$ | mmol/l | concentration |

## Parameters

| Symbol | Id | Name | Value | Unit | Constant |
| --- | --- | --- | --- | --- | :---: |
| $\mathrm{RXBEX}_{\mathrm{k}}$ | `RXBEX_k` | rate of the biliary excretion of the rapamycin metabolites | $0.15$ | 1/min | ✓ |
| $\mathrm{RAPIM}_{\mathrm{k}}$ | `RAPIM_k` | rate rapamycin import | $100$ | 1/min | ✓ |
| $\mathrm{RAP2RX}_{\mathrm{k}}$ | `RAP2RX_k` | rate of rapamycin metabolism by CYP3A4 | $0.0261505045623823$ | 1/min | ✓ |
| $r_{\mathrm{cyp3a5}}$ | `r_cyp3a5` | capacity of CYP3A5 relative to CYP3A4 | $1$ | \- | ✓ |
| $f_{\mathrm{cyp3a4}}$ | `f_cyp3a4` | CYP3A4 activity | $1$ | \- | ✓ |
| $f_{\mathrm{cyp3a5}}$ | `f_cyp3a5` | CYP3A5 activity | $0$ | \- | ✓ |
| $\mathrm{RXEX}_{\mathrm{k}}$ | `RXEX_k` | rate rapamycin metabolites export | $100$ | 1/min | ✓ |

## Reactions

| Rate | Id | Name | Equation |
| --- | --- | --- | --- |
| $v_{\mathrm{RAPIM}}$ | `RAPIM` | rapamycin import (RAPIM) | $\mathrm{rap}_{\mathrm{ext}} \rightleftharpoons \mathrm{rap}$ |
| $v_{\mathrm{RAP2RX}}$ | `RAP2RX` | rapamycin metabolism (CYP3A4/5) | $\mathrm{rap} \longrightarrow \mathrm{rx}$ |
| $v_{\mathrm{RXEX}}$ | `RXEX` | rapamycin metabolites export (RXEX) | $\mathrm{rx} \rightleftharpoons \mathrm{rx}_{\mathrm{ext}}$ |
| $v_{\mathrm{RXBEX}}$ | `RXBEX` | rapamycin metabolites biliary excretion | $\mathrm{rx} \longrightarrow \mathrm{rx}_{\mathrm{lumen}}$ |

The rates of the reactions are:

$$
\begin{aligned}
v_{\mathrm{RAPIM}} &= \mathrm{RAPIM}_{\mathrm{k}} \cdot \mathrm{Vli} \cdot \mathopen{}\left(\mathrm{rap}_{\mathrm{ext}} - \mathrm{rap}\right) \\[1ex]
v_{\mathrm{RAP2RX}} &= \mathopen{}\left(f_{\mathrm{cyp3a4}} + r_{\mathrm{cyp3a5}} \cdot f_{\mathrm{cyp3a5}}\right) \cdot \mathrm{RAP2RX}_{\mathrm{k}} \cdot \mathrm{Vli} \cdot \mathrm{rap} \\[1ex]
v_{\mathrm{RXEX}} &= \mathrm{RXEX}_{\mathrm{k}} \cdot \mathrm{Vli} \cdot \mathopen{}\left(\mathrm{rx} - \mathrm{rx}_{\mathrm{ext}}\right) \\[1ex]
v_{\mathrm{RXBEX}} &= \mathrm{RXBEX}_{\mathrm{k}} \cdot \mathrm{Vli} \cdot \mathrm{rx}
\end{aligned}
$$

## ODE system

The states change in time with the rates of the reactions:

$$
\begin{aligned}
\frac{\mathrm{d} \,\mathrm{rap}_{\mathrm{ext}}}{\mathrm{d} t} &= \frac{-v_{\mathrm{RAPIM}}}{\mathrm{Vext}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{rx}_{\mathrm{ext}}}{\mathrm{d} t} &= \frac{v_{\mathrm{RXEX}}}{\mathrm{Vext}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{rap}}{\mathrm{d} t} &= \frac{v_{\mathrm{RAPIM}} - v_{\mathrm{RAP2RX}}}{\mathrm{Vli}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{rx}}{\mathrm{d} t} &= \frac{v_{\mathrm{RAP2RX}} - v_{\mathrm{RXEX}} - v_{\mathrm{RXBEX}}}{\mathrm{Vli}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{rx}_{\mathrm{lumen}}}{\mathrm{d} t} &= \frac{v_{\mathrm{RXBEX}}}{\mathrm{Vlumen}}
\end{aligned}
$$
