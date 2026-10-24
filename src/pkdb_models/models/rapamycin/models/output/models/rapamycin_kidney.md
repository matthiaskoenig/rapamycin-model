# Model for renal rapamycin excretion.

Model `rapamycin_kidney`, SBML Level 3 Version 2, read from rapamycin\_kidney.xml, written by sbmlode 0.3.0.

Model for the renal excretion of the rapamycin metabolites

Rapamycin is not excreted by the kidneys; its metabolites are excreted in the urine: 2.2 % of the radioactivity of an oral dose of 14C-rapamycin is recovered in the urine (Leung2006).

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
| $\mathrm{Vki}$ | `Vki` | kidney | $0.3$ | l | ✓ |
| $\mathrm{Vurine}$ | `Vurine` | urine | $1$ | l | ✓ |

## Species

| Symbol | Id | Name | Compartment | Value | Unit | Properties |
| --- | --- | --- | --- | --- | --- | --- |
| $\mathrm{rx}_{\mathrm{ext}}$ | `rx_ext` | rapamycin metabolites (blood) | $\mathrm{Vext}$ | $0$ | mmol/l | concentration |
| $\mathrm{rx}_{\mathrm{urine}}$ | `rx_urine` | rapamycin metabolites (urine) | $\mathrm{Vurine}$ |  | mmol | amount |

## Parameters

| Symbol | Id | Name | Value | Unit | Constant |
| --- | --- | --- | --- | --- | :---: |
| $f_{\mathrm{renal\_function}}$ | `f_renal_function` | parameter for renal function | $1$ | \- | ✓ |
| $\mathrm{RXEX}_{\mathrm{k}}$ | `RXEX_k` | rate urinary excretion of rapamycin metabolites | $0.12$ | 1/min | ✓ |

## Initial assignments and assignment rules

The initial assignments set the values at $t = 0$:

$$
\begin{aligned}
\mathrm{rx}_{\mathrm{urine}} &= 0 \cdot \mathrm{Vurine}
\end{aligned}
$$

## Reactions

| Rate | Id | Name | Equation |
| --- | --- | --- | --- |
| $v_{\mathrm{RXEX}}$ | `RXEX` | rapamycin metabolites excretion (RXEX) | $\mathrm{rx}_{\mathrm{ext}} \longrightarrow \mathrm{rx}_{\mathrm{urine}}$ |

The rates of the reactions are:

$$
\begin{aligned}
v_{\mathrm{RXEX}} &= f_{\mathrm{renal\_function}} \cdot \mathrm{RXEX}_{\mathrm{k}} \cdot \mathrm{Vki} \cdot \mathrm{rx}_{\mathrm{ext}}
\end{aligned}
$$

## ODE system

The states change in time with the rates of the reactions:

$$
\begin{aligned}
\frac{\mathrm{d} \,\mathrm{rx}_{\mathrm{ext}}}{\mathrm{d} t} &= \frac{-v_{\mathrm{RXEX}}}{\mathrm{Vext}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{rx}_{\mathrm{urine}}}{\mathrm{d} t} &= v_{\mathrm{RXEX}}
\end{aligned}
$$
