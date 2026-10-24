"""Base class of the rapamycin simulation experiments."""

from collections.abc import Iterable, Mapping
from typing import Any, ClassVar, NamedTuple

import pandas as pd
from pkpdsim.experiment import PKPDSimulationExperiment
from pkpdsim.impairment import RENAL_COLORS_END_STAGE, RENAL_FUNCTION_GFR
from sbmlsim import Q
from sbmlsim.units import Quantity

from pkdb_models.models.rapamycin import MODEL
from pkdb_models.models.rapamycin.rapamycin_pk import calculate_rapamycin_pk


class MolecularWeights(NamedTuple):
    """Molecular weights for the conversion of the data."""

    rap: Quantity


class RapamycinSimulationExperiment(PKPDSimulationExperiment):
    """Base class of the rapamycin simulation experiments."""

    pkpd_model = MODEL
    selections = (
        "time",
        # dosing
        "IVDOSE_rap",
        "PODOSE_rap",
        # venous blood and plasma
        "[Cveblood_rap]",
        "[Cveplasma_rap]",
        "[Cve_rx]",
        "[Cveblood_raptot]",
        "[Cveplasma_raptot]",
        # urine
        "Aurine_rx",
        # feces
        "Afeces_rap",
        "Afeces_rx",
        "Afeces_raptot",
        # cases
        "KI__f_renal_function",
        "f_cirrhosis",
        "GU__fed",
        "GU__f_pg",
        "GU__f_cyp3a4",
        "GU__f_cyp3a5",
        "LI__f_cyp3a4",
        "LI__f_cyp3a5",
    )

    font: ClassVar[dict[str, Any]] = {"weight": "bold", "size": 22}
    scan_font: ClassVar[dict[str, Any]] = {"weight": "bold", "size": 15}
    tick_font_size = 15
    legend_font_size = 9
    suptitle_font_size = 25

    # labels
    label_time = "time"

    # rapamycin
    label_rap = "rapamycin"
    label_rx = "rapamycin metabolites \n"
    label_raptot = "rapamycin and metabolites \n"
    label_rap_plasma = label_rap + " plasma"
    label_rap_blood = label_rap + " blood"
    label_raptot_plasma = label_raptot + " plasma"
    label_raptot_blood = label_raptot + " blood"

    label_rx_urine = label_rx + "urine"
    label_rx_feces = label_rx + "feces"
    label_rap_feces = label_rap + " feces"
    label_raptot_feces = label_raptot + "feces"

    labels: ClassVar[Mapping[str, str]] = {
        "time": "time",
        # rapamycin
        "[Cveplasma_rap]": label_rap_plasma,
        "[Cve_rx]": label_rx + " blood",
        "[Cveplasma_raptot]": label_raptot_plasma,
        "[Cveblood_rap]": label_rap_blood,
        "[Cveblood_raptot]": label_raptot_blood,
        "Aurine_rx": label_rx_urine,
        "Afeces_rap": label_rap_feces,
        "Afeces_rx": label_rx_feces,
        "Afeces_raptot": label_raptot_feces,
    }

    # units
    unit_time = "hr"
    unit_metabolite = "nM"
    unit_metabolite_urine = "µmole"
    unit_metabolite_feces = "µmole"
    unit_rap = unit_metabolite
    unit_rx = unit_metabolite
    unit_raptot = unit_metabolite
    unit_rx_urine = unit_metabolite_urine
    unit_rx_feces = unit_metabolite_feces

    units: ClassVar[Mapping[str, str]] = {
        "time": unit_time,
        "[Cveplasma_rap]": unit_rap,
        "[Cve_rx]": unit_rx,
        "[Cveplasma_raptot]": unit_raptot,
        "[Cveblood_rap]": unit_rap,
        "[Cveblood_raptot]": unit_raptot,
        "Aurine_rx": unit_rx_urine,
        "Afeces_rap": unit_rx_feces,
        "Afeces_rx": unit_rx_feces,
        "Afeces_raptot": unit_rx_feces,
    }
    # ----------- Food effect --------------
    #: the fed state of the intestine (`GU__fed`): the effects of a high fat
    #: meal on the dissolution and the fraction absorbed are parameters of the
    #: model (`GU__FED_f_dissolution`, `GU__FED_f_absorption`); a light
    #: breakfast is assumed to have half the effect
    fasting_map: ClassVar[dict[str, float]] = {
        "not reported": 0.0,  # assuming fasted state if nothing is reported
        "fasted": 0.0,
        "light breakfast": 0.5,
        "fed": 1.0,
    }

    # ----------- Variant activities --------------
    #: the activity of CYP3A4 per allele (`f_cyp3a4`, *1/*1 = 1): *1G has an
    #: increased activity, CL/F of *1G/*1G is 1.66 times the one of *1/*1
    #: (Zhang2017)
    cyp3a4_allele_activity: ClassVar[dict[str, float]] = {
        "*1": 1.0,
        "*1G": 1.66,
    }
    #: the functional CYP3A5 per allele (`f_cyp3a5`): *3 is a null allele (a
    #: splicing defect, no functional protein), *3/*3 is the reference of the
    #: model (0); the capacity of CYP3A5 *1/*1 relative to CYP3A4 is the model
    #: parameter `r_cyp3a5`
    cyp3a5_allele_activity: ClassVar[dict[str, float]] = {
        "*1": 1.0,
        "*3": 0.0,
    }
    cyp3a4_activity: ClassVar[dict[str, float]] = {
        "*1/*1": (cyp3a4_allele_activity["*1"] + cyp3a4_allele_activity["*1"]) / 2,
        "*1/*1G": (cyp3a4_allele_activity["*1"] + cyp3a4_allele_activity["*1G"]) / 2,
        "*1G/*1G": (cyp3a4_allele_activity["*1G"] + cyp3a4_allele_activity["*1G"]) / 2,
    }

    cyp3a5_activity: ClassVar[dict[str, float]] = {
        "*1/*1": (cyp3a5_allele_activity["*1"] + cyp3a5_allele_activity["*1"]) / 2,
        "*1/*3": (cyp3a5_allele_activity["*1"] + cyp3a5_allele_activity["*3"]) / 2,
        "*3/*3": (cyp3a5_allele_activity["*3"] + cyp3a5_allele_activity["*3"]) / 2,
    }

    # ----------- Renal map --------------
    renal_map = RENAL_FUNCTION_GFR
    renal_colors = RENAL_COLORS_END_STAGE

    # molecular weights of rapamycin
    Mr: ClassVar[MolecularWeights] = MolecularWeights(rap=Q(914.1719, "g/mole"))

    # --- Pharmacokinetic parameters ---
    def calculate_rapamycin_pk(
        self, scans: Iterable[str] | None = None
    ) -> dict[str, pd.DataFrame]:
        """Calculate the pharmacokinetic parameters of the scans.

        Args:
            scans: keys of the simulations, all simulations by default.

        Returns:
            The parameters of every scan by simulation key, per draw of the
            fitted parameters in a run over them (`uncertainty`).
        """
        keys = self._simulations if scans is None else scans
        across = None if self.uncertainty is None else self.uncertainty.id
        return {
            key: calculate_rapamycin_pk(
                experiment=self, res=self.results[f"task_{key}"], across=across
            )
            for key in keys
        }
