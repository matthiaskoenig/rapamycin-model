"""
Reusable functionality for multiple simulation experiments.
"""

from collections import namedtuple
from typing import Dict
import pandas as pd

from pkdb_models.models.rapamycin.pk import calculate_rapamycin_pk
from pkdb_models.models.rapamycin import MODEL_PATH
from sbmlsim.experiment import SimulationExperiment
from sbmlsim.model import AbstractModel
from sbmlsim.task import Task

# Constants for conversion
MolecularWeights = namedtuple("MolecularWeights", "rap fuj")


class RapamycinSimulationExperiment(SimulationExperiment):
    """Base class for all SimulationExperiments with rapamycin."""

    font = {"weight": "bold", "size": 22}
    scan_font = {"weight": "bold", "size": 15}
    tick_font_size = 15
    legend_font_size = 9
    suptitle_font_size = 25

    # labels
    label_time = "time"

    #rapamycin
    label_rap = "rapamycin"
    label_rx = "rapamycin metabolites \n"
    label_raptot = "rapamycin and metabolites \n"
    label_rap_plasma = label_rap + " plasma"
    label_rap_blood = label_rap + " blood"
    label_raptot_plasma = label_raptot + " plasma"
    label_raptot_blood = label_raptot + " blood"

    label_rx_urine = label_rx + "urine"
    label_rx_feces = label_rx + "feces"

    labels: Dict[str, str] = {
        "time": "time",
        #rapamycin
        "[Cve_rap]": label_rap_plasma,
        "[Cve_rx]": label_rx + " plasma",
        "[Cve_raptot]": label_raptot_plasma,

        "[Cveblood_rap]": label_rap_blood,
        "[Cveblood_raptot]": label_raptot_blood,

        "Aurine_rx": label_rx_urine,
        "Afeces_rx": label_rx_feces,
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

    units: Dict[str, str] = {
        "time": unit_time,
        "[Cve_rap]": unit_rap,
        "[Cve_rx]": unit_rx,
        "[Cve_raptot]": unit_raptot,
        "[Cveblood_rap]": unit_rap,
        "[Cveblood_raptot]": unit_raptot,
        "Aurine_rx": unit_rx_urine,
        "Afeces_rx": unit_rx_feces,
    }
    # ----------- Food effect --------------
    # ----------- Fasting/food -----
    fasting_map = {
        "not reported": 1.0,  # assuming fasted state if nothing is reported
        "fasted": 1.0,
        "fed": 0.7,
    }
    fasting_colors = {
        "fasted": "black",
        "fed": "tab:red",
    }

    # ----------- Variant activities --------------
    cyp3a4_allele_activity = {
        "*1": 1.0,
        "*1G": 1.2,  # increased activity
    }
    cyp3a4_activity = {
        "*1/*1": (cyp3a4_allele_activity["*1"] + cyp3a4_allele_activity["*1"]) / 2,
        "*1/*1G": (cyp3a4_allele_activity["*1"] + cyp3a4_allele_activity["*1G"]) / 2,
        "*1G/*1G": (cyp3a4_allele_activity["*1G"] + cyp3a4_allele_activity["*1G"]) / 2,
    }
    cyp3a5_allele_activity = {
        "*1": 1.0,
        "*3": 0.5,  # decreased activity
    }
    cyp3a5_activity = {
        "*1/*1": (cyp3a5_allele_activity["*1"] + cyp3a5_allele_activity["*1"]) / 2,
        "*1/*3": (cyp3a5_allele_activity["*1"] + cyp3a5_allele_activity["*3"]) / 2,
        "*3/*3": (cyp3a5_allele_activity["*3"] + cyp3a5_allele_activity["*3"]) / 2,
    }

    # ----------- Renal map --------------
    renal_map = {
        "Normal renal function": 101.0 / 101.0,  # 1.0,
        "Mild renal impairment": 50.0 / 101.0,  # 0.5
        "Moderate renal impairment": 35.0 / 101.0,  # 0.35
        "Severe renal impairment": 20.0 / 101.0,  # 0.20
        "End stage renal disease": 10.5 / 101.0,  # 0.1
    }
    renal_colors = {
        "Normal renal function": "black",
        "Mild renal impairment": "#66c2a4",
        "Moderate renal impairment": "#2ca25f",
        "Severe renal impairment": "#006d2c",
        "End stage renal disease": "#006d5e"
    }

    # ----------- Cirrhosis map --------------
    cirrhosis_map = {
        "Control": 0,
        "Mild cirrhosis": 0.3994897959183674,  # CPT A
        "Moderate cirrhosis": 0.6979591836734694,  # CPT B
        "Severe cirrhosis": 0.8127551020408164,  # CPT C
    }
    cirrhosis_colors = {
        "Control": "black",
        "Mild cirrhosis": "#74a9cf",  # CPT A
        "Moderate cirrhosis": "#2b8cbe",  # CPT B
        "Severe cirrhosis": "#045a8d",  # CPT C
    }

    def models(self) -> Dict[str, AbstractModel]:
        Q_ = self.Q_
        return {
            "model": AbstractModel(
                source=MODEL_PATH,
                language_type=AbstractModel.LanguageType.SBML,
                changes={},
            )
        }

    @staticmethod
    def _default_changes(Q_):
        """Default changes to simulations."""

        changes = {
            #>>> !Optimal parameter 'Kp_rap' within 5% of lower bound! <<<
            'ftissue_rap': Q_(6.385157332090921, 'l/min'),  # [0.01 - 10]
            'Kp_rap': Q_(10.000004054918183, 'dimensionless'),  # [10 - 10000]
            'GU__RAPIM_k': Q_(0.010254172647171959, '1/min'),  # [0.0001 - 1]
            'GU__RAP2RX_k': Q_(0.43373578546152136, '1/min'),  # [0.0001 - 1]
            'GU__RXPG_k': Q_(0.010328378821077363, '1/min'),  # [0.0001 - 1]
            'GU__RXEXC_k': Q_(0.023347411187126337, '1/min'),  # [0.0001 - 1]
            'KI__RXEX_k': Q_(0.1693055659116988, '1/min'),  # [0.1 - 10]
            'LI__RAP2RX_Vmax': Q_(0.002858932837010464, 'mmole/min/l'),  # [0.0001 - 1]
            'LI__RXBEX_k': Q_(2.2698133102257085e-06, '1/min'),  # [1e-07 - 0.001]
        }

        return changes

    def default_changes(self: SimulationExperiment) -> Dict:
        """Default changes to simulations."""
        return RapamycinSimulationExperiment._default_changes(Q_=self.Q_)

    def tasks(self) -> Dict[str, Task]:
        if self.simulations():
            return {
                f"task_{key}": Task(model="model", simulation=key)
                for key in self.simulations()
            }
        return {}

    def data(self) -> Dict:
        self.add_selections_data(
            selections=[
                "time",

                # dosing
                "IVDOSE_rap",
                "PODOSE_rap",

                # venous plasma
                "[Cve_rap]",
                "[Cve_rx]",
                "[Cve_raptot]",
                "[Cveblood_rap]",
                "[Cveblood_raptot]",

                # urine
                "Aurine_rx",

                # feces
                "Afeces_rx",

                # cases
                'KI__f_renal_function',
                'f_cirrhosis',
                "GU__f_oatp",
                "GU__f_cyp3a4",
                "GU__f_cyp3a5",
                "LI__f_cyp3a4",
                "LI__f_cyp3a5",
            ]
        )
        return {}

    @property
    def Mr(self):
        return MolecularWeights(
            rap=self.Q_(914.1719, "g/mole"),
            fuj=self.Q_(804, "g/mole")
        )

# --- Pharmacokinetic parameters ---
    pk_labels = {
        "auc": "AUCend",
        "aucinf": "AUC",
        "cl": "Total clearance",
        "cl_renal": "Renal clearance",
        "cl_hepatic": "Hepatic clearance",
        "cmax": "Cmax",
        "thalf": "Half-life",
        "kel": "kel",
        "vd": "vd",
    }

    pk_units = {
        "auc": "µmole/l*hr",
        "aucinf": "µmole/l*hr",
        "cl": "ml/min",
        "cl_renal": "ml/min",
        "cl_hepatic": "ml/min",
        "cmax": "µmole/l",
        "thalf": "hr",
        "kel": "1/hr",
        "vd": "l",
    }

    def calculate_rapamycin_pk(self, scans: list = []) -> Dict[str, pd.DataFrame]:
       """Calculate pk parameters for simulations (scans)"""
       pk_dfs = {}
       if scans:
           for sim_key in scans:
               xres = self.results[f"task_{sim_key}"]
               df = calculate_rapamycin_pk(experiment=self, xres=xres)
               pk_dfs[sim_key] = df
       else:
           for sim_key in self._simulations.keys():
               xres = self.results[f"task_{sim_key}"]
               df = calculate_rapamycin_pk(experiment=self, xres=xres)
               pk_dfs[sim_key] = df
       return pk_dfs
