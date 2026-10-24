"""The simulation experiments of the rapamycin model and their groups.

They are the registry of `MODEL`, the experiments its commands simulate.
"""

from sbmlsim.experiment import SimulationExperiment

from pkdb_models.models.rapamycin import MODEL
from pkdb_models.models.rapamycin.experiments.misc import DoseDependencyExperiment
from pkdb_models.models.rapamycin.experiments.scans import RapamycinParameterScan
from pkdb_models.models.rapamycin.experiments.studies import (
    BasaDenes2019,
    Bottiger2001,
    Brattstrom2000,
    Kelly1999,
    KorthBradley2012,
    Leelahavanichkul2005,
    Leung2006,
    Renders2007,
    Rogers2008,
    Tortorici2013,
    Tortorici2014,
    Wang2014,
    Zha2022,
    Zhang2017,
    Zimmerman1997,
    Zimmerman1999,
    Zimmerman2003,
    Zimmerman2005,
    Zimmerman2008,
)

EXPERIMENTS: dict[str, list[type[SimulationExperiment]]] = {
    "studies": [
        BasaDenes2019,
        Bottiger2001,
        Brattstrom2000,
        Kelly1999,
        KorthBradley2012,
        Leelahavanichkul2005,
        Leung2006,
        Renders2007,
        Rogers2008,
        Tortorici2013,
        Tortorici2014,
        Wang2014,
        Zha2022,
        Zhang2017,
        Zimmerman1997,
        Zimmerman1999,
        Zimmerman2003,
        Zimmerman2005,
        Zimmerman2008,
    ],
    "dose": [
        BasaDenes2019,
        Brattstrom2000,
        Kelly1999,
        Zimmerman1997,
    ],
    "cyp3a4_5": [
        Renders2007,
        Zhang2017,
    ],
    "food": [
        BasaDenes2019,
        Zimmerman1999,
    ],
    "hepatic_impairment": [
        Zimmerman2005,
        Zimmerman2008,
    ],
    "renal_impairment": [
        Rogers2008,
    ],
    "renal_transplant": [
        Wang2014,
    ],
    "misc": [
        DoseDependencyExperiment,
    ],
    "scans": [
        RapamycinParameterScan,
    ],
}
EXPERIMENTS["all"] = EXPERIMENTS["studies"] + EXPERIMENTS["misc"] + EXPERIMENTS["scans"]


if __name__ == "__main__":
    MODEL.run_simulation_experiments(EXPERIMENTS["all"], MODEL.simulation_dir / "all")
