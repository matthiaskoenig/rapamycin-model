"""Rapamycin: PBPK model of rapamycin (sirolimus).

The paths of the model and its model object `MODEL`, which the commands
`run_simulation`, `run_fit`, `run_report` and `run_identifiability` of pkpdsim
run by its package `pkdb_models.models.rapamycin`.
"""

from pathlib import Path

from pkpdsim.model import FigureSettings, PKPDModel, SimulatorSettings
from pkpdsim.release.settings import Author, ReleaseSettings

from pkdb_models.models.rapamycin.tolerances import ABSOLUTE_TOLERANCE

RAPAMYCIN_PATH = Path(__file__).parent

MODEL_BASE_PATH = RAPAMYCIN_PATH / "models" / "output" / "models"
MODEL_PATH = MODEL_BASE_PATH / "rapamycin_body_flat.xml"

DATA_PATH_BASE = RAPAMYCIN_PATH / "data"
DATA_PATH_RAPAMYCIN = DATA_PATH_BASE / "rapamycin"
DATA_PATHS = [
    DATA_PATH_RAPAMYCIN,
]

MODEL = PKPDModel(
    name="rapamycin",
    package=__name__,
    title="Rapamycin PBPK model",
    path=RAPAMYCIN_PATH,
    model_path=MODEL_PATH,
    data_paths=tuple(DATA_PATHS),
    archive_path=MODEL_BASE_PATH.parent / "rapamycin_model.omex",
    factory=".models.factory:create_rapamycin_models",
    experiments=".simulations:EXPERIMENTS",
    fit_definitions=".fitting.fitting:FIT_DEFINITIONS",
    simulator=SimulatorSettings(
        absolute_tolerance=ABSOLUTE_TOLERANCE, relative_tolerance=1e-10
    ),
    figures=FigureSettings(formats=("svg",), fig_dpi=300, legend_fontsize=10),
    fits=("PK_20261017_094949__24b18",),
    fit_identifiability=True,
    # the release of the model repository
    release=ReleaseSettings(
        title="Physiologically based pharmacokinetic (PBPK) model of rapamycin",
        repository="matthiaskoenig/rapamycin-model",
        doi="10.5281/zenodo.17091694",
        authors=(
            Author(
                "Jesionek, Monika",
                "0009-0003-3110-3440",
                "Jagiellonian University Medical College in Krakow, Poland",
            ),
            Author(
                "König, Matthias",
                "0000-0003-1725-179X",
                (
                    "Humboldt-Universität zu Berlin, Faculty of Life Sciences, Department of Biology, Institute of Theoretical Biology, Systems Medicine of the Liver; University of Stuttgart, Institute of Structural Mechanics and Dynamics in Aerospace Engineering; University of Southern Denmark, Clinical Pharmacology, Pharmacy and Environmental Medicine, Department of Public Health; University Hospital Schleswig-Holstein, Campus Lübeck, First Department of Medicine, Lübeck, Germany"
                ),
            ),
        ),
        keywords=(
            "model",
            "rapamycin",
            "sirolimus",
            "mTOR inhibitor",
            "pharmacokinetics",
            "PBPK",
            "physiologically based pharmacokinetic model",
            "SBML",
        ),
        docker_image="matthiaskoenig/rapamycin",
        funding=(
            "Matthias König was supported by the Federal Ministry of Research, Technology and Space (BMFTR, Germany) within ATLAS by grant number 031L0304B and by the German Research Foundation (DFG) within the Research Unit Program FOR 5151 QuaLiPerF (Quantifying Liver Perfusion-Function Relationship in Complex Resection - A Systems Medicine Approach) by grant number 436883643 and by grant number 465194077 (Priority Programme SPP 2311, Subproject SimLivA). This work was supported by the BMBF-funded de.NBI Cloud within the German Network for Bioinformatics Infrastructure (de.NBI) (031A537B, 031A533A, 031A538A, 031A533B, 031A535A, 031A537C, 031A534A, 031A532B)."
        ),
        description=(
            "This repository provides the rapamycin (sirolimus) physiologically"
            " based pharmacokinetic (PBPK) model: its oral absorption, its"
            " metabolism by CYP3A4 and CYP3A5 in the intestine and the liver,"
            " its distribution in whole blood and the tissues and the excretion"
            " of its metabolites."
        ),
        first_year=2025,
    ),
)
"""The rapamycin model: its commands, see `pkpdsim.model`."""
