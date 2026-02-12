"""Run all simulation experiments."""
import shutil
from pathlib import Path

from sbmlutils.console import console

from pkdb_models.models.rapamycin.helpers import run_experiments
from pkdb_models.models.rapamycin.experiments.studies import *
from pkdb_models.models.rapamycin.experiments.misc import *
import pkdb_models.models.rapamycin as rapamycin

from sbmlutils import log
from sbmlutils.console import console


logger = log.get_logger(__name__)

EXPERIMENTS = {
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
    "scan": [
    ]

}
EXPERIMENTS["all"] = EXPERIMENTS["studies"] + EXPERIMENTS["misc"]


def run_simulation_experiments(
        selected: str = None,
        experiment_classes: list = None,
        output_dir: Path = None
) -> None:
    """Run simulation experiments."""

    # Determine which experiments to run
    if experiment_classes is not None:
        experiments_to_run = experiment_classes
        if output_dir is None:
            output_dir = rapamycin.RESULTS_PATH_SIMULATION / "custom_selection"
    elif selected:
        # Using the 'selected' parameter
        if selected not in EXPERIMENTS:
            console.rule(style="red bold")
            console.print(
                f"[red]Error: Unknown group '{selected}'. Valid groups: {', '.join(EXPERIMENTS.keys())}[/red]"
            )
            console.rule(style="red bold")
            return
        experiments_to_run = EXPERIMENTS[selected]
        if output_dir is None:
            output_dir = rapamycin.RESULTS_PATH_SIMULATION / selected
    else:
        console.print("\n[red bold]Error: No experiments specified![/red bold]")
        console.print("[yellow]Use selected='all' or selected='studies' or provide experiment_classes=[...][/yellow]\n")
        return

    # Run the experiments
    run_experiments(experiment_classes=experiments_to_run, output_dir=output_dir)

    # Collect figures into one folder
    figures_dir = output_dir / "_figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    for f in output_dir.glob("**/*.png"):
        if f.parent == figures_dir:
            continue
        try:
            shutil.copy2(f, figures_dir / f.name)
        except Exception as err:
            print(f"file {f.name} in {f.parent} fails, skipping. Error: {err}")
    console.print(f"Figures copied to: file://{figures_dir}", style="info")


if __name__ == "__main__":
    """
    # Run experiments

    # selected = "renal_impairment"
    # selected = "renal_impairment"
    # selected = "cardiac_impairment"
    # selected = "hepatic_impairment"
    # selected = "pharmacodynamics"
    # selected = "all"
    # selected = "scan"
    """

    run_simulation_experiments(selected="all")
