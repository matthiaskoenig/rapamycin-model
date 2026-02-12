
from typing import Dict

from sbmlsim.plot import Axis, Figure, Plot
from sbmlsim.simulation import Timecourse, TimecourseSim

from pkdb_models.models.rapamycin.experiments.base_experiment import (
    RapamycinSimulationExperiment,
)
from pkdb_models.models.rapamycin.helpers import run_experiments


class DoseDependencyExperiment(RapamycinSimulationExperiment):
    """Tests po application."""

    routes = {
        "rap": ["IV", "PO"],
    }
    doses = [0, 10, 20, 40, 80]  # [mg]
    colors = ["black", "tab:orange", "tab:blue", "tab:red", "tab:green"]

    def simulations(self) -> Dict[str, TimecourseSim]:
        Q_ = self.Q_
        tcsims = {}

        for substance, routes in self.routes.items():
            for route in routes:
                for dose in self.doses:
                    tcsims[f"rap_{substance}_{route}_{dose}"] = TimecourseSim(
                        Timecourse(
                            start=0,
                            end=1 * 24 * 60,  # [min]
                            steps=1000,
                            changes={
                                **self.default_changes(),
                                f"{route}DOSE_{substance}": Q_(dose, "mg"),
                            },
                        )
                    )


        return tcsims

    def figures(self) -> Dict[str, Figure]:
        return {
            **self.figure_pk(),
        }

    def figure_pk(self) -> Dict[str, Figure]:
        figures = {}
        for substance, routes in self.routes.items():
            for route in routes:

                fig = Figure(
                    experiment=self,
                    sid=f"Fig_dose_dependency_pk_{substance}_{route}",
                    num_rows=3,
                    num_cols=2,
                    name=f"Dose Dependency {substance}_{route}",
                )
                plots = fig.create_plots(xaxis=Axis("time", unit="hr"), legend=True)
                sids = [
                    # plasma
                    "[Cve_rap]",
                    "[Cve_rx]",

                    "[Cveblood_rap]",
                    "[Cveblood_raptot]",

                    # urine,
                    "Aurine_rx",

                    # feces
                    "Afeces_rx",
                ]
                for ksid, sid in enumerate(sids):
                    if sid:
                        plots[ksid].set_yaxis(label=self.labels[sid], unit=self.units[sid])

                for ksid, sid in enumerate(sids):
                    if sid:
                        for kval, dose in enumerate(self.doses):
                            plots[ksid].add_data(
                                task=f"task_rap_{substance}_{route}_{dose}",
                                xid="time",
                                yid=sid,
                                label=f"{dose} mg",
                                color=self.colors[kval],
                            )

                figures[fig.sid] = fig

        return figures


if __name__ == "__main__":
    run_experiments(DoseDependencyExperiment, output_dir=DoseDependencyExperiment.__name__)