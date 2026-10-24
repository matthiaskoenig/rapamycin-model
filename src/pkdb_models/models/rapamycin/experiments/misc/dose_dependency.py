"""Simulation experiment DoseDependencyExperiment."""

from typing import ClassVar, override

import matplotlib
import numpy as np
from matplotlib.colors import to_hex
from sbmlsim import Q
from sbmlsim.plot import Figure
from sbmlsim.simulation import Simulation

from pkdb_models.models.rapamycin import MODEL
from pkdb_models.models.rapamycin.experiments.base_experiment import (
    RapamycinSimulationExperiment,
)


class DoseDependencyExperiment(RapamycinSimulationExperiment):
    """Dose dependency of single intravenous and oral doses of rapamycin.

    The doses of the studies (2 to 40 mg) are followed for 14 days, about five
    terminal half-lives (60 hr): rapamycin in whole blood, the main readout of
    the studies, and in plasma, its metabolites in blood and the excretion of
    the metabolites into the urine and of rapamycin and its metabolites into
    the feces. In a run over the draws of the fitted parameters
    (`run_simulation pkdb_models.models.rapamycin --uncertainty`) every curve
    is the median of the draws with the band of their quantiles.
    """

    #: the routes of the application
    routes: ClassVar[dict[str, str]] = {"PO": "oral", "IV": "intravenous"}
    #: the single doses of the studies [mg]
    doses: ClassVar[list[float]] = [2, 5, 10, 15, 20, 40]
    #: the color of every dose, light to dark with the dose
    colors: ClassVar[dict[float, str]] = {
        dose: to_hex(matplotlib.colormaps["YlGnBu"](value))
        for dose, value in zip(
            doses, np.linspace(0.35, 1.0, num=len(doses)), strict=True
        )
    }
    #: the end of the simulations [hr], 14 days
    tend = 14 * 24
    #: the output times [hr]: every 0.05 hr in the first day (the absorption
    #: and the maximum), every 0.5 hr after it; one grid for the simulations
    #: and the draws of the fitted parameters
    times: ClassVar[tuple[float, ...]] = tuple(
        np.unique(
            np.concatenate([np.linspace(0, 24, 481), np.linspace(24, tend, 625)])
        ).tolist()
    )
    #: the concentrations, in the first day (linear) and over the 14 days
    #: (logarithmic), with the lower end of their logarithmic axis [nM]: a
    #: decade below the concentration at 14 days after 2 mg oral, which keeps
    #: the bands of the draws of the fitted parameters on the axis
    concentrations: ClassVar[dict[str, float]] = {
        "[Cveblood_rap]": 1e-3,
        "[Cveplasma_rap]": 1e-5,
        "[Cve_rx]": 1e-4,
    }
    #: the excreted amounts over the 14 days
    amounts: ClassVar[tuple[str, ...]] = ("Aurine_rx", "Afeces_rap", "Afeces_rx")
    #: the end of the timecourses of the first day [hr]
    tend_day = 24
    selections = ("time", *concentrations, *amounts)

    @override
    def simulations(self) -> dict[str, Simulation]:
        return {
            self.to_sid(f"rap_{route}_{dose}"): Simulation(
                time_unit="hr",
                end=self.tend,
                times=self.times,
                preinit_changes={
                    **self.default_changes(),
                    f"{route}DOSE_rap": Q(dose, "mg"),
                },
            )
            for route in self.routes
            for dose in self.doses
        }

    @override
    def figures(self) -> dict[str, Figure]:
        figures: dict[str, Figure] = {}
        for route, route_label in self.routes.items():
            fig = Figure(
                experiment=self,
                sid=f"Fig_dose_dependency_rap_{route}",
                num_rows=3,
                num_cols=3,
                name=f"Dose dependency rapamycin {route_label}",
                # the two-line y labels keep off the panel to their left
                settings={"legend_scope": "row", "fig_subplots_wspace": 0.4},
            )
            plots = fig.create_plots(legend=True)
            panels: list[tuple[str, float, float | None]] = [
                *((sid, self.tend_day, None) for sid in self.concentrations),
                *((sid, self.tend, ymin) for sid, ymin in self.concentrations.items()),
                *((sid, self.tend, None) for sid in self.amounts),
            ]
            for plot, (sid, tend, ymin) in zip(plots, panels, strict=True):
                plot.set_xaxis(self.label_time, unit=self.unit_time, min=0, max=tend)
                if ymin is None:
                    # no negative ticks below the curves, e.g. the zero
                    # rapamycin in the feces after an intravenous dose
                    plot.set_yaxis(self.labels[sid], unit=self.units[sid], min=0)
                else:
                    plot.set_yaxis(
                        self.labels[sid], unit=self.units[sid], scale="log", min=ymin
                    )
                for dose in self.doses:
                    plot.add_data(
                        task=self.to_sid(f"task_rap_{route}_{dose}"),
                        xid="time",
                        yid=sid,
                        label=f"{dose:g} mg {route}",
                        color=self.colors[dose],
                    )

            figures |= self.figures_by_sid(fig)

        return figures


if __name__ == "__main__":
    MODEL.run_experiments(
        DoseDependencyExperiment, output_dir=DoseDependencyExperiment.__name__
    )
