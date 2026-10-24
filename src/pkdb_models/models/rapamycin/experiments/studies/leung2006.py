"""Simulation experiment of the study Leung2006."""

from typing import ClassVar, override

from pkpdsim.metadata import (
    ApplicationForm,
    Coadministration,
    Dosing,
    Fasting,
    Health,
    PKPDMappingMetaData,
    Route,
    Tissue,
)
from sbmlsim import Q
from sbmlsim.data import DataSet
from sbmlsim.fit import FitData, FitMapping
from sbmlsim.plot import Axis, Figure
from sbmlsim.simulation import Simulation

from pkdb_models.models.rapamycin import MODEL
from pkdb_models.models.rapamycin.experiments.base_experiment import (
    RapamycinSimulationExperiment,
)


class Leung2006(RapamycinSimulationExperiment):
    """Simulation experiment of Leung2006.

    Mass balance of a single oral dose of [14C]rapamycin (oral solution) in six
    healthy fasted men. Fig2 is rapamycin in whole blood and plasma (LC-MS/MS),
    Fig1 the total radioactivity in whole blood and plasma (rapamycin and its
    metabolites in rapamycin equivalents). The excretion of the radioactivity
    (91.0 % of the dose in the feces, 2.2 % in the urine over 15 days) is
    reported in the text only and not curated.
    """

    #: the nominal dose was 40 mg, the actual dose 42 mg (Methods, Tab3); the
    #: time courses are curated with the intervention RAP40
    dose: ClassVar[float] = 42.0

    tissues: ClassVar[list[str]] = ["blood", "plasma"]
    #: observables of rapamycin (Fig2) and of the total radioactivity (Fig1)
    sids: ClassVar[dict[str, dict[str, str]]] = {
        "rap": {
            "blood": "[Cveblood_rap]",
            "plasma": "[Cveplasma_rap]",
        },
        "raptot": {
            "blood": "[Cveblood_raptot]",
            "plasma": "[Cveplasma_raptot]",
        },
    }
    fig_ids: ClassVar[dict[str, str]] = {
        "rap": "Fig2",
        "raptot": "Fig1",
    }
    colors: ClassVar[dict[str, str]] = {
        "blood": "tab:red",
        "plasma": "tab:blue",
    }

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        for substance, fig_id in self.fig_ids.items():
            df = self.load_dataframe(fig_id)
            for tissue in self.tissues:
                df_sub = df[
                    (df["measurement_type"] == "concentration")
                    & (df["tissue"] == tissue)
                ]
                dset = DataSet.from_df(df_sub.copy(), self.ureg)
                # the total radioactivity is in rapamycin equivalents
                dset.unit_conversion("mean", 1 / self.Mr.rap)
                dsets[f"{fig_id}_{substance}_{tissue}"] = dset

        return dsets

    @override
    def simulations(self) -> dict[str, Simulation]:
        return {
            "rap42": Simulation(
                time_unit="hr",
                end=150,
                preinit_changes={
                    **self.default_changes(),
                    # physiological changes
                    "BW": Q(73, "kg"),  # six men
                    # dose
                    "PODOSE_rap": Q(self.dose, "mg"),
                },
            )
        }

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:
        mappings = {}
        for substance, fig_id in self.fig_ids.items():
            for tissue in self.tissues:
                mappings[f"fm_{substance}_{fig_id}_{tissue}"] = FitMapping(
                    self,
                    reference=FitData(
                        self,
                        dataset=f"{fig_id}_{substance}_{tissue}",
                        xid="time",
                        yid="mean",
                        yid_sd="mean_sd",
                        count="count",
                    ),
                    observable=FitData(
                        self,
                        task="task_rap42",
                        xid="time",
                        yid=self.sids[substance][tissue],
                    ),
                    metadata=PKPDMappingMetaData(
                        tissue=Tissue.BLOOD if tissue == "blood" else Tissue.PLASMA,
                        route=Route.PO,
                        application_form=ApplicationForm.SOLUTION,
                        dosing=Dosing.SINGLE,
                        health=Health.HEALTHY,
                        fasting=Fasting.FASTED,
                        coadministration=Coadministration.NONE,
                    ),
                )

        return mappings

    @override
    def figures(self) -> dict[str, Figure]:
        labels = {
            "rap": self.label_rap,
            "raptot": self.label_raptot,
        }
        figures = []
        for substance, fig_id in self.fig_ids.items():
            fig = Figure(
                experiment=self,
                sid=fig_id,
                num_cols=2,
                name=f"{self.__class__.__name__} ({fig_id})",
            )
            plots = fig.create_plots(
                xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
            )
            for k, tissue in enumerate(self.tissues):
                plots[k].set_yaxis(f"{labels[substance]} {tissue}", unit=self.unit_rap)
                # simulation
                plots[k].add_data(
                    task="task_rap42",
                    xid="time",
                    yid=self.sids[substance][tissue],
                    label=f"{self.dose:g} mg",
                    color=self.colors[tissue],
                )
                # data
                plots[k].add_data(
                    dataset=f"{fig_id}_{substance}_{tissue}",
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                    label=f"{self.dose:g} mg",
                    color=self.colors[tissue],
                )
            figures.append(fig)

        return self.figures_by_sid(*figures)


if __name__ == "__main__":
    MODEL.run_experiments(Leung2006, output_dir=Leung2006.__name__)
