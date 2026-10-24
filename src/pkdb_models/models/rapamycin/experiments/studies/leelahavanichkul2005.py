"""Simulation experiment of the study Leelahavanichkul2005."""

from typing import ClassVar, override

import matplotlib.pyplot as plt
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


class Leelahavanichkul2005(RapamycinSimulationExperiment):
    """Simulation experiment of Leelahavanichkul2005.

    Healthy Thai subjects (n = 12, 64.73 kg, 160 cm) took a single dose of 6 mg
    of the oral solution of sirolimus with 240 ml orange juice (the vehicle of
    the label, no interaction) after an overnight fast, fasting until 4 h after
    the dose. The mean and SD per time point (Table 3) are mapped, the
    individual curves of the same subjects (Fig1) are plotted.
    """

    interventions: ClassVar[list[str]] = ["RAP6"]
    colors: ClassVar[dict[str, str]] = {
        "RAP6": "black",
    }
    subjects: ClassVar[list[str]] = [f"S{k + 1}" for k in range(12)]

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig1", "Tab3A"]:
            df = self.load_dataframe(fig_id)
            for label, df_label in df.groupby("label"):
                label = str(label)
                dset = DataSet.from_df(df_label, self.ureg)
                if label.startswith("rapamycin_"):
                    if fig_id == "Fig1":
                        dset.unit_conversion("value", 1 / self.Mr.rap)
                    elif fig_id == "Tab3A":
                        dset.unit_conversion("mean", 1 / self.Mr.rap)

                dsets[label] = dset

        return dsets

    @override
    def simulations(self) -> dict[str, Simulation]:
        simulations = {}
        for intervention in self.interventions:
            simulations[f"rap_{intervention}"] = Simulation(
                time_unit="hr",
                end=30,
                preinit_changes={
                    **self.default_changes(),
                    # physiological changes
                    "BW": Q(64.73, "kg"),
                    "HEIGHT": Q(160, "cm"),
                    # dose
                    "PODOSE_rap": Q(6, "mg"),
                    "GU__fed": Q(self.fasting_map["fasted"], "dimensionless"),
                },
            )

        return simulations

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:

        mappings = {}
        for intervention in self.interventions:
            mappings[f"fm_rap_{intervention}"] = FitMapping(
                self,
                reference=FitData(
                    self,
                    dataset=f"rapamycin_{intervention}",
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                ),
                observable=FitData(
                    self,
                    task=f"task_rap_{intervention}",
                    xid="time",
                    yid="[Cveblood_rap]",
                ),
                metadata=PKPDMappingMetaData(
                    tissue=Tissue.BLOOD,
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
        figures = {}

        cmap = plt.get_cmap("copper")
        subject_colors = {
            subj: cmap(i / (len(self.subjects) - 1))
            for i, subj in enumerate(self.subjects)
        }

        # individual data
        fig = Figure(
            experiment=self,
            sid="Fig1",
            name=f"{self.__class__.__name__} (individuals)",
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )
        plots[0].set_yaxis(self.label_rap_blood, unit=self.unit_rap, scale="linear")
        for intervention in self.interventions:
            # simulation, once for all subjects
            plots[0].add_data(
                task=f"task_rap_{intervention}",
                xid="time",
                yid="[Cveblood_rap]",
                label="6 mg",
                color=self.colors[intervention],
            )
            # data of every subject
            for subject in self.subjects:
                plots[0].add_data(
                    dataset=f"rapamycin_{intervention}_{subject}",
                    xid="time",
                    yid="value",
                    count="count",
                    yid_sd=None,
                    label=subject,
                    color=subject_colors[subject],
                )
        figures |= self.figures_by_sid(fig)

        # mean data
        fig3a = Figure(
            experiment=self,
            sid="Tab3A",
            name=f"{self.__class__.__name__} (mean)",
        )
        plots3a = fig3a.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time),
            legend=True,
        )

        plots3a[0].set_yaxis(self.label_rap_blood, unit=self.unit_rap)

        for intervention in self.interventions:
            # simulation
            plots3a[0].add_data(
                task=f"task_rap_{intervention}",
                xid="time",
                yid="[Cveblood_rap]",
                label="6 mg",
                color=self.colors[intervention],
            )
            plots3a[0].add_data(
                dataset=f"rapamycin_{intervention}",
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label="6 mg",
                color=self.colors[intervention],
            )

        figures |= self.figures_by_sid(fig3a)

        return figures


if __name__ == "__main__":
    MODEL.run_experiments(
        Leelahavanichkul2005, output_dir=Leelahavanichkul2005.__name__
    )
