"""Simulation experiment of the study Brattstrom2000."""

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


class Brattstrom2000(RapamycinSimulationExperiment):
    """Simulation experiment of Brattstrom2000.

    Healthy white men (n = 6 per dose) took a single dose of 0.3, 1, 3, 5 or
    8 mg/m^2 of the oral solution of sirolimus after an overnight fast, fasting
    until 4 h after the dose (Fig1, whole blood, mean and SE; the SD of the
    data is SE * sqrt(6)). The doses are per body surface area, the absolute
    dose is the dose times the mean body surface area of the group (`bsas`).
    """

    interventions: ClassVar[list[str]] = ["RAP03", "RAP1", "RAP3", "RAP5", "RAP8"]
    colors: ClassVar[dict[str, str]] = {
        "RAP03": "black",
        "RAP1": "#FFE299",
        "RAP3": "#FFB347",
        "RAP5": "#FF944D",
        "RAP8": "#D2691E",
    }
    doses: ClassVar[dict[str, float]] = {  # [mg/m^2]
        "RAP03": 0.3,
        "RAP1": 1,
        "RAP3": 3,
        "RAP5": 5,
        "RAP8": 8,
    }
    bsas: ClassVar[dict[str, float]] = {  # [m^2]
        "RAP03": 1.95,
        "RAP1": 2.06,
        "RAP3": 1.94,
        "RAP5": 1.95,
        "RAP8": 1.98,
    }
    bodyweights: ClassVar[dict[str, float]] = {  # [kg], Table 1
        "RAP03": 74.2,
        "RAP1": 82.3,
        "RAP3": 74.5,
        "RAP5": 75.7,
        "RAP8": 74.8,
    }

    def dose_label(self, intervention: str) -> str:
        """Get the label of a dose, per body surface area and absolute."""
        dose = self.doses[intervention]
        return f"{dose:g} mg/m² ({dose * self.bsas[intervention]:.2g} mg)"

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig1"]:
            df = self.load_dataframe(fig_id)
            for label, df_label in df.groupby("label"):
                label = str(label)
                dset = DataSet.from_df(df_label, self.ureg)

                # unit conversion
                if label.startswith("rapamycin_"):
                    dset.unit_conversion("mean", 1 / self.Mr.rap)

                dsets[label] = dset

        return dsets

    @override
    def simulations(self) -> dict[str, Simulation]:
        simulations = {}
        for intervention in self.interventions:
            simulations[f"rap_{intervention}"] = Simulation(
                time_unit="hr",
                end=350,
                preinit_changes={
                    **self.default_changes(),
                    # physiological changes
                    "BW": Q(self.bodyweights[intervention], "kg"),
                    # dose
                    "PODOSE_rap": Q(self.doses[intervention], "mg/m^2")
                    * Q(self.bsas[intervention], "m^2"),
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

        fig = Figure(
            experiment=self,
            sid="Fig1",
            name=f"{self.__class__.__name__}",
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )
        plots[0].set_yaxis(
            self.label_rap_blood, unit=self.unit_rap, scale="log", min=0.2
        )

        fig_12 = Figure(
            experiment=self,
            sid="Fig1_12h",
            name=f"{self.__class__.__name__} (12h)",
        )
        plots_12 = fig_12.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )
        plots_12[0].set_yaxis(self.label_rap_blood, unit=self.unit_rap)
        plots_12[0].set_xaxis(self.label_time, unit=self.unit_time, min=0, max=12)

        for intervention in self.interventions:
            # simulation
            plots[0].add_data(
                task=f"task_rap_{intervention}",
                xid="time",
                yid="[Cveblood_rap]",
                label=self.dose_label(intervention),
                color=self.colors[intervention],
            )
            # data
            plots[0].add_data(
                dataset=f"rapamycin_{intervention}",
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label=self.dose_label(intervention),
                color=self.colors[intervention],
            )

            plots_12[0].add_data(
                task=f"task_rap_{intervention}",
                xid="time",
                yid="[Cveblood_rap]",
                label=self.dose_label(intervention),
                color=self.colors[intervention],
            )
            plots_12[0].add_data(
                dataset=f"rapamycin_{intervention}",
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label=self.dose_label(intervention),
                color=self.colors[intervention],
            )

        return self.figures_by_sid(fig, fig_12)


if __name__ == "__main__":
    MODEL.run_experiments(Brattstrom2000, output_dir=Brattstrom2000.__name__)
