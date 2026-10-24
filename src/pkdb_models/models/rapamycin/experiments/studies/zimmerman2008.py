"""Simulation experiment of the study Zimmerman2008."""

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


class Zimmerman2008(RapamycinSimulationExperiment):
    """Simulation experiment of Zimmerman2008.

    Nine subjects with hepatic impairment of Child-Pugh grade C and nine
    matched healthy controls took 15 mg of the oral solution after an overnight
    fast of at least 10 h. Fig1 shows the mean and the standard error (n=9,
    "Mean ± SE"). Child-Pugh C is simulated as severe cirrhosis.
    """

    interventions: ClassVar[list[str]] = ["RAP15"]
    groups: ClassVar[list[str]] = [
        "CONTROL",
        "GRADE_C",
    ]
    colors: ClassVar[dict[str, str]] = {
        "CONTROL": "black",
        "GRADE_C": RapamycinSimulationExperiment.cirrhosis_colors["Severe cirrhosis"],
    }
    labels_group: ClassVar[dict[str, str]] = {
        "CONTROL": "15 mg healthy",
        "GRADE_C": "15 mg Child-Pugh C",
    }
    bodyweights: ClassVar[dict[str, float]] = {  # [kg]
        "CONTROL": 64.3,
        "GRADE_C": 69.2,
    }
    heights: ClassVar[dict[str, float]] = {  # [cm]
        "CONTROL": 167,
        "GRADE_C": 169,
    }
    f_cirrhosis: ClassVar[dict[str, float]] = {
        "CONTROL": RapamycinSimulationExperiment.cirrhosis_map["Control"],
        "GRADE_C": RapamycinSimulationExperiment.cirrhosis_map["Severe cirrhosis"],
    }

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

                dsets[self.to_sid(label)] = dset

        return dsets

    @override
    def simulations(self) -> dict[str, Simulation]:
        simulations = {}
        for group in self.groups:
            simulations[f"rap_RAP15_{group}"] = Simulation(
                time_unit="hr",
                end=320,
                # every 6 min, the absorption peak is within the first hour
                preinit_changes={
                    **self.default_changes(),
                    # physiological changes
                    "BW": Q(self.bodyweights[group], "kg"),
                    "HEIGHT": Q(self.heights[group], "cm"),
                    # dose
                    "PODOSE_rap": Q(15, "mg"),
                    # fasted
                    "GU__fed": Q(self.fasting_map["fasted"], "dimensionless"),
                    # hepatic impairment
                    "f_cirrhosis": Q(self.f_cirrhosis[group], "dimensionless"),
                },
            )

        return simulations

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:

        mappings = {}

        for group in self.groups:
            mappings[f"fm_rap_RAP15_{group}"] = FitMapping(
                self,
                reference=FitData(
                    self,
                    dataset=f"rapamycin_RAP15_{group}",
                    xid="time",
                    yid="mean",
                    yid_se="mean_se",
                    count="count",
                ),
                observable=FitData(
                    self,
                    task=f"task_rap_RAP15_{group}",
                    xid="time",
                    yid="[Cveblood_rap]",
                ),
                metadata=PKPDMappingMetaData(
                    tissue=Tissue.BLOOD,
                    route=Route.PO,
                    application_form=ApplicationForm.SOLUTION,
                    dosing=Dosing.SINGLE,
                    health=Health.HEPATIC_IMPAIRMENT
                    if group == "GRADE_C"
                    else Health.HEALTHY,
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
        plots[0].set_yaxis(self.label_rap_blood, unit=self.unit_rap, scale="linear")

        for group in self.groups:
            # simulation
            plots[0].add_data(
                task=f"task_rap_RAP15_{group}",
                xid="time",
                yid="[Cveblood_rap]",
                label=self.labels_group[group],
                color=self.colors[group],
            )
            # data
            plots[0].add_data(
                dataset=f"rapamycin_RAP15_{group}",
                xid="time",
                yid="mean",
                yid_se="mean_se",
                count="count",
                label=self.labels_group[group],
                color=self.colors[group],
            )

        return self.figures_by_sid(fig)


if __name__ == "__main__":
    MODEL.run_experiments(Zimmerman2008, output_dir=Zimmerman2008.__name__)
