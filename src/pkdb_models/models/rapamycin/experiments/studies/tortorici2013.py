"""Simulation experiment of the study Tortorici2013."""

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


class Tortorici2013(RapamycinSimulationExperiment):
    """Simulation experiment of Tortorici2013 (rapamycin only).

    Four way crossover in healthy fasted men: rapamycin 15 mg alone, together
    with tacrolimus 10 mg, and 4 h after tacrolimus 10 mg (Fig1, mean and SE,
    n = 27). Tacrolimus does not change the pharmacokinetics of rapamycin
    (geometric mean ratio of AUCinf 97 % and 107 %, of Cmax 95 % and 120 %),
    so the arms with tacrolimus are simulated without an interaction.
    """

    dose: ClassVar[float] = 15.0
    interventions: ClassVar[list[str]] = ["RAP15", "TAC10", "TAC10_4H"]
    dataset_ids: ClassVar[dict[str, str]] = {
        "RAP15": "rapamycin_RAP15",
        "TAC10": "rapamycin_RAP15_TAC10",
        "TAC10_4H": "rapamycin_RAP15_TAC10_4H",
    }
    names: ClassVar[dict[str, str]] = {
        "RAP15": "15 mg",
        "TAC10": "15 mg + tacrolimus 10 mg",
        "TAC10_4H": "15 mg 4 h after tacrolimus 10 mg",
    }
    colors: ClassVar[dict[str, str]] = {
        "RAP15": "black",
        "TAC10": "#e3256b",
        "TAC10_4H": "#d65282",
    }

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig1"]:
            df = self.load_dataframe(fig_id)
            for label, df_label in df.groupby("label"):
                label = str(label)
                if not label.startswith("rapamycin_"):
                    continue
                dset = DataSet.from_df(df_label, self.ureg)

                dset.unit_conversion("mean", 1 / self.Mr.rap)
                dsets[label] = dset

        return dsets

    @override
    def simulations(self) -> dict[str, Simulation]:
        simulations = {}
        for intervention in self.interventions:
            simulations[f"rap_{intervention}"] = Simulation(
                time_unit="hr",
                end=150,
                preinit_changes={
                    **self.default_changes(),
                    "BW": Q(71.4, "kg"),
                    "PODOSE_rap": Q(self.dose, "mg"),
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
                    dataset=self.dataset_ids[intervention],
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd" if intervention == "RAP15" else None,
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
                    application_form=ApplicationForm.NR,
                    dosing=Dosing.SINGLE,
                    health=Health.HEALTHY,
                    fasting=Fasting.FASTED,
                    # tacrolimus without an effect on rapamycin
                    coadministration=Coadministration.NONE
                    if intervention == "RAP15"
                    else Coadministration.TACROLIMUS,
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

        for intervention in self.interventions:
            # simulation
            plots[0].add_data(
                task=f"task_rap_{intervention}",
                xid="time",
                yid="[Cveblood_rap]",
                label=self.names[intervention],
                color=self.colors[intervention],
            )
            # data
            plots[0].add_data(
                dataset=self.dataset_ids[intervention],
                xid="time",
                yid="mean",
                yid_sd="mean_sd" if intervention == "RAP15" else None,
                count="count",
                label=self.names[intervention],
                color=self.colors[intervention],
            )

        return self.figures_by_sid(fig)


if __name__ == "__main__":
    MODEL.run_experiments(Tortorici2013, output_dir=Tortorici2013.__name__)
