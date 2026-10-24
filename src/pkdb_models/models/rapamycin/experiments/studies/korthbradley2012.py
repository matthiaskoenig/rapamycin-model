"""Simulation experiment of the study KorthBradley2012."""

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


class KorthBradley2012(RapamycinSimulationExperiment):
    """Simulation experiment of KorthBradley2012.

    Healthy subjects (29 men, 1 woman; 14 black, 16 white; 79 kg) took a single
    dose of 5 mg sirolimus as 5 x 1 mg Rapamune tablets (reference) or as
    10 x 0.5 mg tablets of two prototypes (nonshellac, shellac) in a randomized
    crossover after an overnight fast of at least 10 h (Fig1, 0-144 h, and
    Fig2, 0-6 h, mean and SE). The formulations are bioequivalent and are
    simulated alike. Fig1 is mapped, Fig2 (the early phase on a linear axis)
    is plotted.
    """

    interventions: ClassVar[list[str]] = ["RAP_T5", "RAP_NS5", "RAP_S5"]
    colors: ClassVar[dict[str, str]] = {
        "RAP_T5": "black",
        "RAP_NS5": "#FFB98C",
        "RAP_S5": "#D28F66",
    }
    labels_intervention: ClassVar[dict[str, str]] = {
        "RAP_T5": "5 mg Rapamune tablets",
        "RAP_NS5": "5 mg nonshellac tablets",
        "RAP_S5": "5 mg shellac tablets",
    }

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig1", "Fig2"]:
            df = self.load_dataframe(fig_id)
            for label, df_label in df.groupby("label"):
                label = str(label)
                dset = DataSet.from_df(df_label, self.ureg)

                # unit conversion
                if label.startswith("rapamycin"):
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
                    # physiological changes
                    "BW": Q(79, "kg"),
                    # dose
                    "PODOSE_rap": Q(5, "mg"),
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
                    application_form=ApplicationForm.TABLET,
                    dosing=Dosing.SINGLE,
                    health=Health.HEALTHY,
                    fasting=Fasting.FASTED,
                    coadministration=Coadministration.NONE,
                ),
            )

        return mappings

    @override
    def figures(self) -> dict[str, Figure]:

        fig1 = Figure(
            experiment=self,
            sid="Fig1",
            num_cols=2,
            name=f"{self.__class__.__name__}",
        )
        plots = fig1.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )
        plots[0].set_xaxis(self.label_time, unit=self.unit_time, min=0, max=12)
        plots[0].set_yaxis(self.label_rap_blood, unit=self.unit_rap, scale="linear")
        plots[1].set_yaxis(
            self.label_rap_blood, unit=self.unit_rap, scale="log", min=0.2
        )

        for intervention in self.interventions:
            label = self.labels_intervention[intervention]
            color = self.colors[intervention]
            for plot in plots:
                plot.add_data(
                    task=f"task_rap_{intervention}",
                    xid="time",
                    yid="[Cveblood_rap]",
                    label=label,
                    color=color,
                )
            # data: Fig2 (0-6 h, linear axis) and Fig1 (0-144 h, log axis)
            for plot, dataset in [
                (plots[0], f"rapamycin2_{intervention}"),
                (plots[1], f"rapamycin_{intervention}"),
            ]:
                plot.add_data(
                    dataset=dataset,
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                    label=label,
                    color=color,
                )
        return self.figures_by_sid(fig1)


if __name__ == "__main__":
    MODEL.run_experiments(KorthBradley2012, output_dir=KorthBradley2012.__name__)
