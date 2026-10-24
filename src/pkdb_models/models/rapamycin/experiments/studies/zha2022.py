"""Simulation experiment of the study Zha2022."""

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


class Zha2022(RapamycinSimulationExperiment):
    """Simulation experiment of Zha2022.

    Healthy subjects (sirolimus arm: 11 with pharmacokinetics, 77.2 kg) took a
    single oral dose of 2 mg sirolimus alone (period 1) and of 0.5 mg on day 15
    of the 3D regimen (ombitasvir/paritaprevir/ritonavir 25/150/100 mg once
    daily and dasabuvir 250 mg twice daily, days 1 to 34, period 2). The
    ritonavir of the regimen, a strong mechanism-based inhibitor of CYP3A,
    increased the dose-normalized AUCinf 38-fold and the half-life from 96 to
    249 h (Table 1). All drugs were taken "under non-fasting conditions"
    (Supplement S1) without a reported meal, i.e., the fed state is not
    reported. The formulation of sirolimus is not reported.
    """

    interventions: ClassVar[list[str]] = ["RAP2", "RAP05_3D"]
    colors: ClassVar[dict[str, str]] = {
        "RAP2": "black",
        "RAP05_3D": "#e3256b",
    }
    labels_intervention: ClassVar[dict[str, str]] = {
        "RAP2": "2 mg",
        "RAP05_3D": "0.5 mg + 3D regimen (ritonavir)",
    }
    doses: ClassVar[dict[str, float]] = {
        "RAP2": 2,
        "RAP05_3D": 0.5,
    }
    #: the residual activity of CYP3A4 and CYP3A5 in the gut wall and the liver
    #: under ritonavir. Ritonavir (100 mg once daily, at steady state from 14
    #: days before to 19 days after the dose of sirolimus) is a strong
    #: mechanism-based inactivator of CYP3A4 and CYP3A5 in the gut and the liver,
    #: which inhibits CYP3A by more than 90 % (oral midazolam AUC ratios of 20 to
    #: 30 with 100 mg ritonavir), a residual activity of 5 %. Its effect in the
    #: study, a dose-normalized AUCinf ratio of 38 with a half-life of 96 -> 249 h,
    #: is on the first pass and on the elimination, so the inhibition acts in the
    #: gut and the liver. P-glycoprotein is not changed. The value is not fitted.
    f_cyp3a_ritonavir: ClassVar[float] = 0.05
    f_cyp3a: ClassVar[dict[str, float]] = {
        "RAP2": 1.0,
        "RAP05_3D": f_cyp3a_ritonavir,
    }

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig2"]:
            df = self.load_dataframe(fig_id)
            for label, df_label in df.groupby("label"):
                label = str(label)
                dset = DataSet.from_df(df_label, self.ureg)
                dset.unit_conversion("mean", 1 / self.Mr.rap)
                dsets[self.to_sid(label)] = dset

        return dsets

    @override
    def simulations(self) -> dict[str, Simulation]:
        simulations = {}
        for intervention in self.interventions:
            simulations[f"rap_{intervention}"] = Simulation(
                time_unit="hr",
                end=500,
                preinit_changes={
                    **self.default_changes(),
                    # physiological changes (sirolimus arm, Table S1)
                    "BW": Q(77.2, "kg"),
                    "HEIGHT": Q(174, "cm"),
                    "PODOSE_rap": Q(self.doses[intervention], "mg"),
                    # fed state not reported (non-fasting without a meal)
                    "GU__fed": Q(self.fasting_map["not reported"], "dimensionless"),
                    # CYP3A4 and CYP3A5 activity (inhibition by ritonavir), the
                    # CYP3A5 genotype is not reported (reference *3/*3)
                    **{
                        f"{organ}__f_{enzyme}": Q(
                            activity * self.f_cyp3a[intervention], "dimensionless"
                        )
                        for organ in ["GU", "LI"]
                        for enzyme, activity in [
                            ("cyp3a4", self.cyp3a4_activity["*1/*1"]),
                            ("cyp3a5", self.cyp3a5_activity["*3/*3"]),
                        ]
                    },
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
                    application_form=ApplicationForm.NR,
                    dosing=Dosing.SINGLE,
                    health=Health.HEALTHY,
                    fasting=Fasting.NR,
                    coadministration=Coadministration.NONE
                    if intervention == "RAP2"
                    else Coadministration.RITONAVIR,
                ),
            )

        return mappings

    @override
    def figures(self) -> dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig2",
            name=f"{self.__class__.__name__}",
        )
        plots = fig.create_plots(xaxis=Axis(self.label_time, unit=self.unit_time))
        plots[0].set_yaxis(self.label_rap_blood, unit=self.unit_rap)
        for intervention in self.interventions:
            # simulation
            plots[0].add_data(
                task=f"task_rap_{intervention}",
                xid="time",
                yid="[Cveblood_rap]",
                label=self.labels_intervention[intervention],
                color=self.colors[intervention],
            )
            # data
            plots[0].add_data(
                dataset=f"rapamycin_{intervention}",
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label=self.labels_intervention[intervention],
                color=self.colors[intervention],
            )

        return self.figures_by_sid(fig)


if __name__ == "__main__":
    MODEL.run_experiments(Zha2022, output_dir=Zha2022.__name__)
