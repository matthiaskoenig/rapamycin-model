"""Simulation experiment of the study Bottiger2001."""

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


class Bottiger2001(RapamycinSimulationExperiment):
    """Simulation experiment of Bottiger2001.

    Healthy subjects (12 men of 82 kg, 6 women of 64 kg) took a single dose of
    10 mg of the oral solution of sirolimus alone or 1 h after a single dose of
    120 mg diltiazem (two 60 mg tablets) in a randomized crossover, fasted from
    the evening before until 4 h after the dose (Fig1, mean and SD). Diltiazem
    increased the AUC of sirolimus 1.60-fold and its Cmax 1.43-fold, the
    half-life did not increase (79 h alone, 67 h with diltiazem, Table I):
    diltiazem inhibits the first pass of sirolimus (intestinal CYP3A and
    P-glycoprotein), not its elimination.
    """

    interventions: ClassVar[list[str]] = ["RAP10", "RAP10_DIL120"]
    colors: ClassVar[dict[str, str]] = {
        "RAP10": "black",
        "RAP10_DIL120": "#e3256b",
    }
    labels_intervention: ClassVar[dict[str, str]] = {
        "RAP10": "10 mg",
        "RAP10_DIL120": "10 mg 1 h after diltiazem 120 mg",
    }
    #: the residual activity of CYP3A4, CYP3A5 and P-glycoprotein of the gut wall
    #: under diltiazem; the liver is not inhibited (half-life not prolonged,
    #: 79 h alone and 67 h with diltiazem). Diltiazem (and N-desmethyldiltiazem)
    #: is a mechanism based inhibitor of CYP3A4 and an inhibitor of
    #: P-glycoprotein; the AUC ratio of sirolimus was 1.60 (90 % CI 1.35-1.90).
    #: Relative to the near complete inhibition of the first pass by
    #: cyclosporine taken together with sirolimus (AUC ratio 3.30, residual
    #: activity 10 %, Zimmerman2003) the effect is 0.39 on a log scale
    #: (ln 1.60 / ln 3.30), a residual activity of 0.1^0.39 = 40 %. The value is
    #: not fitted.
    f_gut_diltiazem: ClassVar[dict[str, float]] = {
        "RAP10": 1.0,
        "RAP10_DIL120": 0.4,
    }

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig1"]:
            df = self.load_dataframe(fig_id)
            for label, df_label in df.groupby("label"):
                label = str(label)
                dset = DataSet.from_df(df_label, self.ureg)

                if label.startswith("rapamycin_"):
                    dset.unit_conversion("mean", 1 / self.Mr.rap)

                dsets[label] = dset

        return dsets

    @override
    def simulations(self) -> dict[str, Simulation]:
        simulations = {}
        for intervention in self.interventions:
            f_gut = self.f_gut_diltiazem[intervention]
            simulations[f"rap_{intervention}"] = Simulation(
                time_unit="hr",
                end=150,
                preinit_changes={
                    **self.default_changes(),
                    # mean of 12 men (82 kg) and 6 women (64 kg)
                    "BW": Q((82 * 12 + 64 * 6) / 18, "kg"),
                    "PODOSE_rap": Q(10, "mg"),
                    "GU__fed": Q(self.fasting_map["fasted"], "dimensionless"),
                    # gut CYP3A4, CYP3A5 (genotype not reported, reference *3/*3)
                    # and P-glycoprotein (inhibition by diltiazem)
                    "GU__f_cyp3a4": Q(
                        self.cyp3a4_activity["*1/*1"] * f_gut, "dimensionless"
                    ),
                    "GU__f_cyp3a5": Q(
                        self.cyp3a5_activity["*3/*3"] * f_gut, "dimensionless"
                    ),
                    "GU__f_pg": Q(f_gut, "dimensionless"),
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
                    coadministration=Coadministration.NONE
                    if intervention == "RAP10"
                    else Coadministration.DILTIAZEM,
                ),
            )

        return mappings

    @override
    def figures(self) -> dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig1",
            num_cols=2,
            name=f"{self.__class__.__name__}",
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )
        plots[0].set_xaxis(self.label_time, unit=self.unit_time, min=0, max=24)
        plots[0].set_yaxis(self.label_rap_blood, unit=self.unit_rap, scale="linear")
        plots[1].set_yaxis(
            self.label_rap_blood, unit=self.unit_rap, scale="log", min=0.5
        )

        for plot in plots:
            for intervention in self.interventions:
                plot.add_data(
                    task=f"task_rap_{intervention}",
                    xid="time",
                    yid="[Cveblood_rap]",
                    label=self.labels_intervention[intervention],
                    color=self.colors[intervention],
                )
                plot.add_data(
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
    MODEL.run_experiments(Bottiger2001, output_dir=Bottiger2001.__name__)
