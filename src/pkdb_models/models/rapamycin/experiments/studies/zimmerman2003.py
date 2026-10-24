"""Simulation experiment of the study Zimmerman2003."""

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


class Zimmerman2003(RapamycinSimulationExperiment):
    """Simulation experiment of Zimmerman2003.

    Healthy subjects (21 completed, 74.6 kg) took 10 mg of the oral solution of
    sirolimus alone, together with 300 mg cyclosporine (Neoral) or 4 h after
    the cyclosporine (randomized crossover), after an overnight fast and a clear
    liquid breakfast. Cyclosporine increased the AUC of sirolimus 3.30-fold
    when taken together and 1.80-fold when taken 4 h before, the half-life and
    the mean residence time were unchanged (Table I): cyclosporine inhibits the
    first pass of sirolimus (intestinal CYP3A and P-glycoprotein), not its
    elimination.
    """

    interventions: ClassVar[list[str]] = ["RAP10", "RAP10_CSA300", "RAP10_CSA300B"]
    colors: ClassVar[dict[str, str]] = {
        "RAP10": "black",
        "RAP10_CSA300": "#e3256b",
        "RAP10_CSA300B": "#d65282",
    }
    labels_intervention: ClassVar[dict[str, str]] = {
        "RAP10": "10 mg",
        "RAP10_CSA300": "10 mg + cyclosporine 300 mg",
        "RAP10_CSA300B": "10 mg 4 h after cyclosporine 300 mg",
    }
    #: the residual activity of CYP3A4, CYP3A5 and P-glycoprotein of the gut wall
    #: under cyclosporine; the liver is not inhibited (half-life unchanged).
    #: Taken together, the gut wall is exposed to the high luminal concentration
    #: of the 300 mg dose (about 1 mg/mL, far above the inhibitory constants of
    #: cyclosporine for CYP3A4 and P-glycoprotein in the low micromolar range),
    #: a near complete inhibition of the first pass, a residual activity of
    #: 10 %. Taken 4 h before sirolimus, cyclosporine is absorbed (tmax 1.5 h)
    #: and only its blood concentration reaches the gut wall: the AUC ratio of
    #: 1.80 is half of the joint effect (3.30) on a log scale, a residual
    #: activity of 30 %. The values are not fitted.
    f_gut_cyclosporine: ClassVar[dict[str, float]] = {
        "RAP10": 1.0,
        "RAP10_CSA300": 0.1,
        "RAP10_CSA300B": 0.3,
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
                # the predose samples (time 0) are below the limit of
                # quantification
                df_label = df_label[df_label.time > 0]
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
                end=150,
                preinit_changes={
                    **self.default_changes(),
                    "BW": Q(74.6, "kg"),
                    "PODOSE_rap": Q(10, "mg"),
                    # fasted (clear liquid breakfast)
                    "GU__fed": Q(self.fasting_map["fasted"], "dimensionless"),
                    # gut CYP3A4, CYP3A5 (genotype not reported, reference *3/*3)
                    # and P-glycoprotein (inhibition by cyclosporine)
                    "GU__f_cyp3a4": Q(
                        self.cyp3a4_activity["*1/*1"]
                        * self.f_gut_cyclosporine[intervention],
                        "dimensionless",
                    ),
                    "GU__f_cyp3a5": Q(
                        self.cyp3a5_activity["*3/*3"]
                        * self.f_gut_cyclosporine[intervention],
                        "dimensionless",
                    ),
                    "GU__f_pg": Q(
                        self.f_gut_cyclosporine[intervention], "dimensionless"
                    ),
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
                    else Coadministration.CYCLOSPORINE,
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
    MODEL.run_experiments(Zimmerman2003, output_dir=Zimmerman2003.__name__)
