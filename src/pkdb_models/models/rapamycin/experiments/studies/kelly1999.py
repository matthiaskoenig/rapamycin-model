"""Simulation experiment of the study Kelly1999."""

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


class Kelly1999(RapamycinSimulationExperiment):
    """Simulation experiment of Kelly1999.

    Stable renal transplant recipients (n = 21, 75.8 kg) on the oral solution of
    sirolimus once daily for more than a year (mean 3.6 mg/day) and on
    cyclosporine (Neoral, twice daily, the morning dose taken together with
    sirolimus) were converted mg to mg to the 1 mg tablet. The 12 h profiles at
    steady state after the last dose of the solution and after 2 (n = 12), 4
    and 8 weeks of tablets are shown (Fig1, mean and SD; the liquid profile is
    the same in the three panels of the figure, the panel of week 4 is
    labelled week 6). Fasting is not reported.

    The doses are simulated once daily for 8 weeks, i.e., at steady state (more
    than 5 half-lives); the profile is the one of the last dose (time 0). The
    model has no difference between the solution and the tablet, both are
    simulated alike.
    """

    #: the liquid profile and the tablet profiles at 2, 4 and 8 weeks
    interventions: ClassVar[list[str]] = [
        "RAP_L, CSA",
        "RAP_T2, CSA",
        "RAP_T4, CSA",
        "RAP_T8, CSA",
    ]
    formulations: ClassVar[dict[str, str]] = {
        "RAP_L, CSA": "L",
        "RAP_T2, CSA": "T",
        "RAP_T4, CSA": "T",
        "RAP_T8, CSA": "T",
    }
    colors: ClassVar[dict[str, str]] = {
        "L": "black",
        "T": "#D2691E",
    }
    labels_formulation: ClassVar[dict[str, str]] = {
        "L": "3.6 mg solution",
        "T": "3.6 mg tablet",
    }
    #: the panels of the figure: the liquid profile and the tablet of the week
    panels: ClassVar[dict[str, str]] = {
        "week 2": "RAP_T2, CSA",
        "week 4": "RAP_T4, CSA",
        "week 8": "RAP_T8, CSA",
    }
    dose = 3.6  # [mg/day]
    n_days = 56  # once daily doses until steady state

    #: the residual activity of CYP3A4, CYP3A5 and P-glycoprotein of the gut wall
    #: under cyclosporine. The morning dose of cyclosporine (about 120 mg of the
    #: 240.5 mg/day) was taken together with sirolimus: the gut wall is exposed
    #: to its luminal concentration, far above the inhibitory constants of
    #: cyclosporine for CYP3A4 and P-glycoprotein, as with 300 mg cyclosporine
    #: taken together with sirolimus in Zimmerman2003 (AUC ratio 3.30, half-life
    #: unchanged), a residual activity of 10 %. The liver is not inhibited. The
    #: value is not fitted.
    f_gut_cyclosporine: ClassVar[float] = 0.1

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
        # once daily doses, the data is reported from the last dose
        dose_times = [24 * k for k in range(self.n_days)]
        end = 24 * (self.n_days - 1) + 14
        for formulation in sorted(set(self.formulations.values())):
            simulations[f"rap_{formulation}"] = Simulation(
                time_unit="hr",
                end=end,
                # every 10 min
                time_shift=-24 * (self.n_days - 1),
                preinit_changes={
                    **self.default_changes(),
                    "BW": Q(75.8, "kg"),
                    "GU__fed": Q(self.fasting_map["not reported"], "dimensionless"),
                    # gut CYP3A4, CYP3A5 (genotype not reported, reference *3/*3)
                    # and P-glycoprotein (inhibition by cyclosporine)
                    "GU__f_cyp3a4": Q(
                        self.cyp3a4_activity["*1/*1"] * self.f_gut_cyclosporine,
                        "dimensionless",
                    ),
                    "GU__f_cyp3a5": Q(
                        self.cyp3a5_activity["*3/*3"] * self.f_gut_cyclosporine,
                        "dimensionless",
                    ),
                    "GU__f_pg": Q(self.f_gut_cyclosporine, "dimensionless"),
                },
                changes=[
                    self.add_doses(dose_times, {"PODOSE_rap": Q(self.dose, "mg")})
                ],
            )
        return simulations

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:

        mappings = {}

        for intervention in self.interventions:
            formulation = self.formulations[intervention]
            mappings[self.to_sid(f"fm_rap_{intervention}")] = FitMapping(
                self,
                reference=FitData(
                    self,
                    dataset=self.to_sid(f"rapamycin_{intervention}"),
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                ),
                observable=FitData(
                    self,
                    task=f"task_rap_{formulation}",
                    xid="time",
                    yid="[Cveblood_rap]",
                ),
                metadata=PKPDMappingMetaData(
                    tissue=Tissue.BLOOD,
                    route=Route.PO,
                    application_form=ApplicationForm.SOLUTION
                    if formulation == "L"
                    else ApplicationForm.TABLET,
                    dosing=Dosing.MULTIPLE,
                    health=Health.RENAL_TRANSPLANT,
                    fasting=Fasting.NR,
                    coadministration=Coadministration.CYCLOSPORINE,
                ),
            )
        return mappings

    @override
    def figures(self) -> dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig1",
            num_cols=len(self.panels),
            name=f"{self.__class__.__name__}",
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time, min=-1, max=13),
            legend=True,
        )
        for kp, (title, tablet) in enumerate(self.panels.items()):
            plots[kp].set_title(title)
            plots[kp].set_yaxis(self.label_rap_blood, unit=self.unit_rap)
            for intervention in ["RAP_L, CSA", tablet]:
                formulation = self.formulations[intervention]
                # simulation (the same for both formulations, the tablet dashed)
                plots[kp].add_data(
                    task=f"task_rap_{formulation}",
                    xid="time",
                    yid="[Cveblood_rap]",
                    label=self.labels_formulation[formulation],
                    color=self.colors[formulation],
                    linestyle="solid" if formulation == "L" else "dashed",
                )
                # data
                plots[kp].add_data(
                    dataset=self.to_sid(f"rapamycin_{intervention}"),
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                    label=self.labels_formulation[formulation],
                    color=self.colors[formulation],
                )

        return self.figures_by_sid(fig)


if __name__ == "__main__":
    MODEL.run_experiments(Kelly1999, output_dir=Kelly1999.__name__)
