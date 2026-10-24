"""Simulation experiment of the study Zimmerman1997."""

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


class Zimmerman1997(RapamycinSimulationExperiment):
    """Simulation experiment of Zimmerman1997.

    Stable renal transplant patients (sirolimus group: n=30, 83.2 kg, BSA 1.89
    m^2) on cyclosporine at steady state, prednisone and azathioprine took an
    oral solution of sirolimus per body surface area twice daily for 13 days
    and once on the morning of day 14 (27 doses, Fig1 shows 0.5, 1.5, 2.5, 3.5
    and 6.5 mg/m^2). On the days of the profiles the morning dose was taken
    about 1 h after a light low-fat liquid breakfast, the evening dose 2 to 3 h
    after the evening meal; cyclosporine was taken 3 h after sirolimus. The
    troughs were sampled before the evening dose of the days 2, 4, 6, 9 and 13.
    """

    interventions: ClassVar[list[str]] = [
        "RAP05",
        "RAP1_5",
        "RAP2_5",
        "RAP3_5",
        "RAP6_5",
    ]
    colors: ClassVar[dict[str, str]] = {
        "RAP05": "black",
        "RAP1_5": "#FFE299",
        "RAP2_5": "#FFB347",
        "RAP3_5": "#FF944D",
        "RAP6_5": "#D2691E",
    }
    doses: ClassVar[dict[str, float]] = {  # [mg/m^2]
        "RAP05": 0.5,
        "RAP1_5": 1.5,
        "RAP2_5": 2.5,
        "RAP3_5": 3.5,
        "RAP6_5": 6.5,
    }
    bsa = 1.89  # [m^2]

    #: the residual activity of CYP3A4, CYP3A5 and P-glycoprotein of the gut wall
    #: under cyclosporine. Cyclosporine inhibits the intestinal CYP3A and the
    #: P-glycoprotein efflux, i.e., the first pass of sirolimus, but not its
    #: elimination: the half-life is unchanged with cyclosporine (Zimmerman2003,
    #: t1/2 ratio 0.97 and 1.10). Here cyclosporine was taken at steady state 3 h
    #: after sirolimus, not with it: the inhibition of the first pass is
    #: incomplete, as with cyclosporine 4 h before sirolimus in Zimmerman2003
    #: (AUC ratio 1.80 vs 3.30 when taken together, half of the effect on a log
    #: scale), a residual activity of 30 %. The liver is not inhibited. The value
    #: is not fitted.
    f_gut_cyclosporine: ClassVar[float] = 0.3

    def dose_label(self, intervention: str) -> str:
        """Get the label of a dose, per body surface area and absolute."""
        dose = self.doses[intervention]
        return f"{dose:g} mg/m² ({dose * self.bsa:.2g} mg) bid"

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
        # Doses were administered twice daily for 13 days and once on the morning of study day 14.
        # 27 doses, every 12 hours, the data is reported from the first dose
        dose_times = [12 * k for k in range(27)]
        for intervention in self.interventions:
            # the doses are per body surface area, see `bsa`
            dose = self.doses[intervention]
            simulations[f"rap_{intervention}"] = Simulation(
                time_unit="hr",
                end=12 * 26 + 200,
                # every 3 min
                preinit_changes={
                    **self.default_changes(),
                    # physiological changes
                    "BW": Q(83.2, "kg"),
                    # light breakfast (morning) and evening meal (evening dose)
                    "GU__fed": Q(self.fasting_map["light breakfast"], "dimensionless"),
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
                    # dose
                    self.add_doses(
                        dose_times,
                        {"PODOSE_rap": Q(dose, "mg/m^2") * Q(self.bsa, "m^2")},
                    ),
                ],
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
                    dosing=Dosing.MULTIPLE,
                    health=Health.RENAL_TRANSPLANT,
                    fasting=Fasting.LIGHT_BREAKFAST,
                    coadministration=Coadministration.CYCLOSPORINE,
                ),
            )
        return mappings

    @override
    def figures(self) -> dict[str, Figure]:

        fig = Figure(
            experiment=self,
            sid="Fig1",
            name=f"{self.__class__.__name__}",
            # the legend of ten entries right of the panel, not on the data
            settings={"legend_scope": "row"},
        )
        plots = fig.create_plots(xaxis=Axis(self.label_time, unit="day"), legend=True)
        plots[0].set_yaxis(self.label_rap_blood, unit=self.unit_rap, scale="linear")

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

        return self.figures_by_sid(fig)


if __name__ == "__main__":
    MODEL.run_experiments(Zimmerman1997, output_dir=Zimmerman1997.__name__)
