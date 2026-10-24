"""Simulation experiment of the study Renders2007."""

from typing import ClassVar, override

from pkpdsim.metadata import (
    ApplicationForm,
    Coadministration,
    Dosing,
    Fasting,
    Genotype,
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


class Renders2007(RapamycinSimulationExperiment):
    """Simulation experiment of Renders2007.

    Steady state profiles of rapamycin in whole blood of stable kidney
    transplant recipients by CYP3A5 genotype (Fig4, 7 *3/*3 and 3 *1/*3, no
    *1/*1 in the profile cohort). The assay is described for plasma, but the
    concentrations are whole blood levels. The patients received
    mycophenolate mofetil and/or prednisolone, some verapamil or diltiazem
    (CYP3A inhibitors), which is not reported per patient and not simulated.
    """

    groups: ClassVar[list[str]] = ["CYP3A5_1_3", "CYP3A5_3_3"]
    colors: ClassVar[dict[str, str]] = {
        "CYP3A5_1_3": "#a50f15",
        "CYP3A5_3_3": "#fb6a4a",
    }
    genotypes: ClassVar[dict[str, Genotype]] = {
        "CYP3A5_1_3": Genotype.CYP3A5_1_3,
        "CYP3A5_3_3": Genotype.CYP3A5_3_3,
    }
    f_cyp3a5: ClassVar[dict[str, float]] = {
        "CYP3A5_1_3": RapamycinSimulationExperiment.cyp3a5_activity["*1/*3"],
        "CYP3A5_3_3": RapamycinSimulationExperiment.cyp3a5_activity["*3/*3"],
    }
    #: legend names of the groups
    legends: ClassVar[dict[str, str]] = {
        "CYP3A5_1_3": "CYP3A5 *1/*3",
        "CYP3A5_3_3": "CYP3A5 *3/*3",
    }
    #: the once daily dose of the profile cohort is not reported per genotype
    #: (2.9 +- 2.26 mg overall; 6.5 and 2.31 mg are the doses of the trough
    #: cohort of 20 patients); it is estimated as clearance times AUC0-24 of
    #: Tab3, read as L/h (13.6 L/h * 163.4 ng*h/ml = 2.2 mg for *3/*3,
    #: 17.92 L/h * 95.2 ng*h/ml = 1.7 mg for *1/*3)
    doses: ClassVar[dict[str, float]] = {
        "CYP3A5_1_3": 1.7,
        "CYP3A5_3_3": 2.2,
    }

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig4"]:
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
        # 14 doses, every 24 hours, the data is reported from the last dose
        dose_times = [24 * k for k in range(14)]
        for group in self.groups:
            simulations[f"rap_{group}"] = Simulation(
                time_unit="hr",
                end=24 * 13 + 25,
                # every 6 min
                time_shift=-24 * 13,
                preinit_changes={
                    **self.default_changes(),
                    "BW": Q(70, "kg"),  # no bodyweight in article
                    # CYP3A5 activity
                    "GU__f_cyp3a5": Q(self.f_cyp3a5[group], "dimensionless"),
                    "LI__f_cyp3a5": Q(self.f_cyp3a5[group], "dimensionless"),
                },
                changes=[
                    self.add_doses(
                        dose_times, {"PODOSE_rap": Q(self.doses[group], "mg")}
                    ),
                ],
            )

        return simulations

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:

        mappings = {}

        for group in self.groups:
            mappings[f"fm_rap_{group}"] = FitMapping(
                self,
                reference=FitData(
                    self,
                    dataset=f"rapamycin_{group}",
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                ),
                observable=FitData(
                    self,
                    task=f"task_rap_{group}",
                    xid="time",
                    yid="[Cveblood_rap]",
                ),
                metadata=PKPDMappingMetaData(
                    tissue=Tissue.BLOOD,
                    route=Route.PO,
                    application_form=ApplicationForm.NR,
                    dosing=Dosing.MULTIPLE,
                    health=Health.RENAL_TRANSPLANT,
                    fasting=Fasting.NR,
                    # mycophenolate mofetil, prednisolone, partly verapamil or
                    # diltiazem
                    coadministration=Coadministration.COMEDICATION,
                    genotype=self.genotypes[group],
                ),
            )

        return mappings

    @override
    def figures(self) -> dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig4",
            name=f"{self.__class__.__name__}",
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time, min=-5, max=25),
            legend=True,
        )
        plots[0].set_yaxis(self.label_rap_blood, unit=self.unit_rap, scale="linear")

        for group in self.groups:
            # simulation
            plots[0].add_data(
                task=f"task_rap_{group}",
                xid="time",
                yid="[Cveblood_rap]",
                label=f"{self.legends[group]} ({self.doses[group]:g} mg)",
                color=self.colors[group],
            )
            # data
            plots[0].add_data(
                dataset=f"rapamycin_{group}",
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label=f"{self.legends[group]} ({self.doses[group]:g} mg)",
                color=self.colors[group],
            )

        return self.figures_by_sid(fig)


if __name__ == "__main__":
    MODEL.run_experiments(Renders2007, output_dir=Renders2007.__name__)
