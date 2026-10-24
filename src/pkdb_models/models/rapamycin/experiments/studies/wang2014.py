"""Simulation experiment of the study Wang2014."""

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


class Wang2014(RapamycinSimulationExperiment):
    """Simulation experiment of Wang2014.

    Steady state profiles of rapamycin in whole blood of Chinese renal
    transplant recipients (1 mg tablets once daily, individual doses of 1 to
    4 mg), without (n = 15) and with cyclosporine (n = 9) (Fig1, mean and SE).
    The comedication was prednisone (23 of 24), mycophenolate mofetil (21 of
    24) and partly tacrolimus; the CYP3A5 genotypes are mixed (3 *1/*1, 8
    *1/*3, 13 *3/*3). The dose was taken with or without food as usual.
    Cyclosporine (low doses, mostly 4 to 6 h apart) shows no effect on
    rapamycin (dose normalised Cmax 15.9 vs 18.3 ng/ml, AUCtau 231 vs 246
    ng*h/ml, CL/F 11.1 vs 9.6 L/h), so no interaction is simulated.
    """

    interventions: ClassVar[list[str]] = ["RAP", "RAP_CSA"]
    names: ClassVar[dict[str, str]] = {
        "RAP": "without cyclosporine",
        "RAP_CSA": "with cyclosporine",
    }
    colors: ClassVar[dict[str, str]] = {
        "RAP": "black",
        "RAP_CSA": "#e3256b",
    }
    #: the mean dose [mg] is not reported; the dose normalised (to 2 mg)
    #: values of Tab1 give it as value / normalised value * 2 mg (a ratio of
    #: means): without cyclosporine 1.53 (Cmax), 1.54 (AUCtau) and 1.53 mg
    #: (Ctrough), with cyclosporine 1.80, 1.88 and 1.86 mg
    doses: ClassVar[dict[str, float]] = {
        "RAP": 1.5,
        "RAP_CSA": 1.85,
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

                dsets[label] = dset

        return dsets

    @override
    def simulations(self) -> dict[str, Simulation]:
        simulations = {}
        # 15 doses, every 24 hours, the data is reported from the last dose
        dose_times = [24 * k for k in range(15)]
        for intervention in self.interventions:
            simulations[f"rap_{intervention}"] = Simulation(
                time_unit="hr",
                end=24 * 14 + 25,
                # every 10 min
                time_shift=-24 * 14,
                preinit_changes={
                    **self.default_changes(),
                    # physiological changes
                    "BW": Q(70, "kg"),  # no weight, only bmi
                },
                changes=[
                    # dose
                    self.add_doses(
                        dose_times,
                        {"PODOSE_rap": Q(self.doses[intervention], "mg")},
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
                    application_form=ApplicationForm.TABLET,
                    dosing=Dosing.MULTIPLE,
                    health=Health.RENAL_TRANSPLANT,
                    # with or without food
                    fasting=Fasting.NR,
                    # prednisone, mycophenolate mofetil, partly tacrolimus
                    coadministration=Coadministration.CYCLOSPORINE
                    if intervention == "RAP_CSA"
                    else Coadministration.COMEDICATION,
                ),
            )

        return mappings

    @override
    def figures(self) -> dict[str, Figure]:

        fig1 = Figure(
            experiment=self,
            sid="Fig1",
            name=f"{self.__class__.__name__}",
        )
        plots = fig1.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time, min=-25, max=25),
            legend=True,
        )
        plots[0].set_yaxis(self.label_rap_blood, unit=self.unit_rap, scale="linear")

        for intervention in self.interventions:
            # simulation
            plots[0].add_data(
                task=f"task_rap_{intervention}",
                xid="time",
                yid="[Cveblood_rap]",
                label=f"{self.names[intervention]} ({self.doses[intervention]:g} mg)",
                color=self.colors[intervention],
            )
            # data
            plots[0].add_data(
                dataset=f"rapamycin_{intervention}",
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label=f"{self.names[intervention]} ({self.doses[intervention]:g} mg)",
                color=self.colors[intervention],
            )
        return self.figures_by_sid(fig1)


if __name__ == "__main__":
    MODEL.run_experiments(Wang2014, output_dir=Wang2014.__name__)
