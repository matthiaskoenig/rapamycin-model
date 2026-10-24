"""Simulation experiment of the study Rogers2008."""

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


class Rogers2008(RapamycinSimulationExperiment):
    """Simulation experiment of Rogers2008.

    Single oral dose of 6 mg rapamycin (tablet) after an 8 h fast, together
    with tacrolimus 4 mg and mycophenolate mofetil 1 g (at 0 h and 12 h), in
    four dialysis dependent patients with end stage renal disease after a
    Roux-en-Y gastric bypass (P5 to P8, individual whole blood profiles; P5
    lost the intravenous access at about 13.5 h). The gastric bypass is not
    represented in the model.
    """

    groups: ClassVar[list[str]] = ["P5", "P6", "P7", "P8"]
    colors: ClassVar[dict[str, str]] = {
        "P5": "tab:blue",
        "P6": "tab:red",
        "P7": "tab:green",
        "P8": "tab:orange",
    }
    markers: ClassVar[dict[str, str]] = {
        "P5": "s",
        "P6": "o",
        "P7": "^",
        "P8": "*",
    }
    #: body weight [kg] of the patients (Tab1)
    bodyweights: ClassVar[dict[str, float]] = {
        "P5": 99,
        "P6": 93,
        "P7": 92.2,
        "P8": 91,
    }
    #: body mass index [kg/m^2] of the patients (Tab1), the height is
    #: calculated from the body weight and the body mass index
    bmis: ClassVar[dict[str, float]] = {
        "P5": 32.0,
        "P6": 33.1,
        "P7": 27.9,
        "P8": 31.5,
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
                    dset.unit_conversion("value", 1 / self.Mr.rap)

                dsets[label] = dset

        return dsets

    @override
    def simulations(self) -> dict[str, Simulation]:
        simulations = {}
        for group in self.groups:
            simulations[f"rap_{group}"] = Simulation(
                time_unit="hr",
                end=24,
                preinit_changes={
                    **self.default_changes(),
                    # physiological changes
                    "BW": Q(self.bodyweights[group], "kg"),
                    "HEIGHT": Q(
                        (self.bodyweights[group] / self.bmis[group]) ** 0.5, "m"
                    ),
                    # dose
                    "PODOSE_rap": Q(6, "mg"),
                    # end stage renal disease, dialysis dependent
                    "KI__f_renal_function": Q(
                        self.renal_map["End stage renal disease"],
                        "dimensionless",
                    ),
                },
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
                    yid="value",
                    yid_sd=None,
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
                    application_form=ApplicationForm.TABLET,
                    dosing=Dosing.SINGLE,
                    health=Health.RENAL_IMPAIRMENT,
                    fasting=Fasting.FASTED,
                    coadministration=Coadministration.TACROLIMUS_MYCOPHENOLATE_MOFETIL,
                ),
            )

        return mappings

    @override
    def figures(self) -> dict[str, Figure]:

        fig1 = Figure(
            experiment=self,
            sid="Fig1",
            name=f"{self.__class__.__name__} (6 mg, tablet)",
        )
        plots = fig1.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )
        plots[0].set_yaxis(self.label_rap_blood, unit=self.unit_rap, scale="linear")

        for group in self.groups:
            # simulation
            plots[0].add_data(
                task=f"task_rap_{group}",
                xid="time",
                yid="[Cveblood_rap]",
                label=group,
                color=self.colors[group],
            )
            # data
            plots[0].add_data(
                dataset=f"rapamycin_{group}",
                xid="time",
                yid="value",
                yid_sd=None,
                count="count",
                label=group,
                color=self.colors[group],
                marker=self.markers[group],
            )
        return self.figures_by_sid(fig1)


if __name__ == "__main__":
    MODEL.run_experiments(Rogers2008, output_dir=Rogers2008.__name__)
