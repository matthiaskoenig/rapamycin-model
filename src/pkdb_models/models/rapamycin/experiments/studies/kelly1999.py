from typing import Dict

from sbmlsim.data import DataSet, load_pkdb_dataframe
from sbmlsim.fit import FitMapping, FitData

from pkdb_models.models.rapamycin.experiments.base_experiment import (
    RapamycinSimulationExperiment,
)
from pkdb_models.models.rapamycin.experiments.metadata import (
    RapamycinMappingMetaData,
    Tissue, Route, Dosing, ApplicationForm, Health, Fasting, Coadministration
)
from sbmlsim.plot import Axis, Figure
from sbmlsim.simulation import Timecourse, TimecourseSim

from pkdb_models.models.rapamycin.helpers import run_experiments


class Kelly1999(RapamycinSimulationExperiment):
    """Simulation experiment of Kelly1999."""

    # FIXME: different times after transplantation
    interventions = [
        "RAP_L2, CSA",
        "RAP_L4, CSA",
        "RAP_L8, CSA",
        "RAP_T2, CSA",
        "RAP_T4, CSA",
        "RAP_T8, CSA",
    ]

    colors = {
        "RAP_L2, CSA": "black",
        "RAP_L4, CSA": "#FFB98C",
        "RAP_L8, CSA": "#D28F66",
        "RAP_T2, CSA": "black",
        "RAP_T4, CSA": "#FF944D",
        "RAP_T8, CSA": "#D2691E",
    }

    f_cyp3a4_cyclosporine = 0.6

    def datasets(self) -> Dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig1"]:
            df = load_pkdb_dataframe(f"{self.sid}_{fig_id}", data_path=self.data_path)
            for label, df_label in df.groupby("label"):
                dset = DataSet.from_df(df_label, self.ureg)

                # unit conversion
                if label.startswith("rapamycin_"):
                   dset.unit_conversion("mean", 1 / self.Mr.rap)

                dsets[label] = dset

        # console.print(dsets)

        #console.print(dsets.keys())
        return dsets

    def simulations(self) -> Dict[str, TimecourseSim]:
        Q_ = self.Q_
        tcsims = {}

        tc0 = Timecourse(
            start=0,
            end=24*60,
            steps=100,
            changes={
                **self.default_changes(),
                "BW": Q_(75.8, "kg"),
                "PODOSE_rap": Q_(3.6, "mg"),

                # CYP3A4 activity (inhibition by cyclosporine)
                "GU__f_cyp3a4": Q_(self.f_cyp3a4_cyclosporine, "dimensionless"),
                "LI__f_cyp3a4": Q_(self.f_cyp3a4_cyclosporine, "dimensionless"),
            }
        )
        tc1 = Timecourse(
            start=0,
            end=24 * 60,
            steps=50,
            changes={
                "PODOSE_rap": Q_(3.6, "mg")
            }
        )
        tc2 = Timecourse(
            start=0,
            end=14 * 60,
            steps=200,
            changes={
                "PODOSE_rap": Q_(3.6, "mg")
            }
        )

        for intervention in self.interventions:
            if "2" in intervention:
                duration = 2
            elif "4" in intervention:
                duration = 4
            elif "8" in intervention:
                duration = 8

            tcsims[f"rap_{intervention}"] = TimecourseSim(
                [tc0] + [tc1] * (7 * duration -2) + [tc2],
                time_offset= -(7*duration -1) * 24 * 60,
            )
        return tcsims

    def fit_mappings(self) -> Dict[str, FitMapping]:

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
                    count="count"
                ),
                observable=FitData(
                    self, task=f"task_rap_{intervention}", xid="time",yid="[Cveblood_rap]",
                ),
                metadata=RapamycinMappingMetaData(
                    tissue=Tissue.BLOOD,
                    route=Route.PO,
                    application_form=ApplicationForm.SOLUTION if "_L" in intervention else ApplicationForm.TABLET,
                    dosing=Dosing.MULTIPLE,
                    health=Health.RENAL_TRANSPLANT,
                    fasting=Fasting.FASTED,
                    coadministration=Coadministration.CYCLOSPORINE
                )
            )
        return mappings

    def figures(self) -> Dict[str, Figure]:
        figures = {}
        groups = {
            "L": [i for i in self.interventions if "_L" in i],
            "T": [i for i in self.interventions if "_T" in i],
        }
        for group_name, interventions in groups.items():
            fig = Figure(
                experiment=self,
                sid=f"Fig1_{group_name}",
                name=f"{self.__class__.__name__}",
            )
            plots = fig.create_plots(xaxis=Axis(self.label_time, unit="day", min=-1.5, max=0.7), legend=True)
            plots[0].set_yaxis(self.label_rap_blood, unit=self.unit_rap)

            for intervention in interventions:
                # simulation
                plots[0].add_data(
                    task=f"task_rap_{intervention}",
                    xid="time",
                    yid="[Cveblood_rap]",
                    label=intervention,
                    color=self.colors[intervention],
                )
                # data
                plots[0].add_data(
                    dataset=f"rapamycin_{intervention}",
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                    label=intervention,
                    color=self.colors[intervention],
                )
            figures[fig.sid] = fig

        return figures


if __name__ == "__main__":
    run_experiments(Kelly1999, output_dir=Kelly1999.__name__)
