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


class Zimmerman1997(RapamycinSimulationExperiment):
    """Simulation experiment of Zimmerman1997."""

    interventions = ["RAP05","RAP1_5", "RAP2_5", "RAP3_5", "RAP6_5"]
    colors = {
        "RAP05": "black",
        "RAP1_5": "#FFE299",
        "RAP2_5": "#FFB347",
        "RAP3_5": "#FF944D",
        "RAP6_5": "#D2691E",
    }
    doses = {  # [mg/m^2]
        "RAP05": 0.5,
        "RAP1_5": 1.5,
        "RAP2_5": 2.5,
        "RAP3_5": 3.5,
        "RAP6_5": 6.5,
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
        for intervention in self.interventions:
            dose = self.doses[intervention] * 1.78  # scaling with body surface area
            tc0 = Timecourse(
                start=0,
                end=12 * 60,  # [min]
                steps=200,
                changes={
                    **self.default_changes(),
                    # physiological changes
                    "BW": Q_(83.2, "kg"),
                    # dose
                    "PODOSE_rap": Q_(dose, "mg"),

                    # CYP3A4 activity (inhibition by cyclosporine)
                    "GU__f_cyp3a4": Q_(self.f_cyp3a4_cyclosporine, "dimensionless"),
                    "LI__f_cyp3a4": Q_(self.f_cyp3a4_cyclosporine, "dimensionless"),
                },
            )
            tc1 = Timecourse(
                start=0,
                end=12 * 60,
                steps=200,
                changes={
                    "PODOSE_rap": Q_(dose, "mg"),
                },
            )
            tc2 = Timecourse(
                start=0,
                end=200 * 60,
                steps=2000,
                changes={
                    "PODOSE_rap": Q_(dose, "mg"),
                },
            )
            # Doses were administered twice daily for 13 days and once on the morning of study day 14.
            tcsims[f"rap_{intervention}"] = TimecourseSim(
                [tc0] + [tc1 for _ in range(25)] + [tc2]
            )

        #console.print(tcsims.keys())
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
                    application_form=ApplicationForm.SOLUTION,
                    dosing=Dosing.MULTIPLE,
                    health=Health.RENAL_TRANSPLANT,
                    fasting=Fasting.FASTED,
                    coadministration=Coadministration.CYCLOSPORINE and Coadministration.PREDNISONE
                )
            )
        return mappings

    def figures(self) -> Dict[str, Figure]:

        fig = Figure(
            experiment=self,
            sid="Fig1",
            name=f"{self.__class__.__name__}",
        )
        plots = fig.create_plots(xaxis=Axis(self.label_time, unit="day", max=30), legend=True)
        plots[0].set_yaxis(self.label_rap_blood, unit=self.unit_rap, scale="linear")

        for intervention in self.interventions:
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

        return {
            fig.sid: fig,
        }


if __name__ == "__main__":
    run_experiments(Zimmerman1997, output_dir=Zimmerman1997.__name__)
