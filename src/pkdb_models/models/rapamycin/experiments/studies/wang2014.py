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


class Wang2014(RapamycinSimulationExperiment):
    """Simulation experiment of Wang2014."""

    interventions = ["RAP", "RAP_CSA"]
    colors = {
        "RAP": "black",
        "RAP_CSA": "#e3256b",
    }
    f_cyp3a4 = {
        "RAP": 1.0,
        "RAP_CSA": 0.6,
    }

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

            tc0 = Timecourse(
                start=0,
                end=24* 60,  # [min]
                steps=50,
                changes={
                    **self.default_changes(),

                    # physiological changes
                    "BW": Q_(70, "kg"),  # no weight, only bmi

                    # dose
                    "PODOSE_rap": Q_(2.5, "mg"),

                    # endstage renal disease
                    #"KI__f_renal_function": Q_(self.renal_map["Mild renal impairment"], "dimensionless"),

                    # CYP3A4 activity (inhibition by cyclosporine)
                    "GU__f_cyp3a4": Q_(self.f_cyp3a4[intervention], "dimensionless"),
                    "LI__f_cyp3a4": Q_(self.f_cyp3a4[intervention], "dimensionless"),
                },
            )
            tc1 = Timecourse(
                start=0,
                end=24* 60,  # [min]
                steps=50,
                changes={
                    "PODOSE_rap": Q_(2.5, "mg"),
                },
            )
            tc2 = Timecourse(
                start=0,
                end=25 * 60,  # [min]
                steps=300,
                changes={
                    "PODOSE_rap": Q_(2.5, "mg"),
                },
            )

            tcsims[f"rap_{intervention}"] = TimecourseSim(
                [tc0] + [tc1]*13 + [tc2],
                time_offset=-24*60*14
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
                    self, task=f"task_rap_{intervention}", xid="time", yid="[Cveblood_rap]",
                ),
                metadata=RapamycinMappingMetaData(
                    tissue=Tissue.BLOOD,
                    route=Route.PO,
                    application_form=ApplicationForm.TABLET,
                    dosing=Dosing.MULTIPLE,
                    health=Health.RENAL_TRANSPLANT,
                    fasting=Fasting.FED,
                    coadministration=Coadministration.CYCLOSPORINE if "CSA" in intervention else Coadministration.NONE
                )
            )

        return mappings

    def figures(self) -> Dict[str, Figure]:

        fig1 = Figure(
            experiment=self,
            sid="Fig1",
            name=f"{self.__class__.__name__}",
        )
        plots = fig1.create_plots(xaxis=Axis(self.label_time, unit=self.unit_time, min=-25, max=25), legend=True)
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
            fig1.sid: fig1
        }


if __name__ == "__main__":
    run_experiments(Wang2014, output_dir=Wang2014.__name__)
