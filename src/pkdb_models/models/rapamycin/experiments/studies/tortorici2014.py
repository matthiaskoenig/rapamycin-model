"""Simulation experiment of the study Tortorici2014."""

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


class Tortorici2014(RapamycinSimulationExperiment):
    """Simulation experiment of Tortorici2014.

    Rapamycin 20 mg oral solution in healthy fasted men, alone and on day 9 of
    rifampin 600 mg once daily for 14 days (Fig1, mean and SD, n = 14).
    Rifampin induces CYP3A4 and P-glycoprotein. The effect is mostly on the
    first pass: the AUCinf decreases to 18 % (geometric mean ratio, 90 % CI
    16-21 %), the Cmax to 29 % (26-32 %), while the half-life is almost
    unchanged (92 %), i.e. the systemic (hepatic) clearance changes little.
    The induction is at steady state on day 9 and maintained over the
    sampling (rifampin until day 14), so the multiple dosing of rifampin is
    represented by constant induction factors.
    """

    dose: ClassVar[float] = 20.0
    interventions: ClassVar[list[str]] = ["RAP20", "RAP20_RIF600"]
    dataset_ids: ClassVar[dict[str, str]] = {
        "RAP20": "rapamycin_RAP20",
        "RAP20_RIF600": "rapamycin_RAP20_RIF600",
    }
    names: ClassVar[dict[str, str]] = {
        "RAP20": "20 mg",
        "RAP20_RIF600": "20 mg + rifampin 600 mg/day",
    }
    colors: ClassVar[dict[str, str]] = {
        "RAP20": "black",
        "RAP20_RIF600": "#e3256b",
    }
    #: induction by rifampin (not fitted): the intestinal CYP3A4 (several
    #: fold in duodenal biopsies, Kolars1992) and P-glycoprotein (3.5 fold,
    #: Greiner1999) are strongly induced, the hepatic CYP3A4 only a little as
    #: the half-life of rapamycin is unchanged (92 %). The intestinal CYP3A4
    #: factor of 8 reproduces the AUC ratio of 0.18 with the hepatic factor
    #: 1.1 and the P-glycoprotein factor 3.5 (AUC ratio 0.18, Cmax ratio 0.18,
    #: half-life ratio 0.94 with the parameters of the model before the
    #: refit); CYP3A5 is not induced as the subjects are taken as *3/*3.
    induction: ClassVar[dict[str, dict[str, float]]] = {
        "RAP20": {
            "GU__f_cyp3a4": 1.0,
            "LI__f_cyp3a4": 1.0,
            "GU__f_pg": 1.0,
        },
        "RAP20_RIF600": {
            "GU__f_cyp3a4": 8.0,
            "LI__f_cyp3a4": 1.1,
            "GU__f_pg": 3.5,
        },
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
                dset = DataSet.from_df(df_label, self.ureg)

                dset.unit_conversion("mean", 1 / self.Mr.rap)
                dsets[label] = dset

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
                    "BW": Q(77.1, "kg"),
                    "PODOSE_rap": Q(self.dose, "mg"),
                    # CYP3A4 and P-glycoprotein activity (induction by rifampin)
                    **{
                        pid: Q(value, "dimensionless")
                        for pid, value in self.induction[intervention].items()
                    },
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
                    dataset=self.dataset_ids[intervention],
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
                    if intervention == "RAP20"
                    else Coadministration.RIFAMPICIN,
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
                label=self.names[intervention],
                color=self.colors[intervention],
            )
            # data
            plots[0].add_data(
                dataset=self.dataset_ids[intervention],
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label=self.names[intervention],
                color=self.colors[intervention],
            )

        return self.figures_by_sid(fig)


if __name__ == "__main__":
    MODEL.run_experiments(Tortorici2014, output_dir=Tortorici2014.__name__)
