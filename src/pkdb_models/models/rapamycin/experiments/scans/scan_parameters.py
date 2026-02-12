"""Parameter scans rapamycin."""
from typing import Dict

import matplotlib
import matplotlib.axes
import matplotlib.cm as cm
from matplotlib.colors import LinearSegmentedColormap, Normalize, LogNorm
from matplotlib.ticker import ScalarFormatter
from matplotlib.lines import Line2D
import numpy as np
from sbmlsim.simulation import Timecourse, TimecourseSim, ScanSim, Dimension
from sbmlsim.plot.serialization_matplotlib import FigureMPL, MatplotlibFigureSerializer
from sbmlsim.plot.serialization_matplotlib import plt
from sbmlutils.console import console

from pkdb_models.models.rapamycin.experiments.base_experiment import (
    RapamycinSimulationExperiment,
)
from pkdb_models.models.rapamycin.helpers import run_experiments


class RapamycinParameterScan(RapamycinSimulationExperiment):
    """Scan the effect of parameters on pharmacokinetics."""

    tend = 5 * 24 * 60
    steps = 2000
    dose_rap = 2  # [mg]

    num_points = 10
    scan_map = {
        "dose_scan": {
            "parameter": "PODOSE_rap",
            "default": 2,
            "range": np.sort(
                np.append(np.linspace(0.1, 20, num=num_points), [2])
            ),  # [10^-1=0.1, 10^1=10]
            "scale": "linear",
            "colormap": "PuRd",
            "units": "mg",
            "label": "rapamycin dose [mg]",
        },
        "food_scan": {
            "parameter": "GU__f_oatp",
            "default": 1.0,
            "range": np.sort(
                np.append(np.logspace(-1, 1, num=num_points), [1.0])
            ),  # [10^-1=0.1, 10^1=10]
            "scale": "log",
            "colormap": "seismic_r",
            "units": "dimensionless",
            "label": "absorption activity [-]",
        },
        "hepatic_scan": {
            "parameter": "f_cirrhosis",
            "default": 0.0,
            "range": np.linspace(0, 0.9, num=num_points),
            "scale": "linear",
            "units": "dimensionless",
            "label": "cirrhosis degree [-]",
        },
        "renal_scan": {
            "parameter": "KI__f_renal_function",
            "default": 1.0,
            "range": np.sort(
                np.append(np.logspace(-1, 1, num=num_points), [1.0])
            ),  # [10^-1=0.1, 10^1=10]
            "scale": "log",
            "units": "dimensionless",
            "label": "renal function [-]",
        },
        "cyp3a4_gu_scan": {
            "parameter": "GU__f_cyp3a4",
            "default": 1.0,
            "range": np.sort(
                np.append(np.logspace(-1, 1, num=num_points), [1.0])
            ),  # [10^-1=0.1, 10^1=10]
            "scale": "log",
            "colormap": "seismic_r",
            "units": "dimensionless",
            "label": "CYP3A4 gut activity [-]",
        },
        "cyp3a5_gu_scan": {
            "parameter": "GU__f_cyp3a5",
            "default": 1.0,
            "range": np.sort(
                np.append(np.logspace(-1, 1, num=num_points), [1.0])
            ),  # [10^-1=0.1, 10^1=10]
            "scale": "log",
            "colormap": "seismic_r",
            "units": "dimensionless",
            "label": "CYP3A5 gut activity [-]",
        },
        "cyp3a4_li_scan": {
            "parameter": "LI__f_cyp3a4",
            "default": 1.0,
            "range": np.sort(
                np.append(np.logspace(-1, 1, num=num_points), [1.0])
            ),  # [10^-1=0.1, 10^1=10]
            "scale": "log",
            "colormap": "seismic_r",
            "units": "dimensionless",
            "label": "CYP3A4 liver activity [-]",
        },
        "cyp3a5_li_scan": {
            "parameter": "LI__f_cyp3a5",
            "default": 1.0,
            "range": np.sort(
                np.append(np.logspace(-1, 1, num=num_points), [1.0])
            ),  # [10^-1=0.1, 10^1=10]
            "scale": "log",
            "colormap": "seismic_r",
            "units": "dimensionless",
            "label": "CYP3A5 liver activity [-]",
        },
    }

    def simulations(self) -> Dict[str, ScanSim]:
        Q_ = self.Q_
        tcscans = {}

        for scan_key, scan_data in self.scan_map.items():
            tcscans[f"scan_po_{scan_key}"] = ScanSim(
                simulation=TimecourseSim(
                    Timecourse(
                        start=0,
                        end=self.tend,
                        steps=self.steps,
                        changes={
                            **self.default_changes(),
                            "PODOSE_rap": Q_(self.dose_rap, "mg"),
                        },
                    )
                ),
                dimensions=[
                    Dimension(
                        "dim_scan",
                        changes={
                            scan_data["parameter"]: Q_(
                                scan_data["range"], scan_data["units"]
                            )
                        },
                    ),
                ],
            )

        return tcscans

    def figures_mpl(self) -> Dict[str, FigureMPL]:
        """Matplotlib figures."""
        # calculate pharmacokinetic parameters
        self.pk_dfs = self.calculate_rapamycin_pk()

        return {
            **self.figures_mpl_timecourses(),
            **self.figures_mpl_pharmacokinetics(),
        }

    def figures_mpl_timecourses(self) -> Dict[str, FigureMPL]:
        """Timecourse plots for key variables depending on degree of renal impairment."""

        sids = [
            "[Cve_rap]",
            "[Cve_rx]",
            "[Cveblood_rap]",
            "Aurine_rx",
            "Afeces_rx",
        ]

        figures = {}
        for scan_key, scan_data in self.scan_map.items():
            range_vals = scan_data["range"]
            rmin, rmax = range_vals[0], range_vals[-1]

            if scan_key == "hepatic_scan":
                cmap = LinearSegmentedColormap.from_list(
                    "cirrhosis_blues",
                    [
                        self.cirrhosis_colors["Mild cirrhosis"],
                        self.cirrhosis_colors["Moderate cirrhosis"],
                        self.cirrhosis_colors["Severe cirrhosis"],
                    ],
                )
                norm = Normalize(vmin=0.0, vmax=0.9, clip=False)

            elif scan_key == "renal_scan":
                vmin, vmax = 0.1, 10.0

                def _pos(v):
                    return (np.log10(v) - np.log10(vmin)) / (np.log10(vmax) - np.log10(vmin))

                severe = self.renal_colors["Severe renal impairment"]       # ~0.1
                moderate = self.renal_colors["Moderate renal impairment"]   # ~0.32
                mild = self.renal_colors["Mild renal impairment"]           # ~0.69
                high_ext = "#e5f5e0"  # very light green for >1.0 extension

                cmap = LinearSegmentedColormap.from_list(
                    "renal_greens",
                    [
                        (_pos(0.1), severe),
                        (_pos(0.32), moderate),
                        (_pos(0.69), mild),
                        (_pos(10.0), high_ext),
                    ],
                )
                norm = LogNorm(vmin=vmin, vmax=vmax, clip=False)

            else:
                cmap_str = scan_data["colormap"]
                cmap = matplotlib.colormaps.get_cmap(cmap_str)
                if scan_data["scale"] == "linear":
                    norm = Normalize(vmin=rmin, vmax=rmax, clip=False)
                else:
                    norm = LogNorm(vmin=rmin, vmax=rmax, clip=False)

            f, axes = plt.subplots(
                nrows=1,
                ncols=len(sids),
                figsize=(6 * len(sids), 6),
                dpi=150,
                layout="constrained"
            )

            ymax = {}
            for kcol, sid in enumerate(sids):
                ymax[sid] = 0.0
                ax = axes[kcol]

                # get data
                Q_ = self.Q_
                xres = self.results[f"task_scan_po_{scan_key}"]

                # scanned dimension
                scandim = xres._redop_dims()[0]
                parameter_id = scan_data["parameter"]
                par_vec = Q_(xres[parameter_id].values[0], xres.uinfo[parameter_id])
                t_vec = xres.dim_mean("time").to(self.units["time"])

                # placeholders for default/reference curve
                t_vec_default = None
                c_vec_default = None

                for k_par, par in enumerate(par_vec):
                    c_vec = Q_(
                        xres[sid].sel({scandim: k_par}).values,
                        xres.uinfo[sid],
                    ).to(self.units[sid])

                    # update ymax
                    cmax = np.nanmax(c_vec.magnitude)
                    if cmax > ymax[sid]:
                        ymax[sid] = cmax

                    linewidth = 2.0
                    if np.isclose(scan_data["default"], par.magnitude):
                        color = "black"
                        t_vec_default = t_vec
                        c_vec_default = c_vec
                    else:
                        color = cmap(norm(par.magnitude))

                    ax.plot(
                        t_vec.magnitude,
                        c_vec.magnitude,
                        color=color,
                        linewidth=linewidth,
                    )

                # plot the reference line in black (if present)
                if t_vec_default is not None and c_vec_default is not None:
                    ax.plot(
                        t_vec_default.magnitude,
                        c_vec_default.magnitude,
                        color="black",
                        linewidth=2.0,
                    )

                ax.set_xlabel(
                    f"{self.label_time} [{self.units['time']}]",
                    fontdict=self.font,
                )
                ax.tick_params(axis="x", labelsize=self.tick_font_size)
                ax.tick_params(axis="y", labelsize=self.tick_font_size)

                ax.set_ylabel(
                    f"{self.labels[sid]} [{self.units[sid]}]",
                    fontdict=self.font,
                )

            cb_ax = f.add_axes(rect=[0.08, 0.85, 0.1, 0.08])
            cb_ax.set_in_layout(True)

            cbar = f.colorbar(
                cm.ScalarMappable(norm=norm, cmap=cmap),
                cax=cb_ax,
                orientation="horizontal",
            )

            if scan_key == "hepatic_scan":
                ticks = [0.0, 0.9]
                cbar.set_ticks(ticks)
                cbar.set_ticklabels(ticks, **{"size": 15, "weight": "medium"})
                cbar.ax.axvline(x=scan_data["default"], color="black", linewidth=2)

            elif scan_key == "renal_scan":
                ticks = [0.1, 10.0]
                cbar.set_ticks(ticks)
                cbar.set_ticklabels(ticks, **{"size": 15, "weight": "medium"})
                cbar.ax.axvline(x=scan_data["default"], color="black", linewidth=2)

            else:
                ticks = [rmin, rmax]
                if scan_data["default"] not in ticks:
                    ticks.append(scan_data["default"])
                    ticks = sorted(ticks)
                cbar.set_ticks(ticks)
                cbar.set_ticklabels(ticks, **{"size": 15, "weight": "medium"})
                cbar.ax.axvline(x=scan_data["default"], color="black", linewidth=2)

            cbar.ax.set_xlabel(
                scan_data["label"], **{"size": 15, "weight": "bold"}
            )
            console.print(f"{scan_key} colorbar ticks: {ticks}")

            figures[f"timecourse__{scan_key}"] = f

        return figures

    def figures_mpl_pharmacokinetics(self):
        """Visualize dependency of pharmacokinetics parameters."""
        Q_ = self.Q_
        figures = {}

        parameters_info = {
            "rap": ["aucinf", "cmax", "kel", "thalf"],
            "rx": ["aucinf", "cmax", "kel", "thalf"],
        }

        # analyte line/edge colours
        series_colors = {"rap": "black", "rx": "darkgrey"}

        # marker sizing
        marker_size_pts = 11.0
        scatter_s = marker_size_pts ** 2
        edge_width = 1.8
        line_width = 2.0

        for scan_key, scan_data in self.scan_map.items():
            range_vals = scan_data["range"]
            rmin, rmax = float(range_vals[0]), float(range_vals[-1])

            if scan_key == "hepatic_scan":
                cmap = LinearSegmentedColormap.from_list(
                    "cirrhosis_blues",
                    [
                        self.cirrhosis_colors["Mild cirrhosis"],
                        self.cirrhosis_colors["Moderate cirrhosis"],
                        self.cirrhosis_colors["Severe cirrhosis"],
                    ],
                )
                norm = Normalize(vmin=0.0, vmax=0.9, clip=False)

            elif scan_key == "renal_scan":
                vmin, vmax = 0.1, 10.0

                def _pos(v):
                    return (np.log10(v) - np.log10(vmin)) / (np.log10(vmax) - np.log10(vmin))

                severe = self.renal_colors["Severe renal impairment"]
                moderate = self.renal_colors["Moderate renal impairment"]
                mild = self.renal_colors["Mild renal impairment"]
                high_ext = "#e5f5e0"
                cmap = LinearSegmentedColormap.from_list(
                    "renal_greens",
                    [
                        (_pos(0.1), severe),
                        (_pos(0.32), moderate),
                        (_pos(0.69), mild),
                        (_pos(10.0), high_ext),
                    ],
                )
                norm = LogNorm(vmin=vmin, vmax=vmax, clip=False)

            else:
                cmap_str = scan_data["colormap"]
                cmap = matplotlib.colormaps.get_cmap(cmap_str)
                if scan_data["scale"] == "linear":
                    norm = Normalize(vmin=rmin, vmax=rmax, clip=False)
                else:
                    norm = LogNorm(vmin=rmin, vmax=rmax, clip=False)

            metrics = parameters_info["rap"]
            f, axes = plt.subplots(
                nrows=1, ncols=len(metrics), figsize=(6 * len(metrics), 6), dpi=150, layout="constrained"
            )
            axes = axes.flatten()

            sim_key = f"scan_po_{scan_key}"
            xres = self.results[f"task_{sim_key}"]
            df_all = self.pk_dfs[sim_key]

            parameter_id = scan_data["parameter"]
            x_q = Q_(xres[parameter_id].values[0], xres.uinfo[parameter_id])
            x_vals = np.asarray([float(v) for v in x_q.magnitude])

            for k, pk_key in enumerate(metrics):
                ax = axes[k]
                ax.axvline(x=scan_data["default"], color="grey", linestyle="--", linewidth=1.2, zorder=1)

                for substance in ["rap", "rx"]:
                    df = df_all[df_all.substance == substance].copy()
                    yq = Q_(df[pk_key].to_numpy(), df[f"{pk_key}_unit"].values[0]).to(self.pk_units[pk_key])
                    y = yq.magnitude

                    ax.plot(
                        x_vals,
                        y,
                        linestyle="-",
                        linewidth=line_width,
                        color=series_colors[substance],
                        alpha=0.9,
                        zorder=2,
                    )

                    facecols = [cmap(norm(x)) for x in x_vals]
                    ax.scatter(
                        x_vals,
                        y,
                        s=scatter_s,
                        facecolors=facecols,
                        edgecolors=series_colors[substance],
                        linewidths=edge_width,
                        zorder=3,
                    )

                # axes formatting
                ax.tick_params(axis="x", labelsize=self.tick_font_size)
                ax.tick_params(axis="y", labelsize=self.tick_font_size)
                ax.set_xlabel(scan_data["label"], fontdict=self.scan_font)
                ax.set_ylabel(f"{self.pk_labels[pk_key]} [{self.pk_units[pk_key]}]", fontdict=self.scan_font)

                if scan_data["scale"] == "log":
                    ax.set_xscale("log")
                    ax.xaxis.set_major_formatter(ScalarFormatter())
                    ax.ticklabel_format(style="plain", axis="x")

                if k == 0:
                    legend_handles = [
                        Line2D([0], [0], marker="o", linestyle="", markerfacecolor="white",
                               markeredgecolor=series_colors["rap"], markeredgewidth=edge_width,
                               markersize=marker_size_pts, label="rap"),
                        Line2D([0], [0], marker="o", linestyle="", markerfacecolor="white",
                               markeredgecolor=series_colors["rx"], markeredgewidth=edge_width,
                               markersize=marker_size_pts, label="rx"),
                    ]
                    ax.legend(handles=legend_handles,
                              loc="upper left", frameon=False,
                              fontsize=self.legend_font_size)

            figures[f"pk_{scan_key}"] = f

        return figures


if __name__ == "__main__":
    run_experiments(RapamycinParameterScan, output_dir=RapamycinParameterScan.__name__)
