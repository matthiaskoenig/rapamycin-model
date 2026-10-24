#set page(paper: "a4", margin: 2cm)
#set text(size: 10pt)
#set heading(numbering: "1.1")
#set math.equation(numbering: "(1)")
// a table longer than a page continues on the next page, its
// caption above it and its header repeated
#show figure.where(kind: table): set block(breakable: true)
#show figure.where(kind: table): set figure.caption(position: top)

#align(center, text(size: 16pt, weight: "bold")[Rapamycin PBPK model: fit PK])

= Methods: Rapamycin PBPK model: fit PK

== Parameters and data

9 parameters were optimized: Ka\_dis \[1\/hr\] dissolution rapamycin (#raw("GU__Ka_dis_rap")), rate of rapamycin metabolism by CYP3A4 (#raw("LI__RAP2RX_k")), tissue\/blood partition coefficient rap (#raw("Kp_rap")), tissue distribution rap (#raw("ftissue_rap")), bound\/free partition coefficient rap (tissue) (#raw("Kpdeep_rap")), rate of the slow binding rap (tissue) (#raw("kdeep_rap")), fed effect on the fraction absorbed (#raw("GU__FED_f_absorption")), fed effect on the dissolution and gastric emptying (#raw("GU__FED_f_dissolution")) and capacity of CYP3A5 relative to CYP3A4 (#raw("r_cyp3a5")).

A fit mapping is a curve of a study, the data and the observable of the model it is compared with. The fit has 22 fit mappings of training data with 308 data points, whose residuals enter the objective; and 32 fit mappings of validation data with 479 data points, which were simulated and evaluated with the fitted parameters but not fitted. The fit was set up and run with sbmlsim #cite(<Koenig2026_sbmlsim>).

== Objective function

The parameters were estimated by minimizing the objective function

$ F(arrow(p)) = 1/2 sum_k sum_i rho(w_(i,k) r_(i,k)(arrow(p))^2) $

where $k$ runs over the fit mappings of the training data and $i$ over their data points, $r_(i,k)$ is the residual, $w_(i,k) = w_k w^p_(i,k)$ the weight of the data point and $rho$ the loss function, applied to the squared weighted residual $z = w_(i,k) r_(i,k)^2$. The loss function is linear, $rho(z) = z$, i.e. a standard least squares problem.

The residuals are normalized by the mean $overline(y)_k$ of the data of their curve, which makes curves of different quantities and magnitudes comparable:

$ r_(i,k) = (y_(i,k) - f_k(x_(i,k), arrow(p))) / overline(y)_k $

The weight of curve $k$ is $w_k = n_k / N_k$, the weight $n_k$ of its fit mapping (the weight its fit mapping collection assigns, 1 by default, or, for a collection which uses the weights of its mappings, the weight of the mapping, by default the number of subjects of its data) divided by the number $N_k$ of its data points.

A data point is weighted with the inverse of its coefficient of variation, $w^p_(i,k) = |y_(i,k)| / sigma_(i,k)$, with $sigma_(i,k)$ the standard deviation of the data point, or its standard error for data without a standard deviation. An error of zero is missing, a missing error of a curve with errors is replaced by the largest error of the curve, and the data points of a curve without errors have the weight $w^p_(i,k) = 2$, the weight of a coefficient of variation of 0.5.

== Optimization

The optimizer searched the parameters in log10 space within their bounds; bounds, values and intervals are given in the units of the fit parameters, the units of the table of the parameters.

The objective was minimized in 3 runs with differential evolution #cite(<Storn1997>) of SciPy #cite(<Virtanen2020>) within the bounds of the parameters, with the strategy #raw("best1bin"), a population of 15 times the number of parameters, at most 1000 generations, the relative tolerance 0.01 of convergence, the recombination constant 0.7 and a mutation dithered between 0.5 and 1. The initial population of a run was drawn by latin hypercube sampling #cite(<McKay1979>), with a seed of every run derived from the seed 1234 of the fit. The best member of the final population was polished with the L-BFGS-B method of SciPy.

The convergence of the runs was assessed from the traces of their cost and the waterfall plot of their final costs, and the parameter set with the lowest cost was used.

The models were simulated with libRoadRunner #cite(<Welsh2023>) and its integrator CVODE #cite(<Hindmarsh2005>) with the relative tolerance $10^(-6)$ and the absolute tolerances $10^(-13)$ of an amount, $10^(-13)$ of a concentration and $10^(-6)$ of the other states, which for a species in a concentration is multiplied by the initial volume of its compartment. It ran with the variable step size of libRoadRunner, i.e. with an output at every step of the integrator. A data point at the time of a change after the start of a simulation, e.g. a dose or the reset of a collection, was compared with the state before the change, except the first data point of a curve with data at later times, which starts a collection and was compared with the state after the change.

== Goodness of fit

The goodness of fit is given per kind of data (training data and validation data) by metrics over its $n$ data points with the residuals $e_i = y_i - f(x_i, hat(arrow(p)))$ at the fitted parameters $hat(arrow(p))$, the mean squared error, its root and the root mean squared error normalized by the mean $overline(y)_(k(i))$ of the data of the curve of a point:

$ "MSE" = 1/n sum_(i=1)^n e_i^2, quad "RMSE" = sqrt("MSE"), quad "NRMSE" = sqrt(1/n sum_(i=1)^n (e_i / overline(y)_(k(i)))^2) $

The root mean square of the residuals of the objective, i.e. what the fit minimizes, and the coefficient of determination:

$ "RMSE"_w = sqrt(1/n sum_(i=1)^n rho(w_i r_i^2)), quad R^2 = 1 - (sum_i e_i^2) / (sum_i (y_i - overline(y))^2) $

The predictions of a non-linear model are no linear regression of the data, so $R^2$ can be negative: the prediction is worse than the mean of the data. The information criteria of Akaike #cite(<Akaike1974>) and Schwarz #cite(<Schwarz1978>) with the number $n_p$ of fitted parameters are those of a least squares fit with normally distributed residuals, up to an additive constant, so only differences between models fitted to the same data are meaningful:

$ "AIC" = n ln("MSE") + 2 n_p, quad "BIC" = n ln("MSE") + n_p ln(n) $

The percentage of the data points whose prediction $f_i = f(x_i, hat(arrow(p)))$ lies within twofold of the data is taken over the $n_(+)$ data points with positive data $y_i > 0$ and a positive prediction $f_i > 0$; the others have no ratio and are not counted, and a kind of data without such a point has no percentage:

$ P_2 = 100 / n_(+) abs(\{ i : 1/2 <= f_i / y_i <= 2 \}) $

The goodness of fit plot shows the predictions against the data on logarithmic axes, with the lines of twofold $f = 2y$ and $f = y/2$. The Bland-Altman plot #cite(<Bland1986>) shows the ratio of prediction and data as $d_i = log_10(f_i / y_i)$ against the geometric mean $sqrt(f_i y_i)$ for the data points with positive values, with the bias $overline(d)$ and the limits of agreement $overline(d) plus.minus 1.96 s_d$, with the sample standard deviation $s_d$ of $d_i$, of the training data; they are read as the fold factors $10^(overline(d))$ and $10^(overline(d) plus.minus 1.96 s_d)$.

== Identifiability

The identifiability of the parameters was analysed with the profile likelihood #cite(<Raue2009>) #cite(<Kreutz2013>). The profile of a parameter is the lowest objective with the parameter fixed:

$ "PL"(p_j) = min_(p_(l != j)) F(arrow(p)) $

At every point of a scan the other parameters were re-optimized from the previous point with the least squares optimizer of SciPy with the relative step 0.05. A scan ran from the optimum in both directions in the space of the optimizer, with steps in decades of the parameter, starting at 0.1. A step which raised the objective by more than 0.2 of the distance between the objective at the optimum and the threshold was reduced by the factor 2 and repeated, down to 0.01, and the step after one which raised it by less than that divided by the factor was enlarged by the factor, up to 1. A scan stopped at the threshold, at the bound of the parameter, at a point which failed or after 50 points.

The confidence intervals are the values at which the profile crosses the threshold of the likelihood ratio test

$ F_"thr" = F_"min" + chi^2_alpha ("df") / 2 $

with $F_"min"$ the lowest objective of the fit and of all profiles and $chi^2_alpha ("df")$ the quantile of the $chi^2$ distribution at the confidence level $alpha = 0.95$ with df \= 1 degrees of freedom, i.e. 1.92 above the lowest objective. The threshold has this confidence level if $F$ is the negative log-likelihood of the data up to a constant, i.e. if the weighted residuals are independent and normally distributed with unit variance; otherwise the level is approximate. A crossing is interpolated between two points of the profile, on the logarithm of a parameter on a logarithmic scale. A parameter is identifiable if its profile reaches the threshold within the scan on both sides of the optimum, i.e. its confidence interval is finite; practically non-identifiable towards the lower bound if its profile reaches the threshold above the optimum but not below it within the scan, which ends at the lower bound or after 50 points; practically non-identifiable towards the upper bound if its profile reaches the threshold below the optimum but not above it within the scan, which ends at the upper bound or after 50 points; practically non-identifiable if its profile reaches the threshold on neither side within the scans, which end at the bounds or after 50 points, and is not flat; and structurally non-identifiable if its profile is flat on both sides. A side of a profile is flat if it rises by less than 0.05 of the distance between the lowest objective and the threshold, and a profile with a point which failed or did not converge is not classified. The profile likelihood is preferred over the intervals of the Fisher information, its intervals are invariant under a transformation of the parameters and can be asymmetric #cite(<Wieland2021>).

The Fisher information is the local alternative to the profiles: it is the Gauss-Newton approximation of the curvature of the objective at the optimum, from the Jacobian $J$ of the weighted residuals of the training data with respect to the parameters $arrow(p)$ in the space of the optimizer, computed by central finite differences with the step $h_j = epsilon max(|p_j|, 1), space epsilon = 10^(-6)$ of the parameter $p_j$:

$ I = J^top J, quad "cov" = sigma^2 I^(-1), quad sigma^2 = (2 F) / (n - n_p) $

with the number $n$ of data points of the training data and the number $n_p$ of parameters. The information has rank 9 of 9, the number of its eigenvalues above $10^(-8)$ of the largest one, and the condition number 1020; the covariance of an information of a lower rank is computed with the pseudo-inverse. The Wald intervals are $hat(p)_j plus.minus t_((1+alpha)/2, n-n_p) "SE"_j$ with the standard errors $"SE"_j = sqrt("cov"_(j j))$ at the confidence level $alpha = 0.95$, transformed back into the units of the fit parameters.

== Provenance

The fit PK\_20261017\_094949\_\_24b18 was run on 2026-10-17 07:49:51 UTC with sbmlsim 0.10.0 #cite(<Koenig2026_sbmlsim>), libRoadRunner 2.10.0 #cite(<Welsh2023>), SciPy 1.18.1 #cite(<Virtanen2020>), NumPy 2.5.3 #cite(<Harris2020>) and Python 3.14.4 (seed 1234). The models and data were those of commit #raw("138cd5e26b655379e4c34971f8cb649aa6cb81cd") of #raw("/home/mkoenig/git/pkdb_models-rapamycin"); their SHA-256 checksums are listed in #ref(<tab:checksums:PK_20261017_094949__24b18>).

#figure(
  table(
  columns: 2,
  align: (left, left,),
  table.header([File], [SHA-256]),
  raw("pkdb_models/models/rapamycin/data/rapamycin/BasaDenes2019/.BasaDenes2019_Fig1A.tsv"),
  raw("13006e53327e..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/BasaDenes2019/.BasaDenes2019_Fig1C.tsv"),
  raw("75f1a06f66cc..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/BasaDenes2019/.BasaDenes2019_Tab2.tsv"),
  raw("538a078117f8..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/BasaDenes2019/.BasaDenes2019_TabGroups.tsv"),
  raw("131175414a11..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/BasaDenes2019/.BasaDenes2019_TabIndividuals.tsv"),
  raw("f924ac035979..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/BasaDenes2019/.BasaDenes2019_TabInterventions.tsv"),
  raw("49f81db9d585..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Bottiger2001/.Bottiger2001_Fig1.tsv"),
  raw("fb4e1d2fb87a..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Bottiger2001/.Bottiger2001_Fig2.tsv"),
  raw("07207cba6fe4..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Bottiger2001/.Bottiger2001_Fig2A.tsv"),
  raw("08691b58bbb4..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Bottiger2001/.Bottiger2001_Fig2B.tsv"),
  raw("fbb5bbb441c8..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Bottiger2001/.Bottiger2001_Fig3.tsv"),
  raw("b015ed4f7389..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Bottiger2001/.Bottiger2001_Tab1.tsv"),
  raw("0a4267b84efb..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Bottiger2001/.Bottiger2001_Tab2.tsv"),
  raw("3d61fb9188dc..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Bottiger2001/.Bottiger2001_TabGroups.tsv"),
  raw("bf66a2f415f4..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Bottiger2001/.Bottiger2001_TabIndividuals.tsv"),
  raw("f2a01715e7f8..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Bottiger2001/.Bottiger2001_TabInterventions.tsv"),
  raw("d97a0a2d7ed9..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Brattstrom2000/.Brattstrom2000_Fig1.tsv"),
  raw("c580613148d3..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Brattstrom2000/.Brattstrom2000_Fig2.tsv"),
  raw("6777cec976ed..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Brattstrom2000/.Brattstrom2000_Fig3.tsv"),
  raw("c680268839fa..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Brattstrom2000/.Brattstrom2000_Tab3.tsv"),
  raw("fb4bd0572679..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Brattstrom2000/.Brattstrom2000_Tab4.tsv"),
  raw("89ad59cb83db..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Brattstrom2000/.Brattstrom2000_TabGroups.tsv"),
  raw("e583f7ad4731..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Brattstrom2000/.Brattstrom2000_TabIndividuals.tsv"),
  raw("811aaf8d038b..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Brattstrom2000/.Brattstrom2000_TabInterventions.tsv"),
  raw("294a2affd173..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Kelly1999/.Kelly1999_Fig1.tsv"),
  raw("67a785017573..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Kelly1999/.Kelly1999_Tab2.tsv"),
  raw("ee4a399e466e..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Kelly1999/.Kelly1999_Tab3.tsv"),
  raw("85ae5e50f51f..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Kelly1999/.Kelly1999_TabGroups.tsv"),
  raw("1ab69cc2c5d0..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Kelly1999/.Kelly1999_TabInterventions.tsv"),
  raw("2662d51a8cf5..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/KorthBradley2012/.KorthBradley2012_Fig1.tsv"),
  raw("65b0c8650ac5..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/KorthBradley2012/.KorthBradley2012_Fig2.tsv"),
  raw("924012e8580d..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/KorthBradley2012/.KorthBradley2012_Tab1.tsv"),
  raw("5243aeeb4d12..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/KorthBradley2012/.KorthBradley2012_TabGroups.tsv"),
  raw("98eae14b8b4c..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/KorthBradley2012/.KorthBradley2012_TabInterventions.tsv"),
  raw("ffd8b9f8fc89..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Leelahavanichkul2005/.Leelahavanichkul2005_Fig1.tsv"),
  raw("a1d90dac9290..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Leelahavanichkul2005/.Leelahavanichkul2005_Tab3.tsv"),
  raw("9c0869686e57..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Leelahavanichkul2005/.Leelahavanichkul2005_Tab3A.tsv"),
  raw("9228a4518ca7..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Leelahavanichkul2005/.Leelahavanichkul2005_TabGroups.tsv"),
  raw("fce6b139031e..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Leelahavanichkul2005/.Leelahavanichkul2005_TabIndividuals.tsv"),
  raw("5a7e18c59123..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Leelahavanichkul2005/.Leelahavanichkul2005_TabInterventions.tsv"),
  raw("a583b930eb2b..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Leung2006/.Leung2006_Fig1.tsv"),
  raw("cb15f2cb0453..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Leung2006/.Leung2006_Fig2.tsv"),
  raw("5a1c52fcf55a..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Leung2006/.Leung2006_Tab1.tsv"),
  raw("102ab99e42a9..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Leung2006/.Leung2006_Tab1A.tsv"),
  raw("b697031a63f1..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Leung2006/.Leung2006_Tab2.tsv"),
  raw("e3d0ffa70b99..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Leung2006/.Leung2006_Tab2A.tsv"),
  raw("ba474077480b..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Leung2006/.Leung2006_Tab3.tsv"),
  raw("7912c555985f..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Leung2006/.Leung2006_Tab3A.tsv"),
  raw("eb0a81143bf3..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Leung2006/.Leung2006_TabGroups.tsv"),
  raw("46bd4793ec75..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Leung2006/.Leung2006_TabIndividuals.tsv"),
  raw("9a8fe0cfb5ce..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Leung2006/.Leung2006_TabInterventions.tsv"),
  raw("6d7aa16fd67a..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Renders2007/.Renders2007_Fig3.tsv"),
  raw("0df0bddfe28e..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Renders2007/.Renders2007_Fig4.tsv"),
  raw("565d75d1f4d0..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Renders2007/.Renders2007_Tab2.tsv"),
  raw("53794ef66237..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Renders2007/.Renders2007_Tab3.tsv"),
  raw("40e2060550a0..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Renders2007/.Renders2007_TabGroups.tsv"),
  raw("794762a06fe0..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Renders2007/.Renders2007_TabInterventions.tsv"),
  raw("dd93d598090e..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Rogers2008/.Rogers2008_Fig1.tsv"),
  raw("c16175688c37..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Rogers2008/.Rogers2008_Fig2.tsv"),
  raw("dde8f3a111ff..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Rogers2008/.Rogers2008_Fig3.tsv"),
  raw("b7f8ad4e6d3b..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Rogers2008/.Rogers2008_Fig4.tsv"),
  raw("2a5a833be9c7..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Rogers2008/.Rogers2008_Tab1.tsv"),
  raw("e5a58b02c2ad..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Rogers2008/.Rogers2008_Tab2.tsv"),
  raw("a2390d3901ab..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Rogers2008/.Rogers2008_Tab3.tsv"),
  raw("d4350ffebca8..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Rogers2008/.Rogers2008_Tab4.tsv"),
  raw("070098eba6bf..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Rogers2008/.Rogers2008_TabGroups.tsv"),
  raw("3ad97e0734d5..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Rogers2008/.Rogers2008_TabInterventions.tsv"),
  raw("bb78c72e1407..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Tortorici2013/.Tortorici2013_Fig1.tsv"),
  raw("d8520a4e20c3..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Tortorici2013/.Tortorici2013_Fig2.tsv"),
  raw("03a53753cdb0..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Tortorici2013/.Tortorici2013_Tab1.tsv"),
  raw("3f03d7aecac4..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Tortorici2013/.Tortorici2013_Tab2.tsv"),
  raw("615dd9564327..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Tortorici2013/.Tortorici2013_TabGroups.tsv"),
  raw("ac44e6798412..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Tortorici2013/.Tortorici2013_TabInterventions.tsv"),
  raw("4f7df7473aa5..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Tortorici2014/.Tortorici2014_Fig1.tsv"),
  raw("0c8f186cc2ac..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Tortorici2014/.Tortorici2014_Tab1.tsv"),
  raw("abba1042c8e5..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Tortorici2014/.Tortorici2014_TabGroups.tsv"),
  raw("89f51cae990c..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Tortorici2014/.Tortorici2014_TabInterventions.tsv"),
  raw("2aee1a655c60..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Wang2014/.Wang2014_Fig1.tsv"),
  raw("2ace696344df..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Wang2014/.Wang2014_Fig3.tsv"),
  raw("d7a7b60f6514..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Wang2014/.Wang2014_Tab1.tsv"),
  raw("5e7a49654cdd..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Wang2014/.Wang2014_TabGroups.tsv"),
  raw("d158d890c392..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Wang2014/.Wang2014_TabIndividuals.tsv"),
  raw("a8be4c361011..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Wang2014/.Wang2014_TabInterventions.tsv"),
  raw("386c7de063b1..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Zha2022/.Zha2022_Fig2.tsv"),
  raw("f34b0903cc21..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Zha2022/.Zha2022_Fig3.tsv"),
  raw("432e4727302b..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Zha2022/.Zha2022_Tab1.tsv"),
  raw("b0b5c0c3120a..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Zha2022/.Zha2022_Tab2.tsv"),
  raw("8ddf2a7f9e9a..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Zha2022/.Zha2022_TabGroups.tsv"),
  raw("35f121debc24..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Zha2022/.Zha2022_TabInterventions.tsv"),
  raw("058004c07f72..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Zhang2017/.Zhang2017_Fig1A.tsv"),
  raw("67b022e9c0f7..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Zhang2017/.Zhang2017_Fig1B.tsv"),
  raw("17bb9402f369..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Zhang2017/.Zhang2017_Tab4.tsv"),
  raw("aeebed617f32..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Zhang2017/.Zhang2017_Tab6.tsv"),
  raw("9edd0b2873fe..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Zhang2017/.Zhang2017_TabGroups.tsv"),
  raw("abb9ad0260c6..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Zhang2017/.Zhang2017_TabInterventions.tsv"),
  raw("4b53a9b9f25b..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Zimmerman1997/.Zimmerman1997_Fig1.tsv"),
  raw("5c57542e6fdc..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Zimmerman1997/.Zimmerman1997_Fig2.tsv"),
  raw("482ccc06e7c0..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Zimmerman1997/.Zimmerman1997_Fig3.tsv"),
  raw("d16982579e20..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Zimmerman1997/.Zimmerman1997_Fig4.tsv"),
  raw("1fed866e9a85..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Zimmerman1997/.Zimmerman1997_Tab2.tsv"),
  raw("4495d595e5c7..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Zimmerman1997/.Zimmerman1997_TabGroups.tsv"),
  raw("4fc5ebaae762..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Zimmerman1997/.Zimmerman1997_TabIndividuals.tsv"),
  raw("2bb7f6fce1c4..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Zimmerman1997/.Zimmerman1997_TabInterventions.tsv"),
  raw("794b07c4f898..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Zimmerman1999/.Zimmerman1999_Fig1.tsv"),
  raw("67be37dab38d..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Zimmerman1999/.Zimmerman1999_Fig2.tsv"),
  raw("84c5e98cbe15..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Zimmerman1999/.Zimmerman1999_Tab1.tsv"),
  raw("2f18356303a3..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Zimmerman1999/.Zimmerman1999_TabGroups.tsv"),
  raw("99647fb116b8..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Zimmerman1999/.Zimmerman1999_TabIndividuals.tsv"),
  raw("cc5759fd82ea..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Zimmerman1999/.Zimmerman1999_TabInterventions.tsv"),
  raw("c8f5ecd371ef..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Zimmerman2003/.Zimmerman2003_Fig1.tsv"),
  raw("ddecd6ef2916..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Zimmerman2003/.Zimmerman2003_Fig2.tsv"),
  raw("0a8154d390c6..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Zimmerman2003/.Zimmerman2003_Fig3.tsv"),
  raw("f8d4e5d95e71..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Zimmerman2003/.Zimmerman2003_Fig4.tsv"),
  raw("b7aa09394894..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Zimmerman2003/.Zimmerman2003_Tab1.tsv"),
  raw("1fe7c403d91a..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Zimmerman2003/.Zimmerman2003_Tab2.tsv"),
  raw("64824d228af5..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Zimmerman2003/.Zimmerman2003_TabGroups.tsv"),
  raw("a2d068c934a5..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Zimmerman2003/.Zimmerman2003_TabIndividuals.tsv"),
  raw("3e5fd6d3d543..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Zimmerman2003/.Zimmerman2003_TabInterventions.tsv"),
  raw("a7402cfb0c65..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Zimmerman2005/.Zimmerman2005_Tab2.tsv"),
  raw("403f1434f473..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Zimmerman2005/.Zimmerman2005_TabGroups.tsv"),
  raw("4356c1f9ddb6..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Zimmerman2005/.Zimmerman2005_TabInterventions.tsv"),
  raw("fd4763f3c6e5..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Zimmerman2008/.Zimmerman2008_Fig1.tsv"),
  raw("3b0b3f20b169..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Zimmerman2008/.Zimmerman2008_Tab2.tsv"),
  raw("1f91affc6218..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Zimmerman2008/.Zimmerman2008_TabGroups.tsv"),
  raw("f5cf2f65304c..."),
  raw("pkdb_models/models/rapamycin/data/rapamycin/Zimmerman2008/.Zimmerman2008_TabInterventions.tsv"),
  raw("38ec39cadeaa..."),
  raw("pkdb_models/models/rapamycin/models/output/models/rapamycin_body_flat.xml"),
  raw("147d4aa9d8ed..."),
  ),
  caption: [SHA-256 checksums of the models and data of the fit PK\_20261017\_094949\_\_24b18, their first 12 characters; the full checksums are in the provenance of the result (optimization\_result.json).],
) <tab:checksums:PK_20261017_094949__24b18>

= Results: Rapamycin PBPK model: fit PK

== Parameters

#figure(
  table(
  columns: 6,
  align: (left, left, right, left, right, left,),
  table.header([Parameter], [Description], [Value \[95% CI\]], [Unit], [Bounds], [Identifiable]),
  raw("GU__Ka_dis_rap"),
  [Ka\_dis \[1\/hr\] dissolution rapamycin],
  [3.65 \[2.13, 7.88\]],
  [1\/hr],
  [\[0.1, 20\]],
  [yes],
  raw("LI__RAP2RX_k"),
  [rate of rapamycin metabolism by CYP3A4],
  [0.0262 \[LB, 0.0393\]],
  [1\/min],
  [\[0.001, 1\]],
  [no (LB)],
  raw("Kp_rap"),
  [tissue\/blood partition coefficient rap],
  [0.648 \[0.551, 0.82\]],
  [dimensionless],
  [\[0.3, 30\]],
  [yes],
  raw("ftissue_rap"),
  [tissue distribution rap],
  [3.72 \[1.66, UB\]],
  [l\/min],
  [\[0.001, 10\]],
  [no (UB)],
  raw("Kpdeep_rap"),
  [bound\/free partition coefficient rap (tissue)],
  [3.5 \[1.67, 19.8\]],
  [dimensionless],
  [\[0.1, 100\]],
  [yes],
  raw("kdeep_rap"),
  [rate of the slow binding rap (tissue)],
  [0.000648 \[0.000122, 0.00143\]],
  [1\/min],
  [\[1e-05, 0.1\]],
  [yes],
  raw("GU__FED_f_absorption"),
  [fed effect on the fraction absorbed],
  [1.2 \[1.07, 1.35\]],
  [dimensionless],
  [\[1, 1.9\]],
  [yes],
  raw("GU__FED_f_dissolution"),
  [fed effect on the dissolution and gastric emptying],
  [0.171 \[0.082, 0.278\]],
  [dimensionless],
  [\[0.01, 1\]],
  [yes],
  raw("r_cyp3a5"),
  [capacity of CYP3A5 relative to CYP3A4],
  [0.458 \[LB, 1.18\]],
  [dimensionless],
  [\[0.01, 5\]],
  [no (LB)],
  ),
  caption: [Parameters of the fit PK\_20261017\_094949\_\_24b18, the parameter set PK\_20261017\_094949\_\_24b18\_0, cost 150, with their bounds in the optimization. The 95% confidence intervals are those of the profile likelihood, LB and UB an interval which is open towards the lower or the upper bound, and the Wald intervals of the Fisher information for a parameter without a profile. Identifiable is the classification of the profile: yes, no (LB) or no (UB) if the profile does not reach the threshold towards the lower or the upper bound, no if it reaches it on neither side, structural if it is flat.],
) <tab:parameters:PK_20261017_094949__24b18>

== Data

#figure(
  table(
  columns: 5,
  align: (left, right, right, right, right,),
  table.header([Study], [Training], [Validation], [Outlier], [Total]),
  [BasaDenes2019],
  [-],
  [5 (70)],
  [-],
  [5 (70)],
  [Bottiger2001],
  [1 (15)],
  [1 (15)],
  [-],
  [2 (30)],
  [Brattstrom2000],
  [5 (73)],
  [-],
  [-],
  [5 (73)],
  [Kelly1999],
  [-],
  [4 (32)],
  [-],
  [4 (32)],
  [KorthBradley2012],
  [3 (39)],
  [-],
  [-],
  [3 (39)],
  [Leelahavanichkul2005],
  [1 (11)],
  [-],
  [-],
  [1 (11)],
  [Leung2006],
  [2 (16)],
  [2 (20)],
  [-],
  [4 (36)],
  [Renders2007],
  [-],
  [2 (26)],
  [-],
  [2 (26)],
  [Rogers2008],
  [-],
  [4 (52)],
  [-],
  [4 (52)],
  [Tortorici2013],
  [1 (15)],
  [2 (30)],
  [-],
  [3 (45)],
  [Tortorici2014],
  [1 (14)],
  [1 (14)],
  [-],
  [2 (28)],
  [Wang2014],
  [-],
  [2 (26)],
  [-],
  [2 (26)],
  [Zha2022],
  [1 (22)],
  [1 (31)],
  [-],
  [2 (53)],
  [Zhang2017],
  [3 (42)],
  [-],
  [-],
  [3 (42)],
  [Zimmerman1997],
  [-],
  [5 (119)],
  [-],
  [5 (119)],
  [Zimmerman1999],
  [2 (32)],
  [-],
  [-],
  [2 (32)],
  [Zimmerman2003],
  [1 (13)],
  [2 (28)],
  [-],
  [3 (41)],
  [Zimmerman2008],
  [1 (16)],
  [1 (16)],
  [-],
  [2 (32)],
  [Total],
  [22 (308)],
  [32 (479)],
  [-],
  [54 (787)],
  ),
  caption: [Data of the fit PK\_20261017\_094949\_\_24b18: the number of fit mappings and, in parentheses, of data points per study and kind. The training data enters the cost, the validation data and the outliers are evaluated only.],
) <tab:data:PK_20261017_094949__24b18>

== Goodness of fit

#figure(
  table(
  columns: 11,
  align: (left, right, right, right, right, right, right, right, right, right, right,),
  table.header([Kind], [n], [cost], [MSE], [RMSE], [NRMSE], [RMSE\_w], [R²], [2-fold \[%\]], [AIC], [BIC]),
  [training],
  [308],
  [150],
  [5.58e-11],
  [7.47e-06],
  [0.572],
  [0.986],
  [0.86],
  [87.3],
  [-7253.5],
  [-7220.0],
  [validation],
  [479],
  [-],
  [5.17e-10],
  [2.27e-05],
  [0.758],
  [1.14],
  [0.694],
  [68.6],
  [-10224.4],
  [-10186.8],
  ),
  caption: [Goodness of fit of the fit PK\_20261017\_094949\_\_24b18 per kind of data, the parameter set PK\_20261017\_094949\_\_24b18\_0. The cost is defined on the training data only, 2-fold is the percentage of the data points with a positive value and prediction whose prediction is within twofold of the data.],
) <tab:metrics:PK_20261017_094949__24b18>

== Figures

#figure(
  image("plots/traces.svg", width: 100%),
  caption: [Cost of the optimizers over their steps, fit PK\_20261017\_094949\_\_24b18. The cost of every optimization run over its steps. Runs which end at the same cost found the same optimum; a run which stops much higher is stuck in a local minimum.],
) <fig:traces:PK_20261017_094949__24b18>

#figure(
  image("plots/waterfall.svg", width: 100%),
  caption: [Cost of the optimization runs, ordered, fit PK\_20261017\_094949\_\_24b18. The final cost of the runs, ordered. A flat plateau at the left is the global optimum found repeatedly, which is the evidence that the multistart converged; steps to the right are local minima.],
) <fig:waterfall:PK_20261017_094949__24b18>

#figure(
  image("plots/goodness_of_fit.svg", width: 100%),
  caption: [Prediction against the measured data points, per kind of fit mapping, with the agreement of the training data, the lines of twofold and the metrics of the panel, fit PK\_20261017\_094949\_\_24b18. Within twofold of the data are 269 of the 308 training data points (87.3%) and 327 of the 477 validation data points (68.6%); only the data points with a positive value and prediction are counted. Prediction against measurement, one point per data point, with a panel per kind of fit mapping. The points scatter around the diagonal when the model describes the data; a systematic deviation from it is a systematic error of the model. The band is the agreement of the training data, i.e. the same band the Bland-Altman plot draws, here as lines parallel to the diagonal. The dark red dotted lines are twofold, prediction \= 2 x measurement and prediction \= measurement \/ 2; the caption states the points between them of every kind. A point which is zero or negative is not drawn on the logarithmic axes and not counted, the title of a panel then states the points drawn of all. The validation panel shows how the fit describes data it was not fitted on, the outlier panel the data the fit dropped. The box of a panel carries its R², its normalized RMSE, the RMSE of the residuals of the cost and the percentage of its points within twofold.],
) <fig:goodness_of_fit:PK_20261017_094949__24b18>

#figure(
  image("plots/bland_altman.svg", width: 100%),
  caption: [Agreement of prediction and measurement as a ratio, per kind of fit mapping, with the agreement of the training data and the agreement of the panel, fit PK\_20261017\_094949\_\_24b18. The ratio of prediction and measurement over the geometric mean of the two, i.e. the goodness of fit with the diagonal turned into the horizontal. The band is the same, the bias and the limits of agreement #raw("bias ± 1.96 SD") of the training data as fold factors. A bias away from 1 is a systematic over- or underprediction, wide limits are a large scatter, and a trend over the mean is a model which describes the large or the small values better. The box of a panel carries the bias and the SD of its own points as fold factors and the share of them inside the limits of agreement of the training data.],
) <fig:bland_altman:PK_20261017_094949__24b18>

#figure(
  image("plots/profiles.svg", width: 100%),
  caption: [Profile likelihood of every parameter: the cost with the parameter fixed and the other parameters optimized, the threshold of the confidence level and the confidence interval, fit PK\_20261017\_094949\_\_24b18. One profile per parameter. A profile which crosses the dashed threshold on both sides gives a finite confidence interval; one which stays below it up to a bound is practically non-identifiable; a flat profile is structurally non-identifiable, the parameter is compensated by the others.],
) <fig:profiles:PK_20261017_094949__24b18>

== Provenance

#figure(
  table(
  columns: 2,
  align: (left, left,),
  table.header([Item], [Value]),
  [Fit],
  raw("PK_20261017_094949__24b18"),
  [Created],
  [2026-10-17 07:49:51 UTC],
  [Seed],
  [1234],
  [Python],
  [3.14.4],
  [sbmlsim],
  [0.10.0],
  [libRoadRunner],
  [2.10.0],
  [SciPy],
  [1.18.1],
  [NumPy],
  [2.5.3],
  [Repository],
  raw("/home/mkoenig/git/pkdb_models-rapamycin"),
  [Commit],
  raw("138cd5e26b655379e4c34971f8cb649aa6cb81cd"),
  [Uncommitted changes],
  [no],
  ),
  caption: [Provenance of the fit PK\_20261017\_094949\_\_24b18: the versions it ran with and the revision of its models and data; the SHA-256 checksums of the models and data are listed in #ref(<tab:checksums:PK_20261017_094949__24b18>).],
) <tab:provenance:PK_20261017_094949__24b18>

#bibliography("references.bib", style: "references.csl")
