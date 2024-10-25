{smcl}
{* *! version 0.1, 1 July 2024}{...}
{cmd:help for dmlpath}{right:Geoffrey T. Wodtke}
{hline}

{title:Title}

{p2colset 5 18 18 2}{...}
{p2col : {cmd:dmlpath} {hline 2}} analysis of path-specific effects using de-biased machine learning {p_end}
{p2colreset}{...}


{title:Syntax}

{p 8 18 2}
{cmd:dmlpath} {depvar} {help indepvars:mvars} {ifin} {cmd:,} 
{opt model(string)}
{opt dvar(varname)} 
{opt d(real)} 
{opt dstar(real)} 
{opt cvars(varlist)} 
{opt xfits(integer)} 
{opt seed(integer)}
{opt censor(numlist)}
[{it:{help rforest##options:rforest_options}}]
[{it:{help lasso2##options:lasso2_options}}]

{phang}{opt model(string)} - this specifies which machine learning algorithm to implement. Options are rforest and lasso. 
For model(rforest), random forests are used to predict the nuisance terms. 
For model(lasso), LASSO models with all two-way interactions are used to predict the nuisance terms.

{phang}{opt depvar} - this specifies the outcome variable.

{phang}{opt mvars} - this specifies the mediators in causal order, beggining with the first in the hypothesized causal sequence
and ending with the last. Up to 5 causally ordered mediators are permitted.

{phang}{opt dvar(varname)} - this specifies the treatment (exposure) variable, which must be binary and coded 0/1.

{phang}{opt d(real)} - this specifies the reference level of treatment.

{phang}{opt dstar(real)} - this specifies the alternative level of treatment. Together, (d - dstar) defines
the treatment contrast of interest.

{title:Options}

{phang}{opt cvars(varlist)} - this option specifies the list of baseline covariates to be included in the analysis. Categorical 
variables need to be coded as a series of dummy variables before being entered as covariates.

{phang}{opt xfits(integer)} - this option specifies the number of sample partitions to use for cross-fitting (the default is 5).

{phang}{opt seed(integer)} - this option specifies the seed for cross-fitting and model training.

{phang}{opt censor(numlist)} - this option specifies that the inverse probability weights used in the robust estimating equations are censored at the percentiles supplied in {numlist}. For example,
censor(1 99) censors the weights at their 1st and 99th percentiles.

{phang}{it:{help rforest##options:rforest_options}} - all {help rforest} options are available when using model(rforest). 

{phang}{it:{help lasso2##options:lasso2_options}} - all {help lasso2} options are available when using model(lasso). {p_end}

{title:Description}

{pstd}{cmd:dmlpath} performs analysis of path-specific effects using de-biased machine learning (DML), and it computes inferential statistics using
analytic standard errors and a normal approximation. Currently, the command supports implementation of DML with random forests and LASSO models. 
It requires prior installation of the {cmd:rforest} and {cmd:lassopack} modules. {p_end}

{pstd}More specifically, {cmd:dmlpath} estimates path-specific effects using the type mr2 esimtator, as described in Chapter 6, Section 6.4, 
of Wodtke and Zhou "Causal Mediation Analysis." If there are K causally ordered mediators, {cmd:dmlpath} provides estimates for the total effect 
and then for K+1 path-specific effects: the direct effect of the exposure on the outcome that does not operate through any of the mediators, 
and then a separate path-specific effect operating through each of the K mediators, net of the mediators that precede it in causal order. If 
only a single mediator is specified, {cmd:dmlpath} computes and reports DML estimates of conventional natural direct and indirect 
effects through a univariate mediator. {p_end}

{title:Examples}

{pstd}Setup{p_end}
{phang2}{cmd:. use nlsy79.dta} {p_end}

{pstd} K=2 causally ordered mediators with the LASSO and weights censored at their 1st and 99th percentiles: {p_end}
 
{phang2}{cmd:. dmlpath std_cesd_age40 ever_unemp_age3539 log_faminc_adj_age3539, model(lasso) dvar(att22) cvars(female black hispan paredu parprof parinc_prank famsize afqt3) d(1) dstar(0) censor(1 99)} {p_end}

{pstd} K=3 causally ordered mediators with the LASSO and censored weights: {p_end}
 
{phang2}{cmd:. dmlpath std_cesd_age40 cesd_1992 ever_unemp_age3539 log_faminc_adj_age3539, model(lasso) dvar(att22) cvars(female black hispan paredu parprof parinc_prank famsize afqt3) d(1) dstar(0) censor(1 99)} {p_end}

{pstd} K=2 causally ordered mediators with random forests and custom hyperparameters: {p_end}
 
{phang2}{cmd:. dmlpath std_cesd_age40 ever_unemp_age3539 log_faminc_adj_age3539, model(rforest) iter(200) lsize(20) dvar(att22) cvars(female black hispan paredu parprof parinc_prank famsize afqt3) d(1) dstar(0) censor(1 99)} {p_end}

{title:Saved results}

{pstd}{cmd:dmlpath} saves the following results in {cmd:e()}:

{synoptset 15 tabbed}{...}
{p2col 5 15 19 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}matrix containing the total and path-specific effect estimates{p_end}


{title:Author}

{pstd}Geoffrey T. Wodtke {break}
Department of Sociology{break}
University of Chicago{p_end}

{phang}Email: wodtke@uchicago.edu


{title:References}

{pstd}Wodtke, GT and X Zhou. Causal Mediation Analysis. In preparation. {p_end}

{title:Also see}

{psee}
Help: {manhelp regress R}, {manhelp logit R}, {manhelp bootstrap R}
{p_end}
