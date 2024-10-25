# dmlpath: A Stata Module for Analysis of Path-Specific Effects using De-Biased Machine Learning

## Overview

**dmlpath** is a Stata module designed to estimate path-specific effects using de-biased machine learning (DML). It provides estimates for the total effect and path-specific effects operating through multiple mediators in a causal sequence, leveraging machine learning algorithms for improved robustness.

## Syntax

```stata
dmlpath depvar mvars [if] [in], model(string) dvar(varname) d(real) dstar(real) [options]
```

### Required Arguments

- **depvar**: Specifies the outcome variable.
- **mvars**: Specifies the mediators in causal order, starting with the first in the hypothesized causal sequence and ending with the last. Up to 5 causally ordered mediators are allowed.
- **model(string)**: Specifies the machine learning algorithm to use for estimating nuisance functions. Options are `rforest` (random forests) and `lasso` (LASSO with two-way interactions).
- **dvar(varname)**: Specifies the treatment (exposure) variable, which must be binary (0/1).
- **d(real)**: Specifies the reference level of treatment.
- **dstar(real)**: Specifies the alternative level of treatment. Together, `d` and `dstar` define the treatment contrast of interest.

### Options

- **cvars(varlist)**: Specifies the list of baseline covariates to be included in the analysis. Categorical variables must be coded as dummy variables.
- **xfits(integer)**: Specifies the number of sample partitions to use for cross-fitting (default is 5).
- **seed(integer)**: Sets the seed for cross-fitting and model training.
- **censor(numlist)**: Specifies that the inverse probability weights used in the robust estimating equations are censored at the percentiles provided in `numlist`.
- **rforest_options**: All options available in the `rforest` module can be specified when using `model(rforest)`.
- **lasso2_options**: All options available in the `lasso2` module can be specified when using `model(lasso)`.

## Description

`dmlpath` performs analysis of path-specific effects using de-biased machine learning (DML) and computes inferential statistics using analytic standard errors and a normal approximation. The command supports DML implementations using either random forests or LASSO models and requires prior installation of the **rforest** and **lassopack** modules.

More specifically, `dmlpath` estimates path-specific effects using the type MR2 estimator, as described in Chapter 6, Section 6.4, of Wodtke and Zhou's *Causal Mediation Analysis*. If there are K causally ordered mediators, `dmlpath` provides estimates for the total effect and K+1 path-specific effects: the direct effect of the exposure on the outcome that does not operate through any of the mediators, and a separate path-specific effect operating through each of the K mediators, net of the mediators that precede it in causal order. If only a single mediator is specified, `dmlpath` computes and reports DML estimates of conventional natural direct and indirect effects through a univariate mediator.

## Examples

### Example 1: K=2 causally ordered mediators with LASSO and weights censored at 1st and 99th percentiles

```stata
use nlsy79.dta
dmlpath std_cesd_age40 ever_unemp_age3539 log_faminc_adj_age3539, model(lasso) dvar(att22) cvars(female black hispan paredu parprof parinc_prank famsize afqt3) d(1) dstar(0) censor(1 99)
```

### Example 2: K=3 causally ordered mediators with LASSO and censored weights

```stata
dmlpath std_cesd_age40 cesd_1992 ever_unemp_age3539 log_faminc_adj_age3539, model(lasso) dvar(att22) cvars(female black hispan paredu parprof parinc_prank famsize afqt3) d(1) dstar(0) censor(1 99)
```

### Example 3: K=2 causally ordered mediators with random forests and custom hyperparameters

```stata
dmlpath std_cesd_age40 ever_unemp_age3539 log_faminc_adj_age3539, model(rforest) iter(200) lsize(20) dvar(att22) cvars(female black hispan paredu parprof parinc_prank famsize afqt3) d(1) dstar(0) censor(1 99)
```

## Saved Results

The following results are saved in `e()`:

- **Matrices:**
  - `e(b)`: Matrix containing the total and path-specific effect estimates.

## Author

**Geoffrey T. Wodtke**  
Department of Sociology  
University of Chicago  
Email: [wodtke@uchicago.edu](mailto:wodtke@uchicago.edu)

## References

- Wodtke, GT, and X. Zhou. *Causal Mediation Analysis*. In preparation.

## See Also

- `rforest`, `lasso2`, `lassologit`
  
