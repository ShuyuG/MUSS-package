# MUSS-package
Spike and Slab Variable Selector under Matrix Uncertainty

In high-dimensional sparse regression with measurement error in variables, assuming the errors are Normally distributed with mean 0 and errors for each variable have a specific variance, then 'MUSS' can select variables for such cases. 'MUSS' selects variables by spike and slab priors for the regression coefficients. It works under EM framework, where the unknown true design matrix and indicators of spike-slab priors are treated as latent variables. To avoid the effect of inappropriate choices of spike and slab scale parameters on variable selection, the final output coefficients are obtained following a decreasing sequence of spike parameters and the path can be displayed by functions from 'MUSS' package.

## Installation
'MUSS' package can be installed by

```devtools::install_github("ShuyuG/MUSS-package")``` 

## Vignette
[Vignette](https://github.com/ShuyuG/MUSS-package/blob/master/inst/doc/Vignette.pdf)

## Manual
[Manual](https://github.com/ShuyuG/MUSS-package/blob/master/inst/doc/MUSS_1.0.0.pdf)

## Siulation Results


## Acknowledgement
Thanks [Bo Y.-C. Ning](https://bo-ning.github.io/) for contributing original idea of 'MUSS' and discussion.
