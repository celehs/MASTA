## Overview

Large clinical datasets derived from insurance claims and electronic health record (EHR) systems are valuable sources for precision medicine research. These datasets can be used to develop models for personalized prediction of risk or treatment response. Efficiently deriving
prediction models using real world data, however, faces practical and methodological challenges. For example, event time information such as time to cancer progression is not readily available. Simple estimates of the event time based on billing or procedure codes may poorly approximate the true event time.

This package implements a two-step semi-supervised learning method on predicting event time with longitudinal encounter records (PETLER). In step I, we employ a functional principal component analysis approach to efficiently approximate the intensity functions of individual encounters using the observed point processes in the unlabeled set. In step II, we fit a penalized proportional odds model to the event time outcomes with features derived in step I in the labeled data where the non-parametric baseline function is approximated using B-splines. 

## Installation

If `devtools` is not installed, uncomment the code below and install it from CRAN.

``` r
# install.packages("devtools")
```

Install development version from GitHub:

``` r
devtools::install_github("celehs/PETLER")
```

## Getting Started

WORK IN PROGRESS...

- [STEP 1](https://celehs.github.io/PETLER/demo/step1.html): Feature Engineering

- [STEP 2](https://celehs.github.io/PETLER/demo/step2.html): Model Training & Validation

The data used for this demo can be found [HERE](https://github.com/celehs/PETLER/tree/master/demo). 

Click [HERE](https://celehs.github.io/PETLER/reference/petler.html) to view the documentation for the main function `petler`.

## References

- Liang, L., Uno, H., Ma, Y., Cai, T. __Robust Approach to Event Time Annotation
Using Longitudinal Medical Encounters__. _Working Paper_.

- Wu, S., Müller, H., Zhang, Z. (2013). __Functional Data Analysis for Point Processes with Rare Events__. _Statistica Sinica_, 23:1-23. <https://doi.org/10.5705/ss.2010.162>
