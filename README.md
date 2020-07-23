__Multi-modal Automated Survival Time Annotation (MASTA) Algorithm with Longitudinal Electronic Health Records__

## Overview

Large clinical datasets derived from insurance claims and electronic
health record (EHR) systems are valuable sources for precision medicine
research. These datasets can be used to develop models for personalized
prediction of risk or treatment response. Efficiently deriving
prediction models using real world data, however, faces practical and
methodological challenges. For example, event time information such as
time to cancer progression is not readily available. Simple estimates of
the event time based on billing or procedure codes may poorly
approximate the true event time.

This package implements a two-step semi-supervised learning method on
predicting event time with longitudinal encounter records. In
step I, we employ a functional principal component analysis approach to
efficiently approximate the intensity functions of individual encounters
using the observed point processes in the unlabeled set. In step II, we
fit a penalized proportional odds model to the event time outcomes with
features derived in step I in the labeled data where the non-parametric
baseline function is approximated using B-splines.

![](https://github.com/celehs/MASTA/raw/master/flowchart/flowchart.png)

## Installation

Install development version from GitHub.

``` r
# install.packages("remotes")
remotes::install_github("celehs/MASTA")
```

Load the package into R.

``` r
library(MATA)
```

## Documentation

- [Input Data for the MASTA Algorithm](https://celehs.github.io/MASTA/articles/data.html)

- [Step I. Feature Extraction with Functional PCA](https://celehs.github.io/MASTA/articles/step1.html)

- [Step II. Model Estimation with B-Splines](https://celehs.github.io/MASTA/articles/step2.html)

- [Acceleration with Multicore Parallel Computing](https://celehs.github.io/MASTA/articles/multicore.html)

## References

- Liang, L., Uno, H., Ma, Y., Cai, T. __Robust Approach to Event Time Annotation Using Longitudinal Medical Encounters__. _Working Paper_.

- Wu, S., Müller, H., Zhang, Z. (2013). __Functional Data Analysis for Point Processes with Rare Events__. _Statistica Sinica_, 23:1-23. <https://doi.org/10.5705/ss.2010.162>
