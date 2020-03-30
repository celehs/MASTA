**Predicting Event Time with Longitudinal Encounter Records**

Overview
--------

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
predicting event time with longitudinal encounter records (PETLER). In
step I, we employ a functional principal component analysis approach to
efficiently approximate the intensity functions of individual encounters
using the observed point processes in the unlabeled set. In step II, we
fit a penalized proportional odds model to the event time outcomes with
features derived in step I in the labeled data where the non-parametric
baseline function is approximated using B-splines.

Installation
------------

Install development version from GitHub.

``` r
# install.packages("remotes")
remotes::install_github("celehs/PETLER")
```

Load the package into R.

``` r
library(PETLER)
```

    ## Loading required package: survival
    
Getting Started
---------------

- Data <https://celehs.github.io/PETLER/articles/data.html>

- Step I <https://celehs.github.io/PETLER/articles/step1.html>

- Step II <https://celehs.github.io/PETLER/articles/step2.html>

- Multicore <https://celehs.github.io/PETLER/articles/scalability.html>

- Reference <https://celehs.github.io/PETLER/reference/index.html>

References
----------

-   Liang, L., Uno, H., Ma, Y., Cai, T. **Robust Approach to Event Time
    Annotation Using Longitudinal Medical Encounters**. *Working Paper*.

-   Wu, S., Müller, H., Zhang, Z. (2013). **Functional Data Analysis for
    Point Processes with Rare Events**. *Statistica Sinica*, 23:1-23.
    <a href="https://doi.org/10.5705/ss.2010.162" class="uri">https://doi.org/10.5705/ss.2010.162</a>
