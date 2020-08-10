Multi-modal Automated Survival Time Annotation
================

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
predicting event time with longitudinal encounter records. In step I, we
employ a functional principal component analysis approach to efficiently
approximate the intensity functions of individual encounters using the
observed point processes in the unlabeled set. In step II, we fit a
penalized proportional odds model to the event time outcomes with
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
library(MASTA)
```

    ## Loading required package: survival

## Documentation

  - [Input Data for the MASTA
    Algorithm](https://celehs.github.io/MASTA/articles/data.html)

  - [Step I. Feature Extraction with Functional
    PCA](https://celehs.github.io/MASTA/articles/step1.html)

  - [Step II. Model Estimation with
    B-Splines](https://celehs.github.io/MASTA/articles/step2.html)

  - [Acceleration with Multicore Parallel
    Computing](https://celehs.github.io/MASTA/articles/multicore.html)

## References

  - Liang, L., Uno, H., Ma, Y., Cai, T. **Robust Approach to Event Time
    Annotation Using Longitudinal Medical Encounters**. *Working Paper*.

  - Wu, S., Müller, H., Zhang, Z. (2013). **Functional Data Analysis for
    Point Processes with Rare Events**. *Statistica Sinica*, 23:1-23.
    <https://doi.org/10.5705/ss.2010.162>

## Work in Progress

``` r
data_org2
```

    ## $TrainLong
    ## $TrainLong[[1]]
    ##           id  fu_time month pred1
    ##     1:     4 14.85010     1     3
    ##     2:     4 14.85010     2     1
    ##     3:     5 80.65708     4     3
    ##     4:     5 80.65708     5     4
    ##     5:     5 80.65708     6     4
    ##    ---                           
    ## 72913: 20599 13.93018     8     1
    ## 72914: 20599 13.93018    12     1
    ## 72915: 20600 24.70637    12     2
    ## 72916: 20600 24.70637    13     1
    ## 72917: 20600 24.70637    16     1
    ## 
    ## $TrainLong[[2]]
    ##           id  fu_time month pred2
    ##     1:     1 49.41273    11     1
    ##     2:     5 80.65708     1     1
    ##     3:     5 80.65708     2     1
    ##     4:     5 80.65708     4     1
    ##     5:     6 42.64476     2     1
    ##    ---                           
    ## 57206: 20600 24.70637     1     1
    ## 57207: 20600 24.70637    14     1
    ## 57208: 20600 24.70637    16     1
    ## 57209: 20600 24.70637    18     1
    ## 57210: 20600 24.70637    19     1
    ## 
    ## $TrainLong[[3]]
    ##           id   fu_time month pred3
    ##     1:     6 42.644764    37     2
    ##     2:     6 42.644764    38     2
    ##     3:     6 42.644764    39     2
    ##     4:     6 42.644764    40     2
    ##     5:     6 42.644764    41     2
    ##    ---                            
    ## 13892: 20597 14.652977    14     1
    ## 13893: 20598  8.147844     8     2
    ## 13894: 20599 13.930185    13     2
    ## 13895: 20600 24.706366    23     1
    ## 13896: 20600 24.706366    24     2
    ## 
    ## 
    ## $ValidLong
    ## $ValidLong[[1]]
    ##          id  fu_time month pred1
    ##    1: 90001 71.72074     9     1
    ##    2: 90001 71.72074    52     1
    ##    3: 90001 71.72074    54     1
    ##    4: 90003 14.94867     5     2
    ##    5: 90003 14.94867     6     3
    ##   ---                           
    ## 1358: 90496 35.58111     2     3
    ## 1359: 90496 35.58111     3     3
    ## 1360: 90496 35.58111     4     2
    ## 1361: 90496 35.58111     5     3
    ## 1362: 90496 35.58111    34     1
    ## 
    ## $ValidLong[[2]]
    ##         id   fu_time month pred2
    ##   1: 90001  71.72074    44     1
    ##   2: 90002  70.37372    21     1
    ##   3: 90002  70.37372    30     1
    ##   4: 90003  14.94867     1     1
    ##   5: 90003  14.94867    10     1
    ##  ---                            
    ## 783: 90492 122.34908   121     1
    ## 784: 90494  76.78029    26     1
    ## 785: 90494  76.78029    38     1
    ## 786: 90498  24.60780    15     1
    ## 787: 90499  27.17043    21     1
    ## 
    ## $ValidLong[[3]]
    ##         id   fu_time month pred3
    ##   1: 90003 14.948665    14     1
    ##   2: 90010  7.622177     6     1
    ##   3: 90027 41.659138    38     2
    ##   4: 90027 41.659138    39     1
    ##   5: 90032 14.357290    13     2
    ##  ---                            
    ## 243: 90487 81.051335    77    17
    ## 244: 90487 81.051335    78    16
    ## 245: 90487 81.051335    79    15
    ## 246: 90487 81.051335    80    17
    ## 247: 90487 81.051335    81     1
