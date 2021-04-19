Sample Data
================

``` r
library(tidyverse)
```

    ## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.0 ──

    ## ✓ ggplot2 3.3.3     ✓ purrr   0.3.4
    ## ✓ tibble  3.0.5     ✓ dplyr   1.0.3
    ## ✓ tidyr   1.1.2     ✓ stringr 1.4.0
    ## ✓ readr   1.4.0     ✓ forcats 0.5.0

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## x tidyr::chop()   masks uno1misc::chop()
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()

``` r
TrainCode <- read_csv("data_org/TrainCode.csv")
```

    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   id = col_double(),
    ##   fu_time = col_double(),
    ##   month = col_double(),
    ##   pred1 = col_double(),
    ##   pred2 = col_double(),
    ##   pred3 = col_double()
    ## )

``` r
ValidCode <- read_csv("data_org/ValidCode.csv")
```

    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   id = col_double(),
    ##   fu_time = col_double(),
    ##   month = col_double(),
    ##   pred1 = col_double(),
    ##   pred2 = col_double(),
    ##   pred3 = col_double()
    ## )

``` r
TrainSurv <- read_csv("data_org/TrainSurv.csv")
```

    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   id = col_double(),
    ##   event_ind = col_double(),
    ##   event_time = col_double(),
    ##   fu_time = col_double(),
    ##   cov_1 = col_double(),
    ##   cov_2 = col_double(),
    ##   cov_3 = col_double()
    ## )

``` r
ValidSurv <- read_csv("data_org/ValidSurv.csv")
```

    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   id = col_double(),
    ##   event_ind = col_double(),
    ##   event_time = col_double(),
    ##   fu_time = col_double(),
    ##   cov_1 = col_double(),
    ##   cov_2 = col_double(),
    ##   cov_3 = col_double()
    ## )

``` r
f <- function(data, var) {
  data2 <- data %>% 
    mutate(code = substr(var, 5, 5)) %>%
    filter(.data[[var]] > 0) %>%
    select(code, id, time = month, .data[[var]]) %>%  
    slice(rep(1:n(), times = .data[[var]])) %>% 
    select(-.data[[var]])
  return(data2)
}
```

``` r
follow_up_time <- bind_rows(
  bind_cols(TrainCode[, 1:2], train_valid = 1), 
  bind_cols(ValidCode[, 1:2], train_valid = 2)) %>%
  distinct() 
follow_up_time
```

    ## # A tibble: 21,100 x 3
    ##       id fu_time train_valid
    ##    <dbl>   <dbl>       <dbl>
    ##  1     1   49.4            1
    ##  2     2   13.9            1
    ##  3     3   12.6            1
    ##  4     4   14.9            1
    ##  5     5   80.7            1
    ##  6     6   42.6            1
    ##  7     7   13.7            1
    ##  8     8   21.0            1
    ##  9     9    5.16           1
    ## 10    10    6.60           1
    ## # … with 21,090 more rows

``` r
longitudinal <- bind_rows(
  f(TrainCode, "pred1"),
  f(TrainCode, "pred2"),
  f(TrainCode, "pred3"),  
  f(ValidCode, "pred1"),
  f(ValidCode, "pred2"),
  f(ValidCode, "pred3")) %>%
  arrange(code, id, time)
longitudinal
```

    ## # A tibble: 258,117 x 3
    ##    code     id  time
    ##    <chr> <dbl> <dbl>
    ##  1 1         4     1
    ##  2 1         4     1
    ##  3 1         4     1
    ##  4 1         4     2
    ##  5 1         5     4
    ##  6 1         5     4
    ##  7 1         5     4
    ##  8 1         5     5
    ##  9 1         5     5
    ## 10 1         5     5
    ## # … with 258,107 more rows

``` r
survival <- bind_rows(
  TrainSurv[,-4], 
  ValidSurv[,-4])
survival
```

    ## # A tibble: 1,100 x 6
    ##       id event_ind event_time cov_1 cov_2 cov_3
    ##    <dbl>     <dbl>      <dbl> <dbl> <dbl> <dbl>
    ##  1     1         1       9.36    79     1     0
    ##  2     2         0      13.9     81     0     0
    ##  3     3         0      12.6     55     1     1
    ##  4     4         0      14.9     72     1     0
    ##  5     5         0      80.7     83     1     1
    ##  6     6         1      15.7     47     1     0
    ##  7     7         0      13.7     86     0     1
    ##  8     8         0      21.0     40     0     0
    ##  9     9         0       5.16    85     1     0
    ## 10    10         0       6.60    58     0     1
    ## # … with 1,090 more rows

``` r
write_csv(follow_up_time, "follow_up_time.csv")
write_csv(longitudinal, "longitudinal.csv")
write_csv(survival, "survival.csv")
```

``` r
usethis::use_data(follow_up_time, overwrite = TRUE)
```

    ## ✓ Setting active project to '/Users/hajime/GitHub/MASTA'

    ## ✓ Saving 'follow_up_time' to 'data/follow_up_time.rda'

    ## ● Document your data (see 'https://r-pkgs.org/data.html')

``` r
usethis::use_data(longitudinal, overwrite = TRUE)
```

    ## ✓ Saving 'longitudinal' to 'data/longitudinal.rda'
    ## ● Document your data (see 'https://r-pkgs.org/data.html')

``` r
usethis::use_data(survival, overwrite = TRUE)
```

    ## ✓ Saving 'survival' to 'data/survival.rda'
    ## ● Document your data (see 'https://r-pkgs.org/data.html')

``` r
proc.time()
```

    ##    user  system elapsed 
    ##   3.400   0.225   3.740
