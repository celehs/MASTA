## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE)

## ---- echo=FALSE, fig.width=5, fig.height=3-----------------------------------
DiagrammeR::grViz("
digraph flowchart {
  node [fontname = Helvetica, shape = rectangle]        
  tab1 [label = '@@1']
  tab2 [label = '@@2']
  tab3 [label = '@@3']
  tab4 [label = '@@4']
  tab5 [label = '@@5']
  tab1 -> tab2;
  tab1 -> tab3;
  tab2 -> tab4;
  tab2 -> tab5;
}
[1]: 'Input Data'
[2]: 'Labeled (n = 1,100)'
[3]: 'Unlabeled (N = 20,000)'
[4]: 'Training (600)'
[5]: 'Validation (500)'
")

## ---- eval=TRUE---------------------------------------------------------------
library(MASTA)

## ---- eval=FALSE--------------------------------------------------------------
#  ?longitudinal

## -----------------------------------------------------------------------------
head(longitudinal)
table(longitudinal$code)

## ---- eval=FALSE--------------------------------------------------------------
#  ?follow_up_time

## -----------------------------------------------------------------------------
head(follow_up_time)

## ---- eval=FALSE--------------------------------------------------------------
#  ?survival

## -----------------------------------------------------------------------------
head(survival)

