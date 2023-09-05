# medNet

## A New Framework for Exploratory Network Mediator Analysis in Omics Data

![Setup of medNet](Overview.png)
We propose a new framework for predictive mediation analysis that can handle tens of thousands can-
didate mediators, model nonlinear relations, and incorporate existing knowledge. The goal of our method is to sensitively detect potential mediators and mediator networks, which can be validated by further experiments.Instead of using the traditional mediation model which is very rigid, we focus on the prediction accuracy estimated by cross-validation as a criterion for mediator selection. medNet is implemented as an open-source R package.

Installation
------------
You can install the released version of medNet from Github with the following code.

## Dependencies 
* R version >= 4.1.0.
* R packages: grDevices, igraph, abind, caret, doParallel,
        e1071, energy, foreach, igraph, parallel, pROC, randomForest, foreach, doParallel, parallel
        
``` r
# install devtools if necessary
install.packages('devtools')

# install the CARD package
devtools::install_github('EddieFua/medNet')

# load package
library(medNet)

```
