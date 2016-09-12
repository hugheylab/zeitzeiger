# zeitzeiger

`zeitzeiger` is a package for regularized supervised learning on high-dimensional data from an oscillatory system. `zeitzeiger` can make accurate predictions, identify major patterns and important features, and detect when the oscillator is perturbed.

For details about the method and to see how we used it to analyze circadian gene expression in mice, please check out [Hughey et al. (2016)](http://dx.doi.org/10.1093/nar/gkw030) and the [accompanying results](http://dx.doi.org/10.5061/dryad.hn8gp).

To see how we used `zeitzeiger` to analyze circadian gene expression from human blood, please check out [Hughey (2016)](http://dx.doi.org/10.1101/066126) and the [accompanying results](http://bit.ly/2a6M2HT).

## Installation
```R
source('https://bioconductor.org/biocLite.R')
install.packages('devtools')
devtools::install_github('jakejh/zeitzeiger', repos=BiocInstaller::biocinstallRepos(), build_vignettes=TRUE, dependencies=TRUE)
```

## Getting started
See `vignette('introduction', package='zeitzeiger')` and `vignette('association', package='zeitzeiger')`.
