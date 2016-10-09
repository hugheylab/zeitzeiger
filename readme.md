# zeitzeiger

`zeitzeiger` is a package for regularized supervised learning on high-dimensional data from an oscillatory system. `zeitzeiger` can quantify rhythmic behavior, make accurate predictions, identify major patterns and important features, and detect when the oscillator is perturbed.

For details about the method and to see how we used it to analyze circadian gene expression in mice, check out [Hughey et al. (2016)](http://dx.doi.org/10.1093/nar/gkw030) and the [accompanying results](http://dx.doi.org/10.5061/dryad.hn8gp).

To see how we used `zeitzeiger` to analyze the phasing of circadian clocks in humans and other mammals, check out [Hughey and Butte (2016)](http://dx.doi.org/10.1177/0748730416668049) and the [accompanying results]( http://dx.doi.org/10.5061/dryad.g928q).

To see how we used `zeitzeiger` to predict circadian time from gene expression from human blood, check out [Hughey (2016)](http://dx.doi.org/10.1101/066126) and the [accompanying results](http://bit.ly/2a6M2HT).

## Installation
Install the package in your local version of R:
```R
source('https://bioconductor.org/biocLite.R')
install.packages('devtools')
devtools::install_github('jakejh/zeitzeiger', repos=BiocInstaller::biocinstallRepos(), build_vignettes=TRUE, dependencies=TRUE)
```

Or use a pre-built [docker image](https://hub.docker.com/r/jakejh/hugheyverse), which has all dependencies already installed:
```
docker pull jakejh/hugheyverse
```

## Getting started
See `vignette('introduction', package='zeitzeiger')` and `vignette('association', package='zeitzeiger')`.
