# zeitzeiger
`zeitzeiger` is a package for regularized supervised learning on high-dimensional data from an oscillatory system. `zeitzeiger` can quantify rhythmic behavior, make accurate predictions, identify major patterns and important features, and detect when the oscillator is perturbed.

For details about the method and to see how we used it to analyze circadian gene expression in mice, check out [Hughey et al. (2016)](https://dx.doi.org/10.1093/nar/gkw030) and the [accompanying results](https://dx.doi.org/10.5061/dryad.hn8gp).

To see how we used `zeitzeiger` to analyze the phasing of circadian clocks in humans and other mammals, check out [Hughey and Butte (2016)](https://dx.doi.org/10.1177/0748730416668049) and the [accompanying results](https://dx.doi.org/10.5061/dryad.g928q).

To see how we used `zeitzeiger` to predict circadian time from gene expression in human blood, check out [Hughey (2017)](https://dx.doi.org/10.1186/s13073-017-0406-4) and the [accompanying results](https://dx.doi.org/10.6084/m9.figshare.3756375.v1).

## Install using drat
```R
install.packages('drat')
setRepositories(ind=1:5)
drat::addRepo('hugheylab')
install.packages('zeitzeiger', type='source')
```
You can update the package by calling `drat::addRepo('hugheylab')`, then `update.packages`.

## Install using devtools
```R
install.packages('devtools')
setRepositories(ind=1:5)
devtools::install_github('hugheylab/zeitzeiger')
```
You can update the package using these same three lines.

## Install using docker
You can use a pre-built [docker image](https://hub.docker.com/r/hugheylab/hugheyverse), which has all dependencies already installed:
```
docker pull hugheylab/hugheyverse
```

## Getting started
See `vignette('introduction', package='zeitzeiger')` and `vignette('association', package='zeitzeiger')`.
