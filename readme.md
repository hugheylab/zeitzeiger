# zeitzeiger

`zeitzeiger` is a package for regularized supervised learning on high-dimensional data from an oscillatory system. `zeitzeiger` can make accurate predictions, identify major patterns and important features, and detect when the oscillator is perturbed. For more details about the method and how we used it to analyze circadian gene expression, please see [Hughey et al. (2016)](https://isitchristmas.com/).

## Installation
```R
install.packages('devtools')
devtools::install_github('hadley/ggplot2') # used in the vignette
devtools::install_github('jakejh/zeitzeiger', build_vignettes=TRUE, dependencies=TRUE)
```

## Getting started
See `vignette('introduction', package='zeitzeiger')`.
