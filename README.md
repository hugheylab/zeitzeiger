# zeitzeiger

[![check-deploy](https://github.com/hugheylab/zeitzeiger/workflows/check-deploy/badge.svg)](https://github.com/hugheylab/zeitzeiger/actions)
[![Netlify Status](https://api.netlify.com/api/v1/badges/80a704ab-1728-4e6a-a17d-a7cce4e6f198/deploy-status)](https://app.netlify.com/sites/dazzling-wilson-37374e/deploys)
[![drat version](https://raw.githubusercontent.com/hugheylab/drat/gh-pages/badges/zeitzeiger_drat_badge.svg)](https://github.com/hugheylab/drat/tree/gh-pages/src/contrib)

`zeitzeiger` is a package for regularized supervised learning on high-dimensional data from an oscillatory system. `zeitzeiger` can quantify rhythmic behavior, make accurate predictions, identify major patterns and important features, and detect when the oscillator is perturbed.

Update (Nov 2018): ZeitZeiger now uses [limma](https://doi.org/doi:10.18129/B9.bioc.limma) internally, which makes training (previously the slowest step by far) about 250x faster. Additional optimizations have made calculation of the sparse principal components about 6x faster, and prediction about 20% faster.

For details about the method and to see how we used it to analyze circadian gene expression in mice, check out [Hughey et al. (2016)](https://doi.org/10.1093/nar/gkw030) and the [accompanying results](https://doi.org/10.5061/dryad.hn8gp).

To see how we used `zeitzeiger` to analyze the phasing of circadian clocks in humans and other mammals, check out [Hughey and Butte (2016)](https://doi.org/10.1177/0748730416668049) and the [accompanying results](https://doi.org/10.5061/dryad.g928q).

To see how we used `zeitzeiger` to predict circadian time from gene expression in human blood, check out [Hughey (2017)](https://doi.org/10.1186/s13073-017-0406-4) and the [accompanying results](https://doi.org/10.6084/m9.figshare.3756375.v1).

## Installation

1. Install [`BiocManager`](https://cran.r-project.org/package=BiocManager).

    ```r
    if (!requireNamespace('BiocManager', quietly = TRUE))
      install.packages('BiocManager')
    ```

1. If you use RStudio, go to Tools → Global Options... → Packages → Add... (under Secondary repositories), then enter:

    - Name: hugheylab
    - Url: https://hugheylab.github.io/drat/

    You only have to do this once. Then you can install or update the package by entering:

    ```r
    BiocManager::install('zeitzeiger')
    ```

    Alternatively, you can install or update the package by entering:

    ```r
    BiocManager::install('zeitzeiger', site_repository = 'https://hugheylab.github.io/drat/')
    ```

## Usage

For an introduction to the package, read the [vignette](https://zeitzeiger.hugheylab.org/articles/introduction.html). For detailed help on specific functions, check out the [reference documentation](https://zeitzeiger.hugheylab.org/reference/index.html).
