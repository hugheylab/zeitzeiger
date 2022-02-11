## R CMD check results

### Local check
`devtools::check()` result:

  0 errors ✓ | 0 warnings ✓ | 0 notes ✓

### Online check
`devtools::check_rhub()` result:

  > checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Jake Hughey <jakejhughey@gmail.com>'
  
  New submission
  
  Possibly misspelled words in DESCRIPTION:
    Hughey (8:53, 9:32)
    al (8:63)
    et (8:60)

> checking for detritus in the temp directory ... NOTE
  Found the following files/directories:
    'lastMiKTeXException'

  0 errors ✓ | 0 warnings ✓ | 2 notes x

Notes:
  - This is a the first time this package was sent to CRAN, and the words are names and correctly spelled terms, so ignore this.
  - The second note only occurs on the Windows Server rhub environment, and from what I have seen about these types of notes they do not occur when building and checking on CRAN.

You can also see the results of R CMD check on Windows, Linux, and MacOS by going to the GitHub Actions run labeled `check-deploy` [here](https://github.com/hugheylab/zeitzeiger/actions).

## Downstream dependencies
There are no downstream dependencies for zeitzeiger.
