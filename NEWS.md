# modeest 2.4.0 (2019-11-17)

## Breaking changes

* Defunct `hrm()` due to difficulties in installing automatically Bioconductor's 'genefilter' dependency. This function shall be reintroduced in a future version.

## Performance

* Remove dependency to the 'genefilter' package.


# modeest 2.3.3

## Features

* Simplify and improve the code of `parzen()`, `vieu()`, `tsybakov()`, `meanshift()`, `mlv()`. 

* Improve documentation. 


# modeest 2.3.2

## Breaking changes

* `mlv()` no longer returns a list, it returns a vector of 
values (usually one single value) for the sake of simplicity and homogeneity with functions such as `mean()` or `median()`. 

* Move `mfv()` and `mfv1()` to package `statip` (and reexport them).

* Remove `discrete()`. Use `mfv()` or `mfv1()` instead. 


# modeest 2.3.0 

## Features

* Move hidden kernel related functions to package `statip`. 

* Deprecate `discrete()`; it will be removed in a future version of the package. Start using `mfv()` instead. 

* Add `mfv1()`, which always returns a length 1 value (so that `mfv1(x)==mfv(x)[[1L]]`). 


# modeest 2.2

## Bug fixes

* Thank you to W. H. Beasley who pointed out a slight mistake in the 
calculation of Bickel's skewness in `mlv.integer()`. Now the skewness is set at `NA` in case of multiple modes. 

* As documented under `?as.numeric`, the function `as.numeric.mlv()` was not correct, and is now replaced by `as.double.mlv()`. 

## Features

* Add the meanshift mode estimator. 


# modeest 2.1

## Breaking changes

* Remove function `symstbMode()`. 

## Bug fixes

* Thank you to C. Lepoittevin and K. Fijorek who pointed out a misuse of 
`ifelse` in the function `hsm()`. This has been corrected, so now `hsm` works correctly.

## Features

* Add functions `fiskMode()`, `gompertzMode()`, `koenkerMode()`, `kumarMode()`, 
`laplaceMode()`, `paralogisticMode()`, `paretoMode()`, `rayleighMode()`. 


# modeest 1.09

## Breaking changes

* Remove methods for the Chernoff distribution, because of their lack of efficiency. 

* Remove the DIP statistic. 

* In function `tsybakov()`, the argument `djeddour` is renamed `dmp`. 

* In functions `parzen()` and `mlv.density()`, the argument `biau` is renamed `abc`. 

## Features

* Add the Asselin de Beauville mode estimator. 

* Add function `as.numeric.mlv()`. 
