# modeest 2.3.2

* BREAKING CHANGE: `mlv()` no longer returns a list, it returns a vector of 
values (usually one single value) for the sake of simplicity and homogeneity 
with functions such as `mean()` or `median()`. 
* `discrete()` is now removed. Use `mfv()` or `mfv1()` instead. 


# modeest 2.3.0 

* Hidden kernel related functions are transferred to package `statip`. 
* `discrete()` is now deprecated and will be removed in a future version of the 
package. Use `mfv()` instead. 
* `mfv1()` is a new function that always returns a length 1 value (so that `mfv1(x)==mfv(x)[[1L]]`). 


# modeest 2.2

* Thank you to W. H. Beasley who pointed out a slight mistake in the 
calculation of Bickel's skewness in `mlv.integer()`. Now the skewness is set at 
`NA` in case of multiple modes. 
* The meanshift mode estimator is added. 
* As documented under `?as.numeric`, the function `as.numeric.mlv()` was not 
correct, and is now replaced by `as.double.mlv()`. 


# modeest 2.1

* Thank you to C. Lepoittevin and K. Fijorek who pointed out a misuse of 
`ifelse` in the function `hsm()`. This has been corrected, so now `hsm` works 
correctly. 
* Functions `fiskMode()`, `gompertzMode()`, `koenkerMode()`, `kumarMode()`, 
`laplaceMode()`, `paralogisticMode()`, `paretoMode()`, `rayleighMode()` are 
added. 
* Function `symstbMode()` is removed. 
* Package modeest now needs to load `stabledist` and `VGAM` (in addition to 
`stats`, `evd`, and `fBasics`). 
  

# modeest 1.09

* The Asselin de Beauville mode estimator has been added. 
* A function `as.numeric.mlv()` is provided. 
* Methods for the Chernoff distribution are provisionally suppressed, because 
of their lack of efficiency. 
* The DIP statistic is no longer provided. 
* In function `tsybakov()`, the argument `djeddour` is renamed `dmp`. 
* In functions `parzen()` and `mlv.density()`, the argument `biau` is renamed `abc`. 
