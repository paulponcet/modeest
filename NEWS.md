# modeest 2.5.0 (2026-09-20)

* Lint various parts of the code.


# modeest 2.3.4 (2019-07-15)

* Fix bug in `betaMode()`. Results were incorrect when `shape1` or `shape2` was 
less than `1`. 


# modeest 2.3.3

* Update documentation. 

* Simplify and improve the code of `parzen()`, `vieu()`, 
`tsybakov()`, `meanshift()`, `mlv()`. 


# modeest 2.3.2

* BREAKING CHANGE: `mlv()` no longer returns a list, it returns a vector of 
values (usually one single value) for the sake of simplicity and homogeneity 
with functions such as `mean()` or `median()`. 

* Move `mfv()` and `mfv1()` to package `statip` and reexport then.

* Remove `discrete()`. Use `mfv()` or `mfv1()` instead. 


# modeest 2.3.0 

* Move hidden kernel related functions to package `statip`. 

* `discrete()` is now deprecated and will be removed in a future version of the 
package. Use `mfv()` instead. 

* `mfv1()` is a new function that always returns a length 1 value (so that `mfv1(x)==mfv(x)[[1L]]`). 


# modeest 2.2

* Thank you to W. H. Beasley who pointed out a slight mistake in the 
calculation of Bickel's skewness in `mlv.integer()`. Now the skewness is set at 
`NA` in case of multiple modes. 

* Add meanshift mode estimator. 

* As documented under `?as.numeric`, the function `as.numeric.mlv()` was not 
correct, and is now replaced by `as.double.mlv()`. 


# modeest 2.1

* Thank you to C. Lepoittevin and K. Fijorek who pointed out a misuse of 
`ifelse` in the function `hsm()`. This has been corrected, so now `hsm` works 
correctly. 

* Add functions `fiskMode()`, `gompertzMode()`, `koenkerMode()`, `kumarMode()`, 
`laplaceMode()`, `paralogisticMode()`, `paretoMode()`, `rayleighMode()`. 

* Remove function `symstbMode()`. 


# modeest 1.09

* Add the Asselin de Beauville mode estimator. 

* Add function `as.numeric.mlv()`. 

* Methods for the Chernoff distribution are provisionally suppressed, because 
of their lack of efficiency. 

* Remove the DIP statistic. 

* In function `tsybakov()`, the argument `djeddour` is renamed `dmp`. 

* In functions `parzen()` and `mlv.density()`, the argument `biau` is renamed `abc`. 
