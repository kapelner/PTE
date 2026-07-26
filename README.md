Personalized Treatment Evaluator
===

Code for the R package PTE as well as code to duplicate the figures found in 

Kapelner, A, Bleich, J, Cohen, ZD, DeRubeis, RJ and Berk, R (2014) Inference for Treatment Regime Models in Personalized Medicine, arXiv

## Installing for maximum speed

By default the package compiles portably (no machine-specific CPU flags), which is required for CRAN and safe to redistribute. If you're installing from source on the machine you'll actually run it on and want the fastest possible build, opt into native-CPU optimization:

```
PTE_NATIVE_SPEED=1 R CMD INSTALL .
```

This compiles with `-march=native -mtune=native -O3 -flto`, tying the resulting binary to the exact CPU that compiled it -- don't use this for a build you plan to share or distribute.




