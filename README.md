## SelectMixedModels

The package **SelectMixedModels** provides useful functions for model selection in linear mixed effects model (LMM).

This package follows idea and algorithm in the package [**MixedModels**](https://github.com/dmbates/MixedModels.jl) by Prof. D. Bates.

Here is a list of functions avaliable in  **SelectMixedModels**.

* trace of Hat matrix in LMM
* Conditional AIC

### An example 

An example of LMM  is a model fit to the `Dyestuff` data ([lme4 package for R](https://github.com/lme4/lme4))

In `R`, `cAIC4` package provides conditional AIC.

```R
> library(lme4)
> library(Matrix)
> library(cAIC4) 
> fm1 <- lmer(Yield ~ 1|Batch, Dyestuff ,REML=F)
> system.time(print(cAIC(fm1)))
$loglikelihood
[1] -157.4164

$df
[1] 6.282338

$reducedModel
NULL

$new
[1] FALSE

$caic
[1] 327.3974

   user  system elapsed 
  0.029   0.005   0.032 
```

In `julia`, `conAIC` type for conditional AIC is defined as

```Julia
type conAIC <: AIC
    condlike :: float64         #conditional likelihood
    tracehat :: Float64         # trace of Hat matrix
    method   :: ASCIIString     # method
    corterm  :: Vector          # correction term
    value    :: float64         # conditional AIC
end
```

A LMM model should be fitted with function `lmmg`, which forcely fit the model with `PLGGeneral` type in the packge `MixedModels`. Then use the function `conAIC` to compute conditional AIC.


```julia
julia> using MixedModels, SelectMixedModels, RDatasets

julia> ds = dataset("lme4", "Dyestuff");

julia> fm1 = fit(lmmg(Yield ~ 1|Batch, ds))
Linear mixed model fit by maximum likelihood
Formula: Yield ~ 1 | Batch

 logLik: -163.663530, deviance: 327.327060

 Variance components:
                Variance    Std.Dev.
 Batch        1388.333335   37.260345
 Residual     2451.250000   49.510100
 Number of obs: 30; levels of grouping factors: 6

  Fixed-effects parameters:
             Estimate Std.Error z value
(Intercept)    1527.5   17.6946  86.326


julia> @time conAIC(fm1)
elapsed time: 0.00011262 seconds (6680 bytes allocated)
conAIC(-157.41636044100218,4.6951603612889725,"VB",[1.0728,0.0766284],327.2054787453064)

```

### An Example with large data

R function `calAIC.lmer` for a simple conditional AIC is also avaliable in this package. 

My laptop is Ubuntu 14.04 64bit, 4GB memory, Intel® Core™ i5-3317U CPU @ 1.70GHz × 4.


```R
> fm4 <- lmer(score ~  (1|lea), Chem97)
> fm4
Linear mixed model fit by REML ['lmerMod']
Formula: score ~ (1 | lea)
   Data: Chem97
REML criterion at convergence: 161967
Random effects:
 Groups   Name        Std.Dev.
 lea      (Intercept) 0.7132  
 Residual             3.2768  
Number of obs: 31022, groups:  lea, 131
Fixed Effects:
(Intercept)  
      5.684  
> system.time(print(cAIC(fm4)))
Error in print(cAIC(fm4)) : 
  error in evaluating the argument 'x' in selecting a method for function 'print': Error: cannot allocate vector of size 7.2 Gb
Timing stopped at: 0.228 0 0.226 
> system.time(print(calAIC.lmer(fm4)))
[1] -80780.6168    113.4699 161788.1733
   user  system elapsed 
  0.041   0.000   0.041 
```

Also In `Julia`!

```Julia
julia> chem = dataset("mlmRev", "Chem97");

julia> fm4 = fit(lmmg(Score ~  (1|Lea), chem))
Linear mixed model fit by maximum likelihood
Formula: Score ~ 1 | Lea

 logLik: -80981.739482, deviance: 161963.478964

 Variance components:
                Variance    Std.Dev.
 Lea            0.503081    0.709282
 Residual      10.737177    3.276763
 Number of obs: 31022; levels of grouping factors: 131

  Fixed-effects parameters:
             Estimate Std.Error z value
(Intercept)   5.68449 0.0669684 84.8833


julia> @time conAIC(fm4)
elapsed time: 0.05279801 seconds (3231464 bytes allocated)
conAIC(-80780.87350003979,112.30490263302406,"VB",[1.00006,6.44787e-5],161788.3715450998)
```