# Firm-Credit-Risk
Credit risk measured by intensity model with frailty

> ## Intensity Model without Frailty
> Function `intensity.model()` performs the default intensity model based on *Duffie (2007)*.
>> Note that the model do not introduce the frailty variable.  
> ## Intensity Model with Industrial Frailty
> Function `industrial.frailty()` performs the default intensity model with a industrial frailty based on *Chava et al. (2009)*.  
> ## Fisher Dispersion Test (FDT)
> Function `FDT()` performs Fisher dispersion test used to check default dependency under the model assumption. The concept of FDT is kind of like chi-square test.  
>> You could also set the argument `fdt=TRUE` in `intensity.model()` and `industrial.frailty()` to perform FDT.
