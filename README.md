# barcoo
## RNL and NL models for barcoo

### adultjuv
### Netlogo:
Sleepy IBM_v.6.2.nlogo = Sleepy IBM_v.6.1.1_two strategies_shadedens.nlogo (as of 19-7-17)

### R:
RNL_adultjuv.R = RNL_new trans model_with DEB_1.6.2.1_juvenile.R (as of 19-7-17)

### barcoo arguments in RNL_adultjuv.R
```{r}
strat<-as.numeric(args[2]) # movement strategy (0 = opt; 1 = sat)
resource<-as.numeric(args[3]) # resource density (0 = low; 1 = high)
juv <- as.numeric(args[4]) # life stage (0 = adult; 1 = juvenile). search 
```

