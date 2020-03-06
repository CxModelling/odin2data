# odin2data
A model fitting toolkit for simulation models with the odin package


## Install package
In R, you can install from github
```r
devtools::install_github("TimeWz667/odin2data", upgrade = FALSE)
```

Load the package before use it
```r
library(odin2data)
```

## Construct a model
```r
### Load the model file
f <- system.file("models/SIR.txt", package = "odin2data")
test_m <- odin::odin(f, verbose=F)


### Generate a prior and set up it as a list
r_prior <- function() {
  list(beta = runif(1, 0.1, 5), gamma = runif(1, 0.1, 0.3))
}

d_prior <- function(pars) {
  dunif(pars$beta, 0.1, 5, log = T) + dunif(pars$gamma, 0.1, 0.3, log = T)
}

times = seq(0, 10, 0.2)
y0 <- c(995, 5, 0)

### Compile all elements as a simulation model
sim <- odin2data::compile_model(d_prior, r_prior, y0, ts_sim = times, m_sim = test_m)

```

## Prepare data
```r
test_data <- data.frame(
  t = 1:5,
  incidence = c(20, 49, 109, 184, 206) / 1000
)

### Compile the model with data
lf <- odin2data::compile_model_likefree(test_data, sim)
```

## Fit model to data

odin2data provides three types of fitting algorithm currently

- **ABC** 
- **ABC-SMC**
- **ABC-PMC**


```r
post <- odin2data::fit(lf, 1000, method = "abc", target_acc = 0.01)

post <- odin2data::fit(lf, 1000, method = "abcsmc", max_round = 20, alpha = 0.80)

post <- odin2data::fit(lf, 1000, method = "abcpmc", max_round = 20)

summary(post)
```


## Academic contact

Chu-Chang Ku,

Health Economic and Decision Science, University of Sheffield, UK

Email: c.ku@sheffield.ac.uk


## License

MIT (c) 2020 Chu-Chang Ku

