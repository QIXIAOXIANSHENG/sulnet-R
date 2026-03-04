library(devtools)
clean_dll()
document()
load_all()

y = sample(0:1, 10, replace = TRUE)
x = matrix(rnorm(10*5), ncol=5)
sul = uniFit(x,y,family = "binomial")
ul = uniInfo(x,y,"binomial", loo = TRUE)
