## ---- include = FALSE--------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)


## ----setup, include = FALSE--------------------------------------------
set.seed(1341)
library(causalOT)


## ----------------------------------------------------------------------
n <- 5
d <- 3
x <- matrix(stats::rnorm(n*d), nrow = n)
w <- stats::runif(n)
w <- w/sum(w)

m <- Measure(x = x, weights = w)


## ----------------------------------------------------------------------
 m$x
 m$weights


## ----------------------------------------------------------------------
m <- Measure(x = x, weights = w, 
             probability.measure = TRUE,
             adapt = "none")


## ----------------------------------------------------------------------
target.data <- matrix(rnorm(n*d),n,d)
target.values <- colMeans(target.data)

m <- Measure(x = x, weights = w, 
             probability.measure = TRUE,
             adapt = "weights",
             target.values = target.values)


## ---- eval = FALSE-----------------------------------------------------
## m$balance_functions # to view the balance functions
## m$balance_target # to view the target values


## ----------------------------------------------------------------------
all.equal(as.numeric(m$balance_target), target.values)
all.equal(as.matrix(m$balance_functions), x)

sds <- apply(x,2,sd)
all.equal(as.numeric(m$balance_target), target.values/sds)
all.equal(as.matrix(m$balance_functions), sweep(x,2,sds,"/"))


## ----------------------------------------------------------------------
m


## ----measure-----------------------------------------------------------
m_target <- Measure(x = matrix(rnorm(n*2*d), n*2,d))
m_source <- Measure(x = x, weights = w, adapt = "weights")


## ----problem-----------------------------------------------------------
otp <- OTProblem(m_source, m_target)


## ----options, warning = FALSE------------------------------------------
otp$setup_arguments(
  lambda = NULL, # penalty values of the optimal transport (OT) distances to try
  delta = NULL, # constraint values to try for balancing functinos
  grid.length = 7L, # number of values of lambda and delta to try
  # if none are provided
  cost.function = NULL, # the ground cost to use between covariates
  # default is the Euclidean distance
  p = 2, # power to raise the cost by
  cost.online = "auto", #Should cost be calculated "online" or "tensorized" (stored in memory). "auto" will try to decide for you
  debias = TRUE, # use Sinkhorn divergences (debias = TRUE), i.e. debiased Sinkhorn distances, 
  # or use the Sinkhorn distances (debias = FALSE)
  diameter = NULL, # the diameter of the covariate space if known
  ot_niter = 1000L, # the number of iterations to run when solving OT distances
  ot_tol = 0.001 # the tolerance for convergance of OT distances
)


## ----solve, cache = TRUE-----------------------------------------------
otp$solve(
  niter = 1000L, # maximum number of iterations
  tol = 1e-5, # tolerance for convergence
  optimizer = "torch", # which optimizer to use "torch" or "frank-wolfe"
  torch_optim = torch::optim_lbfgs, # torch optimizer to use if required
  torch_scheduler = torch::lr_reduce_on_plateau, # torch scheduler to use if required
  torch_args = list(line_search_fn = "strong_wolfe"), # args passed to the torch functions,
  osqp_args = NULL, #arguments passed to the osqp solver used for "frank-wolfe" and balance functions
  quick.balance.function = TRUE # if balance functions are also present, should an approximate value of the hyperparameter "delta" be found first
)


## ---- cache = TRUE, echo = FALSE---------------------------------------
cbind(adapted = as.numeric(m_source$weights), 
      original = w)


## ----hyper, cache = TRUE-----------------------------------------------
otp$choose_hyperparameters(
  n_boot_lambda = 100L, #Number of bootstrap iterations to choose lambda
  n_boot_delta = 1000L, #Number of bootstrap iterations to choose delta
  lambda_bootstrap = Inf # penalty parameter to use for OT distances
)


## ---- cache = TRUE-----------------------------------------------------
otp$selected_lambda


## ---- cache = TRUE-----------------------------------------------------
as.numeric(m_source$weights)


## ----------------------------------------------------------------------
otp$loss


## ---- eval = FALSE-----------------------------------------------------
## m_target <- Measure(x = matrix(rnorm(n*2*d), n*2,d))
## m_source <- Measure(x = x, weights = w, adapt = "weights")


## ---- eval = FALSE-----------------------------------------------------
## otp <- OTProblem(m_source, m_target)


## ---- eval = FALSE-----------------------------------------------------
## otp$setup_arguments(
##   lambda = NULL, # penalty values of the optimal transport (OT) distances to try
##   delta = NULL, # constraint values to try for balancing functinos
##   grid.length = 7L, # number of values of lambda and delta to try
##   # if none are provided
##   cost.function = NULL, # the ground cost to use between covariates
##   # default is the Euclidean distance
##   p = 2, # power to raise the cost by
##   cost.online = "auto", #Should cost be calculated "online" or "tensorized" (stored in memory). "auto" will try to decide for you
##   debias = TRUE, # use Sinkhorn divergences (debias = TRUE), i.e. debiased Sinkhorn distances,
##   # or use the Sinkhorn distances (debias = FALSE)
##   diameter = NULL, # the diameter of the covariate space if known
##   ot_niter = 1000L, # the number of iterations to run when solving OT distances
##   ot_tol = 0.001 # the tolerance for convergance of OT distances
## )


## ---- eval = FALSE-----------------------------------------------------
## otp$solve(
##   niter = 1000L, # maximum number of iterations
##   tol = 1e-5, # tolerance for convergence
##   optimizer = "torch", # which optimizer to use "torch" or "frank-wolfe"
##   torch_optim = torch::optim_lbfgs, # torch optimizer to use if required
##   torch_scheduler = torch::lr_reduce_on_plateau, # torch scheduler to use if required
##   torch_args = list(line_search_fn = "strong_wolfe"), # args passed to the torch functions,
##   osqp_args = NULL, #arguments passed to the osqp solver used for "frank-wolfe" and balance functions
##   quick.balance.function = TRUE # if balance functions are also present, should an approximate value of the hyperparameter "delta" be found first
## )


## ---- eval = FALSE-----------------------------------------------------
## otp$choose_hyperparameters(
##   n_boot_lambda = 100L, #Number of bootstrap iterations to choose lambda
##   n_boot_delta = 1000L, #Number of bootstrap iterations to choose delta
##   lambda_bootstrap = Inf # penalty parameter to use for OT distances
## )


## ----------------------------------------------------------------------
as.numeric(m_source$weights)


## ----------------------------------------------------------------------
nrow <- 100
ncol <- 2
a <- Measure(x = matrix(rnorm(nrow*ncol,mean=c(0.1,0.1)) + 0.1,nrow,ncol,byrow = TRUE), adapt = "weights")
b <- Measure(x = matrix(rnorm(nrow*ncol,mean=c(-0.1,-0.1),sd=0.25),nrow,ncol,byrow = TRUE), adapt = "weights")
c <- Measure(x = matrix(rnorm(nrow*ncol,mean=c(0.1,-0.1)),nrow,ncol,byrow = TRUE), adapt = "weights")
d <- Measure(x = matrix(rnorm(nrow*ncol,mean= c(-0.1,0.1),sd=0.25),nrow,ncol,byrow = TRUE), adapt = "weights")

overall <- Measure(x = torch::torch_vstack(lapply(list(a,b,c,d), function(meas) meas$x)), 
                   adapt = "none")

overall_ot <- OTProblem(a,overall) + OTProblem(b, overall) +
  OTProblem(c, overall) + OTProblem(d, overall)

## ---- eval = FALSE-----------------------------------------------------
## overall_ot$setup_arguments()
## overall_ot$solve()
## overall_ot$choose_hyperparameters()
## 


## ----------------------------------------------------------------------
overall_ot


## ---- eval = FALSE-----------------------------------------------------
## source_measures <- list(a,b,c,d)
## meas <- x_temp <- NULL
## z_temp <- c(rep(1, nrow*4), rep(0,nrow))
## wt  <- list()
## for(i in seq_along(source_measures)) {
##   meas <- source_measures[[i]]
##   x_temp <- as.matrix(torch::torch_vstack(list(overall$x,meas$x)))
##   wt[[i]] <- calc_weight(x = x_temp,
##                         z = z_temp,
##                         estimand = "ATT",
##                         method = "COT")
## }
## 


## ----------------------------------------------------------------------
target.values <- 
  as.numeric(a$x$mean(1) + b$x$mean(1) + 
             c$x$mean(1) + d$x$mean(1))/4
a_t <- Measure(x = a$x, adapt = "weights",
               target.values = target.values)
b_t <- Measure(x = a$x, adapt = "weights",
               target.values = target.values)
c_t <- Measure(x = a$x, adapt = "weights",
               target.values = target.values)
d_t <- Measure(x = a$x, adapt = "weights",
               target.values = target.values)

all.target.measures <- list(a_t, b_t, c_t, d_t)


## ----warning=FALSE, cache = TRUE---------------------------------------
ot_targ <- NULL
for(meas in all.target.measures) {
  ot_targ <- OTProblem(meas, meas)
  ot_targ$setup_arguments(lambda = 100)
  ot_targ$solve(torch_optim = torch::optim_lbfgs,
                torch_args = list(line_search_fn = "strong_wolfe"))
}


## ---- cache = TRUE-----------------------------------------------------
final.bal <- as.numeric(a_t$x$mT()$matmul(a_t$weights$detach()))
original  <- as.numeric(a_t$x$mean(1))
rbind(original,
      `final balance` = final.bal, 
      `target values` = target.values)


## ----------------------------------------------------------------------
pseudo <- Measure(x = matrix(rnorm(nrow*4*ncol), nrow*4, ncol),
                  adapt = "x")


## ----------------------------------------------------------------------
pseudo_a <- pseudo$detach()
pseudo_b <- pseudo$detach()
pseudo_c <- pseudo$detach()
pseudo_d <- pseudo$detach()

pseudo_a$requires_grad <- pseudo_b$requires_grad <-
pseudo_c$requires_grad <- pseudo_d$requires_grad <- "x"
ota <- OTProblem(a$detach(), # don't update a
                 pseudo_a)

otb <- OTProblem(b$detach(), # don't update b
                 pseudo_b)

otc <- OTProblem(c$detach(), # don't update c
                 pseudo_c)

otd <- OTProblem(d$detach(), # don't update c
                 pseudo_d)



## ----------------------------------------------------------------------
ota$setup_arguments(lambda = .1)
otb$setup_arguments(lambda = .1)
otc$setup_arguments(lambda = .1)
otd$setup_arguments(lambda = .1)


## ----------------------------------------------------------------------
 opt <- torch::optim_rmsprop(pseudo$x)
 sched <- torch::lr_multiplicative(opt, lr_lambda = function(epoch) {0.99})


## ---- eval = FALSE-----------------------------------------------------
## 
## #optimization loop
## for (i in 1:100) {
##   # zero grad of main optimizer
##     opt$zero_grad()
##   # get gradients at each site
##     ota$loss$backward()
##     otb$loss$backward()
##     otc$loss$backward()
##     otd$loss$backward()
##   # pass grads back to main site
##     pseudo$grad <- pseudo_a$grad + pseudo_b$grad +
##       pseudo_c$grad + pseudo_d$grad
## 
##   # update pseudo data at main site
##     opt$step()
## 
##   # zero site gradients
##     torch::with_no_grad({
##       pseudo_a$grad$copy_(0.0)
##       pseudo_b$grad$copy_(0.0)
##       pseudo_c$grad$copy_(0.0)
##       pseudo_d$grad$copy_(0.0)
##       })
##   # update scheduler
##     sched$step()
## }
## 


## ---- eval = FALSE-----------------------------------------------------
## pseudo_a$x <- pseudo_b$x <- pseudo_c$x <- pseudo_d$x <-
##    pseudo$x
## 
## ota_w <- OTProblem(a, pseudo_a$detach())
## otb_w <- OTProblem(b, pseudo_b$detach())
## otc_w <- OTProblem(c, pseudo_c$detach())
## otd_w <- OTProblem(d, pseudo_d$detach())
## 
## 
## ota_w$setup_arguments()
## ota_w$solve(torch_args = list(line_search_fn = "strong_wolfe"))
## ota_w$choose_hyperparameters()
## 
## otb_w$setup_arguments()
## otb_w$solve(torch_args = list(line_search_fn = "strong_wolfe"))
## otb_w$choose_hyperparameters()
## 
## 
## otc_w$setup_arguments()
## otc_w$solve(torch_args = list(line_search_fn = "strong_wolfe"))
## otc_w$choose_hyperparameters()
## 
## otd_w$setup_arguments()
## otd_w$solve(torch_args = list(line_search_fn = "strong_wolfe"))
## otd_w$choose_hyperparameters()
## 


## ----------------------------------------------------------------------
a$weights <- a$init_weights
b$weights <- b$init_weights
c$weights <- c$init_weights
d$weights <- d$init_weights


## ----------------------------------------------------------------------
pseudo <- Measure(x = matrix(rnorm(nrow*4*ncol), nrow*4, ncol),
                  adapt = "x")


## ----------------------------------------------------------------------
pseudo_a <- pseudo$detach()
pseudo_b <- pseudo$detach()
pseudo_c <- pseudo$detach()
pseudo_d <- pseudo$detach()

pseudo_a$requires_grad <- pseudo_b$requires_grad <- 
pseudo_c$requires_grad <- pseudo_d$requires_grad <- "x"
ota <- OTProblem(a$detach(), # don't update a
                 pseudo_a)

otb <- OTProblem(b$detach(), # don't update b
                 pseudo_b)

otc <- OTProblem(c$detach(), # don't update c
                 pseudo_c)

otd <- OTProblem(d$detach(), # don't update c
                 pseudo_d)



## ----------------------------------------------------------------------
ota$setup_arguments(lambda = .1)
otb$setup_arguments(lambda = .1)
otc$setup_arguments(lambda = .1)
otd$setup_arguments(lambda = .1)


## ---- eval = FALSE-----------------------------------------------------
## # run separately at each site
## ota$solve(torch_optim = torch::optim_rmsprop)
## otb$solve(torch_optim = torch::optim_rmsprop)
## otc$solve(torch_optim = torch::optim_rmsprop)
## otd$solve(torch_optim = torch::optim_rmsprop)
## 
## 


## ----echo = FALSE, fig.show='hold', fig.path = "oop-", cache = TRUE----
otb$solve(torch_optim = torch::optim_rmsprop)
plot(as.matrix(b$x),xlim = c(-3,3), ylim = c(-3,3),
     xlab = expression(X[1]),
     ylab = expression(X[2]),
     main = "Before")
points(as.matrix(pseudo$x), col = rgb(red = 0, green = 0, blue = 1, alpha = 0.25))

plot(as.matrix(b$x),xlim = c(-3,3), ylim = c(-3,3),
     xlab = expression(X[1]),
     ylab = expression(X[2]),
     main = "After")
points(as.matrix(pseudo_b$x), col = rgb(red = 0, green = 0, blue = 1, alpha = 0.25))


## ---- eval = FALSE-----------------------------------------------------
## # send back to the main site and create overall problem
## ot_overall <-
##   OTProblem(pseudo_a$detach(),
##             pseudo) +
##   OTProblem(pseudo_b$detach(),
##             pseudo) +
##   OTProblem(pseudo_c$detach(),
##             pseudo) +
##   OTProblem(pseudo_d$detach(),
##             pseudo)
## 
## ot_overall$setup_arguments(lambda = 0.1)
## 
## 
## ot_overall$solve(torch_optim = torch::optim_rmsprop)
## 


## ---- eval = FALSE-----------------------------------------------------
## # pass pseudo to each site then setup the problems again
## ota2 <- OTProblem(a,
##                  pseudo$detach())
## 
## otb2 <- OTProblem(b, # don't update b
##                  pseudo$detach())
## 
## otc2 <- OTProblem(c,
##                  pseudo$detach())
## 
## otd2 <- OTProblem(d,
##                  pseudo$detach())
## 
## all.problems <- list(ota2,
##                      otb2,
##                      otc2,
##                      otd2)
## 
## # then we optimize the weights at each site separately.
## for (prob in all.problems) {
##   prob$setup_arguments()
##   prob$solve(
##     torch_optim = torch::optim_lbfgs,
##     torch_args = list(line_search_fn = "strong_wolfe")
##     )
##   prob$choose_hyperparameters()
## }
## 


## ----------------------------------------------------------------------
x_1 <- matrix(rnorm(128*2),128) + 
  matrix(c(-0.1,-0.1), 128, 2,byrow = TRUE)
x_2 <- matrix(rnorm(256*2), 256) + 
  matrix(c(0.1,0.1), 256, 2,byrow = TRUE)

target.data <- matrix(rnorm(512*2), 512, 2) * 0.5 + 
  matrix(c(0.1,-0.1), 512, 2, byrow = TRUE)
constructor.formula <- formula("~ 0 + . + I(V1^2) + I(V2^2)")
target.values <- colMeans(model.matrix(constructor.formula,
                                       as.data.frame(target.data)))

m_1 <- Measure(x = x_1, adapt = "weights",
               balance.functions = model.matrix(constructor.formula, 
                                                as.data.frame(x_1)),
               target.values = target.values)
m_2 <- Measure(x = x_2, adapt = "weights",
               balance.functions = model.matrix(constructor.formula, 
                                                as.data.frame(x_2)),
               target.values = target.values)

ot_binary <- OTProblem(m_1, m_2)


## ---- cache = TRUE-----------------------------------------------------
ot_binary$setup_arguments() 

ot_binary$solve(torch_optim = torch::optim_lbfgs,
                torch_args = list(line_search_fn = "strong_wolfe"))

ot_binary$choose_hyperparameters()



## ----------------------------------------------------------------------
info <- ot_binary$info()
names(info)


## ----------------------------------------------------------------------
info$balance.function.differences


## ---- cache = TRUE-----------------------------------------------------
c(initial = ot_distance(m_1$x, m_2$x, 
            a = m_1$init_weights, b = m_2$init_weights, penalty = 1),
final = ot_distance(m_1$x, m_2$x, 
            a = m_1$weights, b = m_2$weights, penalty = 1))


