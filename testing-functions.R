rm(list = ls())
import::here(fit.gaussian.process,
             generating.models,
             markov.decision.process,
             .from = '~/Dropbox/Shiny_Apps/MOGPDP/helpers.R')
test.time.series.length <- 100
test.data <- matrix(NA,nrow = 101)
test.data[1] <- 2
for (t in 1:(test.time.series.length)) {
  test.data[t + 1,] = generating.models(test.data[t,], r_ricker = 2.5, k_ricker = 0.5, 0) * 
    exp(rnorm(n = 1, mean = 0, sd = 0.0001 ^ 0.5))
}

plot(test.data, type = 'l')
test.output = fit.gaussian.process(test.data)
test.mdp.output <- markov.decision.process(test.output)
data.frame(hist = test.output$hist, 
           obs = test.output$obs)
data.frame(pred.hist = test.output$state.grid,
           pred.mean = exp(test.output$MM[[1]][,1]),
           pred.low.q = exp(test.output$MM[[1]][,1] - 2 * sqrt(diag(test.output$sig) + test.output$var[[1]])),
           pred.high.q = exp(test.output$MM[[1]][,1] + 2 * sqrt(diag(test.output$sig) + test.output$var[[1]])))
lines(test.output$state.grid,exp(test.output$MM[[1]][,1]),col = 'red')
dim(test.output$MM[[1]])

         