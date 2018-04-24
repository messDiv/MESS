## function to make the mainland meta community with a phylo, traits, and abundances
## required packages:
##    ape,
##    TreeSim,
##    sads
## arguments:
#' @param Jm the number of individuals in the meta community
#' @param S the number of species in the meta community
#' @param lambda the speciation rate
#' @param deathFrac the proportional extinction rate
#' @param sigma2 the rate of brownian motion

makeMeta <- function(Jm, S, lambda, deathFrac, sigma2) {
    ## the tree
    tre <- TreeSim::sim.bd.taxa(S, numbsim = 1, lambda = lambda, mu = lambda * deathFrac, 
                                complete = FALSE)[[1]]
    
    ## the traits
    trt <- ape::rTraitCont(tre, sigma = sqrt(sigma2))
    
    ## parameters for the log-series
    nBar <- Jm / S
    p <- 1 - 1/nBar
    b <- -log(p)
    
    ## the abundances
    abund <- sads::rls(length(trt), length(trt), 0.01)
    
    ## return it all in a list
    return(list(phylo = tre, traits = trt, abundance = abund))
}
