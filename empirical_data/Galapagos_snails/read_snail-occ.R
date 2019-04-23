library(readxl)
library(reshape2)

setwd('~/Dropbox/Research/continuousity/MESS')

sheets <- 1:17
snailOcc <- lapply(sheets, function(i) {
    thisDat <- read_excel('empirical_data/Galapagos_snails/Incidence_GalapagosSnails.xlsx', 
                          sheet = i)
    names(thisDat)[names(thisDat) == 'X__1'] <- 'name2'
    names(thisDat)[names(thisDat) == 'elevation'] <- 'ele'
    
    i <- which(names(thisDat) == 'ele')
    
    info <- thisDat[, 1:i]
    if(!('name2' %in% names(info))) info$name2 <- NA
    info <- info[, c('name', 'name2', 'island', 'lon', 'lat', 'ele')]
    
    thisDat <- cbind(info, thisDat[, -(1:i)])
    thisDat <- melt(thisDat, id.vars = names(info))
    names(thisDat)[7:8] <- c('species', 'occ')
    
    occ <- rep(0, nrow(thisDat))
    occ[!is.na(thisDat$occ)] <- 1
    thisDat$occ <- occ
    thisDat <- thisDat[thisDat$occ == 1, ]
    
    return(thisDat)
})

bla <- do.call(rbind, snailOcc)
