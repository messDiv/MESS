library(readxl)
library(reshape2)

setwd('~/Dropbox/Research/continuousity/MESS')

# read and clean occurrence data ----

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

snailOcc <- do.call(rbind, snailOcc)

# get the two most species rich islands
islandDiv <- tapply(snailOcc$species, snailOcc$island, function(x) length(unique(x)))
snailOccMaxRich <- snailOcc[snailOcc$island %in% names(sort(islandDiv, TRUE))[1:2], ]

# compress data down to just number of occs per species per island
x <- aggregate(list(num_occ = snailOccMaxRich[, 'occ']), 
               snailOccMaxRich[, c('island', 'species')], 
               sum)


# read and clean mt data ----

mt <- read_excel('empirical_data/galapagos_snails/MtDNAspecimenList.xlsx')
mt <- mt[mt$Island %in% snailOccMaxRich$island, ]
mt <- as.data.frame(mt)
names(mt) <- c('phyloID', 'island', 'pop', 'species')
mt <- mt[, c('island', 'species', 'phyloID')]

# record whether we have spp by island matches
mt$matched2occ <- paste(mt$island, mt$species) %in% paste(x$island, x$species)
x$matched2mt <- paste(x$island, x$species) %in% paste(mt$island, mt$species)

# get number of mt seq per spp per island
numMt <- sapply(paste(x$island, x$species), function(i) {
    sum(paste(mt$island, mt$species) == i)
})
x$num_mtSeq <- numMt
