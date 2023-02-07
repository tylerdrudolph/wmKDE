#########################################################################################################################################################
## Conduct model validation of an RSF model based on training and validation datasets
## Implements procedure published by Howlin, Erickson & Nielson (2003): "A validation technique for assessing predictive abilities of resource selection functions" (sample validation procedure)
##
## Written by Tyler D. Rudolph, M.Sc.
## tylerdrudolph@gmail.com
##
## Date: Aug. 27, 2021
##
## Arguments:
# trainset = GPS relocations serving to estimate the reference (e.g. mean weighted population) kernel (training relocations)
# validset = GPS (e.g. new) relocations not serving to estimate the reference kernel (validation relocations)
# kmat = R object of class "RasterLayer" corresponding to the reference kernel, rescaled between 0 and 1 (max probability of occurrence)
#
# plot = logical argument whether or not to plot linear model of observed vs. predicted
# id = unique identifier associated with validation iteration
##
## Output = list containing three elements:
# 1) validation.tab = summary of validation output
# 2) deviance.tab = distribution of observed counts vs. expected (scaled RSF) values
# 3) linear.model = lm object resulting from linear regression of observed vs. expected values
##
## COMMENTS:

# validateKDE = function(trainset, validset, kernmod, plot=FALSE, id=NA) {
#
#   # Estimate model
#   #sampmod = coxph(formula, trainset, robust=TRUE)
#
#   # Ensure kernel densities have been rescaled between 0 and 1 (** seeking correct approach)
#   #sampleRandom(kernmod, size=100) * 100
#
#   # Derive predicted values
#   #trainset$RSFexp = exp(predict(sampmod, trainset, type="lp", reference="sample"))
#   #validset$RSFexp = exp(predict(sampmod, newdata=validset, type="lp", reference="sample"))
#   trainset$kpred <- extract(kernmod, trainset)
#   validset$kpred <- extract(kernmod, validset)
#
#   # Scale predicted relative probabilities of selection (see Howlin, Erickson, & Neilson 2003 in Resource Selection Methods and Applications)
#   #trainset = trainset[trainset$response==0,] # remove used points
#   #trainset$RSFexp.scaled = (trainset$RSFexp/sum(trainset$RSFexp))*nrow(validset)
#   #validset$RSFexp.scaled = (validset$RSFexp/sum(trainset$RSFexp))*nrow(validset)
#   trainset$kpred.scaled <- (trainset$kpred / sum(trainset$kpred)) * nrow(validset)
#   validset$kpred.scaled <- (validset$kpred / sum(trainset$kpred)) * nrow(validset)
#
#   # Create bins and divide scaled predicted relative probabilities of selection (kpred.scaled) for training set
#   #nbins <- ceiling(nrow(trainset)/nrow(validset))
#   #cutblocks <- cut(trainset$RSFexp.scaled, breaks=nbins, ordered.result=TRUE)
#   nbins <- ceiling(nrow(trainset) / nrow(validset))
#   cutblocks <- cut(trainset$kpred.scaled, breaks=nbins, ordered.result=T)
#
#   # Sum scaled predicted values by bin (equivalent to expected values)
#   #devtab = tapply(trainset$RSFexp.scaled[trainset$response==0], list(cutblocks[trainset$response==0]), function(x) sum(x, na.rm=T))
#   #devtab = data.frame(bins = names(devtab), expected = devtab)
#   #devtab[is.na(devtab$expected),"expected"] = 0
#   #row.names(devtab) = 1:nrow(devtab)
#   devtab = tapply(trainset$kpred.scaled, list(cutblocks), function(x) sum(x, na.rm=T))
#   devtab = data.frame(bins = names(devtab), expected = devtab)
#   devtab[is.na(devtab$expected),"expected"] = 0
#   row.names(devtab) = 1:nrow(devtab)
#
#   # Extract bin values
#   binvals = do.call(rbind, lapply(strsplit(as.character(devtab$bins), ","), function(x) {
#     x[1] = strsplit(x[1],"(", fixed=TRUE)[[1]][2]
#     x[2] = strsplit(x[2],"]", fixed=TRUE)[[1]][1]
#     return(as.numeric(x))
#   }))
#
#   # Sum of sums should be equal to number of observations in validation set
#   if(round(sum(devtab$expected))!=nrow(validset)) stop("incorrect distribution of expected values")
#
#   # Extend bin values to ensure inclusion of maximum RSF values observed in validation dataset
#   #if(max(binvals) < max(validset$RSFexp.scaled)) binvals[which.max(binvals)] = ceiling(max(validset$RSFexp.scaled))
#   if(max(binvals) < max(validset$kpred.scaled)) binvals[which.max(binvals)] = max(validset$kpred.scaled)
#
#   # Divide validation set data into bins
#   #devtab$observed = apply(binvals, 1, function(bins) sum(validset$RSFexp.scaled > bins[1] & validset$RSFexp.scaled <= bins[2]))
#   devtab$observed = apply(binvals, 1, function(bins) sum(validset$kpred.scaled > bins[1] & validset$kpred.scaled <= bins[2]))
#
#   # Ensure conformity
#   if(round(sum(devtab$expected))!=round(sum(devtab$observed))) stop("error in count distribution")
#
#   # Fit linear model
#   diagmod = lm(observed~expected, devtab)
#   civals = confint(diagmod)
#
#   if(plot) {
#     plot(observed~expected, devtab)
#     lines(fitted(diagmod)~devtab$expected, pch=19, lty="dotted", main=id)
#   }
#
#   # Stock output
#   output = data.frame(id, nbins, lower.ci=civals[2,1], beta=diagmod$coef[2], upper.ci=civals[2,2], pval.lm=summary(diagmod)$coefficients[2,4], spearman.cor=cor(devtab$observed, devtab$expected, method="spearman"))#, pval.ks=ks.test(devtab$expected,devtab$observed)$p.value)
#
#   return(list(validation.tab = output, deviance.tab=devtab, linear.model=diagmod))
#
# }
#
#
# ## Implementation (working example in development)
# setwd('/mnt/DATA/R/PACKAGES/mwKDE')
# #data(cartab)
#
# ## Load training set (< 2020 relocation data)
# load('/mnt/DATA/MFFP/Analyses_GSzor/mwKDE/data/cartab.RData')
# source('/mnt/DATA/MFFP/Analyses_GSzor/mwKDE/R/01_fonctions_de_base.R')
#
# ## Load functions
# sapply(dir('R'), function(script) source(str_c('R/', script)))
#
# ## Rarify and convert to sf object
# trainset <- subsetDB(rarifyGPS(cartab, date.heure=cartab$date.heure, id=cartab$IDAnimal), pop='Assinica')
#
# ## Estimate population kernel
# kde.list <- mwKDE(xy=trainset[,c('x','y')],
#       spatres=1000,
#       spwt=F, #wts=trainset$poids,
#       avg=TRUE,
#       id=trainset$IDAnimal,
#       bw.global=TRUE,
#       zscale=T,
#       export=F,
#       obj.ret = T)
#
# ## Rescale kernel densities to between 0 (lowest probability) and 1 (highest probability)
# kernmod <- kde.list$mwKern
#   kernmod$fhat <- (100 - fhat2confin(kernmod$fhat)) / 100
#   kernmod <- UD2rast(kernmod)
#
# ## Import validation set (>= 2020 collar telemetry relocations)
# validset <- filter(readRDS('/mnt/DATA/MFFP/Analyses_GSzor/mwKDE/data/vdata.rds'), Pop=='Assinica')
#
# ## Validate fit of new data to kernel model
# vout  <- lapply(unique(validset$IDAnimal), function(i) validateKDE(trainset, filter(validset, IDAnimal==i), kernmod, plot=T, id=i))
#
# do.call(rbind, lapply(vout, function(x) x$validation.tab))

