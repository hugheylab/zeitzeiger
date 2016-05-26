zeitzeigerPredictGivenDensity = function(xTest, xFitMean, xFitVar, beta, timeRange=seq(0, 1, 0.01)) {
	negLoglike = -zeitzeigerLikelihood(xTest, xFitMean, xFitVar, beta, timeRange, logArg=TRUE)
	mleFit = list()
	for (ii in 1:nrow(xTest)) {
		xTestNow = xTest[ii,, drop=FALSE]
		timeStart = timeRange[which.min(negLoglike[ii,])]
		names(timeStart) = 'time'
		negLoglikeFunc = function(time) -sum(fx(xTestNow, time, xFitMean, xFitVar, logArg=TRUE) * beta)
		bbmle::parnames(negLoglikeFunc) = names(timeStart)
		warnOrig = getOption('warn')
		options(warn=-1)
		mleFit[[ii]] = bbmle::mle2(negLoglikeFunc, start=timeStart, vecpar=TRUE, method='L-BFGS-B', lower=0, upper=1)
		options(warn=warnOrig)}
	timePred = sapply(mleFit, function(x) x@coef)
	return(list(timeDepLike=exp(-negLoglike), mleFit=mleFit, timePred=timePred))}


#' Predict corresponding time for test observations.
#'
#' \code{zeitzeigerPredict} predicts the value of the periodic variable
#' for test observations, given training data and SPCs. This is a wrapper
#' around \code{bbmle::mle2}.
#'
#' @param xTrain Matrix of measurements for training data, observations in rows
#' and features in columns.
#' @param timeTrain Vector of values of the periodic variable for training observations,
#' where 0 corresponds to the lowest possible value and 1 corresponds to the highest
#' possible value.
#' @param xTest Matrix of measurements for test data, observations in rows
#' and features in columns.
#' @param spcResult Result from \code{zeitzeigerSpc}.
#' @param fitMeanArgs List of arguments to pass to \code{bigspline} for fitting mean of each SPC.
#' @param constVar Logical indicating whether to assume constant variance as a function
#' of the periodic variable.
#' @param fitVarArgs List of arguments to pass to \code{bigspline} for fitting variance of each SPC.
#' Unused if \code{constVar==TRUE}.
#' @param nSpc Vector of the number of SPCs to use for prediction. If \code{NA} (default),
#' \code{nSpc} will become \code{1:K}, where \code{K} is the number of SPCs in \code{spcResult}.
#' Each value in \code{nSpc} will correspond to one prediction for each test observation.
#' A value of 2 means that the prediction will be based on the first 2 SPCs.
#' @param betaSv Logical indicating whether to use the singular values of the SPCs
#' as weights in the likelihood calculation.
#' @param timeRange Vector of values of the periodic variable at which to calculate likelihood.
#' The time with the highest likelihood is used as the initial value for the MLE optimizer.
#'
#' @return
#' \item{timeDepLike}{3-D array of likelihood, with dimensions for each test observation,
#' each element of \code{nSpc}, and each element of \code{timeRange}.}
#' \item{mleFit}{List (for each element in \code{nSpc}) of lists (for each test observation)
#' of \code{mle2} objects.}
#' \item{timePred}{Matrix of predicted times for test observations by values of \code{nSpc}.}
#'
#' @export
zeitzeigerPredict = function(xTrain, timeTrain, xTest, spcResult, fitMeanArgs=list(rparm=NA), constVar=TRUE,
									  fitVarArgs=list(rparm=NA), nSpc=NA, betaSv=FALSE, timeRange=seq(0, 1, 0.01)) {
	zTrain = xTrain %*% spcResult$v
	zTest = xTest %*% spcResult$v

	zFitResult = zeitzeigerFit(zTrain, timeTrain, fitMeanArgs)
	zFitMean = zFitResult$xFitMean
	zFitResid = zFitResult$xFitResid
	zFitVar = zeitzeigerFitVar(timeTrain, zFitResid, constVar, fitVarArgs)

	if (length(nSpc)==1 && is.na(nSpc)) {
		nSpc = 1:length(spcResult$d)
	} else {
		if (!all(nSpc %in% 1:length(spcResult$d)) || anyDuplicated(nSpc)>0) {
			stop('nSpc must be unique integers in 1 to the number of singular vectors.')}}

	timePred = matrix(NA, nrow=nrow(xTest), ncol=length(nSpc))
	timeDepLike = array(NA, dim=c(nrow(xTest), length(nSpc), length(timeRange)))
	mleFit = list() # list (for each nSpc) of lists (for each observation)
	for (ii in 1:length(nSpc)) {
		if (betaSv) {
			beta = spcResult$d[1:nSpc[ii]]
			# beta = spcResult$d[1:nSpc[ii]]^2 / sum(spcResult$d[1:nSpc[ii]]^2)
		} else {
			beta = rep_len(1, nSpc[ii])}
		predResult = zeitzeigerPredictGivenDensity(zTest[,1:nSpc[ii], drop=FALSE], zFitMean[1:nSpc[ii]], zFitVar[1:nSpc[ii]], beta)
		timeDepLike[,ii,] = predResult$timeDepLike
		mleFit[[ii]] = predResult$mleFit
		timePred[,ii] = predResult$timePred}
	return(list(timeDepLike=timeDepLike, mleFit=mleFit, timePred=timePred))}


#' Train and test a ZeitZeiger predictor.
#'
#' \code{zeitzeiger} sequentially calls \code{zeitzeigerFit}, \code{zeitzeigerSpc},
#' and \code{zeitzeigerPredict}.
#'
#' @param xTrain Matrix of measurements for training data, observations in rows
#' and features in columns.
#' @param timeTrain Vector of values of the periodic variable for training observations,
#' where 0 corresponds to the lowest possible value and 1 corresponds to the highest
#' possible value.
#' @param xTest Matrix of measurements for test data, observations in rows
#' and features in columns.
#' @param fitMeanArgs List of arguments to pass to \code{bigspline} for fitting mean of each SPC.
#' @param constVar Logical indicating whether to assume constant variance as a function
#' of the periodic variable.
#' @param fitVarArgs List of arguments to pass to \code{bigspline} for fitting variance of each SPC.
#' Unused if \code{constVar==TRUE}.
#' @param nTime Number of time-points by which to discretize the time-dependent
#' behavior of each feature. Corresponds to the number of rows in the matrix for
#' which the SPCs will be calculated.
#' @param useSpc Logical indicating whether to use \code{SPC} (default) or \code{svd}.
#' @param sumabsv L1-constraint on the SPCs, passed to \code{SPC}.
#' @param orth Logical indicating whether to require left singular vectors
#' be orthogonal to each other, passed to \code{SPC}.
#' @param nSpc Vector of the number of SPCs to use for prediction. If \code{NA} (default),
#' \code{nSpc} will become \code{1:K}, where \code{K} is the number of SPCs in \code{spcResult}.
#' Each value in \code{nSpc} will correspond to one prediction for each test observation.
#' A value of 2 means that the prediction will be based on the first 2 SPCs.
#' @param betaSv Logical indicating whether to use the singular values of the SPCs
#' as weights in the likelihood calculation.
#' @param timeRange Vector of values of the periodic variable at which to calculate likelihood.
#' The time with the highest likelihood is used as the initial value for the
#' MLE optimizer.
#'
#' @return
#' \item{fitResult}{Result from \code{zeitzeigerFit}}
#' \item{spcResult}{Result from \code{zeitzeigerSpc}}
#' \item{predResult}{Result from \code{zeitzeigerPredict}}
#'
#' @export
zeitzeiger = function(xTrain, timeTrain, xTest, fitMeanArgs=list(rparm=NA), constVar=TRUE, fitVarArgs=list(rparm=NA), nTime=10,
							 useSpc=TRUE, sumabsv=2, orth=TRUE, nSpc=2, betaSv=FALSE, timeRange=seq(0, 1, 0.01)) {
	fitResult = zeitzeigerFit(xTrain, timeTrain, fitMeanArgs)
	spcResult = zeitzeigerSpc(fitResult$xFitMean, fitResult$xFitResid, nTime, useSpc, sumabsv, orth)
	predResult = zeitzeigerPredict(xTrain, timeTrain, xTest, spcResult, fitMeanArgs, constVar, fitVarArgs, nSpc, betaSv,
											 timeRange)
	return(list(fitResult=fitResult, spcResult=spcResult, predResult=predResult))}


#' Train and test a ZeitZeiger predictor, accounting for batch effects.
#'
#' \code{zeitzeigerBatch} trains and tests a predictor on multiple datasets
#' independently, using \code{ComBat} to correct for batch effects prior
#' to running \code{zeitzeiger}. This function requires the \code{metapredict}
#' package.
#'
#' @param ematList Named list of matrices of measurements, one for each dataset,
#' some of which will be for training, others for testing. Each matrix should
#' have rownames corresponding to sample names and colnames corresponding to
#' feature names.
#' @param trainStudyNames Character vector of names in \code{ematList} corresponding
#' to datasets for training.
#' @param sampleMetadata data.frame containing relevant information for each sample
#' across all datasets.
#' @param studyColname Name of column in \code{sampleMetdata} that contains
#' information about which dataset each sample belongs to.
#' @param batchColname Name of column in \code{sampleMetdata} that contains
#' information about which dataset each sample belongs to. This should correspond
#' to the names of \code{ematList}, and will often be the same as
#' \code{studyColname}, but doesn't have to be.
#' @param timeColname Name of column in \code{sampleMetdata} that contains
#' the values of the periodic variable.
#' @param fitMeanArgs List of arguments to pass to \code{bigspline} for fitting mean of each SPC.
#' @param constVar Logical indicating whether to assume constant variance as a function
#' of the periodic variable.
#' @param fitVarArgs List of arguments to pass to \code{bigspline} for fitting variance of each SPC.
#' Unused if \code{constVar==TRUE}.
#' @param nTime Number of time-points by which to discretize the time-dependent
#' behavior of each feature. Corresponds to the number of rows in the matrix for
#' which the SPCs will be calculated.
#' @param useSpc Logical indicating whether to use \code{SPC} (default) or \code{svd}.
#' @param sumabsv L1-constraint on the SPCs, passed to \code{SPC}.
#' @param orth Logical indicating whether to require left singular vectors
#' be orthogonal to each other, passed to \code{SPC}.
#' @param nSpc Vector of the number of SPCs to use for prediction. If \code{NA} (default),
#' \code{nSpc} will become \code{1:K}, where \code{K} is the number of SPCs in \code{spcResult}.
#' Each value in \code{nSpc} will correspond to one prediction for each test observation.
#' A value of 2 means that the prediction will be based on the first 2 SPCs.
#' @param betaSv Logical indicating whether to use the singular values of the SPCs
#' as weights in the likelihood calculation.
#' @param timeRange Vector of values of the periodic variable at which to calculate likelihood.
#' The time with the highest likelihood is used as the initial value for the
#' MLE optimizer.
#' @param covariateName Name of column(s) in \code{sampleMetadata} containing
#' information about other covariates for \code{ComBat}, besides \code{batchColname}.
#' If \code{NA} (default), then there are no other covariates.
#' @param featuresExclude Named list of character vectors corresponding to features
#' to exclude from being used for prediction for the respective test datasets.
#'
#' @return
#' \item{spcResultList}{List of results from \code{zeitzeigerSpc}, one for each test dataset.}
#' \item{timeDepLike}{3-D array of likelihood, with dimensions for each test observation
#' (across all datasets), each element of \code{nSpc}, and each element of \code{timeRange}.}
#' \item{mleFit}{List (for each element in \code{nSpc}) of lists (for each test observation)
#' of \code{mle2} objects.}
#' \item{timePred}{Matrix of predicted times for test observations by values of \code{nSpc}.}
#'
#' @export
zeitzeigerBatch = function(ematList, trainStudyNames, sampleMetadata, studyColname, batchColname, timeColname,
									fitMeanArgs=list(rparm=NA), constVar=TRUE, fitVarArgs=list(rparm=NA), nTime=10, useSpc=TRUE,
									sumabsv=2, orth=TRUE, nSpc=2, betaSv=FALSE, timeRange=seq(0, 1, 0.01), covariateName=NA,
									featuresExclude=NULL) {
	if (!requireNamespace('metapredict', quietly=TRUE)) {
		stop('This function requires the metapredict package. Please see https://github.com/jakejh/metapredict.', call.=FALSE)}

	testStudyNames = names(ematList)[!(names(ematList) %in% trainStudyNames)]

	batchResult = foreach(testStudyName=testStudyNames) %dopar% {
		ematListNow = ematList[c(trainStudyNames, testStudyName)]
		if (!is.null(featuresExclude) && !is.na(featuresExclude[[testStudyName]])) {
			ematListNow = lapply(ematListNow, function(emat) emat[!(rownames(emat) %in% featuresExclude[[testStudyName]]),])}

		ematMerged = metapredict::mergeStudyData(ematListNow, sampleMetadata, batchColname, covariateName)

		idxTrain = sampleMetadata[colnames(ematMerged), studyColname] %in% trainStudyNames
		xTrain = t(ematMerged[,idxTrain])
		xTest = t(ematMerged[,!idxTrain])
		timeTrain = sampleMetadata[colnames(ematMerged)[idxTrain], timeColname]

		fitResult = zeitzeigerFit(xTrain, timeTrain, fitMeanArgs)
		spcResult = zeitzeigerSpc(fitResult$xFitMean, fitResult$xFitResid, nTime, useSpc, sumabsv, orth)
		predResult = zeitzeigerPredict(xTrain, timeTrain, xTest, spcResult, fitMeanArgs, constVar, fitVarArgs, nSpc, betaSv, timeRange)

		dimnames(predResult$timeDepLike)[[1]] = rownames(xTest)
		for (ii in 1:length(nSpc)) {
			names(predResult$mleFit[[ii]]) = rownames(xTest)}
		rownames(predResult$timePred) = rownames(xTest)

		return(list(spcResult=spcResult, predResult=predResult))}

	spcResultList = lapply(batchResult, function(x) x$spcResult)
	names(spcResultList) = testStudyNames

	timeDepLikeList = lapply(batchResult, function(x) x$predResult$timeDepLike)
	timeDepLike = do.call(abind::abind, list(timeDepLikeList, along=1))

	mleFitList = lapply(batchResult, function(x) x$predResult$mleFit) # list (by batch) of lists (by nSpc) of lists (by obs)
	mleFit = list() # list (by nSpc) of lists (by obs)
	for (ii in 1:length(nSpc)) {
		mleFit[[ii]] = do.call(c, lapply(mleFitList, function(x) x[[ii]]))}

	timePred = do.call(rbind, lapply(batchResult, function(x) x$predResult$timePred))

	return(list(spcResultList=spcResultList, timeDepLike=timeDepLike, mleFit=mleFit, timePred=timePred))}
