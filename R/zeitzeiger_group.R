fxGroup = function(x, timeDiff, time, xFitMean, xFitVar, logArg=FALSE) {
	# dim(x): c(m, p)
	# length(timeDiff): m
	# length(time): n
	timeRep = rowSums(expand.grid(time, timeDiff)) %% 1
	xRep = x[rep(1:nrow(x), each=length(time)),, drop=FALSE]
	logLike = matrix(data=NA, nrow=length(time), ncol=ncol(x))

	for (jj in 1:ncol(x)) {
		xPredMean = predict(xFitMean[[jj]], newdata=timeRep)
		xPredSd = sqrt(predict(xFitVar[[jj]], newdata=timeRep)) # safe if constVar==TRUE
		logLikeTmp = dnorm(xRep[,jj], mean=xPredMean, sd=xPredSd, log=TRUE)
		logLike[,jj] = rowSums(matrix(logLikeTmp, nrow=length(time)))}

	if (logArg) {
		return(logLike)
	} else {
		return(exp(logLike))}}


zeitzeigerPredictGivenDensityGroup = function(xTest, groupTest, xFitMean, xFitVar, beta, timeRange=seq(0, 1, 0.01)) {
	groups = groupTest[['group']]
	groupsUnique = sort(unique(groups))
	negLogLike = matrix(NA, nrow=length(groupsUnique), ncol=length(timeRange))
	mleFit = list()

	for (ii in 1:length(groupsUnique)) {
		xTestNow = xTest[groups==groupsUnique[ii],, drop=FALSE]
		timeDiffNow = groupTest[groups==groupsUnique[ii],][['timeDiff']]
		negLogLike[ii,] = -fxGroup(xTestNow, timeDiffNow, timeRange, xFitMean, xFitVar, logArg=TRUE) %*% matrix(beta)

		timeStart = timeRange[which.min(negLogLike[ii,])]
		names(timeStart) = 'time'
		negLogLikeFunc = function(time) -sum(fxGroup(xTestNow, timeDiffNow, time, xFitMean, xFitVar) * beta)
		bbmle::parnames(negLogLikeFunc) = names(timeStart)
		warnOrig = getOption('warn')
		options(warn=-1)
		mleFit[[ii]] = bbmle::mle2(negLogLikeFunc, start=timeStart, vecpar=TRUE, method='L-BFGS-B', lower=0, upper=1)
		options(warn=warnOrig)}

	timePred = sapply(mleFit, function(x) x@coef)
	return(list(timeDepLike=exp(-negLogLike), mleFit=mleFit, timePred=timePred))}


#' Predict corresponding time for groups of test observations.
#'
#' \code{zeitzeigerPredictGroup} predicts the value of the periodic variable
#' for each group of test observations, where the amount of time between each
#' observation in a group is known. This function calls \code{bbmle::mle2}.
#'
#' @param xTrain Matrix of measurements for training data, observations in rows
#' and features in columns.
#' @param timeTrain Vector of values of the periodic variable for training observations,
#' where 0 corresponds to the lowest possible value and 1 corresponds to the highest
#' possible value.
#' @param xTest Matrix of measurements for test data, observations in rows
#' and features in columns.
#' @param groupTest data.frame with one row per observation in \code{xTest}, and columns
#' for \code{group} and \code{timeDiff}. Observations in the same group should have the
#' same value of \code{group}. Within each group, the value of \code{timeDiff}
#' should correspond to the amount of time between that observation and a reference time.
#' Typically, \code{timeDiff} will equal zero for one observation per group.
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
#' The time with the highest likelihood is used as the initial value for the
#' MLE optimizer.
#'
#' @return A list with the following elements, where the groups will be sorted by their names.
#' \item{timeDepLike}{3-D array of likelihood, with dimensions for each group of test observations,
#' each element of \code{nSpc}, and each element of \code{timeRange}.}
#' \item{mleFit}{List (for each element in \code{nSpc}) of lists (for each group of test observations)
#' of \code{mle2} objects.}
#' \item{timePred}{Matrix of predicted times for each group of test observations by values of \code{nSpc}.}
#'
#' @export
zeitzeigerPredictGroup = function(xTrain, timeTrain, xTest, groupTest, spcResult, fitMeanArgs=list(rparm=NA), constVar=TRUE,
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

	timePred = matrix(NA, nrow=length(unique(groupTest[['group']])), length(nSpc))
	timeDepLike = array(NA, dim=c(length(unique(groupTest[['group']])), length(nSpc), length(timeRange)))
	mleFit = list() # list (for each nSpc) of lists (for each group)
	for (ii in 1:length(nSpc)) {
		if (betaSv) {
			beta = spcResult$d[1:nSpc[ii]]
		} else {
			beta = rep_len(1, nSpc[ii])}
		predResult = zeitzeigerPredictGivenDensityGroup(zTest[,1:nSpc[ii], drop=FALSE], groupTest, zFitMean[1:nSpc[ii]],
																		zFitVar[1:nSpc[ii]], beta)
		timeDepLike[,ii,] = predResult$timeDepLike
		mleFit[[ii]] = predResult$mleFit
		timePred[,ii] = predResult$timePred}
	return(list(timeDepLike=timeDepLike, mleFit=mleFit, timePred=timePred))}
