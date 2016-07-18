#' Estimate time-dependent mean on cross-validation.
#'
#' \code{zeitzeigerFitCv} calls \code{zeitzeigerFit} for each fold of cross-validation.
#' By default, if a parallel backend is registered, this function processes the folds in parallel.
#'
#' @param x Matrix of measurements, with observations in rows and features in columns.
#' @param time Vector of values of the periodic variable for the observations, where 0
#' corresponds to the lowest possible value and 1 corresponds to the highest possible value.
#' @param foldid Vector of values indicating which fold each observation is in.
#' @param fitMeanArgs List of arguments to pass to \code{bigspline}.
#' @param dopar Logical indicating whether to process the folds in parallel. Before calling
#' this function, use \code{doParallel::registerDoParallel} to register the parallel backend.
#'
#' @return A list consisting of the result from \code{zeitzeigerFit} for each fold.
#'
#' @seealso \code{\link{zeitzeigerFit}}, \code{\link{zeitzeigerSpcCv}}, \code{\link{zeitzeigerPredictCv}}
#'
#' @export
zeitzeigerFitCv = function(x, time, foldid, fitMeanArgs=list(rparm=NA), dopar=TRUE) {
	foldidUnique = sort(unique(foldid))
	doOp = ifelse(dopar, `%dopar%`, `%do%`)
	fitResultList = doOp(foreach(foldidNow=foldidUnique), {
		idxTrain = foldid!=foldidNow
		return(zeitzeigerFit(x[idxTrain,], time[idxTrain], fitMeanArgs))})
	return(fitResultList)}


#' Calculate sparse principal components of time-dependent variation on cross-validation.
#'
#' \code{zeitzeigerSpcCv} calls \code{zeitzeigerSpc} for each fold of cross-validation.
#' By default, if a parallel backend is registered, this function processes the folds in parallel.
#'
#' @param fitResultList Result from \code{zeitzeigerFitCv}.
#' @param nTime Number of time-points by which to discretize the time-dependent
#' behavior of each feature. Corresponds to the number of rows in the matrix for
#' which the SPCs will be calculated.
#' @param useSpc Logical indicating whether to use \code{SPC} (default) or \code{svd}.
#' @param sumabsv L1-constraint on the SPCs, passed to \code{SPC}.
#' @param orth Logical indicating whether to require left singular vectors
#' be orthogonal to each other, passed to \code{SPC}.
#' @param dopar Logical indicating whether to process the folds in parallel. Before calling
#' this function, use \code{doParallel::registerDoParallel} to register the parallel backend.
#'
#' @return A list consisting of the result from \code{zeitzeigerSpc} for each fold.
#'
#' @seealso \code{\link{zeitzeigerSpc}}, \code{\link{zeitzeigerFitCv}}, \code{\link{zeitzeigerPredictCv}}
#'
#' @export
zeitzeigerSpcCv = function(fitResultList, nTime=10, useSpc=TRUE, sumabsv=1, orth=TRUE, dopar=TRUE) {
	doOp = ifelse(dopar, `%dopar%`, `%do%`)
	spcResultList = doOp(foreach(fitResult=fitResultList), {
		zeitzeigerSpc(fitResult$xFitMean, fitResult$xFitResid, nTime, useSpc, sumabsv, orth)})
	return(spcResultList)}


#' Predict corresponding time for observations on cross-validation.
#'
#' \code{zeitzeigerPredictCv} calls \code{zeitzeigerPredict} for each fold of cross-validation.
#' By default, if a parallel backend is registered, this function processes the folds in parallel.
#'
#' @param x Matrix of measurements, observations in rows and features in columns.
#' @param time Vector of values of the periodic variable for observations, where 0
#' corresponds to the lowest possible value and 1 corresponds to the highest
#' possible value.
#' @param foldid Vector of values indicating which fold each observation is in.
#' @param spcResultList Result from \code{zeitzeigerSpcCv}.
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
#' @param dopar Logical indicating whether to process the folds in parallel. Before calling
#' this function, use \code{doParallel::registerDoParallel} to register the parallel backend.
#'
#' @return A list of the same structure as \code{zeitzeigerPredict}, combining the results
#' from each fold of cross-validation.
#' \item{timeDepLike}{3-D array of likelihood, with dimensions for each observation,
#' each element of \code{nSpc}, and each element of \code{timeRange}.}
#' \item{mleFit}{List (for each element in \code{nSpc}) of lists (for each observation)
#' of \code{mle2} objects.}
#' \item{timePred}{Matrix of predicted times for observations by values of \code{nSpc}.}
#'
#' @seealso \code{\link{zeitzeigerPredict}}, \code{\link{zeitzeigerFitCv}}, \code{\link{zeitzeigerSpcCv}}
#'
#' @export
zeitzeigerPredictCv = function(x, time, foldid, spcResultList, fitMeanArgs=list(rparm=NA), constVar=TRUE,
										 fitVarArgs=list(rparm=NA), nSpc=NA, betaSv=FALSE, timeRange=seq(0, 1, 0.01),
										 dopar=TRUE) {
	foldidUnique = sort(unique(foldid))
	doOp = ifelse(dopar, `%dopar%`, `%do%`)

	predResultList = doOp(foreach(foldidNow=foldidUnique, spcResult=spcResultList), {
		idxTrain = foldid!=foldidNow
		xTrain = x[idxTrain,, drop=FALSE]
		xTest = x[!idxTrain,, drop=FALSE]
		zeitzeigerPredict(xTrain, time[idxTrain], xTest, spcResult, fitMeanArgs, constVar,
								fitVarArgs, nSpc, betaSv, timeRange)})

	nSpcLen = dim(predResultList[[1]]$timePred)[2]
	timeDepLike = array(NA, dim=c(nrow(x), nSpcLen, length(timeRange)))
	timePred = matrix(NA, nrow=nrow(x), ncol=nSpcLen)
	mleFit = vector('list', nSpcLen)
	mleFit = lapply(mleFit, function(a) vector('list', nrow(x)))
	for (foldidNow in foldidUnique) {
		timeDepLike[foldid==foldidNow,,] = predResultList[[which(foldidUnique==foldidNow)]]$timeDepLike
		for (ii in 1:nSpcLen) {
			mleFit[[ii]][foldid==foldidNow] = predResultList[[which(foldidUnique==foldidNow)]]$mleFit[[ii]]}
		timePred[foldid==foldidNow,] = predResultList[[which(foldidUnique==foldidNow)]]$timePred}

	return(list(timeDepLike=timeDepLike, mleFit=mleFit, timePred=timePred))}


#' Predict corresponding time for groups of observations on cross-validation.
#'
#' \code{zeitzeigerPredictGroupCv} calls \code{zeitzeigerPredictGroup} for each fold of
#' cross-validation. Thus, each fold is equivalent to a group. By default, if a parallel
#' backend is registered, this function processes the folds in parallel.
#'
#' @param x Matrix of measurements, observations in rows and features in columns.
#' @param time Vector of values of the periodic variable for observations, where 0
#' corresponds to the lowest possible value and 1 corresponds to the highest
#' possible value.
#' @param foldid Vector of values indicating which fold each observation is in.
#' @param spcResultList Result from \code{zeitzeigerSpcCv}.
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
#' @param dopar Logical indicating whether to process the folds in parallel. Before calling
#' this function, use \code{doParallel::registerDoParallel} to register the parallel backend.
#'
#' @return A list of the same structure as \code{zeitzeigerPredictGroup}, combining the
#' results from each fold of cross-validation. Folds (i.e, groups) will be sorted by foldid.
#' \item{timeDepLike}{3-D array of likelihood, with dimensions for each fold,
#' each element of \code{nSpc}, and each element of \code{timeRange}.}
#' \item{mleFit}{List (for each element in \code{nSpc}) of lists (for each fold)
#' of \code{mle2} objects.}
#' \item{timePred}{Matrix of predicted times for folds by values of \code{nSpc}.}
#'
#' @seealso \code{\link{zeitzeigerFitCv}}, \code{\link{zeitzeigerSpcCv}}, \code{\link{zeitzeigerPredictGroup}}
#'
#' @export
zeitzeigerPredictGroupCv = function(x, time, foldid, spcResultList, fitMeanArgs=list(rparm=NA),
												constVar=TRUE, fitVarArgs=list(rparm=NA), nSpc=NA, betaSv=FALSE,
												timeRange=seq(0, 1, 0.01), dopar=TRUE) {
	foldidUnique = sort(unique(foldid))
	doOp = ifelse(dopar, `%dopar%`, `%do%`)

	groupDf = data.frame(time = time, group = foldid, stringsAsFactors=FALSE)
	groupDf = groupDf %>%
		inner_join(groupDf %>% group_by(group) %>% summarize(timeMin = min(time)), by='group') %>%
		mutate(timeDiff = time - timeMin)

	predResultList = doOp(foreach(foldidNow=foldidUnique, spcResult=spcResultList), {
		idxTrain = foldid!=foldidNow
		xTrain = x[idxTrain,, drop=FALSE]
		xTest = x[!idxTrain,, drop=FALSE]
		groupTest = groupDf[!idxTrain,]
		zeitzeigerPredictGroup(xTrain, time[idxTrain], xTest, groupTest, spcResult, fitMeanArgs,
									  constVar, fitVarArgs, nSpc, betaSv, timeRange)})

	timeDepLike = abind::abind(lapply(predResultList, function(a) a$timeDepLike), along=1)
	mleFit = lapply(predResultFit, function(a) a$mleFit)
	timePred = do.call(rbind, lapply(predResultList, function(a) a$timePred))

	return(list(timeDepLike=timeDepLike, mleFit=mleFit, timePred=timePred))}
