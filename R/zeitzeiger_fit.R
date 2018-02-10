#' @importFrom bigsplines bigspline
#' @importFrom bigsplines predict.bigspline
#' @importFrom foreach foreach
#' @importFrom foreach "%dopar%"
#' @importFrom foreach "%do%"
#' @importFrom magrittr "%>%"
NULL

globalVariables(c('fitResult', 'foldidNow', 'group', 'ii', 'jj', 'predResultFit', 'spcResult',
						'testStudyName', 'timeMin')) # because codetools doesn't understand foreach


#' Estimate time-dependent mean.
#'
#' \code{zeitzeigerFit} estimates the mean of each feature as a function of the
#' periodic variable, using a periodic smoothing spline.
#'
#' @param x Matrix of measurements, with observations in rows and features in columns.
#' Missing values are allowed.
#' @param time Vector of values of the periodic variable for the observations, where 0
#' corresponds to the lowest possible value and 1 corresponds to the highest possible value.
#' @param fitMeanArgs List of arguments to pass to \code{\link[bigsplines]{bigspline}}.
#' @param dopar Logical indicating whether to process features in parallel.
#' Use \code{\link[doParallel]{registerDoParallel}} to register the parallel backend.
#'
#' @return
#' \item{xFitMean}{List of results from \code{bigspline}. Length is number of columns in \code{x}.}
#' \item{xFitResid}{Matrix of residuals, same dimensions as \code{x}.}
#'
#' @seealso \code{\link[bigsplines]{bigspline}}, \code{\link{zeitzeigerSpc}}, \code{\link{zeitzeigerPredict}},
#' \code{\link{zeitzeigerProp}}
#'
#' @export
zeitzeigerFit = function(x, time, fitMeanArgs=list(rparm=NA, nknots=3), dopar=FALSE) {
	doOp = ifelse(dopar, `%dopar%`, `%do%`)
	idx = !is.na(x)
	resultList = doOp(foreach(jj=1:ncol(x)), {
		xFitMean = do.call(bigspline, c(list(time[idx[,jj]], x[idx[,jj], jj], type='per', xmin=0, xmax=1),
												  fitMeanArgs))
		xFitResid = predict.bigspline(xFitMean, newdata=time) - x[,jj]
		list(xFitMean, xFitResid)})
	xFitMean = lapply(resultList, function(a) a[[1]])
	xFitResid = do.call(cbind, lapply(resultList, function(a) a[[2]]))
	return(list(xFitMean=xFitMean, xFitResid=xFitResid))}


#' Estimate peaks and troughs.
#'
#' \code{zeitzeigerExtrema} estimates the extremum (peak or trough) for each feature
#' by using \code{\link[stats]{optimize}} and the periodic spline fit. This function
#' is now deprecated, and it is recommended to instead use \code{\link{zeitzeigerProp}}.
#'
#' @param fitResult Output of \code{zeitzeigerFit}.
#' @param maximum Logical indicating whether to find maximum or minimum.
#' @param dopar Logical indicating whether to process features in parallel.
#' Use \code{\link[doParallel]{registerDoParallel}} to register the parallel backend.
#'
#' @return Matrix with a row for each feature and columns for location and value.
#'
#' @seealso \code{\link{zeitzeigerFit}}, \code{\link{zeitzeigerProp}}
#'
#' @export
zeitzeigerExtrema = function(fitResult, maximum=TRUE, dopar=TRUE) {
	doOp = ifelse(dopar, `%dopar%`, `%do%`)
	extrema = doOp(foreach(ii=1:length(fitResult$xFitMean), .combine=rbind), {
		f = function(time) predict.bigspline(fitResult$xFitMean[[ii]], newdata=time)
		optResult = stats::optimize(f, interval=c(0, 1), maximum=maximum)
		matrix(do.call(c, optResult), nrow=1, dimnames=list(NULL, c('location', 'value')))})}


#' Calculate the signal-to-noise of the periodic spline fits.
#'
#' \code{zeitzeigerSnr} calculates the signal-to-noise ratio of the spline fit
#' for each feature. The SNR is calculated as the peak-to-trough difference,
#' divided by the square root of the mean of the squared residuals. This function
#' is now deprecated, and it is recommended to instead use \code{\link{zeitzeigerProp}}.
#'
#' @param fitResult Output of \code{zeitzeigerFit}.
#' @param dopar Logical indicating whether to process features in parallel.
#' Use \code{\link[doParallel]{registerDoParallel}} to register the parallel backend.
#'
#' @return Vector of signal-to-noise values.
#'
#' @seealso \code{\link{zeitzeigerFit}}, \code{\link{zeitzeigerProp}}
#'
#' @export
zeitzeigerSnr = function(fitResult, dopar=TRUE) {
	maxVal = zeitzeigerExtrema(fitResult, dopar=dopar)[,'value']
	minVal = zeitzeigerExtrema(fitResult, maximum=FALSE, dopar=dopar)[,'value']
	return((maxVal - minVal) / sqrt(colMeans(fitResult$xFitResid^2)))}


#' Calculate the rhythmic properties of each feature.
#'
#' \code{zeitzeigerProp} calculates the rhythmic properties of each feature's
#' spline fit: location and value of peak, location and value of trough,
#' amplitude measured peak to trough, and signal-to-noise ratio (amplitude
#' divided by the square root of the mean of the squared residuals).
#'
#' @param fitResult Output of \code{zeitzeigerFit}.
#' @param dopar Logical indicating whether to process features in parallel.
#' Use \code{\link[doParallel]{registerDoParallel}} to register the parallel backend.
#'
#' @return Data frame with a row for each feature.
#'
#' @seealso \code{\link{zeitzeigerFit}}
#'
#' @export
zeitzeigerProp = function(fitResult, dopar=TRUE) {
	peak = zeitzeigerExtrema(fitResult, dopar=dopar)
	trough = zeitzeigerExtrema(fitResult, maximum=FALSE, dopar=dopar)
	amp = peak[,'value'] - trough[,'value']
	snr = amp / sqrt(colMeans(fitResult$xFitResid^2))
	return(data.frame(snr = snr, amp = amp, peakLoc = peak[,'location'], peakVal = peak[,'value'],
							troughLoc = trough[,'location'], troughVal = trough[,'value']))}


#' Estimate significance of periodicity by permutation testing.
#'
#' \code{zeitzeigerSig} estimates the statistical significance of the periodic
#' smoothing spline fit. At each permutation, the time vector is scrambled and then
#' zeitzeigerFit is used to fit a periodic smoothing spline for each feature as a
#' function of time. The p-value for each feature is calculated based on the
#' of permutations that had a signal-to-noise ratio at least as large as the
#' observed signal-to-noise ratio, adjusted by the method of Phipson and Smyth (2010).
#' Make sure to first register the parallel backend using \code{registerDoParallel}.
#' For genome-scale data, this function will be very slow.
#'
#' @param x Matrix of measurements, with observations in rows and features in columns.
#' Missing values are allowed.
#' @param time Vector of values of the periodic variable for the observations, where 0
#' corresponds to the lowest possible value and 1 corresponds to the highest possible value.
#' @param fitMeanArgs List of arguments to pass to \code{bigspline}.
#' @param nIter Number of permutations.
#' @param dopar Logical indicating whether to process features in parallel.
#' Use \code{\link[doParallel]{registerDoParallel}} to register the parallel backend.
#'
#' @return Vector of p-values.
#'
#' @seealso \code{\link{zeitzeigerFit}}
#'
#' @export
zeitzeigerSig = function(x, time, fitMeanArgs=list(rparm=NA, nknots=3), nIter=200, dopar=TRUE) {
	doOp = ifelse(dopar, `%dopar%`, `%do%`)
	timeIdx = do.call(rbind, lapply(1:nIter, function(x) sample.int(length(time))))
	snrRand = doOp(foreach(ii=1:nIter, .combine=rbind), {
		timeRand = time[timeIdx[ii,]]
		fitResult = zeitzeigerFit(x, timeRand, fitMeanArgs)
		zeitzeigerSnr(fitResult)})
	fitResult = zeitzeigerFit(x, time, fitMeanArgs)
	snr = zeitzeigerSnr(fitResult, dopar=dopar)
	snrMat = matrix(rep(snr, nIter), nrow=nIter, byrow=TRUE)
	# assume x values within a feature are unique, find number of unique permutations of time vector
	totalNperm = factorial(length(time)) / prod(factorial(table(time)))
	sig = statmod::permp(colSums(snrMat <= snrRand), nIter, total.nperm=totalNperm, twosided=FALSE)
	return(sig)}


zeitzeigerFitVar = function(time, xFitResid, constVar=TRUE, fitVarArgs=list(rparm=NA)) {
	# length(time): n
	# dim(xFitResid): c(n, p)
	warnOrig = getOption('warn')
	options(warn=-1)
	xFitVar = list()
	if (constVar) {
		sigmaAll = colMeans(xFitResid^2, na.rm=TRUE)
		for (jj in 1:ncol(xFitResid)) {
			xFitVar[[jj]] = bigspline(c(0, 0.3, 0.7), rep(sigmaAll[jj], 3), type='per', xmin=0, xmax=1,
											  rparm=NA, nknots=3)}
	} else {
		idx = !is.na(xFitResid)
		for (jj in 1:ncol(xFitResid)) {
			# todo: fix so that variance can't be less than zero
			xFitVar[[jj]] = do.call(bigspline, c(list(time[idx[,jj]], xFitResid[idx[,jj], jj]^2, type='per',
																	xmin=0, xmax=1), fitVarArgs))}}
	options(warn=warnOrig)
	return(xFitVar)}


#' Calculate sparse principal components of time-dependent variation.
#'
#' \code{zeitzeigerSpc} calculates the sparse principal components (SPCs),
#' given the time-dependent means and the residuals from \code{zeitzeigerFit}.
#' This function calls \code{PMA::SPC}.
#'
#' @param xFitMean List of bigsplines, length is number of features.
#' @param xFitResid Matrix of residuals, dimensions are observations by features.
#' @param nTime Number of time-points by which to discretize the time-dependent
#' behavior of each feature. Corresponds to the number of rows in the matrix for
#' which the SPCs will be calculated.
#' @param useSpc Logical indicating whether to use \code{SPC} (default) or \code{svd}.
#' @param sumabsv L1-constraint on the SPCs, passed to \code{SPC}.
#' @param orth Logical indicating whether to require left singular vectors
#' be orthogonal to each other, passed to \code{SPC}.
#' @param ... Other arguments passed to \code{SPC}.
#'
#' @return Result from \code{SPC}, unless \code{useSpc==FALSE}, then result from \code{svd}.
#'
#' @seealso \code{\link{zeitzeigerFit}}, \code{\link{zeitzeigerPredict}}
#'
#' @export
zeitzeigerSpc = function(xFitMean, xFitResid, nTime=10, useSpc=TRUE, sumabsv=1, orth=TRUE, ...) {
	timeRange = seq(0, 1 - 1/nTime, 1/nTime)
	xMean = matrix(data=NA, nrow=length(timeRange), ncol=length(xFitMean))

	for (jj in 1:ncol(xMean)) {
		xMean[,jj] = predict.bigspline(xFitMean[[jj]], newdata=timeRange)}
	xMeanScaled = scale(xMean, center=TRUE, scale=FALSE)
	z = xMeanScaled %*% diag(1/sqrt(colMeans(xFitResid^2)))

	if (useSpc) {
		spcResult = PMA::SPC(z, sumabsv=sumabsv, K=min(dim(z)), orth=orth, trace=FALSE, compute.pve=FALSE, ...)
	} else {
		spcResult = svd(z)}
	return(spcResult)}


fx = function(x, time, xFitMean, xFitVar, logArg=FALSE) {
	# dim(x): c(n, p) or c(1, p)
	# length(time): n
	if (is.vector(x)) {
		x = matrix(x, nrow=1)}
	like = matrix(data=NA, nrow=length(time), ncol=ncol(x))
	for (jj in 1:ncol(x)) {
		xPredMean = predict.bigspline(xFitMean[[jj]], newdata=time)
		xPredSd = sqrt(pmax(predict.bigspline(xFitVar[[jj]], newdata=time), 0.0025))
		like[,jj] = stats::dnorm(x[,jj], mean=xPredMean, sd=xPredSd, log=logArg)}
	return(like)}


# #' Calculate time-dependent likelihood.
# #'
# #' Given a matrix of test observations, the estimated time-dependent means
# #' and variances of the features, and a vector of times, \code{zeitzeigerLikelihood}
# #' calculates the likelihood of each time for each test observation. The calculation
# #' assumes that conditioned on the periodic variable, the densities of the features
# #' are normally distributed. The calculation also assumes that the features are
# #' independent.
# #'
# #' @param xTest Matrix of measurements, with observations in rows and features in columns.
# #' @param xFitMean List of bigsplines for time-dependent mean, length is number of features.
# #' @param xFitVar List of bigsplines for time-dependent variance, length is number of features.
# #' @param beta Vector of coefficients for weighted likelihood. If \code{NA} (default),
# #' then each feature is weighted equally.
# #' @param timeRange Vector of values of the periodic variable at which to calculate likelihood.
# #' @param logArg Logical indicating whether to return likelihood (default) or log-likelihood.
# #'
# #' @return Matrix with observations in rows and times in columns.
# #'
# #' @export
zeitzeigerLikelihood = function(xTest, xFitMean, xFitVar, beta=NA, timeRange=seq(0, 1, 0.01), logArg=FALSE) {
	if (is.na(beta[1])) {
		beta = rep_len(1, length(xFitMean))}
	betaMat = matrix(rep.int(beta, length(timeRange)), nrow=length(timeRange), byrow=TRUE)
	loglike = matrix(NA, nrow=nrow(xTest), ncol=length(timeRange))
	for (ii in 1:nrow(xTest)) {
		xTestNow = xTest[ii,, drop=FALSE]
		loglikeTmp = fx(xTestNow, timeRange, xFitMean, xFitVar, logArg=TRUE)
		loglike[ii,] = rowSums(loglikeTmp * betaMat)}
	if (logArg) {
		return(loglike)
	} else {
		return(exp(loglike))}}


#' Calculate difference between values of a periodic variable.
#'
#' \code{calcTimeDiff} calculates the difference between values
#' of a periodic variable in a sensible way, making the difference
#' as close to zero as possible.
#'
#' @param time1 Vector.
#' @param time2 Vector (same length as \code{time1}) or matrix
#' (number of rows equal to length of \code{time1}). If \code{time2} is a
#' matrix, \code{time1} is expanded to have the same number of columns.
#' @param timeMax Maximum value of the periodic variable, i.e., the value
#' that is equivalent to zero. Typically, all values in \code{time1} and
#' \code{time2} should be between 0 and \code{timeMax}.
#'
#' @return Vector or matrix corresponding to \code{time2 - time1}.
#'
#' @export
calcTimeDiff = function(time1, time2, timeMax=1) {
	time1 = time1 / timeMax
	time2 = time2 / timeMax
	if (is.vector(time2)) {
		d = unname(time2 - time1)
		d2 = cbind(d, d-1, d+1)
		d3 = apply(d2, MARGIN=1, function(x) x[which.min(abs(x))])
		d4 = sapply(d3, function(x) ifelse(length(x)==0, NA, x))
	} else {
		d4 = time2
		for (jj in 1:ncol(time2)) {
			d = unname(time2[,jj] - time1)
			d2 = cbind(d, d-1, d+1)
			d3 = apply(d2, MARGIN=1, function(x) x[which.min(abs(x))])
			d4[,jj] = sapply(d3, function(x) ifelse(length(x)==0, NA, x))}}
	return(d4 * timeMax)}
