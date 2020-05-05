#' Estimate peaks and troughs
#'
#' DEPRECATED: We recommend instead using
#' [limorhyde](https://github.com/hugheylab/limorhyde), which has support for
#' cosinor and periodic splines, in combination with methods such as
#' [limma](https://doi.org/doi:10.18129/B9.bioc.limma).
#'
#' Estimate the extremum (peak or trough) for each feature by using
#' [stats::optimize()] and the periodic spline fit. Directly calling
#' this function is deprecated. Instead, use [zeitzeigerProp()].
#'
#' @param fitResult Output of [zeitzeigerFit()].
#' @param maximum Logical indicating whether to find maximum or minimum.
#' @param dopar Logical indicating whether to process features in parallel.
#'   Use [doParallel::registerDoParallel()] to register the parallel backend.
#'
#' @return Matrix with a row for each feature and columns for location and value.
#'
#' @seealso [zeitzeigerFit()], [zeitzeigerProp()]
#'
#' @export
zeitzeigerExtrema = function(fitResult, maximum = TRUE, dopar = TRUE) {
  doOp = ifelse(dopar, `%dopar%`, `%do%`)
  extrema = doOp(foreach(ii = 1:nrow(fitResult$xFitMean), .combine = rbind), {
    f = function(time) predictIntensity(fitResult$xFitMean[ii, , drop = FALSE], time)[, 1]
    optResult = stats::optimize(f, interval = c(0, 1), maximum = maximum)
    matrix(do.call(c, optResult), nrow = 1,
           dimnames = list(NULL, c('location', 'value')))})}


#' Calculate the signal-to-noise of the periodic spline fits
#'
#' DEPRECATED: We recommend instead using
#' [`limorhyde`](https://github.com/hugheylab/limorhyde), which has support for
#' cosinor and periodic splines, in combination with methods such as
#' [`limma`](https://doi.org/doi:10.18129/B9.bioc.limma).
#'
#' Calculate the signal-to-noise ratio of the spline fit for each feature. The
#' SNR is calculated as the peak-to-trough difference, divided by the square
#' root of the mean of the squared residuals.
#'
#' @param fitResult Output of [zeitzeigerFit()].
#' @param dopar Logical indicating whether to process features in parallel.
#'   Use [doParallel::registerDoParallel()] to register the parallel backend.
#'
#' @return Vector of signal-to-noise values.
#'
#' @seealso [zeitzeigerFit()], [zeitzeigerProp()]
#'
#' @export
zeitzeigerSnr = function(fitResult, dopar = TRUE) {
  maxVal = zeitzeigerExtrema(fitResult, dopar = dopar)[, 'value']
  minVal = zeitzeigerExtrema(fitResult, maximum = FALSE, dopar = dopar)[, 'value']
  return((maxVal - minVal) / fitResult$xFitResid)}


#' Calculate the rhythmic properties of each feature
#'
#' DEPRECATED: We recommend instead using
#' [limorhyde](https://github.com/hugheylab/limorhyde), which has support for
#' cosinor and periodic splines, in combination with methods such as
#' [limma](https://doi.org/doi:10.18129/B9.bioc.limma).
#'
#' Calculate the rhythmic properties of each feature's spline fit: location and
#' value of peak, location and value of trough, amplitude measured peak to
#' trough, and signal-to-noise ratio (amplitude divided by the square root of
#' the mean of the squared residuals).
#'
#' @param fitResult Output of [zeitzeigerFit()].
#' @param dopar Logical indicating whether to process features in parallel.
#'   Use [doParallel::registerDoParallel()] to register the parallel backend.
#'
#' @return `data.frame` with a row for each feature.
#'
#' @seealso [zeitzeigerFit()]
#'
#' @export
zeitzeigerProp = function(fitResult, dopar = TRUE) {
  peak = zeitzeigerExtrema(fitResult, dopar = dopar)
  trough = zeitzeigerExtrema(fitResult, maximum = FALSE, dopar = dopar)
  amp = peak[, 'value'] - trough[, 'value']
  snr = amp / fitResult$xFitResid
  d = data.frame(snr = snr, amp = amp,
                 peakLoc = peak[, 'location'], peakVal = peak[, 'value'],
                 troughLoc = trough[, 'location'], troughVal = trough[, 'value'])
  return(d)}


#' Estimate significance of periodicity by permutation testing
#'
#' DEPRECATED: We recommend instead using
#' [limorhyde](https://github.com/hugheylab/limorhyde), which has support for
#' cosinor and periodic splines, in combination with methods such as
#' [limma](https://doi.org/doi:10.18129/B9.bioc.limma).
#'
#' Estimate the statistical significance of the periodic smoothing spline fit.
#' At each permutation, the time vector is scrambled and then [zeitzeigerFit()]
#' is used to fit a periodic smoothing spline for each feature as a function of
#' time. The p-value for each feature is calculated based on the of permutations
#' that had a signal-to-noise ratio at least as large as the observed
#' signal-to-noise ratio, adjusted by the method of
#' [Phipson and Smyth (2010)](https://doi.org/10.2202/1544-6115.1585). Make sure
#' to first register the parallel backend using
#' [doParallel::registerDoParallel()]. For genome-scale data, this will be slow.
#'
#' @param x Matrix of measurements, with observations in rows and features in
#'   columns. Missing values are allowed.
#' @param time Vector of values of the periodic variable for the observations,
#'   where 0 corresponds to the lowest possible value and 1 corresponds to the
#'   highest possible value.
#' @param nKnots Number of internal knots to use for the periodic smoothing
#'   spline.
#' @param nIter Number of permutations.
#' @param dopar Logical indicating whether to process features in parallel. Use
#'   [doParallel::registerDoParallel()] to register the parallel backend.
#'
#' @return Vector of p-values.
#'
#' @seealso [zeitzeigerFit()]
#'
#' @export
zeitzeigerSig = function(x, time, nKnots = 3, nIter = 200, dopar = TRUE) {
  doOp = ifelse(dopar, `%dopar%`, `%do%`)
  timeIdx = do.call(rbind, lapply(1:nIter, function(x) sample.int(length(time))))
  snrRand = doOp(foreach(ii = 1:nIter, .combine = rbind), {
    timeRand = time[timeIdx[ii, ]]
    fitResult = zeitzeigerFit(x, timeRand, nKnots)
    zeitzeigerSnr(fitResult)})
  fitResult = zeitzeigerFit(x, time, nKnots)
  snr = zeitzeigerSnr(fitResult, dopar = dopar)
  snrMat = matrix(rep(snr, nIter), nrow = nIter, byrow = TRUE)
  # assume x values within a feature are unique,
  # find number of unique permutations of time vector
  totalNperm = factorial(length(time)) / prod(factorial(table(time)))
  sig = statmod::permp(colSums(snrMat <= snrRand), nIter,
                       total.nperm = totalNperm, twosided = FALSE)
  return(sig)}
