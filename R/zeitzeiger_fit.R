#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#' @importFrom foreach %do%
#' @importFrom magrittr %>%
NULL

globalVariables(c('fitResult', 'foldidNow', 'group', 'ii', 'jj', 'predResultFit', 'spcResult',
                  'testStudyName', 'timeMin')) # because codetools doesn't understand foreach


#' Fit a periodic spline for each feature
#'
#' `zeitzeigerFit` fits a periodic smoothing spline to the measurements for each
#' feature as a function of the periodic variable.
#'
#' @param x Matrix of measurements, with observations in rows and features in columns.
#' Missing values are allowed.
#' @param time Vector of values of the periodic variable for the observations, where 0
#' corresponds to the lowest possible value and 1 corresponds to the highest possible value.
#' @param nKnots Number of internal knots to use for the periodic smoothing spline.
#'
#' @return
#' \item{xFitMean}{Matrix of coefficients, where rows correspond to features and
#' columns correspond to variables in the fit.}
#' \item{xFitResid}{Vector of root mean square of residuals, same length as `x`.}
#'
#' @seealso `\link{zeitzeigerSpc}`, `\link{zeitzeigerPredict}`
#'
#' @export
zeitzeigerFit = function(x, time, nKnots = 3) {
  knots = seq(0, 1 - 1 / nKnots, length = nKnots)
  basis = getSplineBasis(time, knots)
  fit = limma::lmFit(t(x), basis)
  xFitMean = stats::coef(fit)
  xFitResid = fit$sigma
  return(list(xFitMean = xFitMean, xFitResid = xFitResid))}


# modified from the pbs package
# Jake Hughey, Nov 2018
getPeriodBsplineKnots = function(knots, degree = 3) {
  nKnots = length(knots)

  closeKnots1 = rep(0, degree)
  for (i in 1:degree) {
    closeKnots1[i] = knots[nKnots] + knots[i + 1] - knots[1]}

  closeKnots2 = rep(0, degree)
  for (i in 1:degree){
    closeKnots2[i] = knots[1] - knots[nKnots] + knots[nKnots - i]}

  closeKnots = c(closeKnots2, knots, closeKnots1)
  return(closeKnots)}


# modified from the pbs package
# Jake Hughey, Nov 2018
getSplineBasis = function(x, knots, period = 1, degree = 3) {
  # basis = cbind(1, cos(x / period * 2 * pi), sin(x / period * 2 * pi))
  x = x %% period
  degree = as.integer(degree)
  knotsFull = getPeriodBsplineKnots(c(0, knots, period), degree)

  basisTmp = splines::spline.des(knotsFull, x, 1 + degree)$design
  basisLeft = basisTmp[, 1:degree, drop = FALSE]
  basisRight = basisTmp[, (ncol(basisTmp) - degree + 1):ncol(basisTmp),
                        drop = FALSE]

  idx = -c(1:degree, (ncol(basisTmp) - degree + 1):ncol(basisTmp))
  basis = cbind(basisTmp[, idx, drop = FALSE], basisLeft + basisRight)
  return(basis)}


#' Calculate time-dependent mean
#'
#' `predictIntensity` calculates the expected value of each.
#'
#' @param fitCoef Matrix of coefficients from the spline fits, where rows
#' correspond to features and columns correspond to variables in the model.
#' @param time Vector of values of the periodic variable for the observations, where 0
#' corresponds to the lowest possible value and 1 corresponds to the highest possible value.
#' @param period Period for the periodic variable.
#' @param knots Optional vector of knots. This argument is designed for internal use.
#'
#' @return Matrix of predicted measurements, where rows correspond to time-points
#' and columns correspond to features.
#'
#' @seealso `\link{zeitzeigerFit}`
#'
#' @export
predictIntensity = function(fitCoef, time, period = 1, knots = NULL) {
  if (is.null(knots)) {
    nKnots = ncol(fitCoef) - 1
    knots = seq(0, period - period / nKnots, length = nKnots)}
  basis = getSplineBasis(time, knots, period)
  y = basis %*% t(fitCoef)
  return(y)}


#' Calculate sparse principal components of time-dependent variation
#'
#' `zeitzeigerSpc` calculates the sparse principal components (SPCs),
#' given the time-dependent means and the residuals from `zeitzeigerFit`.
#' This function calls `PMA::SPC`.
#'
#' @param xFitMean List of bigsplines, length is number of features.
#' @param xFitResid Matrix of residuals, dimensions are observations by features.
#' @param nTime Number of time-points by which to discretize the time-dependent
#' behavior of each feature. Corresponds to the number of rows in the matrix for
#' which the SPCs will be calculated.
#' @param useSpc Logical indicating whether to use `SPC` (default) or `svd`.
#' @param sumabsv L1-constraint on the SPCs, passed to `SPC`.
#' @param orth Logical indicating whether to require left singular vectors
#' be orthogonal to each other, passed to `SPC`.
#' @param ... Other arguments passed to `SPC`.
#'
#' @return Result from `SPC`, unless `useSpc` is `FALSE`, then result from `svd`.
#'
#' @seealso `\link{zeitzeigerFit}`, `\link{zeitzeigerPredict}`
#'
#' @export
zeitzeigerSpc = function(xFitMean, xFitResid, nTime = 10, useSpc = TRUE,
                         sumabsv = 1, orth = TRUE, ...) {
  timeRange = seq(0, 1 - 1 / nTime, 1 / nTime)
  xMean = predictIntensity(xFitMean, timeRange)
  xMeanScaled = scale(xMean, center = TRUE, scale = FALSE)
  z = xMeanScaled / matrix(rep(xFitResid, nTime), nrow = nTime, byrow = TRUE)
  if (useSpc) {
    spcResult = PMA::SPC(z, sumabsv = sumabsv, K = min(dim(z)), orth = orth,
                         trace = FALSE, compute.pve = FALSE, ...)
  } else {
    spcResult = svd(z)}
  return(spcResult)}


fx = function(x, time, xFitMean, xFitResid, knots, logArg = FALSE) {
  # dim(x): c(n, p) or c(1, p)
  # length(time): n
  if (is.vector(x)) {
    x = matrix(x, nrow = 1)}

  like = matrix(data = 0, nrow = length(time), ncol = ncol(x))
  for (jj in 1:ncol(x)) {
    xPredMean = predictIntensity(xFitMean[jj, , drop = FALSE], time, knots = knots)[, 1]
    like[, jj] = stats::dnorm(x[, jj], mean = xPredMean, sd = xFitResid[jj], log = logArg)}
  return(like)}


# #' Calculate time-dependent likelihood
# #'
# #' Given a matrix of test observations, the estimated time-dependent means
# #' and variances of the features, and a vector of times, `zeitzeigerLikelihood`
# #' calculates the likelihood of each time for each test observation. The calculation
# #' assumes that conditioned on the periodic variable, the densities of the features
# #' are normally distributed. The calculation also assumes that the features are
# #' independent.
# #'
# #' @param xTest Matrix of measurements, with observations in rows and features in columns.
# #' @param xFitMean Matrix of coefficients from the spline fits.
# #' @param xFitResid Vector of root mean square residuals.
# #' @param timeRange Vector of values of the periodic variable at which to calculate likelihood.
# #' @param logArg Logical indicating whether to return likelihood (default) or log-likelihood.
# #'
# #' @return Matrix with observations in rows and times in columns.
# #'
# #' @export
zeitzeigerLikelihood = function(xTest, xFitMean, xFitResid, knots,
                                timeRange = seq(0, 1 - 0.01, 0.01),
                                logArg = FALSE) {
  loglike = matrix(NA, nrow = nrow(xTest), ncol = length(timeRange))
  for (ii in 1:nrow(xTest)) {
    xTestNow = xTest[ii, , drop=FALSE]
    loglikeTmp = fx(xTestNow, timeRange, xFitMean, xFitResid, knots, logArg = TRUE)
    loglike[ii, ] = rowSums(loglikeTmp)}
  if (logArg) {
    return(loglike)
  } else {
    return(exp(loglike))}}


#' Calculate circular difference
#'
#' `getCircDiff` calculates the circular difference.
#'
#' @param x Numeric vector or matrix.
#' @param y Numeric vector or matrix.
#' @param period Period of the periodic variable.
#' @param towardZero If `TRUE`, returned values will be between `-period / 2`
#' and `period / 2`. If `FALSE`, returned values will be between 0 and `period`.
#'
#' @return Vector or matrix corresponding to `x - y`.
#'
#' @export
getCircDiff = function(x, y, period = 1, towardZero = TRUE) {
  d = (x - y) %% period
  if (towardZero) {
    d = ifelse(d <= (period / 2), d, d - period)}
  return(d)}
