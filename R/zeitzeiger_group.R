fxGroup = function(x, timeDiff, time, xFitMean, xFitResid, knots, logArg = FALSE) {
  # dim(x): c(m, p)
  # length(timeDiff): m
  # length(time): n
  timeRep = rowSums(expand.grid(time, timeDiff)) %% 1
  xRep = x[rep(1:nrow(x), each=length(time)), , drop = FALSE]
  logLike = matrix(data = NA, nrow = length(time), ncol = ncol(x))

  for (jj in 1:ncol(x)) {
    xPredMean = predictIntensity(
      xFitMean[jj, , drop = FALSE], timeRep, knots = knots)[, 1]
    logLikeTmp = stats::dnorm(
      xRep[, jj], mean = xPredMean, sd = xFitResid[jj], log = TRUE)
    logLike[, jj] = rowSums(matrix(logLikeTmp, nrow = length(time)))}

  if (logArg) {
    return(logLike)
  } else {
    return(exp(logLike))}}


zeitzeigerPredictGivenDensityGroup = function(
  xTest, groupTest, xFitMean, xFitResid, knots, timeRange) {
  groups = groupTest$group
  groupsUnique = sort(unique(groups))
  negLogLike = matrix(NA, nrow = length(groupsUnique), ncol = length(timeRange))
  mleFit = list()

  for (ii in 1:length(groupsUnique)) {
    xTestNow = xTest[groups == groupsUnique[ii], , drop=FALSE]
    timeDiffNow = groupTest[groups == groupsUnique[ii], ]$timeDiff
    negLogLike[ii, ] = -fxGroup(xTestNow, timeDiffNow, timeRange, xFitMean,
                                xFitResid, knots, logArg = TRUE)

    timeStart = timeRange[which.min(negLogLike[ii, ])]
    names(timeStart) = 'time'
    negLogLikeFunc = function(time) -sum(fxGroup(xTestNow, timeDiffNow, time,
                                                 xFitMean, xFitResid, knots,
                                                 logArg = TRUE))
    bbmle::parnames(negLogLikeFunc) = names(timeStart)
    warnOrig = getOption('warn')
    options(warn = -1)
    mleFit[[ii]] = bbmle::mle2(negLogLikeFunc, start = timeStart, vecpar = TRUE,
                               method = 'L-BFGS-B', lower = 0, upper = 1)
    options(warn = warnOrig)}

  timePred = sapply(mleFit, function(x) x@coef)
  return(list(timeDepLike = exp(-negLogLike), mleFit = mleFit, timePred = timePred))}


#' Predict corresponding time for groups of test observations
#'
#' Predict the value of the periodic variable for each group of test
#' observations, where the amount of time between each observation in a group is
#' known.
#'
#' @param xTrain Matrix of measurements for training data, observations in rows
#'   and features in columns.
#' @param timeTrain Vector of values of the periodic variable for training
#'   observations, where 0 corresponds to the lowest possible value and 1
#'   corresponds to the highest possible value.
#' @param xTest Matrix of measurements for test data, observations in rows
#'   and features in columns.
#' @param groupTest data.frame with one row per observation in `xTest`, and
#'   columns for `group` and `timeDiff`. Observations in the same group should
#'   have the same value of `group`. Within each group, the value of `timeDiff`
#'   should correspond to the amount of time between that observation and a
#'   reference time. Typically, `timeDiff` will equal zero for one observation
#'   per group.
#' @param spcResult Output of [zeitzeigerSpc()].
#' @param nKnots Number of internal knots to use for the periodic smoothing
#'   spline.
#' @param nSpc Vector of the number of SPCs to use for prediction. If `NA`
#'   (default), `nSpc` will become `1:K`, where `K` is the number of SPCs in
#'   `spcResult`. Each value in `nSpc` will correspond to one prediction for
#'   each test observation. A value of 2 means that the prediction will be based
#'   on the first 2 SPCs.
#' @param timeRange Vector of values of the periodic variable at which to
#'   calculate likelihood. The time with the highest likelihood is used as the
#'   initial value for the MLE optimizer.
#'
#' @return A list with the following elements, where the groups will be sorted
#'   by their names.
#' \item{timeDepLike}{3-D array of likelihood, with dimensions for each group of
#'   test observations, each element of `nSpc`, and each element of
#'   `timeRange`.}
#' \item{mleFit}{List (for each element in `nSpc`) of lists (for each group of
#'   test observations) of `mle2` objects.}
#' \item{timePred}{Matrix of predicted times for each group of test observations
#'   by values of `nSpc`.}
#'
#' @seealso [zeitzeigerPredict()]
#'
#' @export
zeitzeigerPredictGroup = function(
  xTrain, timeTrain, xTest, groupTest, spcResult, nKnots = 3, nSpc = NA,
  timeRange = seq(0, 1 - 0.01, 0.01)) {

  zTrain = xTrain %*% spcResult$v
  zTest = xTest %*% spcResult$v

  zFitResult = zeitzeigerFit(zTrain, timeTrain, nKnots)
  zFitMean = zFitResult$xFitMean
  zFitResid = zFitResult$xFitResid

  knots = seq(0, 1 - 1 / nKnots, length = nKnots)

  if (length(nSpc)==1 && is.na(nSpc)) {
    nSpc = 1:length(spcResult$d)
  } else {
    if (!all(nSpc %in% 1:length(spcResult$d)) || anyDuplicated(nSpc) > 0) {
      stop('nSpc must be unique integers in 1 to the number of singular vectors.')}}

  timePred = matrix(NA, nrow = length(unique(groupTest$group)), length(nSpc))
  timeDepLike = array(
    NA, dim = c(length(unique(groupTest$group)), length(nSpc), length(timeRange)))
  mleFit = list() # list (for each nSpc) of lists (for each group)
  for (ii in 1:length(nSpc)) {
    predResult = zeitzeigerPredictGivenDensityGroup(
      zTest[, 1:nSpc[ii], drop = FALSE], groupTest,
      zFitMean[1:nSpc[ii], , drop = FALSE], zFitResid[1:nSpc[ii]], knots,
      timeRange)
    timeDepLike[, ii, ] = predResult$timeDepLike
    mleFit[[ii]] = predResult$mleFit
    timePred[, ii] = predResult$timePred}
  return(list(timeDepLike = timeDepLike, mleFit = mleFit, timePred = timePred))}
