#' Fit a periodic spline for each feature on cross-validation
#'
#' Fit a periodic spline for each feature for each fold of cross-validation.
#'
#' @param x Matrix of measurements, with observations in rows and features in
#'   columns.
#' @param time Vector of values of the periodic variable for the observations,
#'   where 0 corresponds to the lowest possible value and 1 corresponds to the
#'   highest possible value.
#' @param foldid Vector of values indicating the fold to which each observation
#'   belongs.
#' @param nKnots Number of internal knots to use for the periodic smoothing
#'   spline.
#'
#' @return A list consisting of the result from [zeitzeigerFit()] for each fold.
#'
#' @seealso [zeitzeigerFit()], [zeitzeigerSpcCv()], [zeitzeigerPredictCv()]
#'
#' @export
zeitzeigerFitCv = function(x, time, foldid, nKnots = 3) {
  foldidNow = NULL
  fitResultList = foreach(foldidNow = sort(unique(foldid))) %do% {
    idxTrain = foldid != foldidNow
    zeitzeigerFit(x[idxTrain, ], time[idxTrain], nKnots)}
  return(fitResultList)}


#' Calculate sparse principal components of time-dependent variation on cross-validation
#'
#' Calculate SPCs for each fold of cross-validation.
#'
#' @param fitResultList Output of [zeitzeigerFitCv()].
#' @param nTime Number of time-points by which to discretize the time-dependent
#'   behavior of each feature. Corresponds to the number of rows in the matrix
#'   for which the SPCs will be calculated.
#' @param useSpc Logical indicating whether to use `SPC` (default) or `svd`.
#' @param sumabsv L1-constraint on the SPCs, passed to `SPC`.
#' @param orth Logical indicating whether to require left singular vectors
#'   be orthogonal to each other, passed to `SPC`.
#' @param dopar Logical indicating whether to process the folds in parallel. Use
#'   [doParallel::registerDoParallel()] to register the parallel backend.
#'
#' @return A list consisting of the result from [zeitzeigerSpc()] for each fold.
#'
#' @seealso [zeitzeigerSpc()], [zeitzeigerFitCv()], [zeitzeigerPredictCv()]
#'
#' @export
zeitzeigerSpcCv = function(
  fitResultList, nTime = 10, useSpc = TRUE, sumabsv = 1, orth = TRUE, dopar = TRUE) {
  fitResult = NULL
  doOp = ifelse(dopar, `%dopar%`, `%do%`)
  spcResultList = doOp(foreach(fitResult = fitResultList), {
    zeitzeigerSpc(fitResult$xFitMean, fitResult$xFitResid, nTime, useSpc, sumabsv, orth)})
  return(spcResultList)}


#' Predict corresponding time for observations on cross-validation
#'
#' Make predictions for each observation for each fold of cross-validation.
#'
#' @param x Matrix of measurements, observations in rows and features in
#'   columns.
#' @param time Vector of values of the periodic variable for observations, where
#'   0 corresponds to the lowest possible value and 1 corresponds to the highest
#'   possible value.
#' @param foldid Vector of values indicating the fold to which each observation
#'   belongs.
#' @param spcResultList Output of [zeitzeigerSpcCv()].
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
#' @param dopar Logical indicating whether to process the folds in parallel.
#' Use [doParallel::registerDoParallel()] to register the parallel backend.
#'
#' @return A list of the same structure as [zeitzeigerPredict()], combining the
#'   results from each fold of cross-validation.
#' \item{timeDepLike}{3-D array of likelihood, with dimensions for each
#'   observation, each element of `nSpc`, and each element of `timeRange`.}
#' \item{mleFit}{List (for each element in `nSpc`) of lists (for each
#'   observation) of `mle2` objects.}
#' \item{timePred}{Matrix of predicted times for observations by values of
#'   `nSpc`.}
#'
#' @seealso [zeitzeigerPredict()], [zeitzeigerFitCv()], [zeitzeigerSpcCv()]
#'
#' @export
zeitzeigerPredictCv = function(
  x, time, foldid, spcResultList, nKnots = 3, nSpc = NA,
  timeRange = seq(0, 1 - 0.01, 0.01), dopar = TRUE) {

  spcResult = NULL
  foldidUnique = sort(unique(foldid))
  doOp = ifelse(dopar, `%dopar%`, `%do%`)

  predResultList = doOp(foreach(foldidNow = foldidUnique, spcResult = spcResultList), {
    idxTrain = foldid != foldidNow
    xTrain = x[idxTrain, , drop = FALSE]
    xTest = x[!idxTrain, , drop = FALSE]
    zeitzeigerPredict(xTrain, time[idxTrain], xTest, spcResult, nKnots, nSpc,
                      timeRange)})

  nSpcLen = dim(predResultList[[1]]$timePred)[2]
  timeDepLike = array(NA, dim = c(nrow(x), nSpcLen, length(timeRange)))
  timePred = matrix(NA, nrow = nrow(x), ncol = nSpcLen)
  mleFit = vector('list', nSpcLen)
  mleFit = lapply(mleFit, function(a) vector('list', nrow(x)))
  for (foldidNow in foldidUnique) {
    timeDepLike[foldid == foldidNow, , ] = predResultList[[which(foldidUnique == foldidNow)]]$timeDepLike
    for (ii in 1:nSpcLen) {
      mleFit[[ii]][foldid == foldidNow] = predResultList[[which(foldidUnique == foldidNow)]]$mleFit[[ii]]}
    timePred[foldid == foldidNow, ] = predResultList[[which(foldidUnique == foldidNow)]]$timePred}

  return(list(timeDepLike = timeDepLike, mleFit = mleFit, timePred = timePred))}


#' Predict corresponding time for groups of observations on cross-validation
#'
#' Predict corresponding time for each group of observations in
#' cross-validation. Thus, each fold is equivalent to a group.
#'
#' @param x Matrix of measurements, observations in rows and features in columns.
#' @param time Vector of values of the periodic variable for observations, where
#'   0 corresponds to the lowest possible value and 1 corresponds to the highest
#'   possible value.
#' @param foldid Vector of values indicating the fold to which each observation
#'   belongs.
#' @param spcResultList Result from [zeitzeigerSpcCv()].
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
#' @param dopar Logical indicating whether to process the folds in parallel. Use
#'   [doParallel::registerDoParallel()] to register the parallel backend.
#'
#' @return A list of the same structure as `zeitzeigerPredictGroup`, combining
#'   the results from each fold of cross-validation. Folds (i.e, groups) will be
#'   sorted by foldid.
#' \item{timeDepLike}{3-D array of likelihood, with dimensions for each fold,
#'   each element of `nSpc`, and each element of `timeRange`.}
#' \item{mleFit}{List (for each element in `nSpc`) of lists (for each fold) of
#'   `mle2` objects.}
#' \item{timePred}{Matrix of predicted times for folds by values of `nSpc`.}
#'
#' @seealso [zeitzeigerFitCv()], [zeitzeigerSpcCv()], [zeitzeigerPredictGroup()]
#'
#' @export
zeitzeigerPredictGroupCv = function(
  x, time, foldid, spcResultList, nKnots = 3, nSpc = NA,
  timeRange = seq(0, 1 - 0.01, 0.01), dopar = TRUE) {

  group = timeDiff = foldidNow = spcResult = predResultFit = NULL
  foldidUnique = sort(unique(foldid))
  doOp = ifelse(dopar, `%dopar%`, `%do%`)

  groupDf = data.table(time = time, group = foldid)
  groupDf[, timeDiff := time - min(time), by = 'group']
  # groupDf = data.frame(time = time, group = foldid, stringsAsFactors = FALSE)
  # groupDf = groupDf %>%
  #   dplyr::inner_join(groupDf %>%
  #                       dplyr::group_by(group) %>%
  #                       dplyr::summarize(timeMin = min(time)), by = 'group') %>%
  #   dplyr::mutate(timeDiff = time - timeMin)

  predResultList = doOp(foreach(foldidNow = foldidUnique, spcResult = spcResultList), {
    idxTrain = foldid != foldidNow
    xTrain = x[idxTrain, , drop = FALSE]
    xTest = x[!idxTrain, , drop = FALSE]
    groupTest = groupDf[!idxTrain, ]
    zeitzeigerPredictGroup(xTrain, time[idxTrain], xTest, groupTest,
                           spcResult, nKnots, nSpc, timeRange)})

  timeDepLike = abind::abind(
    lapply(predResultList, function(a) a$timeDepLike), along = 1)
  mleFit = lapply(predResultFit, function(a) a$mleFit)
  timePred = do.call(rbind, lapply(predResultList, function(a) a$timePred))

  return(list(timeDepLike = timeDepLike, mleFit = mleFit, timePred = timePred))}
