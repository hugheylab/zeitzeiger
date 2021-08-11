zeitzeigerPredictGivenDensity = function(
  xTest, xFitMean, xFitResid, knots, timeRange) {

  negLoglike = -zeitzeigerLikelihood(
    xTest, xFitMean, xFitResid, knots, timeRange, logArg = TRUE)

  mleFit = list()
  for (ii in 1:nrow(xTest)) {
    xTestNow = xTest[ii, , drop = FALSE]
    timeStart = timeRange[which.min(negLoglike[ii, ])]
    names(timeStart) = 'time'
    negLoglikeFunc = function(time) -sum(fx(xTestNow, time, xFitMean, xFitResid,
                                            knots, logArg = TRUE))
    bbmle::parnames(negLoglikeFunc) = names(timeStart)
    warnOrig = getOption('warn')
    options(warn = -1)
    mleFit[[ii]] = bbmle::mle2(negLoglikeFunc, start = timeStart, vecpar = TRUE,
                               method = 'L-BFGS-B', lower = 0, upper = 1)
    options(warn = warnOrig)}
  timePred = sapply(mleFit, function(x) x@coef)
  return(list(timeDepLike = exp(-negLoglike), mleFit = mleFit, timePred = timePred))}


#' Predict corresponding time for test observations
#'
#' Predict the value of the periodic variable for test observations given
#' training data and SPCs.
#'
#' @param xTrain Matrix of measurements for training data, observations in rows
#'   and features in columns.
#' @param timeTrain Vector of values of the periodic variable for training
#'   observations, where 0 corresponds to the lowest possible value and 1
#'   corresponds to the highest possible value.
#' @param xTest Matrix of measurements for test data, observations in rows and
#'   features in columns.
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
#' @return
#' \item{timeDepLike}{3-D array of likelihood, with dimensions for each test
#'   observation, each element of `nSpc`, and each element of `timeRange`.}
#' \item{mleFit}{List (for each element in `nSpc`) of lists (for each test
#'   observation) of `mle2` objects.}
#' \item{timePred}{Matrix of predicted times for test observations by values of
#'   `nSpc`.}
#'
#' @seealso [zeitzeigerFit()], [zeitzeigerSpc()]
#'
#' @export
zeitzeigerPredict = function(
  xTrain, timeTrain, xTest, spcResult, nKnots = 3, nSpc = NA,
  timeRange = seq(0, 1 - 0.01, 0.01)) {

  zTrain = xTrain %*% spcResult$v
  zTest = xTest %*% spcResult$v

  zFitResult = zeitzeigerFit(zTrain, timeTrain, nKnots)
  zFitMean = zFitResult$xFitMean
  zFitResid = zFitResult$xFitResid

  knots = seq(0, 1 - 1 / nKnots, length = nKnots)

  if (length(nSpc) == 1 && is.na(nSpc)) {
    nSpc = 1:length(spcResult$d)
  } else {
    if (!all(nSpc %in% 1:length(spcResult$d)) || anyDuplicated(nSpc) > 0) {
      stop('nSpc must be unique integers in 1 to the number of singular vectors.')}}

  timePred = matrix(NA, nrow = nrow(xTest), ncol = length(nSpc))
  timeDepLike = array(NA, dim = c(nrow(xTest), length(nSpc), length(timeRange)))
  mleFit = list() # list (for each nSpc) of lists (for each observation)
  for (ii in 1:length(nSpc)) {
    predResult = zeitzeigerPredictGivenDensity(
      zTest[, 1:nSpc[ii], drop = FALSE], zFitMean[1:nSpc[ii], , drop = FALSE],
      zFitResid[1:nSpc[ii]], knots, timeRange)
    timeDepLike[, ii, ] = predResult$timeDepLike
    mleFit[[ii]] = predResult$mleFit
    timePred[, ii] = predResult$timePred}
  return(list(timeDepLike = timeDepLike, mleFit = mleFit, timePred = timePred))}


#' Train and test a ZeitZeiger predictor
#'
#' Train and test a ZeitZeiger predictor, calling the necessary functions.
#'
#' @param xTrain Matrix of measurements for training data, observations in rows
#'   and features in columns.
#' @param timeTrain Vector of values of the periodic variable for training
#'   observations, where 0 corresponds to the lowest possible value and 1
#'   corresponds to the highest possible value.
#' @param xTest Matrix of measurements for test data, observations in rows
#' and features in columns.
#' @param nKnots Number of internal knots to use for the periodic smoothing
#'   spline.
#' @param nTime Number of time-points by which to discretize the time-dependent
#'   behavior of each feature. Corresponds to the number of rows in the matri
#'   for which the SPCs will be calculated.
#' @param useSpc Logical indicating whether to use [PMA::SPC()] (default) or
#'   [base::svd()].
#' @param sumabsv L1-constraint on the SPCs, passed to [PMA::SPC()].
#' @param orth Logical indicating whether to require left singular vectors
#'   be orthogonal to each other, passed to [PMA::SPC()].
#' @param nSpc Vector of the number of SPCs to use for prediction. If `NA`
#'   (default), `nSpc` will become `1:K`, where `K` is the number of SPCs in
#'   `spcResult`. Each value in `nSpc` will correspond to one prediction for
#'   each test observation. A value of 2 means that the prediction will be based
#'   on the first 2 SPCs.
#' @param timeRange Vector of values of the periodic variable at which to
#'   calculate likelihood. The time with the highest likelihood is used as the
#'   initial value for the MLE optimizer.
#'
#' @return
#' \item{fitResult}{Output of [zeitzeigerFit()]}
#' \item{spcResult}{Output of [zeitzeigerSpc()]}
#' \item{predResult}{Output of [zeitzeigerPredict()]}
#'
#' @seealso [zeitzeigerFit()], [zeitzeigerSpc()], [zeitzeigerPredict()]
#'
#' @export
zeitzeiger = function(
  xTrain, timeTrain, xTest, nKnots = 3, nTime = 10, useSpc = TRUE, sumabsv = 2,
  orth = TRUE, nSpc = 2, timeRange = seq(0, 1 - 0.01, 0.01)) {
  fitResult = zeitzeigerFit(xTrain, timeTrain, nKnots)
  spcResult = zeitzeigerSpc(
    fitResult$xFitMean, fitResult$xFitResid, nTime, useSpc, sumabsv, orth)
  predResult = zeitzeigerPredict(
    xTrain, timeTrain, xTest, spcResult, nKnots, nSpc, timeRange)
  result = list(
    fitResult = fitResult, spcResult = spcResult, predResult = predResult)
  return(result)}


#' Train and test a ZeitZeiger predictor, accounting for batch effects
#'
#' Train and test a predictor on multiple datasets independently, using
#' [sva::ComBat()] to correct for batch effects prior to running [zeitzeiger()].
#' This function requires the `metapredict` package.
#'
#' @param ematList Named list of matrices of measurements, one for each dataset,
#'   some of which will be for training, others for testing. Each matrix should
#'   have rownames corresponding to sample names and colnames corresponding to
#'   feature names.
#' @param trainStudyNames Character vector of names in `ematList` corresponding
#'   to datasets for training.
#' @param sampleMetadata data.frame containing relevant information for each
#'   sample across all datasets. Must have a column named `sample`.
#' @param studyColname Name of column in `sampleMetdata` that contains
#'   information about which dataset each sample belongs to.
#' @param batchColname Name of column in `sampleMetdata` that contains
#'   information about which dataset each sample belongs to. This should
#'   correspond to the names of `ematList`, and will often be the same as
#'   `studyColname`, but doesn't have to be.
#' @param timeColname Name of column in `sampleMetdata` that contains the values
#'   of the periodic variable.
#' @param nKnots Number of internal knots to use for the periodic smoothing
#'   spline.
#' @param nTime Number of time-points by which to discretize the time-dependent
#'   behavior of each feature. Corresponds to the number of rows in the matrix
#'   for which the SPCs will be calculated.
#' @param useSpc Logical indicating whether to use [PMA::SPC()] (default) or
#'   [base::svd()].
#' @param sumabsv L1-constraint on the SPCs, passed to [PMA::SPC()].
#' @param orth Logical indicating whether to require left singular vectors
#'   be orthogonal to each other, passed to [PMA::SPC()].
#' @param nSpc Vector of the number of SPCs to use for prediction. If `NA`
#'   (default), `nSpc` will become `1:K`, where `K` is the number of SPCs in
#'   `spcResult`. Each value in `nSpc` will correspond to one prediction for
#'   each test observation. A value of 2 means that the prediction will be based
#'   on the first 2 SPCs.
#' @param timeRange Vector of values of the periodic variable at which to
#'   calculate likelihood. The time with the highest likelihood is used as the
#'   initial value for the MLE optimizer.
#' @param covariateName Name of column(s) in `sampleMetadata` containing
#'   information about other covariates for [sva::ComBat()], besides
#'   `batchColname`. If `NA` (default), then there are no other covariates.
#' @param featuresExclude Named list of character vectors corresponding to
#'   features to exclude from being used for prediction for the respective test
#'   datasets.
#' @param dopar Logical indicating whether to process the folds in parallel. Use
#'   [doParallel::registerDoParallel()] to register the parallel backend.

#' @return
#' \item{spcResultList}{List of output from [zeitzeigerSpc()], one for each test
#'   dataset.}
#' \item{timeDepLike}{3-D array of likelihood, with dimensions for each test
#'   observation (across all datasets), each element of `nSpc`, and each element
#'   of `timeRange`.}
#' \item{mleFit}{List (for each element in `nSpc`) of lists (for each test
#'   observation) of `mle2` objects.}
#' \item{timePred}{Matrix of predicted times for test observations by values of
#'   `nSpc`.}
#'
#' @seealso [zeitzeiger()], [metapredict::metapredict()], [sva::ComBat()]
#'
#' @export
zeitzeigerBatch = function(
  ematList, trainStudyNames, sampleMetadata, studyColname, batchColname,
  timeColname, nKnots = 3, nTime = 10, useSpc = TRUE, sumabsv = 2, orth = TRUE,
  nSpc = 2, timeRange = seq(0, 1 - 0.01, 0.01), covariateName = NA,
  featuresExclude = NULL, dopar = TRUE) {

  testStudyName = NULL
  if (!requireNamespace('metapredict', quietly = TRUE)) {
    stop(paste('This function requires the metapredict package.',
               'Please see https://github.com/hugheylab/metapredict.'),
         call. = FALSE)}

  doOp = ifelse(dopar, `%dopar%`, `%do%`)
  testStudyNames = names(ematList)[!(names(ematList) %in% trainStudyNames)]

  batchResult = doOp(foreach(testStudyName = testStudyNames), {
    ematListNow = ematList[c(trainStudyNames, testStudyName)]
    if (!is.null(featuresExclude) && !is.na(featuresExclude[[testStudyName]])) {
      f = function(emat) emat[!(rownames(emat) %in% featuresExclude[[testStudyName]]),]
      ematListNow = lapply(ematListNow, f)}

    ematMerged = metapredict::mergeStudyData(ematListNow, sampleMetadata, batchColname, covariateName)

    sampleMetadataTrain = sampleMetadata[sampleMetadata[[studyColname]] %in% trainStudyNames, , drop = FALSE]
    xTrain = t(ematMerged[, sampleMetadataTrain$sample])
    xTest = t(ematMerged[, !(colnames(ematMerged) %in% sampleMetadataTrain$sample)])
    timeTrain = sampleMetadataTrain[[timeColname]]

    fitResult = zeitzeigerFit(xTrain, timeTrain, nKnots)
    spcResult = zeitzeigerSpc(
      fitResult$xFitMean, fitResult$xFitResid, nTime, useSpc, sumabsv, orth)
    predResult = zeitzeigerPredict(
      xTrain, timeTrain, xTest, spcResult, nKnots, nSpc, timeRange)

    dimnames(predResult$timeDepLike)[[1]] = rownames(xTest)
    for (ii in 1:length(nSpc)) {
      names(predResult$mleFit[[ii]]) = rownames(xTest)}
    rownames(predResult$timePred) = rownames(xTest)

    list(spcResult = spcResult, predResult = predResult)})

  spcResultList = lapply(batchResult, function(x) x$spcResult)
  names(spcResultList) = testStudyNames

  timeDepLikeList = lapply(batchResult, function(x) x$predResult$timeDepLike)
  timeDepLike = do.call(abind::abind, list(timeDepLikeList, along = 1))

  # list (by batch) of lists (by nSpc) of lists (by obs)
  mleFitList = lapply(batchResult, function(x) x$predResult$mleFit)
  mleFit = list() # list (by nSpc) of lists (by obs)
  for (ii in 1:length(nSpc)) {
    mleFit[[ii]] = do.call(c, lapply(mleFitList, function(x) x[[ii]]))}

  timePred = do.call(rbind, lapply(batchResult, function(x) x$predResult$timePred))

  return(list(spcResultList = spcResultList, timeDepLike = timeDepLike,
              mleFit = mleFit, timePred = timePred))}


#' Combine predictions into an ensemble using the log-likelihood
#'
#' Make predictions by finding the maximum of the sum of the log-likelihoods.
#'
#' @param timeDepLike List or 3-D array of time-dependent likelihood from
#'   [zeitzeigerPredict()]. If a list, then each element (for each member of the
#'   ensemble) should be a matrix in which rows correspond to observations and
#'   columns correspond to time-points. If a 3-D array, the three dimensions
#'   should correspond to observations, time-points, and members of the
#'   ensemble.
#' @param timeRange Vector of time-points at which the likelihood was
#'   calculated.
#'
#' @return
#' \item{timeDepLike}{Matrix of likelihood for observations by time-points.}
#' \item{timePred}{Vector of predicted times. Each predicted time will be an
#'   element of timeRange.}
#'
#' @seealso [zeitzeigerPredict()], [zeitzeigerEnsembleMean()]
#'
#' @export
zeitzeigerEnsembleLikelihood = function(timeDepLike, timeRange) {
  if (is.list(timeDepLike)) {
    if (!all(sapply(timeDepLike, is.matrix)) |
        !all(apply(sapply(timeDepLike, dim), 1, function(r) all(r == r[1])))) {
      stop('If timeDepLike is a list, each element must be a matrix of the same size.')}
    timeDepLike = abind::abind(timeDepLike, along=3)

  } else if (!is.array(timeDepLike) | length(dim(timeDepLike)) != 3){
    stop('timeDepLike must be a list or a 3-D array.')}

  loglike = apply(log(timeDepLike), 1:2, sum)
  timePred = timeRange[apply(loglike, 1, which.max)]
  return(list(timeDepLike = exp(loglike), timePred = timePred))}


#' Combine predictions into an ensemble using the circular mean
#'
#' Make predictions by calculating the circular mean of the predictions across
#'   members of the ensemble.
#'
#' @param timePredInput Matrix of predicted times in which rows correspond to
#'   observations and columns correspond to members of the ensemble.
#' @param timeMax Maximum value of the periodic variable, i.e., the value that
#'   is equivalent to zero.
#' @param na.rm Logical indicating whether `NA` values should be removed from
#'   the calculation.
#'
#' @return Matrix with a row for each observation and columns for the predicted
#'   time and the normalized magnitude of the circular mean. The latter can
#'   range from 0 to 1, with 1 indicating perfect agreement among members of the
#'   ensemble.
#'
#' @seealso [zeitzeigerPredict()], [zeitzeigerEnsembleLikelihood()]
#'
#' @export
zeitzeigerEnsembleMean = function(timePredInput, timeMax = 1, na.rm = TRUE) {
  x = cos(timePredInput / timeMax * 2 * pi)
  y = sin(timePredInput / timeMax * 2 * pi)

  xMean = rowMeans(x, na.rm = na.rm)
  yMean = rowMeans(y, na.rm = na.rm)

  timePred = atan2(yMean, xMean) / 2 / pi
  timePred = ifelse(timePred < 0, timePred + 1, timePred) * timeMax
  magnitude = sqrt(xMean^2 + yMean^2)
  return(cbind(timePred, magnitude))}
