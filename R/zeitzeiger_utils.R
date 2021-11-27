mergeStudyData = function(
  ematList, sampleMetadata, batchColname = 'study', covariateName = NA,
  batchCorrection = TRUE, parPrior = TRUE) {

  studyName = NULL
  sampleNames = do.call(c, lapply(ematList, colnames))
  if (!all(sampleNames %in% sampleMetadata$sample)) {
    stop('sampleMetadata must have samples corresponding to the colnames of each matrix in ematList.',
         call. = FALSE)}

  ematList = ematList[sapply(ematList, ncol) > 0]
  geneIds = Reduce(intersect, lapply(ematList, function(x) rownames(x)))
  ematList2 = foreach(studyName = names(ematList)) %do% {
    ematNow = ematList[[studyName]][geneIds, ]}

  if (isTRUE(batchCorrection)) {
    # if both one-color and two-color data is present and data is not scaled beforehand,
    # ComBat can fail catastrophically
    ematListScaled = lapply(ematList2, function(emat) (emat - mean(emat)) / stats::sd(emat))
    ematMerged = do.call(cbind, ematListScaled)
    # sm = mergeDataTable(colnames(ematMerged), sampleMetadata)
    sm = merge(data.frame(sample = colnames(ematMerged), stringsAsFactors = FALSE),
               sampleMetadata, by = 'sample', sort = FALSE)

    if (is.na(covariateName) || length(unique(sm[[covariateName]])) < 2) {
      covariateInfo = stats::model.matrix(~ rep_len(1, ncol(ematMerged)))
    } else {
      covariateInfo = stats::model.matrix(~ sm[[covariateName]])}

    if (length(unique(sm[[batchColname]])) > 1) {
      ematMergedNorm = sva::ComBat(
        ematMerged, batch = as.character(sm[[batchColname]]),
        mod = covariateInfo, par.prior = parPrior)
    } else {
      ematMergedNorm = ematMerged}

    return(ematMergedNorm)
  } else {
    return(do.call(cbind, ematList2))}}
