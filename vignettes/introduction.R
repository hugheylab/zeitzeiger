## ---- echo = FALSE-------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = '#>')

## ---- message=FALSE------------------------------------------------------
library('tidyverse')
library('doParallel')
library('zeitzeiger')

## ------------------------------------------------------------------------
nObs = 100
nFeatures = 200
amp = 0.5

time = rep((0:9)/10, length.out=nObs)
set.seed(42)
x = matrix(rnorm(nObs*nFeatures, sd=0.7), nrow=nObs)

## ------------------------------------------------------------------------
baseSignal = sin(time*2*pi)
x[,1] = x[,1] + 3 * baseSignal
x[,2] = x[,2] - 0.5 * baseSignal
for (ii in 3:10) {
	x[,ii] = x[,ii] + runif(1, -amp, amp) * baseSignal}

baseSignal = cos(time*2*pi)
x[,11] = x[,11] + 2 * baseSignal
x[,12] = x[,12] - baseSignal
for (ii in 13:20) {
	x[,ii] = x[,ii] + runif(1, -amp, amp) * baseSignal}

baseSignal = cos(time*4*pi + pi/6)
x[,21] = x[,21] + baseSignal
x[,22] = x[,22] - baseSignal
for (ii in 23:30) {
	x[,ii] = x[,ii] + runif(1, -amp, amp) * baseSignal}

## ------------------------------------------------------------------------
idxTrain = 1:round(nObs*0.6)

xTrain = x[idxTrain,]
timeTrain = time[idxTrain]

xTest = x[-idxTrain,]
timeTest = time[-idxTrain]

## ------------------------------------------------------------------------
fitResult = zeitzeigerFit(xTrain, timeTrain)
spcResult = zeitzeigerSpc(fitResult$xFitMean, fitResult$xFitResid, sumabsv=1)
predResult = zeitzeigerPredict(xTrain, timeTrain, xTest, spcResult, nSpc=2)

## ------------------------------------------------------------------------
zzFit = zeitzeiger(xTrain, timeTrain, xTest, sumabsv=1, nSpc=2)

## ------------------------------------------------------------------------
timeRange = seq(0, 1, 0.01)
jj = 11

df1 = data.frame(timeTrain, xTrain=xTrain[,jj])
df2 = data.frame(timeRange,
				  xPred = predict(fitResult$xFitMean[[jj]], newdata=timeRange),
				  xObs = 2*cos(timeRange*2*pi)) %>%
	gather(key=key, value=value, -timeRange) %>%
	mutate(key = factor(key, levels=c('xObs', 'xPred'),
							labels=c('Observed mean', 'Predicted mean')))

## ---- fig.width=4, fig.height=3------------------------------------------
ggplot() +
	geom_point(aes(x=timeTrain, y=xTrain), data=df1, shape=1, size=2) +
	geom_line(aes(x=timeRange, y=value, linetype=key), data=df2) +
	labs(x='Time', y=sprintf('Feature %d', jj)) +
	theme_bw() + theme(legend.position=c(0.5, 0.85), legend.title=element_blank())

## ------------------------------------------------------------------------
dfTest = data.frame(timeObs = timeTest,
					 timePred = predResult$timePred,
					 timeError = calcTimeDiff(timeTest, predResult$timePred))

## ---- fig.width=3.5, fig.height=3----------------------------------------
ggplot(dfTest) +
	geom_point(aes(x=timeObs, y=timePred), size=2, shape=1) +
	geom_abline(slope=1, intercept=0, linetype='dashed') +
	scale_x_continuous(limits=c(0, 1)) + scale_y_continuous(limits=c(0, 1)) +
	labs(x='Observed time', y='Predicted time') + theme_bw()

## ---- fig.width=3.5, fig.height=3----------------------------------------
ggplot(dfTest) +
	geom_point(aes(x=timeObs, y=timeError), size=2, shape=1) +
	scale_x_continuous(limits=c(0, 1)) + scale_y_continuous(limits=c(-0.2, 0.2)) +
	labs(x='Observed time', y='Error') + theme_bw()

## ------------------------------------------------------------------------
obs = 2
dfLike = data.frame(timeRange, likelihood = predResult$timeDepLike[obs, 1,])

dfVert = data.frame(type = c('Observed time', 'Predicted time'),
					 xint = c(timeTest[obs], predResult$timePred[obs]))

## ---- fig.width=4, fig.height=3------------------------------------------
ggplot(dfLike) +
	geom_point(aes(x=timeRange, y=likelihood), size=2, shape=1) +
	geom_vline(aes(xintercept=xint, linetype=type), data=dfVert, show.legend=TRUE) +
	scale_linetype_manual(values=c('solid', 'dashed')) +
	labs(x='Time', y='Likelihood') + theme_bw() +
	theme(legend.position=c(0.7, 0.8), legend.title=element_blank())

## ------------------------------------------------------------------------
registerDoParallel(cores=2)

sumabsv = c(1, 1.5, 3)
nSpc = 1:4
nFolds = 10
foldid = sample(rep(1:nFolds, length.out=nObs))

fitResultList = zeitzeigerFitCv(x, time, foldid)

spcResultList = list()
for (ii in 1:length(sumabsv)) {
	spcResultList[[ii]] = zeitzeigerSpcCv(fitResultList, sumabsv=sumabsv[ii])}

predResultList = list()
for (ii in 1:length(sumabsv)) {
	predResultList[[ii]] = zeitzeigerPredictCv(x, time, foldid, spcResultList[[ii]], nSpc=nSpc)}

## ------------------------------------------------------------------------
timePredList = lapply(predResultList, function(a) a$timePred)

cvResult = data.frame(do.call(rbind, timePredList),
					   timeObs = rep(time, length(sumabsv)),
					   sumabsv = rep(sumabsv, each=length(time)),
					   obs = rep(1:nObs, length(sumabsv)),
					   stringsAsFactors=FALSE)

cvResultGath = gather(cvResult, key=nSpc, value=timePred, -obs, -timeObs, -sumabsv)
cvResultGath$nSpc = as.integer(sapply(as.character(cvResultGath$nSpc),
										   function(a) substr(a, 2, nchar(a))))
cvResultGath$sumabsv = factor(cvResultGath$sumabsv)
cvResultGath$timeError = calcTimeDiff(cvResultGath$timeObs, cvResultGath$timePred)

## ------------------------------------------------------------------------
cvResultGathGroup = cvResultGath %>%
	group_by(sumabsv, nSpc) %>%
	summarize(medae = median(abs(timeError)))

## ---- fig.width=3, fig.height=3------------------------------------------
ggplot(cvResultGathGroup) +
	geom_point(aes(x=nSpc, y=medae, shape=sumabsv, color=sumabsv), size=2) +
	labs(x='Number of SPCs', y='Median absolute error') +
	theme_bw() + theme(legend.position=c(0.7, 0.7))

## ---- fig.width=4, fig.height=3------------------------------------------
fitResultFinal = zeitzeigerFit(x, time)
spcResultFinal = zeitzeigerSpc(fitResultFinal$xFitMean, fitResultFinal$xFitResid, sumabsv=1.5)

dfVar = data.frame(spc = 1:length(spcResultFinal$d),
				    propVar = spcResultFinal$d^2 / sum(spcResultFinal$d^2))

## ------------------------------------------------------------------------
ggplot(dfVar) +
	geom_point(aes(x=spc, y=propVar), size=2, shape=1) +
	scale_x_continuous(breaks=seq(1, 10)) +
	labs(x='SPC', y='Proportion of\nvariance explained') + theme_bw()

## ---- fig.width=5, fig.height=5------------------------------------------
z = x %*% spcResultFinal$v[,1:3]
colnames(z) = c('SPC 1', 'SPC 2', 'SPC 3')
zGath = gather(data.frame(z, obs=1:nObs, Time=time, check.names=FALSE),
				key=SPC, value=Abundance, -obs, -Time)

## ---- fig.width=5, fig.height=5------------------------------------------
ggplot(zGath) +
	facet_grid(SPC ~ ., scales='free_y') +
	geom_point(aes(x=Time, y=Abundance), size=2, shape=1) + theme_bw()

## ------------------------------------------------------------------------
v = data.frame(spcResultFinal$v[,1:3])
colnames(v) = c('SPC 1', 'SPC 2', 'SPC 3')
v = v[apply(v, 1, function(r) any(r!=0)),]
v[v==0] = NA
v = v[do.call(order, v),]
v$feature = rownames(v)
vGath = gather(v, key=spc, value=Coefficient, -feature) %>%
	mutate(feature = factor(feature, levels=rev(v$feature)))

## ---- fig.width=5, fig.height=5------------------------------------------
ggplot(vGath) + facet_wrap(~ spc, nrow=1) +
	geom_bar(aes(x=feature, y=Coefficient), stat='identity') +
	labs(x='Feature') + coord_flip() +
	theme_bw() + theme(panel.spacing=unit(1.2, 'lines'))

