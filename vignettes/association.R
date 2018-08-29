## ---- echo = FALSE-------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = '#>')

## ---- message=FALSE------------------------------------------------------
library('doParallel')
library('dplyr')
library('ggplot2')
library('zeitzeiger')

## ------------------------------------------------------------------------
registerDoParallel(cores=2)
nObs = 50
nTime = 10
nIter = 10
params = expand.grid(snr = c(0, 0.5, 1, 2, 4), phase = c(0, 0.25, 0.5, 0.75))

paramsRep = params[rep(1:nrow(params), each=nIter),]
rownames(paramsRep) = NULL

time = rep_len(seq(0, 1 - 1/nTime, length.out=nTime), length.out=nObs)
set.seed(37203)
x = foreach(ii=1:nrow(params), .combine=cbind) %do% {
	if (params[ii, 'snr']==0) {
		xNow = matrix(rnorm(nObs * nIter), nrow=nObs)
	} else {
		xNoise = matrix(rnorm(nObs * nIter, sd = 2 / params[ii, 'snr']), nrow=nObs)
		xNow = apply(xNoise, 2, function(x) x + cos(2*pi * (time - params[ii, 'phase'])))}}

## ------------------------------------------------------------------------
fitMeanArgs = list(rparm=NA, nknots=3)
fitResult = zeitzeigerFit(x, time, fitMeanArgs, dopar=TRUE)

## ---- fig.width=5, fig.height=3------------------------------------------
timeRange = seq(0, 1, 0.01)
df = data.frame(time = rep(timeRange, times=2),
			condition = rep(c('True signal', 'Spline fit'), each=length(timeRange))) %>%
	mutate(condition = factor(condition, levels=c('True signal', 'Spline fit')))

idx = 16
df[['expr']] = c(cos(2 * pi * (timeRange - paramsRep[idx, 'phase'])),
				predict(fitResult$xFitMean[[idx]], newdata=timeRange))

ggplot(data.frame(time, expr=x[,idx])) +
	geom_line(aes(x=time, y=expr, linetype=condition), data=df) +
	geom_point(aes(x=time, y=expr), shape=1) +
	scale_x_continuous(limits=c(0, 1)) +
	labs(x='Time', y='Expression', title='SNR ≈ 0.5') +
	theme_bw() + theme(legend.title=element_blank())

## ---- fig.width=5, fig.height=3------------------------------------------
idx = 78
df[['expr']] = c(cos(2 * pi * (timeRange - paramsRep[idx, 'phase'])),
				predict(fitResult$xFitMean[[idx]], newdata=timeRange))

ggplot(data.frame(time, expr=x[,idx])) +
	geom_line(aes(x=time, y=expr, linetype=condition), data=df) +
	geom_point(aes(x=time, y=expr), shape=1) +
	scale_x_continuous(limits=c(0, 1)) +
	labs(x='Time', y='Expression', title='SNR ≈ 1') +
	theme_bw() + theme(legend.title=element_blank())

## ---- fig.width=5, fig.height=3------------------------------------------
idx = 135
df[['expr']] = c(cos(2 * pi * (timeRange - paramsRep[idx, 'phase'])),
				predict(fitResult$xFitMean[[idx]], newdata=timeRange))

ggplot(data.frame(time, expr=x[,idx])) +
	geom_line(aes(x=time, y=expr, linetype=condition), data=df) +
	geom_point(aes(x=time, y=expr), shape=1) +
	scale_x_continuous(limits=c(0, 1)) +
	labs(x='Time', y='Expression', title='SNR ≈ 2') +
	theme_bw() + theme(legend.title=element_blank())

## ---- fig.width=5, fig.height=3------------------------------------------
idx = 197
df[['expr']] = c(cos(2 * pi * (timeRange - paramsRep[idx, 'phase'])),
				predict(fitResult$xFitMean[[idx]], newdata=timeRange))

ggplot(data.frame(time, expr=x[,idx])) +
	geom_line(aes(x=time, y=expr, linetype=condition), data=df) +
	geom_point(aes(x=time, y=expr), shape=1) +
	scale_x_continuous(limits=c(0, 1)) +
	labs(x='Time', y='Expression', title='SNR ≈ 4') +
	theme_bw() + theme(legend.title=element_blank())

## ------------------------------------------------------------------------
snrEst = zeitzeigerSnr(fitResult)

## ------------------------------------------------------------------------
phaseEst = zeitzeigerExtrema(fitResult)[,'location']

## ------------------------------------------------------------------------
sig = zeitzeigerSig(x, time, fitMeanArgs, nIter=50)

## ------------------------------------------------------------------------
result = data.frame(paramsRep, snrEst, phaseEst, sig) %>%
	mutate(snrFactor = factor(snr),
			snrChar = paste('SNR =', snr),
			phaseFactor = factor(phase))

## ---- fig.width=4, fig.height=3------------------------------------------
ggplot(result) +
	geom_boxplot(aes(x=snrFactor, y=snrEst), outlier.shape=1) +
	labs(x='Expected SNR', y='Estimated SNR') +
	theme_bw()

## ---- fig.width=4, fig.height=4------------------------------------------
ggplot(result %>% filter(snr > 0)) +
	facet_wrap(~ snrChar, ncol=2) +
	geom_jitter(aes(x=phaseFactor, y=phaseEst), size=2, shape=1, height=0, width=0.2) +
	labs(x='Expected phase', y='Estimated phase') +
	theme_bw()

## ---- fig.width=4, fig.height=3------------------------------------------
ggplot(result) +
	geom_boxplot(aes(x=snrFactor, y=sig), outlier.shape=1) +
	labs(x='Expected SNR', y='Nominal p-value') +
	theme_bw()

