cores = 5

# -----------------------------------------------------------------------------------
# ---- Packages ---------------------------------------------------------------------
# -----------------------------------------------------------------------------------

{

# these are the packages needed to run this script file

# plotting
require(ggplot2)
require(grid)
require(gridExtra)

# modeling
require(forecast)
require(nnet)
require(randomForest)
require(e1071)

# parallel computing
require(foreach)
require(doParallel)

# install them if you haven't installed them before

# install.packages("ggplot2")
# install.packages("grid")
# install.packages("gridExtra")
# install.packages("forecast")
# install.packages("nnet")
# install.packages("foreach")
# install.packages("doParallel")

}

# -----------------------------------------------------------------------------------
# ---- Functions --------------------------------------------------------------------
# -----------------------------------------------------------------------------------

{

# create a sequence of numbers on a log scale

log.seq = function(from, to, length.out)
{
	return(exp(seq(log(from), log(to), length.out = length.out)))
}

# ---- prints the data types of each column in a data frame ----------------------------------------------------------------------------------------------------

types = function(dat)
{
	dat = data.frame(dat)
  
	Column = sapply(1:ncol(dat), function(i)
	(
		colnames(dat)[i]
    ))
  
	Data_Type = sapply(1:ncol(dat), function(i)
    (
		class(dat[,i])
    ))
  
	results = data.frame(cbind(Column, Data_Type))	
	results
}

# ---- plots 6 residual plots ---------------------------------------------------------------------------------------------------------------------------------

residplots = function(actual, fit, binwidth = NULL, from = NULL, to = NULL, by = NULL, histlabel.y = -10, n = NULL, basefont = 20)
{
	require(ggplot2)
	
	residual = actual - fit 
	DF = data.frame("actual" = actual, "fit" = fit, "residual" = residual)
	
    rvfPlot = ggplot(DF, aes(x = fit, y = residual)) + 
			  geom_point(na.rm = TRUE) +
			  stat_smooth(method = "loess", se = FALSE, na.rm = TRUE, color = "blue") +
			  geom_hline(yintercept = 0, col = "red", linetype = "dashed") +
			  xlab("Fitted values") +
			  ylab("Residuals") +
			  ggtitle("Residual vs Fitted Plot") + 
			  theme_bw(base_size = basefont) +
			  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
    
	ggqq = function(x, distribution = "norm", ..., conf = 0.95, probs = c(0.25, 0.75), note = TRUE, alpha = 0.33, main = "", xlab = "\nTheoretical Quantiles", ylab = "Empirical Quantiles\n")
	{
		# compute the sample quantiles and theoretical quantiles
		q.function = eval(parse(text = paste0("q", distribution)))
  		d.function = eval(parse(text = paste0("d", distribution)))
  		x = na.omit(x)
  		ord = order(x)
  		n = length(x)
  		P = ppoints(length(x))
  		DF = data.frame(ord.x = x[ord], z = q.function(P, ...))
  
  		# compute the quantile line
  		Q.x = quantile(DF$ord.x, c(probs[1], probs[2]))
  		Q.z = q.function(c(probs[1], probs[2]), ...)
  		b = diff(Q.x) / diff(Q.z)
  		coef = c(Q.x[1] - (b * Q.z[1]), b)
  
  		# compute the confidence interval band
  		zz = qnorm(1 - (1 - conf) / 2)
  		SE = (coef[2] / d.function(DF$z, ...)) * sqrt(P * (1 - P) / n)
  		fit.value = coef[1] + (coef[2] * DF$z)
  		DF$upper = fit.value + (zz * SE)
  		DF$lower = fit.value - (zz * SE)
  
  		# plot the qqplot
  		p = ggplot(DF, aes(x = z, y = ord.x)) + 
    		geom_point(color = "black", alpha = alpha) +
    		geom_abline(intercept = coef[1], slope = coef[2], size = 1, color = "blue") +
    		geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15) +
    		coord_cartesian(ylim = c(min(DF$ord.x), max(DF$ord.x))) + 
    		labs(x = xlab, y = ylab) +
    		theme_bw(base_size = basefont) +
			theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
						
  		# conditional additions
  		if(main != "")(p = p + ggtitle(main))
  		
  		return(p)
	}

    qqPlot = ggqq(residual, 
				  alpha = 1,				  
				  main = "Normal Q-Q Plot", 
				  xlab = "Theoretical Quantiles", 
				  ylab = "Residuals")
    
    rvtPlot = ggplot(data.frame("x" = 1:length(DF$residual), "y" = DF$residual), aes(x = x, y = y)) + 
			  geom_line(na.rm = TRUE) +
			  geom_hline(yintercept = 0, col = "red", linetype = "dashed") +
			  xlab("Obs. Number") +
			  ylab("Residuals") +
			  ggtitle("Residual Time Series") + 
			  theme_bw(base_size = basefont) +
			  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
        
	variogramDF = function(x)
	{
		n = length(x) - 2
	
		num = sapply(1:n, function(k)
						  sapply(1:(length(x) - k), function(i)
													x[i + k] - x[i]))
	
		num = sapply(1:length(num), function(j)
									var(num[[j]]))

		den = var(sapply(1:(length(x) - 1), function(i)
											x[i + 1] - x[i]))
	
		val = num / den
	
		DF = data.frame("Lag" = 1:n, "Variogram" = val)
	
		return(DF)
	}

	DFv = variogramDF(x = DF$residual)

	varioPlot = ggplot(DFv, aes(x = Lag, y = Variogram)) + 
				geom_point() +
				geom_line(color = "blue") +
				xlab("Lag") +
				ylab("Variogram") +
				ggtitle("Variogram of Residuals") + 
				theme_bw(base_size = basefont) +
				theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
	
	test = t.test(DF$residual)
	
	CI = data.frame("x" = test$estimate, 
					"LCB" = test$conf.int[1], 
					"UCB" = test$conf.int[2], 
					row.names = 1)
	
	histPlot = ggplot(DF, aes(x = residual)) +
			   geom_histogram(color = "white", fill = "black", binwidth = binwidth) +
			   geom_segment(data = CI, aes(x = LCB, xend = LCB, y = 0, yend = Inf), color = "blue") +
			   geom_segment(data = CI, aes(x = UCB, xend = UCB, y = 0, yend = Inf), color = "blue") +
			   annotate("text", x = CI$x, y = histlabel.y, 
						label = "T-Test C.I.", size = 5, 
						color = "blue", fontface = 2) + 
			    ggtitle("Residual Histogram") +
			   labs(x = "Residuals", y = "Frequency") +
	 		   theme_bw(base_size = basefont) +
			   theme(legend.key.size = unit(.25, "in"),
					 legend.position = "bottom", plot.title = element_text(hjust = 0.5))
	
	if(class(from) != "NULL" & class(to) != "NULL" & class(by) != "NULL") (histPlot = histPlot + scale_x_continuous(breaks = seq(from = from, to = to, by = by)))
	
	ggacf = function(x, n = NULL, conf.level = 0.95, main = "ACF Plot", xlab = "Lag", ylab = "Autocorrelation", basefont = 20) 
	{
		if(class(n) == "NULL")
		{
			n = length(x) - 2
		}
	
		ciline = qnorm((1 - conf.level) / 2) / sqrt(length(x))
		bacf = acf(x, lag.max = n, plot = FALSE)
		bacfdf = with(bacf, data.frame(lag, acf))
		bacfdf = bacfdf[-1,]
		
		p = ggplot(bacfdf, aes(x = lag, y = acf)) + 
			geom_bar(stat = "identity", position = "dodge", fill = "black") +
			geom_hline(yintercept = -ciline, color = "blue", size = 1) +
			geom_hline(yintercept = ciline, color = "blue", size = 1) +
			geom_hline(yintercept = 0, color = "red", size = 1) +
			labs(x = xlab, y = ylab) +
			ggtitle(main) +
			theme_bw(base_size = basefont) +
			theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

		return(p)
	}

	acfPlot = ggacf(x = DF$residual, main = "ACF Plot of Residuals", basefont = basefont, n = n)
	
    return(list("rvfPlot" = rvfPlot, 
				"qqPlot" = qqPlot, 
				"rvtPlot" = rvtPlot, 
				"varioPlot" = varioPlot, 
				"histPlot" = histPlot, 
				"acfPlot" = acfPlot))
}

# ---- plots the fitted v. actuals of a predictive model -----------------------------------------------------------------------------------------------------------------------------------------------

modplot = function(mod, xaxis = NULL, y.levels = NULL, actuals = NULL, n.ahead = NULL, level = c(80, 95), newdata = NULL, limits = NULL, nnet.y.min = 0, nnet.y.max = 1, xlab = "Observation", ylab = "Value", main = "Fitted v. Actual\nPrediction Plot", basefont = 20)
{
	require(ggplot2)
	require(forecast)
	# require(nnetpredint)
	# require(dplyr)
	
	if(class(actuals) == "NULL" & "svm" %in% class(mod) & !("residuals" %in% names(mod)))
	{
		p = "Input a vector for 'actuals' to plot an SVM Classification model"
		
	} else
	{
		
		# build fits
		
		if("HoltWinters" %in% class(mod))
		{
			fits = as.numeric(fitted(mod)[,1])
			fits = c(rep(NA, length(mod$x) - length(fits)), fits)
			
		} else if("nnet" %in% class(mod))
		{
			fits = (as.numeric(fitted(mod)) * (nnet.y.max - nnet.y.min)) + nnet.y.min
			
		} else if("randomForest" %in% class(mod))
		{
			fits = as.numeric(mod$predicted)
			
		} else
		{
			fits = as.numeric(fitted(mod))
		}
		
		# build actuals
		
		if(class(actuals) == "NULL")
		{
			if("HoltWinters" %in% class(mod))
			{
				actuals = as.numeric(mod$x)
				
			} else if("nnet" %in% class(mod))
			{
				actuals = ((as.numeric(resid(mod)) + as.numeric(fitted(mod))) * (nnet.y.max - nnet.y.min)) + nnet.y.min
				
			} else if("randomForest" %in% class(mod)) 
			{
				actuals = as.numeric(mod$y)
				
			} else
			{
				actuals = as.numeric(resid(mod)) + fits
			}
			
			if(class(n.ahead) != "NULL")
			{
				actuals = c(actuals, rep(NA, n.ahead))
			}

		} else
		{
			actuals = as.numeric(actuals)
		}
		
		# build xaxis
		
		if(class(xaxis) == "NULL")
		{
			if(class(n.ahead) == "NULL")
			{
				xaxis = 1:length(fits)
				
			} else
			{
				xaxis = 1:(length(fits) + n.ahead)
			}	
		}
		
		# initialize upper and lower prediction interval values
		
		if(class(mod)[1] == "lm")
		{
			upper1 = as.numeric(predict(mod, interval = "confidence", level = 0.95)[,3])
			upper2 = as.numeric(predict(mod, interval = "confidence", level = 0.80)[,3])
			lower1 = as.numeric(predict(mod, interval = "confidence", level = 0.95)[,2])
			lower2 = as.numeric(predict(mod, interval = "confidence", level = 0.80)[,2])
			
		} else
		{
			upper1 = as.numeric(rep(NA, length(fits)))
			upper2 = upper1
			lower1 = upper1
			lower2 = upper1
		}
		
		# initialize a color vector to seperate training data from testing data
		
		color = factor(rep(0, length(fits)), levels = c(0, 1))
		
		# compute predictions if n.ahead is specified 
			# also update fits, upper and lower prediction interval values, and color 
		
		if(class(n.ahead) != "NULL")
		{
			if(class(mod)[1] == "lm")
			{
				predictions = data.frame(forecast(mod, newdata = newdata, h = n.ahead, level = level))
				
			} else if(class(mod)[1] == "glm")
			{
				predictions = data.frame(predict(mod, newdata = newdata, n.ahead = n.ahead, se.fit = TRUE))
				
				predictions = data.frame("fit" = predictions$fit, 
										 "lower2" = predictions$fit - predictions$se.fit,
										 "upper2" = predictions$fit + predictions$se.fit,
										 "lower1" = predictions$fit - (3 * predictions$se.fit),
										 "upper1" = predictions$fit + (3 * predictions$se.fit))
				
			} else if(class(mod)[1] == "ARIMA")
			{
				predictions = data.frame(forecast(mod, h = n.ahead, level = level, xreg = newdata))
				
			} else if(class(mod)[1] == "nnet")
			{				
				predictions = data.frame("fit" = predict(mod, newdata = newdata, n.ahead = n.ahead), 
										 "lower2" = rep(NA, n.ahead),
										 "upper2" = rep(NA, n.ahead),
										 "lower1" = rep(NA, n.ahead),
										 "upper1" = rep(NA, n.ahead))
										 
				predictions$fit = (predictions$fit * (nnet.y.max - nnet.y.min)) + nnet.y.min
				predictions$lower2 = (predictions$lower2 * (nnet.y.max - nnet.y.min)) + nnet.y.min
				predictions$upper2 = (predictions$upper2 * (nnet.y.max - nnet.y.min)) + nnet.y.min
				predictions$lower1 = (predictions$lower1 * (nnet.y.max - nnet.y.min)) + nnet.y.min
				predictions$upper1 = (predictions$upper1 * (nnet.y.max - nnet.y.min)) + nnet.y.min
				
			} else if("randomForest" %in% class(mod))
			{
				predictions = predict(mod, newdata = newdata, predict.all = TRUE)
				
				se = function(x)
				{
					result = sd(x) / sqrt(length(x))
					return(result)
				}
				
				se.fit = as.numeric(apply(X = predictions$individual, MARGIN = 1, FUN = se))
				fit = as.numeric(predictions$aggregate)
				
				predictions = data.frame("fit" = fit, 
										 "lower2" = fit - se.fit,
										 "upper2" = fit + se.fit,
										 "lower1" = fit - (3 * se.fit),
										 "upper1" = fit + (3 * se.fit))
				
			} else if("svm" %in% class(mod))
			{
				predictions = data.frame("fit" = predict(mod, newdata = newdata, n.ahead = n.ahead), 
										 "lower2" = rep(NA, n.ahead),
										 "upper2" = rep(NA, n.ahead),
										 "lower1" = rep(NA, n.ahead),
										 "upper1" = rep(NA, n.ahead))
				
			} else
			{
				predictions = data.frame(forecast(mod, h = n.ahead, level = level))
			}
			
			fits = as.numeric(c(fits, predictions[,1]))
			
			upper1 = as.numeric(c(upper1, predictions[,5]))
			lower1 = as.numeric(c(lower1, predictions[,4]))
			upper2 = as.numeric(c(upper2, predictions[,3]))
			lower2 = as.numeric(c(lower2, predictions[,2]))
			
			color = c(as.numeric(color) * 0, rep(1, length(predictions[,1])))
			color = factor(color, levels = c(0, 1))
		} 

		# build DF
		
		DF = data.frame("Observation" = xaxis, 
						"Actuals" = actuals, 
						"Fitted" = fits,
						"Upper1" = upper1,
						"Lower1" = lower1,
						"Upper2" = upper2,
						"Lower2" = lower2,
						"Color" = color)
			
		# subset DF if limits was specified
		
		if(class(limits) != "NULL")
		{
			DF = DF[-c(which(DF[,1] < limits[1]), which(DF[,1] > limits[2])),]
		}
		
		# plot DF
		
		if(class(y.levels) != "NULL")
		{
			if("glm" %in% class(mod))
			{
				Offset = 1
			} else
			{
				Offset = 0
			}
			
			p = ggplot(DF, aes(x = Observation, y = Actuals, color = Color)) + 
				scale_color_manual(values = c("black", "red"), drop = TRUE, limits = levels(DF$Color)) +
				geom_point(na.rm = TRUE) + 
				geom_line(aes(x = Observation, y = Fitted), color = "blue", na.rm = TRUE) + 
				geom_ribbon(aes(ymin = Lower1, ymax = Upper1), alpha = .1, color = NA) +
				geom_ribbon(aes(ymin = Lower2, ymax = Upper2), alpha = .125, color = NA) +
				scale_y_continuous(breaks = (1:length(y.levels) - Offset), labels = y.levels, limits = (c(0.5, length(y.levels) + 0.5) - Offset)) +
				labs(x = xlab, y = ylab) +
				ggtitle(main) + 
				theme_bw(base_size = basefont) +
				theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
				
		} else
		{
			p = ggplot(DF, aes(x = Observation, y = Actuals, color = Color)) + 
				scale_color_manual(values = c("black", "red"), drop = TRUE, limits = levels(DF$Color)) +
				geom_point(na.rm = TRUE) + 
				geom_line(aes(x = Observation, y = Fitted), color = "blue", na.rm = TRUE) + 
				geom_ribbon(aes(ymin = Lower1, ymax = Upper1), alpha = .1, color = NA) +
				geom_ribbon(aes(ymin = Lower2, ymax = Upper2), alpha = .125, color = NA) +
				labs(x = xlab, y = ylab) +
				ggtitle(main) + 
				theme_bw(base_size = basefont) +
				theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

		}
	}
	
	return(p)
}

}

# -----------------------------------------------------------------------------------
# ---- Import the Data --------------------------------------------------------------
# -----------------------------------------------------------------------------------

{

# make sure the working directory is where crime_weekly.csv is located

getwd()
setwd("F:/Documents/Working/Forecasting Methods/Independent Study/PHI Crime")

# load the data

crime = read.csv("crime_weekly.csv")
head(crime)
types(crime)

# the YEAR.WEEK column is read in as a numeric data type
# lets convert this to a factor that has levels ordered with respect to time

# these nested ifelse statements converts each entry in YEAR.WEEK to a character data type without losing the values for WEEK
# because the YEAR.WEEK column is a numeric data type, the WEEK portion is treated as decimal values
	# week values such as "00" and "10" for week 0 and week 10 will be converted into "" and "1", so we paste "00" and "0" to get "00" and "10"

YEAR.WEEK = ifelse(nchar(crime$YEAR.WEEK) == 4, 
						paste0(crime$YEAR.WEEK, ".00"), 
							ifelse(nchar(crime$YEAR.WEEK) == 6, 
										paste0(crime$YEAR.WEEK, "0"), 
											crime$YEAR.WEEK))

# convert the above character vector into a factor data type to replace the YEAR.WEEK column in crime 
	# the YEAR.WEEK column is already ordered with respect to time so we can use the unique values in the vector to specify the order of the levels
	# unique() returns the first observed values, removing any duplicates that follow

crime$YEAR.WEEK = factor(YEAR.WEEK, levels = unique(YEAR.WEEK))

rm(YEAR.WEEK)

}

# -----------------------------------------------------------------------------------
# ---- Linear Models ----------------------------------------------------------------
# -----------------------------------------------------------------------------------

{

# ---- Central Philly --------------------------------------------------------------

{

# extract the CP subset from crime

DF = crime[which(crime$Area == "CP"),]

# build sine and cosine variables to capture the 54 week seasonality

DF$SIN = sin((2 * pi * (1:nrow(DF))) / 54)
DF$COS = cos((2 * pi * (1:nrow(DF))) / 54)

# ---- fitted model ----

# create the linear model

mod = lm(CRIMES ~ Feb + Mar + Apr + May + Jun + Jul + Aug + Sep + Oct + Nov + Dec + DAYS + UNEMPLOYMENT + UTILITY.GAS + SIN + COS, data = DF)

# plot the model 

plot_lm = modplot(mod = mod,
					xlab = "Observation Number", 
					ylab = "Crimes", 
					main = "Crimes Per Week: Central Philly\nMultiple Linear Regression: Fitted Model")

plot_lm

# plot the residuals

plot_resid = residplots(actual = DF$CRIMES, 
						fit = as.numeric(fitted(mod)), 
						histlabel.y = -3, 
						basefont = 20)

grid.arrange(plot_resid[[1]], 
			 plot_resid[[2]], 
			 plot_resid[[3]], 
			 plot_resid[[4]], 
			 plot_resid[[5]], 
			 plot_resid[[6]], 
			 ncol = 3)		  


# build a table of stats

ttest = data.frame("t.pval" = t.test(DF$CRIMES - as.numeric(fitted(mod)))$p.value)

stats = data.frame(accuracy(f = as.numeric(fitted(mod)), 
							x = DF$CRIMES), 
					row.names = 1)

table_lm = cbind(ttest, stats)
cnames = colnames(table_lm)

table_lm[,2:ncol(table_lm)] = sapply(2:ncol(table_lm), function(i) table_lm[,i] = as.numeric(table_lm[,i]))
colnames(table_lm) = cnames

table_lm

# lets save these results

LM_CP_fit_plot = plot_lm
LM_CP_fit_mod = mod

LM_CP_fit_table = table_lm
rownames(LM_CP_fit_table) = "LM"

rm(table_lm, stats, ttest, plot_resid, plot_lm, mod)

# ---- forecasting model ----

# create the linear model

mod = lm(CRIMES ~ Feb + Mar + Apr + May + Jun + Jul + Aug + Sep + Oct + Nov + Dec + DAYS + UNEMPLOYMENT + UTILITY.GAS + SIN + COS, data = DF[1:478,])

# plot the model 

plot_lm = modplot(mod = mod,
					newdata = DF[479:565, c(4:17, 19:20)],
					n.ahead = 87,
					actuals = DF$CRIMES,
					limits = c(479, 565),
					xlab = "Observation Number", 
					ylab = "Crimes", 
					main = "Crimes Per Week: Central Philly\nMultiple Linear Regression: Forecasting Model")

plot_lm

# plot the residuals

plot_resid = residplots(actual = DF$CRIMES[479:565], 
						fit = as.numeric(predict(mod, newdata = DF[479:565, c(4:17, 19:20)])), 
						histlabel.y = -0.5, 
						basefont = 20)

grid.arrange(plot_resid[[1]], 
			 plot_resid[[2]], 
			 plot_resid[[3]], 
			 plot_resid[[4]], 
			 plot_resid[[5]], 
			 plot_resid[[6]], 
			 ncol = 3)		  

# build a table of stats

ttest = data.frame("t.pval" = t.test(DF$CRIMES[479:565] - as.numeric(predict(mod, newdata = DF[479:565, c(4:17, 19:20)])))$p.value)

stats = data.frame(accuracy(f = as.numeric(predict(mod, newdata = DF[479:565, c(4:17, 19:20)])), 
							x = DF$CRIMES[479:565]), 
					row.names = 1)

table_lm = cbind(ttest, stats)
cnames = colnames(table_lm)

table_lm[,2:ncol(table_lm)] = sapply(2:ncol(table_lm), function(i) table_lm[,i] = as.numeric(table_lm[,i]))
colnames(table_lm) = cnames

table_lm

# lets save these results

LM_CP_forecast_plot = plot_lm
LM_CP_forecast_mod = mod

LM_CP_forecast_table = table_lm
rownames(LM_CP_forecast_table) = "LM"

rm(table_lm, stats, ttest, plot_resid, plot_lm, mod)

}

# ---- West Philly --------------------------------------------------------------

{

# extract the WP subset from crime

DF = crime[which(crime$Area == "WP"),]

# build sine and cosine variables to capture the 54 week seasonality

DF$SIN = sin((2 * pi * (1:nrow(DF))) / 54)
DF$COS = cos((2 * pi * (1:nrow(DF))) / 54)

# ---- fitted model ----

# create the linear model

mod = lm(CRIMES ~ Feb + Mar + Apr + May + Jun + Jul + Aug + Sep + Oct + Nov + Dec + DAYS + UNEMPLOYMENT + UTILITY.GAS + SIN + COS, data = DF)

# plot the model 

plot_lm = modplot(mod = mod,
					xlab = "Observation Number", 
					ylab = "Crimes", 
					main = "Crimes Per Week: West Philly\nMultiple Linear Regression: Fitted Model")

plot_lm

# plot the residuals

plot_resid = residplots(actual = DF$CRIMES, 
						fit = as.numeric(fitted(mod)), 
						histlabel.y = -3, 
						basefont = 20)

grid.arrange(plot_resid[[1]], 
			 plot_resid[[2]], 
			 plot_resid[[3]], 
			 plot_resid[[4]], 
			 plot_resid[[5]], 
			 plot_resid[[6]], 
			 ncol = 3)		  


# build a table of stats

ttest = data.frame("t.pval" = t.test(DF$CRIMES - as.numeric(fitted(mod)))$p.value)

stats = data.frame(accuracy(f = as.numeric(fitted(mod)), 
							x = DF$CRIMES), 
					row.names = 1)

table_lm = cbind(ttest, stats)
cnames = colnames(table_lm)

table_lm[,2:ncol(table_lm)] = sapply(2:ncol(table_lm), function(i) table_lm[,i] = as.numeric(table_lm[,i]))
colnames(table_lm) = cnames

table_lm

# lets save these results

LM_WP_fit_plot = plot_lm
LM_WP_fit_mod = mod

LM_WP_fit_table = table_lm
rownames(LM_WP_fit_table) = "LM"

rm(table_lm, stats, ttest, plot_resid, plot_lm, mod)

# ---- forecasting model ----

# create the linear model

mod = lm(CRIMES ~ Feb + Mar + Apr + May + Jun + Jul + Aug + Sep + Oct + Nov + Dec + DAYS + UNEMPLOYMENT + UTILITY.GAS + SIN + COS, data = DF[1:478,])

# plot the model 

plot_lm = modplot(mod = mod,
					newdata = DF[479:565, c(4:17, 19:20)],
					n.ahead = 87,
					actuals = DF$CRIMES,
					limits = c(479, 565),
					xlab = "Observation Number", 
					ylab = "Crimes", 
					main = "Crimes Per Week: West Philly\nMultiple Linear Regression: Forecasting Model")

plot_lm

# plot the residuals

plot_resid = residplots(actual = DF$CRIMES[479:565], 
						fit = as.numeric(predict(mod, newdata = DF[479:565, c(4:17, 19:20)])), 
						histlabel.y = -0.5, 
						basefont = 20)

grid.arrange(plot_resid[[1]], 
			 plot_resid[[2]], 
			 plot_resid[[3]], 
			 plot_resid[[4]], 
			 plot_resid[[5]], 
			 plot_resid[[6]], 
			 ncol = 3)		  

# build a table of stats

ttest = data.frame("t.pval" = t.test(DF$CRIMES[479:565] - as.numeric(predict(mod, newdata = DF[479:565, c(4:17, 19:20)])))$p.value)

stats = data.frame(accuracy(f = as.numeric(predict(mod, newdata = DF[479:565, c(4:17, 19:20)])), 
							x = DF$CRIMES[479:565]), 
					row.names = 1)

table_lm = cbind(ttest, stats)
cnames = colnames(table_lm)

table_lm[,2:ncol(table_lm)] = sapply(2:ncol(table_lm), function(i) table_lm[,i] = as.numeric(table_lm[,i]))
colnames(table_lm) = cnames

table_lm

# lets save these results

LM_WP_forecast_plot = plot_lm
LM_WP_forecast_mod = mod

LM_WP_forecast_table = table_lm
rownames(LM_WP_forecast_table) = "LM"

rm(table_lm, stats, ttest, plot_resid, plot_lm, mod)

}

# ---- North Philly --------------------------------------------------------------

{

# extract the NP subset from crime

DF = crime[which(crime$Area == "NP"),]

# build sine and cosine variables to capture the 54 week seasonality

DF$SIN = sin((2 * pi * (1:nrow(DF))) / 54)
DF$COS = cos((2 * pi * (1:nrow(DF))) / 54)

# ---- fitted model ----

# create the linear model

mod = lm(CRIMES ~ Feb + Mar + Apr + May + Jun + Jul + Aug + Sep + Oct + Nov + Dec + DAYS + UNEMPLOYMENT + UTILITY.GAS + SIN + COS, data = DF)

# plot the model 

plot_lm = modplot(mod = mod,
					xlab = "Observation Number", 
					ylab = "Crimes", 
					main = "Crimes Per Week: North Philly\nMultiple Linear Regression: Fitted Model")

plot_lm

# plot the residuals

plot_resid = residplots(actual = DF$CRIMES, 
						fit = as.numeric(fitted(mod)), 
						histlabel.y = -3, 
						basefont = 20)

grid.arrange(plot_resid[[1]], 
			 plot_resid[[2]], 
			 plot_resid[[3]], 
			 plot_resid[[4]], 
			 plot_resid[[5]], 
			 plot_resid[[6]], 
			 ncol = 3)		  


# build a table of stats

ttest = data.frame("t.pval" = t.test(DF$CRIMES - as.numeric(fitted(mod)))$p.value)

stats = data.frame(accuracy(f = as.numeric(fitted(mod)), 
							x = DF$CRIMES), 
					row.names = 1)

table_lm = cbind(ttest, stats)
cnames = colnames(table_lm)

table_lm[,2:ncol(table_lm)] = sapply(2:ncol(table_lm), function(i) table_lm[,i] = as.numeric(table_lm[,i]))
colnames(table_lm) = cnames

table_lm

# lets save these results

LM_NP_fit_plot = plot_lm
LM_NP_fit_mod = mod

LM_NP_fit_table = table_lm
rownames(LM_NP_fit_table) = "LM"

rm(table_lm, stats, ttest, plot_resid, plot_lm, mod)

# ---- forecasting model ----

# create the linear model

mod = lm(CRIMES ~ Feb + Mar + Apr + May + Jun + Jul + Aug + Sep + Oct + Nov + Dec + DAYS + UNEMPLOYMENT + UTILITY.GAS + SIN + COS, data = DF[1:478,])

# plot the model 

plot_lm = modplot(mod = mod,
					newdata = DF[479:565, c(4:17, 19:20)],
					n.ahead = 87,
					actuals = DF$CRIMES,
					limits = c(479, 565),
					xlab = "Observation Number", 
					ylab = "Crimes", 
					main = "Crimes Per Week: North Philly\nMultiple Linear Regression: Forecasting Model")

plot_lm

# plot the residuals

plot_resid = residplots(actual = DF$CRIMES[479:565], 
						fit = as.numeric(predict(mod, newdata = DF[479:565, c(4:17, 19:20)])), 
						histlabel.y = -0.5, 
						basefont = 20)

grid.arrange(plot_resid[[1]], 
			 plot_resid[[2]], 
			 plot_resid[[3]], 
			 plot_resid[[4]], 
			 plot_resid[[5]], 
			 plot_resid[[6]], 
			 ncol = 3)		  

# build a table of stats

ttest = data.frame("t.pval" = t.test(DF$CRIMES[479:565] - as.numeric(predict(mod, newdata = DF[479:565, c(4:17, 19:20)])))$p.value)

stats = data.frame(accuracy(f = as.numeric(predict(mod, newdata = DF[479:565, c(4:17, 19:20)])), 
							x = DF$CRIMES[479:565]), 
					row.names = 1)

table_lm = cbind(ttest, stats)
cnames = colnames(table_lm)

table_lm[,2:ncol(table_lm)] = sapply(2:ncol(table_lm), function(i) table_lm[,i] = as.numeric(table_lm[,i]))
colnames(table_lm) = cnames

table_lm

# lets save these results

LM_NP_forecast_plot = plot_lm
LM_NP_forecast_mod = mod

LM_NP_forecast_table = table_lm
rownames(LM_NP_forecast_table) = "LM"

rm(table_lm, stats, ttest, plot_resid, plot_lm, mod)

}

}

# -----------------------------------------------------------------------------------
# ---- HoltWinters Models -----------------------------------------------------------
# -----------------------------------------------------------------------------------

{

# ---- Central Philly --------------------------------------------------------------

# ---- fitted model ----

# extract the CP subset from crime

DF = crime[which(crime$Area == "CP"),]

# build the best fitted model that i found

HW_CP_fit_mod = HoltWinters(ts(DF$CRIMES, frequency = 54), 
							alpha = 0.15, 
							beta = 0.01, 
							gamma = 0.58)

# extract the fitted values 
# extract the actual values 
	# we are using the entries 55:565 because the model is using the first season of values to build seasonal coefficients

fit = as.numeric(fitted(HW_CP_fit_mod)[,1])
actual = DF$CRIMES[55:565]

# compute the performance metrics of this model
	# t.test()$p.value to see if the residuals have an expected value of zero
		# we want p.value > 0.05 to fail to reject the null: mu == 0
	# accuracy() output to see model fitness

ttest = t.test(actual - fit)$p.value
stats = accuracy(f = fit, x = actual)

HW_CP_fit_table = data.frame("t.pval" = ttest, stats)
rownames(HW_CP_fit_table) = "HW"

HW_CP_fit_table

# remove objects we no longer need

rm(fit, actual, ttest, stats)

# plot the fits on the actuals

HW_CP_fit_plot = modplot(mod = HW_CP_fit_mod,
							xlab = "Observation Number", 
							ylab = "Crimes", 
							main = "Crimes Per Week: Central Philly\nHolt-Winters Smoothing: Fitted Model")

HW_CP_fit_plot

# ---- forecasting model ----

# build the best forecasting model that i found
	# we are using entries 1:478 because this is the training data

HW_CP_forecast_mod = HoltWinters(ts(DF$CRIMES[1:478], frequency = 54), 
									alpha = 0.11, 
									beta = 0.08, 
									gamma = 0.4)

# extract the fitted values
# extract the actual values 
	# we are using the entries 479:565 because this is the forecasting range
	
fit = as.numeric(predict(HW_CP_forecast_mod, n.ahead = 87))
actual = DF$CRIMES[479:565]

# compute the performance metrics of this model
	# t.test()$p.value to see if the residuals have an expected value of zero
		# we want p.value > 0.05 to fail to reject the null: mu == 0
	# accuracy() output to see model fitness

ttest = t.test(actual - fit)$p.value
stats = accuracy(f = fit, x = actual)

HW_CP_forecast_table = data.frame("t.pval" = ttest, stats)
rownames(HW_CP_forecast_table) = "HW"

HW_CP_forecast_table

# remove objects we no longer need

rm(fit, actual, ttest, stats)

# plot the fits on the actuals

HW_CP_forecast_plot = modplot(mod = HW_CP_forecast_mod,
								n.ahead = 87,
								actuals = DF$CRIMES,
								limits = c(479, 565),
								level = c(80, 95),
								xlab = "Observation Number", 
								ylab = "Crimes", 
								main = "Crimes Per Week: Central Philly\nHolt-Winters Smoothing: Forecasting Model")

HW_CP_forecast_plot

# ---- West Philly --------------------------------------------------------------

# ---- fitted model ----

# extract the WP subset from crime

DF = crime[which(crime$Area == "WP"),]

# build the best fitted model that i found

HW_WP_fit_mod = HoltWinters(ts(DF$CRIMES, frequency = 54), 
							alpha = 0.3, 
							beta = 0, 
							gamma = 0.51)

# extract the fitted values
# extract the actual values 
	# we are using the entries 55:565 because the model is using the first season of values to build seasonal coefficients

fit = as.numeric(fitted(HW_WP_fit_mod)[,1])
actual = DF$CRIMES[55:565]

# compute the performance metrics of this model
	# t.test()$p.value to see if the residuals have an expected value of zero
		# we want p.value > 0.05 to fail to reject the null: mu == 0
	# accuracy() output to see model fitness

ttest = t.test(actual - fit)$p.value
stats = accuracy(f = fit, x = actual)

HW_WP_fit_table = data.frame("t.pval" = ttest, stats)
rownames(HW_WP_fit_table) = "HW"

HW_WP_fit_table

# remove objects we no longer need

rm(fit, actual, ttest, stats)

# plot the fits on the actuals

HW_WP_fit_plot = modplot(mod = HW_WP_fit_mod,
							xlab = "Observation Number", 
							ylab = "Crimes", 
							main = "Crimes Per Week: West Philly\nHolt-Winters Smoothing: Fitted Model")

HW_WP_fit_plot

# ---- forecasting model ----

# build the best forecasting model that i found
	# we are using entries 1:478 because this is the training data
	
HW_WP_forecast_mod = HoltWinters(ts(DF$CRIMES[1:478], frequency = 54), 
									alpha = 0.02, 
									beta = 0.17, 
									gamma = 0.3)

# extract the fitted values
# extract the actual values 
	# we are using the entries 479:565 because this is the forecasting range
	
fit = as.numeric(predict(HW_WP_forecast_mod, n.ahead = 87))
actual = DF$CRIMES[479:565]

# compute the performance metrics of this model
	# t.test()$p.value to see if the residuals have an expected value of zero
		# we want p.value > 0.05 to fail to reject the null: mu == 0
	# accuracy() output to see model fitness

ttest = t.test(actual - fit)$p.value
stats = accuracy(f = fit, x = actual)

HW_WP_forecast_table = data.frame("t.pval" = ttest, stats)
rownames(HW_WP_forecast_table) = "HW"

HW_WP_forecast_table

# remove objects we no longer need

rm(fit, actual, ttest, stats)

# plot the fits on the actuals

HW_WP_forecast_plot = modplot(mod = HW_WP_forecast_mod,
								n.ahead = 87,
								actuals = DF$CRIMES,
								limits = c(479, 565),
								level = c(80, 95),
								xlab = "Observation Number", 
								ylab = "Crimes", 
								main = "Crimes Per Week: West Philly\nHolt-Winters Smoothing: Forecasting Model")

HW_WP_forecast_plot

# ---- North Philly --------------------------------------------------------------

# ---- fitted model ----

# extract the NP subset from crime

DF = crime[which(crime$Area == "NP"),]

# build the best fitted model that i found

HW_NP_fit_mod = HoltWinters(ts(DF$CRIMES, frequency = 54), 
							alpha = 0.12, 
							beta = 0, 
							gamma = 0.6)

# extract the fitted values
# extract the actual values 
	# we are using the entries 55:565 because the model is using the first season of values to build seasonal coefficients

fit = as.numeric(fitted(HW_NP_fit_mod)[,1])
actual = DF$CRIMES[55:565]

# compute the performance metrics of this model
	# t.test()$p.value to see if the residuals have an expected value of zero
		# we want p.value > 0.05 to fail to reject the null: mu == 0
	# accuracy() output to see model fitness

ttest = t.test(actual - fit)$p.value
stats = accuracy(f = fit, x = actual)

HW_NP_fit_table = data.frame("t.pval" = ttest, stats)
rownames(HW_NP_fit_table) = "HW"

HW_NP_fit_table

# remove objects we no longer need

rm(fit, actual, ttest, stats)

# plot the fits on the actuals

HW_NP_fit_plot = modplot(mod = HW_NP_fit_mod,
							xlab = "Observation Number", 
							ylab = "Crimes", 
							main = "Crimes Per Week: North Philly\nHolt-Winters Smoothing: Fitted Model")

HW_NP_fit_plot

# ---- forecasting model ----

# build the best forecasting model that i found
	# we are using entries 1:478 because this is the training data
	
HW_NP_forecast_mod = HoltWinters(ts(DF$CRIMES[1:478], frequency = 54), 
									alpha = 0.13, 
									beta = 0.05, 
									gamma = 0.42)

# extract the fitted values
# extract the actual values 
	# we are using the entries 479:565 because this is the forecasting range
	
fit = as.numeric(predict(HW_NP_forecast_mod, n.ahead = 87))
actual = DF$CRIMES[479:565]

# compute the performance metrics of this model
	# t.test()$p.value to see if the residuals have an expected value of zero
		# we want p.value > 0.05 to fail to reject the null: mu == 0
	# accuracy() output to see model fitness

ttest = t.test(actual - fit)$p.value
stats = accuracy(f = fit, x = actual)

HW_NP_forecast_table = data.frame("t.pval" = ttest, stats)
rownames(HW_NP_forecast_table) = "HW"

HW_NP_forecast_table

# remove objects we no longer need

rm(fit, actual, ttest, stats)

# plot the fits on the actuals

HW_NP_forecast_plot = modplot(mod = HW_NP_forecast_mod,
								n.ahead = 87,
								actuals = DF$CRIMES,
								limits = c(479, 565),
								level = c(80, 95),
								xlab = "Observation Number", 
								ylab = "Crimes", 
								main = "Crimes Per Week: North Philly\nHolt-Winters Smoothing: Forecasting Model")

HW_NP_forecast_plot

}

# -----------------------------------------------------------------------------------
# ---- ARIMA Models -----------------------------------------------------------------
# -----------------------------------------------------------------------------------

{

# ---- Central Philly --------------------------------------------------------------

# ---- fitted model ----

# extract the CP subset from crime

DF = crime[which(crime$Area == "CP"),]

# build the best fitted model that i found

ARIMA_CP_fit_mod = Arima(ts(DF$CRIMES, frequency = 54), 
							order = c(2, 1, 1), 
							seasonal = list(order = c(0, 0, 1), period = 54),
							optim.control = list(maxit = 1000))

# extract the fitted values 
# extract the actual values 
	# we are using the entries 55:565 because when finding this model, it was compared against other models that used the seasonal difference parameter D
	# when the parameter D is used, the first season of values are used to build seasonal coefficients
	# so for sake of consistency with the experiment, we will exclude the first season when evaluating the performance metrics of this model
	
fit = as.numeric(fitted(ARIMA_CP_fit_mod))[55:565]
actual = DF$CRIMES[55:565]

# compute the performance metrics of this model
	# t.test()$p.value to see if the residuals have an expected value of zero
		# we want p.value > 0.05 to fail to reject the null: mu == 0
	# accuracy() output to see model fitness

ttest = t.test(actual - fit)$p.value
stats = accuracy(f = fit, x = actual)

ARIMA_CP_fit_table = data.frame("t.pval" = ttest, stats)
rownames(ARIMA_CP_fit_table) = "ARIMA"

ARIMA_CP_fit_table

# remove objects we no longer need

rm(fit, actual, ttest, stats)

# plot the fits on the actuals

ARIMA_CP_fit_plot = modplot(mod = ARIMA_CP_fit_mod,
							xlab = "Observation Number", 
							ylab = "Crimes", 
							main = "Crimes Per Week: Central Philly\nARIMA: Fitted Model")

ARIMA_CP_fit_plot

# ---- forecasting model ----

# build the best forecasting model that i found
	# we are using entries 1:478 because this is the training data
	
ARIMA_CP_forecast_mod = Arima(ts(DF$CRIMES[1:478], frequency = 54), 
								order = c(2, 1, 1), 
								seasonal = list(order = c(0, 1, 1), period = 54),
								optim.control = list(maxit = 1000))

# extract the fitted values
# extract the actual values 
	# we are using the entries 479:565 because this is the forecasting range
	# we need to specify se.fit = FALSE or else we will get a list of fits and se.fits
	
fit = as.numeric(predict(ARIMA_CP_forecast_mod, n.ahead = 87, se.fit = FALSE))
actual = DF$CRIMES[479:565]

# compute the performance metrics of this model
	# t.test()$p.value to see if the residuals have an expected value of zero
		# we want p.value > 0.05 to fail to reject the null: mu == 0
	# accuracy() output to see model fitness

ttest = t.test(actual - fit)$p.value
stats = accuracy(f = fit, x = actual)

ARIMA_CP_forecast_table = data.frame("t.pval" = ttest, stats)
rownames(ARIMA_CP_forecast_table) = "ARIMA"

ARIMA_CP_forecast_table

# remove objects we no longer need

rm(fit, actual, ttest, stats)

# plot the fits on the actuals

ARIMA_CP_forecast_plot = modplot(mod = ARIMA_CP_forecast_mod,
									n.ahead = 87,
									actuals = DF$CRIMES,
									limits = c(479, 565),
									level = c(80, 95),
									xlab = "Observation Number", 
									ylab = "Crimes", 
									main = "Crimes Per Week: Central Philly\nARIMA: Forecasting Model")

ARIMA_CP_forecast_plot

# ---- West Philly --------------------------------------------------------------

# ---- fitted model ----

# extract the WP subset from crime

DF = crime[which(crime$Area == "WP"),]

# build the best fitted model that i found

ARIMA_WP_fit_mod = Arima(ts(DF$CRIMES, frequency = 54), 
							order = c(1, 0, 1), 
							seasonal = list(order = c(0, 0, 1), period = 54),
							optim.control = list(maxit = 1000))

# extract the fitted values 
# extract the actual values 
	# we are using the entries 55:565 because when finding this model, it was compared against other models that used the seasonal difference parameter D
	# when the parameter D is used, the first season of values are used to build seasonal coefficients
	# so for sake of consistency with the experiment, we will exclude the first season when evaluating the performance metrics of this model
	
fit = as.numeric(fitted(ARIMA_WP_fit_mod))[55:565]
actual = DF$CRIMES[55:565]

# compute the performance metrics of this model
	# t.test()$p.value to see if the residuals have an expected value of zero
		# we want p.value > 0.05 to fail to reject the null: mu == 0
	# accuracy() output to see model fitness

ttest = t.test(actual - fit)$p.value
stats = accuracy(f = fit, x = actual)

ARIMA_WP_fit_table = data.frame("t.pval" = ttest, stats)
rownames(ARIMA_WP_fit_table) = "ARIMA"

ARIMA_WP_fit_table

# remove objects we no longer need

rm(fit, actual, ttest, stats)

# plot the fits on the actuals

ARIMA_WP_fit_plot = modplot(mod = ARIMA_WP_fit_mod,
							xlab = "Observation Number", 
							ylab = "Crimes", 
							main = "Crimes Per Week: West Philly\nARIMA: Fitted Model")

ARIMA_WP_fit_plot

# ---- forecasting model ----

# build the best forecasting model that i found
	# we are using entries 1:478 because this is the training data

ARIMA_WP_forecast_mod = Arima(ts(DF$CRIMES[1:478], frequency = 54), 
								order = c(3, 1, 1), 
								seasonal = list(order = c(0, 1, 1), period = 54),
								optim.control = list(maxit = 1000))

# extract the fitted values
# extract the actual values 
	# we are using the entries 479:565 because this is the forecasting range
	# we need to specify se.fit = FALSE or else we will get a list of fits and se.fits
	
fit = as.numeric(predict(ARIMA_WP_forecast_mod, n.ahead = 87, se.fit = FALSE))
actual = DF$CRIMES[479:565]

# compute the performance metrics of this model
	# t.test()$p.value to see if the residuals have an expected value of zero
		# we want p.value > 0.05 to fail to reject the null: mu == 0
	# accuracy() output to see model fitness

ttest = t.test(actual - fit)$p.value
stats = accuracy(f = fit, x = actual)

ARIMA_WP_forecast_table = data.frame("t.pval" = ttest, stats)
rownames(ARIMA_WP_forecast_table) = "ARIMA"

ARIMA_WP_forecast_table

# remove objects we no longer need

rm(fit, actual, ttest, stats)

# plot the fits on the actuals

ARIMA_WP_forecast_plot = modplot(mod = ARIMA_WP_forecast_mod,
									n.ahead = 87,
									actuals = DF$CRIMES,
									limits = c(479, 565),
									level = c(80, 95),
									xlab = "Observation Number", 
									ylab = "Crimes", 
									main = "Crimes Per Week: West Philly\nARIMA: Forecasting Model")

ARIMA_WP_forecast_plot

# ---- North Philly --------------------------------------------------------------

# ---- fitted model ----

# extract the NP subset from crime

DF = crime[which(crime$Area == "NP"),]

# build the best fitted model that i found

ARIMA_NP_fit_mod = Arima(ts(DF$CRIMES, frequency = 54), 
							order = c(2, 1, 2), 
							seasonal = list(order = c(1, 0, 0), period = 54),
							optim.control = list(maxit = 1000))

# extract the fitted values 
# extract the actual values 
	# we are using the entries 55:565 because when finding this model, it was compared against other models that used the seasonal difference parameter D
	# when the parameter D is used, the first season of values are used to build seasonal coefficients
	# so for sake of consistency with the experiment, we will exclude the first season when evaluating the performance metrics of this model
	
fit = as.numeric(fitted(ARIMA_NP_fit_mod))[55:565]
actual = DF$CRIMES[55:565]

# compute the performance metrics of this model
	# t.test()$p.value to see if the residuals have an expected value of zero
		# we want p.value > 0.05 to fail to reject the null: mu == 0
	# accuracy() output to see model fitness

ttest = t.test(actual - fit)$p.value
stats = accuracy(f = fit, x = actual)

ARIMA_NP_fit_table = data.frame("t.pval" = ttest, stats)
rownames(ARIMA_NP_fit_table) = "ARIMA"

ARIMA_NP_fit_table

# remove objects we no longer need

rm(fit, actual, ttest, stats)

# plot the fits on the actuals

ARIMA_NP_fit_plot = modplot(mod = ARIMA_NP_fit_mod,
							xlab = "Observation Number", 
							ylab = "Crimes", 
							main = "Crimes Per Week: North Philly\nARIMA: Fitted Model")

ARIMA_NP_fit_plot

# ---- forecasting model ----

# build the best forecasting model that i found
	# we are using entries 1:478 because this is the training data

ARIMA_NP_forecast_mod = Arima(ts(DF$CRIMES[1:478], frequency = 54), 
								order = c(2, 1, 1), 
								seasonal = list(order = c(0, 1, 1), period = 54),
								optim.control = list(maxit = 1000))

# extract the fitted values
# extract the actual values 
	# we are using the entries 479:565 because this is the forecasting range
	# we need to specify se.fit = FALSE or else we will get a list of fits and se.fits
	
fit = as.numeric(predict(ARIMA_NP_forecast_mod, n.ahead = 87, se.fit = FALSE))
actual = DF$CRIMES[479:565]

# compute the performance metrics of this model
	# t.test()$p.value to see if the residuals have an expected value of zero
		# we want p.value > 0.05 to fail to reject the null: mu == 0
	# accuracy() output to see model fitness

ttest = t.test(actual - fit)$p.value
stats = accuracy(f = fit, x = actual)

ARIMA_NP_forecast_table = data.frame("t.pval" = ttest, stats)
rownames(ARIMA_NP_forecast_table) = "ARIMA"

ARIMA_NP_forecast_table

# remove objects we no longer need

rm(fit, actual, ttest, stats)

# plot the fits on the actuals

ARIMA_NP_forecast_plot = modplot(mod = ARIMA_NP_forecast_mod,
									n.ahead = 87,
									actuals = DF$CRIMES,
									limits = c(479, 565),
									level = c(80, 95),
									xlab = "Observation Number", 
									ylab = "Crimes", 
									main = "Crimes Per Week: North Philly\nARIMA: Forecasting Model")

ARIMA_NP_forecast_plot

}

# -----------------------------------------------------------------------------------
# ---- Dynamic Regression Models ----------------------------------------------------
# -----------------------------------------------------------------------------------

{

# ---- Central Philly --------------------------------------------------------------

# ---- fitted model ----

# extract the CP subset from crime

DF = crime[which(crime$Area == "CP"),]

# build the best fitted model that i found

DR_CP_fit_mod = Arima(ts(DF$CRIMES, frequency = 54), 
						order = c(0, 1, 1), 
						seasonal = list(order = c(0, 0, 0), period = 54),
						xreg = DF[,c("SUMMER", "DAYS")],
						optim.control = list(maxit = 1000))

# extract the fitted values 
# extract the actual values 
	# we are using the entries 55:565 because when finding this model, it was compared against other models that used the seasonal difference parameter D
	# when the parameter D is used, the first season of values are used to build seasonal coefficients
	# so for sake of consistency with the experiment, we will exclude the first season when evaluating the performance metrics of this model
	
fit = as.numeric(fitted(DR_CP_fit_mod))[55:565]
actual = DF$CRIMES[55:565]

# compute the performance metrics of this model
	# t.test()$p.value to see if the residuals have an expected value of zero
		# we want p.value > 0.05 to fail to reject the null: mu == 0
	# accuracy() output to see model fitness

ttest = t.test(actual - fit)$p.value
stats = accuracy(f = fit, x = actual)

DR_CP_fit_table = data.frame("t.pval" = ttest, stats)
rownames(DR_CP_fit_table) = "DR"

DR_CP_fit_table

# remove objects we no longer need

rm(fit, actual, ttest, stats)

# plot the fits on the actuals

DR_CP_fit_plot = modplot(mod = DR_CP_fit_mod,
							xlab = "Observation Number", 
							ylab = "Crimes", 
							main = "Crimes Per Week: Central Philly\nDynamic Regression: Fitted Model")

DR_CP_fit_plot

# ---- forecasting model ----

# build the best forecasting model that i found
	# we are using entries 1:478 because this is the training data

DR_CP_forecast_mod = Arima(ts(DF$CRIMES[1:478], frequency = 54), 
						order = c(1, 1, 0), 
						seasonal = list(order = c(0, 1, 1), period = 54),
						xreg = DF[1:478, c("SUMMER", "DAYS")],
						optim.control = list(maxit = 1000))

# extract the fitted values 
# extract the actual values 
	# we are using the entries 479:565 because this is the forecasting range
	# we need to specify se.fit = FALSE or else we will get a list of fits and se.fits
	
fit = as.numeric(predict(DR_CP_forecast_mod, n.ahead = 87, se.fit = FALSE, newxreg = DF[479:565, c("SUMMER", "DAYS")]))
actual = DF$CRIMES[479:565]

# compute the performance metrics of this model
	# t.test()$p.value to see if the residuals have an expected value of zero
		# we want p.value > 0.05 to fail to reject the null: mu == 0
	# accuracy() output to see model fitness

ttest = t.test(actual - fit)$p.value
stats = accuracy(f = fit, x = actual)

DR_CP_forecast_table = data.frame("t.pval" = ttest, stats)
rownames(DR_CP_forecast_table) = "DR"

DR_CP_forecast_table

# remove objects we no longer need

rm(fit, actual, ttest, stats)

# plot the fits on the actuals

DR_CP_forecast_plot = modplot(mod = DR_CP_forecast_mod,
								n.ahead = 87,
								actuals = DF$CRIMES,
								limits = c(479, 565),
								level = c(80, 95),
								newdata = DF[479:565, c("SUMMER", "DAYS")],
								xlab = "Observation Number", 
								ylab = "Crimes", 
								main = "Crimes Per Week: Central Philly\nDynamic Regression: Forecasting Model")

DR_CP_forecast_plot

# ---- West Philly --------------------------------------------------------------

# ---- fitted model ----

# extract the WP subset from crime

DF = crime[which(crime$Area == "WP"),]

# build the best fitted model that i found

DR_WP_fit_mod = Arima(ts(DF$CRIMES, frequency = 54), 
						order = c(0, 1, 1), 
						seasonal = list(order = c(0, 1, 1), period = 54),
						xreg = DF[,c("SUMMER", "DAYS")],
						optim.control = list(maxit = 1000))

# extract the fitted values 
# extract the actual values 
	# we are using the entries 55:565 because the model is using the first season of values to build the seasonal difference parameter

fit = as.numeric(fitted(DR_WP_fit_mod))[55:565]
actual = DF$CRIMES[55:565]

# compute the performance metrics of this model
	# t.test()$p.value to see if the residuals have an expected value of zero
		# we want p.value > 0.05 to fail to reject the null: mu == 0
	# accuracy() output to see model fitness

ttest = t.test(actual - fit)$p.value
stats = accuracy(f = fit, x = actual)

DR_WP_fit_table = data.frame("t.pval" = ttest, stats)
rownames(DR_WP_fit_table) = "DR"

DR_WP_fit_table

# remove objects we no longer need

rm(fit, actual, ttest, stats)

# plot the fits on the actuals

DR_WP_fit_plot = modplot(mod = DR_WP_fit_mod,
							xlab = "Observation Number", 
							ylab = "Crimes", 
							main = "Crimes Per Week: West Philly\nDynamic Regression: Fitted Model")

DR_WP_fit_plot

# ---- forecasting model ----

# build the best forecasting model that i found
	# we are using entries 1:478 because this is the training data

DR_WP_forecast_mod = Arima(ts(DF$CRIMES[1:478], frequency = 54), 
						order = c(0, 1, 1), 
						seasonal = list(order = c(0, 1, 1), period = 54),
						xreg = DF[1:478, c("SUMMER", "DAYS")],
						optim.control = list(maxit = 1000))

# extract the fitted values 
# extract the actual values 
	# we are using the entries 479:565 because this is the forecasting range
	# we need to specify se.fit = FALSE or else we will get a list of fits and se.fits
	
fit = as.numeric(predict(DR_WP_forecast_mod, n.ahead = 87, se.fit = FALSE, newxreg = DF[479:565, c("SUMMER", "DAYS")]))
actual = DF$CRIMES[479:565]

# compute the performance metrics of this model
	# t.test()$p.value to see if the residuals have an expected value of zero
		# we want p.value > 0.05 to fail to reject the null: mu == 0
	# accuracy() output to see model fitness

ttest = t.test(actual - fit)$p.value
stats = accuracy(f = fit, x = actual)

DR_WP_forecast_table = data.frame("t.pval" = ttest, stats)
rownames(DR_WP_forecast_table) = "DR"

DR_WP_forecast_table

# remove objects we no longer need

rm(fit, actual, ttest, stats)

# plot the fits on the actuals

DR_WP_forecast_plot = modplot(mod = DR_WP_forecast_mod,
								n.ahead = 87,
								actuals = DF$CRIMES,
								limits = c(479, 565),
								level = c(80, 95),
								newdata = DF[479:565, c("SUMMER", "DAYS")],
								xlab = "Observation Number", 
								ylab = "Crimes", 
								main = "Crimes Per Week: West Philly\nDynamic Regression: Forecasting Model")

DR_WP_forecast_plot

# ---- North Philly --------------------------------------------------------------

# ---- fitted model ----

# extract the NP subset from crime

DF = crime[which(crime$Area == "NP"),]

# build the best fitted model that i found

DR_NP_fit_mod = Arima(ts(DF$CRIMES, frequency = 54), 
						order = c(2, 1, 0), 
						seasonal = list(order = c(0, 1, 1), period = 54),
						xreg = DF[,c("SUMMER", "DAYS")],
						optim.control = list(maxit = 1000))

# extract the fitted values 
# extract the actual values 
	# we are using the entries 55:565 because the model is using the first season of values to build the seasonal difference parameter

fit = as.numeric(fitted(DR_NP_fit_mod))[55:565]
actual = DF$CRIMES[55:565]

# compute the performance metrics of this model
	# t.test()$p.value to see if the residuals have an expected value of zero
		# we want p.value > 0.05 to fail to reject the null: mu == 0
	# accuracy() output to see model fitness

ttest = t.test(actual - fit)$p.value
stats = accuracy(f = fit, x = actual)

DR_NP_fit_table = data.frame("t.pval" = ttest, stats)
rownames(DR_NP_fit_table) = "DR"

DR_NP_fit_table

# remove objects we no longer need

rm(fit, actual, ttest, stats)

# plot the fits on the actuals

DR_NP_fit_plot = modplot(mod = DR_NP_fit_mod,
							xlab = "Observation Number", 
							ylab = "Crimes", 
							main = "Crimes Per Week: North Philly\nDynamic Regression: Fitted Model")

DR_NP_fit_plot

# ---- forecasting model ----

# build the best forecasting model that i found
	# we are using entries 1:478 because this is the training data
	
DR_NP_forecast_mod = Arima(ts(DF$CRIMES[1:478], frequency = 54), 
						order = c(1, 1, 1), 
						seasonal = list(order = c(0, 1, 1), period = 54),
						xreg = DF[1:478, c("SUMMER", "DAYS")],
						optim.control = list(maxit = 1000))

# extract the fitted values 
# extract the actual values 
	# we are using the entries 479:565 because this is the forecasting range
	# we need to specify se.fit = FALSE or else we will get a list of fits and se.fits
	
fit = as.numeric(predict(DR_NP_forecast_mod, n.ahead = 87, se.fit = FALSE, newxreg = DF[479:565, c("SUMMER", "DAYS")]))
actual = DF$CRIMES[479:565]

# compute the performance metrics of this model
	# t.test()$p.value to see if the residuals have an expected value of zero
		# we want p.value > 0.05 to fail to reject the null: mu == 0
	# accuracy() output to see model fitness

ttest = t.test(actual - fit)$p.value
stats = accuracy(f = fit, x = actual)

DR_NP_forecast_table = data.frame("t.pval" = ttest, stats)
rownames(DR_NP_forecast_table) = "DR"

DR_NP_forecast_table

# remove objects we no longer need

rm(fit, actual, ttest, stats)

# plot the fits on the actuals

DR_NP_forecast_plot = modplot(mod = DR_NP_forecast_mod,
								n.ahead = 87,
								actuals = DF$CRIMES,
								limits = c(479, 565),
								level = c(80, 95),
								newdata = DF[479:565, c("SUMMER", "DAYS")],
								xlab = "Observation Number", 
								ylab = "Crimes", 
								main = "Crimes Per Week: North Philly\nDynamic Regression: Forecasting Model")

DR_NP_forecast_plot

}

# -----------------------------------------------------------------------------------
# ---- Neural Network Models --------------------------------------------------------
# -----------------------------------------------------------------------------------

{

# ---- Central Philly --------------------------------------------------------------

# ---- fitted model ----

# extract the CP subset from crime

DF = crime[which(crime$Area == "CP"),]

# keep a copy of the response variable
# normalize all variables that will be used in the model to a 0-1 scale
	# the month columns are already on a 0-1 scale

CRIMES = DF$CRIMES

for(i in 15:18)
{
	DF[,i] = (DF[,i] - min(DF[,i], na.rm = TRUE)) / (max(DF[,i], na.rm = TRUE) - min(DF[,i], na.rm = TRUE))
}

rm(i)

# build the best fitted model that i found

NNET_CP_fit_mod = nnet(y = DF$CRIMES,
						x = DF[,4:17],
						size = 30,
						decay = 0,
						maxit = 2000,
						linout = TRUE,
						trace = FALSE)

# extract the fitted values 
# extract the actual values 
	# we are un-normalizing the fits so we can see the model performance in terms of the original scale

fit = (as.numeric(fitted(NNET_CP_fit_mod)) * (max(CRIMES) - min(CRIMES))) + min(CRIMES)
actual = CRIMES

# compute the performance metrics of this model
	# t.test()$p.value to see if the residuals have an expected value of zero
		# we want p.value > 0.05 to fail to reject the null: mu == 0
	# accuracy() output to see model fitness

ttest = t.test(actual - fit)$p.value
stats = accuracy(f = fit, x = actual)

NNET_CP_fit_table = data.frame("t.pval" = ttest, stats)
rownames(NNET_CP_fit_table) = "NNET"

NNET_CP_fit_table

# remove objects we no longer need

rm(fit, actual, ttest, stats)

# plot the fits on the actuals

NNET_CP_fit_plot = modplot(mod = NNET_CP_fit_mod,
							xlab = "Observation Number", 
							ylab = "Crimes", 
							main = "Crimes Per Week: Central Philly\nNeural Network: Fitted Model")

NNET_CP_fit_plot

# ---- forecasting model ----

# build the best forecasting model that i found
	# we are using entries 1:478 because this is the training data

NNET_CP_forecast_mod = nnet(y = DF$CRIMES[1:478],
							x = DF[1:478, 4:17],
							size = 40,
							decay = 0.05,
							maxit = 2000,
							linout = TRUE,
							trace = FALSE)

# extract the fitted values 
# extract the actual values 
	# we are un-normalizing the fits so we can see the model performance in terms of the original scale
	# we are using the entries 479:565 because this is the forecasting range

fit = (as.numeric(predict(NNET_CP_forecast_mod, newdata = DF[479:565, 4:17], n.ahead = 87)) * (max(CRIMES) - min(CRIMES))) + min(CRIMES)
actual = CRIMES[479:565]

# compute the performance metrics of this model
	# t.test()$p.value to see if the residuals have an expected value of zero
		# we want p.value > 0.05 to fail to reject the null: mu == 0
	# accuracy() output to see model fitness

ttest = t.test(actual - fit)$p.value
stats = accuracy(f = fit, x = actual)

NNET_CP_forecast_table = data.frame("t.pval" = ttest, stats)
rownames(NNET_CP_forecast_table) = "NNET"

NNET_CP_forecast_table

# remove objects we no longer need

rm(fit, actual, ttest, stats)

# plot the fits on the actuals

NNET_CP_forecast_plot = modplot(mod = NNET_CP_forecast_mod,
									n.ahead = 87,
									actuals = CRIMES,
									limits = c(479, 565),
									level = c(80, 95),
									newdata = DF[479:565, 4:17],
									nnet.y.min = min(CRIMES), 
									nnet.y.max = max(CRIMES),
									xlab = "Observation Number", 
									ylab = "Crimes", 
									main = "Crimes Per Week: Central Philly\nNeural Network: Forecasting Model")

NNET_CP_forecast_plot

# ---- West Philly --------------------------------------------------------------

# ---- fitted model ----

# extract the WP subset from crime

DF = crime[which(crime$Area == "WP"),]

# keep a copy of the response variable
# normalize all variables that will be used in the model to a 0-1 scale
	# the month columns are already on a 0-1 scale
	
CRIMES = DF$CRIMES

for(i in 15:18)
{
	DF[,i] = (DF[,i] - min(DF[,i], na.rm = TRUE)) / (max(DF[,i], na.rm = TRUE) - min(DF[,i], na.rm = TRUE))
}

rm(i)

# build the best fitted model that i found

NNET_WP_fit_mod = nnet(y = DF$CRIMES,
						x = DF[,4:17],
						size = 30,
						decay = 0,
						maxit = 2000,
						linout = TRUE,
						trace = FALSE)

# extract the fitted values 
# extract the actual values 
	# we are un-normalizing the fits so we can see the model performance in terms of the original scale

fit = (as.numeric(fitted(NNET_WP_fit_mod)) * (max(CRIMES) - min(CRIMES))) + min(CRIMES)
actual = CRIMES

# compute the performance metrics of this model
	# t.test()$p.value to see if the residuals have an expected value of zero
		# we want p.value > 0.05 to fail to reject the null: mu == 0
	# accuracy() output to see model fitness

ttest = t.test(actual - fit)$p.value
stats = accuracy(f = fit, x = actual)

NNET_WP_fit_table = data.frame("t.pval" = ttest, stats)
rownames(NNET_WP_fit_table) = "NNET"

NNET_WP_fit_table

# remove objects we no longer need

rm(fit, actual, ttest, stats)

# plot the fits on the actuals

NNET_WP_fit_plot = modplot(mod = NNET_WP_fit_mod,
							xlab = "Observation Number", 
							ylab = "Crimes", 
							main = "Crimes Per Week: West Philly\nNeural Network: Fitted Model")

NNET_WP_fit_plot

# ---- forecasting model ----

# build the best forecasting model that i found
	# we are using entries 1:478 because this is the training data

NNET_WP_forecast_mod = nnet(y = DF$CRIMES[1:478],
							x = DF[1:478, 4:17],
							size = 50,
							decay = 0.2,
							maxit = 2000,
							linout = TRUE,
							trace = FALSE)

# extract the fitted values 
# extract the actual values 
	# we are un-normalizing the fits so we can see the model performance in terms of the original scale
	# we are using the entries 479:565 because this is the forecasting range

fit = (as.numeric(predict(NNET_WP_forecast_mod, newdata = DF[479:565, 4:17], n.ahead = 87)) * (max(CRIMES) - min(CRIMES))) + min(CRIMES)
actual = CRIMES[479:565]

# compute the performance metrics of this model
	# t.test()$p.value to see if the residuals have an expected value of zero
		# we want p.value > 0.05 to fail to reject the null: mu == 0
	# accuracy() output to see model fitness

ttest = t.test(actual - fit)$p.value
stats = accuracy(f = fit, x = actual)

NNET_WP_forecast_table = data.frame("t.pval" = ttest, stats)
rownames(NNET_WP_forecast_table) = "NNET"

NNET_WP_forecast_table

# remove objects we no longer need

rm(fit, actual, ttest, stats)

# plot the fits on the actuals

NNET_WP_forecast_plot = modplot(mod = NNET_WP_forecast_mod,
									n.ahead = 87,
									actuals = CRIMES,
									limits = c(479, 565),
									level = c(80, 95),
									newdata = DF[479:565, 4:17],
									nnet.y.min = min(CRIMES), 
									nnet.y.max = max(CRIMES),
									xlab = "Observation Number", 
									ylab = "Crimes", 
									main = "Crimes Per Week: West Philly\nNeural Network: Forecasting Model")

NNET_WP_forecast_plot

# ---- North Philly --------------------------------------------------------------

# ---- fitted model ----

# extract the NP subset from crime

DF = crime[which(crime$Area == "NP"),]

# keep a copy of the response variable
# normalize all variables that will be used in the model to a 0-1 scale
	# the month columns are already on a 0-1 scale
	
CRIMES = DF$CRIMES

for(i in 15:18)
{
	DF[,i] = (DF[,i] - min(DF[,i], na.rm = TRUE)) / (max(DF[,i], na.rm = TRUE) - min(DF[,i], na.rm = TRUE))
}

rm(i)

# build the best fitted model that i found

NNET_NP_fit_mod = nnet(y = DF$CRIMES,
						x = DF[,4:17],
						size = 30,
						decay = 0,
						maxit = 2000,
						linout = TRUE,
						trace = FALSE)

# extract the fitted values 
# extract the actual values 
	# we are un-normalizing the fits so we can see the model performance in terms of the original scale

fit = (as.numeric(fitted(NNET_NP_fit_mod)) * (max(CRIMES) - min(CRIMES))) + min(CRIMES)
actual = CRIMES

# compute the performance metrics of this model
	# t.test()$p.value to see if the residuals have an expected value of zero
		# we want p.value > 0.05 to fail to reject the null: mu == 0
	# accuracy() output to see model fitness

ttest = t.test(actual - fit)$p.value
stats = accuracy(f = fit, x = actual)

NNET_NP_fit_table = data.frame("t.pval" = ttest, stats)
rownames(NNET_NP_fit_table) = "NNET"

NNET_NP_fit_table

# remove objects we no longer need

rm(fit, actual, ttest, stats)

# plot the fits on the actuals

NNET_NP_fit_plot = modplot(mod = NNET_NP_fit_mod,
							xlab = "Observation Number", 
							ylab = "Crimes", 
							main = "Crimes Per Week: North Philly\nNeural Network: Fitted Model")

NNET_NP_fit_plot

# ---- forecasting model ----

# build the best forecasting model that i found
	# we are using entries 1:478 because this is the training data

NNET_NP_forecast_mod = nnet(y = DF$CRIMES[1:478],
							x = DF[1:478, 4:17],
							size = 50,
							decay = 0.05,
							maxit = 2000,
							linout = TRUE,
							trace = FALSE)

# extract the fitted values 
# extract the actual values 
	# we are un-normalizing the fits so we can see the model performance in terms of the original scale
	# we are using the entries 479:565 because this is the forecasting range

fit = (as.numeric(predict(NNET_NP_forecast_mod, newdata = DF[479:565, 4:17], n.ahead = 87)) * (max(CRIMES) - min(CRIMES))) + min(CRIMES)
actual = CRIMES[479:565]

# compute the performance metrics of this model
	# t.test()$p.value to see if the residuals have an expected value of zero
		# we want p.value > 0.05 to fail to reject the null: mu == 0
	# accuracy() output to see model fitness

ttest = t.test(actual - fit)$p.value
stats = accuracy(f = fit, x = actual)

NNET_NP_forecast_table = data.frame("t.pval" = ttest, stats)
rownames(NNET_NP_forecast_table) = "NNET"

NNET_NP_forecast_table

# remove objects we no longer need

rm(fit, actual, ttest, stats)

# plot the fits on the actuals

NNET_NP_forecast_plot = modplot(mod = NNET_NP_forecast_mod,
									n.ahead = 87,
									actuals = CRIMES,
									limits = c(479, 565),
									level = c(80, 95),
									newdata = DF[479:565, 4:17],
									nnet.y.min = min(CRIMES), 
									nnet.y.max = max(CRIMES),
									xlab = "Observation Number", 
									ylab = "Crimes", 
									main = "Crimes Per Week: North Philly\nNeural Network: Forecasting Model")

NNET_NP_forecast_plot

rm(DF, CRIMES)

}

# -----------------------------------------------------------------------------------
# ---- Random Forest Models ---------------------------------------------------------
# -----------------------------------------------------------------------------------

{

# ---- Central Philly --------------------------------------------------------------

{

# ---- fitted model ----

# extract the CP subset from crime

DF = crime[which(crime$Area == "CP"),]

# lets look at the automated model

set.seed(42)

rfor = randomForest(formula = CRIMES ~ Feb + Mar + Apr + May + Jun + Jul + Aug + Sep + Oct + Nov + Dec + DAYS + UNEMPLOYMENT + UTILITY.GAS, data = DF)

# randomForest input parameters of interest:
	# ntree = 500 ~ number of decision trees to create
	# mtry = floor(p/3) (where p is number of regressors) ~ the number of regressors to randomly choose, from which, one is randomly chosen to then determine how to recursively split a parent node into two child nodes
	# nodesize = 5 ~ minimum size of terminal nodes (ie. the minimum number of data points that can be grouped together in any node of a tree)

# lets create a DOE to test a range of values for each of the above randomForest parameters

DOE = expand.grid(ntree = seq(from = 500, to = 1000, by = 100),
				  mtry = seq(from = floor(14 / 4), to = floor(14 / 2), by = 1),
				  nodesize = seq(from = 3, to = 15, by = 1))

DOE = rbind(c(500, 4, 5), DOE)
DOE = DOE[!duplicated(DOE),]
rownames(DOE) = 1:nrow(DOE)

head(DOE)
tail(DOE)
nrow(DOE)

rm(rfor)

# lets compute the residuals for each random forest model in the DOE

registerDoParallel(cores = cores) 

rfors = foreach(i = 1:nrow(DOE)) %dopar%
{
	require(foreach)
	require(randomForest)
	
	set.seed(42)
	
	DFs = data.frame("actual" = DF$CRIMES,
					"fit" = as.numeric(randomForest(formula = CRIMES ~ Feb + Mar + Apr + May + Jun + Jul + Aug + Sep + Oct + Nov + Dec + DAYS + UNEMPLOYMENT + UTILITY.GAS,
													data = DF,
													ntree = DOE$ntree[i],
													nodesize = DOE$nodesize[i])$predicted))
													
	return(DFs)
}

registerDoSEQ()

# values greater than 0.05 are ideal 

ttests = sapply(1:nrow(DOE), function(i) 
							 t.test(rfors[[i]]$actual - rfors[[i]]$fit)$p.value)

# values close to zero are ideal 

stats = data.frame(t(sapply(1:nrow(DOE), function(i)
										 data.frame(accuracy(f = rfors[[i]]$fit, 
															 x = rfors[[i]]$actual), 
													row.names = 1))))

# lets add these test results to the DOE

DOE$t.pval = ttests
DOE = cbind(DOE, stats)
cnames = colnames(DOE)

DOE[,5:9] = sapply(5:9, function(i) DOE[,i] = as.numeric(DOE[,i]))
colnames(DOE) = cnames

head(DOE)

rm(rfors, ttests, stats, cnames)

# lets filter out the top models

# lets plot the histograms of these reponses to find filter limits

par(mfrow = c(3, 2))
lapply(4:9, function(i) hist(DOE[,i], main = colnames(DOE)[i]))

head(DOE)

# apply filters

x = DOE[which(DOE$t.pval >= 0.72 & 
			  DOE$ME >= -1 &
			  DOE$MAPE <= 10 & 
			  DOE$RMSE <= 62),]

par(mfrow = c(3, 2))
lapply(4:9, function(i) hist(x[,i], main = colnames(x)[i]))

x

# lets compare the auto fit with the chosen fit from the DOE

DOE[c(1, 168),]

set.seed(42)

rfor = lapply(c(1, 168), function(i)
			   randomForest(formula = CRIMES ~ Feb + Mar + Apr + May + Jun + Jul + Aug + Sep + Oct + Nov + Dec + DAYS + UNEMPLOYMENT + UTILITY.GAS,
							data = DF,
							ntree = DOE$ntree[i],
							nodesize = DOE$nodesize[i]))

# lets plot the model & residuals for the auto fit and chosen fit 

	# auto fit
	
plot_rfor1 = modplot(mod = rfor[[1]],
					  xlab = "Observation Number", 
					  ylab = "Crimes", 
					  main = "Crimes Per Week: Central Philly\nRandom Forest: Fitted Model")

plot_rfor1

plot_resid1 = residplots(actual = DF$CRIMES, 
						fit = as.numeric(rfor[[1]]$predicted), 
						histlabel.y = -3, 
						basefont = 20)

grid.arrange(plot_resid1[[1]], 
			 plot_resid1[[2]], 
			 plot_resid1[[3]], 
			 plot_resid1[[4]], 
			 plot_resid1[[5]], 
			 plot_resid1[[6]], 
			 ncol = 3)		  

	# chosen fit

plot_rfor2 = modplot(mod = rfor[[2]],
					  xlab = "Observation Number", 
					  ylab = "Crimes", 
					  main = "Crimes Per Week: Central Philly\nRandom Forest: Fitted Model")

plot_rfor2

plot_resid2 = residplots(actual = DF$CRIMES, 
						fit = as.numeric(rfor[[2]]$predicted), 
						histlabel.y = -3, 
						basefont = 20)

grid.arrange(plot_resid2[[1]], 
			 plot_resid2[[2]], 
			 plot_resid2[[3]], 
			 plot_resid2[[4]], 
			 plot_resid2[[5]], 
			 plot_resid2[[6]], 
			 ncol = 3)		  

RF_CP_fit_plot = plot_rfor2
RF_CP_fit_mod = rfor[[2]]

RF_CP_fit_table = DOE[168, 4:9]
rownames(RF_CP_fit_table) = "RF"

rm(x, DOE, rfor, plot_rfor1, plot_resid1, plot_rfor2, plot_resid2)

# ---- forecasting model ----

# build the best forecasting model that i found
	# we are using entries 1:478 because this is the training data

# lets look at the automated model

set.seed(42)

rfor = randomForest(formula = CRIMES ~ Feb + Mar + Apr + May + Jun + Jul + Aug + Sep + Oct + Nov + Dec + DAYS + UNEMPLOYMENT + UTILITY.GAS, data = DF[1:478,])

# randomForest input parameters of interest:
	# ntree = 500 ~ number of trees to create (ie. the number of row-column subsets to create) 
	# nodesize = 5 ~ size of terminal nodes (ie. the minimum number of data points that can be clustered together in any node of a tree)

# lets create a DOE to test a range of values for each of the above randomForest parameters

DOE = expand.grid(ntree = seq(from = 500, to = 1000, by = 100),
				  mtry = seq(from = floor(14 / 4), to = floor(14 / 2), by = 1),
				  nodesize = seq(from = 3, to = 15, by = 1))

DOE = rbind(c(500, 4, 5), DOE)
DOE = DOE[!duplicated(DOE),]
rownames(DOE) = 1:nrow(DOE)

head(DOE)
tail(DOE)
nrow(DOE)

rm(rfor)

# lets compute the residuals for each random forest model in the DOE

registerDoParallel(cores = cores) 

rfors = foreach(i = 1:nrow(DOE)) %dopar%
{
	require(foreach)
	require(randomForest)

	set.seed(42)
	
	DFs = data.frame("actual" = DF$CRIMES[479:565],
					"fit" = as.numeric(predict(randomForest(formula = CRIMES ~ Feb + Mar + Apr + May + Jun + Jul + Aug + Sep + Oct + Nov + Dec + DAYS + UNEMPLOYMENT + UTILITY.GAS,
															data = DF[1:478,],
															ntree = DOE$ntree[i],
															nodesize = DOE$nodesize[i]),
												newdata = DF[479:565, 4:17])))
													
	return(DFs)
}

registerDoSEQ()

# values greater than 0.05 are ideal 

ttests = sapply(1:nrow(DOE), function(i) 
							 t.test(rfors[[i]]$actual - rfors[[i]]$fit)$p.value)

# values close to zero are ideal 

stats = data.frame(t(sapply(1:nrow(DOE), function(i)
										 data.frame(accuracy(f = rfors[[i]]$fit, 
															 x = rfors[[i]]$actual), 
													row.names = 1))))

# lets add these test results to the DOE

DOE$t.pval = ttests
DOE = cbind(DOE, stats)
cnames = colnames(DOE)

DOE[,5:9] = sapply(5:9, function(i) DOE[,i] = as.numeric(DOE[,i]))
colnames(DOE) = cnames

head(DOE)

rm(rfors, ttests, stats, cnames)

# lets filter out the top models

# lets plot the histograms of these reponses to find filter limits

par(mfrow = c(2, 3))
lapply(4:9, function(i) hist(DOE[,i], main = colnames(DOE)[i]))

head(DOE)

# apply filters

x = DOE[which(DOE$t.pval >= 0.3),]

par(mfrow = c(2, 3))
lapply(4:9, function(i) hist(x[,i], main = colnames(x)[i]))

x

# lets compare the auto fit with the chosen fit from the DOE

DOE[c(1, 133),]

set.seed(42)

rfor = lapply(c(1, 133), function(i)
			   randomForest(formula = CRIMES ~ Feb + Mar + Apr + May + Jun + Jul + Aug + Sep + Oct + Nov + Dec + DAYS + UNEMPLOYMENT + UTILITY.GAS,
							data = DF[1:478,],
							ntree = DOE$ntree[i],
							nodesize = DOE$nodesize[i]))

# lets plot the model & residuals for the auto fit and chosen fit 

	# auto fit
	
plot_rfor1 = modplot(mod = rfor[[1]],
					  newdata = DF[479:565,4:17],
					  n.ahead = 87,
					  actuals = DF$CRIMES,
					  limits = c(479, 565),
					  xlab = "Observation Number", 
					  ylab = "Crimes", 
					  main = "Crimes Per Week: Central Philly\nRandom Forest: Forecasting Model")

plot_rfor1

plot_resid1 = residplots(actual = DF$CRIMES[479:565], 
						fit = as.numeric(predict(rfor[[1]], newdata = DF[479:565, 4:17])), 
						histlabel.y = -0.5, 
						basefont = 20)

grid.arrange(plot_resid1[[1]], 
			 plot_resid1[[2]], 
			 plot_resid1[[3]], 
			 plot_resid1[[4]], 
			 plot_resid1[[5]], 
			 plot_resid1[[6]], 
			 ncol = 3)		  

	# chosen fit

plot_rfor2 = modplot(mod = rfor[[2]],
					  newdata = DF[479:565,4:17],
					  n.ahead = 87,
					  actuals = DF$CRIMES,
					  limits = c(479, 565),
					  xlab = "Observation Number", 
					  ylab = "Crimes", 
					  main = "Crimes Per Week: Central Philly\nRandom Forest: Forecasting Model")

plot_rfor2

plot_resid2 = residplots(actual = DF$CRIMES[479:565], 
						fit = as.numeric(predict(rfor[[2]], newdata = DF[479:565, 4:17])), 
						histlabel.y = -0.5, 
						basefont = 20)

grid.arrange(plot_resid2[[1]], 
			 plot_resid2[[2]], 
			 plot_resid2[[3]], 
			 plot_resid2[[4]], 
			 plot_resid2[[5]], 
			 plot_resid2[[6]], 
			 ncol = 3)		  

RF_CP_forecast_plot = plot_rfor2
RF_CP_forecast_mod = rfor[[2]]

RF_CP_forecast_table = DOE[133, 4:9]
rownames(RF_CP_forecast_table) = "RF"

rm(x, DOE, rfor, plot_rfor1, plot_resid1, plot_rfor2, plot_resid2)

}

# ---- West Philly --------------------------------------------------------------

{

# ---- fitted model ----

# extract the WP subset from crime

DF = crime[which(crime$Area == "WP"),]

# lets look at the automated model

set.seed(42)

rfor = randomForest(formula = CRIMES ~ Feb + Mar + Apr + May + Jun + Jul + Aug + Sep + Oct + Nov + Dec + DAYS + UNEMPLOYMENT + UTILITY.GAS, data = DF)

# randomForest input parameters of interest:
	# ntree = 500 ~ number of decision trees to create
	# mtry = floor(p/3) (where p is number of regressors) ~ the number of regressors to randomly choose, from which, one is randomly chosen to then determine how to recursively split a parent node into two child nodes
	# nodesize = 5 ~ minimum size of terminal nodes (ie. the minimum number of data points that can be grouped together in any node of a tree)

# lets create a DOE to test a range of values for each of the above randomForest parameters

DOE = expand.grid(ntree = seq(from = 500, to = 1000, by = 100),
				  mtry = seq(from = floor(14 / 4), to = floor(14 / 2), by = 1),
				  nodesize = seq(from = 3, to = 15, by = 1))

DOE = rbind(c(500, 4, 5), DOE)
DOE = DOE[!duplicated(DOE),]
rownames(DOE) = 1:nrow(DOE)

head(DOE)
tail(DOE)
nrow(DOE)

rm(rfor)

# lets compute the residuals for each random forest model in the DOE

registerDoParallel(cores = cores) 

rfors = foreach(i = 1:nrow(DOE)) %dopar%
{
	require(foreach)
	require(randomForest)
	
	set.seed(42)
	
	DFs = data.frame("actual" = DF$CRIMES,
					"fit" = as.numeric(randomForest(formula = CRIMES ~ Feb + Mar + Apr + May + Jun + Jul + Aug + Sep + Oct + Nov + Dec + DAYS + UNEMPLOYMENT + UTILITY.GAS,
													data = DF,
													ntree = DOE$ntree[i],
													nodesize = DOE$nodesize[i])$predicted))
													
	return(DFs)
}

registerDoSEQ()

# values greater than 0.05 are ideal 

ttests = sapply(1:nrow(DOE), function(i) 
							 t.test(rfors[[i]]$actual - rfors[[i]]$fit)$p.value)

# values close to zero are ideal 

stats = data.frame(t(sapply(1:nrow(DOE), function(i)
										 data.frame(accuracy(f = rfors[[i]]$fit, 
															 x = rfors[[i]]$actual), 
													row.names = 1))))

# lets add these test results to the DOE

DOE$t.pval = ttests
DOE = cbind(DOE, stats)
cnames = colnames(DOE)

DOE[,5:9] = sapply(5:9, function(i) DOE[,i] = as.numeric(DOE[,i]))
colnames(DOE) = cnames

head(DOE)

rm(rfors, ttests, stats, cnames)

# lets filter out the top models

# lets plot the histograms of these reponses to find filter limits

par(mfrow = c(3, 2))
lapply(4:9, function(i) hist(DOE[,i], main = colnames(DOE)[i]))

head(DOE)

# apply filters

x = DOE[which(DOE$t.pval >= 0.95 & 
			  DOE$ME >= -1 &
			  DOE$MAPE <= 10 & 
			  DOE$RMSE <= 62),]

par(mfrow = c(3, 2))
lapply(4:9, function(i) hist(x[,i], main = colnames(x)[i]))

x

# lets compare the auto fit with the chosen fit from the DOE

DOE[c(1, 283),]

set.seed(42)

rfor = lapply(c(1, 283), function(i)
			   randomForest(formula = CRIMES ~ Feb + Mar + Apr + May + Jun + Jul + Aug + Sep + Oct + Nov + Dec + DAYS + UNEMPLOYMENT + UTILITY.GAS,
							data = DF,
							ntree = DOE$ntree[i],
							nodesize = DOE$nodesize[i]))

# lets plot the model & residuals for the auto fit and chosen fit 

	# auto fit
	
plot_rfor1 = modplot(mod = rfor[[1]],
					  xlab = "Observation Number", 
					  ylab = "Crimes", 
					  main = "Crimes Per Week: West Philly\nRandom Forest: Fitted Model")

plot_rfor1

plot_resid1 = residplots(actual = DF$CRIMES, 
						fit = as.numeric(rfor[[1]]$predicted), 
						histlabel.y = -3, 
						basefont = 20)

grid.arrange(plot_resid1[[1]], 
			 plot_resid1[[2]], 
			 plot_resid1[[3]], 
			 plot_resid1[[4]], 
			 plot_resid1[[5]], 
			 plot_resid1[[6]], 
			 ncol = 3)		  

	# chosen fit

plot_rfor2 = modplot(mod = rfor[[2]],
					  xlab = "Observation Number", 
					  ylab = "Crimes", 
					  main = "Crimes Per Week: West Philly\nRandom Forest: Fitted Model")

plot_rfor2

plot_resid2 = residplots(actual = DF$CRIMES, 
						fit = as.numeric(rfor[[2]]$predicted), 
						histlabel.y = -3, 
						basefont = 20)

grid.arrange(plot_resid2[[1]], 
			 plot_resid2[[2]], 
			 plot_resid2[[3]], 
			 plot_resid2[[4]], 
			 plot_resid2[[5]], 
			 plot_resid2[[6]], 
			 ncol = 3)		  

RF_WP_fit_plot = plot_rfor2
RF_WP_fit_mod = rfor[[2]]

RF_WP_fit_table = DOE[283, 4:9]
rownames(RF_WP_fit_table) = "RF"

rm(x, DOE, rfor, plot_rfor1, plot_resid1, plot_rfor2, plot_resid2)

# ---- forecasting model ----

# build the best forecasting model that i found
	# we are using entries 1:478 because this is the training data

# lets look at the automated model

set.seed(42)

rfor = randomForest(formula = CRIMES ~ Feb + Mar + Apr + May + Jun + Jul + Aug + Sep + Oct + Nov + Dec + DAYS + UNEMPLOYMENT + UTILITY.GAS, data = DF[1:478,])

# randomForest input parameters of interest:
	# ntree = 500 ~ number of trees to create (ie. the number of row-column subsets to create) 
	# nodesize = 5 ~ size of terminal nodes (ie. the minimum number of data points that can be clustered together in any node of a tree)

# lets create a DOE to test a range of values for each of the above randomForest parameters

DOE = expand.grid(ntree = seq(from = 500, to = 1000, by = 100),
				  mtry = seq(from = floor(14 / 4), to = floor(14 / 2), by = 1),
				  nodesize = seq(from = 3, to = 15, by = 1))

DOE = rbind(c(500, 4, 5), DOE)
DOE = DOE[!duplicated(DOE),]
rownames(DOE) = 1:nrow(DOE)

head(DOE)
tail(DOE)
nrow(DOE)

rm(rfor)

# lets compute the residuals for each random forest model in the DOE

registerDoParallel(cores = cores) 

rfors = foreach(i = 1:nrow(DOE)) %dopar%
{
	require(foreach)
	require(randomForest)

	set.seed(42)
	
	DFs = data.frame("actual" = DF$CRIMES[479:565],
					"fit" = as.numeric(predict(randomForest(formula = CRIMES ~ Feb + Mar + Apr + May + Jun + Jul + Aug + Sep + Oct + Nov + Dec + DAYS + UNEMPLOYMENT + UTILITY.GAS,
															data = DF[1:478,],
															ntree = DOE$ntree[i],
															nodesize = DOE$nodesize[i]),
												newdata = DF[479:565, 4:17])))
													
	return(DFs)
}

registerDoSEQ()

# values greater than 0.05 are ideal 

ttests = sapply(1:nrow(DOE), function(i) 
							 t.test(rfors[[i]]$actual - rfors[[i]]$fit)$p.value)

# values close to zero are ideal 

stats = data.frame(t(sapply(1:nrow(DOE), function(i)
										 data.frame(accuracy(f = rfors[[i]]$fit, 
															 x = rfors[[i]]$actual), 
													row.names = 1))))

# lets add these test results to the DOE

DOE$t.pval = ttests
DOE = cbind(DOE, stats)
cnames = colnames(DOE)

DOE[,5:9] = sapply(5:9, function(i) DOE[,i] = as.numeric(DOE[,i]))
colnames(DOE) = cnames

head(DOE)

rm(rfors, ttests, stats, cnames)

# lets filter out the top models

# lets plot the histograms of these reponses to find filter limits

par(mfrow = c(2, 3))
lapply(4:9, function(i) hist(DOE[,i], main = colnames(DOE)[i]))

head(DOE)

# apply filters

x = DOE[which(DOE$MAPE <= 16 &
			  DOE$ME >= -40),]

par(mfrow = c(2, 3))
lapply(4:9, function(i) hist(x[,i], main = colnames(x)[i]))

x

# lets compare the auto fit with the chosen fit from the DOE

DOE[c(1, 194),]

set.seed(42)

rfor = lapply(c(1, 194), function(i)
			   randomForest(formula = CRIMES ~ Feb + Mar + Apr + May + Jun + Jul + Aug + Sep + Oct + Nov + Dec + DAYS + UNEMPLOYMENT + UTILITY.GAS,
							data = DF[1:478,],
							ntree = DOE$ntree[i],
							nodesize = DOE$nodesize[i]))

# lets plot the model & residuals for the auto fit and chosen fit 

	# auto fit
	
plot_rfor1 = modplot(mod = rfor[[1]],
					  newdata = DF[479:565,4:17],
					  n.ahead = 87,
					  actuals = DF$CRIMES,
					  limits = c(479, 565),
					  xlab = "Observation Number", 
					  ylab = "Crimes", 
					  main = "Crimes Per Week: West Philly\nRandom Forest: Forecasting Model")

plot_rfor1

plot_resid1 = residplots(actual = DF$CRIMES[479:565], 
						fit = as.numeric(predict(rfor[[1]], newdata = DF[479:565, 4:17])), 
						histlabel.y = -0.5, 
						basefont = 20)

grid.arrange(plot_resid1[[1]], 
			 plot_resid1[[2]], 
			 plot_resid1[[3]], 
			 plot_resid1[[4]], 
			 plot_resid1[[5]], 
			 plot_resid1[[6]], 
			 ncol = 3)		  

	# chosen fit

plot_rfor2 = modplot(mod = rfor[[2]],
					  newdata = DF[479:565,4:17],
					  n.ahead = 87,
					  actuals = DF$CRIMES,
					  limits = c(479, 565),
					  xlab = "Observation Number", 
					  ylab = "Crimes", 
					  main = "Crimes Per Week: West Philly\nRandom Forest: Forecasting Model")

plot_rfor2

plot_resid2 = residplots(actual = DF$CRIMES[479:565], 
						fit = as.numeric(predict(rfor[[2]], newdata = DF[479:565, 4:17])), 
						histlabel.y = -0.5, 
						basefont = 20)

grid.arrange(plot_resid2[[1]], 
			 plot_resid2[[2]], 
			 plot_resid2[[3]], 
			 plot_resid2[[4]], 
			 plot_resid2[[5]], 
			 plot_resid2[[6]], 
			 ncol = 3)		  

RF_WP_forecast_plot = plot_rfor2
RF_WP_forecast_mod = rfor[[2]]

RF_WP_forecast_table = DOE[194, 4:9]
rownames(RF_WP_forecast_table) = "RF"

rm(x, DOE, rfor, plot_rfor1, plot_resid1, plot_rfor2, plot_resid2)

}

# ---- North Philly --------------------------------------------------------------

{

# ---- fitted model ----

# extract the NP subset from crime

DF = crime[which(crime$Area == "NP"),]

# lets look at the automated model

set.seed(42)

rfor = randomForest(formula = CRIMES ~ Feb + Mar + Apr + May + Jun + Jul + Aug + Sep + Oct + Nov + Dec + DAYS + UNEMPLOYMENT + UTILITY.GAS, data = DF)

# randomForest input parameters of interest:
	# ntree = 500 ~ number of decision trees to create
	# mtry = floor(p/3) (where p is number of regressors) ~ the number of regressors to randomly choose, from which, one is randomly chosen to then determine how to recursively split a parent node into two child nodes
	# nodesize = 5 ~ minimum size of terminal nodes (ie. the minimum number of data points that can be grouped together in any node of a tree)

# lets create a DOE to test a range of values for each of the above randomForest parameters

DOE = expand.grid(ntree = seq(from = 500, to = 1000, by = 100),
				  mtry = seq(from = floor(14 / 4), to = floor(14 / 2), by = 1),
				  nodesize = seq(from = 3, to = 15, by = 1))

DOE = rbind(c(500, 4, 5), DOE)
DOE = DOE[!duplicated(DOE),]
rownames(DOE) = 1:nrow(DOE)

head(DOE)
tail(DOE)
nrow(DOE)

rm(rfor)

# lets compute the residuals for each random forest model in the DOE

registerDoParallel(cores = cores) 

rfors = foreach(i = 1:nrow(DOE)) %dopar%
{
	require(foreach)
	require(randomForest)
	
	set.seed(42)
	
	DFs = data.frame("actual" = DF$CRIMES,
					"fit" = as.numeric(randomForest(formula = CRIMES ~ Feb + Mar + Apr + May + Jun + Jul + Aug + Sep + Oct + Nov + Dec + DAYS + UNEMPLOYMENT + UTILITY.GAS,
													data = DF,
													ntree = DOE$ntree[i],
													nodesize = DOE$nodesize[i])$predicted))
													
	return(DFs)
}

registerDoSEQ()

# values greater than 0.05 are ideal 

ttests = sapply(1:nrow(DOE), function(i) 
							 t.test(rfors[[i]]$actual - rfors[[i]]$fit)$p.value)

# values close to zero are ideal 

stats = data.frame(t(sapply(1:nrow(DOE), function(i)
										 data.frame(accuracy(f = rfors[[i]]$fit, 
															 x = rfors[[i]]$actual), 
													row.names = 1))))

# lets add these test results to the DOE

DOE$t.pval = ttests
DOE = cbind(DOE, stats)
cnames = colnames(DOE)

DOE[,5:9] = sapply(5:9, function(i) DOE[,i] = as.numeric(DOE[,i]))
colnames(DOE) = cnames

head(DOE)

rm(rfors, ttests, stats, cnames)

# lets filter out the top models

# lets plot the histograms of these reponses to find filter limits

par(mfrow = c(3, 2))
lapply(4:9, function(i) hist(DOE[,i], main = colnames(DOE)[i]))

head(DOE)

# apply filters

x = DOE[which(DOE$t.pval >= 0.8 & 
			  DOE$ME >= -1 &
			  DOE$MAPE <= 9),]

par(mfrow = c(3, 2))
lapply(4:9, function(i) hist(x[,i], main = colnames(x)[i]))

x

# lets compare the auto fit with the chosen fit from the DOE

DOE[c(1, 313),]

set.seed(42)

rfor = lapply(c(1, 313), function(i)
			   randomForest(formula = CRIMES ~ Feb + Mar + Apr + May + Jun + Jul + Aug + Sep + Oct + Nov + Dec + DAYS + UNEMPLOYMENT + UTILITY.GAS,
							data = DF,
							ntree = DOE$ntree[i],
							nodesize = DOE$nodesize[i]))

# lets plot the model & residuals for the auto fit and chosen fit 

	# auto fit
	
plot_rfor1 = modplot(mod = rfor[[1]],
					  xlab = "Observation Number", 
					  ylab = "Crimes", 
					  main = "Crimes Per Week: North Philly\nRandom Forest: Fitted Model")

plot_rfor1

plot_resid1 = residplots(actual = DF$CRIMES, 
						fit = as.numeric(rfor[[1]]$predicted), 
						histlabel.y = -3, 
						basefont = 20)

grid.arrange(plot_resid1[[1]], 
			 plot_resid1[[2]], 
			 plot_resid1[[3]], 
			 plot_resid1[[4]], 
			 plot_resid1[[5]], 
			 plot_resid1[[6]], 
			 ncol = 3)		  

	# chosen fit

plot_rfor2 = modplot(mod = rfor[[2]],
					  xlab = "Observation Number", 
					  ylab = "Crimes", 
					  main = "Crimes Per Week: North Philly\nRandom Forest: Fitted Model")

plot_rfor2

plot_resid2 = residplots(actual = DF$CRIMES, 
						fit = as.numeric(rfor[[2]]$predicted), 
						histlabel.y = -3, 
						basefont = 20)

grid.arrange(plot_resid2[[1]], 
			 plot_resid2[[2]], 
			 plot_resid2[[3]], 
			 plot_resid2[[4]], 
			 plot_resid2[[5]], 
			 plot_resid2[[6]], 
			 ncol = 3)		  

RF_NP_fit_plot = plot_rfor2
RF_NP_fit_mod = rfor[[2]]

RF_NP_fit_table = DOE[313, 4:9]
rownames(RF_NP_fit_table) = "RF"

rm(x, DOE, rfor, plot_rfor1, plot_resid1, plot_rfor2, plot_resid2)

# ---- forecasting model ----

# build the best forecasting model that i found
	# we are using entries 1:478 because this is the training data

# lets look at the automated model

set.seed(42)

rfor = randomForest(formula = CRIMES ~ Feb + Mar + Apr + May + Jun + Jul + Aug + Sep + Oct + Nov + Dec + DAYS + UNEMPLOYMENT + UTILITY.GAS, data = DF[1:478,])

# randomForest input parameters of interest:
	# ntree = 500 ~ number of trees to create (ie. the number of row-column subsets to create) 
	# nodesize = 5 ~ size of terminal nodes (ie. the minimum number of data points that can be clustered together in any node of a tree)

# lets create a DOE to test a range of values for each of the above randomForest parameters

DOE = expand.grid(ntree = seq(from = 500, to = 1000, by = 100),
				  mtry = seq(from = floor(14 / 4), to = floor(14 / 2), by = 1),
				  nodesize = seq(from = 3, to = 15, by = 1))

DOE = rbind(c(500, 4, 5), DOE)
DOE = DOE[!duplicated(DOE),]
rownames(DOE) = 1:nrow(DOE)

head(DOE)
tail(DOE)
nrow(DOE)

rm(rfor)

# lets compute the residuals for each random forest model in the DOE

registerDoParallel(cores = cores) 

rfors = foreach(i = 1:nrow(DOE)) %dopar%
{
	require(foreach)
	require(randomForest)

	set.seed(42)
	
	DFs = data.frame("actual" = DF$CRIMES[479:565],
					"fit" = as.numeric(predict(randomForest(formula = CRIMES ~ Feb + Mar + Apr + May + Jun + Jul + Aug + Sep + Oct + Nov + Dec + DAYS + UNEMPLOYMENT + UTILITY.GAS,
															data = DF[1:478,],
															ntree = DOE$ntree[i],
															nodesize = DOE$nodesize[i]),
												newdata = DF[479:565, 4:17])))
													
	return(DFs)
}

registerDoSEQ()

# values greater than 0.05 are ideal 

ttests = sapply(1:nrow(DOE), function(i) 
							 t.test(rfors[[i]]$actual - rfors[[i]]$fit)$p.value)

# values close to zero are ideal 

stats = data.frame(t(sapply(1:nrow(DOE), function(i)
										 data.frame(accuracy(f = rfors[[i]]$fit, 
															 x = rfors[[i]]$actual), 
													row.names = 1))))

# lets add these test results to the DOE

DOE$t.pval = ttests
DOE = cbind(DOE, stats)
cnames = colnames(DOE)

DOE[,5:9] = sapply(5:9, function(i) DOE[,i] = as.numeric(DOE[,i]))
colnames(DOE) = cnames

head(DOE)

rm(rfors, ttests, stats, cnames)

# lets filter out the top models

# lets plot the histograms of these reponses to find filter limits

par(mfrow = c(2, 3))
lapply(4:9, function(i) hist(DOE[,i], main = colnames(DOE)[i]))

head(DOE)

# apply filters

x = DOE[which(DOE$t.pval >= 0.2 &
			  DOE$MAPE <= 10 &
			  DOE$ME >= -25),]

par(mfrow = c(2, 3))
lapply(4:9, function(i) hist(x[,i], main = colnames(x)[i]))

x

# lets compare the auto fit with the chosen fit from the DOE

DOE[c(1, 253),]

set.seed(42)

rfor = lapply(c(1, 253), function(i)
			   randomForest(formula = CRIMES ~ Feb + Mar + Apr + May + Jun + Jul + Aug + Sep + Oct + Nov + Dec + DAYS + UNEMPLOYMENT + UTILITY.GAS,
							data = DF[1:478,],
							ntree = DOE$ntree[i],
							nodesize = DOE$nodesize[i]))

# lets plot the model & residuals for the auto fit and chosen fit 

	# auto fit
	
plot_rfor1 = modplot(mod = rfor[[1]],
					  newdata = DF[479:565,4:17],
					  n.ahead = 87,
					  actuals = DF$CRIMES,
					  limits = c(479, 565),
					  xlab = "Observation Number", 
					  ylab = "Crimes", 
					  main = "Crimes Per Week: North Philly\nRandom Forest: Forecasting Model")

plot_rfor1

plot_resid1 = residplots(actual = DF$CRIMES[479:565], 
						fit = as.numeric(predict(rfor[[1]], newdata = DF[479:565, 4:17])), 
						histlabel.y = -0.5, 
						basefont = 20)

grid.arrange(plot_resid1[[1]], 
			 plot_resid1[[2]], 
			 plot_resid1[[3]], 
			 plot_resid1[[4]], 
			 plot_resid1[[5]], 
			 plot_resid1[[6]], 
			 ncol = 3)		  

	# chosen fit

plot_rfor2 = modplot(mod = rfor[[2]],
					  newdata = DF[479:565,4:17],
					  n.ahead = 87,
					  actuals = DF$CRIMES,
					  limits = c(479, 565),
					  xlab = "Observation Number", 
					  ylab = "Crimes", 
					  main = "Crimes Per Week: North Philly\nRandom Forest: Forecasting Model")

plot_rfor2

plot_resid2 = residplots(actual = DF$CRIMES[479:565], 
						fit = as.numeric(predict(rfor[[2]], newdata = DF[479:565, 4:17])), 
						histlabel.y = -0.5, 
						basefont = 20)

grid.arrange(plot_resid2[[1]], 
			 plot_resid2[[2]], 
			 plot_resid2[[3]], 
			 plot_resid2[[4]], 
			 plot_resid2[[5]], 
			 plot_resid2[[6]], 
			 ncol = 3)		  

RF_NP_forecast_plot = plot_rfor2
RF_NP_forecast_mod = rfor[[2]]

RF_NP_forecast_table = DOE[253, 4:9]
rownames(RF_NP_forecast_table) = "RF"

rm(x, DOE, rfor, plot_rfor1, plot_resid1, plot_rfor2, plot_resid2)

}

}

# -----------------------------------------------------------------------------------
# ---- Support Vector Machine Models ------------------------------------------------
# -----------------------------------------------------------------------------------

{

# we will focus on nu-regression, aka v-regression
# this is becuase the parameter 'nu' takes on values [0, 1] so it makes tuning simplier than esp-regression which uses the 'epsilon' parameter which takes on values >= 0

# we will focus on the radial basis as our kernel function for projecting our data into a higher demensional feature space becuase:
	# this seems to be the popular method
	# this only requires one parameter: 'gamma'
	# this may end up being less computationally intensive becuase we are not actually putting our data into higher demsional space, we are estimating it accurately
	# this is the default method of choice in the e1071 package

# important information regarding our parameters:
	# cost:
		# The penalty factor. If it is too large, we have a high penalty for nonseparable points and we may store many support vectors and overfit. If it is too small, we may have underfitting.
		# In practice, a logarithmic grid from 10^(-3) to 10^3 is usually sufficient for tuning 'cost'
	# nu:
		# It is used to determine the proportion of the number of support vectors you desire to keep in your solution with respect to the total number of samples in the dataset.
		# If this takes on a large value, then all of your data points could become support vectors, which is overfitting
	# gamma:
		# If gamma is too large, the radius of the area of influence of the support vectors only includes the support vector itself and no amount of regularization with 'cost' will be able to prevent overfitting.
		# When gamma is very small, the model is too constrained and cannot capture the complexity or “shape” of the data. The region of influence of any selected support vector would include the whole training set.
		# grid.py tool from the libSVM package checks values of gamma from 2^-15 to 2^3 
		
# the default values of our parameters are:
	# cost = 1
	# nu = 0.5
	# gamma = 1 / (number of regressors)

# ---- Central Philly --------------------------------------------------------------

{

# ---- fitted model ----

# extract the CP subset from crime

DF = crime[which(crime$Area == "CP"),]

# lets create a DOE to test a range of values for each of the svm parameters of interest

DOE = expand.grid(cost = log.seq(from = 10^(-3), to = 10^3, length.out = 10),
				  gamma = seq(from = 2^(-15), to = 2^3, length.out = 10),
				  nu = seq(from = 0.3, to = 0.7, by = 0.1))

DOE = rbind(c(1, 1/14, 0.5), DOE)
DOE = DOE[!duplicated(DOE),]
rownames(DOE) = 1:nrow(DOE)

head(DOE)
tail(DOE)
nrow(DOE)

# lets compute the residuals for each svm model in the DOE

registerDoParallel(cores = cores) 

svms = foreach(i = 1:nrow(DOE)) %dopar%
{
	require(foreach)
	require(e1071)
	
	DFs = data.frame("actual" = DF$CRIMES,
					"fit" = as.numeric(fitted(svm(formula = CRIMES ~ Feb + Mar + Apr + May + Jun + Jul + Aug + Sep + Oct + Nov + Dec + DAYS + UNEMPLOYMENT + UTILITY.GAS,
													data = DF,
													type = "nu-regression",
													cost = DOE$cost[i],
													gamma = DOE$gamma[i],
													nu = DOE$nu[i]))))
													
	return(DFs)
}

registerDoSEQ()

# values greater than 0.05 are ideal 

ttests = sapply(1:nrow(DOE), function(i) 
							 t.test(svms[[i]]$actual - svms[[i]]$fit)$p.value)

# values close to zero are ideal 

stats = data.frame(t(sapply(1:nrow(DOE), function(i)
										 data.frame(accuracy(f = svms[[i]]$fit, 
															 x = svms[[i]]$actual), 
													row.names = 1))))

# lets add these test results to the DOE

DOE$t.pval = ttests
DOE = cbind(DOE, stats)
cnames = colnames(DOE)

DOE[,5:9] = sapply(5:9, function(i) DOE[,i] = as.numeric(DOE[,i]))
colnames(DOE) = cnames

head(DOE)

rm(svms, ttests, stats, cnames)

# lets filter out the top models

# lets plot the histograms of these reponses to find filter limits

par(mfrow = c(3, 2))
lapply(4:9, function(i) hist(DOE[,i], main = colnames(DOE)[i]))

head(DOE)

# apply filters

x = DOE[which(DOE$t.pval >= 0.9 & 
			  DOE$MAPE <= 5.2),]

par(mfrow = c(3, 2))
lapply(4:9, function(i) hist(x[,i], main = colnames(x)[i]))

x

# lets compare the auto fit with the chosen fit from the DOE

DOE[c(1, 440),]

svms = lapply(c(1, 440), function(i)
			   svm(formula = CRIMES ~ Feb + Mar + Apr + May + Jun + Jul + Aug + Sep + Oct + Nov + Dec + DAYS + UNEMPLOYMENT + UTILITY.GAS,
					data = DF,
					type = "nu-regression",
					cost = DOE$cost[i],
					gamma = DOE$gamma[i],
					nu = DOE$nu[i]))

# lets plot the model & residuals for the auto fit and chosen fit 

	# auto fit
	
plot_svm1 = modplot(mod = svms[[1]],
					  xlab = "Observation Number", 
					  ylab = "Crimes", 
					  main = "Crimes Per Week: Central Philly\nSupport Vector Machine: Fitted Model")

plot_svm1

plot_resid1 = residplots(actual = DF$CRIMES, 
						fit = as.numeric(fitted(svms[[1]])), 
						histlabel.y = -3, 
						basefont = 20)

grid.arrange(plot_resid1[[1]], 
			 plot_resid1[[2]], 
			 plot_resid1[[3]], 
			 plot_resid1[[4]], 
			 plot_resid1[[5]], 
			 plot_resid1[[6]], 
			 ncol = 3)		  

	# chosen fit

plot_svm2 = modplot(mod = svms[[2]],
					  xlab = "Observation Number", 
					  ylab = "Crimes", 
					  main = "Crimes Per Week: Central Philly\nSupport Vector Machine: Fitted Model")

plot_svm2

plot_resid2 = residplots(actual = DF$CRIMES, 
						fit = as.numeric(fitted(svms[[2]])), 
						histlabel.y = -3, 
						basefont = 20)

grid.arrange(plot_resid2[[1]], 
			 plot_resid2[[2]], 
			 plot_resid2[[3]], 
			 plot_resid2[[4]], 
			 plot_resid2[[5]], 
			 plot_resid2[[6]], 
			 ncol = 3)		  

SVM_CP_fit_plot = plot_svm2
SVM_CP_fit_mod = svms[[2]]

SVM_CP_fit_table = DOE[440, 4:9]
rownames(SVM_CP_fit_table) = "SVM"

rm(x, DOE, svms, plot_svm1, plot_resid1, plot_svm2, plot_resid2)

# ---- forecasting model ----

# lets create a DOE to test a range of values for each of the svm parameters of interest

DOE = expand.grid(cost = log.seq(from = 10^(-3), to = 10^3, length.out = 10),
				  gamma = seq(from = 2^(-15), to = 2^3, length.out = 10),
				  nu = seq(from = 0.3, to = 0.7, by = 0.1))

DOE = rbind(c(1, 1/14, 0.5), DOE)
DOE = DOE[!duplicated(DOE),]
rownames(DOE) = 1:nrow(DOE)

head(DOE)
tail(DOE)
nrow(DOE)

# lets compute the residuals for each svm model in the DOE

registerDoParallel(cores = cores) 

svms = foreach(i = 1:nrow(DOE)) %dopar%
{
	require(foreach)
	require(e1071)
	
	DFs = data.frame("actual" = DF$CRIMES[479:565],
					"fit" = as.numeric(predict(svm(formula = CRIMES ~ Feb + Mar + Apr + May + Jun + Jul + Aug + Sep + Oct + Nov + Dec + DAYS + UNEMPLOYMENT + UTILITY.GAS,
													data = DF[1:478,],
													type = "nu-regression",
													cost = DOE$cost[i],
													gamma = DOE$gamma[i],
													nu = DOE$nu[i]),
												newdata = DF[479:565, 4:17])))
													
	return(DFs)
}

registerDoSEQ()

# values greater than 0.05 are ideal 

ttests = sapply(1:nrow(DOE), function(i) 
							 t.test(svms[[i]]$actual - svms[[i]]$fit)$p.value)

# values close to zero are ideal 

stats = data.frame(t(sapply(1:nrow(DOE), function(i)
										 data.frame(accuracy(f = svms[[i]]$fit, 
															 x = svms[[i]]$actual), 
													row.names = 1))))

# lets add these test results to the DOE

DOE$t.pval = ttests
DOE = cbind(DOE, stats)
cnames = colnames(DOE)

DOE[,5:9] = sapply(5:9, function(i) DOE[,i] = as.numeric(DOE[,i]))
colnames(DOE) = cnames

head(DOE)

rm(svms, ttests, stats, cnames)

# lets filter out the top models

# lets plot the histograms of these reponses to find filter limits

par(mfrow = c(3, 2))
lapply(4:9, function(i) hist(DOE[,i], main = colnames(DOE)[i]))

head(DOE)

# apply filters

x = DOE[which(DOE$MAPE <= 15 & 
				DOE$ME >= -40),]

par(mfrow = c(3, 2))
lapply(4:9, function(i) hist(x[,i], main = colnames(x)[i]))

x

# lets compare the auto fit with the chosen fit from the DOE

DOE[c(1, 311),]

svms = lapply(c(1, 311), function(i)
			   svm(formula = CRIMES ~ Feb + Mar + Apr + May + Jun + Jul + Aug + Sep + Oct + Nov + Dec + DAYS + UNEMPLOYMENT + UTILITY.GAS,
					data = DF[1:478,],
					type = "nu-regression",
					cost = DOE$cost[i],
					gamma = DOE$gamma[i],
					nu = DOE$nu[i]))

# lets plot the model & residuals for the auto fit and chosen fit 

	# auto fit
	
plot_svm1 = modplot(mod = svms[[1]],
					  newdata = DF[479:565,4:17],
					  n.ahead = 87,
					  actuals = DF$CRIMES,
					  limits = c(479, 565),
					  xlab = "Observation Number", 
					  ylab = "Crimes", 
					  main = "Crimes Per Week: Central Philly\nSupport Vector Machine: Forecasting Model")

plot_svm1

plot_resid1 = residplots(actual = DF$CRIMES[479:565], 
						fit = as.numeric(predict(svms[[1]], newdata = DF[479:565, 4:17])), 
						histlabel.y = -0.5, 
						basefont = 20)

grid.arrange(plot_resid1[[1]], 
			 plot_resid1[[2]], 
			 plot_resid1[[3]], 
			 plot_resid1[[4]], 
			 plot_resid1[[5]], 
			 plot_resid1[[6]], 
			 ncol = 3)		  

	# chosen fit

plot_svm2 = modplot(mod = svms[[2]],
					newdata = DF[479:565,4:17],
					n.ahead = 87,
					actuals = DF$CRIMES,
					limits = c(479, 565),
					xlab = "Observation Number", 
					ylab = "Crimes", 
					main = "Crimes Per Week: Central Philly\nSupport Vector Machine: Forecasting Model")

plot_svm2

plot_resid2 = residplots(actual = DF$CRIMES[479:565], 
						fit = as.numeric(predict(svms[[2]], newdata = DF[479:565, 4:17])), 
						histlabel.y = -0.5, 
						basefont = 20)

grid.arrange(plot_resid2[[1]], 
			 plot_resid2[[2]], 
			 plot_resid2[[3]], 
			 plot_resid2[[4]], 
			 plot_resid2[[5]], 
			 plot_resid2[[6]], 
			 ncol = 3)		  

SVM_CP_forecast_plot = plot_svm2
SVM_CP_forecast_mod = svms[[2]]

SVM_CP_forecast_table = DOE[311, 4:9]
rownames(SVM_CP_forecast_table) = "SVM"

rm(x, DOE, svms, plot_svm1, plot_resid1, plot_svm2, plot_resid2)

}

# ---- West Philly --------------------------------------------------------------

{

# ---- fitted model ----

# extract the WP subset from crime

DF = crime[which(crime$Area == "WP"),]

# lets create a DOE to test a range of values for each of the svm parameters of interest

DOE = expand.grid(cost = log.seq(from = 10^(-3), to = 10^3, length.out = 10),
				  gamma = seq(from = 2^(-15), to = 2^3, length.out = 10),
				  nu = seq(from = 0.3, to = 0.7, by = 0.1))

DOE = rbind(c(1, 1/14, 0.5), DOE)
DOE = DOE[!duplicated(DOE),]
rownames(DOE) = 1:nrow(DOE)

head(DOE)
tail(DOE)
nrow(DOE)

# lets compute the residuals for each svm model in the DOE

registerDoParallel(cores = cores) 

svms = foreach(i = 1:nrow(DOE)) %dopar%
{
	require(foreach)
	require(e1071)
	
	DFs = data.frame("actual" = DF$CRIMES,
					"fit" = as.numeric(fitted(svm(formula = CRIMES ~ Feb + Mar + Apr + May + Jun + Jul + Aug + Sep + Oct + Nov + Dec + DAYS + UNEMPLOYMENT + UTILITY.GAS,
													data = DF,
													type = "nu-regression",
													cost = DOE$cost[i],
													gamma = DOE$gamma[i],
													nu = DOE$nu[i]))))
													
	return(DFs)
}

registerDoSEQ()

# values greater than 0.05 are ideal 

ttests = sapply(1:nrow(DOE), function(i) 
							 t.test(svms[[i]]$actual - svms[[i]]$fit)$p.value)

# values close to zero are ideal 

stats = data.frame(t(sapply(1:nrow(DOE), function(i)
										 data.frame(accuracy(f = svms[[i]]$fit, 
															 x = svms[[i]]$actual), 
													row.names = 1))))

# lets add these test results to the DOE

DOE$t.pval = ttests
DOE = cbind(DOE, stats)
cnames = colnames(DOE)

DOE[,5:9] = sapply(5:9, function(i) DOE[,i] = as.numeric(DOE[,i]))
colnames(DOE) = cnames

head(DOE)

rm(svms, ttests, stats, cnames)

# lets filter out the top models

# lets plot the histograms of these reponses to find filter limits

par(mfrow = c(3, 2))
lapply(4:9, function(i) hist(DOE[,i], main = colnames(DOE)[i]))

head(DOE)

# apply filters

x = DOE[which(DOE$t.pval >= 0.9 & 
			  DOE$MAPE <= 6),]

par(mfrow = c(3, 2))
lapply(4:9, function(i) hist(x[,i], main = colnames(x)[i]))

x

# lets compare the auto fit with the chosen fit from the DOE

DOE[c(1, 320),]

svms = lapply(c(1, 320), function(i)
			   svm(formula = CRIMES ~ Feb + Mar + Apr + May + Jun + Jul + Aug + Sep + Oct + Nov + Dec + DAYS + UNEMPLOYMENT + UTILITY.GAS,
					data = DF,
					type = "nu-regression",
					cost = DOE$cost[i],
					gamma = DOE$gamma[i],
					nu = DOE$nu[i]))

# lets plot the model & residuals for the auto fit and chosen fit 

	# auto fit
	
plot_svm1 = modplot(mod = svms[[1]],
					  xlab = "Observation Number", 
					  ylab = "Crimes", 
					  main = "Crimes Per Week: West Philly\nSupport Vector Machine: Fitted Model")

plot_svm1

plot_resid1 = residplots(actual = DF$CRIMES, 
						fit = as.numeric(fitted(svms[[1]])), 
						histlabel.y = -3, 
						basefont = 20)

grid.arrange(plot_resid1[[1]], 
			 plot_resid1[[2]], 
			 plot_resid1[[3]], 
			 plot_resid1[[4]], 
			 plot_resid1[[5]], 
			 plot_resid1[[6]], 
			 ncol = 3)		  

	# chosen fit

plot_svm2 = modplot(mod = svms[[2]],
					  xlab = "Observation Number", 
					  ylab = "Crimes", 
					  main = "Crimes Per Week: West Philly\nSupport Vector Machine: Fitted Model")

plot_svm2

plot_resid2 = residplots(actual = DF$CRIMES, 
						fit = as.numeric(fitted(svms[[2]])), 
						histlabel.y = -3, 
						basefont = 20)

grid.arrange(plot_resid2[[1]], 
			 plot_resid2[[2]], 
			 plot_resid2[[3]], 
			 plot_resid2[[4]], 
			 plot_resid2[[5]], 
			 plot_resid2[[6]], 
			 ncol = 3)		  

SVM_WP_fit_plot = plot_svm2
SVM_WP_fit_mod = svms[[2]]

SVM_WP_fit_table = DOE[320, 4:9]
rownames(SVM_WP_fit_table) = "SVM"

rm(x, DOE, svms, plot_svm1, plot_resid1, plot_svm2, plot_resid2)

# ---- forecasting model ----

# lets create a DOE to test a range of values for each of the svm parameters of interest

DOE = expand.grid(cost = log.seq(from = 10^(-3), to = 10^3, length.out = 10),
				  gamma = seq(from = 2^(-15), to = 2^3, length.out = 10),
				  nu = seq(from = 0.3, to = 0.7, by = 0.1))

DOE = rbind(c(1, 1/14, 0.5), DOE)
DOE = DOE[!duplicated(DOE),]
rownames(DOE) = 1:nrow(DOE)

head(DOE)
tail(DOE)
nrow(DOE)

# lets compute the residuals for each svm model in the DOE

registerDoParallel(cores = cores) 

svms = foreach(i = 1:nrow(DOE)) %dopar%
{
	require(foreach)
	require(e1071)
	
	DFs = data.frame("actual" = DF$CRIMES[479:565],
					"fit" = as.numeric(predict(svm(formula = CRIMES ~ Feb + Mar + Apr + May + Jun + Jul + Aug + Sep + Oct + Nov + Dec + DAYS + UNEMPLOYMENT + UTILITY.GAS,
													data = DF[1:478,],
													type = "nu-regression",
													cost = DOE$cost[i],
													gamma = DOE$gamma[i],
													nu = DOE$nu[i]),
												newdata = DF[479:565, 4:17])))
													
	return(DFs)
}

registerDoSEQ()

# values greater than 0.05 are ideal 

ttests = sapply(1:nrow(DOE), function(i) 
							 t.test(svms[[i]]$actual - svms[[i]]$fit)$p.value)

# values close to zero are ideal 

stats = data.frame(t(sapply(1:nrow(DOE), function(i)
										 data.frame(accuracy(f = svms[[i]]$fit, 
															 x = svms[[i]]$actual), 
													row.names = 1))))

# lets add these test results to the DOE

DOE$t.pval = ttests
DOE = cbind(DOE, stats)
cnames = colnames(DOE)

DOE[,5:9] = sapply(5:9, function(i) DOE[,i] = as.numeric(DOE[,i]))
colnames(DOE) = cnames

head(DOE)

rm(svms, ttests, stats, cnames)

# lets filter out the top models

# lets plot the histograms of these reponses to find filter limits

par(mfrow = c(3, 2))
lapply(4:9, function(i) hist(DOE[,i], main = colnames(DOE)[i]))

head(DOE)

# apply filters

x = DOE[which(DOE$MAPE <= 17.5 & 
				DOE$t.pval >= 0.6),]

par(mfrow = c(3, 2))
lapply(4:9, function(i) hist(x[,i], main = colnames(x)[i]))

x

# lets compare the auto fit with the chosen fit from the DOE

DOE[c(1, 18),]

svms = lapply(c(1, 18), function(i)
			   svm(formula = CRIMES ~ Feb + Mar + Apr + May + Jun + Jul + Aug + Sep + Oct + Nov + Dec + DAYS + UNEMPLOYMENT + UTILITY.GAS,
					data = DF[1:478,],
					type = "nu-regression",
					cost = DOE$cost[i],
					gamma = DOE$gamma[i],
					nu = DOE$nu[i]))

# lets plot the model & residuals for the auto fit and chosen fit 

	# auto fit
	
plot_svm1 = modplot(mod = svms[[1]],
					  newdata = DF[479:565,4:17],
					  n.ahead = 87,
					  actuals = DF$CRIMES,
					  limits = c(479, 565),
					  xlab = "Observation Number", 
					  ylab = "Crimes", 
					  main = "Crimes Per Week: West Philly\nSupport Vector Machine: Forecasting Model")

plot_svm1

plot_resid1 = residplots(actual = DF$CRIMES[479:565], 
						fit = as.numeric(predict(svms[[1]], newdata = DF[479:565, 4:17])), 
						histlabel.y = -0.5, 
						basefont = 20)

grid.arrange(plot_resid1[[1]], 
			 plot_resid1[[2]], 
			 plot_resid1[[3]], 
			 plot_resid1[[4]], 
			 plot_resid1[[5]], 
			 plot_resid1[[6]], 
			 ncol = 3)		  

	# chosen fit

plot_svm2 = modplot(mod = svms[[2]],
					newdata = DF[479:565,4:17],
					n.ahead = 87,
					actuals = DF$CRIMES,
					limits = c(479, 565),
					xlab = "Observation Number", 
					ylab = "Crimes", 
					main = "Crimes Per Week: West Philly\nSupport Vector Machine: Forecasting Model")

plot_svm2

plot_resid2 = residplots(actual = DF$CRIMES[479:565], 
						fit = as.numeric(predict(svms[[2]], newdata = DF[479:565, 4:17])), 
						histlabel.y = -0.5, 
						basefont = 20)

grid.arrange(plot_resid2[[1]], 
			 plot_resid2[[2]], 
			 plot_resid2[[3]], 
			 plot_resid2[[4]], 
			 plot_resid2[[5]], 
			 plot_resid2[[6]], 
			 ncol = 3)		  

SVM_WP_forecast_plot = plot_svm2
SVM_WP_forecast_mod = svms[[2]]

SVM_WP_forecast_table = DOE[18, 4:9]
rownames(SVM_WP_forecast_table) = "SVM"

rm(x, DOE, svms, plot_svm1, plot_resid1, plot_svm2, plot_resid2)

}

# ---- North Philly --------------------------------------------------------------

{

# ---- fitted model ----

# extract the NP subset from crime

DF = crime[which(crime$Area == "NP"),]

# lets create a DOE to test a range of values for each of the svm parameters of interest

DOE = expand.grid(cost = log.seq(from = 10^(-3), to = 10^3, length.out = 10),
				  gamma = seq(from = 2^(-15), to = 2^3, length.out = 10),
				  nu = seq(from = 0.3, to = 0.7, by = 0.1))

DOE = rbind(c(1, 1/14, 0.5), DOE)
DOE = DOE[!duplicated(DOE),]
rownames(DOE) = 1:nrow(DOE)

head(DOE)
tail(DOE)
nrow(DOE)

# lets compute the residuals for each svm model in the DOE

registerDoParallel(cores = cores) 

svms = foreach(i = 1:nrow(DOE)) %dopar%
{
	require(foreach)
	require(e1071)
	
	DFs = data.frame("actual" = DF$CRIMES,
					"fit" = as.numeric(fitted(svm(formula = CRIMES ~ Feb + Mar + Apr + May + Jun + Jul + Aug + Sep + Oct + Nov + Dec + DAYS + UNEMPLOYMENT + UTILITY.GAS,
													data = DF,
													type = "nu-regression",
													cost = DOE$cost[i],
													gamma = DOE$gamma[i],
													nu = DOE$nu[i]))))
													
	return(DFs)
}

registerDoSEQ()

# values greater than 0.05 are ideal 

ttests = sapply(1:nrow(DOE), function(i) 
							 t.test(svms[[i]]$actual - svms[[i]]$fit)$p.value)

# values close to zero are ideal 

stats = data.frame(t(sapply(1:nrow(DOE), function(i)
										 data.frame(accuracy(f = svms[[i]]$fit, 
															 x = svms[[i]]$actual), 
													row.names = 1))))

# lets add these test results to the DOE

DOE$t.pval = ttests
DOE = cbind(DOE, stats)
cnames = colnames(DOE)

DOE[,5:9] = sapply(5:9, function(i) DOE[,i] = as.numeric(DOE[,i]))
colnames(DOE) = cnames

head(DOE)

rm(svms, ttests, stats, cnames)

# lets filter out the top models

# lets plot the histograms of these reponses to find filter limits

par(mfrow = c(3, 2))
lapply(4:9, function(i) hist(DOE[,i], main = colnames(DOE)[i]))

head(DOE)

# apply filters

x = DOE[which(DOE$t.pval >= 0.9 & 
			  DOE$MAPE <= 6 &
			  abs(DOE$ME) <= 0.1 &
			  DOE$MAE <= 71),]

par(mfrow = c(3, 2))
lapply(4:9, function(i) hist(x[,i], main = colnames(x)[i]))

x

# lets compare the auto fit with the chosen fit from the DOE

DOE[c(1, 340),]

svms = lapply(c(1, 340), function(i)
			   svm(formula = CRIMES ~ Feb + Mar + Apr + May + Jun + Jul + Aug + Sep + Oct + Nov + Dec + DAYS + UNEMPLOYMENT + UTILITY.GAS,
					data = DF,
					type = "nu-regression",
					cost = DOE$cost[i],
					gamma = DOE$gamma[i],
					nu = DOE$nu[i]))

# lets plot the model & residuals for the auto fit and chosen fit 

	# auto fit
	
plot_svm1 = modplot(mod = svms[[1]],
					  xlab = "Observation Number", 
					  ylab = "Crimes", 
					  main = "Crimes Per Week: North Philly\nSupport Vector Machine: Fitted Model")

plot_svm1

plot_resid1 = residplots(actual = DF$CRIMES, 
						fit = as.numeric(fitted(svms[[1]])), 
						histlabel.y = -3, 
						basefont = 20)

grid.arrange(plot_resid1[[1]], 
			 plot_resid1[[2]], 
			 plot_resid1[[3]], 
			 plot_resid1[[4]], 
			 plot_resid1[[5]], 
			 plot_resid1[[6]], 
			 ncol = 3)		  

	# chosen fit

plot_svm2 = modplot(mod = svms[[2]],
					  xlab = "Observation Number", 
					  ylab = "Crimes", 
					  main = "Crimes Per Week: North Philly\nSupport Vector Machine: Fitted Model")

plot_svm2

plot_resid2 = residplots(actual = DF$CRIMES, 
						fit = as.numeric(fitted(svms[[2]])), 
						histlabel.y = -3, 
						basefont = 20)

grid.arrange(plot_resid2[[1]], 
			 plot_resid2[[2]], 
			 plot_resid2[[3]], 
			 plot_resid2[[4]], 
			 plot_resid2[[5]], 
			 plot_resid2[[6]], 
			 ncol = 3)		  

SVM_NP_fit_plot = plot_svm2
SVM_NP_fit_mod = svms[[2]]

SVM_NP_fit_table = DOE[340, 4:9]
rownames(SVM_NP_fit_table) = "SVM"

rm(x, DOE, svms, plot_svm1, plot_resid1, plot_svm2, plot_resid2)

# ---- forecasting model ----

# lets create a DOE to test a range of values for each of the svm parameters of interest

DOE = expand.grid(cost = log.seq(from = 10^(-3), to = 10^3, length.out = 10),
				  gamma = seq(from = 2^(-15), to = 2^3, length.out = 10),
				  nu = seq(from = 0.3, to = 0.7, by = 0.1))

DOE = rbind(c(1, 1/14, 0.5), DOE)
DOE = DOE[!duplicated(DOE),]
rownames(DOE) = 1:nrow(DOE)

head(DOE)
tail(DOE)
nrow(DOE)

# lets compute the residuals for each svm model in the DOE

registerDoParallel(cores = cores) 

svms = foreach(i = 1:nrow(DOE)) %dopar%
{
	require(foreach)
	require(e1071)
	
	DFs = data.frame("actual" = DF$CRIMES[479:565],
					"fit" = as.numeric(predict(svm(formula = CRIMES ~ Feb + Mar + Apr + May + Jun + Jul + Aug + Sep + Oct + Nov + Dec + DAYS + UNEMPLOYMENT + UTILITY.GAS,
													data = DF[1:478,],
													type = "nu-regression",
													cost = DOE$cost[i],
													gamma = DOE$gamma[i],
													nu = DOE$nu[i]),
												newdata = DF[479:565, 4:17])))
													
	return(DFs)
}

registerDoSEQ()

# values greater than 0.05 are ideal 

ttests = sapply(1:nrow(DOE), function(i) 
							 t.test(svms[[i]]$actual - svms[[i]]$fit)$p.value)

# values close to zero are ideal 

stats = data.frame(t(sapply(1:nrow(DOE), function(i)
										 data.frame(accuracy(f = svms[[i]]$fit, 
															 x = svms[[i]]$actual), 
													row.names = 1))))

# lets add these test results to the DOE

DOE$t.pval = ttests
DOE = cbind(DOE, stats)
cnames = colnames(DOE)

DOE[,5:9] = sapply(5:9, function(i) DOE[,i] = as.numeric(DOE[,i]))
colnames(DOE) = cnames

head(DOE)

rm(svms, ttests, stats, cnames)

# lets filter out the top models

# lets plot the histograms of these reponses to find filter limits

par(mfrow = c(3, 2))
lapply(4:9, function(i) hist(DOE[,i], main = colnames(DOE)[i]))

head(DOE)

# apply filters

x = DOE[which(DOE$MAPE <= 15 &
			  DOE$t.pval >= 0.39),]

par(mfrow = c(3, 2))
lapply(4:9, function(i) hist(x[,i], main = colnames(x)[i]))

x

# lets compare the auto fit with the chosen fit from the DOE

DOE[c(1, 218),]

svms = lapply(c(1, 218), function(i)
			   svm(formula = CRIMES ~ Feb + Mar + Apr + May + Jun + Jul + Aug + Sep + Oct + Nov + Dec + DAYS + UNEMPLOYMENT + UTILITY.GAS,
					data = DF[1:478,],
					type = "nu-regression",
					cost = DOE$cost[i],
					gamma = DOE$gamma[i],
					nu = DOE$nu[i]))

# lets plot the model & residuals for the auto fit and chosen fit 

	# auto fit
	
plot_svm1 = modplot(mod = svms[[1]],
					  newdata = DF[479:565,4:17],
					  n.ahead = 87,
					  actuals = DF$CRIMES,
					  limits = c(479, 565),
					  xlab = "Observation Number", 
					  ylab = "Crimes", 
					  main = "Crimes Per Week: North Philly\nSupport Vector Machine: Forecasting Model")

plot_svm1

plot_resid1 = residplots(actual = DF$CRIMES[479:565], 
						fit = as.numeric(predict(svms[[1]], newdata = DF[479:565, 4:17])), 
						histlabel.y = -0.5, 
						basefont = 20)

grid.arrange(plot_resid1[[1]], 
			 plot_resid1[[2]], 
			 plot_resid1[[3]], 
			 plot_resid1[[4]], 
			 plot_resid1[[5]], 
			 plot_resid1[[6]], 
			 ncol = 3)		  

	# chosen fit

plot_svm2 = modplot(mod = svms[[2]],
					newdata = DF[479:565,4:17],
					n.ahead = 87,
					actuals = DF$CRIMES,
					limits = c(479, 565),
					xlab = "Observation Number", 
					ylab = "Crimes", 
					main = "Crimes Per Week: North Philly\nSupport Vector Machine: Forecasting Model")

plot_svm2

plot_resid2 = residplots(actual = DF$CRIMES[479:565], 
						fit = as.numeric(predict(svms[[2]], newdata = DF[479:565, 4:17])), 
						histlabel.y = -0.5, 
						basefont = 20)

grid.arrange(plot_resid2[[1]], 
			 plot_resid2[[2]], 
			 plot_resid2[[3]], 
			 plot_resid2[[4]], 
			 plot_resid2[[5]], 
			 plot_resid2[[6]], 
			 ncol = 3)		  

SVM_NP_forecast_plot = plot_svm1
SVM_NP_forecast_mod = svms[[1]]

SVM_NP_forecast_table = DOE[1, 4:9]
rownames(SVM_NP_forecast_table) = "SVM"

rm(x, DOE, svms, plot_svm1, plot_resid1, plot_svm2, plot_resid2)

}

}

# -----------------------------------------------------------------------------------
# ---- Summary of Models ------------------------------------------------------------
# -----------------------------------------------------------------------------------

{

# ---- Central Philly --------------------------------------------------------------

# ---- tabular summaries ----

CP_fit_table = rbind(LM_CP_fit_table, HW_CP_fit_table, ARIMA_CP_fit_table, DR_CP_fit_table, NNET_CP_fit_table, RF_CP_fit_table, SVM_CP_fit_table)

CP_forecast_table = rbind(LM_CP_forecast_table, HW_CP_forecast_table, ARIMA_CP_forecast_table, DR_CP_forecast_table, NNET_CP_forecast_table, RF_CP_forecast_table, SVM_CP_forecast_table)

# performance metrics of all fitted models
 
CP_fit_table

# performance metrics of all forecasting models

CP_forecast_table

# ---- graphic summary ----

CP_forecast_plots = arrangeGrob(LM_CP_forecast_plot, HW_CP_forecast_plot, ncol = 2)

# top two forecasts

grid.draw(CP_forecast_plots)

# ---- West Philly --------------------------------------------------------------

# ---- tabular summaries ----

WP_fit_table = rbind(LM_WP_fit_table, HW_WP_fit_table, ARIMA_WP_fit_table, DR_WP_fit_table, NNET_WP_fit_table, RF_WP_fit_table, SVM_WP_fit_table)

WP_forecast_table = rbind(LM_WP_forecast_table, HW_WP_forecast_table, ARIMA_WP_forecast_table, DR_WP_forecast_table, NNET_WP_forecast_table, RF_WP_forecast_table, SVM_WP_forecast_table)

# performance metrics of all fitted models

WP_fit_table

# performance metrics of all forecasting models

WP_forecast_table

# ---- graphic summary ----

WP_forecast_plots = arrangeGrob(HW_WP_forecast_plot, DR_WP_forecast_plot, ncol = 2)

# top two forecasts

grid.draw(WP_forecast_plots)

# ---- North Philly --------------------------------------------------------------

# ---- tabular summaries ----

NP_fit_table = rbind(LM_NP_fit_table, HW_NP_fit_table, ARIMA_NP_fit_table, DR_NP_fit_table, NNET_NP_fit_table, RF_NP_fit_table, SVM_NP_fit_table)

NP_forecast_table = rbind(LM_NP_forecast_table, HW_NP_forecast_table, ARIMA_NP_forecast_table, DR_NP_forecast_table, NNET_NP_forecast_table, RF_NP_forecast_table, SVM_NP_forecast_table)

# performance metrics of all fitted models
 
NP_fit_table

# performance metrics of all forecasting models

NP_forecast_table

# ---- graphic summary ----

NP_forecast_plots = arrangeGrob(HW_NP_forecast_plot, DR_NP_forecast_plot, ncol = 2)

# top two forecasts

grid.draw(NP_forecast_plots)

}

CP_fit_table
WP_fit_table
NP_fit_table

CP_forecast_table
WP_forecast_table
NP_forecast_table

LM_CP_forecast_plot
RF_CP_forecast_plot
SVM_CP_forecast_plot

DR_WP_forecast_plot
RF_WP_forecast_plot
SVM_WP_forecast_plot

DR_NP_forecast_plot
RF_NP_forecast_plot
SVM_NP_forecast_plot




