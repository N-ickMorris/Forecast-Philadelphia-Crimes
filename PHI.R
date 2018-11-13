# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ---- Packages ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

{

require(data.table)
require(ggmap)
require(ggplot2)
require(RColorBrewer)
require(gridExtra)
require(plyr)
require(dplyr)
require(forecast)
require(car)
require(GGally)
require(grid)
require(EnvStats)
require(nortest)
require(fields)
plot(1:5)

}

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ---- Functions ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

{

# ---------------------------------------------------------------------------------------------------------------------------------------------------

types = function(dat)
{
	dat = data.frame(dat)
  
	Column = sapply(1:NCOL(dat), function(i)
	(
		colnames(dat)[i]
    ))
  
	Data_Type = sapply(1:NCOL(dat), function(i)
    (
		class(dat[,i])
    ))
  
	results = data.frame(cbind(Column, Data_Type))	
	results
}

# ---------------------------------------------------------------------------------------------------------------------------------------------------

ggacf = function(x, n = NULL, partial = FALSE, conf.level = 0.95, main = "ACF Plot", xlab = "Lag", ylab = "Autocorrelation", basefont = 20) 
{
	require(ggplot2)
	
	if(class(n) == "NULL")
	{
		n = length(x) - 2
	}
	
	ciline = qnorm((1 - conf.level) / 2) / sqrt(length(x))
	
	if(partial == TRUE)
	{
		bacf = pacf(x, lag.max = n, plot = FALSE)
	} else
	{
		bacf = acf(x, lag.max = n, plot = FALSE)
	}
	
	bacfdf = with(bacf, data.frame(lag, acf))
	
	if(partial == FALSE)
	{
		bacfdf = bacfdf[-1,]
	}
	
	acfplot = ggplot(bacfdf, aes(x = lag, y = acf)) + 
			  geom_bar(stat = "identity", position = "dodge", fill = "black") +
			  geom_hline(yintercept = -ciline, color = "blue", size = 1) +
			  geom_hline(yintercept = ciline, color = "blue", size = 1) +
			  geom_hline(yintercept = 0, color = "red", size = 1) +
			  labs(x = xlab, y = ylab) +
			  ggtitle(main) +
			  theme_light(base_size = basefont) +
			  theme(legend.position = "none")

	return(acfplot)
}

# ---------------------------------------------------------------------------------------------------------------------------------------------------

variogramDF = function(x, n = NULL)
{
	if(class(n) == "NULL")(n = length(x) - 2)
	
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

# ---------------------------------------------------------------------------------------------------------------------------------------------------

statslm = function(model)
{
    # Calculate the Predictive Residuals of 'model'
    PR = residuals(model)/(1 - lm.influence(model)$hat)
    
    # Calculate the Predicted Residual Sum of Squares
    PRESS=sum(PR^2)
    
    # Calculate the Total Sum of Squares
    TSS=sum(anova(model)$"Sum Sq")
    
    # Summary
    return(list("Prediction" = c("R^2-Pred" = round(1 - (PRESS/TSS), 4)),
				"Fitness" = c("AIC" = AIC(model),
							  "BIC" = BIC(model))))
}

# ---------------------------------------------------------------------------------------------------------------------------------------------------

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
			  theme_light(base_size = basefont) +
			  theme(legend.position = "none")
    
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
    		theme_light(base_size = basefont) +
			theme(legend.position = "none")
						
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
			  theme_light(base_size = basefont) +
			  theme(legend.position = "none")
        
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
				theme_light(base_size = basefont) +
				theme(legend.position = "none")
	
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
	 		   theme_light(base_size = basefont) +
			   theme(legend.key.size = unit(.25, "in"),
					 legend.position = "bottom")
	
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
			theme_light(base_size = basefont) +
			theme(legend.position = "none")

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

# ---------------------------------------------------------------------------------------------------------------------------------------------------

testslm = function(LM, LAG = 1)
{
    tab1 = t.test(LM$residuals)
    
	require(car)
    tab2 = ncvTest(LM)
	
	if(length(LM$residuals) > 5000)
	{
		require(nortest)
		tab3 = ad.test(x = LM$residuals)
	} else
	{
		tab3 = shapiro.test(LM$residuals)
	}
	
    tab4 = durbinWatsonTest(LM, LAG)
    tab5 = summary(LM)
    
	if(length(LM[[1]][which(names(LM[[1]]) != "(Intercept)")]) > 1)
	{
		tab6 = vif(LM)
		
		Results = list("Assumption 1 - Residuals have an Expected Value of Zero - Zero Within CI" = tab1, 
					   "Assumption 2 - Residuals have Constant Variance - P-Value > 0.05" = tab2,
					   "Assumption 3 - Residuals are Normally Distributed - P-Value > 0.05" = tab3,
					   "Assumption 4 - Residuals are Uncorrelated - P-value > 0.05" = tab4,
					   "Assumption 5 - The Relationship between the Response and Regressors is Correct - All P-Value < 0.1" = tab5,
					   "Assumption 6 - The Regressors are Independent - All vif < 10" = tab6)
	} else
	{
		Results = list("Assumption 1 - Residuals have an Expected Value of Zero - Zero Within CI" = tab1,
					   "Assumption 2 - Residuals have Constant Variance - P-Value > 0.05" = tab2,
					   "Assumption 3 - Residuals are Normally Distributed - P-Value > 0.05" = tab3,
					   "Assumption 4 - Residuals are Uncorrelated - P-value > 0.05" = tab4,
					   "Assumption 5 - The Relationship between the Response and Regressors is Correct - All P-Value < 0.1" = tab5)
	}
	
    return(Results)
}

# ---------------------------------------------------------------------------------------------------------------------------------------------------

termplots = function(dat, model, basefont = 15)
{
	require(ggplot2)
	
	vars = names(model$coefficients)[which(names(model$coefficients) %in% names(dat))]
	pterms = predict(model, type = "terms")
	pres = apply(pterms, 2, function(i) i + resid(model))
	DF = lapply(1:length(vars), function(i) data.frame("x" = model$model[,(i + 1)], "y" = pres[,i]))
	ablines = lapply(1:length(vars), function(i) lm(pres[,i] ~ model$model[,(i + 1)]))
	
	p = lapply(1:length(DF), function(i)
										ggplot(DF[[i]], aes(x = x, y = y)) + 
										geom_point(color = "black") +
										geom_abline(intercept = ablines[[i]]$coefficients[[1]], slope = ablines[[i]]$coefficients[[2]], size = 2, color = "cornflowerblue") +
										labs(x = vars[i], y = "Partial Residuals") +
										ggtitle(paste0("Partial Residuals of ", vars[i])) +
										theme_light(base_size = basefont))
	
	names(p) = vars
	
	return(p)
}

# ---------------------------------------------------------------------------------------------------------------------------------------------------

predplot = function(obj, xaxis = NULL, actuals = NULL, n.ahead = NULL, level = c(80, 95), newdata = NULL, xreg = NULL, limits = NULL, xlab = "Observation", ylab = "Value", main = "Fitted v. Actual\nPrediction Plot", basefont = 20)
{
	require(ggplot2)
	require(forecast)
	
	# build fits

	if(class(obj)[1] == "HoltWinters")
	{
		fits = as.numeric(fitted(obj)[,1])
		fits = c(rep(NA, length(obj$x) - length(fits)), fits)
	} else
	{
		fits = as.numeric(fitted(obj))
	}

	# build actuals

	if(class(actuals) == "NULL")
	{
		if(class(obj)[1] == "HoltWinters")
		{
			actuals = as.numeric(obj$x)
		} else
		{
			actuals = as.numeric(resid(obj)) + fits
		}
		
		if(class(n.ahead) != "NULL")
		{
			actuals = c(actuals, rep(NA, n.ahead))
		}
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
	
	upper1 = as.numeric(rep(NA, length(fits)))
	upper2 = upper1
	lower1 = upper1
	lower2 = upper1
		
	# initialize a color vector to seperate training data from testing data
		
	color = factor(rep(0, length(fits)), levels = c(0, 1))
	
	# compute predictions if n.ahead is specified 
		# also update fits, upper and lower prediction interval values, and color 
		
	if(class(n.ahead) != "NULL")
	{
		if(class(obj)[1] == "lm")
		{
			predictions = data.frame(forecast(obj, newdata = newdata, h = n.ahead, level = level))
		} else if(class(obj)[1] == "glm")
		{
			predictions = data.frame(predict(obj, newdata = newdata, n.ahead = n.ahead, se.fit = TRUE))
			
			predictions = data.frame("fit" = predictions$fit, 
									 "lower2" = predictions$fit - predictions$se.fit,
									 "upper2" = predictions$fit + predictions$se.fit,
									 "lower1" = predictions$fit - (3 * predictions$se.fit),
									 "upper1" = predictions$fit + (3 * predictions$se.fit))
			
		} else if(class(obj)[1] == "ARIMA")
		{
			predictions = data.frame(forecast(obj, h = n.ahead, level = level, xreg = xreg))
		} else
		{
			predictions = data.frame(forecast(obj, h = n.ahead, level = level))
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
	
	p = ggplot(DF, aes(x = Observation, y = Actuals, color = Color)) + 
		scale_color_manual(values = c("black", "red"), drop = TRUE, limits = levels(DF$Color)) +
		geom_point(na.rm = TRUE) + 
		geom_line(aes(x = Observation, y = Fitted), color = "blue", na.rm = TRUE) + 
		geom_ribbon(aes(ymin = Lower1, ymax = Upper1), alpha = .1, color = NA) +
		geom_ribbon(aes(ymin = Lower2, ymax = Upper2), alpha = .125, color = NA) +
		labs(x = xlab, y = ylab) +
		ggtitle(main) + 
		theme_light(base_size = basefont) +
		theme(legend.position = "none")

	return(p)
}

# ---------------------------------------------------------------------------------------------------------------------------------------------------

Arima.if = function(...) 
{
	return(tryCatch(Arima(...), error = function(e) NA))
}

arima.if = function(...) 
{
	return(tryCatch(arima(...), error = function(e) NA))
}

predict.if = function(...)
{
	return(tryCatch(predict(...), error = function(e) NA))
}

fitted.if = function(...)
{
	return(tryCatch(fitted(...), error = function(e) NA))
}

forecast.if = function(...)
{
	return(tryCatch(forecast(...), error = function(e) NA))
}

# ---------------------------------------------------------------------------------------------------------------------------------------------------


}

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ---- Importing Data ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

{

# PHI crime data: crime up through 2016-08-18
  # https://www.kaggle.com/mchirico/philadelphiacrimedata
  # choose: crime.csv

PHI = data.table(read.csv(file.choose(), header = TRUE))

# lets evaluate the data types of the columns

head(PHI)
types(PHI)

# lets evaluate each column for NA's
# this will show us the proportion of NA's in each of the columns

na = sapply(names(PHI), function(i) NROW(na.omit(PHI, cols = i, invert = TRUE)) / NROW(PHI))
names(na) = colnames(PHI)
na = data.frame("Percent.NA" = na)
na

# given that there is a small proportion of NA's in 4 columns, lets remove the rows with NA's 

PHI = data.table(na.omit(PHI))

# lets look at the daylight saving hours in this data set
# this shows us how R sequences datetimes based on daylight savings

daylightsavings = data.frame("FORWARD" = seq.POSIXt(from = as.POSIXct("2006-04-02 00:00:01"), to = as.POSIXct("2006-04-02 04:00:01"), by = "hour"), 
							 "BACKWARD" = seq.POSIXt(from = as.POSIXct("2006-10-29 00:00:01"), to = as.POSIXct("2006-10-29 02:00:01"), by = "hour"))

daylightsavings

# lets create a dataframe of date-hours that shouldn't be present in the data set based on FORWARD daylight savings

forward = data.frame("DATE" = as.Date(c("2006-04-02", "2007-03-11", "2008-03-09", "2009-03-08", 
										"2010-03-14", "2011-03-13", "2012-03-11", "2013-03-10", 
										"2014-03-09", "2015-03-08", "2016-03-13")), 
					 "HOUR" = rep(2, 11))

forward

# lets update the PHI data to not have any of the date-hours in the dataframe forward
# observations that have such a date-hour (2am) will be moved up to the next hour (3am)

PHI[, Dispatch_Date := as.Date(Dispatch_Date)]
PHI[, Hour := as.numeric(Hour)]
PHI_forward = data.table(PHI[Dispatch_Date %in% forward[,1] & Hour %in% forward[,2]])
PHI_forward

PHI[Dispatch_Date %in% forward$DATE & Hour %in% forward$HOUR, Hour := Hour + 1]

# verify the change

PHI[Dispatch_Date %in% forward[,1] & Hour %in% forward[,2]]

rm(daylightsavings, forward)

# lets transform the dispatch_date_time column to a POSIXct object
# we are just concerned with hourly dispatches, so lets just combine the dispatch_date and hour columns together to recreate dispatch_date_time

PHI[,Dispatch_Date_Time := as.POSIXct(paste0(PHI[,Dispatch_Date], 
											 if_else(PHI[,Hour] < 10, " 0", " "), 
											 PHI[,Hour], 
											 ":00:01"))]

# lets add a day column, week column, and year column

PHI[, Day := factor(strftime(PHI[,Dispatch_Date_Time], format = "%A"), 
					levels = c("Sunday", "Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday"))]

PHI[, Week := factor(strftime(PHI[,Dispatch_Date_Time], format = "%W"), 
					 levels = c("00", "01", "02", "03", "04", "05", "06", "07", "08", "09", as.character(10:53)))]

PHI[, Year := factor(strftime(PHI[,Dispatch_Date_Time], format = "%Y"), 
					 levels = as.character(2006:2016))]
					 
# change the following columns into proper data types

PHI[,Dc_Dist := as.factor(Dc_Dist)]
PHI[,Hour := as.factor(Hour)]
PHI[,UCR_General := as.factor(UCR_General)]
PHI[,Police_Districts := as.factor(Police_Districts)]
PHI[, Month := factor(strftime(PHI[,Dispatch_Date_Time], format = "%b"), 
					  levels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))]

# remove the following columns

PHI[,c("Dispatch_Time", "Dc_Key") := NULL]

# reorder the data by dispatch_date_time

setorder(PHI, cols = "Dispatch_Date_Time")

# the last day in the data set isn't a full day so lets remove that date

PHI = data.table(PHI[Dispatch_Date < "2016-08-18"])

# lets update the PHI data by spliting its 1am data into two 1am sets for the days corresponding to daylight savings moving backward

backward = data.frame("DATETIME" = as.POSIXct(c("2006-10-29 01:00:01", "2007-11-04 01:00:01",
												"2008-11-02 01:00:01", "2009-11-01 01:00:01",
												"2010-11-07 01:00:01", "2011-11-06 01:00:01",
												"2012-11-04 01:00:01", "2013-11-03 01:00:01",
												"2014-11-02 01:00:01", "2015-11-01 01:00:01")))

PHI = data.frame(PHI)												

x = lapply(1:NROW(backward), function(i) PHI[which(PHI[,3] == backward[i,1]),])
y = lapply(1:NROW(backward), function(i) x[[i]][1:round(NROW(x[[i]]) / 2, 0),])
z = lapply(1:NROW(backward), function(i) x[[i]][(round(NROW(x[[i]]) / 2, 0) + 1):NROW(x[[i]]),])
z = lapply(1:NROW(backward), function(i) data.frame(z[[i]][,1:2], "Dispatch_Date_Time" = z[[i]][,3] + 3601, z[[i]][,4:NCOL(z[[i]])]))
a = lapply(1:NROW(backward), function(i) rbind(y[[i]], z[[i]]))
b = rbind.fill(a)

PHI_backward = PHI[which(PHI[,3] %in% backward[,1]),]
PHI_backward

PHI[which(PHI[,3] %in% backward[,1]),] = b

rm(backward, x, y, z, a, b)

PHI = data.table(PHI)

# lets put these results in lists to organize our workspace

cleaning = list("na" = na, "PHI_backward" = PHI_backward, "PHI_forward" = PHI_forward)

rm(na, PHI_backward, PHI_forward)

# lets look at all of the levels for each of the factor variables
	# we won't look at Location_Block becuase this has too many levels to evaluate

levels(PHI[,Dc_Dist])
levels(PHI[,Psa])
levels(PHI[,Hour])
levels(PHI[,UCR_General])
levels(PHI[,Text_General_Code])
levels(PHI[,Police_Districts])
levels(PHI[,Month])
levels(PHI[,Day])
levels(PHI[,Week])
levels(PHI[,Year])

# there is a "" level in the factor Text_General_Code, lets remove that

PHI[,Text_General_Code := droplevels(Text_General_Code)]

}

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ---- Heat Maps ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

{

# PHI is a dataframe
# Lat and Lon are columns of the data.frame that contain longitude and latitude numbers

phillymap = get_map(location = c(lon = median(PHI$Lon), lat = median(PHI$Lat)), zoom = 11, color = "bw", maptype = "roadmap")

# 2006 - 2016

DF = data.frame(PHI)

heat = ggmap(phillymap, extent = "device") + 
	   stat_bin2d(data = DF, aes(x = Lon, y = Lat), alpha = 2/3, bins = 40) +
	   scale_fill_gradientn(colours = colorRampPalette(brewer.pal(n = 9, name = "YlOrRd"))(3), labels = c("5000", "15000", "25000"), breaks = c(5000, 15000, 25000)) +
	   scale_x_continuous(limits = c(min(PHI$Lon), max(PHI$Lon))) +
	   scale_y_continuous(limits = c(min(PHI$Lat), max(PHI$Lat))) +
	   labs(x = "Longitude", y = "Latitude", fill = "Crimes") +
	   ggtitle("Heat Map of Crimes\n2006 - 2016") + 
	   theme_classic(base_size = 20) +
	   theme(legend.key.size = unit(2/5, "in"))

heat
  
# 2006

DF = data.frame(PHI[Year == "2006"])

heat06 = ggmap(phillymap, extent = "device") + 
		 stat_bin2d(data = DF, aes(x = Lon, y = Lat), alpha = 2/3, bins = 40) +
		 scale_fill_gradientn(colours = colorRampPalette(brewer.pal(n = 9, name = "YlOrRd"))(3), labels = c("250", "1000", "2000"), breaks = c(250, 1000, 2000)) +
		 scale_x_continuous(limits = c(min(PHI$Lon), max(PHI$Lon))) +
		 scale_y_continuous(limits = c(min(PHI$Lat), max(PHI$Lat))) +
		 labs(x = "Longitude", y = "Latitude", fill = "Crimes") +
		 ggtitle("Heat Map of Crimes\n2006") + 
		 theme_classic(base_size = 20) +
		 theme(legend.key.size = unit(2/5, "in"))

heat06

# 2007

DF = data.frame(PHI[Year == "2007"])

heat07 = ggmap(phillymap, extent = "device") + 
		 stat_bin2d(data = DF, aes(x = Lon, y = Lat), alpha = 2/3, bins = 40) +
		 scale_fill_gradientn(colours = colorRampPalette(brewer.pal(n = 9, name = "YlOrRd"))(3), labels = c("250", "1000", "2000"), breaks = c(250, 1000, 2000)) +
		 scale_x_continuous(limits = c(min(PHI$Lon), max(PHI$Lon))) +
		 scale_y_continuous(limits = c(min(PHI$Lat), max(PHI$Lat))) +
		 labs(x = "Longitude", y = "Latitude", fill = "Crimes") +
		 ggtitle("Heat Map of Crimes\n2007") + 
		 theme_classic(base_size = 20) +
		 theme(legend.key.size = unit(2/5, "in"))

heat07

# 2008

DF = data.frame(PHI[Year == "2008"])

heat08 = ggmap(phillymap, extent = "device") + 
		 stat_bin2d(data = DF, aes(x = Lon, y = Lat), alpha = 2/3, bins = 40) +
		 scale_fill_gradientn(colours = colorRampPalette(brewer.pal(n = 9, name = "YlOrRd"))(3), labels = c("250", "1000", "2000"), breaks = c(250, 1000, 2000)) +
		 scale_x_continuous(limits = c(min(PHI$Lon), max(PHI$Lon))) +
		 scale_y_continuous(limits = c(min(PHI$Lat), max(PHI$Lat))) +
		 labs(x = "Longitude", y = "Latitude", fill = "Crimes") +
		 ggtitle("Heat Map of Crimes\n2008") + 
		 theme_classic(base_size = 20) +
		 theme(legend.key.size = unit(2/5, "in"))

heat08

# 2009

DF = data.frame(PHI[Year == "2009"])

heat09 = ggmap(phillymap, extent = "device") + 
		 stat_bin2d(data = DF, aes(x = Lon, y = Lat), alpha = 2/3, bins = 40) +
		 scale_fill_gradientn(colours = colorRampPalette(brewer.pal(n = 9, name = "YlOrRd"))(3), labels = c("250", "1000", "2000"), breaks = c(250, 1000, 2000)) +
		 scale_x_continuous(limits = c(min(PHI$Lon), max(PHI$Lon))) +
		 scale_y_continuous(limits = c(min(PHI$Lat), max(PHI$Lat))) +
		 labs(x = "Longitude", y = "Latitude", fill = "Crimes") +
		 ggtitle("Heat Map of Crimes\n2009") + 
		 theme_classic(base_size = 20) +
		 theme(legend.key.size = unit(2/5, "in"))

heat09

# 2010

DF = data.frame(PHI[Year == "2010"])

heat10 = ggmap(phillymap, extent = "device") + 
		 stat_bin2d(data = DF, aes(x = Lon, y = Lat), alpha = 2/3, bins = 40) +
		 scale_fill_gradientn(colours = colorRampPalette(brewer.pal(n = 9, name = "YlOrRd"))(3), labels = c("250", "1000", "2000"), breaks = c(250, 1000, 2000)) +
		 scale_x_continuous(limits = c(min(PHI$Lon), max(PHI$Lon))) +
		 scale_y_continuous(limits = c(min(PHI$Lat), max(PHI$Lat))) +
		 labs(x = "Longitude", y = "Latitude", fill = "Crimes") +
		 ggtitle("Heat Map of Crimes\n2010") + 
		 theme_classic(base_size = 20) +
		 theme(legend.key.size = unit(2/5, "in"))

heat10

# 2011

DF = data.frame(PHI[Year == "2011"])

heat11 = ggmap(phillymap, extent = "device") + 
		 stat_bin2d(data = DF, aes(x = Lon, y = Lat), alpha = 2/3, bins = 40) +
		 scale_fill_gradientn(colours = colorRampPalette(brewer.pal(n = 9, name = "YlOrRd"))(3), labels = c("250", "1000", "2000"), breaks = c(250, 1000, 2000)) +
		 scale_x_continuous(limits = c(min(PHI$Lon), max(PHI$Lon))) +
		 scale_y_continuous(limits = c(min(PHI$Lat), max(PHI$Lat))) +
		 labs(x = "Longitude", y = "Latitude", fill = "Crimes") +
		 ggtitle("Heat Map of Crimes\n2011") + 
		 theme_classic(base_size = 20) +
		 theme(legend.key.size = unit(2/5, "in"))

heat11

# 2012

DF = data.frame(PHI[Year == "2012"])

heat12 = ggmap(phillymap, extent = "device") + 
		 stat_bin2d(data = DF, aes(x = Lon, y = Lat), alpha = 2/3, bins = 40) +
		 scale_fill_gradientn(colours = colorRampPalette(brewer.pal(n = 9, name = "YlOrRd"))(3), labels = c("250", "1000", "2000"), breaks = c(250, 1000, 2000)) +
		 scale_x_continuous(limits = c(min(PHI$Lon), max(PHI$Lon))) +
		 scale_y_continuous(limits = c(min(PHI$Lat), max(PHI$Lat))) +
		 labs(x = "Longitude", y = "Latitude", fill = "Crimes") +
		 ggtitle("Heat Map of Crimes\n2012") + 
		 theme_classic(base_size = 20) +
		 theme(legend.key.size = unit(2/5, "in"))

heat12

# 2013

DF = data.frame(PHI[Year == "2013"])

heat13 = ggmap(phillymap, extent = "device") + 
		 stat_bin2d(data = DF, aes(x = Lon, y = Lat), alpha = 2/3, bins = 40) +
		 scale_fill_gradientn(colours = colorRampPalette(brewer.pal(n = 9, name = "YlOrRd"))(3), labels = c("250", "1000", "2000"), breaks = c(250, 1000, 2000)) +
		 scale_x_continuous(limits = c(min(PHI$Lon), max(PHI$Lon))) +
		 scale_y_continuous(limits = c(min(PHI$Lat), max(PHI$Lat))) +
		 labs(x = "Longitude", y = "Latitude", fill = "Crimes") +
		 ggtitle("Heat Map of Crimes\n2013") + 
		 theme_classic(base_size = 20) +
		 theme(legend.key.size = unit(2/5, "in"))

heat13

# 2014

DF = data.frame(PHI[Year == "2014"])

heat14 = ggmap(phillymap, extent = "device") + 
		 stat_bin2d(data = DF, aes(x = Lon, y = Lat), alpha = 2/3, bins = 40) +
		 scale_fill_gradientn(colours = colorRampPalette(brewer.pal(n = 9, name = "YlOrRd"))(3), labels = c("250", "1000", "2000"), breaks = c(250, 1000, 2000)) +
		 scale_x_continuous(limits = c(min(PHI$Lon), max(PHI$Lon))) +
		 scale_y_continuous(limits = c(min(PHI$Lat), max(PHI$Lat))) +
		 labs(x = "Longitude", y = "Latitude", fill = "Crimes") +
		 ggtitle("Heat Map of Crimes\n2014") + 
		 theme_classic(base_size = 20) +
		 theme(legend.key.size = unit(2/5, "in"))

heat14

# 2015

DF = data.frame(PHI[Year == "2015"])

heat15 = ggmap(phillymap, extent = "device") + 
		 stat_bin2d(data = DF, aes(x = Lon, y = Lat), alpha = 2/3, bins = 40) +
		 scale_fill_gradientn(colours = colorRampPalette(brewer.pal(n = 9, name = "YlOrRd"))(3), labels = c("250", "1000", "2000"), breaks = c(250, 1000, 2000)) +
		 scale_x_continuous(limits = c(min(PHI$Lon), max(PHI$Lon))) +
		 scale_y_continuous(limits = c(min(PHI$Lat), max(PHI$Lat))) +
		 labs(x = "Longitude", y = "Latitude", fill = "Crimes") +
		 ggtitle("Heat Map of Crimes\n2015") + 
		 theme_classic(base_size = 20) +
		 theme(legend.key.size = unit(2/5, "in"))

heat15

# define area1

cutA1 = heat + 
		geom_vline(xintercept = -75.1733, color = "blue") + 
		geom_vline(xintercept = -75.1562, color = "blue") + 
		geom_hline(yintercept = 39.9475, color = "blue") + 
		geom_hline(yintercept = 39.9553, color = "blue")

cutA1

# define area2

cutA2 = heat + 
		geom_vline(xintercept = -75.1329, color = "blue") + 
		geom_vline(xintercept = -75.1083, color = "blue") + 
		geom_hline(yintercept = 39.9875, color = "blue") + 
		geom_hline(yintercept = 40.0008, color = "blue")

cutA2

# define West Philly 1

cutWP1 = heat +  
		 geom_vline(xintercept = -75.2532, color = "blue") + 
		 geom_vline(xintercept = -75.2129, color = "blue") + 
		 geom_hline(yintercept = 39.9152, color = "blue") + 
		 geom_hline(yintercept = 39.9475, color = "blue")

cutWP1

# define West Philly 2

cutWP2 = heat + 
		 geom_vline(xintercept = -75.2532, color = "blue") + 
		 geom_vline(xintercept = -75.1877, color = "blue") +
		 geom_hline(yintercept = 39.9819, color = "blue") +
		 geom_hline(yintercept = 39.9475, color = "blue")

cutWP2

# define Central Philly 1

cutCP1 = heat + 
		 geom_vline(xintercept = -75.1970, color = "blue") + 
		 geom_vline(xintercept = -75.1322, color = "blue") + 
		 geom_hline(yintercept = 39.9152, color = "blue") + 
		 geom_hline(yintercept = 39.9418, color = "blue")

cutCP1

# define Central Philly 2

cutCP2 = heat + 
		 geom_vline(xintercept = -75.1877, color = "blue") +
		 geom_vline(xintercept = -75.1322, color = "blue") + 
		 geom_hline(yintercept = 39.9677, color = "blue") +
		 geom_hline(yintercept = 39.9418, color = "blue")

cutCP2

# define North Philly 

cutNP = heat + 
		geom_vline(xintercept = -75.1877, color = "blue") +
		geom_vline(xintercept = -75.0352, color = "blue") + 
		geom_hline(yintercept = 39.9677, color = "blue") +
		geom_hline(yintercept = 40.0470, color = "blue")

cutNP

# create an Area1 column for West Philly, North Philly, and Central Philly

PHI[, Area1 := factor(if_else(-75.1970 <= PHI[,Lon] & PHI[,Lon] <= -75.1322 & 
							  39.9152 <= PHI[,Lat] & PHI[,Lat] <= 39.9418,
						      "CP",
				      if_else(-75.1877 <= PHI[,Lon] & PHI[,Lon] <= -75.1322 & 
						      39.9418 < PHI[,Lat] & PHI[,Lat] <= 39.9677,
						      "CP",
				      if_else(-75.2532 <= PHI[,Lon] & PHI[,Lon] <= -75.2129 & 
						      39.9152 <= PHI[,Lat] & PHI[,Lat] <= 39.9475,
						      "WP",
				      if_else(-75.2532 <= PHI[,Lon] & PHI[,Lon] < -75.1877 & 
						      39.9475 < PHI[,Lat] & PHI[,Lat] <= 39.9819,
						      "WP",
				      if_else(-75.1877 <= PHI[,Lon] & PHI[,Lon] <= -75.0352 & 
						      39.9677 < PHI[,Lat] & PHI[,Lat] <= 40.0470,
						      "NP", "Other"))))),
				      levels = c("CP", "WP", "NP", "Other"))]

# verify the change worked

DF = data.frame(PHI[Area1 == "CP"])

heatCP = ggmap(phillymap, extent = "device") + 
		 stat_bin2d(data = DF, aes(x = Lon, y = Lat), alpha = 2/3, bins = 40) +
		 scale_fill_gradientn(colours = colorRampPalette(brewer.pal(n = 9, name = "YlOrRd"))(3), labels = c("5000", "15000", "25000"), breaks = c(5000, 15000, 25000)) +
		 scale_x_continuous(limits = c(min(PHI$Lon), max(PHI$Lon))) +
		 scale_y_continuous(limits = c(min(PHI$Lat), max(PHI$Lat))) +
		 labs(x = "Longitude", y = "Latitude", fill = "Crimes") +
		 ggtitle("Heat Map of Central Philly Crimes\n2006 - 2016") + 
		 theme_classic(base_size = 20) +
		 theme(legend.key.size = unit(2/5, "in"))

heatCP

DF = data.frame(PHI[Area1 == "WP"])

heatWP = ggmap(phillymap, extent = "device") + 
		 stat_bin2d(data = DF, aes(x = Lon, y = Lat), alpha = 2/3, bins = 40) +
		 scale_fill_gradientn(colours = colorRampPalette(brewer.pal(n = 9, name = "YlOrRd"))(3), labels = c("5000", "15000", "25000"), breaks = c(5000, 15000, 25000)) +
		 scale_x_continuous(limits = c(min(PHI$Lon), max(PHI$Lon))) +
		 scale_y_continuous(limits = c(min(PHI$Lat), max(PHI$Lat))) +
		 labs(x = "Longitude", y = "Latitude", fill = "Crimes") +
		 ggtitle("Heat Map of West Philly Crimes\n2006 - 2016") + 
		 theme_classic(base_size = 20) +
		 theme(legend.key.size = unit(2/5, "in"))

heatWP

DF = data.frame(PHI[Area1 == "NP"])

heatNP = ggmap(phillymap, extent = "device") + 
		 stat_bin2d(data = DF, aes(x = Lon, y = Lat), alpha = 2/3, bins = 40) +
		 scale_fill_gradientn(colours = colorRampPalette(brewer.pal(n = 9, name = "YlOrRd"))(3), labels = c("5000", "15000", "25000"), breaks = c(5000, 15000, 25000)) +
		 scale_x_continuous(limits = c(min(PHI$Lon), max(PHI$Lon))) +
		 scale_y_continuous(limits = c(min(PHI$Lat), max(PHI$Lat))) +
		 labs(x = "Longitude", y = "Latitude", fill = "Crimes") +
		 ggtitle("Heat Map of North Philly Crimes\n2006 - 2016") + 
		 theme_classic(base_size = 20) +
		 theme(legend.key.size = unit(2/5, "in"))

heatNP

# create an Area2 column for A1 and A2

PHI[, Area2 := factor(if_else(-75.1733 <= PHI[,Lon] & PHI[,Lon] <= -75.1562 & 
						      39.9475 <= PHI[,Lat] & PHI[,Lat] <= 39.9553,
						      "A1",
				      if_else(-75.1329 <= PHI[,Lon] & PHI[,Lon] <= -75.1083 & 
						      39.9875 <= PHI[,Lat] & PHI[,Lat] <= 40.0008,
						      "A2", "Other")),
				      levels = c("A1", "A2", "Other"))]

# verify the change worked

DF = data.frame(PHI[Area2 == "A1"])

heatA1 = ggmap(phillymap, extent = "device") + 
		 stat_bin2d(data = DF, aes(x = Lon, y = Lat), alpha = 2/3, bins = 40) +
		 scale_fill_gradientn(colours = colorRampPalette(brewer.pal(n = 9, name = "YlOrRd"))(3), labels = c("5000", "15000", "25000"), breaks = c(5000, 15000, 25000)) +
		 scale_x_continuous(limits = c(min(PHI$Lon), max(PHI$Lon))) +
		 scale_y_continuous(limits = c(min(PHI$Lat), max(PHI$Lat))) +
		 labs(x = "Longitude", y = "Latitude", fill = "Crimes") +
		 ggtitle("Heat Map of A1 Crimes\n2006 - 2016") + 
		 theme_classic(base_size = 20) +
		 theme(legend.key.size = unit(2/5, "in"))

heatA1

DF = data.frame(PHI[Area2 == "A2"])

heatA2 = ggmap(phillymap, extent = "device") + 
		 stat_bin2d(data = DF, aes(x = Lon, y = Lat), alpha = 2/3, bins = 40) +
		 scale_fill_gradientn(colours = colorRampPalette(brewer.pal(n = 9, name = "YlOrRd"))(3), labels = c("5000", "15000", "25000"), breaks = c(5000, 15000, 25000)) +
		 scale_x_continuous(limits = c(min(PHI$Lon), max(PHI$Lon))) +
		 scale_y_continuous(limits = c(min(PHI$Lat), max(PHI$Lat))) +
		 labs(x = "Longitude", y = "Latitude", fill = "Crimes") +
		 ggtitle("Heat Map of A2 Crimes\n2006 - 2016") + 
		 theme_classic(base_size = 20) +
		 theme(legend.key.size = unit(2/5, "in"))

heatA2

# create a list for these heat maps 

heatmaps = list("heat" = heat,
				"heat06" = heat06,
				"heat07" = heat07,
				"heat08" = heat08,
				"heat09" = heat09,
				"heat10" = heat10,
				"heat11" = heat11,
				"heat12" = heat12,
				"heat13" = heat13,
				"heat14" = heat14,
				"heat15" = heat15,
				"cutA1" = cutA1,
				"cutA2" = cutA2,
				"cutCP1" = cutCP1,
				"cutCP2" = cutCP2,
				"cutNP" = cutNP,
				"cutWP1" = cutWP1,
				"cutWP2" = cutWP2,
				"heatA1" = heatA1,
				"heatA2" = heatA2,
				"heatCP" = heatCP,
				"heatNP" = heatNP,
				"heatWP" = heatWP)

rm(heat, 
   heat06, 
   heat07, 
   heat08, 
   heat09, 
   heat10, 
   heat11, 
   heat12, 
   heat13, 
   heat14, 
   heat15, 
   cutA1, 
   cutA2, 
   cutCP1, 
   cutCP2, 
   cutNP, 
   cutWP1, 
   cutWP2,
   heatA1, 
   heatA2, 
   heatCP, 
   heatNP, 
   heatWP,
   DF, 
   phillymap)

}

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ---- Crimes -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

{

# lets create a table of crimes per hour ------------------------------------------------------------------------------------------------------------------------------------

crimes_hourly = PHI[,.("CRIMES" = .N), by = .("DATETIME" = Dispatch_Date_Time)]

# lets make sure that there is one hourly datetime of data everyday between the min and max datetimes

datetimes = data.frame("DATETIME" = seq.POSIXt(from = min(crimes_hourly[,DATETIME]), 
											   to = max(crimes_hourly[,DATETIME]), 
											   by = "hour"))

# lets update datetimes to correspond to the two sets of 1am's for days that have daylight savings moving backward

backward = data.frame("DATETIME" = as.POSIXct(c("2006-10-29 01:00:01", "2007-11-04 01:00:01",
												"2008-11-02 01:00:01", "2009-11-01 01:00:01",
												"2010-11-07 01:00:01", "2011-11-06 01:00:01",
												"2012-11-04 01:00:01", "2013-11-03 01:00:01",
												"2014-11-02 01:00:01", "2015-11-01 01:00:01")))

x = which(datetimes[,1] %in% backward[,1]) + 1
datetimes[x,] = datetimes[x,] + 1
datetimes = data.table(datetimes)

rm(backward, x)

# lets do a right join to find out if there are missing datetimes

setkey(datetimes, DATETIME)
setkey(crimes_hourly, DATETIME)
crimes_hourly = data.table(crimes_hourly[datetimes])
increment1 = na.omit(crimes_hourly, cols = "CRIMES", invert = TRUE)
increment1

# based on the missing datetimes lets replicate the data in PHI for the previous week corresponding to the NA's in x, and then use this replicated data to fill in these NA's in x

dup = na.omit(crimes_hourly, cols = "CRIMES", invert = TRUE)[,DATETIME] - (3600 * 24 * 7)
dup = PHI[Dispatch_Date_Time %in% dup]
dup[, Dispatch_Date_Time := Dispatch_Date_Time + (3600 * 24 * 7)]

dup[, Dispatch_Date := as.Date(strftime(dup[,Dispatch_Date_Time], format = "%Y-%m-%d"))]

dup[, Day := factor(strftime(dup[,Dispatch_Date_Time], format = "%A"), 
					levels = c("Sunday", "Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday"))]

dup[, Week := factor(strftime(dup[,Dispatch_Date_Time], format = "%W"), 
					 levels = c("00", "01", "02", "03", "04", "05", "06", "07", "08", "09", as.character(10:53)))]

dup[, Month := factor(strftime(dup[,Dispatch_Date_Time], format = "%b"), 
					  levels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))]
				  
dup[, Year := factor(strftime(dup[,Dispatch_Date_Time], format = "%Y"), 
					 levels = as.character(2006:2016))]
					 
PHI = rbindlist(list(PHI, dup))
setorder(PHI, cols = "Dispatch_Date_Time")

# lets re-create crimes_hourly and see if all NA's have been removed

crimes_hourly = PHI[,.("CRIMES" = .N), by = .("DATETIME" = Dispatch_Date_Time)]

setkey(crimes_hourly, DATETIME)
crimes_hourly = data.table(crimes_hourly[datetimes])
increment2 = na.omit(crimes_hourly, cols = "CRIMES", invert = TRUE)
increment2

# based on the missing datetimes lets replicate the data in PHI for the previous week corresponding to the NA's in x, and then use this replicated data to fill in these NA's in x

dup = na.omit(crimes_hourly, cols = "CRIMES", invert = TRUE)[,DATETIME] - (3600 * 24 * 7)
dup = PHI[Dispatch_Date_Time %in% dup]
dup[, Dispatch_Date_Time := Dispatch_Date_Time + (3600 * 24 * 7)]

dup[, Dispatch_Date := as.Date(strftime(dup[,Dispatch_Date_Time], format = "%Y-%m-%d"))]

dup[, Day := factor(strftime(dup[,Dispatch_Date_Time], format = "%A"), 
					levels = c("Sunday", "Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday"))]

dup[, Week := factor(strftime(dup[,Dispatch_Date_Time], format = "%W"), 
					 levels = c("00", "01", "02", "03", "04", "05", "06", "07", "08", "09", as.character(10:53)))]

dup[, Month := factor(strftime(dup[,Dispatch_Date_Time], format = "%b"), 
					  levels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))]

dup[, Year := factor(strftime(dup[,Dispatch_Date_Time], format = "%Y"), 
					 levels = as.character(2006:2016))]
					 
PHI = rbindlist(list(PHI, dup))
setorder(PHI, cols = "Dispatch_Date_Time")

# lets re-create crimes_hourly and see if all NA's have been removed

crimes_hourly = PHI[,.("CRIMES" = .N), by = .("DATETIME" = Dispatch_Date_Time)]

setkey(crimes_hourly, DATETIME)
crimes_hourly = data.table(crimes_hourly[datetimes])
na.omit(crimes_hourly, cols = "CRIMES", invert = TRUE)

# all NA's have been removed and there is one hourly datetime from the min to max datetime in PHI

increments = list("increment1" = increment1, "increment2" = increment2)
cleaning$increments = increments

rm(dup, datetimes, increment1, increment2, increments)

# lets add a hour column to crimes_hourly for graphing purposes

crimes_hourly = PHI[,.("CRIMES" = .N), by = .("DATETIME" = Dispatch_Date_Time, "HOUR" = Hour)]

# lets graph crimes_hourly

DF = data.frame(crimes_hourly)

plot_crimes_hourly = ggplot(data = DF, aes(x = DATETIME, y = CRIMES)) +
					 geom_point(na.rm = TRUE) +
				   	 labs(x = "Hour", y = "Crimes") +
				   	 ggtitle("Crimes Per Hour") +
				   	 theme_light(base_size = 20)

plot_crimes_hourly

plot_crimes_hourly_outlier = plot_crimes_hourly
cleaning$plot_crimes_hourly_outlier = plot_crimes_hourly_outlier

rm(plot_crimes_hourly_outlier)

# lets check out that outlier up top

crimes_hourly[CRIMES == max(CRIMES)][,DATETIME]
print(PHI[Dispatch_Date_Time == crimes_hourly[CRIMES == max(CRIMES)][,DATETIME]], nrow = 177)

# there is alot of non-residential burgleries at 7000 BLOCK LINDBERGH BLVD
# this is a public storage center, perhaps alot of units were stolen from
# but lets remove some of these observations from the data becuase they are anomolies

# this is the average number of crimes at 8am for the whole dataset

newpoint = round(crimes_hourly[HOUR == "8", mean(CRIMES)], 0)
newpoint

# out of all 7000 BLOCK LINDBERGH BLVD crimes at this datetime, only keep newpoint-many

dat1 = data.table(PHI[Dispatch_Date_Time == crimes_hourly[CRIMES == max(CRIMES)][,DATETIME] & Location_Block == "7000 BLOCK LINDBERGH BLVD"][1:newpoint,])

# keep all NON-7000 BLOCK LINDBERGH BLVD crimes at this datetime

dat2 = data.table(PHI[Dispatch_Date_Time == crimes_hourly[CRIMES == max(CRIMES)][,DATETIME] & Location_Block != "7000 BLOCK LINDBERGH BLVD"])

# keep all crimes not at this datetime

PHI = data.table(PHI[Dispatch_Date_Time != crimes_hourly[CRIMES == max(CRIMES)][,DATETIME]])

# combine new data

PHI = rbindlist(list(PHI, dat1, dat2))
setorder(PHI, cols = "Dispatch_Date_Time")

rm(newpoint, dat1, dat2)

# lets rebuild crimes_hourly and plot_crimes_hourly

crimes_hourly = PHI[,.("CRIMES" = .N), by = .("DATETIME" = Dispatch_Date_Time, "HOUR" = Hour)]

DF = data.frame(crimes_hourly)

plot_crimes_hourly = ggplot(data = DF, aes(x = DATETIME, y = CRIMES)) +
					 geom_point(na.rm = TRUE) +
				   	 labs(x = "Hour", y = "Crimes") +
				   	 ggtitle("Crimes Per Hour") +
				   	 theme_light(base_size = 20)

plot_crimes_hourly

# lets zoom in on one week to see how crimes per hour behave

DF = data.frame(crimes_hourly[1:168,])

plot_crimes_hourly_zoom = ggplot(data = DF, aes(x = DATETIME, y = CRIMES)) +
						  geom_point(na.rm = TRUE) +
						  labs(x = "Hour", y = "Crimes") +
						  ggtitle("Crimes Per Hour\n2006-01-01 00:00:01 EST - 2006-01-07 23:00:01 EST") +
						  theme_light(base_size = 20)
						  
plot_crimes_hourly_zoom

plot_crimes_hourly_zoom_color = ggplot(data = DF, aes(x = DATETIME, y = CRIMES, color = HOUR)) +
								geom_point(na.rm = TRUE) +
								scale_color_manual(values = colorRampPalette(brewer.pal(n = 9, name = "YlOrRd"))(24)) +
								labs(x = "Hour", y = "Crimes") +
								ggtitle("Crimes Per Hour\n2006-01-01 00:00:01 EST - 2006-01-07 23:00:01 EST") +
								theme_dark(base_size = 20) + 
								theme(legend.position = "top",
									  legend.key.size = unit(.5, "in")) + 
								guides(color = guide_legend(nrow = 2, override.aes = list(size = 10)))

plot_crimes_hourly_zoom_color

# lets put these results in lists to organize our workspace

crime_tables = list("crimes_hourly" = crimes_hourly)
crime_plots = list("plot_crimes_hourly" = plot_crimes_hourly, 
				   "plot_crimes_hourly_zoom" = plot_crimes_hourly_zoom, 
				   "plot_crimes_hourly_zoom_color" = plot_crimes_hourly_zoom_color)

rm(crimes_hourly, plot_crimes_hourly, plot_crimes_hourly_zoom, plot_crimes_hourly_zoom_color)

# lets create a table of crimes per day ------------------------------------------------------------------------------------------------------------------------------------

crimes_daily = PHI[,.("CRIMES" = .N), by = .("DATE" = Dispatch_Date)]

# lets plot crimes_daily

DF = data.frame(crimes_daily)

plot_crimes_daily = ggplot(data = DF, aes(x = DATE, y = CRIMES)) +
					geom_point(na.rm = TRUE) +
				   	labs(x = "Day", y = "Crimes") +
				   	ggtitle("Crimes Per Day") +
				   	theme_light(base_size = 20)

plot_crimes_daily

# lets add a month column to crimes_daily and plot it

crimes_daily = PHI[,.("CRIMES" = .N), by = .("DATE" = Dispatch_Date, "MONTH" = Month)]

DF = data.frame(crimes_daily)

plot_crimes_daily_color = ggplot(data = DF, aes(x = DATE, y = CRIMES, color = MONTH)) +
						  geom_point(na.rm = TRUE) +
				   		  scale_color_manual(values = colorRampPalette(brewer.pal(n = 9, name = "YlOrRd"))(12)) +		
						  labs(x = "Day", y = "Crimes") +
				   		  ggtitle("Crimes Per Day") +
				   		  theme_dark(base_size = 20) + 
						  theme(legend.position = "top",
								legend.key.size = unit(.5, "in")) + 
						  guides(color = guide_legend(nrow = 2, override.aes = list(size = 10)))

plot_crimes_daily_color

# lets zoom in on three months to see how crimes per day behave

DF = data.frame(crimes_daily[1:90,])

plot_crimes_daily_zoom = ggplot(data = DF, aes(x = DATE, y = CRIMES)) +
						  geom_point(na.rm = TRUE) +
						  scale_x_date(limits = c(DF$DATE[1], DF$DATE[90])) +
						  labs(x = "Day", y = "Crimes") +
						  ggtitle("Crimes Per Day\n2006-01-01 - 2006-03-31") +
						  theme_light(base_size = 20)
						  
plot_crimes_daily_zoom

# lets add a week column to crimes_daily, zoom in to three months, and plot it

crimes_daily = PHI[,.("CRIMES" = .N), by = .("DATE" = Dispatch_Date, "MONTH" = Month, "WEEK" = Week)]

DF = data.frame(crimes_daily[1:90,])

plot_crimes_daily_zoom_color = ggplot(data = DF, aes(x = DATE, y = CRIMES, color = WEEK)) +
								geom_point(na.rm = TRUE) +
								scale_color_manual(values = colorRampPalette(brewer.pal(n = 12, name = "Paired"))(14)) +		
								labs(x = "Day", y = "Crimes") +
								ggtitle("Crimes Per Day\n2006-01-01 - 2006-03-31") +
								theme_gray(base_size = 20) + 
								theme(legend.position = "top",
									  legend.key.size = unit(.5, "in")) + 
								guides(color = guide_legend(nrow = 1, override.aes = list(size = 10)))
						  
plot_crimes_daily_zoom_color

# lets put these results in lists to organize our workspace

crime_tables$crimes_daily = crimes_daily
crime_plots$plot_crimes_daily = plot_crimes_daily
crime_plots$plot_crimes_daily_color = plot_crimes_daily_color
crime_plots$plot_crimes_daily_zoom = plot_crimes_daily_zoom
crime_plots$plot_crimes_daily_zoom_color = plot_crimes_daily_zoom_color

rm(crimes_daily, plot_crimes_daily, plot_crimes_daily_color, plot_crimes_daily_zoom, plot_crimes_daily_zoom_color)

# lets create a table of crimes per week ------------------------------------------------------------------------------------------------------------------------------------

crimes_weekly = PHI[,.("CRIMES" = .N), by = .("YEAR" = Year, "WEEK" = Week)]
crimes_weekly[, YEAR.WEEK := factor(paste0(crimes_weekly[,YEAR], ".", crimes_weekly[,WEEK]),
							 levels = paste0(crimes_weekly[,YEAR], ".", crimes_weekly[,WEEK]))]
crimes_weekly[, YEAR := NULL]
setcolorder(crimes_weekly, c("YEAR.WEEK", "WEEK", "CRIMES"))

# lets plot crimes_weekly

DF = data.frame(crimes_weekly)

plot_crimes_weekly = ggplot(data = DF, aes(x = YEAR.WEEK, y = CRIMES)) +
					 geom_point(na.rm = TRUE) +
					 scale_x_discrete(breaks = c("2006.01", "2007.01", "2008.01", "2009.01", "2010.01", "2011.01", "2012.01", "2013.01", "2014.01", "2015.01", "2016.01")) +
				   	 labs(x = "Year.Week", y = "Crimes") +
				   	 ggtitle("Crimes Per Week") +
				   	 theme_light(base_size = 20) + 
					 theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_crimes_weekly

# lets color crime_weekly by week

plot_crimes_weekly_color = ggplot(data = DF, aes(x = YEAR.WEEK, y = CRIMES, color = WEEK)) +
						   geom_point(na.rm = TRUE) +
						   scale_x_discrete(breaks = c("2006.01", "2007.01", "2008.01", "2009.01", "2010.01", "2011.01", "2012.01", "2013.01", "2014.01", "2015.01", "2016.01")) +
				   		   scale_color_manual(values = colorRampPalette(brewer.pal(n = 9, name = "YlOrRd"))(54)) +		
						   labs(x = "Year.Week", y = "Crimes") +
				   		   ggtitle("Crimes Per Week") +
				   		   theme_dark(base_size = 20) + 
						   theme(legend.position = "top",
								 legend.key.size = unit(.5, "in"),
						   	   	 axis.text.x = element_text(angle = 45, hjust = 1)) + 
						   guides(color = guide_legend(nrow = 3, override.aes = list(size = 10)))

plot_crimes_weekly_color

# lets put these results in lists to organize our workspace

crime_tables$crimes_weekly = crimes_weekly
crime_plots$plot_crimes_weekly = plot_crimes_weekly
crime_plots$plot_crimes_weekly_color = plot_crimes_weekly_color

rm(crimes_weekly, plot_crimes_weekly, plot_crimes_weekly_color)

# lets create a table of crimes per month ------------------------------------------------------------------------------------------------------------------------------------

crimes_monthly = PHI[,.("CRIMES" = .N), by = .("YEAR" = Year, "MONTH" = Month)]
crimes_monthly[, YEAR.MONTH := factor(paste0(crimes_monthly[,YEAR], ".", crimes_monthly[,MONTH]),
							 levels = paste0(crimes_monthly[,YEAR], ".", crimes_monthly[,MONTH]))]
crimes_monthly[, YEAR := NULL]
setcolorder(crimes_monthly, c("YEAR.MONTH", "MONTH", "CRIMES"))

# lets plot crimes_monthly

DF = data.frame(crimes_monthly)

plot_crimes_monthly = ggplot(data = DF, aes(x = YEAR.MONTH, y = CRIMES)) +
					  geom_point(na.rm = TRUE) +
					  scale_x_discrete(breaks = c("2006.Jan", "2007.Jan", "2008.Jan", "2009.Jan", "2010.Jan", "2011.Jan", "2012.Jan", "2013.Jan", "2014.Jan", "2015.Jan", "2016.Jan")) +
				   	  labs(x = "Year.Month", y = "Crimes") +
				   	  ggtitle("Crimes Per Month") +
				   	  theme_light(base_size = 20) + 
					  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_crimes_monthly

# lets color crimes_monthly by month

plot_crimes_monthly_color = ggplot(data = DF, aes(x = YEAR.MONTH, y = CRIMES, color = MONTH)) +
						    geom_point(na.rm = TRUE) +
						    scale_x_discrete(breaks = c("2006.Jan", "2007.Jan", "2008.Jan", "2009.Jan", "2010.Jan", "2011.Jan", "2012.Jan", "2013.Jan", "2014.Jan", "2015.Jan", "2016.Jan")) +
				   		    scale_color_manual(values = colorRampPalette(brewer.pal(n = 9, name = "YlOrRd"))(12)) +		
						    labs(x = "Year.Month", y = "Crimes") +
				   		    ggtitle("Crimes Per Month") +
				   		    theme_dark(base_size = 20) + 
						    theme(legend.position = "top",
								  legend.key.size = unit(.5, "in"),
						   	   	  axis.text.x = element_text(angle = 45, hjust = 1)) + 
						    guides(color = guide_legend(nrow = 1, override.aes = list(size = 10)))

plot_crimes_monthly_color

# lets put these results in lists to organize our workspace

crime_tables$crimes_monthly = crimes_monthly
crime_plots$plot_crimes_monthly = plot_crimes_monthly
crime_plots$plot_crimes_monthly_color = plot_crimes_monthly_color

rm(crimes_monthly, plot_crimes_monthly, plot_crimes_monthly_color)

# lets create a table of crimes per day of week ------------------------------------------------------------------------------------------------------------------------------------

crimes_dow = PHI[,.("CRIMES" = .N), by = .("DATE" = Dispatch_Date, "DAY" = Day)]
crimes_dow[, CRIMES := as.numeric(CRIMES)]
crimes_dow = crimes_dow[, .("MEAN" = mean(CRIMES), "MEDIAN" = median(CRIMES), "SD" = sd(CRIMES)), 
						  by = .("DAY" = DAY)]

# lets plot crimes_dow

DF = data.frame(crimes_dow)

plot_crimes_dow = ggplot(data = DF, aes(x = DAY, y = MEAN)) +
				  geom_bar(stat = "identity") +
				  labs(x = "Day", y = "Crimes") +
				  ggtitle("Average Crimes Per Day of Week") +
				  theme_light(base_size = 20)

plot_crimes_dow

# lets put these results in lists to organize our workspace

crime_tables$crimes_dow = crimes_dow
crime_plots$plot_crimes_dow = plot_crimes_dow

rm(crimes_dow, plot_crimes_dow)

}

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ---- Crimes By Dc_Dist -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

{

# lets create a table of crimes per hour ------------------------------------------------------------------------------------------------------------------------------------

crimes_hourly = PHI[,.("CRIMES" = .N), by = .("DC" = Dc_Dist, "DATETIME" = Dispatch_Date_Time, "HOUR" = Hour)]

# lets plot crimes_hourly

DF = data.frame(crimes_hourly[DC %in% levels(DC)[1:6]])

plot_crimes_hourly_set1 = ggplot(data = DF, aes(x = DATETIME, y = CRIMES, color = DC)) +
						  geom_point(na.rm = TRUE) +
						  facet_wrap(~DC, scales = "free_y") +
						  labs(x = "Hour", y = "Crimes") +
						  ggtitle("Crimes Per Hour By DC") +
						  theme_light(base_size = 15) + 
						  theme(legend.position = "top",
								legend.key.size = unit(.25, "in")) + 
						  guides(color = guide_legend(nrow = 1, override.aes = list(size = 10)))

plot_crimes_hourly_set1

DF = data.frame(crimes_hourly[DC %in% levels(DC)[7:12]])

plot_crimes_hourly_set2 = ggplot(data = DF, aes(x = DATETIME, y = CRIMES, color = DC)) +
						  geom_point(na.rm = TRUE) +
						  facet_wrap(~DC, scales = "free_y") +
						  labs(x = "Hour", y = "Crimes") +
						  ggtitle("Crimes Per Hour By DC") +
						  theme_light(base_size = 15) + 
						  theme(legend.position = "top",
								legend.key.size = unit(.25, "in")) + 
						  guides(color = guide_legend(nrow = 1, override.aes = list(size = 10)))

plot_crimes_hourly_set2

DF = data.frame(crimes_hourly[DC %in% levels(DC)[13:18]])

plot_crimes_hourly_set3 = ggplot(data = DF, aes(x = DATETIME, y = CRIMES, color = DC)) +
						  geom_point(na.rm = TRUE) +
						  facet_wrap(~DC, scales = "free_y") +
						  labs(x = "Hour", y = "Crimes") +
						  ggtitle("Crimes Per Hour By DC") +
						  theme_light(base_size = 15) + 
						  theme(legend.position = "top",
								legend.key.size = unit(.25, "in")) + 
						  guides(color = guide_legend(nrow = 1, override.aes = list(size = 10)))

plot_crimes_hourly_set3

DF = data.frame(crimes_hourly[DC %in% levels(DC)[19:24]])

plot_crimes_hourly_set4 = ggplot(data = DF, aes(x = DATETIME, y = CRIMES, color = DC)) +
						  geom_point(na.rm = TRUE) +
						  facet_wrap(~DC, scales = "free_y") +
						  labs(x = "Hour", y = "Crimes") +
						  ggtitle("Crimes Per Hour By DC") +
						  theme_light(base_size = 15) + 
						  theme(legend.position = "top",
								legend.key.size = unit(.25, "in")) + 
						  guides(color = guide_legend(nrow = 1, override.aes = list(size = 10)))

plot_crimes_hourly_set4

DF = data.frame(crimes_hourly[DC %in% levels(DC)[25]])

plot_crimes_hourly_set5 = ggplot(data = DF, aes(x = DATETIME, y = CRIMES, color = DC)) +
						  geom_point(na.rm = TRUE) +
						  facet_wrap(~DC, scales = "free_y") +
						  labs(x = "Hour", y = "Crimes") +
						  ggtitle("Crimes Per Hour By DC") +
						  theme_light(base_size = 15) + 
						  theme(legend.position = "top",
								legend.key.size = unit(.25, "in")) + 
						  guides(color = guide_legend(nrow = 1, override.aes = list(size = 10)))

plot_crimes_hourly_set5

# lets put these results in lists to organize our workspace

crime_DC_tables = list("crimes_hourly" = crimes_hourly)
crime_DC_plots = list("plot_crimes_hourly_set1" = plot_crimes_hourly_set1, 
					  "plot_crimes_hourly_set2" = plot_crimes_hourly_set2,
					  "plot_crimes_hourly_set3" = plot_crimes_hourly_set3,
					  "plot_crimes_hourly_set4" = plot_crimes_hourly_set4,
					  "plot_crimes_hourly_set5" = plot_crimes_hourly_set5)

rm(crimes_hourly, plot_crimes_hourly_set1, plot_crimes_hourly_set2, plot_crimes_hourly_set3, plot_crimes_hourly_set4, plot_crimes_hourly_set5)

# lets create a table of crimes per day ------------------------------------------------------------------------------------------------------------------------------------

crimes_daily = PHI[,.("CRIMES" = .N), by = .("DC" = Dc_Dist, "DATE" = Dispatch_Date)]

# lets plot crimes_daily

DF = data.frame(crimes_daily[DC %in% levels(DC)[1:6]])

plot_crimes_daily_set1 = ggplot(data = DF, aes(x = DATE, y = CRIMES, color = DC)) +
						 geom_point(na.rm = TRUE) +
						 facet_wrap(~DC, scales = "free_y") +
						 labs(x = "Day", y = "Crimes") +
						 ggtitle("Crimes Per Day By DC") +
						 theme_light(base_size = 15) + 
						 theme(legend.position = "top",
						       legend.key.size = unit(.25, "in")) + 
						 guides(color = guide_legend(nrow = 1, override.aes = list(size = 10)))

plot_crimes_daily_set1

DF = data.frame(crimes_daily[DC %in% levels(DC)[7:12]])

plot_crimes_daily_set2 = ggplot(data = DF, aes(x = DATE, y = CRIMES, color = DC)) +
						 geom_point(na.rm = TRUE) +
						 facet_wrap(~DC, scales = "free_y") +
						 labs(x = "Day", y = "Crimes") +
						 ggtitle("Crimes Per Day By DC") +
						 theme_light(base_size = 15) + 
						 theme(legend.position = "top",
						       legend.key.size = unit(.25, "in")) + 
						 guides(color = guide_legend(nrow = 1, override.aes = list(size = 10)))

plot_crimes_daily_set2

DF = data.frame(crimes_daily[DC %in% levels(DC)[13:18]])

plot_crimes_daily_set3 = ggplot(data = DF, aes(x = DATE, y = CRIMES, color = DC)) +
						 geom_point(na.rm = TRUE) +
						 facet_wrap(~DC, scales = "free_y") +
						 labs(x = "Day", y = "Crimes") +
						 ggtitle("Crimes Per Day By DC") +
						 theme_light(base_size = 15) + 
						 theme(legend.position = "top",
						       legend.key.size = unit(.25, "in")) + 
						 guides(color = guide_legend(nrow = 1, override.aes = list(size = 10)))

plot_crimes_daily_set3

DF = data.frame(crimes_daily[DC %in% levels(DC)[19:24]])

plot_crimes_daily_set4 = ggplot(data = DF, aes(x = DATE, y = CRIMES, color = DC)) +
						 geom_point(na.rm = TRUE) +
						 facet_wrap(~DC, scales = "free_y") +
						 labs(x = "Day", y = "Crimes") +
						 ggtitle("Crimes Per Day By DC") +
						 theme_light(base_size = 15) + 
						 theme(legend.position = "top",
						       legend.key.size = unit(.25, "in")) + 
						 guides(color = guide_legend(nrow = 1, override.aes = list(size = 10)))

plot_crimes_daily_set4

DF = data.frame(crimes_daily[DC %in% levels(DC)[25]])

plot_crimes_daily_set5 = ggplot(data = DF, aes(x = DATE, y = CRIMES, color = DC)) +
						 geom_point(na.rm = TRUE) +
						 facet_wrap(~DC, scales = "free_y") +
						 labs(x = "Day", y = "Crimes") +
						 ggtitle("Crimes Per Day By DC") +
						 theme_light(base_size = 15) + 
						 theme(legend.position = "top",
						       legend.key.size = unit(.25, "in")) + 
						 guides(color = guide_legend(nrow = 1, override.aes = list(size = 10)))

plot_crimes_daily_set5

# lets put these results in lists to organize our workspace

crime_DC_tables$crimes_daily = crimes_daily
crime_DC_plots$plot_crimes_daily_set1 = plot_crimes_daily_set1
crime_DC_plots$plot_crimes_daily_set2 = plot_crimes_daily_set2
crime_DC_plots$plot_crimes_daily_set3 = plot_crimes_daily_set3
crime_DC_plots$plot_crimes_daily_set4 = plot_crimes_daily_set4
crime_DC_plots$plot_crimes_daily_set5 = plot_crimes_daily_set5

rm(crimes_daily, plot_crimes_daily_set1, plot_crimes_daily_set2, plot_crimes_daily_set3, plot_crimes_daily_set4, plot_crimes_daily_set5)

# lets create a table of crimes per week ------------------------------------------------------------------------------------------------------------------------------------

crimes_weekly = PHI[,.("CRIMES" = .N), by = .("DC" = Dc_Dist, "YEAR" = Year, "WEEK" = Week)]
crimes_weekly[, YEAR.WEEK := factor(paste0(crimes_weekly[,YEAR], ".", crimes_weekly[,WEEK]),
							 levels = paste0(crimes_weekly[,YEAR], ".", crimes_weekly[,WEEK]))]
crimes_weekly[, c("YEAR", "WEEK") := NULL]
setcolorder(crimes_weekly, c("DC", "YEAR.WEEK", "CRIMES"))

# lets plot crimes_weekly

DF = data.frame(crimes_weekly[DC %in% levels(DC)[1:6]])

plot_crimes_weekly_set1 = ggplot(data = DF, aes(x = YEAR.WEEK, y = CRIMES, color = DC)) +
						  geom_point(na.rm = TRUE) +
						  scale_x_discrete(breaks = c("2006.01", "2007.01", "2008.01", "2009.01", "2010.01", "2011.01", "2012.01", "2013.01", "2014.01", "2015.01", "2016.01")) +
						  facet_wrap(~DC, scales = "free_y") +
						  labs(x = "Year.Week", y = "Crimes") +
						  ggtitle("Crimes Per Week By DC") +
						  theme_light(base_size = 15) + 
						  theme(legend.position = "top",
								legend.key.size = unit(.25, "in"),
								axis.text.x = element_text(angle = 45, hjust = 1)) + 
						  guides(color = guide_legend(nrow = 1, override.aes = list(size = 10)))

plot_crimes_weekly_set1

DF = data.frame(crimes_weekly[DC %in% levels(DC)[7:12]])

plot_crimes_weekly_set2 = ggplot(data = DF, aes(x = YEAR.WEEK, y = CRIMES, color = DC)) +
						  geom_point(na.rm = TRUE) +
						  scale_x_discrete(breaks = c("2006.01", "2007.01", "2008.01", "2009.01", "2010.01", "2011.01", "2012.01", "2013.01", "2014.01", "2015.01", "2016.01")) +
						  facet_wrap(~DC, scales = "free_y") +
						  labs(x = "Year.Week", y = "Crimes") +
						  ggtitle("Crimes Per Week By DC") +
						  theme_light(base_size = 15) + 
						  theme(legend.position = "top",
								legend.key.size = unit(.25, "in"),
								axis.text.x = element_text(angle = 45, hjust = 1)) + 
						  guides(color = guide_legend(nrow = 1, override.aes = list(size = 10)))

plot_crimes_weekly_set2

DF = data.frame(crimes_weekly[DC %in% levels(DC)[13:18]])

plot_crimes_weekly_set3 = ggplot(data = DF, aes(x = YEAR.WEEK, y = CRIMES, color = DC)) +
						  geom_point(na.rm = TRUE) +
						  scale_x_discrete(breaks = c("2006.01", "2007.01", "2008.01", "2009.01", "2010.01", "2011.01", "2012.01", "2013.01", "2014.01", "2015.01", "2016.01")) +
						  facet_wrap(~DC, scales = "free_y") +
						  labs(x = "Year.Week", y = "Crimes") +
						  ggtitle("Crimes Per Week By DC") +
						  theme_light(base_size = 15) + 
						  theme(legend.position = "top",
								legend.key.size = unit(.25, "in"),
								axis.text.x = element_text(angle = 45, hjust = 1)) + 
						  guides(color = guide_legend(nrow = 1, override.aes = list(size = 10)))

plot_crimes_weekly_set3

DF = data.frame(crimes_weekly[DC %in% levels(DC)[19:24]])

plot_crimes_weekly_set4 = ggplot(data = DF, aes(x = YEAR.WEEK, y = CRIMES, color = DC)) +
						  geom_point(na.rm = TRUE) +
						  scale_x_discrete(breaks = c("2006.01", "2007.01", "2008.01", "2009.01", "2010.01", "2011.01", "2012.01", "2013.01", "2014.01", "2015.01", "2016.01")) +
						  facet_wrap(~DC, scales = "free_y") +
						  labs(x = "Year.Week", y = "Crimes") +
						  ggtitle("Crimes Per Week By DC") +
						  theme_light(base_size = 15) + 
						  theme(legend.position = "top",
								legend.key.size = unit(.25, "in"),
								axis.text.x = element_text(angle = 45, hjust = 1)) + 
						  guides(color = guide_legend(nrow = 1, override.aes = list(size = 10)))

plot_crimes_weekly_set4

DF = data.frame(crimes_weekly[DC %in% levels(DC)[25]])

plot_crimes_weekly_set5 = ggplot(data = DF, aes(x = YEAR.WEEK, y = CRIMES, color = DC)) +
						  geom_point(na.rm = TRUE) +
						  scale_x_discrete(breaks = c("2006.01", "2007.01", "2008.01", "2009.01", "2010.01", "2011.01", "2012.01", "2013.01", "2014.01", "2015.01", "2016.01")) +
						  facet_wrap(~DC, scales = "free_y") +
						  labs(x = "Year.Week", y = "Crimes") +
						  ggtitle("Crimes Per Week By DC") +
						  theme_light(base_size = 15) + 
						  theme(legend.position = "top",
								legend.key.size = unit(.25, "in"),
								axis.text.x = element_text(angle = 45, hjust = 1)) + 
						  guides(color = guide_legend(nrow = 1, override.aes = list(size = 10)))

plot_crimes_weekly_set5

# lets put these results in lists to organize our workspace

crime_DC_tables$crimes_weekly = crimes_weekly
crime_DC_plots$plot_crimes_weekly_set1 = plot_crimes_weekly_set1
crime_DC_plots$plot_crimes_weekly_set2 = plot_crimes_weekly_set2
crime_DC_plots$plot_crimes_weekly_set3 = plot_crimes_weekly_set3
crime_DC_plots$plot_crimes_weekly_set4 = plot_crimes_weekly_set4
crime_DC_plots$plot_crimes_weekly_set5 = plot_crimes_weekly_set5

rm(crimes_weekly, plot_crimes_weekly_set1, plot_crimes_weekly_set2, plot_crimes_weekly_set3, plot_crimes_weekly_set4, plot_crimes_weekly_set5)

# lets create a table of crimes per month ------------------------------------------------------------------------------------------------------------------------------------

crimes_monthly = PHI[,.("CRIMES" = .N), by = .("DC" = Dc_Dist, "YEAR" = Year, "MONTH" = Month)]
crimes_monthly[, YEAR.MONTH := factor(paste0(crimes_monthly[,YEAR], ".", crimes_monthly[,MONTH]),
							 levels = paste0(crimes_monthly[,YEAR], ".", crimes_monthly[,MONTH]))]
crimes_monthly[, c("YEAR", "MONTH") := NULL]
setcolorder(crimes_monthly, c("DC", "YEAR.MONTH", "CRIMES"))

# lets plot crimes_monthly

DF = data.frame(crimes_monthly[DC %in% levels(DC)[1:6]])

plot_crimes_monthly_set1 = ggplot(data = DF, aes(x = YEAR.MONTH, y = CRIMES, color = DC)) +
						   geom_point(na.rm = TRUE) +
						   scale_x_discrete(breaks = c("2006.Jan", "2007.Jan", "2008.Jan", "2009.Jan", "2010.Jan", "2011.Jan", "2012.Jan", "2013.Jan", "2014.Jan", "2015.Jan", "2016.Jan")) +
						   facet_wrap(~DC, scales = "free_y") +
						   labs(x = "Year.Month", y = "Crimes") +
						   ggtitle("Crimes Per Month By DC") +
						   theme_light(base_size = 15) + 
						   theme(legend.position = "top",
								 legend.key.size = unit(.25, "in"),
								 axis.text.x = element_text(angle = 45, hjust = 1)) + 
						   guides(color = guide_legend(nrow = 1, override.aes = list(size = 10)))

plot_crimes_monthly_set1

DF = data.frame(crimes_monthly[DC %in% levels(DC)[7:12]])

plot_crimes_monthly_set2 = ggplot(data = DF, aes(x = YEAR.MONTH, y = CRIMES, color = DC)) +
						   geom_point(na.rm = TRUE) +
						   scale_x_discrete(breaks = c("2006.Jan", "2007.Jan", "2008.Jan", "2009.Jan", "2010.Jan", "2011.Jan", "2012.Jan", "2013.Jan", "2014.Jan", "2015.Jan", "2016.Jan")) +
						   facet_wrap(~DC, scales = "free_y") +
						   labs(x = "Year.Month", y = "Crimes") +
						   ggtitle("Crimes Per Month By DC") +
						   theme_light(base_size = 15) + 
						   theme(legend.position = "top",
								 legend.key.size = unit(.25, "in"),
								 axis.text.x = element_text(angle = 45, hjust = 1)) + 
						   guides(color = guide_legend(nrow = 1, override.aes = list(size = 10)))

plot_crimes_monthly_set2

DF = data.frame(crimes_monthly[DC %in% levels(DC)[13:18]])

plot_crimes_monthly_set3 = ggplot(data = DF, aes(x = YEAR.MONTH, y = CRIMES, color = DC)) +
						   geom_point(na.rm = TRUE) +
						   scale_x_discrete(breaks = c("2006.Jan", "2007.Jan", "2008.Jan", "2009.Jan", "2010.Jan", "2011.Jan", "2012.Jan", "2013.Jan", "2014.Jan", "2015.Jan", "2016.Jan")) +
						   facet_wrap(~DC, scales = "free_y") +
						   labs(x = "Year.Month", y = "Crimes") +
						   ggtitle("Crimes Per Month By DC") +
						   theme_light(base_size = 15) + 
						   theme(legend.position = "top",
								 legend.key.size = unit(.25, "in"),
								 axis.text.x = element_text(angle = 45, hjust = 1)) + 
						   guides(color = guide_legend(nrow = 1, override.aes = list(size = 10)))

plot_crimes_monthly_set3

DF = data.frame(crimes_monthly[DC %in% levels(DC)[19:24]])

plot_crimes_monthly_set4 = ggplot(data = DF, aes(x = YEAR.MONTH, y = CRIMES, color = DC)) +
						   geom_point(na.rm = TRUE) +
						   scale_x_discrete(breaks = c("2006.Jan", "2007.Jan", "2008.Jan", "2009.Jan", "2010.Jan", "2011.Jan", "2012.Jan", "2013.Jan", "2014.Jan", "2015.Jan", "2016.Jan")) +
						   facet_wrap(~DC, scales = "free_y") +
						   labs(x = "Year.Month", y = "Crimes") +
						   ggtitle("Crimes Per Month By DC") +
						   theme_light(base_size = 15) + 
						   theme(legend.position = "top",
								 legend.key.size = unit(.25, "in"),
								 axis.text.x = element_text(angle = 45, hjust = 1)) + 
						   guides(color = guide_legend(nrow = 1, override.aes = list(size = 10)))

plot_crimes_monthly_set4

DF = data.frame(crimes_monthly[DC %in% levels(DC)[25]])

plot_crimes_monthly_set5 = ggplot(data = DF, aes(x = YEAR.MONTH, y = CRIMES, color = DC)) +
						   geom_point(na.rm = TRUE) +
						   scale_x_discrete(breaks = c("2006.Jan", "2007.Jan", "2008.Jan", "2009.Jan", "2010.Jan", "2011.Jan", "2012.Jan", "2013.Jan", "2014.Jan", "2015.Jan", "2016.Jan")) +
						   facet_wrap(~DC, scales = "free_y") +
						   labs(x = "Year.Month", y = "Crimes") +
						   ggtitle("Crimes Per Month By DC") +
						   theme_light(base_size = 15) + 
						   theme(legend.position = "top",
								 legend.key.size = unit(.25, "in"),
								 axis.text.x = element_text(angle = 45, hjust = 1)) + 
						   guides(color = guide_legend(nrow = 1, override.aes = list(size = 10)))

plot_crimes_monthly_set5

# lets put these results in lists to organize our workspace

crime_DC_tables$crimes_monthly = crimes_monthly
crime_DC_plots$plot_crimes_monthly_set1 = plot_crimes_monthly_set1
crime_DC_plots$plot_crimes_monthly_set2 = plot_crimes_monthly_set2
crime_DC_plots$plot_crimes_monthly_set3 = plot_crimes_monthly_set3
crime_DC_plots$plot_crimes_monthly_set4 = plot_crimes_monthly_set4
crime_DC_plots$plot_crimes_monthly_set5 = plot_crimes_monthly_set5

rm(crimes_monthly, plot_crimes_monthly_set1, plot_crimes_monthly_set2, plot_crimes_monthly_set3, plot_crimes_monthly_set4, plot_crimes_monthly_set5)

}

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ---- Crimes By Psa -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

{

# lets create a table of crimes per hour ------------------------------------------------------------------------------------------------------------------------------------

crimes_hourly = PHI[,.("CRIMES" = .N), by = .("PSA" = Psa, "DATETIME" = Dispatch_Date_Time, "HOUR" = Hour)]

# lets plot crimes_hourly

DF = data.frame(crimes_hourly[PSA %in% levels(PSA)[1:6]])

plot_crimes_hourly_set1 = ggplot(data = DF, aes(x = DATETIME, y = CRIMES, color = PSA)) +
						  geom_point(na.rm = TRUE) +
						  facet_wrap(~PSA, scales = "free_y") +
						  labs(x = "Hour", y = "Crimes") +
						  ggtitle("Crimes Per Hour By PSA") +
						  theme_light(base_size = 15) + 
						  theme(legend.position = "top",
								legend.key.size = unit(.25, "in")) + 
						  guides(color = guide_legend(nrow = 1, override.aes = list(size = 10)))

plot_crimes_hourly_set1

DF = data.frame(crimes_hourly[PSA %in% levels(PSA)[7:12]])

plot_crimes_hourly_set2 = ggplot(data = DF, aes(x = DATETIME, y = CRIMES, color = PSA)) +
						  geom_point(na.rm = TRUE) +
						  facet_wrap(~PSA, scales = "free_y") +
						  labs(x = "Hour", y = "Crimes") +
						  ggtitle("Crimes Per Hour By PSA") +
						  theme_light(base_size = 15) + 
						  theme(legend.position = "top",
								legend.key.size = unit(.25, "in")) + 
						  guides(color = guide_legend(nrow = 1, override.aes = list(size = 10)))

plot_crimes_hourly_set2

DF = data.frame(crimes_hourly[PSA %in% levels(PSA)[13:18]])

plot_crimes_hourly_set3 = ggplot(data = DF, aes(x = DATETIME, y = CRIMES, color = PSA)) +
						  geom_point(na.rm = TRUE) +
						  facet_wrap(~PSA, scales = "free_y") +
						  labs(x = "Hour", y = "Crimes") +
						  ggtitle("Crimes Per Hour By PSA") +
						  theme_light(base_size = 15) + 
						  theme(legend.position = "top",
								legend.key.size = unit(.25, "in")) + 
						  guides(color = guide_legend(nrow = 1, override.aes = list(size = 10)))

plot_crimes_hourly_set3

DF = data.frame(crimes_hourly[PSA %in% levels(PSA)[19:24]])

plot_crimes_hourly_set4 = ggplot(data = DF, aes(x = DATETIME, y = CRIMES, color = PSA)) +
						  geom_point(na.rm = TRUE) +
						  facet_wrap(~PSA, scales = "free_y") +
						  labs(x = "Hour", y = "Crimes") +
						  ggtitle("Crimes Per Hour By PSA") +
						  theme_light(base_size = 15) + 
						  theme(legend.position = "top",
								legend.key.size = unit(.25, "in")) + 
						  guides(color = guide_legend(nrow = 1, override.aes = list(size = 10)))

plot_crimes_hourly_set4

DF = data.frame(crimes_hourly[PSA %in% levels(PSA)[25:30]])

plot_crimes_hourly_set5 = ggplot(data = DF, aes(x = DATETIME, y = CRIMES, color = PSA)) +
						  geom_point(na.rm = TRUE) +
						  facet_wrap(~PSA, scales = "free_y") +
						  labs(x = "Hour", y = "Crimes") +
						  ggtitle("Crimes Per Hour By PSA") +
						  theme_light(base_size = 15) + 
						  theme(legend.position = "top",
								legend.key.size = unit(.25, "in")) + 
						  guides(color = guide_legend(nrow = 1, override.aes = list(size = 10)))

plot_crimes_hourly_set5

# lets put these results in lists to organize our workspace

crime_PSA_tables = list("crimes_hourly" = crimes_hourly)
crime_PSA_plots = list("plot_crimes_hourly_set1" = plot_crimes_hourly_set1, 
					  "plot_crimes_hourly_set2" = plot_crimes_hourly_set2,
					  "plot_crimes_hourly_set3" = plot_crimes_hourly_set3,
					  "plot_crimes_hourly_set4" = plot_crimes_hourly_set4,
					  "plot_crimes_hourly_set5" = plot_crimes_hourly_set5)

rm(crimes_hourly, plot_crimes_hourly_set1, plot_crimes_hourly_set2, plot_crimes_hourly_set3, plot_crimes_hourly_set4, plot_crimes_hourly_set5)

# lets create a table of crimes per day ------------------------------------------------------------------------------------------------------------------------------------

crimes_daily = PHI[,.("CRIMES" = .N), by = .("PSA" = Psa, "DATE" = Dispatch_Date)]

# lets plot crimes_daily

DF = data.frame(crimes_daily[PSA %in% levels(PSA)[1:6]])

plot_crimes_daily_set1 = ggplot(data = DF, aes(x = DATE, y = CRIMES, color = PSA)) +
						 geom_point(na.rm = TRUE) +
						 facet_wrap(~PSA, scales = "free_y") +
						 labs(x = "Day", y = "Crimes") +
						 ggtitle("Crimes Per Day By PSA") +
						 theme_light(base_size = 15) + 
						 theme(legend.position = "top",
						       legend.key.size = unit(.25, "in")) + 
						 guides(color = guide_legend(nrow = 1, override.aes = list(size = 10)))

plot_crimes_daily_set1

DF = data.frame(crimes_daily[PSA %in% levels(PSA)[7:12]])

plot_crimes_daily_set2 = ggplot(data = DF, aes(x = DATE, y = CRIMES, color = PSA)) +
						 geom_point(na.rm = TRUE) +
						 facet_wrap(~PSA, scales = "free_y") +
						 labs(x = "Day", y = "Crimes") +
						 ggtitle("Crimes Per Day By PSA") +
						 theme_light(base_size = 15) + 
						 theme(legend.position = "top",
						       legend.key.size = unit(.25, "in")) + 
						 guides(color = guide_legend(nrow = 1, override.aes = list(size = 10)))

plot_crimes_daily_set2

DF = data.frame(crimes_daily[PSA %in% levels(PSA)[13:18]])

plot_crimes_daily_set3 = ggplot(data = DF, aes(x = DATE, y = CRIMES, color = PSA)) +
						 geom_point(na.rm = TRUE) +
						 facet_wrap(~PSA, scales = "free_y") +
						 labs(x = "Day", y = "Crimes") +
						 ggtitle("Crimes Per Day By PSA") +
						 theme_light(base_size = 15) + 
						 theme(legend.position = "top",
						       legend.key.size = unit(.25, "in")) + 
						 guides(color = guide_legend(nrow = 1, override.aes = list(size = 10)))

plot_crimes_daily_set3

DF = data.frame(crimes_daily[PSA %in% levels(PSA)[19:24]])

plot_crimes_daily_set4 = ggplot(data = DF, aes(x = DATE, y = CRIMES, color = PSA)) +
						 geom_point(na.rm = TRUE) +
						 facet_wrap(~PSA, scales = "free_y") +
						 labs(x = "Day", y = "Crimes") +
						 ggtitle("Crimes Per Day By PSA") +
						 theme_light(base_size = 15) + 
						 theme(legend.position = "top",
						       legend.key.size = unit(.25, "in")) + 
						 guides(color = guide_legend(nrow = 1, override.aes = list(size = 10)))

plot_crimes_daily_set4

DF = data.frame(crimes_daily[PSA %in% levels(PSA)[25:30]])

plot_crimes_daily_set5 = ggplot(data = DF, aes(x = DATE, y = CRIMES, color = PSA)) +
						 geom_point(na.rm = TRUE) +
						 facet_wrap(~PSA, scales = "free_y") +
						 labs(x = "Day", y = "Crimes") +
						 ggtitle("Crimes Per Day By PSA") +
						 theme_light(base_size = 15) + 
						 theme(legend.position = "top",
						       legend.key.size = unit(.25, "in")) + 
						 guides(color = guide_legend(nrow = 1, override.aes = list(size = 10)))

plot_crimes_daily_set5

# lets put these results in lists to organize our workspace

crime_PSA_tables$crimes_daily = crimes_daily
crime_PSA_plots$plot_crimes_daily_set1 = plot_crimes_daily_set1
crime_PSA_plots$plot_crimes_daily_set2 = plot_crimes_daily_set2
crime_PSA_plots$plot_crimes_daily_set3 = plot_crimes_daily_set3
crime_PSA_plots$plot_crimes_daily_set4 = plot_crimes_daily_set4
crime_PSA_plots$plot_crimes_daily_set5 = plot_crimes_daily_set5

rm(crimes_daily, plot_crimes_daily_set1, plot_crimes_daily_set2, plot_crimes_daily_set3, plot_crimes_daily_set4, plot_crimes_daily_set5)

# lets create a table of crimes per week ------------------------------------------------------------------------------------------------------------------------------------

crimes_weekly = PHI[,.("CRIMES" = .N), by = .("PSA" = Psa, "YEAR" = Year, "WEEK" = Week)]
crimes_weekly[, YEAR.WEEK := factor(paste0(crimes_weekly[,YEAR], ".", crimes_weekly[,WEEK]),
							 levels = paste0(crimes_weekly[,YEAR], ".", crimes_weekly[,WEEK]))]
crimes_weekly[, c("YEAR", "WEEK") := NULL]
setcolorder(crimes_weekly, c("PSA", "YEAR.WEEK", "CRIMES"))

# lets plot crimes_weekly

DF = data.frame(crimes_weekly[PSA %in% levels(PSA)[1:6]])

plot_crimes_weekly_set1 = ggplot(data = DF, aes(x = YEAR.WEEK, y = CRIMES, color = PSA)) +
						  geom_point(na.rm = TRUE) +
						  scale_x_discrete(breaks = c("2006.01", "2007.01", "2008.01", "2009.01", "2010.01", "2011.01", "2012.01", "2013.01", "2014.01", "2015.01", "2016.01")) +
						  facet_wrap(~PSA, scales = "free_y") +
						  labs(x = "Year.Week", y = "Crimes") +
						  ggtitle("Crimes Per Week By PSA") +
						  theme_light(base_size = 15) + 
						  theme(legend.position = "top",
								legend.key.size = unit(.25, "in"),
								axis.text.x = element_text(angle = 45, hjust = 1)) + 
						  guides(color = guide_legend(nrow = 1, override.aes = list(size = 10)))

plot_crimes_weekly_set1

DF = data.frame(crimes_weekly[PSA %in% levels(PSA)[7:12]])

plot_crimes_weekly_set2 = ggplot(data = DF, aes(x = YEAR.WEEK, y = CRIMES, color = PSA)) +
						  geom_point(na.rm = TRUE) +
						  scale_x_discrete(breaks = c("2006.01", "2007.01", "2008.01", "2009.01", "2010.01", "2011.01", "2012.01", "2013.01", "2014.01", "2015.01", "2016.01")) +
						  facet_wrap(~PSA, scales = "free_y") +
						  labs(x = "Year.Week", y = "Crimes") +
						  ggtitle("Crimes Per Week By PSA") +
						  theme_light(base_size = 15) + 
						  theme(legend.position = "top",
								legend.key.size = unit(.25, "in"),
								axis.text.x = element_text(angle = 45, hjust = 1)) + 
						  guides(color = guide_legend(nrow = 1, override.aes = list(size = 10)))

plot_crimes_weekly_set2

DF = data.frame(crimes_weekly[PSA %in% levels(PSA)[13:18]])

plot_crimes_weekly_set3 = ggplot(data = DF, aes(x = YEAR.WEEK, y = CRIMES, color = PSA)) +
						  geom_point(na.rm = TRUE) +
						  scale_x_discrete(breaks = c("2006.01", "2007.01", "2008.01", "2009.01", "2010.01", "2011.01", "2012.01", "2013.01", "2014.01", "2015.01", "2016.01")) +
						  facet_wrap(~PSA, scales = "free_y") +
						  labs(x = "Year.Week", y = "Crimes") +
						  ggtitle("Crimes Per Week By PSA") +
						  theme_light(base_size = 15) + 
						  theme(legend.position = "top",
								legend.key.size = unit(.25, "in"),
								axis.text.x = element_text(angle = 45, hjust = 1)) + 
						  guides(color = guide_legend(nrow = 1, override.aes = list(size = 10)))

plot_crimes_weekly_set3

DF = data.frame(crimes_weekly[PSA %in% levels(PSA)[19:24]])

plot_crimes_weekly_set4 = ggplot(data = DF, aes(x = YEAR.WEEK, y = CRIMES, color = PSA)) +
						  geom_point(na.rm = TRUE) +
						  scale_x_discrete(breaks = c("2006.01", "2007.01", "2008.01", "2009.01", "2010.01", "2011.01", "2012.01", "2013.01", "2014.01", "2015.01", "2016.01")) +
						  facet_wrap(~PSA, scales = "free_y") +
						  labs(x = "Year.Week", y = "Crimes") +
						  ggtitle("Crimes Per Week By PSA") +
						  theme_light(base_size = 15) + 
						  theme(legend.position = "top",
								legend.key.size = unit(.25, "in"),
								axis.text.x = element_text(angle = 45, hjust = 1)) + 
						  guides(color = guide_legend(nrow = 1, override.aes = list(size = 10)))

plot_crimes_weekly_set4

DF = data.frame(crimes_weekly[PSA %in% levels(PSA)[25:30]])

plot_crimes_weekly_set5 = ggplot(data = DF, aes(x = YEAR.WEEK, y = CRIMES, color = PSA)) +
						  geom_point(na.rm = TRUE) +
						  scale_x_discrete(breaks = c("2006.01", "2007.01", "2008.01", "2009.01", "2010.01", "2011.01", "2012.01", "2013.01", "2014.01", "2015.01", "2016.01")) +
						  facet_wrap(~PSA, scales = "free_y") +
						  labs(x = "Year.Week", y = "Crimes") +
						  ggtitle("Crimes Per Week By PSA") +
						  theme_light(base_size = 15) + 
						  theme(legend.position = "top",
								legend.key.size = unit(.25, "in"),
								axis.text.x = element_text(angle = 45, hjust = 1)) + 
						  guides(color = guide_legend(nrow = 1, override.aes = list(size = 10)))

plot_crimes_weekly_set5

# lets put these results in lists to organize our workspace

crime_PSA_tables$crimes_weekly = crimes_weekly
crime_PSA_plots$plot_crimes_weekly_set1 = plot_crimes_weekly_set1
crime_PSA_plots$plot_crimes_weekly_set2 = plot_crimes_weekly_set2
crime_PSA_plots$plot_crimes_weekly_set3 = plot_crimes_weekly_set3
crime_PSA_plots$plot_crimes_weekly_set4 = plot_crimes_weekly_set4
crime_PSA_plots$plot_crimes_weekly_set5 = plot_crimes_weekly_set5

rm(crimes_weekly, plot_crimes_weekly_set1, plot_crimes_weekly_set2, plot_crimes_weekly_set3, plot_crimes_weekly_set4, plot_crimes_weekly_set5)

# lets create a table of crimes per month ------------------------------------------------------------------------------------------------------------------------------------

crimes_monthly = PHI[,.("CRIMES" = .N), by = .("PSA" = Psa, "YEAR" = Year, "MONTH" = Month)]
crimes_monthly[, YEAR.MONTH := factor(paste0(crimes_monthly[,YEAR], ".", crimes_monthly[,MONTH]),
							 levels = paste0(crimes_monthly[,YEAR], ".", crimes_monthly[,MONTH]))]
crimes_monthly[, c("YEAR", "MONTH") := NULL]
setcolorder(crimes_monthly, c("PSA", "YEAR.MONTH", "CRIMES"))

# lets plot crimes_monthly

DF = data.frame(crimes_monthly[PSA %in% levels(PSA)[1:6]])

plot_crimes_monthly_set1 = ggplot(data = DF, aes(x = YEAR.MONTH, y = CRIMES, color = PSA)) +
						   geom_point(na.rm = TRUE) +
						   scale_x_discrete(breaks = c("2006.Jan", "2007.Jan", "2008.Jan", "2009.Jan", "2010.Jan", "2011.Jan", "2012.Jan", "2013.Jan", "2014.Jan", "2015.Jan", "2016.Jan")) +
						   facet_wrap(~PSA, scales = "free_y") +
						   labs(x = "Year.Month", y = "Crimes") +
						   ggtitle("Crimes Per Month By PSA") +
						   theme_light(base_size = 15) + 
						   theme(legend.position = "top",
								 legend.key.size = unit(.25, "in"),
								 axis.text.x = element_text(angle = 45, hjust = 1)) + 
						   guides(color = guide_legend(nrow = 1, override.aes = list(size = 10)))

plot_crimes_monthly_set1

DF = data.frame(crimes_monthly[PSA %in% levels(PSA)[7:12]])

plot_crimes_monthly_set2 = ggplot(data = DF, aes(x = YEAR.MONTH, y = CRIMES, color = PSA)) +
						   geom_point(na.rm = TRUE) +
						   scale_x_discrete(breaks = c("2006.Jan", "2007.Jan", "2008.Jan", "2009.Jan", "2010.Jan", "2011.Jan", "2012.Jan", "2013.Jan", "2014.Jan", "2015.Jan", "2016.Jan")) +
						   facet_wrap(~PSA, scales = "free_y") +
						   labs(x = "Year.Month", y = "Crimes") +
						   ggtitle("Crimes Per Month By PSA") +
						   theme_light(base_size = 15) + 
						   theme(legend.position = "top",
								 legend.key.size = unit(.25, "in"),
								 axis.text.x = element_text(angle = 45, hjust = 1)) + 
						   guides(color = guide_legend(nrow = 1, override.aes = list(size = 10)))

plot_crimes_monthly_set2

DF = data.frame(crimes_monthly[PSA %in% levels(PSA)[13:18]])

plot_crimes_monthly_set3 = ggplot(data = DF, aes(x = YEAR.MONTH, y = CRIMES, color = PSA)) +
						   geom_point(na.rm = TRUE) +
						   scale_x_discrete(breaks = c("2006.Jan", "2007.Jan", "2008.Jan", "2009.Jan", "2010.Jan", "2011.Jan", "2012.Jan", "2013.Jan", "2014.Jan", "2015.Jan", "2016.Jan")) +
						   facet_wrap(~PSA, scales = "free_y") +
						   labs(x = "Year.Month", y = "Crimes") +
						   ggtitle("Crimes Per Month By PSA") +
						   theme_light(base_size = 15) + 
						   theme(legend.position = "top",
								 legend.key.size = unit(.25, "in"),
								 axis.text.x = element_text(angle = 45, hjust = 1)) + 
						   guides(color = guide_legend(nrow = 1, override.aes = list(size = 10)))

plot_crimes_monthly_set3

DF = data.frame(crimes_monthly[PSA %in% levels(PSA)[19:24]])

plot_crimes_monthly_set4 = ggplot(data = DF, aes(x = YEAR.MONTH, y = CRIMES, color = PSA)) +
						   geom_point(na.rm = TRUE) +
						   scale_x_discrete(breaks = c("2006.Jan", "2007.Jan", "2008.Jan", "2009.Jan", "2010.Jan", "2011.Jan", "2012.Jan", "2013.Jan", "2014.Jan", "2015.Jan", "2016.Jan")) +
						   facet_wrap(~PSA, scales = "free_y") +
						   labs(x = "Year.Month", y = "Crimes") +
						   ggtitle("Crimes Per Month By PSA") +
						   theme_light(base_size = 15) + 
						   theme(legend.position = "top",
								 legend.key.size = unit(.25, "in"),
								 axis.text.x = element_text(angle = 45, hjust = 1)) + 
						   guides(color = guide_legend(nrow = 1, override.aes = list(size = 10)))

plot_crimes_monthly_set4

DF = data.frame(crimes_monthly[PSA %in% levels(PSA)[25:30]])

plot_crimes_monthly_set5 = ggplot(data = DF, aes(x = YEAR.MONTH, y = CRIMES, color = PSA)) +
						   geom_point(na.rm = TRUE) +
						   scale_x_discrete(breaks = c("2006.Jan", "2007.Jan", "2008.Jan", "2009.Jan", "2010.Jan", "2011.Jan", "2012.Jan", "2013.Jan", "2014.Jan", "2015.Jan", "2016.Jan")) +
						   facet_wrap(~PSA, scales = "free_y") +
						   labs(x = "Year.Month", y = "Crimes") +
						   ggtitle("Crimes Per Month By PSA") +
						   theme_light(base_size = 15) + 
						   theme(legend.position = "top",
								 legend.key.size = unit(.25, "in"),
								 axis.text.x = element_text(angle = 45, hjust = 1)) + 
						   guides(color = guide_legend(nrow = 1, override.aes = list(size = 10)))

plot_crimes_monthly_set5

# lets put these results in lists to organize our workspace

crime_PSA_tables$crimes_monthly = crimes_monthly
crime_PSA_plots$plot_crimes_monthly_set1 = plot_crimes_monthly_set1
crime_PSA_plots$plot_crimes_monthly_set2 = plot_crimes_monthly_set2
crime_PSA_plots$plot_crimes_monthly_set3 = plot_crimes_monthly_set3
crime_PSA_plots$plot_crimes_monthly_set4 = plot_crimes_monthly_set4
crime_PSA_plots$plot_crimes_monthly_set5 = plot_crimes_monthly_set5

rm(crimes_monthly, plot_crimes_monthly_set1, plot_crimes_monthly_set2, plot_crimes_monthly_set3, plot_crimes_monthly_set4, plot_crimes_monthly_set5)

}

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ---- Crimes By UCR_General -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

{

# lets create a table of crimes per hour ------------------------------------------------------------------------------------------------------------------------------------

crimes_hourly = PHI[,.("CRIMES" = .N), by = .("UCR" = UCR_General, "DATETIME" = Dispatch_Date_Time, "HOUR" = Hour)]

# lets plot crimes_hourly

DF = data.frame(crimes_hourly[UCR %in% levels(UCR)[1:6]])

plot_crimes_hourly_set1 = ggplot(data = DF, aes(x = DATETIME, y = CRIMES, color = UCR)) +
						  geom_point(na.rm = TRUE) +
						  facet_wrap(~UCR, scales = "free_y") +
						  labs(x = "Hour", y = "Crimes") +
						  ggtitle("Crimes Per Hour By UCR") +
						  theme_light(base_size = 15) + 
						  theme(legend.position = "top",
								legend.key.size = unit(.25, "in")) + 
						  guides(color = guide_legend(nrow = 1, override.aes = list(size = 10)))

plot_crimes_hourly_set1

DF = data.frame(crimes_hourly[UCR %in% levels(UCR)[7:12]])

plot_crimes_hourly_set2 = ggplot(data = DF, aes(x = DATETIME, y = CRIMES, color = UCR)) +
						  geom_point(na.rm = TRUE) +
						  facet_wrap(~UCR, scales = "free_y") +
						  labs(x = "Hour", y = "Crimes") +
						  ggtitle("Crimes Per Hour By UCR") +
						  theme_light(base_size = 15) + 
						  theme(legend.position = "top",
								legend.key.size = unit(.25, "in")) + 
						  guides(color = guide_legend(nrow = 1, override.aes = list(size = 10)))

plot_crimes_hourly_set2

DF = data.frame(crimes_hourly[UCR %in% levels(UCR)[13:18]])

plot_crimes_hourly_set3 = ggplot(data = DF, aes(x = DATETIME, y = CRIMES, color = UCR)) +
						  geom_point(na.rm = TRUE) +
						  facet_wrap(~UCR, scales = "free_y") +
						  labs(x = "Hour", y = "Crimes") +
						  ggtitle("Crimes Per Hour By UCR") +
						  theme_light(base_size = 15) + 
						  theme(legend.position = "top",
								legend.key.size = unit(.25, "in")) + 
						  guides(color = guide_legend(nrow = 1, override.aes = list(size = 10)))

plot_crimes_hourly_set3

DF = data.frame(crimes_hourly[UCR %in% levels(UCR)[19:24]])

plot_crimes_hourly_set4 = ggplot(data = DF, aes(x = DATETIME, y = CRIMES, color = UCR)) +
						  geom_point(na.rm = TRUE) +
						  facet_wrap(~UCR, scales = "free_y") +
						  labs(x = "Hour", y = "Crimes") +
						  ggtitle("Crimes Per Hour By UCR") +
						  theme_light(base_size = 15) + 
						  theme(legend.position = "top",
								legend.key.size = unit(.25, "in")) + 
						  guides(color = guide_legend(nrow = 1, override.aes = list(size = 10)))

plot_crimes_hourly_set4

DF = data.frame(crimes_hourly[UCR %in% levels(UCR)[25:26]])

plot_crimes_hourly_set5 = ggplot(data = DF, aes(x = DATETIME, y = CRIMES, color = UCR)) +
						  geom_point(na.rm = TRUE) +
						  facet_wrap(~UCR, scales = "free_y") +
						  labs(x = "Hour", y = "Crimes") +
						  ggtitle("Crimes Per Hour By UCR") +
						  theme_light(base_size = 15) + 
						  theme(legend.position = "top",
								legend.key.size = unit(.25, "in")) + 
						  guides(color = guide_legend(nrow = 1, override.aes = list(size = 10)))

plot_crimes_hourly_set5

# lets put these results in lists to organize our workspace

crime_UCR_tables = list("crimes_hourly" = crimes_hourly)
crime_UCR_plots = list("plot_crimes_hourly_set1" = plot_crimes_hourly_set1, 
					  "plot_crimes_hourly_set2" = plot_crimes_hourly_set2,
					  "plot_crimes_hourly_set3" = plot_crimes_hourly_set3,
					  "plot_crimes_hourly_set4" = plot_crimes_hourly_set4,
					  "plot_crimes_hourly_set5" = plot_crimes_hourly_set5)

rm(crimes_hourly, plot_crimes_hourly_set1, plot_crimes_hourly_set2, plot_crimes_hourly_set3, plot_crimes_hourly_set4, plot_crimes_hourly_set5)

# lets create a table of crimes per day ------------------------------------------------------------------------------------------------------------------------------------

crimes_daily = PHI[,.("CRIMES" = .N), by = .("UCR" = UCR_General, "DATE" = Dispatch_Date)]

# lets plot crimes_daily

DF = data.frame(crimes_daily[UCR %in% levels(UCR)[1:6]])

plot_crimes_daily_set1 = ggplot(data = DF, aes(x = DATE, y = CRIMES, color = UCR)) +
						 geom_point(na.rm = TRUE) +
						 facet_wrap(~UCR, scales = "free_y") +
						 labs(x = "Day", y = "Crimes") +
						 ggtitle("Crimes Per Day By UCR") +
						 theme_light(base_size = 15) + 
						 theme(legend.position = "top",
						       legend.key.size = unit(.25, "in")) + 
						 guides(color = guide_legend(nrow = 1, override.aes = list(size = 10)))

plot_crimes_daily_set1

DF = data.frame(crimes_daily[UCR %in% levels(UCR)[7:12]])

plot_crimes_daily_set2 = ggplot(data = DF, aes(x = DATE, y = CRIMES, color = UCR)) +
						 geom_point(na.rm = TRUE) +
						 facet_wrap(~UCR, scales = "free_y") +
						 labs(x = "Day", y = "Crimes") +
						 ggtitle("Crimes Per Day By UCR") +
						 theme_light(base_size = 15) + 
						 theme(legend.position = "top",
						       legend.key.size = unit(.25, "in")) + 
						 guides(color = guide_legend(nrow = 1, override.aes = list(size = 10)))

plot_crimes_daily_set2

DF = data.frame(crimes_daily[UCR %in% levels(UCR)[13:18]])

plot_crimes_daily_set3 = ggplot(data = DF, aes(x = DATE, y = CRIMES, color = UCR)) +
						 geom_point(na.rm = TRUE) +
						 facet_wrap(~UCR, scales = "free_y") +
						 labs(x = "Day", y = "Crimes") +
						 ggtitle("Crimes Per Day By UCR") +
						 theme_light(base_size = 15) + 
						 theme(legend.position = "top",
						       legend.key.size = unit(.25, "in")) + 
						 guides(color = guide_legend(nrow = 1, override.aes = list(size = 10)))

plot_crimes_daily_set3

DF = data.frame(crimes_daily[UCR %in% levels(UCR)[19:24]])

plot_crimes_daily_set4 = ggplot(data = DF, aes(x = DATE, y = CRIMES, color = UCR)) +
						 geom_point(na.rm = TRUE) +
						 facet_wrap(~UCR, scales = "free_y") +
						 labs(x = "Day", y = "Crimes") +
						 ggtitle("Crimes Per Day By UCR") +
						 theme_light(base_size = 15) + 
						 theme(legend.position = "top",
						       legend.key.size = unit(.25, "in")) + 
						 guides(color = guide_legend(nrow = 1, override.aes = list(size = 10)))

plot_crimes_daily_set4

DF = data.frame(crimes_daily[UCR %in% levels(UCR)[25:26]])

plot_crimes_daily_set5 = ggplot(data = DF, aes(x = DATE, y = CRIMES, color = UCR)) +
						 geom_point(na.rm = TRUE) +
						 facet_wrap(~UCR, scales = "free_y") +
						 labs(x = "Day", y = "Crimes") +
						 ggtitle("Crimes Per Day By UCR") +
						 theme_light(base_size = 15) + 
						 theme(legend.position = "top",
						       legend.key.size = unit(.25, "in")) + 
						 guides(color = guide_legend(nrow = 1, override.aes = list(size = 10)))

plot_crimes_daily_set5

# lets put these results in lists to organize our workspace

crime_UCR_tables$crimes_daily = crimes_daily
crime_UCR_plots$plot_crimes_daily_set1 = plot_crimes_daily_set1
crime_UCR_plots$plot_crimes_daily_set2 = plot_crimes_daily_set2
crime_UCR_plots$plot_crimes_daily_set3 = plot_crimes_daily_set3
crime_UCR_plots$plot_crimes_daily_set4 = plot_crimes_daily_set4
crime_UCR_plots$plot_crimes_daily_set5 = plot_crimes_daily_set5

rm(crimes_daily, plot_crimes_daily_set1, plot_crimes_daily_set2, plot_crimes_daily_set3, plot_crimes_daily_set4, plot_crimes_daily_set5)

# lets create a table of crimes per week ------------------------------------------------------------------------------------------------------------------------------------

crimes_weekly = PHI[,.("CRIMES" = .N), by = .("UCR" = UCR_General, "YEAR" = Year, "WEEK" = Week)]
crimes_weekly[, YEAR.WEEK := factor(paste0(crimes_weekly[,YEAR], ".", crimes_weekly[,WEEK]),
							 levels = paste0(crimes_weekly[,YEAR], ".", crimes_weekly[,WEEK]))]
crimes_weekly[, c("YEAR", "WEEK") := NULL]
setcolorder(crimes_weekly, c("UCR", "YEAR.WEEK", "CRIMES"))

# lets plot crimes_weekly

DF = data.frame(crimes_weekly[UCR %in% levels(UCR)[1:6]])

plot_crimes_weekly_set1 = ggplot(data = DF, aes(x = YEAR.WEEK, y = CRIMES, color = UCR)) +
						  geom_point(na.rm = TRUE) +
						  scale_x_discrete(breaks = c("2006.01", "2007.01", "2008.01", "2009.01", "2010.01", "2011.01", "2012.01", "2013.01", "2014.01", "2015.01", "2016.01")) +
						  facet_wrap(~UCR, scales = "free_y") +
						  labs(x = "Year.Week", y = "Crimes") +
						  ggtitle("Crimes Per Week By UCR") +
						  theme_light(base_size = 15) + 
						  theme(legend.position = "top",
								legend.key.size = unit(.25, "in"),
								axis.text.x = element_text(angle = 45, hjust = 1)) + 
						  guides(color = guide_legend(nrow = 1, override.aes = list(size = 10)))

plot_crimes_weekly_set1

DF = data.frame(crimes_weekly[UCR %in% levels(UCR)[7:12]])

plot_crimes_weekly_set2 = ggplot(data = DF, aes(x = YEAR.WEEK, y = CRIMES, color = UCR)) +
						  geom_point(na.rm = TRUE) +
						  scale_x_discrete(breaks = c("2006.01", "2007.01", "2008.01", "2009.01", "2010.01", "2011.01", "2012.01", "2013.01", "2014.01", "2015.01", "2016.01")) +
						  facet_wrap(~UCR, scales = "free_y") +
						  labs(x = "Year.Week", y = "Crimes") +
						  ggtitle("Crimes Per Week By UCR") +
						  theme_light(base_size = 15) + 
						  theme(legend.position = "top",
								legend.key.size = unit(.25, "in"),
								axis.text.x = element_text(angle = 45, hjust = 1)) + 
						  guides(color = guide_legend(nrow = 1, override.aes = list(size = 10)))

plot_crimes_weekly_set2

DF = data.frame(crimes_weekly[UCR %in% levels(UCR)[13:18]])

plot_crimes_weekly_set3 = ggplot(data = DF, aes(x = YEAR.WEEK, y = CRIMES, color = UCR)) +
						  geom_point(na.rm = TRUE) +
						  scale_x_discrete(breaks = c("2006.01", "2007.01", "2008.01", "2009.01", "2010.01", "2011.01", "2012.01", "2013.01", "2014.01", "2015.01", "2016.01")) +
						  facet_wrap(~UCR, scales = "free_y") +
						  labs(x = "Year.Week", y = "Crimes") +
						  ggtitle("Crimes Per Week By UCR") +
						  theme_light(base_size = 15) + 
						  theme(legend.position = "top",
								legend.key.size = unit(.25, "in"),
								axis.text.x = element_text(angle = 45, hjust = 1)) + 
						  guides(color = guide_legend(nrow = 1, override.aes = list(size = 10)))

plot_crimes_weekly_set3

DF = data.frame(crimes_weekly[UCR %in% levels(UCR)[19:24]])

plot_crimes_weekly_set4 = ggplot(data = DF, aes(x = YEAR.WEEK, y = CRIMES, color = UCR)) +
						  geom_point(na.rm = TRUE) +
						  scale_x_discrete(breaks = c("2006.01", "2007.01", "2008.01", "2009.01", "2010.01", "2011.01", "2012.01", "2013.01", "2014.01", "2015.01", "2016.01")) +
						  facet_wrap(~UCR, scales = "free_y") +
						  labs(x = "Year.Week", y = "Crimes") +
						  ggtitle("Crimes Per Week By UCR") +
						  theme_light(base_size = 15) + 
						  theme(legend.position = "top",
								legend.key.size = unit(.25, "in"),
								axis.text.x = element_text(angle = 45, hjust = 1)) + 
						  guides(color = guide_legend(nrow = 1, override.aes = list(size = 10)))

plot_crimes_weekly_set4

DF = data.frame(crimes_weekly[UCR %in% levels(UCR)[25:26]])

plot_crimes_weekly_set5 = ggplot(data = DF, aes(x = YEAR.WEEK, y = CRIMES, color = UCR)) +
						  geom_point(na.rm = TRUE) +
						  scale_x_discrete(breaks = c("2006.01", "2007.01", "2008.01", "2009.01", "2010.01", "2011.01", "2012.01", "2013.01", "2014.01", "2015.01", "2016.01")) +
						  facet_wrap(~UCR, scales = "free_y") +
						  labs(x = "Year.Week", y = "Crimes") +
						  ggtitle("Crimes Per Week By UCR") +
						  theme_light(base_size = 15) + 
						  theme(legend.position = "top",
								legend.key.size = unit(.25, "in"),
								axis.text.x = element_text(angle = 45, hjust = 1)) + 
						  guides(color = guide_legend(nrow = 1, override.aes = list(size = 10)))

plot_crimes_weekly_set5

# lets put these results in lists to organize our workspace

crime_UCR_tables$crimes_weekly = crimes_weekly
crime_UCR_plots$plot_crimes_weekly_set1 = plot_crimes_weekly_set1
crime_UCR_plots$plot_crimes_weekly_set2 = plot_crimes_weekly_set2
crime_UCR_plots$plot_crimes_weekly_set3 = plot_crimes_weekly_set3
crime_UCR_plots$plot_crimes_weekly_set4 = plot_crimes_weekly_set4
crime_UCR_plots$plot_crimes_weekly_set5 = plot_crimes_weekly_set5

rm(crimes_weekly, plot_crimes_weekly_set1, plot_crimes_weekly_set2, plot_crimes_weekly_set3, plot_crimes_weekly_set4, plot_crimes_weekly_set5)

# lets create a table of crimes per month ------------------------------------------------------------------------------------------------------------------------------------

crimes_monthly = PHI[,.("CRIMES" = .N), by = .("UCR" = UCR_General, "YEAR" = Year, "MONTH" = Month)]
crimes_monthly[, YEAR.MONTH := factor(paste0(crimes_monthly[,YEAR], ".", crimes_monthly[,MONTH]),
							 levels = paste0(crimes_monthly[,YEAR], ".", crimes_monthly[,MONTH]))]
crimes_monthly[, c("YEAR", "MONTH") := NULL]
setcolorder(crimes_monthly, c("UCR", "YEAR.MONTH", "CRIMES"))

# lets plot crimes_monthly

DF = data.frame(crimes_monthly[UCR %in% levels(UCR)[1:6]])

plot_crimes_monthly_set1 = ggplot(data = DF, aes(x = YEAR.MONTH, y = CRIMES, color = UCR)) +
						   geom_point(na.rm = TRUE) +
						   scale_x_discrete(breaks = c("2006.Jan", "2007.Jan", "2008.Jan", "2009.Jan", "2010.Jan", "2011.Jan", "2012.Jan", "2013.Jan", "2014.Jan", "2015.Jan", "2016.Jan")) +
						   facet_wrap(~UCR, scales = "free_y") +
						   labs(x = "Year.Month", y = "Crimes") +
						   ggtitle("Crimes Per Month By UCR") +
						   theme_light(base_size = 15) + 
						   theme(legend.position = "top",
								 legend.key.size = unit(.25, "in"),
								 axis.text.x = element_text(angle = 45, hjust = 1)) + 
						   guides(color = guide_legend(nrow = 1, override.aes = list(size = 10)))

plot_crimes_monthly_set1

DF = data.frame(crimes_monthly[UCR %in% levels(UCR)[7:12]])

plot_crimes_monthly_set2 = ggplot(data = DF, aes(x = YEAR.MONTH, y = CRIMES, color = UCR)) +
						   geom_point(na.rm = TRUE) +
						   scale_x_discrete(breaks = c("2006.Jan", "2007.Jan", "2008.Jan", "2009.Jan", "2010.Jan", "2011.Jan", "2012.Jan", "2013.Jan", "2014.Jan", "2015.Jan", "2016.Jan")) +
						   facet_wrap(~UCR, scales = "free_y") +
						   labs(x = "Year.Month", y = "Crimes") +
						   ggtitle("Crimes Per Month By UCR") +
						   theme_light(base_size = 15) + 
						   theme(legend.position = "top",
								 legend.key.size = unit(.25, "in"),
								 axis.text.x = element_text(angle = 45, hjust = 1)) + 
						   guides(color = guide_legend(nrow = 1, override.aes = list(size = 10)))

plot_crimes_monthly_set2

DF = data.frame(crimes_monthly[UCR %in% levels(UCR)[13:18]])

plot_crimes_monthly_set3 = ggplot(data = DF, aes(x = YEAR.MONTH, y = CRIMES, color = UCR)) +
						   geom_point(na.rm = TRUE) +
						   scale_x_discrete(breaks = c("2006.Jan", "2007.Jan", "2008.Jan", "2009.Jan", "2010.Jan", "2011.Jan", "2012.Jan", "2013.Jan", "2014.Jan", "2015.Jan", "2016.Jan")) +
						   facet_wrap(~UCR, scales = "free_y") +
						   labs(x = "Year.Month", y = "Crimes") +
						   ggtitle("Crimes Per Month By UCR") +
						   theme_light(base_size = 15) + 
						   theme(legend.position = "top",
								 legend.key.size = unit(.25, "in"),
								 axis.text.x = element_text(angle = 45, hjust = 1)) + 
						   guides(color = guide_legend(nrow = 1, override.aes = list(size = 10)))

plot_crimes_monthly_set3

DF = data.frame(crimes_monthly[UCR %in% levels(UCR)[19:24]])

plot_crimes_monthly_set4 = ggplot(data = DF, aes(x = YEAR.MONTH, y = CRIMES, color = UCR)) +
						   geom_point(na.rm = TRUE) +
						   scale_x_discrete(breaks = c("2006.Jan", "2007.Jan", "2008.Jan", "2009.Jan", "2010.Jan", "2011.Jan", "2012.Jan", "2013.Jan", "2014.Jan", "2015.Jan", "2016.Jan")) +
						   facet_wrap(~UCR, scales = "free_y") +
						   labs(x = "Year.Month", y = "Crimes") +
						   ggtitle("Crimes Per Month By UCR") +
						   theme_light(base_size = 15) + 
						   theme(legend.position = "top",
								 legend.key.size = unit(.25, "in"),
								 axis.text.x = element_text(angle = 45, hjust = 1)) + 
						   guides(color = guide_legend(nrow = 1, override.aes = list(size = 10)))

plot_crimes_monthly_set4

DF = data.frame(crimes_monthly[UCR %in% levels(UCR)[25:26]])

plot_crimes_monthly_set5 = ggplot(data = DF, aes(x = YEAR.MONTH, y = CRIMES, color = UCR)) +
						   geom_point(na.rm = TRUE) +
						   scale_x_discrete(breaks = c("2006.Jan", "2007.Jan", "2008.Jan", "2009.Jan", "2010.Jan", "2011.Jan", "2012.Jan", "2013.Jan", "2014.Jan", "2015.Jan", "2016.Jan")) +
						   facet_wrap(~UCR, scales = "free_y") +
						   labs(x = "Year.Month", y = "Crimes") +
						   ggtitle("Crimes Per Month By UCR") +
						   theme_light(base_size = 15) + 
						   theme(legend.position = "top",
								 legend.key.size = unit(.25, "in"),
								 axis.text.x = element_text(angle = 45, hjust = 1)) + 
						   guides(color = guide_legend(nrow = 1, override.aes = list(size = 10)))

plot_crimes_monthly_set5

# lets put these results in lists to organize our workspace

crime_UCR_tables$crimes_monthly = crimes_monthly
crime_UCR_plots$plot_crimes_monthly_set1 = plot_crimes_monthly_set1
crime_UCR_plots$plot_crimes_monthly_set2 = plot_crimes_monthly_set2
crime_UCR_plots$plot_crimes_monthly_set3 = plot_crimes_monthly_set3
crime_UCR_plots$plot_crimes_monthly_set4 = plot_crimes_monthly_set4
crime_UCR_plots$plot_crimes_monthly_set5 = plot_crimes_monthly_set5

rm(crimes_monthly, plot_crimes_monthly_set1, plot_crimes_monthly_set2, plot_crimes_monthly_set3, plot_crimes_monthly_set4, plot_crimes_monthly_set5)

}

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ---- Crimes By Text_General_Code -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

{

# lets create a table of crimes per hour ------------------------------------------------------------------------------------------------------------------------------------

crimes_hourly = PHI[,.("CRIMES" = .N), by = .("OFFENSE" = Text_General_Code, "DATETIME" = Dispatch_Date_Time, "HOUR" = Hour)]

# lets plot crimes_hourly

DF = data.frame(crimes_hourly[OFFENSE %in% levels(OFFENSE)[1:6]])

plot_crimes_hourly_set1 = ggplot(data = DF, aes(x = DATETIME, y = CRIMES, color = OFFENSE)) +
						  geom_point(na.rm = TRUE) +
						  facet_wrap(~OFFENSE, scales = "free_y") +
						  labs(x = "Hour", y = "Crimes") +
						  ggtitle("Crimes Per Hour By OFFENSE") +
						  theme_light(base_size = 15) + 
						  theme(legend.position = "top",
								legend.key.size = unit(.25, "in")) + 
						  guides(color = guide_legend(nrow = 2, override.aes = list(size = 10)))

plot_crimes_hourly_set1

DF = data.frame(crimes_hourly[OFFENSE %in% levels(OFFENSE)[7:12]])

plot_crimes_hourly_set2 = ggplot(data = DF, aes(x = DATETIME, y = CRIMES, color = OFFENSE)) +
						  geom_point(na.rm = TRUE) +
						  facet_wrap(~OFFENSE, scales = "free_y") +
						  labs(x = "Hour", y = "Crimes") +
						  ggtitle("Crimes Per Hour By OFFENSE") +
						  theme_light(base_size = 15) + 
						  theme(legend.position = "top",
								legend.key.size = unit(.25, "in")) + 
						  guides(color = guide_legend(nrow = 2, override.aes = list(size = 10)))

plot_crimes_hourly_set2

DF = data.frame(crimes_hourly[OFFENSE %in% levels(OFFENSE)[13:18]])

plot_crimes_hourly_set3 = ggplot(data = DF, aes(x = DATETIME, y = CRIMES, color = OFFENSE)) +
						  geom_point(na.rm = TRUE) +
						  facet_wrap(~OFFENSE, scales = "free_y") +
						  labs(x = "Hour", y = "Crimes") +
						  ggtitle("Crimes Per Hour By OFFENSE") +
						  theme_light(base_size = 15) + 
						  theme(legend.position = "top",
								legend.key.size = unit(.25, "in")) + 
						  guides(color = guide_legend(nrow = 2, override.aes = list(size = 10)))

plot_crimes_hourly_set3

DF = data.frame(crimes_hourly[OFFENSE %in% levels(OFFENSE)[19:24]])

plot_crimes_hourly_set4 = ggplot(data = DF, aes(x = DATETIME, y = CRIMES, color = OFFENSE)) +
						  geom_point(na.rm = TRUE) +
						  facet_wrap(~OFFENSE, scales = "free_y") +
						  labs(x = "Hour", y = "Crimes") +
						  ggtitle("Crimes Per Hour By OFFENSE") +
						  theme_light(base_size = 15) + 
						  theme(legend.position = "top",
								legend.key.size = unit(.25, "in")) + 
						  guides(color = guide_legend(nrow = 2, override.aes = list(size = 10)))

plot_crimes_hourly_set4

DF = data.frame(crimes_hourly[OFFENSE %in% levels(OFFENSE)[25:30]])

plot_crimes_hourly_set5 = ggplot(data = DF, aes(x = DATETIME, y = CRIMES, color = OFFENSE)) +
						  geom_point(na.rm = TRUE) +
						  facet_wrap(~OFFENSE, scales = "free_y") +
						  labs(x = "Hour", y = "Crimes") +
						  ggtitle("Crimes Per Hour By OFFENSE") +
						  theme_light(base_size = 15) + 
						  theme(legend.position = "top",
								legend.key.size = unit(.25, "in")) + 
						  guides(color = guide_legend(nrow = 2, override.aes = list(size = 10)))

plot_crimes_hourly_set5

DF = data.frame(crimes_hourly[OFFENSE %in% levels(OFFENSE)[31:33]])

plot_crimes_hourly_set6 = ggplot(data = DF, aes(x = DATETIME, y = CRIMES, color = OFFENSE)) +
						  geom_point(na.rm = TRUE) +
						  facet_wrap(~OFFENSE, scales = "free_y") +
						  labs(x = "Hour", y = "Crimes") +
						  ggtitle("Crimes Per Hour By OFFENSE") +
						  theme_light(base_size = 15) + 
						  theme(legend.position = "top",
								legend.key.size = unit(.25, "in")) + 
						  guides(color = guide_legend(nrow = 2, override.aes = list(size = 10)))

plot_crimes_hourly_set6

# lets put these results in lists to organize our workspace

crime_OFFENSE_tables = list("crimes_hourly" = crimes_hourly)
crime_OFFENSE_plots = list("plot_crimes_hourly_set1" = plot_crimes_hourly_set1, 
					  "plot_crimes_hourly_set2" = plot_crimes_hourly_set2,
					  "plot_crimes_hourly_set3" = plot_crimes_hourly_set3,
					  "plot_crimes_hourly_set4" = plot_crimes_hourly_set4,
					  "plot_crimes_hourly_set5" = plot_crimes_hourly_set5,
					  "plot_crimes_hourly_set6" = plot_crimes_hourly_set6)

rm(crimes_hourly, plot_crimes_hourly_set1, plot_crimes_hourly_set2, plot_crimes_hourly_set3, plot_crimes_hourly_set4, plot_crimes_hourly_set5, plot_crimes_hourly_set6)

# lets create a table of crimes per day ------------------------------------------------------------------------------------------------------------------------------------

crimes_daily = PHI[,.("CRIMES" = .N), by = .("OFFENSE" = Text_General_Code, "DATE" = Dispatch_Date)]

# lets plot crimes_daily

DF = data.frame(crimes_daily[OFFENSE %in% levels(OFFENSE)[1:6]])

plot_crimes_daily_set1 = ggplot(data = DF, aes(x = DATE, y = CRIMES, color = OFFENSE)) +
						 geom_point(na.rm = TRUE) +
						 facet_wrap(~OFFENSE, scales = "free_y") +
						 labs(x = "Day", y = "Crimes") +
						 ggtitle("Crimes Per Day By OFFENSE") +
						 theme_light(base_size = 15) + 
						 theme(legend.position = "top",
						       legend.key.size = unit(.25, "in")) + 
						 guides(color = guide_legend(nrow = 2, override.aes = list(size = 10)))

plot_crimes_daily_set1

DF = data.frame(crimes_daily[OFFENSE %in% levels(OFFENSE)[7:12]])

plot_crimes_daily_set2 = ggplot(data = DF, aes(x = DATE, y = CRIMES, color = OFFENSE)) +
						 geom_point(na.rm = TRUE) +
						 facet_wrap(~OFFENSE, scales = "free_y") +
						 labs(x = "Day", y = "Crimes") +
						 ggtitle("Crimes Per Day By OFFENSE") +
						 theme_light(base_size = 15) + 
						 theme(legend.position = "top",
						       legend.key.size = unit(.25, "in")) + 
						 guides(color = guide_legend(nrow = 2, override.aes = list(size = 10)))

plot_crimes_daily_set2

DF = data.frame(crimes_daily[OFFENSE %in% levels(OFFENSE)[13:18]])

plot_crimes_daily_set3 = ggplot(data = DF, aes(x = DATE, y = CRIMES, color = OFFENSE)) +
						 geom_point(na.rm = TRUE) +
						 facet_wrap(~OFFENSE, scales = "free_y") +
						 labs(x = "Day", y = "Crimes") +
						 ggtitle("Crimes Per Day By OFFENSE") +
						 theme_light(base_size = 15) + 
						 theme(legend.position = "top",
						       legend.key.size = unit(.25, "in")) + 
						 guides(color = guide_legend(nrow = 2, override.aes = list(size = 10)))

plot_crimes_daily_set3

DF = data.frame(crimes_daily[OFFENSE %in% levels(OFFENSE)[19:24]])

plot_crimes_daily_set4 = ggplot(data = DF, aes(x = DATE, y = CRIMES, color = OFFENSE)) +
						 geom_point(na.rm = TRUE) +
						 facet_wrap(~OFFENSE, scales = "free_y") +
						 labs(x = "Day", y = "Crimes") +
						 ggtitle("Crimes Per Day By OFFENSE") +
						 theme_light(base_size = 15) + 
						 theme(legend.position = "top",
						       legend.key.size = unit(.25, "in")) + 
						 guides(color = guide_legend(nrow = 2, override.aes = list(size = 10)))

plot_crimes_daily_set4

DF = data.frame(crimes_daily[OFFENSE %in% levels(OFFENSE)[25:30]])

plot_crimes_daily_set5 = ggplot(data = DF, aes(x = DATE, y = CRIMES, color = OFFENSE)) +
						 geom_point(na.rm = TRUE) +
						 facet_wrap(~OFFENSE, scales = "free_y") +
						 labs(x = "Day", y = "Crimes") +
						 ggtitle("Crimes Per Day By OFFENSE") +
						 theme_light(base_size = 15) + 
						 theme(legend.position = "top",
						       legend.key.size = unit(.25, "in")) + 
						 guides(color = guide_legend(nrow = 2, override.aes = list(size = 10)))

plot_crimes_daily_set5

DF = data.frame(crimes_daily[OFFENSE %in% levels(OFFENSE)[31:33]])

plot_crimes_daily_set6 = ggplot(data = DF, aes(x = DATE, y = CRIMES, color = OFFENSE)) +
						 geom_point(na.rm = TRUE) +
						 facet_wrap(~OFFENSE, scales = "free_y") +
						 labs(x = "Day", y = "Crimes") +
						 ggtitle("Crimes Per Day By OFFENSE") +
						 theme_light(base_size = 15) + 
						 theme(legend.position = "top",
						       legend.key.size = unit(.25, "in")) + 
						 guides(color = guide_legend(nrow = 2, override.aes = list(size = 10)))

plot_crimes_daily_set6

# lets put these results in lists to organize our workspace

crime_OFFENSE_tables$crimes_daily = crimes_daily
crime_OFFENSE_plots$plot_crimes_daily_set1 = plot_crimes_daily_set1
crime_OFFENSE_plots$plot_crimes_daily_set2 = plot_crimes_daily_set2
crime_OFFENSE_plots$plot_crimes_daily_set3 = plot_crimes_daily_set3
crime_OFFENSE_plots$plot_crimes_daily_set4 = plot_crimes_daily_set4
crime_OFFENSE_plots$plot_crimes_daily_set5 = plot_crimes_daily_set5
crime_OFFENSE_plots$plot_crimes_daily_set6 = plot_crimes_daily_set6

rm(crimes_daily, plot_crimes_daily_set1, plot_crimes_daily_set2, plot_crimes_daily_set3, plot_crimes_daily_set4, plot_crimes_daily_set5, plot_crimes_daily_set6)

# lets create a table of crimes per week ------------------------------------------------------------------------------------------------------------------------------------

crimes_weekly = PHI[,.("CRIMES" = .N), by = .("OFFENSE" = Text_General_Code, "YEAR" = Year, "WEEK" = Week)]
crimes_weekly[, YEAR.WEEK := factor(paste0(crimes_weekly[,YEAR], ".", crimes_weekly[,WEEK]),
							 levels = paste0(crimes_weekly[,YEAR], ".", crimes_weekly[,WEEK]))]
crimes_weekly[, c("YEAR", "WEEK") := NULL]
setcolorder(crimes_weekly, c("OFFENSE", "YEAR.WEEK", "CRIMES"))

# lets plot crimes_weekly

DF = data.frame(crimes_weekly[OFFENSE %in% levels(OFFENSE)[1:6]])

plot_crimes_weekly_set1 = ggplot(data = DF, aes(x = YEAR.WEEK, y = CRIMES, color = OFFENSE)) +
						  geom_point(na.rm = TRUE) +
						  scale_x_discrete(breaks = c("2006.01", "2007.01", "2008.01", "2009.01", "2010.01", "2011.01", "2012.01", "2013.01", "2014.01", "2015.01", "2016.01")) +
						  facet_wrap(~OFFENSE, scales = "free_y") +
						  labs(x = "Year.Week", y = "Crimes") +
						  ggtitle("Crimes Per Week By OFFENSE") +
						  theme_light(base_size = 15) + 
						  theme(legend.position = "top",
								legend.key.size = unit(.25, "in"),
								axis.text.x = element_text(angle = 45, hjust = 1)) + 
						  guides(color = guide_legend(nrow = 2, override.aes = list(size = 10)))

plot_crimes_weekly_set1

DF = data.frame(crimes_weekly[OFFENSE %in% levels(OFFENSE)[7:12]])

plot_crimes_weekly_set2 = ggplot(data = DF, aes(x = YEAR.WEEK, y = CRIMES, color = OFFENSE)) +
						  geom_point(na.rm = TRUE) +
						  scale_x_discrete(breaks = c("2006.01", "2007.01", "2008.01", "2009.01", "2010.01", "2011.01", "2012.01", "2013.01", "2014.01", "2015.01", "2016.01")) +
						  facet_wrap(~OFFENSE, scales = "free_y") +
						  labs(x = "Year.Week", y = "Crimes") +
						  ggtitle("Crimes Per Week By OFFENSE") +
						  theme_light(base_size = 15) + 
						  theme(legend.position = "top",
								legend.key.size = unit(.25, "in"),
								axis.text.x = element_text(angle = 45, hjust = 1)) + 
						  guides(color = guide_legend(nrow = 2, override.aes = list(size = 10)))

plot_crimes_weekly_set2

DF = data.frame(crimes_weekly[OFFENSE %in% levels(OFFENSE)[13:18]])

plot_crimes_weekly_set3 = ggplot(data = DF, aes(x = YEAR.WEEK, y = CRIMES, color = OFFENSE)) +
						  geom_point(na.rm = TRUE) +
						  scale_x_discrete(breaks = c("2006.01", "2007.01", "2008.01", "2009.01", "2010.01", "2011.01", "2012.01", "2013.01", "2014.01", "2015.01", "2016.01")) +
						  facet_wrap(~OFFENSE, scales = "free_y") +
						  labs(x = "Year.Week", y = "Crimes") +
						  ggtitle("Crimes Per Week By OFFENSE") +
						  theme_light(base_size = 15) + 
						  theme(legend.position = "top",
								legend.key.size = unit(.25, "in"),
								axis.text.x = element_text(angle = 45, hjust = 1)) + 
						  guides(color = guide_legend(nrow = 2, override.aes = list(size = 10)))

plot_crimes_weekly_set3

DF = data.frame(crimes_weekly[OFFENSE %in% levels(OFFENSE)[19:24]])

plot_crimes_weekly_set4 = ggplot(data = DF, aes(x = YEAR.WEEK, y = CRIMES, color = OFFENSE)) +
						  geom_point(na.rm = TRUE) +
						  scale_x_discrete(breaks = c("2006.01", "2007.01", "2008.01", "2009.01", "2010.01", "2011.01", "2012.01", "2013.01", "2014.01", "2015.01", "2016.01")) +
						  facet_wrap(~OFFENSE, scales = "free_y") +
						  labs(x = "Year.Week", y = "Crimes") +
						  ggtitle("Crimes Per Week By OFFENSE") +
						  theme_light(base_size = 15) + 
						  theme(legend.position = "top",
								legend.key.size = unit(.25, "in"),
								axis.text.x = element_text(angle = 45, hjust = 1)) + 
						  guides(color = guide_legend(nrow = 2, override.aes = list(size = 10)))

plot_crimes_weekly_set4

DF = data.frame(crimes_weekly[OFFENSE %in% levels(OFFENSE)[25:30]])

plot_crimes_weekly_set5 = ggplot(data = DF, aes(x = YEAR.WEEK, y = CRIMES, color = OFFENSE)) +
						  geom_point(na.rm = TRUE) +
						  scale_x_discrete(breaks = c("2006.01", "2007.01", "2008.01", "2009.01", "2010.01", "2011.01", "2012.01", "2013.01", "2014.01", "2015.01", "2016.01")) +
						  facet_wrap(~OFFENSE, scales = "free_y") +
						  labs(x = "Year.Week", y = "Crimes") +
						  ggtitle("Crimes Per Week By OFFENSE") +
						  theme_light(base_size = 15) + 
						  theme(legend.position = "top",
								legend.key.size = unit(.25, "in"),
								axis.text.x = element_text(angle = 45, hjust = 1)) + 
						  guides(color = guide_legend(nrow = 2, override.aes = list(size = 10)))

plot_crimes_weekly_set5

DF = data.frame(crimes_weekly[OFFENSE %in% levels(OFFENSE)[31:33]])

plot_crimes_weekly_set6 = ggplot(data = DF, aes(x = YEAR.WEEK, y = CRIMES, color = OFFENSE)) +
						  geom_point(na.rm = TRUE) +
						  scale_x_discrete(breaks = c("2006.01", "2007.01", "2008.01", "2009.01", "2010.01", "2011.01", "2012.01", "2013.01", "2014.01", "2015.01", "2016.01")) +
						  facet_wrap(~OFFENSE, scales = "free_y") +
						  labs(x = "Year.Week", y = "Crimes") +
						  ggtitle("Crimes Per Week By OFFENSE") +
						  theme_light(base_size = 15) + 
						  theme(legend.position = "top",
								legend.key.size = unit(.25, "in"),
								axis.text.x = element_text(angle = 45, hjust = 1)) + 
						  guides(color = guide_legend(nrow = 2, override.aes = list(size = 10)))

plot_crimes_weekly_set6

# lets put these results in lists to organize our workspace

crime_OFFENSE_tables$crimes_weekly = crimes_weekly
crime_OFFENSE_plots$plot_crimes_weekly_set1 = plot_crimes_weekly_set1
crime_OFFENSE_plots$plot_crimes_weekly_set2 = plot_crimes_weekly_set2
crime_OFFENSE_plots$plot_crimes_weekly_set3 = plot_crimes_weekly_set3
crime_OFFENSE_plots$plot_crimes_weekly_set4 = plot_crimes_weekly_set4
crime_OFFENSE_plots$plot_crimes_weekly_set5 = plot_crimes_weekly_set5
crime_OFFENSE_plots$plot_crimes_weekly_set6 = plot_crimes_weekly_set6

rm(crimes_weekly, plot_crimes_weekly_set1, plot_crimes_weekly_set2, plot_crimes_weekly_set3, plot_crimes_weekly_set4, plot_crimes_weekly_set5, plot_crimes_weekly_set6)

# lets create a table of crimes per month ------------------------------------------------------------------------------------------------------------------------------------

crimes_monthly = PHI[,.("CRIMES" = .N), by = .("OFFENSE" = Text_General_Code, "YEAR" = Year, "MONTH" = Month)]
crimes_monthly[, YEAR.MONTH := factor(paste0(crimes_monthly[,YEAR], ".", crimes_monthly[,MONTH]),
							 levels = paste0(crimes_monthly[,YEAR], ".", crimes_monthly[,MONTH]))]
crimes_monthly[, c("YEAR", "MONTH") := NULL]
setcolorder(crimes_monthly, c("OFFENSE", "YEAR.MONTH", "CRIMES"))

# lets plot crimes_monthly

DF = data.frame(crimes_monthly[OFFENSE %in% levels(OFFENSE)[1:6]])

plot_crimes_monthly_set1 = ggplot(data = DF, aes(x = YEAR.MONTH, y = CRIMES, color = OFFENSE)) +
						   geom_point(na.rm = TRUE) +
						   scale_x_discrete(breaks = c("2006.Jan", "2007.Jan", "2008.Jan", "2009.Jan", "2010.Jan", "2011.Jan", "2012.Jan", "2013.Jan", "2014.Jan", "2015.Jan", "2016.Jan")) +
						   facet_wrap(~OFFENSE, scales = "free_y") +
						   labs(x = "Year.Month", y = "Crimes") +
						   ggtitle("Crimes Per Month By OFFENSE") +
						   theme_light(base_size = 15) + 
						   theme(legend.position = "top",
								 legend.key.size = unit(.25, "in"),
								 axis.text.x = element_text(angle = 45, hjust = 1)) + 
						   guides(color = guide_legend(nrow = 2, override.aes = list(size = 10)))

plot_crimes_monthly_set1

DF = data.frame(crimes_monthly[OFFENSE %in% levels(OFFENSE)[7:12]])

plot_crimes_monthly_set2 = ggplot(data = DF, aes(x = YEAR.MONTH, y = CRIMES, color = OFFENSE)) +
						   geom_point(na.rm = TRUE) +
						   scale_x_discrete(breaks = c("2006.Jan", "2007.Jan", "2008.Jan", "2009.Jan", "2010.Jan", "2011.Jan", "2012.Jan", "2013.Jan", "2014.Jan", "2015.Jan", "2016.Jan")) +
						   facet_wrap(~OFFENSE, scales = "free_y") +
						   labs(x = "Year.Month", y = "Crimes") +
						   ggtitle("Crimes Per Month By OFFENSE") +
						   theme_light(base_size = 15) + 
						   theme(legend.position = "top",
								 legend.key.size = unit(.25, "in"),
								 axis.text.x = element_text(angle = 45, hjust = 1)) + 
						   guides(color = guide_legend(nrow = 2, override.aes = list(size = 10)))

plot_crimes_monthly_set2

DF = data.frame(crimes_monthly[OFFENSE %in% levels(OFFENSE)[13:18]])

plot_crimes_monthly_set3 = ggplot(data = DF, aes(x = YEAR.MONTH, y = CRIMES, color = OFFENSE)) +
						   geom_point(na.rm = TRUE) +
						   scale_x_discrete(breaks = c("2006.Jan", "2007.Jan", "2008.Jan", "2009.Jan", "2010.Jan", "2011.Jan", "2012.Jan", "2013.Jan", "2014.Jan", "2015.Jan", "2016.Jan")) +
						   facet_wrap(~OFFENSE, scales = "free_y") +
						   labs(x = "Year.Month", y = "Crimes") +
						   ggtitle("Crimes Per Month By OFFENSE") +
						   theme_light(base_size = 15) + 
						   theme(legend.position = "top",
								 legend.key.size = unit(.25, "in"),
								 axis.text.x = element_text(angle = 45, hjust = 1)) + 
						   guides(color = guide_legend(nrow = 2, override.aes = list(size = 10)))

plot_crimes_monthly_set3

DF = data.frame(crimes_monthly[OFFENSE %in% levels(OFFENSE)[19:24]])

plot_crimes_monthly_set4 = ggplot(data = DF, aes(x = YEAR.MONTH, y = CRIMES, color = OFFENSE)) +
						   geom_point(na.rm = TRUE) +
						   scale_x_discrete(breaks = c("2006.Jan", "2007.Jan", "2008.Jan", "2009.Jan", "2010.Jan", "2011.Jan", "2012.Jan", "2013.Jan", "2014.Jan", "2015.Jan", "2016.Jan")) +
						   facet_wrap(~OFFENSE, scales = "free_y") +
						   labs(x = "Year.Month", y = "Crimes") +
						   ggtitle("Crimes Per Month By OFFENSE") +
						   theme_light(base_size = 15) + 
						   theme(legend.position = "top",
								 legend.key.size = unit(.25, "in"),
								 axis.text.x = element_text(angle = 45, hjust = 1)) + 
						   guides(color = guide_legend(nrow = 2, override.aes = list(size = 10)))

plot_crimes_monthly_set4

DF = data.frame(crimes_monthly[OFFENSE %in% levels(OFFENSE)[25:30]])

plot_crimes_monthly_set5 = ggplot(data = DF, aes(x = YEAR.MONTH, y = CRIMES, color = OFFENSE)) +
						   geom_point(na.rm = TRUE) +
						   scale_x_discrete(breaks = c("2006.Jan", "2007.Jan", "2008.Jan", "2009.Jan", "2010.Jan", "2011.Jan", "2012.Jan", "2013.Jan", "2014.Jan", "2015.Jan", "2016.Jan")) +
						   facet_wrap(~OFFENSE, scales = "free_y") +
						   labs(x = "Year.Month", y = "Crimes") +
						   ggtitle("Crimes Per Month By OFFENSE") +
						   theme_light(base_size = 15) + 
						   theme(legend.position = "top",
								 legend.key.size = unit(.25, "in"),
								 axis.text.x = element_text(angle = 45, hjust = 1)) + 
						   guides(color = guide_legend(nrow = 2, override.aes = list(size = 10)))

plot_crimes_monthly_set5

DF = data.frame(crimes_monthly[OFFENSE %in% levels(OFFENSE)[31:33]])

plot_crimes_monthly_set6 = ggplot(data = DF, aes(x = YEAR.MONTH, y = CRIMES, color = OFFENSE)) +
						   geom_point(na.rm = TRUE) +
						   scale_x_discrete(breaks = c("2006.Jan", "2007.Jan", "2008.Jan", "2009.Jan", "2010.Jan", "2011.Jan", "2012.Jan", "2013.Jan", "2014.Jan", "2015.Jan", "2016.Jan")) +
						   facet_wrap(~OFFENSE, scales = "free_y") +
						   labs(x = "Year.Month", y = "Crimes") +
						   ggtitle("Crimes Per Month By OFFENSE") +
						   theme_light(base_size = 15) + 
						   theme(legend.position = "top",
								 legend.key.size = unit(.25, "in"),
								 axis.text.x = element_text(angle = 45, hjust = 1)) + 
						   guides(color = guide_legend(nrow = 2, override.aes = list(size = 10)))

plot_crimes_monthly_set6

# lets put these results in lists to organize our workspace

crime_OFFENSE_tables$crimes_monthly = crimes_monthly
crime_OFFENSE_plots$plot_crimes_monthly_set1 = plot_crimes_monthly_set1
crime_OFFENSE_plots$plot_crimes_monthly_set2 = plot_crimes_monthly_set2
crime_OFFENSE_plots$plot_crimes_monthly_set3 = plot_crimes_monthly_set3
crime_OFFENSE_plots$plot_crimes_monthly_set4 = plot_crimes_monthly_set4
crime_OFFENSE_plots$plot_crimes_monthly_set5 = plot_crimes_monthly_set5
crime_OFFENSE_plots$plot_crimes_monthly_set6 = plot_crimes_monthly_set6

rm(crimes_monthly, plot_crimes_monthly_set1, plot_crimes_monthly_set2, plot_crimes_monthly_set3, plot_crimes_monthly_set4, plot_crimes_monthly_set5, plot_crimes_monthly_set6)

}

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ---- Crimes By Police_Districts -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

{

# lets create a table of crimes per hour ------------------------------------------------------------------------------------------------------------------------------------

crimes_hourly = PHI[,.("CRIMES" = .N), by = .("District" = Police_Districts, "DATETIME" = Dispatch_Date_Time, "HOUR" = Hour)]

# lets plot crimes_hourly

DF = data.frame(crimes_hourly[District %in% levels(District)[1:6]])

plot_crimes_hourly_set1 = ggplot(data = DF, aes(x = DATETIME, y = CRIMES, color = District)) +
						  geom_point(na.rm = TRUE) +
						  facet_wrap(~District, scales = "free_y") +
						  labs(x = "Hour", y = "Crimes") +
						  ggtitle("Crimes Per Hour By District") +
						  theme_light(base_size = 15) + 
						  theme(legend.position = "top",
								legend.key.size = unit(.25, "in")) + 
						  guides(color = guide_legend(nrow = 1, override.aes = list(size = 10)))

plot_crimes_hourly_set1

DF = data.frame(crimes_hourly[District %in% levels(District)[7:12]])

plot_crimes_hourly_set2 = ggplot(data = DF, aes(x = DATETIME, y = CRIMES, color = District)) +
						  geom_point(na.rm = TRUE) +
						  facet_wrap(~District, scales = "free_y") +
						  labs(x = "Hour", y = "Crimes") +
						  ggtitle("Crimes Per Hour By District") +
						  theme_light(base_size = 15) + 
						  theme(legend.position = "top",
								legend.key.size = unit(.25, "in")) + 
						  guides(color = guide_legend(nrow = 1, override.aes = list(size = 10)))

plot_crimes_hourly_set2

DF = data.frame(crimes_hourly[District %in% levels(District)[13:18]])

plot_crimes_hourly_set3 = ggplot(data = DF, aes(x = DATETIME, y = CRIMES, color = District)) +
						  geom_point(na.rm = TRUE) +
						  facet_wrap(~District, scales = "free_y") +
						  labs(x = "Hour", y = "Crimes") +
						  ggtitle("Crimes Per Hour By District") +
						  theme_light(base_size = 15) + 
						  theme(legend.position = "top",
								legend.key.size = unit(.25, "in")) + 
						  guides(color = guide_legend(nrow = 1, override.aes = list(size = 10)))

plot_crimes_hourly_set3

DF = data.frame(crimes_hourly[District %in% levels(District)[19:22]])

plot_crimes_hourly_set4 = ggplot(data = DF, aes(x = DATETIME, y = CRIMES, color = District)) +
						  geom_point(na.rm = TRUE) +
						  facet_wrap(~District, scales = "free_y") +
						  labs(x = "Hour", y = "Crimes") +
						  ggtitle("Crimes Per Hour By District") +
						  theme_light(base_size = 15) + 
						  theme(legend.position = "top",
								legend.key.size = unit(.25, "in")) + 
						  guides(color = guide_legend(nrow = 1, override.aes = list(size = 10)))

plot_crimes_hourly_set4

# lets put these results in lists to organize our workspace

crime_District_tables = list("crimes_hourly" = crimes_hourly)
crime_District_plots = list("plot_crimes_hourly_set1" = plot_crimes_hourly_set1, 
					  "plot_crimes_hourly_set2" = plot_crimes_hourly_set2,
					  "plot_crimes_hourly_set3" = plot_crimes_hourly_set3,
					  "plot_crimes_hourly_set4" = plot_crimes_hourly_set4)

rm(crimes_hourly, plot_crimes_hourly_set1, plot_crimes_hourly_set2, plot_crimes_hourly_set3, plot_crimes_hourly_set4)

# lets create a table of crimes per day ------------------------------------------------------------------------------------------------------------------------------------

crimes_daily = PHI[,.("CRIMES" = .N), by = .("District" = Police_Districts, "DATE" = Dispatch_Date)]

# lets plot crimes_daily

DF = data.frame(crimes_daily[District %in% levels(District)[1:6]])

plot_crimes_daily_set1 = ggplot(data = DF, aes(x = DATE, y = CRIMES, color = District)) +
						 geom_point(na.rm = TRUE) +
						 facet_wrap(~District, scales = "free_y") +
						 labs(x = "Day", y = "Crimes") +
						 ggtitle("Crimes Per Day By District") +
						 theme_light(base_size = 15) + 
						 theme(legend.position = "top",
						       legend.key.size = unit(.25, "in")) + 
						 guides(color = guide_legend(nrow = 1, override.aes = list(size = 10)))

plot_crimes_daily_set1

DF = data.frame(crimes_daily[District %in% levels(District)[7:12]])

plot_crimes_daily_set2 = ggplot(data = DF, aes(x = DATE, y = CRIMES, color = District)) +
						 geom_point(na.rm = TRUE) +
						 facet_wrap(~District, scales = "free_y") +
						 labs(x = "Day", y = "Crimes") +
						 ggtitle("Crimes Per Day By District") +
						 theme_light(base_size = 15) + 
						 theme(legend.position = "top",
						       legend.key.size = unit(.25, "in")) + 
						 guides(color = guide_legend(nrow = 1, override.aes = list(size = 10)))

plot_crimes_daily_set2

DF = data.frame(crimes_daily[District %in% levels(District)[13:18]])

plot_crimes_daily_set3 = ggplot(data = DF, aes(x = DATE, y = CRIMES, color = District)) +
						 geom_point(na.rm = TRUE) +
						 facet_wrap(~District, scales = "free_y") +
						 labs(x = "Day", y = "Crimes") +
						 ggtitle("Crimes Per Day By District") +
						 theme_light(base_size = 15) + 
						 theme(legend.position = "top",
						       legend.key.size = unit(.25, "in")) + 
						 guides(color = guide_legend(nrow = 1, override.aes = list(size = 10)))

plot_crimes_daily_set3

DF = data.frame(crimes_daily[District %in% levels(District)[19:22]])

plot_crimes_daily_set4 = ggplot(data = DF, aes(x = DATE, y = CRIMES, color = District)) +
						 geom_point(na.rm = TRUE) +
						 facet_wrap(~District, scales = "free_y") +
						 labs(x = "Day", y = "Crimes") +
						 ggtitle("Crimes Per Day By District") +
						 theme_light(base_size = 15) + 
						 theme(legend.position = "top",
						       legend.key.size = unit(.25, "in")) + 
						 guides(color = guide_legend(nrow = 1, override.aes = list(size = 10)))

plot_crimes_daily_set4

# lets put these results in lists to organize our workspace

crime_District_tables$crimes_daily = crimes_daily
crime_District_plots$plot_crimes_daily_set1 = plot_crimes_daily_set1
crime_District_plots$plot_crimes_daily_set2 = plot_crimes_daily_set2
crime_District_plots$plot_crimes_daily_set3 = plot_crimes_daily_set3
crime_District_plots$plot_crimes_daily_set4 = plot_crimes_daily_set4

rm(crimes_daily, plot_crimes_daily_set1, plot_crimes_daily_set2, plot_crimes_daily_set3, plot_crimes_daily_set4)

# lets create a table of crimes per week ------------------------------------------------------------------------------------------------------------------------------------

crimes_weekly = PHI[,.("CRIMES" = .N), by = .("District" = Police_Districts, "YEAR" = Year, "WEEK" = Week)]
crimes_weekly[, YEAR.WEEK := factor(paste0(crimes_weekly[,YEAR], ".", crimes_weekly[,WEEK]),
							 levels = paste0(crimes_weekly[,YEAR], ".", crimes_weekly[,WEEK]))]
crimes_weekly[, c("YEAR", "WEEK") := NULL]
setcolorder(crimes_weekly, c("District", "YEAR.WEEK", "CRIMES"))

# lets plot crimes_weekly

DF = data.frame(crimes_weekly[District %in% levels(District)[1:6]])

plot_crimes_weekly_set1 = ggplot(data = DF, aes(x = YEAR.WEEK, y = CRIMES, color = District)) +
						  geom_point(na.rm = TRUE) +
						  scale_x_discrete(breaks = c("2006.01", "2007.01", "2008.01", "2009.01", "2010.01", "2011.01", "2012.01", "2013.01", "2014.01", "2015.01", "2016.01")) +
						  facet_wrap(~District, scales = "free_y") +
						  labs(x = "Year.Week", y = "Crimes") +
						  ggtitle("Crimes Per Week By District") +
						  theme_light(base_size = 15) + 
						  theme(legend.position = "top",
								legend.key.size = unit(.25, "in"),
								axis.text.x = element_text(angle = 45, hjust = 1)) + 
						  guides(color = guide_legend(nrow = 1, override.aes = list(size = 10)))

plot_crimes_weekly_set1

DF = data.frame(crimes_weekly[District %in% levels(District)[7:12]])

plot_crimes_weekly_set2 = ggplot(data = DF, aes(x = YEAR.WEEK, y = CRIMES, color = District)) +
						  geom_point(na.rm = TRUE) +
						  scale_x_discrete(breaks = c("2006.01", "2007.01", "2008.01", "2009.01", "2010.01", "2011.01", "2012.01", "2013.01", "2014.01", "2015.01", "2016.01")) +
						  facet_wrap(~District, scales = "free_y") +
						  labs(x = "Year.Week", y = "Crimes") +
						  ggtitle("Crimes Per Week By District") +
						  theme_light(base_size = 15) + 
						  theme(legend.position = "top",
								legend.key.size = unit(.25, "in"),
								axis.text.x = element_text(angle = 45, hjust = 1)) + 
						  guides(color = guide_legend(nrow = 1, override.aes = list(size = 10)))

plot_crimes_weekly_set2

DF = data.frame(crimes_weekly[District %in% levels(District)[13:18]])

plot_crimes_weekly_set3 = ggplot(data = DF, aes(x = YEAR.WEEK, y = CRIMES, color = District)) +
						  geom_point(na.rm = TRUE) +
						  scale_x_discrete(breaks = c("2006.01", "2007.01", "2008.01", "2009.01", "2010.01", "2011.01", "2012.01", "2013.01", "2014.01", "2015.01", "2016.01")) +
						  facet_wrap(~District, scales = "free_y") +
						  labs(x = "Year.Week", y = "Crimes") +
						  ggtitle("Crimes Per Week By District") +
						  theme_light(base_size = 15) + 
						  theme(legend.position = "top",
								legend.key.size = unit(.25, "in"),
								axis.text.x = element_text(angle = 45, hjust = 1)) + 
						  guides(color = guide_legend(nrow = 1, override.aes = list(size = 10)))

plot_crimes_weekly_set3

DF = data.frame(crimes_weekly[District %in% levels(District)[19:22]])

plot_crimes_weekly_set4 = ggplot(data = DF, aes(x = YEAR.WEEK, y = CRIMES, color = District)) +
						  geom_point(na.rm = TRUE) +
						  scale_x_discrete(breaks = c("2006.01", "2007.01", "2008.01", "2009.01", "2010.01", "2011.01", "2012.01", "2013.01", "2014.01", "2015.01", "2016.01")) +
						  facet_wrap(~District, scales = "free_y") +
						  labs(x = "Year.Week", y = "Crimes") +
						  ggtitle("Crimes Per Week By District") +
						  theme_light(base_size = 15) + 
						  theme(legend.position = "top",
								legend.key.size = unit(.25, "in"),
								axis.text.x = element_text(angle = 45, hjust = 1)) + 
						  guides(color = guide_legend(nrow = 1, override.aes = list(size = 10)))

plot_crimes_weekly_set4

# lets put these results in lists to organize our workspace

crime_District_tables$crimes_weekly = crimes_weekly
crime_District_plots$plot_crimes_weekly_set1 = plot_crimes_weekly_set1
crime_District_plots$plot_crimes_weekly_set2 = plot_crimes_weekly_set2
crime_District_plots$plot_crimes_weekly_set3 = plot_crimes_weekly_set3
crime_District_plots$plot_crimes_weekly_set4 = plot_crimes_weekly_set4

rm(crimes_weekly, plot_crimes_weekly_set1, plot_crimes_weekly_set2, plot_crimes_weekly_set3, plot_crimes_weekly_set4)

# lets create a table of crimes per month ------------------------------------------------------------------------------------------------------------------------------------

crimes_monthly = PHI[,.("CRIMES" = .N), by = .("District" = Police_Districts, "YEAR" = Year, "MONTH" = Month)]
crimes_monthly[, YEAR.MONTH := factor(paste0(crimes_monthly[,YEAR], ".", crimes_monthly[,MONTH]),
							 levels = paste0(crimes_monthly[,YEAR], ".", crimes_monthly[,MONTH]))]
crimes_monthly[, c("YEAR", "MONTH") := NULL]
setcolorder(crimes_monthly, c("District", "YEAR.MONTH", "CRIMES"))

# lets plot crimes_monthly

DF = data.frame(crimes_monthly[District %in% levels(District)[1:6]])

plot_crimes_monthly_set1 = ggplot(data = DF, aes(x = YEAR.MONTH, y = CRIMES, color = District)) +
						   geom_point(na.rm = TRUE) +
						   scale_x_discrete(breaks = c("2006.Jan", "2007.Jan", "2008.Jan", "2009.Jan", "2010.Jan", "2011.Jan", "2012.Jan", "2013.Jan", "2014.Jan", "2015.Jan", "2016.Jan")) +
						   facet_wrap(~District, scales = "free_y") +
						   labs(x = "Year.Month", y = "Crimes") +
						   ggtitle("Crimes Per Month By District") +
						   theme_light(base_size = 15) + 
						   theme(legend.position = "top",
								 legend.key.size = unit(.25, "in"),
								 axis.text.x = element_text(angle = 45, hjust = 1)) + 
						   guides(color = guide_legend(nrow = 1, override.aes = list(size = 10)))

plot_crimes_monthly_set1

DF = data.frame(crimes_monthly[District %in% levels(District)[7:12]])

plot_crimes_monthly_set2 = ggplot(data = DF, aes(x = YEAR.MONTH, y = CRIMES, color = District)) +
						   geom_point(na.rm = TRUE) +
						   scale_x_discrete(breaks = c("2006.Jan", "2007.Jan", "2008.Jan", "2009.Jan", "2010.Jan", "2011.Jan", "2012.Jan", "2013.Jan", "2014.Jan", "2015.Jan", "2016.Jan")) +
						   facet_wrap(~District, scales = "free_y") +
						   labs(x = "Year.Month", y = "Crimes") +
						   ggtitle("Crimes Per Month By District") +
						   theme_light(base_size = 15) + 
						   theme(legend.position = "top",
								 legend.key.size = unit(.25, "in"),
								 axis.text.x = element_text(angle = 45, hjust = 1)) + 
						   guides(color = guide_legend(nrow = 1, override.aes = list(size = 10)))

plot_crimes_monthly_set2

DF = data.frame(crimes_monthly[District %in% levels(District)[13:18]])

plot_crimes_monthly_set3 = ggplot(data = DF, aes(x = YEAR.MONTH, y = CRIMES, color = District)) +
						   geom_point(na.rm = TRUE) +
						   scale_x_discrete(breaks = c("2006.Jan", "2007.Jan", "2008.Jan", "2009.Jan", "2010.Jan", "2011.Jan", "2012.Jan", "2013.Jan", "2014.Jan", "2015.Jan", "2016.Jan")) +
						   facet_wrap(~District, scales = "free_y") +
						   labs(x = "Year.Month", y = "Crimes") +
						   ggtitle("Crimes Per Month By District") +
						   theme_light(base_size = 15) + 
						   theme(legend.position = "top",
								 legend.key.size = unit(.25, "in"),
								 axis.text.x = element_text(angle = 45, hjust = 1)) + 
						   guides(color = guide_legend(nrow = 1, override.aes = list(size = 10)))

plot_crimes_monthly_set3

DF = data.frame(crimes_monthly[District %in% levels(District)[19:22]])

plot_crimes_monthly_set4 = ggplot(data = DF, aes(x = YEAR.MONTH, y = CRIMES, color = District)) +
						   geom_point(na.rm = TRUE) +
						   scale_x_discrete(breaks = c("2006.Jan", "2007.Jan", "2008.Jan", "2009.Jan", "2010.Jan", "2011.Jan", "2012.Jan", "2013.Jan", "2014.Jan", "2015.Jan", "2016.Jan")) +
						   facet_wrap(~District, scales = "free_y") +
						   labs(x = "Year.Month", y = "Crimes") +
						   ggtitle("Crimes Per Month By District") +
						   theme_light(base_size = 15) + 
						   theme(legend.position = "top",
								 legend.key.size = unit(.25, "in"),
								 axis.text.x = element_text(angle = 45, hjust = 1)) + 
						   guides(color = guide_legend(nrow = 1, override.aes = list(size = 10)))

plot_crimes_monthly_set4

# lets put these results in lists to organize our workspace

crime_District_tables$crimes_monthly = crimes_monthly
crime_District_plots$plot_crimes_monthly_set1 = plot_crimes_monthly_set1
crime_District_plots$plot_crimes_monthly_set2 = plot_crimes_monthly_set2
crime_District_plots$plot_crimes_monthly_set3 = plot_crimes_monthly_set3
crime_District_plots$plot_crimes_monthly_set4 = plot_crimes_monthly_set4

rm(crimes_monthly, plot_crimes_monthly_set1, plot_crimes_monthly_set2, plot_crimes_monthly_set3, plot_crimes_monthly_set4)

}

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ---- Crimes By Area1 -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

{

# lets create a table of crimes per hour ------------------------------------------------------------------------------------------------------------------------------------

crimes_hourly = PHI[,.("CRIMES" = .N), by = .("Area1" = Area1, "DATETIME" = Dispatch_Date_Time, "HOUR" = Hour)]

# lets plot crimes_hourly

DF = data.frame(crimes_hourly[Area1 %in% levels(Area1)[1]])

plot_crimes_hourly_CP = ggplot(data = DF, aes(x = DATETIME, y = CRIMES)) +
						  geom_point(na.rm = TRUE) +
						  labs(x = "Hour", y = "Crimes") +
						  ggtitle("Crimes Per Hour\nCentral Philly") +
						  theme_light(base_size = 15)

plot_crimes_hourly_CP

DF = data.frame(crimes_hourly[Area1 %in% levels(Area1)[2]])

plot_crimes_hourly_WP = ggplot(data = DF, aes(x = DATETIME, y = CRIMES)) +
						  geom_point(na.rm = TRUE) +
						  labs(x = "Hour", y = "Crimes") +
						  ggtitle("Crimes Per Hour\nWest Philly") +
						  theme_light(base_size = 15)

plot_crimes_hourly_WP

DF = data.frame(crimes_hourly[Area1 %in% levels(Area1)[3]])

plot_crimes_hourly_NP = ggplot(data = DF, aes(x = DATETIME, y = CRIMES)) +
						  geom_point(na.rm = TRUE) +
						  labs(x = "Hour", y = "Crimes") +
						  ggtitle("Crimes Per Hour\nNorth Philly") +
						  theme_light(base_size = 15)

plot_crimes_hourly_NP

# lets put these results in lists to organize our workspace

crime_Area1_tables = list("crimes_hourly" = crimes_hourly)
crime_Area1_plots = list("plot_crimes_hourly_CP" = plot_crimes_hourly_CP, 
					  "plot_crimes_hourly_WP" = plot_crimes_hourly_WP,
					  "plot_crimes_hourly_NP" = plot_crimes_hourly_NP)

rm(crimes_hourly, plot_crimes_hourly_CP, plot_crimes_hourly_WP, plot_crimes_hourly_NP)

# lets create a table of crimes per day ------------------------------------------------------------------------------------------------------------------------------------

crimes_daily = PHI[,.("CRIMES" = .N), by = .("Area1" = Area1, "DATE" = Dispatch_Date)]

# lets plot crimes_daily

DF = data.frame(crimes_daily[Area1 %in% levels(Area1)[1]])

plot_crimes_daily_CP = ggplot(data = DF, aes(x = DATE, y = CRIMES)) +
						 geom_point(na.rm = TRUE) +
						 labs(x = "Day", y = "Crimes") +
						 ggtitle("Crimes Per Day\nCentral Philly") +
						 theme_light(base_size = 15)

plot_crimes_daily_CP

DF = data.frame(crimes_daily[Area1 %in% levels(Area1)[2]])

plot_crimes_daily_WP = ggplot(data = DF, aes(x = DATE, y = CRIMES)) +
						 geom_point(na.rm = TRUE) +
						 labs(x = "Day", y = "Crimes") +
						 ggtitle("Crimes Per Day\nWest Philly") +
						 theme_light(base_size = 15)

plot_crimes_daily_WP

DF = data.frame(crimes_daily[Area1 %in% levels(Area1)[3]])

plot_crimes_daily_NP = ggplot(data = DF, aes(x = DATE, y = CRIMES)) +
						 geom_point(na.rm = TRUE) +
						 labs(x = "Day", y = "Crimes") +
						 ggtitle("Crimes Per Day\nNorth Philly") +
						 theme_light(base_size = 15)

plot_crimes_daily_NP

# lets put these results in lists to organize our workspace

crime_Area1_tables$crimes_daily = crimes_daily
crime_Area1_plots$plot_crimes_daily_CP = plot_crimes_daily_CP
crime_Area1_plots$plot_crimes_daily_WP = plot_crimes_daily_WP
crime_Area1_plots$plot_crimes_daily_NP = plot_crimes_daily_NP

rm(crimes_daily, plot_crimes_daily_CP, plot_crimes_daily_WP, plot_crimes_daily_NP)

# lets create a table of crimes per week ------------------------------------------------------------------------------------------------------------------------------------

crimes_weekly = PHI[,.("CRIMES" = .N), by = .("Area1" = Area1, "YEAR" = Year, "WEEK" = Week)]
crimes_weekly[, YEAR.WEEK := factor(paste0(crimes_weekly[,YEAR], ".", crimes_weekly[,WEEK]),
							 levels = paste0(crimes_weekly[,YEAR], ".", crimes_weekly[,WEEK]))]
crimes_weekly[, c("YEAR", "WEEK") := NULL]
setcolorder(crimes_weekly, c("Area1", "YEAR.WEEK", "CRIMES"))

# lets plot crimes_weekly

DF = data.frame(crimes_weekly[Area1 %in% levels(Area1)[1]])

plot_crimes_weekly_CP = ggplot(data = DF, aes(x = YEAR.WEEK, y = CRIMES)) +
						  geom_point(na.rm = TRUE) +
						  scale_x_discrete(breaks = c("2006.01", "2007.01", "2008.01", "2009.01", "2010.01", "2011.01", "2012.01", "2013.01", "2014.01", "2015.01", "2016.01")) +
						  labs(x = "Year.Week", y = "Crimes") +
						  ggtitle("Crimes Per Week\nCentral Philly") +
						  theme_light(base_size = 15) + 
						  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_crimes_weekly_CP

DF = data.frame(crimes_weekly[Area1 %in% levels(Area1)[2]])

plot_crimes_weekly_WP = ggplot(data = DF, aes(x = YEAR.WEEK, y = CRIMES)) +
						  geom_point(na.rm = TRUE) +
						  scale_x_discrete(breaks = c("2006.01", "2007.01", "2008.01", "2009.01", "2010.01", "2011.01", "2012.01", "2013.01", "2014.01", "2015.01", "2016.01")) +
						  labs(x = "Year.Week", y = "Crimes") +
						  ggtitle("Crimes Per Week\nWest Philly") +
						  theme_light(base_size = 15) + 
						  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_crimes_weekly_WP

DF = data.frame(crimes_weekly[Area1 %in% levels(Area1)[3]])

plot_crimes_weekly_NP = ggplot(data = DF, aes(x = YEAR.WEEK, y = CRIMES)) +
						  geom_point(na.rm = TRUE) +
						  scale_x_discrete(breaks = c("2006.01", "2007.01", "2008.01", "2009.01", "2010.01", "2011.01", "2012.01", "2013.01", "2014.01", "2015.01", "2016.01")) +
						  labs(x = "Year.Week", y = "Crimes") +
						  ggtitle("Crimes Per Week\nNorth Philly") +
						  theme_light(base_size = 15) + 
						  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_crimes_weekly_NP

# lets put these results in lists to organize our workspace

crime_Area1_tables$crimes_weekly = crimes_weekly
crime_Area1_plots$plot_crimes_weekly_CP = plot_crimes_weekly_CP
crime_Area1_plots$plot_crimes_weekly_WP = plot_crimes_weekly_WP
crime_Area1_plots$plot_crimes_weekly_NP = plot_crimes_weekly_NP

rm(crimes_weekly, plot_crimes_weekly_CP, plot_crimes_weekly_WP, plot_crimes_weekly_NP)

# lets create a table of crimes per month ------------------------------------------------------------------------------------------------------------------------------------

crimes_monthly = PHI[,.("CRIMES" = .N), by = .("Area1" = Area1, "YEAR" = Year, "MONTH" = Month)]
crimes_monthly[, YEAR.MONTH := factor(paste0(crimes_monthly[,YEAR], ".", crimes_monthly[,MONTH]),
							 levels = paste0(crimes_monthly[,YEAR], ".", crimes_monthly[,MONTH]))]
crimes_monthly[, c("YEAR", "MONTH") := NULL]
setcolorder(crimes_monthly, c("Area1", "YEAR.MONTH", "CRIMES"))

# lets plot crimes_monthly

DF = data.frame(crimes_monthly[Area1 %in% levels(Area1)[1]])

plot_crimes_monthly_CP = ggplot(data = DF, aes(x = YEAR.MONTH, y = CRIMES)) +
						   geom_point(na.rm = TRUE) +
						   scale_x_discrete(breaks = c("2006.Jan", "2007.Jan", "2008.Jan", "2009.Jan", "2010.Jan", "2011.Jan", "2012.Jan", "2013.Jan", "2014.Jan", "2015.Jan", "2016.Jan")) +
						   labs(x = "Year.Month", y = "Crimes") +
						   ggtitle("Crimes Per Month\nCentral Philly") +
						   theme_light(base_size = 15) + 
						   theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_crimes_monthly_CP

DF = data.frame(crimes_monthly[Area1 %in% levels(Area1)[2]])

plot_crimes_monthly_WP = ggplot(data = DF, aes(x = YEAR.MONTH, y = CRIMES)) +
						   geom_point(na.rm = TRUE) +
						   scale_x_discrete(breaks = c("2006.Jan", "2007.Jan", "2008.Jan", "2009.Jan", "2010.Jan", "2011.Jan", "2012.Jan", "2013.Jan", "2014.Jan", "2015.Jan", "2016.Jan")) +
						   labs(x = "Year.Month", y = "Crimes") +
						   ggtitle("Crimes Per Month\nWest Philly") +
						   theme_light(base_size = 15) + 
						   theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_crimes_monthly_WP

DF = data.frame(crimes_monthly[Area1 %in% levels(Area1)[3]])

plot_crimes_monthly_NP = ggplot(data = DF, aes(x = YEAR.MONTH, y = CRIMES)) +
						   geom_point(na.rm = TRUE) +
						   scale_x_discrete(breaks = c("2006.Jan", "2007.Jan", "2008.Jan", "2009.Jan", "2010.Jan", "2011.Jan", "2012.Jan", "2013.Jan", "2014.Jan", "2015.Jan", "2016.Jan")) +
						   labs(x = "Year.Month", y = "Crimes") +
						   ggtitle("Crimes Per Month\nNorth Philly") +
						   theme_light(base_size = 15) + 
						   theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_crimes_monthly_NP

# lets put these results in lists to organize our workspace

crime_Area1_tables$crimes_monthly = crimes_monthly
crime_Area1_plots$plot_crimes_monthly_CP = plot_crimes_monthly_CP
crime_Area1_plots$plot_crimes_monthly_WP = plot_crimes_monthly_WP
crime_Area1_plots$plot_crimes_monthly_NP = plot_crimes_monthly_NP

rm(crimes_monthly, plot_crimes_monthly_CP, plot_crimes_monthly_WP, plot_crimes_monthly_NP)

}

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ---- Crimes By Area2 -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

{

# lets create a table of crimes per hour ------------------------------------------------------------------------------------------------------------------------------------

crimes_hourly = PHI[,.("CRIMES" = .N), by = .("Area2" = Area2, "DATETIME" = Dispatch_Date_Time, "HOUR" = Hour)]

# lets plot crimes_hourly

DF = data.frame(crimes_hourly[Area2 %in% levels(Area2)[1]])

plot_crimes_hourly_A1 = ggplot(data = DF, aes(x = DATETIME, y = CRIMES)) +
						  geom_point(na.rm = TRUE) +
						  labs(x = "Hour", y = "Crimes") +
						  ggtitle("Crimes Per Hour\nA1") +
						  theme_light(base_size = 15)

plot_crimes_hourly_A1

DF = data.frame(crimes_hourly[Area2 %in% levels(Area2)[2]])

plot_crimes_hourly_A2 = ggplot(data = DF, aes(x = DATETIME, y = CRIMES)) +
						  geom_point(na.rm = TRUE) +
						  labs(x = "Hour", y = "Crimes") +
						  ggtitle("Crimes Per Hour\nA2") +
						  theme_light(base_size = 15)

plot_crimes_hourly_A2

# lets put these results in lists to organize our workspace

crime_Area2_tables = list("crimes_hourly" = crimes_hourly)
crime_Area2_plots = list("plot_crimes_hourly_A1" = plot_crimes_hourly_A1, 
					  "plot_crimes_hourly_A2" = plot_crimes_hourly_A2)

rm(crimes_hourly, plot_crimes_hourly_A1, plot_crimes_hourly_A2)

# lets create a table of crimes per day ------------------------------------------------------------------------------------------------------------------------------------

crimes_daily = PHI[,.("CRIMES" = .N), by = .("Area2" = Area2, "DATE" = Dispatch_Date)]

# lets plot crimes_daily

DF = data.frame(crimes_daily[Area2 %in% levels(Area2)[1]])

plot_crimes_daily_A1 = ggplot(data = DF, aes(x = DATE, y = CRIMES)) +
						 geom_point(na.rm = TRUE) +
						 labs(x = "Day", y = "Crimes") +
						 ggtitle("Crimes Per Day\nA1") +
						 theme_light(base_size = 15)

plot_crimes_daily_A1

DF = data.frame(crimes_daily[Area2 %in% levels(Area2)[2]])

plot_crimes_daily_A2 = ggplot(data = DF, aes(x = DATE, y = CRIMES)) +
						 geom_point(na.rm = TRUE) +
						 labs(x = "Day", y = "Crimes") +
						 ggtitle("Crimes Per Day\nA2") +
						 theme_light(base_size = 15)

plot_crimes_daily_A2

# lets put these results in lists to organize our workspace

crime_Area2_tables$crimes_daily = crimes_daily
crime_Area2_plots$plot_crimes_daily_A1 = plot_crimes_daily_A1
crime_Area2_plots$plot_crimes_daily_A2 = plot_crimes_daily_A2

rm(crimes_daily, plot_crimes_daily_A1, plot_crimes_daily_A2)

# lets create a table of crimes per week ------------------------------------------------------------------------------------------------------------------------------------

crimes_weekly = PHI[,.("CRIMES" = .N), by = .("Area2" = Area2, "YEAR" = Year, "WEEK" = Week)]
crimes_weekly[, YEAR.WEEK := factor(paste0(crimes_weekly[,YEAR], ".", crimes_weekly[,WEEK]),
							 levels = paste0(crimes_weekly[,YEAR], ".", crimes_weekly[,WEEK]))]
crimes_weekly[, c("YEAR", "WEEK") := NULL]
setcolorder(crimes_weekly, c("Area2", "YEAR.WEEK", "CRIMES"))

# lets plot crimes_weekly

DF = data.frame(crimes_weekly[Area2 %in% levels(Area2)[1]])

plot_crimes_weekly_A1 = ggplot(data = DF, aes(x = YEAR.WEEK, y = CRIMES)) +
						  geom_point(na.rm = TRUE) +
						  scale_x_discrete(breaks = c("2006.01", "2007.01", "2008.01", "2009.01", "2010.01", "2011.01", "2012.01", "2013.01", "2014.01", "2015.01", "2016.01")) +
						  labs(x = "Year.Week", y = "Crimes") +
						  ggtitle("Crimes Per Week\nA1") +
						  theme_light(base_size = 15) + 
						  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_crimes_weekly_A1

DF = data.frame(crimes_weekly[Area2 %in% levels(Area2)[2]])

plot_crimes_weekly_A2 = ggplot(data = DF, aes(x = YEAR.WEEK, y = CRIMES)) +
						  geom_point(na.rm = TRUE) +
						  scale_x_discrete(breaks = c("2006.01", "2007.01", "2008.01", "2009.01", "2010.01", "2011.01", "2012.01", "2013.01", "2014.01", "2015.01", "2016.01")) +
						  labs(x = "Year.Week", y = "Crimes") +
						  ggtitle("Crimes Per Week\nA2") +
						  theme_light(base_size = 15) + 
						  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_crimes_weekly_A2

# lets put these results in lists to organize our workspace

crime_Area2_tables$crimes_weekly = crimes_weekly
crime_Area2_plots$plot_crimes_weekly_A1 = plot_crimes_weekly_A1
crime_Area2_plots$plot_crimes_weekly_A2 = plot_crimes_weekly_A2

rm(crimes_weekly, plot_crimes_weekly_A1, plot_crimes_weekly_A2)

# lets create a table of crimes per month ------------------------------------------------------------------------------------------------------------------------------------

crimes_monthly = PHI[,.("CRIMES" = .N), by = .("Area2" = Area2, "YEAR" = Year, "MONTH" = Month)]
crimes_monthly[, YEAR.MONTH := factor(paste0(crimes_monthly[,YEAR], ".", crimes_monthly[,MONTH]),
							 levels = paste0(crimes_monthly[,YEAR], ".", crimes_monthly[,MONTH]))]
crimes_monthly[, c("YEAR", "MONTH") := NULL]
setcolorder(crimes_monthly, c("Area2", "YEAR.MONTH", "CRIMES"))

# lets plot crimes_monthly

DF = data.frame(crimes_monthly[Area2 %in% levels(Area2)[1]])

plot_crimes_monthly_A1 = ggplot(data = DF, aes(x = YEAR.MONTH, y = CRIMES)) +
						   geom_point(na.rm = TRUE) +
						   scale_x_discrete(breaks = c("2006.Jan", "2007.Jan", "2008.Jan", "2009.Jan", "2010.Jan", "2011.Jan", "2012.Jan", "2013.Jan", "2014.Jan", "2015.Jan", "2016.Jan")) +
						   labs(x = "Year.Month", y = "Crimes") +
						   ggtitle("Crimes Per Month\nA1") +
						   theme_light(base_size = 15) + 
						   theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_crimes_monthly_A1

DF = data.frame(crimes_monthly[Area2 %in% levels(Area2)[2]])

plot_crimes_monthly_A2 = ggplot(data = DF, aes(x = YEAR.MONTH, y = CRIMES)) +
						   geom_point(na.rm = TRUE) +
						   scale_x_discrete(breaks = c("2006.Jan", "2007.Jan", "2008.Jan", "2009.Jan", "2010.Jan", "2011.Jan", "2012.Jan", "2013.Jan", "2014.Jan", "2015.Jan", "2016.Jan")) +
						   labs(x = "Year.Month", y = "Crimes") +
						   ggtitle("Crimes Per Month\nA2") +
						   theme_light(base_size = 15) + 
						   theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_crimes_monthly_A2

# lets put these results in lists to organize our workspace

crime_Area2_tables$crimes_monthly = crimes_monthly
crime_Area2_plots$plot_crimes_monthly_A1 = plot_crimes_monthly_A1
crime_Area2_plots$plot_crimes_monthly_A2 = plot_crimes_monthly_A2

rm(crimes_monthly, plot_crimes_monthly_A1, plot_crimes_monthly_A2)

# lets bundle up all the lists we've created thus far to organize the workspace

tables = list("crime_tables" = crime_tables,
			  "crime_DC_tables" = crime_DC_tables,
			  "crime_PSA_tables" = crime_PSA_tables,
			  "crime_UCR_tables" = crime_UCR_tables,
			  "crime_OFFENSE_tables" = crime_OFFENSE_tables,
			  "crime_District_tables" = crime_District_tables,
			  "crime_Area1_tables" = crime_Area1_tables,
			  "crime_Area2_tables" = crime_Area2_tables)

plots = list("crime_plots" = crime_plots,
			 "crime_DC_plots" = crime_DC_plots,
			 "crime_PSA_plots" = crime_PSA_plots,
			 "crime_UCR_plots" = crime_UCR_plots,
			 "crime_OFFENSE_plots" = crime_OFFENSE_plots,
			 "crime_District_plots" = crime_District_plots,
			 "crime_Area1_plots" = crime_Area1_plots,
			 "crime_Area2_plots" = crime_Area2_plots,
			 "heatmaps" = heatmaps)

rm(crime_tables, crime_DC_tables, crime_PSA_tables, crime_UCR_tables, crime_OFFENSE_tables, crime_District_tables, crime_Area1_tables, crime_Area2_tables,
   crime_plots, crime_DC_plots, crime_PSA_plots, crime_UCR_plots, crime_OFFENSE_plots, crime_District_plots, crime_Area1_plots, crime_Area2_plots, heatmaps)

}

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ---- Finding Trend & Seasonality ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

{

# ---- daily time series ----------------------------------------------------------------------------------------------------------------

# data

plots$crime_plots$plot_crimes_daily_color

# external variables
	
	# Area1

plots$crime_Area1_plots$plot_crimes_daily_CP
plots$crime_Area1_plots$plot_crimes_daily_WP
plots$crime_Area1_plots$plot_crimes_daily_NP

# ---- weekly time series ---------------------------------------------------------------------------------------------------------------

# data

plots$crime_plots$plot_crimes_weekly

# external variables

	# DC
	
plots$crime_DC_plots$plot_crimes_weekly_set1
plots$crime_DC_plots$plot_crimes_weekly_set2
plots$crime_DC_plots$plot_crimes_weekly_set3
plots$crime_DC_plots$plot_crimes_weekly_set4
plots$crime_DC_plots$plot_crimes_weekly_set5

	# PSA

plots$crime_PSA_plots$plot_crimes_weekly_set1
plots$crime_PSA_plots$plot_crimes_weekly_set2
plots$crime_PSA_plots$plot_crimes_weekly_set3
plots$crime_PSA_plots$plot_crimes_weekly_set4
plots$crime_PSA_plots$plot_crimes_weekly_set5

	# UCR

plots$crime_UCR_plots$plot_crimes_weekly_set1
plots$crime_UCR_plots$plot_crimes_weekly_set2
plots$crime_UCR_plots$plot_crimes_weekly_set3
plots$crime_UCR_plots$plot_crimes_weekly_set4
plots$crime_UCR_plots$plot_crimes_weekly_set5

	# OFFENSE

plots$crime_OFFENSE_plots$plot_crimes_weekly_set1
plots$crime_OFFENSE_plots$plot_crimes_weekly_set2
plots$crime_OFFENSE_plots$plot_crimes_weekly_set3
plots$crime_OFFENSE_plots$plot_crimes_weekly_set4
plots$crime_OFFENSE_plots$plot_crimes_weekly_set5
plots$crime_OFFENSE_plots$plot_crimes_weekly_set6

	# District

plots$crime_District_plots$plot_crimes_weekly_set1
plots$crime_District_plots$plot_crimes_weekly_set2
plots$crime_District_plots$plot_crimes_weekly_set3
plots$crime_District_plots$plot_crimes_weekly_set4

	# Area1

plots$crime_Area1_plots$plot_crimes_weekly_CP
plots$crime_Area1_plots$plot_crimes_weekly_WP
plots$crime_Area1_plots$plot_crimes_weekly_NP

#  --- lets try to explain all those really low values in plots$crime_plots$plot_crimes_weekly ------------------------------------------

tables$crime_tables$crimes_weekly[CRIMES <= 2000]

# most of these points correspond to the start/end of a year, these weeks aren't full weeks sometimes, hence the lower values
# there is one point at the end becuase the last day in this dataset is a day in the middle of a week so this is an incomplete week
# lets compute the total number of days for each week in the data

tables$crime_tables$crimes_weekly[, DAYS := sapply(1:NROW(tables$crime_tables$crimes_weekly), function(i)
																							  max(PHI[Year == substr(tables$crime_tables$crimes_weekly[i,YEAR.WEEK], 1, 4) & 
																									  Week == substr(tables$crime_tables$crimes_weekly[i,YEAR.WEEK], 6, 7)][,Dispatch_Date]) - 
																							  min(PHI[Year == substr(tables$crime_tables$crimes_weekly[i,YEAR.WEEK], 1, 4) & 
																									  Week == substr(tables$crime_tables$crimes_weekly[i,YEAR.WEEK], 6, 7)][,Dispatch_Date]) + 1)]

# lets verify the lack of days in these low vlaues 

tables$crime_tables$crimes_weekly[CRIMES <= 2000]

# lets re-plot crimes_weekly

tables$crime_tables$crimes_weekly[, DAYS := factor(DAYS, levels = as.character(1:7))]

DF = data.frame(tables$crime_tables$crimes_weekly)

plot_crimes_weekly_color_days = ggplot(data = DF, aes(x = YEAR.WEEK, y = CRIMES, color = DAYS)) +
								geom_point(na.rm = TRUE) +
								scale_x_discrete(breaks = c("2006.01", "2007.01", "2008.01", "2009.01", "2010.01", "2011.01", "2012.01", "2013.01", "2014.01", "2015.01", "2016.01")) +
								scale_color_manual(values = colorRampPalette(brewer.pal(n = 9, name = "YlOrRd"))(7)) +		
								labs(x = "Year.Week", y = "Crimes") +
								ggtitle("Crimes Per Week") +
								theme_dark(base_size = 20) +
								theme(legend.position = "top",
									  legend.key.size = unit(.5, "in"),
									  axis.text.x = element_text(angle = 45, hjust = 1)) + 
								guides(color = guide_legend(nrow = 1, override.aes = list(size = 10)))

plot_crimes_weekly_color_days

plots$crime_plots$plot_crimes_weekly_color_days = plot_crimes_weekly_color_days

rm(plot_crimes_weekly_color_days)

# lets update all weekly tables with the number of days in each week

setcolorder(tables$crime_tables$crimes_weekly, c(1,2,4,3))

DT = data.table(tables$crime_tables$crimes_weekly[,.(YEAR.WEEK, DAYS)])
setkey(DT, YEAR.WEEK)

setkey(tables$crime_DC_tables$crimes_weekly, YEAR.WEEK)
tables$crime_DC_tables$crimes_weekly = tables$crime_DC_tables$crimes_weekly[DT, nomatch = 0]
setcolorder(tables$crime_DC_tables$crimes_weekly, c(1,2,4,3))

setkey(tables$crime_PSA_tables$crimes_weekly, YEAR.WEEK)
tables$crime_PSA_tables$crimes_weekly = tables$crime_PSA_tables$crimes_weekly[DT, nomatch = 0]
setcolorder(tables$crime_PSA_tables$crimes_weekly, c(1,2,4,3))

setkey(tables$crime_UCR_tables$crimes_weekly, YEAR.WEEK)
tables$crime_UCR_tables$crimes_weekly = tables$crime_UCR_tables$crimes_weekly[DT, nomatch = 0]
setcolorder(tables$crime_UCR_tables$crimes_weekly, c(1,2,4,3))

setkey(tables$crime_OFFENSE_tables$crimes_weekly, YEAR.WEEK)
tables$crime_OFFENSE_tables$crimes_weekly = tables$crime_OFFENSE_tables$crimes_weekly[DT, nomatch = 0]
setcolorder(tables$crime_OFFENSE_tables$crimes_weekly, c(1,2,4,3))

setkey(tables$crime_District_tables$crimes_weekly, YEAR.WEEK)
tables$crime_District_tables$crimes_weekly = tables$crime_District_tables$crimes_weekly[DT, nomatch = 0]
setcolorder(tables$crime_District_tables$crimes_weekly, c(1,2,4,3))

setkey(tables$crime_Area1_tables$crimes_weekly, YEAR.WEEK)
tables$crime_Area1_tables$crimes_weekly = tables$crime_Area1_tables$crimes_weekly[DT, nomatch = 0]
setcolorder(tables$crime_Area1_tables$crimes_weekly, c(1,2,4,3))

setkey(tables$crime_Area2_tables$crimes_weekly, YEAR.WEEK)
tables$crime_Area2_tables$crimes_weekly = tables$crime_Area2_tables$crimes_weekly[DT, nomatch = 0]
setcolorder(tables$crime_Area2_tables$crimes_weekly, c(1,2,4,3))

rm(DT)

}

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ---- Smoothing the Data ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

{

# ---- moving average -------------------------------------------------------------------------------------------------------------------------

# plot a smoother for daily crimes

crimes_daily_MA = data.table(tables$crime_tables$crimes_daily)
crimes_daily_MA[, MA.30 := ma(CRIMES, 30)]

DF = data.frame(crimes_daily_MA)

plot_crimes_daily_MA = ggplot(data = DF, aes(x = DATE, y = CRIMES)) +
					   geom_point(na.rm = TRUE) +
					   geom_line(aes(x = DATE, y = MA.30), color = "cornflowerblue", na.rm = TRUE, size = 1.5) +
				   	   labs(x = "Day", y = "Crimes") +
				   	   ggtitle("Crimes Per Day\nTwo-Sided Moving Average of 30 Days") +
				   	   theme_light(base_size = 20)

plot_crimes_daily_MA

# plot a smoother for weekly crimes

crimes_weekly_MA = data.table(tables$crime_tables$crimes_weekly)
crimes_weekly_MA[, MA.5 := ma(CRIMES, 5)]

DF = data.frame(crimes_weekly_MA)

plot_crimes_weekly_MA = ggplot(data = DF, aes(x = YEAR.WEEK, y = CRIMES, group = 1)) +
						geom_point(na.rm = TRUE) +
						geom_line(aes(x = YEAR.WEEK, y = MA.5), color = "cornflowerblue", na.rm = TRUE, size = 1.5) +
						scale_x_discrete(breaks = c("2006.01", "2007.01", "2008.01", "2009.01", "2010.01", "2011.01", "2012.01", "2013.01", "2014.01", "2015.01", "2016.01")) +
						labs(x = "Year.Week", y = "Crimes") +
						ggtitle("Crimes Per Week\nTwo-Sided Moving Average of 5 Weeks") +
						theme_light(base_size = 20) + 
						theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_crimes_weekly_MA

# store the results in lists

smooth_plots = list("plot_crimes_daily_MA" = plot_crimes_daily_MA,
					"plot_crimes_weekly_MA" = plot_crimes_weekly_MA)

smooth_tables = list("crimes_daily_MA" = crimes_daily_MA,
					 "crimes_weekly_MA" = crimes_weekly_MA)

rm(crimes_daily_MA, crimes_weekly_MA, plot_crimes_daily_MA, plot_crimes_weekly_MA)

# ---- Tukey's Running Median -------------------------------------------------------------------------------------------------------------------------

# plot a smoother for daily crimes

crimes_daily_TRM = data.table(tables$crime_tables$crimes_daily)
crimes_daily_TRM[, TRM := as.numeric(smooth(CRIMES))]

DF = data.frame(crimes_daily_TRM)

plot_crimes_daily_TRM = ggplot(data = DF, aes(x = DATE, y = CRIMES)) +
						geom_point(na.rm = TRUE) +
						geom_line(aes(x = DATE, y = TRM), color = "cornflowerblue", na.rm = TRUE, size = 1) +
						labs(x = "Day", y = "Crimes") +
						ggtitle("Crimes Per Day\nTukey's Running Median") +
						theme_light(base_size = 20)

plot_crimes_daily_TRM

# plot a smoother for weekly crimes

crimes_weekly_TRM = data.table(tables$crime_tables$crimes_weekly)
crimes_weekly_TRM[, TRM := as.numeric(smooth(CRIMES))]

DF = data.frame(crimes_weekly_TRM)

plot_crimes_weekly_TRM = ggplot(data = DF, aes(x = YEAR.WEEK, y = CRIMES, group = 1)) +
						 geom_point(na.rm = TRUE) +
						 geom_line(aes(x = YEAR.WEEK, y = TRM), color = "cornflowerblue", na.rm = TRUE, size = 1) +
						 scale_x_discrete(breaks = c("2006.01", "2007.01", "2008.01", "2009.01", "2010.01", "2011.01", "2012.01", "2013.01", "2014.01", "2015.01", "2016.01")) +
						 labs(x = "Year.Week", y = "Crimes") +
						 ggtitle("Crimes Per Week\nTukey's Running Median") +
						 theme_light(base_size = 20) + 
						 theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_crimes_weekly_TRM

# store the results in lists

smooth_plots$plot_crimes_daily_TRM = plot_crimes_daily_TRM
smooth_plots$plot_crimes_weekly_TRM = plot_crimes_weekly_TRM
smooth_tables$crimes_daily_TRM = crimes_daily_TRM
smooth_tables$crimes_weekly_TRM = crimes_weekly_TRM

rm(plot_crimes_daily_TRM, plot_crimes_weekly_TRM, crimes_daily_TRM, crimes_weekly_TRM)

# ---- LOWESS (Locally Weighted Scatterplot Smoothing) -------------------------------------------------------------------------------------------------------------------------

# plot a smoother for daily crimes

crimes_daily_LOWESS = data.table(tables$crime_tables$crimes_daily)
crimes_daily_LOWESS[, LOWESS := lowess(CRIMES, f = 1/50)$y]

DF = data.frame(crimes_daily_LOWESS)

plot_crimes_daily_LOWESS = ggplot(data = DF, aes(x = DATE, y = CRIMES)) +
					   geom_point(na.rm = TRUE) +
					   geom_line(aes(x = DATE, y = LOWESS), color = "cornflowerblue", na.rm = TRUE, size = 1.5) +
				   	   labs(x = "Day", y = "Crimes") +
				   	   ggtitle("Crimes Per Day\nLocally Weighted Scatterplot Smoothing") +
				   	   theme_light(base_size = 20)

plot_crimes_daily_LOWESS

# plot a smoother for weekly crimes

crimes_weekly_LOWESS = data.table(tables$crime_tables$crimes_weekly)
crimes_weekly_LOWESS[, LOWESS := lowess(CRIMES, f = 1/50)$y]

DF = data.frame(crimes_weekly_LOWESS)

plot_crimes_weekly_LOWESS = ggplot(data = DF, aes(x = YEAR.WEEK, y = CRIMES, group = 1)) +
						geom_point(na.rm = TRUE) +
						geom_line(aes(x = YEAR.WEEK, y = LOWESS), color = "cornflowerblue", na.rm = TRUE, size = 1.5) +
						scale_x_discrete(breaks = c("2006.01", "2007.01", "2008.01", "2009.01", "2010.01", "2011.01", "2012.01", "2013.01", "2014.01", "2015.01", "2016.01")) +
						labs(x = "Year.Week", y = "Crimes") +
						ggtitle("Crimes Per Week\nLocally Weighted Scatterplot Smoothing") +
						theme_light(base_size = 20) + 
						theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_crimes_weekly_LOWESS

# store the results in lists

smooth_plots$plot_crimes_daily_LOWESS = plot_crimes_daily_LOWESS
smooth_plots$plot_crimes_weekly_LOWESS = plot_crimes_weekly_LOWESS
smooth_tables$crimes_daily_LOWESS = crimes_daily_LOWESS
smooth_tables$crimes_weekly_LOWESS = crimes_weekly_LOWESS

plots$smooth_plots = smooth_plots
tables$smooth_tables = smooth_tables

rm(plot_crimes_daily_LOWESS, plot_crimes_weekly_LOWESS, crimes_daily_LOWESS, crimes_weekly_LOWESS, smooth_plots, smooth_tables)

}

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ---- ACF Plots ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

{

plot_crimes_daily_acf = ggacf(x = tables$crime_tables$crimes_daily$CRIMES, 
							  n = 365, 
							  main = "ACF Plot of Crimes Per Day")

plot_crimes_daily_acf

plot_crimes_weekly_acf = ggacf(x = tables$crime_tables$crimes_weekly$CRIMES, 
							   n = 54, 
							   main = "ACF Plot of Crimes Per Week")

plot_crimes_weekly_acf

acf_plots = list("plot_crimes_daily_acf" = plot_crimes_daily_acf,
				 "plot_crimes_weekly_acf" = plot_crimes_weekly_acf)

plots$acf_plots = acf_plots

rm(plot_crimes_daily_acf, plot_crimes_weekly_acf, acf_plots)

}

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ---- Variograms --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

{

DF = variogramDF(x = tables$crime_tables$crimes_daily$CRIMES, n = 365)

plot_crimes_daily_vario = ggplot(DF, aes(x = Lag, y = Variogram)) + 
						  geom_point(size = 2) +
						  geom_line(color = "cornflowerblue", size = 1) +
						  labs(x = "Lag", y = "Variogram") +
						  ggtitle("Variogram of Crimes Per Day") + 
						  theme_light(base_size = 20)

plot_crimes_daily_vario

DF = variogramDF(x = tables$crime_tables$crimes_weekly$CRIMES, n = 54)

plot_crimes_weekly_vario = ggplot(DF, aes(x = Lag, y = Variogram)) + 
						   geom_point(size = 2) +
						   geom_line(color = "cornflowerblue", size = 1) +
						   labs(x = "Lag", y = "Variogram") +
						   ggtitle("Variogram of Crimes Per Week") + 
						   theme_light(base_size = 20)

plot_crimes_weekly_vario

vario_plots = list("plot_crimes_daily_vario" = plot_crimes_daily_vario,
				   "plot_crimes_weekly_vario" = plot_crimes_weekly_vario)

plots$vario_plots = vario_plots

rm(plot_crimes_daily_vario, plot_crimes_weekly_vario, vario_plots)

}

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ---- Removing Seasonality & Trend --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

{

# ---- daily crimes -------------------------------------------------------------------------------------------------------------------------

# lets build a moving average variable (MA) to capture the annual trend in daily crimes

x = as.numeric(ma(tables$crime_tables$crimes_daily$CRIMES, 365))
y = rep(mean(tables$crime_tables$crimes_daily$CRIMES[1:182]), 182)
z = rep(mean(tables$crime_tables$crimes_daily$CRIMES[3701:3882]), 182)
MA.365 = c(y, x[183:3700], z)

rm(x, y, z)

# lets build the linear models data
# SIN and COS are variables that capture the annual seasonality in daily crimes
# MONTH is an external variable to increase prediction accuracy, given the monthly pattern observed in the daily crime plot

DF = data.frame("SIN" = sin((2 * pi * as.numeric(row.names(tables$crime_tables$crimes_daily))) / 365),
				"COS" = cos((2 * pi * as.numeric(row.names(tables$crime_tables$crimes_daily))) / 365),
				"MA.365" = MA.365,
				"MONTH" = tables$crime_tables$crimes_daily$MONTH,
				"CRIMES" = tables$crime_tables$crimes_daily$CRIMES)

rm(MA.365)

# create the linear model

formulation = lm(CRIMES ~ SIN + COS + MA.365 + MONTH, data = DF)

# lets look at the coefficients and model fit

stats1 = summary(formulation)
stats1

stats2 = statslm(formulation)
stats2

# lets plot the model predictions on top of daily crimes

DF2 = data.frame("DATE" = tables$crime_tables$crimes_daily$DATE, 
				 "actual" = tables$crime_tables$crimes_daily$CRIMES, 
				 "fit" = formulation$fitted.values)

plot_crimes_daily_LM = ggplot(data = DF2, aes(x = DATE, y = actual)) +
					   geom_point(na.rm = TRUE) +
					   geom_line(aes(x = DATE, y = fit), color = "cornflowerblue", na.rm = TRUE, size = 1.5) +
				   	   labs(x = "Day", y = "Crimes") +
				   	   ggtitle("Crimes Per Day\nCRIMES = -34.6 - 58.3(SIN) - 48.1(COS) + 1.1(MA.365) + 5.9(MONTHFeb) + 46.9(MONTHMar) + 61.5(MONTHApr) + 35.8(MONTHMay)\n+ 14.8(MONTHJun) - 20.2(MONTHJul) - 24.1(MONTHAug) - 46.3(MONTHSep) - 31.7(MONTHOct) - 28.4(MONTHNov) - 26.6(MONTHDec) + error") +
				   	   theme_light(base_size = 15)

plot_crimes_daily_LM

# lets take a look at the residuals

plot_crimes_daily_resid = residplots(actual = DF2$actual, fit = DF2$fit, n = 365, basefont = 15)

grid.arrange(plot_crimes_daily_resid[[1]], 
			 plot_crimes_daily_resid[[2]], 
			 plot_crimes_daily_resid[[3]], 
			 plot_crimes_daily_resid[[4]], 
			 plot_crimes_daily_resid[[5]], 
			 plot_crimes_daily_resid[[6]], 
			 ncol = 3)

# lets run tests for the 6 linear model assumptions

assumptions = testslm(formulation, 365)
assumptions

# lets store these results in lists

models = list()
LM = list()
daily = list()
tables2 = list()
fit = list()
tests = list()
plots2 = list()

tables2$model_data = data.table(DF)
tables2$model_predictions = data.table(DF2)

fit$formulation = formulation
fit$stats = list("stats1" = stats1, "stats2" = stats2)

tests$assumptions = assumptions

plots2$predictions = plot_crimes_daily_LM
plots2$residual = plot_crimes_daily_resid

daily$tables = tables2
daily$fit = fit
daily$tests = tests
daily$plots = plots2

LM$daily = daily

models$LM = LM

rm(LM, daily, tables2, fit, tests, plots2, DF, DF2, formulation, stats1, stats2, assumptions, plot_crimes_daily_LM, plot_crimes_daily_resid)

# ---- weekly crimes -------------------------------------------------------------------------------------------------------------------------

# lets build a moving average variable (MA) to capture the annual trend in weekly crimes

x = as.numeric(ma(tables$crime_tables$crimes_weekly$CRIMES, 54))
y = rep(mean(tables$crime_tables$crimes_weekly$CRIMES[1:27]), 27)
z = rep(mean(tables$crime_tables$crimes_weekly$CRIMES[539:565]), 27)
MA.54 = c(y, x[28:538], z)

rm(x, y, z)

# lets build the linear models data
# SIN and COS are variables that capture the annual seasonality in weekly crimes
# DAY is an external variable to increase prediction accuracy, given the weekly pattern observed in the weekly crimes plot

DF = data.frame("SIN" = sin((2 * pi * as.numeric(row.names(tables$crime_tables$crimes_weekly))) / 54),
				"COS" = cos((2 * pi * as.numeric(row.names(tables$crime_tables$crimes_weekly))) / 54),
				"MA.54" = MA.54,
				"DAY" = tables$crime_tables$crimes_weekly$DAY,
				"CRIMES" = tables$crime_tables$crimes_weekly$CRIMES)

rm(MA.54)

# create the linear model

formulation = lm(CRIMES ~ SIN + COS + MA.54 + DAY, data = DF)

# lets look at the coefficients and model fit

stats1 = summary(formulation)
stats1

stats2 = statslm(formulation)
stats2

# lets plot the model predictions on top of weekly crimes

DF2 = data.frame("YEAR.WEEK" = tables$crime_tables$crimes_weekly$YEAR.WEEK, 
				 "actual" = tables$crime_tables$crimes_weekly$CRIMES, 
				 "fit" = formulation$fitted.values)

plot_crimes_weekly_LM = ggplot(data = DF2, aes(x = YEAR.WEEK, y = actual, group = 1)) +
						geom_point(na.rm = TRUE) +
						scale_x_discrete(breaks = c("2006.01", "2007.01", "2008.01", "2009.01", "2010.01", "2011.01", "2012.01", "2013.01", "2014.01", "2015.01", "2016.01")) +
						geom_line(aes(x = YEAR.WEEK, y = fit), color = "cornflowerblue", na.rm = TRUE, size = 1.5) +
						labs(x = "Week", y = "Crimes") +
						ggtitle("Crimes Per Week\nCRIMES = -3086.6 + 118.3(SIN) - 377.8(COS) + 1.01(MA.54) + 729.7(DAY2)\n+ 1084.4(DAY3) + 1642.3(DAY4) + 1664.6(DAY5) + 2550.3(DAY6) + 3116.3(DAY7) + error") +
						theme_light(base_size = 15) +
						theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_crimes_weekly_LM

# lets take a look at the residuals

plot_crimes_weekly_resid = residplots(actual = DF2$actual, fit = DF2$fit, n = 54, basefont = 15)

grid.arrange(plot_crimes_weekly_resid[[1]], 
			 plot_crimes_weekly_resid[[2]], 
			 plot_crimes_weekly_resid[[3]], 
			 plot_crimes_weekly_resid[[4]], 
			 plot_crimes_weekly_resid[[5]], 
			 plot_crimes_weekly_resid[[6]], 
			 ncol = 3)

# lets run tests for the 6 linear model assumptions

assumptions = testslm(formulation, 54)
assumptions

# lets store these results in lists

weekly = list()
tables2 = list()
fit = list()
tests = list()
plots2 = list()

tables2$model_data = data.table(DF)
tables2$model_predictions = data.table(DF2)

fit$formulation = formulation
fit$stats = list("stats1" = stats1, "stats2" = stats2)

tests$assumptions = assumptions

plots2$predictions = plot_crimes_weekly_LM
plots2$residual = plot_crimes_weekly_resid

weekly$tables = tables2
weekly$fit = fit
weekly$tests = tests
weekly$plots = plots2

models$LM$weekly = weekly

rm(weekly, tables2, fit, tests, plots2, DF, DF2, formulation, stats1, stats2, assumptions, plot_crimes_weekly_LM, plot_crimes_weekly_resid)

}

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ---- Transformations ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

{

# termplots for the daily crimes linear model

mytermplots = termplots(dat = models$LM$daily$tables$model_data, 
						model = models$LM$daily$fit$formulation)

grid.arrange(mytermplots[[1]], 
			 mytermplots[[2]], 
			 mytermplots[[3]], 
			 ncol = 3)

plot_pairs = ggpairs(models$LM$daily$tables$model_data,
					 columns = colnames(models$LM$daily$tables$model_data)[c(1:3, 5)],
					 lower = list(continuous = wrap(ggally_points, alpha = 1/3)), 
        			 upper = list(continuous = wrap(ggally_cor, size = 10, color = "black")),
        			 diag = list(continuous = wrap(ggally_barDiag, fill = "black", color = "white"))) +
			 theme_light(base_size = 20)

print(plot_pairs, bottomHeightProportion = .25)

models$LM$daily$plots$termplots = mytermplots
models$LM$daily$plots$pairsplot = plot_pairs

# termplots for the weekly crimes linear model

mytermplots = termplots(dat = models$LM$weekly$tables$model_data, 
						model = models$LM$weekly$fit$formulation)

grid.arrange(mytermplots[[1]], 
			 mytermplots[[2]], 
			 mytermplots[[3]], 
			 ncol = 3)

plot_pairs = ggpairs(models$LM$weekly$tables$model_data,
					 columns = colnames(models$LM$weekly$tables$model_data)[c(1:3, 5)],
					 lower = list(continuous = wrap(ggally_points, alpha = 1/3)), 
        			 upper = list(continuous = wrap(ggally_cor, size = 10, color = "black")),
        			 diag = list(continuous = wrap(ggally_barDiag, fill = "black", color = "white"))) +
			 theme_light(base_size = 20)

print(plot_pairs, bottomHeightProportion = .25)

models$LM$weekly$plots$termplots = mytermplots
models$LM$weekly$plots$pairsplot = plot_pairs

rm(mytermplots, plot_pairs)

}

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ---- Exponential Smoothing Decision Matrix -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

{

# Data without trend and without seasonality

DF = c(rnorm(30, 100, 1), rnorm(10, 90, 3), rnorm(10, 80, 3), rnorm(10, 90, 3), rnorm(40, 100, 1))

hw1 = predplot(HoltWinters(ts(DF[1:75], frequency = 1), gamma = FALSE),
			 n.ahead = 25,
			 actuals = DF,
			 xlab = "Observation", 
			 ylab = "Value", 
			 main = "Single Exponential Smoothing\nNo Trend, No Seasonality",
			 basefont = 15)

hw1

# Data with trend and without seasonality

x = seq(60, 113, 17)

DF = unlist(lapply(1:length(x), function(i) rnorm(25, x[i], 5)))

hw2 = predplot(HoltWinters(ts(DF[1:75], frequency = 1), gamma = FALSE),
			 n.ahead = 25,
			 actuals = DF,
			 xlab = "Observation", 
			 ylab = "Value", 
			 main = "Double Exponential Smoothing\nYes Trend, No Seasonality",
			 basefont = 15)

hw2

# Data without trend and with additive seasonality

x = c(rep(c(10 * (7:9), 10 * (9:7)), 4), 70)

DF = unlist(lapply(1:length(x), function(i) rnorm(4, x[i], 5)))

hw3 = predplot(HoltWinters(ts(DF[1:75], frequency = 25), beta = FALSE, seasonal = "additive"),
			 n.ahead = 25,
			 actuals = DF,
			 xlab = "Observation", 
			 ylab = "Value", 
			 main = "Triple Exponential Smoothing\nNo Trend, Yes Seasonality(Additive)",
			 basefont = 15)

hw3

# Data without trend and with multiplicative seasonality

x = c(10 * (7:9), 10 * (9:7),
	  10 * c(6, 8, 10), 10 * c(10, 8, 6),
	  10 * c(5, 8, 11), 10 * c(10, 8, 6),
	  10 * c(4, 8, 12), 10 * c(12, 8, 4),
	  30)

DF = unlist(lapply(1:length(x), function(i) rnorm(4, x[i], 5)))

hw4 = predplot(HoltWinters(ts(DF[1:75], frequency = 25), alpha = .05, beta = FALSE, seasonal = "multiplicative"),
			 n.ahead = 25,
			 actuals = DF,
			 xlab = "Observation", 
			 ylab = "Value", 
			 main = "Triple Exponential Smoothing\nNo Trend, Yes Seasonality(Multiplicative)",
			 basefont = 15)

hw4

# Data with trend and with additive seasonality

x = c(10 * (7:9), 10 * (9:7), 10 * (9:11), 10 * (11:9), 10 * (11:13), 10 * (13:11), 10 * (13:15), 10 * (15:13), 150)

DF = unlist(lapply(1:length(x), function(i) rnorm(4, x[i], 5)))

hw5 = predplot(HoltWinters(ts(DF[1:75], frequency = 25)),
			 n.ahead = 25,
			 actuals = DF,
			 xlab = "Observation", 
			 ylab = "Value", 
			 main = "Triple Exponential Smoothing\nYes Trend, Yes Seasonality(Additive)",
			 basefont = 15)

hw5

# Data with trend and with multiplicative seasonality
 
x = c(10 * (7:9), 10 * (9:7), 
	  10 * c(9, 11, 13), 10 * c(13, 11, 9), 
	  10 * c(13, 16, 19), 10 * c(19, 16, 13), 
	  10 * c(19, 23, 27), 10 * c(27, 23, 19), 
	  270)

DF = unlist(lapply(1:length(x), function(i) rnorm(4, x[i], 6)))

hw6 = predplot(HoltWinters(ts(DF[1:75], frequency = 25), seasonal = "multiplicative"),
			 n.ahead = 25,
			 actuals = DF,
			 xlab = "Observation", 
			 ylab = "Value", 
			 main = "Triple Exponential Smoothing\nYes Trend, Yes Seasonality(Multiplicative)",
			 basefont = 15)

hw6

SmoothingMatrix = arrangeGrob(hw1, 
							  hw3, 
							  hw4,
							  hw2,
							  hw5,
							  hw6,
							  ncol = 3)

grid.draw(SmoothingMatrix)

plots$SmoothingMatrix = SmoothingMatrix

rm(hw1, hw2, hw3, hw4, hw5, hw6, SmoothingMatrix, x, DF)

}

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ---- Exponential Forecasting: Central Philly -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

{

# ---- fit weekly crimes -------------------------------------------------------------------------------------------------------------------------------------------

# lets look at the default parameter values for modeling weekly crimes

DF = data.frame(tables$crime_Area1_tables$crimes_weekly[Area1 == "CP"])

hw = HoltWinters(ts(DF$CRIMES, frequency = 54))
hw$alpha
hw$beta
hw$gamma

# lets create a DOE to test other holtwinter parameter values for modeling weekly crimes

DOE = expand.grid(alpha = seq(0.01, 0.39, 0.01), 
				  beta = seq(0, 0.25, 0.01), 
				  gamma = seq(0.25, 0.75, 0.01))

DOE = rbind(c(hw$alpha, hw$beta, hw$gamma), DOE)

head(DOE)
NROW(DOE)

rm(hw)

# lets compute the residuals for each holt winters model in the DOE

hwobjects = lapply(1:NROW(DOE), function(i)
								 data.frame("actual" = DF$CRIMES[55:565],
											"fit" = as.numeric(HoltWinters(ts(DF$CRIMES, frequency = 54), 
																		   alpha = DOE$alpha[i], 
																		   beta = DOE$beta[i], 
																		   gamma = DOE$gamma[i])$fitted[,1])))

# values greater than 0.05 are ideal 

ttests = sapply(1:NROW(DOE), function(i) 
							 t.test(hwobjects[[i]]$actual - hwobjects[[i]]$fit)$p.value)

# values close to zero are ideal 

stats = data.frame(t(sapply(1:NROW(DOE), function(i)
										 data.frame(accuracy(f = hwobjects[[i]]$fit, 
															 x = hwobjects[[i]]$actual), 
													row.names = 1))))

# lets add these test results to the DOE

DOE$t.pval = ttests
DOE = cbind(DOE, stats)
cnames = colnames(DOE)
DOE[,5:9] = sapply(5:9, function(i) DOE[,i] = as.numeric(DOE[,i]))
colnames(DOE) = cnames

head(DOE)

rm(hwobjects, ttests, stats, cnames)

# lets filter out the top models

# lets plot the histograms of these reponses to find filter limits

par(mfrow = c(2, 3))
lapply(4:9, function(i) hist(DOE[,i], main = colnames(DOE)[i]))

head(DOE)

# apply filters

x = DOE[which(DOE$t.pval > 0.9 & 
			  DOE$MAPE <= 13.94 &
			  DOE$RMSE <= 89),]

par(mfrow = c(2, 3))
lapply(4:9, function(i) hist(x[,i], main = colnames(x)[i]))

x

# lets compare the auto fit (scenario 1) with the chosen fit from the DOE (scenario 33517)

DOE[c(1, 33517),]

hw = lapply(c(1, 33517), function(i)
					   HoltWinters(ts(DF$CRIMES, frequency = 54), 
								   alpha = DOE$alpha[i], 
								   beta = DOE$beta[i], 
								   gamma = DOE$gamma[i]))

# lets plot the model & residuals for the auto fit and chosen fit 

	# auto: scenario 1

plot_crimes_weekly_hw1 = predplot(obj = hw[[1]],
							  xlab = "Observation Number", 
							  ylab = "Crimes", 
							  main = "Crimes Per Week: Central Philly\nHolt-Winters Smoothing: alpha = 0.1389547, beta = 0.001463358, gamma = 0.5037004")

plot_crimes_weekly_hw1

plot_doe_resid1 = residplots(actual = DF$CRIMES[55:565], 
							 fit = as.numeric(hw[[1]]$fitted[,1]), 
							 histlabel.y = -5,
							 n = 54, 
							 basefont = 15)

grid.arrange(plot_doe_resid1[[1]], 
			 plot_doe_resid1[[2]], 
			 plot_doe_resid1[[3]], 
			 plot_doe_resid1[[4]], 
			 plot_doe_resid1[[5]], 
			 plot_doe_resid1[[6]], 
			 ncol = 3)		  
		  
	# chosen: scenario 33517

plot_crimes_weekly_hw2 = predplot(obj = hw[[2]],
							  xlab = "Observation Number", 
							  ylab = "Crimes", 
							  main = "Crimes Per Week: Central Philly\nHolt-Winters Smoothing: alpha = 0.15, beta = 0.01, gamma = 0.58")

plot_crimes_weekly_hw2

plot_doe_resid2 = residplots(actual = DF$CRIMES[55:565], 
							 fit = as.numeric(hw[[2]]$fitted[,1]), 
							 histlabel.y = -5,
							 n = 54, 
							 basefont = 15)

grid.arrange(plot_doe_resid2[[1]], 
			 plot_doe_resid2[[2]], 
			 plot_doe_resid2[[3]], 
			 plot_doe_resid2[[4]], 
			 plot_doe_resid2[[5]], 
			 plot_doe_resid2[[6]], 
			 ncol = 3)		  

# lets store our results into lists

HW = list()
weekly = list()
plots2 = list()
fits = list()
tables2 = list()

fits$crimes_weekly_CP = hw[[2]]
plots2$plot_crimes_weekly_CP = plot_crimes_weekly_hw2
plots2$plot_crimes_weekly_residual_CP = plot_doe_resid2
tables2$model_data_CP = DF
tables2$model_doe_CP = DOE

weekly$fits = fits
weekly$plots = plots2
weekly$tables = tables2

HW$weekly = weekly
models$HW = HW

rm(x, HW, weekly, fits, plots2, tables2, hw, plot_crimes_weekly_hw1, plot_crimes_weekly_hw2, 
   plot_doe_resid1, plot_doe_resid2, DOE, DF)
  
# ---- forecast weekly crimes -------------------------------------------------------------------------------------------------------------------------------------------

# lets look at the default parameter values for modeling weekly crimes

DF = data.frame(tables$crime_Area1_tables$crimes_weekly[Area1 == "CP"])

hw = HoltWinters(ts(DF$CRIMES[1:478], frequency = 54))
hw$alpha
hw$beta
hw$gamma

# lets create a DOE to test other holtwinter parameter values for modeling weekly crimes

DOE = expand.grid(alpha = seq(0.1, 0.35, 0.01), 
				  beta = seq(0, 0.25, 0.01), 
				  gamma = seq(0.25, 0.75, 0.01))

DOE = rbind(c(hw$alpha, hw$beta, hw$gamma), models$HW$weekly$tables$model_doe_CP[1,1:3], DOE)

head(DOE)
NROW(DOE)

rm(hw)

# lets compute the residuals of the forecast for each holt winters model in the DOE

hwobjects = lapply(1:NROW(DOE), function(i)
								 data.frame("actual" = DF$CRIMES[479:565],
											"fit" = as.numeric(predict(HoltWinters(ts(DF$CRIMES[1:478], frequency = 54), 
																				   alpha = DOE$alpha[i], 
																				   beta = DOE$beta[i], 
																				   gamma = DOE$gamma[i]), 
															   n.ahead = 87))))

# values greater than 0.05 are ideal 

ttests = sapply(1:NROW(DOE), function(i) 
							 t.test(hwobjects[[i]]$actual - hwobjects[[i]]$fit)$p.value)

# values close to zero are ideal 

stats = data.frame(t(sapply(1:NROW(DOE), function(i)
										 data.frame(accuracy(f = hwobjects[[i]]$fit, 
															 x = hwobjects[[i]]$actual), 
													row.names = 1))))

# lets add these test results to the DOE

DOE$t.pval = ttests
DOE = cbind(DOE, stats)

cnames = colnames(DOE)
DOE[,5:9] = sapply(5:9, function(i) DOE[,i] = as.numeric(DOE[,i]))
colnames(DOE) = cnames

head(DOE)

rm(hwobjects, ttests, stats, cnames)

# lets filter out the top models

# lets plot the histograms of these reponses to find filter limits

par(mfrow = c(2, 3))
lapply(4:9, function(i) hist(DOE[,i], main = colnames(DOE)[i]))

head(DOE)

# apply filters

x = DOE[which(DOE$t.pval > 0.4 &
			  DOE$MAPE <= 12 &
			  DOE$RMSE <= 76),]

par(mfrow = c(2, 3))
lapply(4:9, function(i) hist(x[,i], main = colnames(x)[i]))

x

# lets compare the auto fits (scenario 1 & scenario 2) with the chosen fit from the DOE (scenario 10352)

DOE[c(1, 2, 10352),]

hw = lapply(c(1, 2, 10352), function(i)
					   HoltWinters(ts(DF$CRIMES[1:478], frequency = 54), 
								   alpha = DOE$alpha[i], 
								   beta = DOE$beta[i], 
								   gamma = DOE$gamma[i]))

# lets plot the model & residuals for each of the final models

	# scenario 1
  
plot_crimes_weekly_hw1 = predplot(obj = hw[[1]],
							  n.ahead = 87,
							  actuals = DF$CRIMES,
							  xlab = "Observation Number", 
							  ylab = "Crimes", 
							  main = "Crimes Per Week: Central Philly\nHolt-Winters Smoothing: alpha = 0.09852259, beta = 0.003108433, gamma = 0.5078421")

plot_crimes_weekly_hw1
	
plot_doe_resid1 = residplots(actual = DF$CRIMES[479:565], 
							 fit = as.numeric(predict(hw[[1]], 
													  n.ahead = 87)), 
							 histlabel.y = -0.5,
							 n = 54,
							 binwidth = 25,							 
							 basefont = 15)

grid.arrange(plot_doe_resid1[[1]], 
			 plot_doe_resid1[[2]], 
			 plot_doe_resid1[[3]], 
			 plot_doe_resid1[[4]], 
			 plot_doe_resid1[[5]], 
			 plot_doe_resid1[[6]], 
			 ncol = 3)		  		  

	# scenario 2

plot_crimes_weekly_hw2 = predplot(obj = hw[[2]],
							  n.ahead = 87,
							  actuals = DF$CRIMES,
							  xlab = "Observation Number", 
							  ylab = "Crimes", 
							  main = "Crimes Per Week: Central Philly\nHolt-Winters Smoothing: alpha = 0.1389547, beta = 0.001463358, gamma = 0.5037004")

plot_crimes_weekly_hw2
	
plot_doe_resid2 = residplots(actual = DF$CRIMES[479:565], 
							 fit = as.numeric(predict(hw[[2]], 
													  n.ahead = 87)), 
							 histlabel.y = -0.5,
							 n = 54,
							 binwidth = 25,
							 basefont = 15)

grid.arrange(plot_doe_resid2[[1]], 
			 plot_doe_resid2[[2]], 
			 plot_doe_resid2[[3]], 
			 plot_doe_resid2[[4]], 
			 plot_doe_resid2[[5]], 
			 plot_doe_resid2[[6]], 
			 ncol = 3)		  
		  
	# scenario 10352

plot_crimes_weekly_hw3 = predplot(obj = hw[[3]],
							  n.ahead = 87,
							  actuals = DF$CRIMES,
							  xlab = "Observation Number", 
							  ylab = "Crimes", 
							  main = "Crimes Per Week: Central Philly\nHolt-Winters Smoothing: alpha = 0.11, beta = 0.08, gamma = 0.40")

plot_crimes_weekly_hw3
	
plot_doe_resid3 = residplots(actual = DF$CRIMES[479:565], 
							 fit = as.numeric(predict(hw[[3]], 
													  n.ahead = 87)), 
							 histlabel.y = -0.5,
							 n = 54, 
							 binwidth = 25,
							 basefont = 15)

grid.arrange(plot_doe_resid3[[1]], 
			 plot_doe_resid3[[2]], 
			 plot_doe_resid3[[3]], 
			 plot_doe_resid3[[4]], 
			 plot_doe_resid3[[5]], 
			 plot_doe_resid3[[6]], 
			 ncol = 3)		  

# lets store our results into lists

models$HW$weekly$fits$crimes_weekly_CP_forecast = hw[[3]]
models$HW$weekly$plots$plot_crimes_weekly_CP_forecast = plot_crimes_weekly_hw3
models$HW$weekly$plots$plot_crimes_weekly_residual_CP_forecast = plot_doe_resid3
models$HW$weekly$tables$model_data_CP_forecast = DF
models$HW$weekly$tables$model_doe_CP_forecast = DOE		  

rm(x, hw, plot_crimes_weekly_hw1, plot_crimes_weekly_hw2, plot_crimes_weekly_hw3, 
   plot_doe_resid1, plot_doe_resid2, plot_doe_resid3, 
   DOE, DF)		  

}		  

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ---- Exponential Forecasting: West Philly -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

{

# ---- fit weekly crimes -------------------------------------------------------------------------------------------------------------------------------------------

# lets look at the default parameter values for modeling weekly crimes

DF = data.frame(tables$crime_Area1_tables$crimes_weekly[Area1 == "WP"])

hw = HoltWinters(ts(DF$CRIMES, frequency = 54))
hw$alpha
hw$beta
hw$gamma

# lets create a DOE to test other holtwinter parameter values for modeling weekly crimes

DOE = expand.grid(alpha = seq(0.01, 0.37, 0.01), 
				  beta = seq(0, 0.25, 0.01), 
				  gamma = seq(0.28, 0.78, 0.01))

DOE = rbind(c(hw$alpha, hw$beta, hw$gamma), DOE)

head(DOE)
NROW(DOE)

rm(hw)

# lets compute the residuals for each holt winters model in the DOE

hwobjects = lapply(1:NROW(DOE), function(i)
								 data.frame("actual" = DF$CRIMES[55:565],
											"fit" = as.numeric(HoltWinters(ts(DF$CRIMES, frequency = 54), 
																		   alpha = DOE$alpha[i], 
																		   beta = DOE$beta[i], 
																		   gamma = DOE$gamma[i])$fitted[,1])))

# values greater than 0.05 are ideal 

ttests = sapply(1:NROW(DOE), function(i) 
							 t.test(hwobjects[[i]]$actual - hwobjects[[i]]$fit)$p.value)

# values close to zero are ideal 

stats = data.frame(t(sapply(1:NROW(DOE), function(i)
										 data.frame(accuracy(f = hwobjects[[i]]$fit, 
															 x = hwobjects[[i]]$actual), 
													row.names = 1))))

# lets add these test results to the DOE

DOE$t.pval = ttests
DOE = cbind(DOE, stats)

cnames = colnames(DOE)
DOE[,5:9] = sapply(5:9, function(i) DOE[,i] = as.numeric(DOE[,i]))
colnames(DOE) = cnames

head(DOE)

rm(hwobjects, ttests, stats, cnames)

# lets filter out the top models

# lets plot the histograms of these reponses to find filter limits

par(mfrow = c(2, 3))
lapply(4:9, function(i) hist(DOE[,i], main = colnames(DOE)[i]))

head(DOE)

# apply filters

x = DOE[which(DOE$t.pval > 0.8 & 
			  DOE$MAPE <= 14.4 &
			  DOE$RMSE <= 90),]

par(mfrow = c(2, 3))
lapply(4:9, function(i) hist(x[,i], main = colnames(x)[i]))

x

# lets compare the auto fit (scenario 1) with the chosen fit from the DOE (scenario 22157)

DOE[c(1, 22157),]

hw = lapply(c(1, 22157), function(i)
					   HoltWinters(ts(DF$CRIMES, frequency = 54), 
								   alpha = DOE$alpha[i], 
								   beta = DOE$beta[i], 
								   gamma = DOE$gamma[i]))

# lets plot the model & residuals for the auto fit and chosen fit 

	# auto: scenario 1

plot_crimes_weekly_hw1 = predplot(obj = hw[[1]],
							  xlab = "Observation Number", 
							  ylab = "Crimes", 
							  main = "Crimes Per Week: West Philly\nHolt-Winters Smoothing: alpha = 0.1210642, beta = 0, gamma = 0.5259288")

plot_crimes_weekly_hw1

plot_doe_resid1 = residplots(actual = DF$CRIMES[55:565], 
							 fit = as.numeric(hw[[1]]$fitted[,1]), 
							 histlabel.y = -5,
							 n = 54, 
							 basefont = 15)

grid.arrange(plot_doe_resid1[[1]], 
			 plot_doe_resid1[[2]], 
			 plot_doe_resid1[[3]], 
			 plot_doe_resid1[[4]], 
			 plot_doe_resid1[[5]], 
			 plot_doe_resid1[[6]], 
			 ncol = 3)		  
		  
	# chosen: scenario 22157

plot_crimes_weekly_hw2 = predplot(obj = hw[[2]],
							  xlab = "Observation Number", 
							  ylab = "Crimes", 
							  main = "Crimes Per Week: West Philly\nHolt-Winters Smoothing: alpha = 0.30, beta = 0, gamma = 0.51")

plot_crimes_weekly_hw2

plot_doe_resid2 = residplots(actual = DF$CRIMES[55:565], 
							 fit = as.numeric(hw[[2]]$fitted[,1]), 
							 histlabel.y = -5,
							 n = 54, 
							 basefont = 15)

grid.arrange(plot_doe_resid2[[1]], 
			 plot_doe_resid2[[2]], 
			 plot_doe_resid2[[3]], 
			 plot_doe_resid2[[4]], 
			 plot_doe_resid2[[5]], 
			 plot_doe_resid2[[6]], 
			 ncol = 3)		  

# lets store our results into lists

models$HW$weekly$fits$crimes_weekly_WP = hw[[2]]
models$HW$weekly$plots$plot_crimes_weekly_WP = plot_crimes_weekly_hw2
models$HW$weekly$plots$plot_crimes_weekly_residual_WP = plot_doe_resid2
models$HW$weekly$tables$model_data_WP = DF
models$HW$weekly$tables$model_doe_WP = DOE		  

rm(x, hw, plot_crimes_weekly_hw1, plot_crimes_weekly_hw2, 
   plot_doe_resid1, plot_doe_resid2, DOE, DF)

# ---- forecast weekly crimes -------------------------------------------------------------------------------------------------------------------------------------------

# lets look at the default parameter values for modeling weekly crimes

DF = data.frame(tables$crime_Area1_tables$crimes_weekly[Area1 == "WP"])

hw = HoltWinters(ts(DF$CRIMES[1:478], frequency = 54))
hw$alpha
hw$beta
hw$gamma

# lets create a DOE to test other holtwinter parameter values for modeling weekly crimes

DOE = expand.grid(alpha = seq(0.01, 0.51, 0.01), 
				  beta = seq(0, 0.25, 0.01), 
				  gamma = seq(0.24, 0.74, 0.01))

DOE = rbind(c(hw$alpha, hw$beta, hw$gamma), models$HW$weekly$tables$model_doe_WP[1,1:3], DOE)

head(DOE)
NROW(DOE)

rm(hw)

# lets compute the residuals of the forecast for each holt winters model in the DOE

hwobjects = lapply(1:NROW(DOE), function(i)
								 data.frame("actual" = DF$CRIMES[479:565],
											"fit" = as.numeric(predict(HoltWinters(ts(DF$CRIMES[1:478], frequency = 54), 
																				   alpha = DOE$alpha[i], 
																				   beta = DOE$beta[i], 
																				   gamma = DOE$gamma[i]), 
															   n.ahead = 87))))

# values greater than 0.05 are ideal 

ttests = sapply(1:NROW(DOE), function(i) 
							 t.test(hwobjects[[i]]$actual - hwobjects[[i]]$fit)$p.value)

# values close to zero are ideal 

stats = data.frame(t(sapply(1:NROW(DOE), function(i)
										 data.frame(accuracy(f = hwobjects[[i]]$fit, 
															 x = hwobjects[[i]]$actual), 
													row.names = 1))))

# lets add these test results to the DOE

DOE$t.pval = ttests
DOE = cbind(DOE, stats)

cnames = colnames(DOE)
DOE[,5:9] = sapply(5:9, function(i) DOE[,i] = as.numeric(DOE[,i]))
colnames(DOE) = cnames

head(DOE)

rm(hwobjects, ttests, stats, cnames)

# lets filter out the top models

# lets plot the histograms of these reponses to find filter limits

par(mfrow = c(2, 3))
lapply(4:9, function(i) hist(DOE[,i], main = colnames(DOE)[i]))

head(DOE)

# apply filters

x = DOE[which(DOE$t.pval > 0.9 &
			  DOE$MAPE <= 12 &
			  DOE$RMSE <= 81 &
			  DOE$MAE <= 59),]

par(mfrow = c(2, 3))
lapply(4:9, function(i) hist(x[,i], main = colnames(x)[i]))

x

# lets compare the auto fits (scenario 1 & scenario 2) with the chosen fit from the DOE (scenario 8827)

DOE[c(1, 2, 8827),]

hw = lapply(c(1, 2, 8827), function(i)
					   HoltWinters(ts(DF$CRIMES[1:478], frequency = 54), 
								   alpha = DOE$alpha[i], 
								   beta = DOE$beta[i], 
								   gamma = DOE$gamma[i]))

# lets plot the model & residuals for each of the final models

	# scenario 1

plot_crimes_weekly_hw1 = predplot(obj = hw[[1]],
							  n.ahead = 87,
							  actuals = DF$CRIMES,
							  xlab = "Observation Number", 
							  ylab = "Crimes", 
							  main = "Crimes Per Week: West Philly\nHolt-Winters Smoothing: alpha = 0.2641452, beta = 0, gamma = 0.4865719")

plot_crimes_weekly_hw1

plot_doe_resid1 = residplots(actual = DF$CRIMES[479:565], 
							 fit = as.numeric(predict(hw[[1]], 
													  n.ahead = 87)), 
							 histlabel.y = -0.5,
							 n = 54, 
							 binwidth = 40,
							 basefont = 15)

grid.arrange(plot_doe_resid1[[1]], 
			 plot_doe_resid1[[2]], 
			 plot_doe_resid1[[3]], 
			 plot_doe_resid1[[4]], 
			 plot_doe_resid1[[5]], 
			 plot_doe_resid1[[6]], 
			 ncol = 3)		  		  

	# scenario 2

plot_crimes_weekly_hw2 = predplot(obj = hw[[2]],
							  n.ahead = 87,
							  actuals = DF$CRIMES,
							  xlab = "Observation Number", 
							  ylab = "Crimes", 
							  main = "Crimes Per Week: West Philly\nHolt-Winters Smoothing: alpha = 0.1210642, beta = 0, gamma = 0.5259288")

plot_crimes_weekly_hw2

plot_doe_resid2 = residplots(actual = DF$CRIMES[479:565], 
							 fit = as.numeric(predict(hw[[2]], 
													  n.ahead = 87)), 
							 histlabel.y = -0.5,
							 n = 54, 
							 binwidth = 40,
							 basefont = 15)

grid.arrange(plot_doe_resid2[[1]], 
			 plot_doe_resid2[[2]], 
			 plot_doe_resid2[[3]], 
			 plot_doe_resid2[[4]], 
			 plot_doe_resid2[[5]], 
			 plot_doe_resid2[[6]], 
			 ncol = 3)		  
		  
	# scenario 8827

plot_crimes_weekly_hw3 = predplot(obj = hw[[3]],
							  n.ahead = 87,
							  actuals = DF$CRIMES,
							  xlab = "Observation Number", 
							  ylab = "Crimes", 
							  main = "Crimes Per Week: West Philly\nHolt-Winters Smoothing: alpha = 0.02, beta = 0.17, gamma = 0.30")

plot_crimes_weekly_hw3

plot_doe_resid3 = residplots(actual = DF$CRIMES[479:565], 
							 fit = as.numeric(predict(hw[[3]], 
													  n.ahead = 87)), 
							 histlabel.y = -0.5,
							 n = 54, 
							 binwidth = 40,
							 basefont = 15)

grid.arrange(plot_doe_resid3[[1]], 
			 plot_doe_resid3[[2]], 
			 plot_doe_resid3[[3]], 
			 plot_doe_resid3[[4]], 
			 plot_doe_resid3[[5]], 
			 plot_doe_resid3[[6]], 
			 ncol = 3)		  

# lets store our results into lists

models$HW$weekly$fits$crimes_weekly_WP_forecast = hw[[3]]
models$HW$weekly$plots$plot_crimes_weekly_WP_forecast = plot_crimes_weekly_hw3
models$HW$weekly$plots$plot_crimes_weekly_residual_WP_forecast = plot_doe_resid3
models$HW$weekly$tables$model_data_WP_forecast = DF
models$HW$weekly$tables$model_doe_WP_forecast = DOE		  

rm(x, hw, plot_crimes_weekly_hw1, plot_crimes_weekly_hw2, plot_crimes_weekly_hw3, 
   plot_doe_resid1, plot_doe_resid2, plot_doe_resid3, 
   DOE, DF)		  

}		  

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ---- Exponential Forecasting: North Philly -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

{

# ---- fit weekly crimes -------------------------------------------------------------------------------------------------------------------------------------------

# lets look at the default parameter values for modeling weekly crimes

DF = data.frame(tables$crime_Area1_tables$crimes_weekly[Area1 == "NP"])

hw = HoltWinters(ts(DF$CRIMES, frequency = 54))
hw$alpha
hw$beta
hw$gamma

# lets create a DOE to test other holtwinter parameter values for modeling weekly crimes

DOE = expand.grid(alpha = seq(0.01, 0.38, 0.01), 
				  beta = seq(0, 0.25, 0.01), 
				  gamma = seq(0.34, 0.84, 0.01))

DOE = rbind(c(hw$alpha, hw$beta, hw$gamma), DOE)

head(DOE)
NROW(DOE)

rm(hw)

# lets compute the residuals for each holt winters model in the DOE

hwobjects = lapply(1:NROW(DOE), function(i)
								 data.frame("actual" = DF$CRIMES[55:565],
											"fit" = as.numeric(HoltWinters(ts(DF$CRIMES, frequency = 54), 
																		   alpha = DOE$alpha[i], 
																		   beta = DOE$beta[i], 
																		   gamma = DOE$gamma[i])$fitted[,1])))

# values greater than 0.05 are ideal 

ttests = sapply(1:NROW(DOE), function(i) 
							 t.test(hwobjects[[i]]$actual - hwobjects[[i]]$fit)$p.value)

# values close to zero are ideal 

stats = data.frame(t(sapply(1:NROW(DOE), function(i)
										 data.frame(accuracy(f = hwobjects[[i]]$fit, 
															 x = hwobjects[[i]]$actual), 
													row.names = 1))))

# lets add these test results to the DOE

DOE$t.pval = ttests
DOE = cbind(DOE, stats)

cnames = colnames(DOE)
DOE[,5:9] = sapply(5:9, function(i) DOE[,i] = as.numeric(DOE[,i]))
colnames(DOE) = cnames

head(DOE)

rm(hwobjects, ttests, stats, cnames)

# lets filter out the top models

# lets plot the histograms of these reponses to find filter limits

par(mfrow = c(2, 3))
lapply(4:9, function(i) hist(DOE[,i], main = colnames(DOE)[i]))

head(DOE)

# apply filters

x = DOE[which(DOE$t.pval > 0.7 & 
			  DOE$MAPE <= 12.58 &
			  DOE$RMSE <= 218.5 &
			  DOE$MAE <= 137.4),]

par(mfrow = c(2, 3))
lapply(4:9, function(i) hist(x[,i], main = colnames(x)[i]))

x

# lets compare the auto fit (scenario 1) with the chosen fit from the DOE (scenario 25701)

DOE[c(1, 25701),]

hw = lapply(c(1, 25701), function(i)
					   HoltWinters(ts(DF$CRIMES, frequency = 54), 
								   alpha = DOE$alpha[i], 
								   beta = DOE$beta[i], 
								   gamma = DOE$gamma[i]))

# lets plots the residuals for the auto fit and chosen fit 

	# auto: scenario 1

plot_crimes_weekly_hw1 = predplot(obj = hw[[1]],
							  xlab = "Observation Number", 
							  ylab = "Crimes", 
							  main = "Crimes Per Week: North Philly\nHolt-Winters Smoothing: alpha = 0.1308771, beta = 0, gamma = 0.5850141")

plot_crimes_weekly_hw1

plot_doe_resid1 = residplots(actual = DF$CRIMES[55:565], 
							 fit = as.numeric(hw[[1]]$fitted[,1]), 
							 histlabel.y = -5,
							 n = 54, 
							 basefont = 15)

grid.arrange(plot_doe_resid1[[1]], 
			 plot_doe_resid1[[2]], 
			 plot_doe_resid1[[3]], 
			 plot_doe_resid1[[4]], 
			 plot_doe_resid1[[5]], 
			 plot_doe_resid1[[6]], 
			 ncol = 3)		  

	# chosen: scenario 6808

plot_crimes_weekly_hw2 = predplot(obj = hw[[2]],
							  xlab = "Observation Number", 
							  ylab = "Crimes", 
							  main = "Crimes Per Week: North Philly\nHolt-Winters Smoothing: alpha = 0.12, beta = 0, gamma = 0.60")

plot_crimes_weekly_hw2

plot_doe_resid2 = residplots(actual = DF$CRIMES[55:565], 
							 fit = as.numeric(hw[[2]]$fitted[,1]), 
							 histlabel.y = -5,
							 n = 54, 
							 basefont = 15)

grid.arrange(plot_doe_resid2[[1]], 
			 plot_doe_resid2[[2]], 
			 plot_doe_resid2[[3]], 
			 plot_doe_resid2[[4]], 
			 plot_doe_resid2[[5]], 
			 plot_doe_resid2[[6]], 
			 ncol = 3)		  

# lets store our results into lists

models$HW$weekly$fits$crimes_weekly_NP = hw[[2]]
models$HW$weekly$plots$plot_crimes_weekly_NP = plot_crimes_weekly_hw2
models$HW$weekly$plots$plot_crimes_weekly_residual_NP = plot_doe_resid2
models$HW$weekly$tables$model_data_NP = DF
models$HW$weekly$tables$model_doe_NP = DOE		  

rm(x, hw, plot_crimes_weekly_hw1, plot_crimes_weekly_hw2, 
   plot_doe_resid1, plot_doe_resid2, DOE, DF)

# ---- forecast weekly crimes -------------------------------------------------------------------------------------------------------------------------------------------

# lets look at the default parameter values for modeling weekly crimes

DF = data.frame(tables$crime_Area1_tables$crimes_weekly[Area1 == "NP"])

hw = HoltWinters(ts(DF$CRIMES[1:478], frequency = 54))
hw$alpha
hw$beta
hw$gamma

# lets create a DOE to test other holtwinter parameter values for modeling weekly crimes

DOE = expand.grid(alpha = seq(0.01, 0.34, 0.01), 
				  beta = seq(0, 0.25, 0.01), 
				  gamma = seq(0.34, 0.84, 0.01))

DOE = rbind(c(hw$alpha, hw$beta, hw$gamma), models$HW$weekly$tables$model_doe_NP[1,1:3], DOE)

head(DOE)
NROW(DOE)

rm(hw)

# lets compute the residuals of the forecast for each holt winters model in the DOE

hwobjects = lapply(1:NROW(DOE), function(i)
								 data.frame("actual" = DF$CRIMES[479:565],
											"fit" = as.numeric(predict(HoltWinters(ts(DF$CRIMES[1:478], frequency = 54), 
																				   alpha = DOE$alpha[i], 
																				   beta = DOE$beta[i], 
																				   gamma = DOE$gamma[i]), 
															   n.ahead = 87))))

# values greater than 0.05 are ideal 

ttests = sapply(1:NROW(DOE), function(i) 
							 t.test(hwobjects[[i]]$actual - hwobjects[[i]]$fit)$p.value)

# values close to zero are ideal 

stats = data.frame(t(sapply(1:NROW(DOE), function(i)
										 data.frame(accuracy(f = hwobjects[[i]]$fit, 
															 x = hwobjects[[i]]$actual), 
													row.names = 1))))

# lets add these test results to the DOE

DOE$t.pval = ttests
DOE = cbind(DOE, stats)

cnames = colnames(DOE)
DOE[,5:9] = sapply(5:9, function(i) DOE[,i] = as.numeric(DOE[,i]))
colnames(DOE) = cnames

head(DOE)

rm(hwobjects, ttests, stats, cnames)

# lets filter out the top models

# lets plot the histograms of these reponses to find filter limits

par(mfrow = c(2, 3))
lapply(4:9, function(i) hist(DOE[,i], main = colnames(DOE)[i]))

head(DOE)

# apply filters

x = DOE[which(DOE$t.pval > 0.6 &
			  DOE$MAPE <= 10 &
			  DOE$RMSE <= 187),]

par(mfrow = c(2, 3))
lapply(4:9, function(i) hist(x[,i], main = colnames(x)[i]))

x

# lets compare the auto fits (scenario 1 & scenario 2) with the chosen fit from the DOE (scenario 7257)

DOE[c(1, 2, 7257),]

hw = lapply(c(1, 2, 7257), function(i)
					   HoltWinters(ts(DF$CRIMES[1:478], frequency = 54), 
								   alpha = DOE$alpha[i], 
								   beta = DOE$beta[i], 
								   gamma = DOE$gamma[i]))

# lets plots the residuals for each of the final models

	# scenario 1

plot_crimes_weekly_hw1 = predplot(obj = hw[[1]],
							  n.ahead = 87,
							  actuals = DF$CRIMES,
							  xlab = "Observation Number", 
							  ylab = "Crimes", 
							  main = "Crimes Per Week: North Philly\nHolt-Winters Smoothing: alpha = 0.0894020, beta = 0, gamma = 0.5868770")

plot_crimes_weekly_hw1

plot_doe_resid1 = residplots(actual = DF$CRIMES[479:565], 
							 fit = as.numeric(predict(hw[[1]], 
													  n.ahead = 87)),  
							 n = 54, 
							 histlabel.y = -0.5,
							 basefont = 15)

grid.arrange(plot_doe_resid1[[1]], 
			 plot_doe_resid1[[2]], 
			 plot_doe_resid1[[3]], 
			 plot_doe_resid1[[4]], 
			 plot_doe_resid1[[5]], 
			 plot_doe_resid1[[6]], 
			 ncol = 3)		  		  

	# scenario 2

plot_crimes_weekly_hw2 = predplot(obj = hw[[2]],
							  n.ahead = 87,
							  actuals = DF$CRIMES,
							  xlab = "Observation Number", 
							  ylab = "Crimes", 
							  main = "Crimes Per Week: North Philly\nHolt-Winters Smoothing: alpha = 0.1308771, beta = 0, gamma = 0.5850141")

plot_crimes_weekly_hw2

plot_doe_resid2 = residplots(actual = DF$CRIMES[479:565], 
							 fit = as.numeric(predict(hw[[2]], 
													  n.ahead = 87)), 
							 n = 54, 
							 histlabel.y = -0.5,
							 basefont = 15)

grid.arrange(plot_doe_resid2[[1]], 
			 plot_doe_resid2[[2]], 
			 plot_doe_resid2[[3]], 
			 plot_doe_resid2[[4]], 
			 plot_doe_resid2[[5]], 
			 plot_doe_resid2[[6]], 
			 ncol = 3)		  
		  
	# scenario 7257

plot_crimes_weekly_hw3 = predplot(obj = hw[[3]],
							  n.ahead = 87,
							  actuals = DF$CRIMES,
							  xlab = "Observation Number", 
							  ylab = "Crimes", 
							  main = "Crimes Per Week: North Philly\nHolt-Winters Smoothing: alpha = 0.13, beta = 0.05, gamma = 0.42")

plot_crimes_weekly_hw3

plot_doe_resid3 = residplots(actual = DF$CRIMES[479:565], 
							 fit = as.numeric(predict(hw[[3]], 
													  n.ahead = 87)), 
							 n = 54, 
							 histlabel.y = -0.5,
							 basefont = 15)

grid.arrange(plot_doe_resid3[[1]], 
			 plot_doe_resid3[[2]], 
			 plot_doe_resid3[[3]], 
			 plot_doe_resid3[[4]], 
			 plot_doe_resid3[[5]], 
			 plot_doe_resid3[[6]], 
			 ncol = 3)		  

# lets store our results into lists

models$HW$weekly$fits$crimes_weekly_NP_forecast = hw[[3]]
models$HW$weekly$plots$plot_crimes_weekly_NP_forecast = plot_crimes_weekly_hw3
models$HW$weekly$plots$plot_crimes_weekly_residual_NP_forecast = plot_doe_resid3
models$HW$weekly$tables$model_data_NP_forecast = DF
models$HW$weekly$tables$model_doe_NP_forecast = DOE		  

rm(x, hw, plot_crimes_weekly_hw1, plot_crimes_weekly_hw2, plot_crimes_weekly_hw3, 
   plot_doe_resid1, plot_doe_resid2, plot_doe_resid3, 
   DOE, DF)		  
  
}	

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ---- ACF & PACF Plots ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

{

# ---- weekly ACF & PACF plots -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	# central philly crimes

		# ACF

models$HW$weekly$plots$plot_crimes_weekly_residual_CP$acfPlot

		# PACF

models$HW$weekly$plots$plot_crimes_weekly_residual_CP$pacfPlot = ggacf(x = as.numeric(resid(models$HW$weekly$fits$crimes_weekly_CP)), partial = TRUE, n = 54, basefont = 15, main = "PACF Plot of Residuals")
models$HW$weekly$plots$plot_crimes_weekly_residual_CP$pacfPlot

grid.arrange(models$HW$weekly$plots$plot_crimes_weekly_residual_CP$acfPlot,
			 models$HW$weekly$plots$plot_crimes_weekly_residual_CP$pacfPlot,
			 ncol = 2)

	# west philly crimes

		# ACF

models$HW$weekly$plots$plot_crimes_weekly_residual_WP$acfPlot

		# PACF

models$HW$weekly$plots$plot_crimes_weekly_residual_WP$pacfPlot = ggacf(x = as.numeric(resid(models$HW$weekly$fits$crimes_weekly_WP)), partial = TRUE, n = 54, basefont = 15, main = "PACF Plot of Residuals")
models$HW$weekly$plots$plot_crimes_weekly_residual_WP$pacfPlot

grid.arrange(models$HW$weekly$plots$plot_crimes_weekly_residual_WP$acfPlot,
			 models$HW$weekly$plots$plot_crimes_weekly_residual_WP$pacfPlot,
			 ncol = 2)

	# north philly crimes

		# ACF

models$HW$weekly$plots$plot_crimes_weekly_residual_NP$acfPlot

		# PACF

models$HW$weekly$plots$plot_crimes_weekly_residual_NP$pacfPlot = ggacf(x = as.numeric(resid(models$HW$weekly$fits$crimes_weekly_NP)), partial = TRUE, n = 54, basefont = 15, main = "PACF Plot of Residuals")
models$HW$weekly$plots$plot_crimes_weekly_residual_NP$pacfPlot

grid.arrange(models$HW$weekly$plots$plot_crimes_weekly_residual_NP$acfPlot,
			 models$HW$weekly$plots$plot_crimes_weekly_residual_NP$pacfPlot,
			 ncol = 2)

}

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ---- ARIMAX Forecasting: Central Philly -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

{

# ---- parameter issues with Arima(): whenever P != 0 --------------------------------------------------------------------------------------------------

DF = data.frame(tables$crime_Area1_tables$crimes_weekly[Area1 == "CP"])

# P != 0, Q != 0

DOE = expand.grid("p" = 1:3, 
				  "d" = 1, 
				  "q" = 1:2,
				  "P" = 1,
				  "D" = 1,
				  "Q" = 1)

DOE

# i = c(1, 4, 6) doesn't work
# i = c(2, 3, 5) does work

lapply(1:NROW(DOE), function(i)
Arima.if(ts(DF$CRIMES, frequency = 54), 
		order = c(DOE$p[i], DOE$d[i], DOE$q[i]), 
		seasonal = list(order = c(DOE$P[i], DOE$D[i], DOE$Q[i]), period = 54),
		optim.control = list(maxit = 1000)))

# P = 0, Q != 0

DOE = expand.grid("p" = 1:3, 
				  "d" = 1, 
				  "q" = 1:2,
				  "P" = 0,
				  "D" = 1,
				  "Q" = 1)

DOE

# i = 1:6 does work

lapply(1:NROW(DOE), function(i)
Arima.if(ts(DF$CRIMES, frequency = 54), 
		order = c(DOE$p[i], DOE$d[i], DOE$q[i]), 
		seasonal = list(order = c(DOE$P[i], DOE$D[i], DOE$Q[i]), period = 54),
		optim.control = list(maxit = 1000)))

# P != 0, Q = 0

DOE = expand.grid("p" = 1:3, 
				  "d" = 1, 
				  "q" = 1:2,
				  "P" = 1,
				  "D" = 1,
				  "Q" = 0)

DOE

# i = c(1, 4, 6) doesn't work
# i = c(2, 3, 5) does work

lapply(1:NROW(DOE), function(i)
Arima.if(ts(DF$CRIMES, frequency = 54), 
		order = c(DOE$p[i], DOE$d[i], DOE$q[i]), 
		seasonal = list(order = c(DOE$P[i], DOE$D[i], DOE$Q[i]), period = 54),
		optim.control = list(maxit = 1000)))

# ---- fit weekly crimes -------------------------------------------------------------------------------------------------------------------------------------------

DF = data.frame(tables$crime_Area1_tables$crimes_weekly[Area1 == "CP"])

# lets look at differenced weekly crimes to find values for d

par(mfrow = c(2,3))

plot(DF$CRIMES, main = "d = 0", type = "b")
lapply(1:5, function(i) plot(diff(DF$CRIMES, differences = i), main = paste0("d = ", i), type ="b"))

d = 1

# lets look at the acf and pacf plots of the differenced weekly crimes to find values of p and q

d = 1
grid.arrange(ggacf(x = diff(DF$CRIMES, differences = d), partial = TRUE, n = 54, basefont = 15, main = paste0("PACF Plot of Central Philly Weekly Crimes\nDifference of Order ", d)),
			 ggacf(x = diff(DF$CRIMES, differences = d), n = 54, basefont = 15, main = paste0("ACF Plot of Central Philly Weekly Crimes\nDifference of Order ", d)),
			 ncol = 2)

p = 1:3
q. = 1:2

# lets look at seasonaly differenced weekly crimes to find values for D

par(mfrow = c(2,3))

plot(DF$CRIMES, main = "d = 0", type = "b")
lapply(1:5, function(i) plot(diff(DF$CRIMES, lag = 54, differences = i), main = paste0("D = ", i), type ="b"))

D = 1

# lets look at the auto arima parameter values for modeling weekly crimes

arima.auto = auto.arima(ts(DF$CRIMES, frequency = 54))
arima.auto

# creating a DOE that varies P will returns an error code: vmmin is not finite
# To address this, we will use an Arima.if() function which performs the same as Arima() but returns NA for error codes
# lets create our DOE based on the acf, pacf, and auto.arima results

DOE = expand.grid("p" = 1:3, 
				  "d" = 1, 
				  "q" = 1:2,
				  "P" = 0:1,
				  "D" = 0:1,
				  "Q" = 0:1)

DOE = rbind(c(1, 1, 1, 0, 0, 1), DOE)
DOE = DOE[!duplicated(DOE),]
rownames(DOE) = 1:NROW(DOE)

head(DOE)
NROW(DOE)

rm(d, D, p, q.)

# lets compute the actual and fitted values for each arima model in the DOE
# given that some models have D = 1, we will evaluate the actual v. fitted after season 1 ie. observations [55:565]
	# this is becuase the the arima model will use the actual values as the fitted values for season 1 to initialize the seasonal parameter(s)

arimas = lapply(1:NROW(DOE), function(i)
								 data.frame("actual" = DF$CRIMES[55:565],
											"fit" = as.numeric(fitted.if(Arima.if(ts(DF$CRIMES, frequency = 54), 
																					order = c(DOE$p[i], DOE$d[i], DOE$q[i]), 
																					seasonal = list(order = c(DOE$P[i], DOE$D[i], DOE$Q[i]), period = 54),
																					optim.control = list(maxit = 1000))))[55:565]))

# filter out invalid DOE scenarios

invalid = sapply(1:NROW(DOE), function(i) all(is.na(arimas[[i]]$fit)))

DOE = DOE[which(invalid == FALSE),]
rownames(DOE) = 1:NROW(DOE)

arimas = arimas[which(invalid == FALSE)]

rm(invalid)

# lower values are ideal

params = sapply(1:NROW(DOE), function(i) 
							 DOE$p[i] + DOE$q[i] + DOE$P[i] + DOE$Q[i])

# values greater than 0.05 are ideal 

ttests = sapply(1:NROW(DOE), function(i) 
							 t.test(arimas[[i]]$actual - arimas[[i]]$fit)$p.value)

# values close to zero are ideal 

stats = data.frame(t(sapply(1:NROW(DOE), function(i)
										 data.frame(accuracy(f = arimas[[i]]$fit, 
															 x = arimas[[i]]$actual), 
													row.names = 1))))

# lets add these results to the DOE

DOE$params = params
DOE$t.pval = ttests
DOE = cbind(DOE, stats)

cnames = colnames(DOE)
DOE[,7:13] = sapply(7:13, function(i) DOE[,i] = as.numeric(DOE[,i]))
colnames(DOE) = cnames

DOE

rm(arimas, params, ttests, stats, cnames)

# lets filter out the top models

# lets plot the histograms of these reponses to find filter limits

par(mfrow = c(2, 4))
lapply(7:13, function(i) hist(DOE[,i], main = colnames(DOE)[i]))

DOE

# apply filters

x = DOE[which(DOE$MAPE <= 14 &
			  DOE$RMSE <= 83 &
			  DOE$params <= 4),]

par(mfrow = c(2, 4))
lapply(7:13, function(i) hist(x[,i], main = colnames(x)[i]))

x

# lets compare the auto fit (scenario 1) with the chosen fit from the DOE (scenario 17)

DOE[c(1, 17),]

# lets build the final models

arimas = lapply(c(1, 17), function(i)
							Arima(ts(DF$CRIMES, frequency = 54), 
									order = c(DOE$p[i], DOE$d[i], DOE$q[i]), 
									seasonal = list(order = c(DOE$P[i], DOE$D[i], DOE$Q[i]), period = 54),					
									optim.control = list(maxit = 1000)))

# lets plot the model & residuals for the auto fit and chosen fit 

	# auto: scenario 1

plot_crimes_weekly_arimas1 = predplot(obj = arimas[[1]],
							  xlab = "Observation Number", 
							  ylab = "Crimes", 
							  main = "Crimes Per Week: Central Philly\nARIMA(1, 1, 1)X(0, 0, 1)54")

plot_crimes_weekly_arimas1

plot_doe_resid1 = residplots(actual = DF$CRIMES,
							 fit = as.numeric(fitted(arimas[[1]])),
							 histlabel.y = -5,
							 n = 54,
							 basefont = 15)

grid.arrange(plot_doe_resid1[[1]], 
			 plot_doe_resid1[[2]], 
			 plot_doe_resid1[[3]], 
			 plot_doe_resid1[[4]], 
			 plot_doe_resid1[[5]], 
			 plot_doe_resid1[[6]], 
			 ncol = 3)		  
		  
	# chosen: scenario 17

plot_crimes_weekly_arimas2 = predplot(obj = arimas[[2]],
							  xlab = "Observation Number", 
							  ylab = "Crimes", 
							  main = "Crimes Per Week: Central Philly\nARIMA(2, 1, 1)X(0, 0, 1)54")

plot_crimes_weekly_arimas2

plot_doe_resid2 = residplots(actual = DF$CRIMES,
							 fit = as.numeric(fitted(arimas[[2]])), 
							 histlabel.y = -5,
							 n = 54, 
							 basefont = 15)

grid.arrange(plot_doe_resid2[[1]], 
			 plot_doe_resid2[[2]], 
			 plot_doe_resid2[[3]], 
			 plot_doe_resid2[[4]], 
			 plot_doe_resid2[[5]], 
			 plot_doe_resid2[[6]], 
			 ncol = 3)		  

# lets store our results into lists

ARIMA = list()
weekly = list()
plots2 = list()
fits = list()
tables2 = list()

fits$crimes_weekly_CP = arimas[[2]]
plots2$plot_crimes_weekly_CP = plot_crimes_weekly_arimas2
plots2$plot_crimes_weekly_residual_CP = plot_doe_resid2
tables2$model_data_CP = DF
tables2$model_doe_CP = DOE

weekly$fits = fits
weekly$plots = plots2
weekly$tables = tables2

ARIMA$weekly = weekly
models$ARIMA = ARIMA

rm(x, ARIMA, weekly, fits, plots2, tables2, arimas, plot_crimes_weekly_arimas1, plot_crimes_weekly_arimas2, 
   plot_doe_resid1, plot_doe_resid2, DOE, DF, arima.auto)

# ---- forecast weekly crimes -------------------------------------------------------------------------------------------------------------------------------------------

DF = data.frame(tables$crime_Area1_tables$crimes_weekly[Area1 == "CP"])

# lets look at differenced weekly crimes to find values for d

par(mfrow = c(2,3))

plot(DF$CRIMES[1:478], main = "d = 0", type ="b")
lapply(1:5, function(i) plot(diff(DF$CRIMES, differences = i), main = paste0("d = ", i), type ="b"))

d = 1

# lets look at the acf and pacf plots of each differenced weekly crimes to find values of p and q

d = 1
grid.arrange(ggacf(x = diff(DF$CRIMES[1:478], differences = d), partial = TRUE, n = 54, basefont = 15, main = paste0("PACF Plot of Central Philly Weekly Crimes\nDifference of Order ", d)),
			 ggacf(x = diff(DF$CRIMES[1:478], differences = d), n = 54, basefont = 15, main = paste0("ACF Plot of Central Philly Weekly Crimes\nDifference of Order ", d)),
			 ncol = 2)

p = 1:3
q. = 1

# lets look at seasonaly differenced weekly crimes to find values for D

par(mfrow = c(2,3))

plot(DF$CRIMES, main = "d = 0", type ="b")
lapply(1:5, function(i) plot(diff(DF$CRIMES[1:478], lag = 54, differences = i), main = paste0("D = ", i), type ="b"))

D = 1

# lets look at the auto arima parameter values for modeling weekly crimes

arima.auto = auto.arima(ts(DF$CRIMES[1:478], frequency = 54))
arima.auto

# creating a DOE that varies P returns an error code: vmmin is not finite
	# I don't know how to fix this so the DOE process for arima modeling will have to ignore values of P (ie. P = 0)

# lets create our DOE based on the acf, pacf, and auto.arima results
# due to how the Arima() function fits P and Q coefficients, this DOE shouldn't be too large (ie. 10 or less scenarios)

DOE = expand.grid("p" = 1:3, 
				  "d" = 1, 
				  "q" = 1,
				  "P" = 0:1,
				  "D" = 0:1,
				  "Q" = 0:1)

DOE = rbind(c(1, 1, 1, 0, 0, 1), DOE)
DOE = DOE[!duplicated(DOE),]
rownames(DOE) = 1:NROW(DOE)

head(DOE)
NROW(DOE)

rm(d, D, p, q.)

# lets compute the actual and fitted values for the testing data of each arima model in the DOE

arimas = lapply(1:NROW(DOE), function(i)
								 data.frame("actual" = DF$CRIMES[479:565],
											"fit" = as.numeric(predict.if(Arima.if(ts(DF$CRIMES[1:478], frequency = 54), 
																			order = c(DOE$p[i], DOE$d[i], DOE$q[i]), 
																			seasonal = list(order = c(DOE$P[i], DOE$D[i], DOE$Q[i]), period = 54),
																			optim.control = list(maxit = 1000)),
																		n.ahead = 87,
																		se.fit = FALSE))))

# filter out invalid DOE scenarios

invalid = sapply(1:NROW(DOE), function(i) all(is.na(arimas[[i]]$fit)))

DOE = DOE[which(invalid == FALSE),]
rownames(DOE) = 1:NROW(DOE)

arimas = arimas[which(invalid == FALSE)]

rm(invalid)

# lower values are ideal

params = sapply(1:NROW(DOE), function(i) 
							 DOE$p[i] + DOE$q[i] + DOE$P[i] + DOE$Q[i])

# values greater than 0.05 are ideal 

ttests = sapply(1:NROW(DOE), function(i) 
							 t.test(arimas[[i]]$actual - arimas[[i]]$fit)$p.value)

# values close to zero are ideal 

stats = data.frame(t(sapply(1:NROW(DOE), function(i)
										 data.frame(accuracy(f = arimas[[i]]$fit, 
															 x = arimas[[i]]$actual), 
													row.names = 1))))

# lets add these results to the DOE

DOE$params = params
DOE$t.pval = ttests
DOE = cbind(DOE, stats)

cnames = colnames(DOE)
DOE[,7:13] = sapply(7:13, function(i) DOE[,i] = as.numeric(DOE[,i]))
colnames(DOE) = cnames

# lets filter out the top models

# lets plot the histograms of these reponses to find filter limits

par(mfrow = c(2, 4))
lapply(7:13, function(i) hist(DOE[,i], main = colnames(DOE)[i]))

DOE

# apply filters

x = DOE[which(DOE$MAPE <= 18),]

par(mfrow = c(2, 4))
lapply(7:13, function(i) hist(x[,i], main = colnames(x)[i]))

x

# lets compare the auto fit (scenario 1) with the chosen fit from the DOE (scenario 13)

DOE[c(1, 13),]

# lets build the final models

arimas = lapply(c(1, 13), function(i)
							Arima(ts(DF$CRIMES[1:478], frequency = 54), 
									order = c(DOE$p[i], DOE$d[i], DOE$q[i]), 
									seasonal = list(order = c(DOE$P[i], DOE$D[i], DOE$Q[i]), period = 54),
									optim.control = list(maxit = 1000)))

# lets plot the model & residuals for the auto fit and chosen fit 

	# auto: scenario 1

plot_crimes_weekly_arimas1 = predplot(obj = arimas[[1]],
							  n.ahead = 87,
							  actuals = DF$CRIMES,
							  xlab = "Observation Number", 
							  ylab = "Crimes", 
							  main = "Crimes Per Week: Central Philly\nARIMA(1, 1, 1)X(0, 0, 1)54")

plot_crimes_weekly_arimas1

plot_doe_resid1 = residplots(actual = DF$CRIMES[479:565],
							 fit = as.numeric(predict(arimas[[1]], n.ahead = 87, se.fit = FALSE)),
							 histlabel.y = -1,
							 n = 54,
							 binwidth = 30,
							 basefont = 15)

grid.arrange(plot_doe_resid1[[1]], 
			 plot_doe_resid1[[2]], 
			 plot_doe_resid1[[3]], 
			 plot_doe_resid1[[4]], 
			 plot_doe_resid1[[5]], 
			 plot_doe_resid1[[6]], 
			 ncol = 3)		  
		  
	# chosen: scenario 13

plot_crimes_weekly_arimas2 = predplot(obj = arimas[[2]],
							  n.ahead = 87,
							  actuals = DF$CRIMES,
							  xlab = "Observation Number", 
							  ylab = "Crimes", 
							  main = "Crimes Per Week: Central Philly\nARIMA(2, 1, 1)X(0, 1, 1)54")

plot_crimes_weekly_arimas2

plot_doe_resid2 = residplots(actual = DF$CRIMES[479:565],
							 fit = as.numeric(predict(arimas[[2]], n.ahead = 87, se.fit = FALSE)),
							 histlabel.y = -0.5,
							 n = 54,
							 binwidth = 30,
							 basefont = 15)

grid.arrange(plot_doe_resid2[[1]], 
			 plot_doe_resid2[[2]], 
			 plot_doe_resid2[[3]], 
			 plot_doe_resid2[[4]], 
			 plot_doe_resid2[[5]], 
			 plot_doe_resid2[[6]], 
			 ncol = 3)		  

# lets store our results into lists

models$ARIMA$weekly$fits$crimes_weekly_CP_forecast = arimas[[2]]
models$ARIMA$weekly$plots$plot_crimes_weekly_CP_forecast = plot_crimes_weekly_arimas2
models$ARIMA$weekly$plots$plot_crimes_weekly_residual_CP_forecast = plot_doe_resid2
models$ARIMA$weekly$tables$model_data_CP_forecast = DF
models$ARIMA$weekly$tables$model_doe_CP_forecast = DOE		  

rm(x, arimas, plot_crimes_weekly_arimas1, plot_crimes_weekly_arimas2, 
   plot_doe_resid1, plot_doe_resid2, DOE, DF, arima.auto)		  

}

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ---- ARIMAX Forecasting: West Philly -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

{

# ---- fit weekly crimes -------------------------------------------------------------------------------------------------------------------------------------------

DF = data.frame(tables$crime_Area1_tables$crimes_weekly[Area1 == "WP"])

# lets look at differenced weekly crimes to find values for d

par(mfrow = c(2,3))

plot(DF$CRIMES, main = "d = 0", type = "b")
lapply(1:5, function(i) plot(diff(DF$CRIMES, differences = i), main = paste0("d = ", i), type ="b"))

d = 1

# lets look at the acf and pacf plots of the differenced weekly crimes to find values of p and q

d = 1
grid.arrange(ggacf(x = diff(DF$CRIMES, differences = d), partial = TRUE, n = 54, basefont = 15, main = paste0("PACF Plot of West Philly Weekly Crimes\nDifference of Order ", d)),
			 ggacf(x = diff(DF$CRIMES, differences = d), n = 54, basefont = 15, main = paste0("ACF Plot of West Philly Weekly Crimes\nDifference of Order ", d)),
			 ncol = 2)

p = 1:3
q. = 1:2

# lets look at seasonaly differenced weekly crimes to find values for D

par(mfrow = c(2,3))

plot(DF$CRIMES, main = "d = 0", type = "b")
lapply(1:5, function(i) plot(diff(DF$CRIMES, lag = 54, differences = i), main = paste0("D = ", i), type ="b"))

D = 1

# lets look at the auto arima parameter values for modeling weekly crimes

arima.auto = auto.arima(ts(DF$CRIMES, frequency = 54))
arima.auto

# creating a DOE that varies P will returns an error code: vmmin is not finite
# To address this, we will use an Arima.if() function which performs the same as Arima() but returns NA for error codes
# lets create our DOE based on the acf, pacf, and auto.arima results

DOE = expand.grid("p" = 1:3, 
				  "d" = 1, 
				  "q" = 1:2,
				  "P" = 0:1,
				  "D" = 0:1,
				  "Q" = 0:1)

DOE = rbind(c(1, 0, 1, 0, 0, 1), DOE)
DOE = DOE[!duplicated(DOE),]
rownames(DOE) = 1:NROW(DOE)

head(DOE)
NROW(DOE)

rm(d, D, p, q.)

# lets compute the actual and fitted values for each arima model in the DOE
# given that some models have D = 1, we will evaluate the actual v. fitted after season 1 ie. observations [55:565]
	# this is becuase the the arima model will use the actual values as the fitted values for season 1 to initialize the seasonal parameter(s)
	
arimas = lapply(1:NROW(DOE), function(i)
								 data.frame("actual" = DF$CRIMES[55:565],
											"fit" = as.numeric(fitted.if(Arima.if(ts(DF$CRIMES, frequency = 54), 
																					order = c(DOE$p[i], DOE$d[i], DOE$q[i]), 
																					seasonal = list(order = c(DOE$P[i], DOE$D[i], DOE$Q[i]), period = 54),
																					optim.control = list(maxit = 1000))))[55:565]))

# filter out invalid DOE scenarios

invalid = sapply(1:NROW(DOE), function(i) all(is.na(arimas[[i]]$fit)))

DOE = DOE[which(invalid == FALSE),]
rownames(DOE) = 1:NROW(DOE)

arimas = arimas[which(invalid == FALSE)]

rm(invalid)

# lower values are ideal

params = sapply(1:NROW(DOE), function(i) 
							 DOE$p[i] + DOE$q[i] + DOE$P[i] + DOE$Q[i])

# values greater than 0.05 are ideal 

ttests = sapply(1:NROW(DOE), function(i) 
							 t.test(arimas[[i]]$actual - arimas[[i]]$fit)$p.value)

# values close to zero are ideal 

stats = data.frame(t(sapply(1:NROW(DOE), function(i)
										 data.frame(accuracy(f = arimas[[i]]$fit, 
															 x = arimas[[i]]$actual), 
													row.names = 1))))

# lets add these results to the DOE

DOE$params = params
DOE$t.pval = ttests
DOE = cbind(DOE, stats)

cnames = colnames(DOE)
DOE[,7:13] = sapply(7:13, function(i) DOE[,i] = as.numeric(DOE[,i]))
colnames(DOE) = cnames

DOE

rm(arimas, params, ttests, stats, cnames)

# lets filter out the top models

# lets plot the histograms of these reponses to find filter limits

par(mfrow = c(2, 4))
lapply(7:13, function(i) hist(DOE[,i], main = colnames(DOE)[i]))

DOE

# apply filters

x = DOE[which(DOE$RMSE <= 85 &
			  DOE$MAE <= 58.6 & 
			  DOE$params <= 5),]

par(mfrow = c(2, 4))
lapply(7:13, function(i) hist(x[,i], main = colnames(x)[i]))

x

# lets compare the auto fit (scenario 1) with the chosen fit from the DOE (scenario 8)

DOE[c(1, 8),]

# lets build the final models

arimas = lapply(c(1), function(i)
							Arima(ts(DF$CRIMES, frequency = 54), 
									order = c(DOE$p[i], DOE$d[i], DOE$q[i]), 
									seasonal = list(order = c(DOE$P[i], DOE$D[i], DOE$Q[i]), period = 54),					
									optim.control = list(maxit = 1000)))

# lets plot the model & residuals for the auto fit and chosen fit 

	# auto: scenario 1

plot_crimes_weekly_arimas1 = predplot(obj = arimas[[1]],
							  xlab = "Observation Number", 
							  ylab = "Crimes", 
							  main = "Crimes Per Week: West Philly\nARIMA(1, 0, 1)X(0, 0, 1)54")

plot_crimes_weekly_arimas1

plot_doe_resid1 = residplots(actual = DF$CRIMES,
							 fit = as.numeric(fitted(arimas[[1]])),
							 histlabel.y = -5,
							 n = 54,
							 basefont = 15)

grid.arrange(plot_doe_resid1[[1]], 
			 plot_doe_resid1[[2]], 
			 plot_doe_resid1[[3]], 
			 plot_doe_resid1[[4]], 
			 plot_doe_resid1[[5]], 
			 plot_doe_resid1[[6]], 
			 ncol = 3)		  
		  
	# chosen: scenario 8

plot_crimes_weekly_arimas2 = predplot(obj = arimas[[2]],
							  xlab = "Observation Number", 
							  ylab = "Crimes", 
							  main = "Crimes Per Week: West Philly\nARIMA(3, 1, 1)X(1, 0, 0)54")

plot_crimes_weekly_arimas2

plot_doe_resid2 = residplots(actual = DF$CRIMES,
							 fit = as.numeric(fitted(arimas[[2]])), 
							 histlabel.y = -5,
							 n = 54, 
							 basefont = 15)

grid.arrange(plot_doe_resid2[[1]], 
			 plot_doe_resid2[[2]], 
			 plot_doe_resid2[[3]], 
			 plot_doe_resid2[[4]], 
			 plot_doe_resid2[[5]], 
			 plot_doe_resid2[[6]], 
			 ncol = 3)		  

# lets store our results into lists

models$ARIMA$weekly$fits$crimes_weekly_WP = arimas[[1]]
models$ARIMA$weekly$plots$plot_crimes_weekly_WP = plot_crimes_weekly_arimas1
models$ARIMA$weekly$plots$plot_crimes_weekly_residual_WP = plot_doe_resid1
models$ARIMA$weekly$tables$model_data_WP = DF
models$ARIMA$weekly$tables$model_doe_WP = DOE		  

rm(x, arimas, plot_crimes_weekly_arimas1, plot_crimes_weekly_arimas2, 
   plot_doe_resid1, plot_doe_resid2, DOE, DF, arima.auto)		  

# ---- forecast weekly crimes -------------------------------------------------------------------------------------------------------------------------------------------

DF = data.frame(tables$crime_Area1_tables$crimes_weekly[Area1 == "WP"])

# lets look at differenced weekly crimes to find values for d

par(mfrow = c(2,3))

plot(DF$CRIMES[1:478], main = "d = 0", type ="b")
lapply(1:5, function(i) plot(diff(DF$CRIMES, differences = i), main = paste0("d = ", i), type ="b"))

d = 1

# lets look at the acf and pacf plots of each differenced weekly crimes to find values of p and q

d = 1
grid.arrange(ggacf(x = diff(DF$CRIMES[1:478], differences = d), partial = TRUE, n = 54, basefont = 15, main = paste0("PACF Plot of West Philly Weekly Crimes\nDifference of Order ", d)),
			 ggacf(x = diff(DF$CRIMES[1:478], differences = d), n = 54, basefont = 15, main = paste0("ACF Plot of West Philly Weekly Crimes\nDifference of Order ", d)),
			 ncol = 2)

p = 1:3
q. = 1

# lets look at seasonaly differenced weekly crimes to find values for D

par(mfrow = c(2,3))

plot(DF$CRIMES, main = "d = 0", type ="b")
lapply(1:5, function(i) plot(diff(DF$CRIMES[1:478], lag = 54, differences = i), main = paste0("D = ", i), type ="b"))

D = 1

# lets look at the auto arima parameter values for modeling weekly crimes

arima.auto = auto.arima(ts(DF$CRIMES[1:478], frequency = 54))
arima.auto

# creating a DOE that varies P will returns an error code: vmmin is not finite
# To address this, we will use an Arima.if() function which performs the same as Arima() but returns NA for error codes
# lets create our DOE based on the acf, pacf, and auto.arima results

DOE = expand.grid("p" = 1:3, 
				  "d" = 1, 
				  "q" = 1,
				  "P" = 0:1,
				  "D" = 0:1,
				  "Q" = 0:1)

DOE = rbind(c(1, 0, 1, 0, 0, 1), DOE)
DOE = DOE[!duplicated(DOE),]
rownames(DOE) = 1:NROW(DOE)

head(DOE)
NROW(DOE)

rm(d, D, p, q.)

# lets compute the actual and fitted values for the testing data of each arima model in the DOE

arimas = lapply(1:NROW(DOE), function(i)
								 data.frame("actual" = DF$CRIMES[479:565],
											"fit" = as.numeric(predict.if(Arima.if(ts(DF$CRIMES[1:478], frequency = 54), 
																			order = c(DOE$p[i], DOE$d[i], DOE$q[i]), 
																			seasonal = list(order = c(DOE$P[i], DOE$D[i], DOE$Q[i]), period = 54),
																			optim.control = list(maxit = 1000)),
																		n.ahead = 87,
																		se.fit = FALSE))))

# filter out invalid DOE scenarios

invalid = sapply(1:NROW(DOE), function(i) all(is.na(arimas[[i]]$fit)))

DOE = DOE[which(invalid == FALSE),]
rownames(DOE) = 1:NROW(DOE)

arimas = arimas[which(invalid == FALSE)]

rm(invalid)

# lower values are ideal

params = sapply(1:NROW(DOE), function(i) 
							 DOE$p[i] + DOE$q[i] + DOE$P[i] + DOE$Q[i])

# values greater than 0.05 are ideal 

ttests = sapply(1:NROW(DOE), function(i) 
							 t.test(arimas[[i]]$actual - arimas[[i]]$fit)$p.value)

# values close to zero are ideal 

stats = data.frame(t(sapply(1:NROW(DOE), function(i)
										 data.frame(accuracy(f = arimas[[i]]$fit, 
															 x = arimas[[i]]$actual), 
													row.names = 1))))

# lets add these results to the DOE

DOE$params = params
DOE$t.pval = ttests
DOE = cbind(DOE, stats)

cnames = colnames(DOE)
DOE[,7:13] = sapply(7:13, function(i) DOE[,i] = as.numeric(DOE[,i]))
colnames(DOE) = cnames

DOE

rm(arimas, params, ttests, stats, cnames)

# lets filter out the top models

# lets plot the histograms of these reponses to find filter limits

par(mfrow = c(2, 4))
lapply(7:13, function(i) hist(DOE[,i], main = colnames(DOE)[i]))

DOE

# apply filters

x = DOE[which(DOE$t.pval >= 0.05 &
			  DOE$MAPE <= 18),]

par(mfrow = c(2, 4))
lapply(7:13, function(i) hist(x[,i], main = colnames(x)[i]))

x

# lets compare the auto fit (scenario 1) with the chosen fit from the DOE (scenario 13)

DOE[c(1, 16),]

# lets build the final models

arimas = lapply(c(1, 16), function(i)
							Arima(ts(DF$CRIMES[1:478], frequency = 54), 
									order = c(DOE$p[i], DOE$d[i], DOE$q[i]), 
									seasonal = list(order = c(DOE$P[i], DOE$D[i], DOE$Q[i]), period = 54),
									optim.control = list(maxit = 1000)))

# lets plot the model & residuals for the auto fit and chosen fit 

	# auto: scenario 1

plot_crimes_weekly_arimas1 = predplot(obj = arimas[[1]],
							  n.ahead = 87,
							  actuals = DF$CRIMES,
							  xlab = "Observation Number", 
							  ylab = "Crimes", 
							  main = "Crimes Per Week: West Philly\nARIMA(1, 0, 1)X(0, 0, 1)54")

plot_crimes_weekly_arimas1

plot_doe_resid1 = residplots(actual = DF$CRIMES[479:565],
							 fit = as.numeric(predict(arimas[[1]], n.ahead = 87, se.fit = FALSE)),
							 histlabel.y = -1,
							 n = 54,
							 binwidth = 40,
							 basefont = 15)

grid.arrange(plot_doe_resid1[[1]], 
			 plot_doe_resid1[[2]], 
			 plot_doe_resid1[[3]], 
			 plot_doe_resid1[[4]], 
			 plot_doe_resid1[[5]], 
			 plot_doe_resid1[[6]], 
			 ncol = 3)		  
		  
	# chosen: scenario 13

plot_crimes_weekly_arimas2 = predplot(obj = arimas[[2]],
							  n.ahead = 87,
							  actuals = DF$CRIMES,
							  xlab = "Observation Number", 
							  ylab = "Crimes", 
							  main = "Crimes Per Week: West Philly\nARIMA(3, 1, 1)X(0, 1, 1)54")

plot_crimes_weekly_arimas2

plot_doe_resid2 = residplots(actual = DF$CRIMES[479:565],
							 fit = as.numeric(predict(arimas[[2]], n.ahead = 87, se.fit = FALSE)),
							 histlabel.y = -0.5,
							 n = 54,
							 binwidth = 40,
							 basefont = 15)

grid.arrange(plot_doe_resid2[[1]], 
			 plot_doe_resid2[[2]], 
			 plot_doe_resid2[[3]], 
			 plot_doe_resid2[[4]], 
			 plot_doe_resid2[[5]], 
			 plot_doe_resid2[[6]], 
			 ncol = 3)		  

# lets store our results into lists

models$ARIMA$weekly$fits$crimes_weekly_WP_forecast = arimas[[2]]
models$ARIMA$weekly$plots$plot_crimes_weekly_WP_forecast = plot_crimes_weekly_arimas2
models$ARIMA$weekly$plots$plot_crimes_weekly_residual_WP_forecast = plot_doe_resid2
models$ARIMA$weekly$tables$model_data_WP_forecast = DF
models$ARIMA$weekly$tables$model_doe_WP_forecast = DOE		  

rm(x, arimas, plot_crimes_weekly_arimas1, plot_crimes_weekly_arimas2, 
   plot_doe_resid1, plot_doe_resid2, DOE, DF, arima.auto)		  

}

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ---- ARIMAX Forecasting: North Philly -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

{

# ---- fit weekly crimes -------------------------------------------------------------------------------------------------------------------------------------------

DF = data.frame(tables$crime_Area1_tables$crimes_weekly[Area1 == "NP"])

# lets look at differenced weekly crimes to find values for d

par(mfrow = c(2,3))

plot(DF$CRIMES, main = "d = 0", type = "b")
lapply(1:5, function(i) plot(diff(DF$CRIMES, differences = i), main = paste0("d = ", i), type ="b"))

d = 1

# lets look at the acf and pacf plots of the differenced weekly crimes to find values of p and q

d = 1
grid.arrange(ggacf(x = diff(DF$CRIMES, differences = d), partial = TRUE, n = 54, basefont = 15, main = paste0("PACF Plot of North Philly Weekly Crimes\nDifference of Order ", d)),
			 ggacf(x = diff(DF$CRIMES, differences = d), n = 54, basefont = 15, main = paste0("ACF Plot of North Philly Weekly Crimes\nDifference of Order ", d)),
			 ncol = 2)

p = 1:3
q. = 1:2

# lets look at seasonaly differenced weekly crimes to find values for D

par(mfrow = c(2,3))

plot(DF$CRIMES, main = "d = 0", type = "b")
lapply(1:5, function(i) plot(diff(DF$CRIMES, lag = 54, differences = i), main = paste0("D = ", i), type ="b"))

D = 1

# lets look at the auto arima parameter values for modeling weekly crimes

arima.auto = auto.arima(ts(DF$CRIMES, frequency = 54))
arima.auto

# creating a DOE that varies P will returns an error code: vmmin is not finite
# To address this, we will use an Arima.if() function which performs the same as Arima() but returns NA for error codes
# lets create our DOE based on the acf, pacf, and auto.arima results

DOE = expand.grid("p" = 1:3, 
				  "d" = 1, 
				  "q" = 1:2,
				  "P" = 0:1,
				  "D" = 0:1,
				  "Q" = 0:1)

DOE = rbind(c(1, 1, 2, 0, 0, 2), DOE)
head(DOE)
NROW(DOE)

rm(d, D, p, q.)

# lets compute the actual and fitted values for each arima model in the DOE
# given that some models have D = 1, we will evaluate the actual v. fitted after season 1 ie. observations [55:565]
	# this is becuase the the arima model will use the actual values as the fitted values for season 1 to initialize the seasonal parameter(s)

arimas = lapply(1:NROW(DOE), function(i)
								 data.frame("actual" = DF$CRIMES[55:565],
											"fit" = as.numeric(fitted.if(Arima.if(ts(DF$CRIMES, frequency = 54), 
																					order = c(DOE$p[i], DOE$d[i], DOE$q[i]), 
																					seasonal = list(order = c(DOE$P[i], DOE$D[i], DOE$Q[i]), period = 54),
																					optim.control = list(maxit = 1000))))[55:565]))

# filter out invalid DOE scenarios

invalid = sapply(1:NROW(DOE), function(i) all(is.na(arimas[[i]]$fit)))

DOE = DOE[which(invalid == FALSE),]
rownames(DOE) = 1:NROW(DOE)

arimas = arimas[which(invalid == FALSE)]

rm(invalid)

# lower values are ideal

params = sapply(1:NROW(DOE), function(i) 
							 DOE$p[i] + DOE$q[i] + DOE$P[i] + DOE$Q[i])

# values greater than 0.05 are ideal 

ttests = sapply(1:NROW(DOE), function(i) 
							 t.test(arimas[[i]]$actual - arimas[[i]]$fit)$p.value)

# values close to zero are ideal 

stats = data.frame(t(sapply(1:NROW(DOE), function(i)
										 data.frame(accuracy(f = arimas[[i]]$fit, 
															 x = arimas[[i]]$actual), 
													row.names = 1))))

# lets add these results to the DOE

DOE$params = params
DOE$t.pval = ttests
DOE = cbind(DOE, stats)

cnames = colnames(DOE)
DOE[,7:13] = sapply(7:13, function(i) DOE[,i] = as.numeric(DOE[,i]))
colnames(DOE) = cnames

DOE

rm(arimas, params, ttests, stats, cnames)

# lets filter out the top models

# lets plot the histograms of these reponses to find filter limits

par(mfrow = c(2, 4))
lapply(7:13, function(i) hist(DOE[,i], main = colnames(DOE)[i]))

DOE

# apply filters

x = DOE[which(DOE$RMSE <= 206 &
			  DOE$params <= 5 &
			  DOE$MAPE <= 12.34),]

par(mfrow = c(2, 4))
lapply(7:13, function(i) hist(x[,i], main = colnames(x)[i]))

x

# lets compare the auto fit (scenario 1) with the chosen fit from the DOE (scenario 10)

DOE[c(1, 10),]

# lets build the final models

arimas = lapply(c(1, 10), function(i)
							Arima(ts(DF$CRIMES, frequency = 54), 
									order = c(DOE$p[i], DOE$d[i], DOE$q[i]), 
									seasonal = list(order = c(DOE$P[i], DOE$D[i], DOE$Q[i]), period = 54),					
									optim.control = list(maxit = 1000)))

# lets plot the model & residuals for the auto fit and chosen fit 

	# auto: scenario 1

plot_crimes_weekly_arimas1 = predplot(obj = arimas[[1]],
							  xlab = "Observation Number", 
							  ylab = "Crimes", 
							  main = "Crimes Per Week: North Philly\nARIMA(1, 1, 2)X(0, 0, 2)54")

plot_crimes_weekly_arimas1

plot_doe_resid1 = residplots(actual = DF$CRIMES,
							 fit = as.numeric(fitted(arimas[[1]])),
							 histlabel.y = -6,
							 n = 54,
							 basefont = 15)

grid.arrange(plot_doe_resid1[[1]], 
			 plot_doe_resid1[[2]], 
			 plot_doe_resid1[[3]], 
			 plot_doe_resid1[[4]], 
			 plot_doe_resid1[[5]], 
			 plot_doe_resid1[[6]], 
			 ncol = 3)		  
		  
	# chosen: scenario 10

plot_crimes_weekly_arimas2 = predplot(obj = arimas[[2]],
							  xlab = "Observation Number", 
							  ylab = "Crimes", 
							  main = "Crimes Per Week: North Philly\nARIMA(2, 1, 2)X(1, 0, 0)54")

plot_crimes_weekly_arimas2

plot_doe_resid2 = residplots(actual = DF$CRIMES,
							 fit = as.numeric(fitted(arimas[[2]])), 
							 histlabel.y = -6,
							 n = 54, 
							 basefont = 15)

grid.arrange(plot_doe_resid2[[1]], 
			 plot_doe_resid2[[2]], 
			 plot_doe_resid2[[3]], 
			 plot_doe_resid2[[4]], 
			 plot_doe_resid2[[5]], 
			 plot_doe_resid2[[6]], 
			 ncol = 3)		  

# lets store our results into lists

models$ARIMA$weekly$fits$crimes_weekly_NP = arimas[[2]]
models$ARIMA$weekly$plots$plot_crimes_weekly_NP = plot_crimes_weekly_arimas2
models$ARIMA$weekly$plots$plot_crimes_weekly_residual_NP = plot_doe_resid2
models$ARIMA$weekly$tables$model_data_NP = DF
models$ARIMA$weekly$tables$model_doe_NP = DOE		  

rm(x, arimas, plot_crimes_weekly_arimas1, plot_crimes_weekly_arimas2, 
   plot_doe_resid1, plot_doe_resid2, DOE, DF, arima.auto)		  

# ---- forecast weekly crimes -------------------------------------------------------------------------------------------------------------------------------------------

DF = data.frame(tables$crime_Area1_tables$crimes_weekly[Area1 == "NP"])

# lets look at differenced weekly crimes to find values for d

par(mfrow = c(2,3))

plot(DF$CRIMES[1:478], main = "d = 0", type ="b")
lapply(1:5, function(i) plot(diff(DF$CRIMES, differences = i), main = paste0("d = ", i), type ="b"))

d = 1

# lets look at the acf and pacf plots of each differenced weekly crimes to find values of p and q

d = 1
grid.arrange(ggacf(x = diff(DF$CRIMES[1:478], differences = d), partial = TRUE, n = 54, basefont = 15, main = paste0("PACF Plot of North Philly Weekly Crimes\nDifference of Order ", d)),
			 ggacf(x = diff(DF$CRIMES[1:478], differences = d), n = 54, basefont = 15, main = paste0("ACF Plot of North Philly Weekly Crimes\nDifference of Order ", d)),
			 ncol = 2)

p = 1:3
q. = 1:2

# lets look at seasonaly differenced weekly crimes to find values for D

par(mfrow = c(2,3))

plot(DF$CRIMES, main = "d = 0", type ="b")
lapply(1:5, function(i) plot(diff(DF$CRIMES[1:478], lag = 54, differences = i), main = paste0("D = ", i), type ="b"))

D = 1

# lets look at the auto arima parameter values for modeling weekly crimes

arima.auto = auto.arima(ts(DF$CRIMES[1:478], frequency = 54))
arima.auto

# creating a DOE that varies P will returns an error code: vmmin is not finite
# To address this, we will use an Arima.if() function which performs the same as Arima() but returns NA for error codes
# lets create our DOE based on the acf, pacf, and auto.arima results

DOE = expand.grid("p" = 1:3, 
				  "d" = 1, 
				  "q" = 1:2,
				  "P" = 0:1,
				  "D" = 0:1,
				  "Q" = 0:1)

DOE = rbind(c(1, 1, 2, 0, 0, 1), DOE)
DOE = DOE[!duplicated(DOE),]
rownames(DOE) = 1:NROW(DOE)

head(DOE)
NROW(DOE)

rm(d, D, p, q.)

# lets compute the actual and fitted values for the testing data of each arima model in the DOE

arimas = lapply(1:NROW(DOE), function(i)
								 data.frame("actual" = DF$CRIMES[479:565],
											"fit" = as.numeric(predict.if(Arima.if(ts(DF$CRIMES[1:478], frequency = 54), 
																			order = c(DOE$p[i], DOE$d[i], DOE$q[i]), 
																			seasonal = list(order = c(DOE$P[i], DOE$D[i], DOE$Q[i]), period = 54),
																			optim.control = list(maxit = 1000)),
																		n.ahead = 87,
																		se.fit = FALSE))))

# filter out invalid DOE scenarios

invalid = sapply(1:NROW(DOE), function(i) all(is.na(arimas[[i]]$fit)))

DOE = DOE[which(invalid == FALSE),]
rownames(DOE) = 1:NROW(DOE)

arimas = arimas[which(invalid == FALSE)]

rm(invalid)

# lower values are ideal

params = sapply(1:NROW(DOE), function(i) 
							 DOE$p[i] + DOE$q[i] + DOE$P[i] + DOE$Q[i])

# values greater than 0.05 are ideal 

ttests = sapply(1:NROW(DOE), function(i) 
							 t.test(arimas[[i]]$actual - arimas[[i]]$fit)$p.value)

# values close to zero are ideal 

stats = data.frame(t(sapply(1:NROW(DOE), function(i)
										 data.frame(accuracy(f = arimas[[i]]$fit, 
															 x = arimas[[i]]$actual), 
													row.names = 1))))

# lets add these results to the DOE

DOE$params = params
DOE$t.pval = ttests
DOE = cbind(DOE, stats)

cnames = colnames(DOE)
DOE[,7:13] = sapply(7:13, function(i) DOE[,i] = as.numeric(DOE[,i]))
colnames(DOE) = cnames

DOE

rm(arimas, params, ttests, stats, cnames)

# lets filter out the top models

# lets plot the histograms of these reponses to find filter limits

par(mfrow = c(2, 4))
lapply(7:13, function(i) hist(DOE[,i], main = colnames(DOE)[i]))

DOE

# apply filters

x = DOE[which(DOE$MAPE <= 15 &
			  DOE$RMSE <= 222 &
			  DOE$MAE <= 156 &
			  DOE$params <= 5),]

par(mfrow = c(2, 4))
lapply(7:13, function(i) hist(x[,i], main = colnames(x)[i]))

x

# lets compare the auto fit (scenario 1) with the chosen fit from the DOE (scenario 31)

DOE[c(1, 31),]

# lets build the final models

arimas = lapply(c(1, 31), function(i)
							Arima(ts(DF$CRIMES[1:478], frequency = 54), 
									order = c(DOE$p[i], DOE$d[i], DOE$q[i]), 
									seasonal = list(order = c(DOE$P[i], DOE$D[i], DOE$Q[i]), period = 54),
									optim.control = list(maxit = 1000)))

# lets plot the model & residuals for the auto fit and chosen fit 

	# auto: scenario 1

plot_crimes_weekly_arimas1 = predplot(obj = arimas[[1]],
							  n.ahead = 87,
							  actuals = DF$CRIMES,
							  xlab = "Observation Number", 
							  ylab = "Crimes", 
							  main = "Crimes Per Week: North Philly\nARIMA(1, 1, 2)X(0, 0, 1)54")

plot_crimes_weekly_arimas1

plot_doe_resid1 = residplots(actual = DF$CRIMES[479:565],
							 fit = as.numeric(predict(arimas[[1]], n.ahead = 87, se.fit = FALSE)),
							 histlabel.y = -0.5,
							 n = 54,
							 binwidth = 60,
							 basefont = 15)

grid.arrange(plot_doe_resid1[[1]], 
			 plot_doe_resid1[[2]], 
			 plot_doe_resid1[[3]], 
			 plot_doe_resid1[[4]], 
			 plot_doe_resid1[[5]], 
			 plot_doe_resid1[[6]], 
			 ncol = 3)		  
		  
	# chosen: scenario 31

plot_crimes_weekly_arimas2 = predplot(obj = arimas[[2]],
							  n.ahead = 87,
							  actuals = DF$CRIMES,
							  xlab = "Observation Number", 
							  ylab = "Crimes", 
							  main = "Crimes Per Week: North Philly\nARIMA(2, 1, 1)X(0, 1, 1)54")

plot_crimes_weekly_arimas2

plot_doe_resid2 = residplots(actual = DF$CRIMES[479:565],
							 fit = as.numeric(predict(arimas[[2]], n.ahead = 87, se.fit = FALSE)),
							 histlabel.y = -0.5,
							 n = 54,
							 binwidth = 65,
							 basefont = 15)

grid.arrange(plot_doe_resid2[[1]], 
			 plot_doe_resid2[[2]], 
			 plot_doe_resid2[[3]], 
			 plot_doe_resid2[[4]], 
			 plot_doe_resid2[[5]], 
			 plot_doe_resid2[[6]], 
			 ncol = 3)		  

# lets store our results into lists

models$ARIMA$weekly$fits$crimes_weekly_NP_forecast = arimas[[2]]
models$ARIMA$weekly$plots$plot_crimes_weekly_NP_forecast = plot_crimes_weekly_arimas2
models$ARIMA$weekly$plots$plot_crimes_weekly_residual_NP_forecast = plot_doe_resid2
models$ARIMA$weekly$tables$model_data_NP_forecast = DF
models$ARIMA$weekly$tables$model_doe_NP_forecast = DOE		  

rm(x, arimas, plot_crimes_weekly_arimas1, plot_crimes_weekly_arimas2, 
   plot_doe_resid1, plot_doe_resid2, DOE, DF, arima.auto)		  

}

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ---- Week 11 Presentation Tables -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

{

# CP ----------------------------------------------------------------

	# exp fit: 51714 + 1
	
DOE = expand.grid(alpha = seq(0.01, 0.39, 0.01), 
				  beta = seq(0, 0.25, 0.01), 
				  gamma = seq(0.25, 0.75, 0.01))

	# arimax fit: 48 + 1

DOE = expand.grid("p" = 1:3, 
				  "d" = 1, 
				  "q" = 1:2,
				  "P" = 0:1,
				  "D" = 0:1,
				  "Q" = 0:1)


	# exp forecast: 34476 + 2
	
DOE = expand.grid(alpha = seq(0.1, 0.35, 0.01), 
				  beta = seq(0, 0.25, 0.01), 
				  gamma = seq(0.25, 0.75, 0.01))

	# arimax forecast: 24 + 1

DOE = expand.grid("p" = 1:3, 
				  "d" = 1, 
				  "q" = 1,
				  "P" = 0:1,
				  "D" = 0:1,
				  "Q" = 0:1)


# WP ----------------------------------------------------------------

	# exp fit: 49062 + 1
	
DOE = expand.grid(alpha = seq(0.01, 0.37, 0.01), 
				  beta = seq(0, 0.25, 0.01), 
				  gamma = seq(0.28, 0.78, 0.01))

	# arimax fit: 48 + 1

DOE = expand.grid("p" = 1:3, 
				  "d" = 1, 
				  "q" = 1:2,
				  "P" = 0:1,
				  "D" = 0:1,
				  "Q" = 0:1)

	# exp forecast: 67626 + 2
	
DOE = expand.grid(alpha = seq(0.01, 0.51, 0.01), 
				  beta = seq(0, 0.25, 0.01), 
				  gamma = seq(0.24, 0.74, 0.01))

	# arimax forecast: 24 + 1

DOE = expand.grid("p" = 1:3, 
				  "d" = 1, 
				  "q" = 1,
				  "P" = 0:1,
				  "D" = 0:1,
				  "Q" = 0:1)

# NP ----------------------------------------------------------------

	# exp fit: 50388 + 1
	
DOE = expand.grid(alpha = seq(0.01, 0.38, 0.01), 
				  beta = seq(0, 0.25, 0.01), 
				  gamma = seq(0.34, 0.84, 0.01))

	# arimax fit: 48 + 1

DOE = expand.grid("p" = 1:3, 
				  "d" = 1, 
				  "q" = 1:2,
				  "P" = 0:1,
				  "D" = 0:1,
				  "Q" = 0:1)

	# exp forecast: 45084 + 2
	
DOE = expand.grid(alpha = seq(0.01, 0.34, 0.01), 
				  beta = seq(0, 0.25, 0.01), 
				  gamma = seq(0.34, 0.84, 0.01))

	# arimax forecast: 48 + 1

DOE = expand.grid("p" = 1:3, 
				  "d" = 1, 
				  "q" = 1:2,
				  "P" = 0:1,
				  "D" = 0:1,
				  "Q" = 0:1)

}

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ---- Week 11 Presentation Graphics -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

{

# fit

i = 3

x1 = 500
y1 = 2250
x2 = 500
y2 = 2250

hw = c(models$HW$weekly$tables$model_doe_CP$MAPE[33517],
		models$HW$weekly$tables$model_doe_WP$MAPE[22157],
		models$HW$weekly$tables$model_doe_NP$MAPE[25701])[i]
		
arimax = c(models$ARIMA$weekly$tables$model_doe_CP$MAPE[17],
			models$ARIMA$weekly$tables$model_doe_WP$MAPE[8],
			models$ARIMA$weekly$tables$model_doe_NP$MAPE[10])[i]

p1 = models$HW$weekly$plots$plot_crimes_weekly_NP + annotate("text", size = 5, x = x1, y = y1, label = paste0("MAPE = ", round(hw,1), "%")) + theme_light(base_size = 14) + theme(legend.position = "none")
p2 = models$ARIMA$weekly$plots$plot_crimes_weekly_NP + annotate("text", size = 5, x = x2, y = y2, label = paste0("MAPE = ", round(arimax,1), "%")) + theme_light(base_size = 14) + theme(legend.position = "none")

grid.arrange(p1, p2, ncol = 1)

# forecast

i = 3

x1 = 500
y1 = 3000
x2 = 500
y2 = 3000

hw = c(models$HW$weekly$tables$model_doe_CP_forecast$MAPE[10352],
		models$HW$weekly$tables$model_doe_WP_forecast$MAPE[8827],
		models$HW$weekly$tables$model_doe_NP_forecast$MAPE[7257])[i]

arimax = c(models$ARIMA$weekly$tables$model_doe_CP_forecast$MAPE[13],
			models$ARIMA$weekly$tables$model_doe_WP_forecast$MAPE[16],
			models$ARIMA$weekly$tables$model_doe_NP_forecast$MAPE[31])[i]

p1 = models$HW$weekly$plots$plot_crimes_weekly_NP_forecast + scale_x_continuous(limits = c(479, 565)) + annotate("text", size = 5, x = x1, y = y1, label = paste0("MAPE = ", round(hw,1), "%")) + theme_light(base_size = 14) + theme(legend.position = "none")
p2 = models$ARIMA$weekly$plots$plot_crimes_weekly_NP_forecast + scale_x_continuous(limits = c(479, 565)) + annotate("text", size = 5, x = x2, y = y2, label = paste0("MAPE = ", round(arimax,1), "%")) + theme_light(base_size = 14) + theme(legend.position = "none")

grid.arrange(p1, p2, ncol = 1)

######################################################

rm(x1, x2, y1, y2, p1, p2, arimax, hw, i)

}

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ---- Dynamic Regression Forecasting: Central Philly -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

{

# ---- how to model with factor regression variables --------------------------------------------------------------------------------------------

# create some artifical data

dat = data.frame(value = rpois(49, 3000), 
				 day = factor(c("Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat"), 
								  levels = c("Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat")))

head(dat)

# create matrix of numeric predictors

xreg = data.frame(model.matrix(~dat$day))

# rename columns

colnames(xreg) = c("(Int)", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat")
head(xreg)

# remove intercept

xreg = xreg[,-1]
head(xreg)

# variable to be modelled

value.ts = ts(dat$value, frequency = 7)

# find ARIMAX model

mod = auto.arima(value.ts, xreg = xreg)
mod

rm(mod, value.ts, xreg, dat)

# ---- updating weekly crimes ---------------------------------------------------------------------------------------------------------------

# lets recreate the weekly crimes dataset with SUMMER as an external variable

dat = PHI[,.("CRIMES" = .N), by = .("Area1" = Area1, "YEAR" = Year, "WEEK" = Week)]

dat[, YEAR.WEEK := factor(paste0(dat[,YEAR], ".", dat[,WEEK]),
							 levels = paste0(dat[,YEAR], ".", dat[,WEEK]))]

summers = data.frame("YEAR" = as.numeric(2006:2016), 
					 "START" = as.numeric(strftime(c("2006-06-21",
													 "2007-06-21",
													 "2008-06-20",
													 "2009-06-21",
													 "2010-06-21",
													 "2011-06-21",
													 "2012-06-20",
													 "2013-06-21",
													 "2014-06-21",
													 "2015-06-21",
													 "2016-06-20"), format = "%W")), 
					 "END" = as.numeric(strftime(c("2006-09-23",
												   "2007-09-23",
												   "2008-09-22",
												   "2009-09-22",
												   "2010-09-22",
												   "2011-09-23",
												   "2012-09-22",
												   "2013-09-22",
												   "2014-09-22",
												   "2015-09-23",
												   "2016-09-22"), format = "%W")))

dat[, SUMMER := factor(as.numeric(sapply(1:NROW(dat), function(i)
													(as.numeric(dat$WEEK[i]) - 1) %in% (summers$START[which(as.numeric(dat$YEAR[i]) + 2005 == summers$YEAR)]:summers$END[which(as.numeric(dat$YEAR[i]) + 2005 == summers$YEAR)]))), levels = 0:1)]

dat[, DAYS := factor(sapply(1:NROW(dat), function(i)
									max(PHI[Year == substr(dat[i,YEAR.WEEK], 1, 4) & 
											Week == substr(dat[i,YEAR.WEEK], 6, 7)][,Dispatch_Date]) - 
									min(PHI[Year == substr(dat[i,YEAR.WEEK], 1, 4) & 
											Week == substr(dat[i,YEAR.WEEK], 6, 7)][,Dispatch_Date]) + 1), levels = 1:7)]

dat[, c("YEAR", "WEEK") := NULL]

setcolorder(dat, c(1, 3, 4, 5, 2))

tables$crime_Area1_tables$crimes_weekly = dat
rm(dat, summers)

# ---- fit weekly crimes -------------------------------------------------------------------------------------------------------------------------------------------

# lets set up our data set and external variables

DF = data.frame(tables$crime_Area1_tables$crimes_weekly[Area1 == "CP"])

# lets look at plots of the dataset

mycolors = colorRampPalette(brewer.pal(n = 8, name = "Dark2"))(9)

p1 = ggplot(data = DF, aes(x = DAYS, y = CRIMES, fill = DAYS, color = DAYS)) +
	 geom_violin(alpha = 1/4, width = 5/4) +
	 geom_jitter(alpha = 1/2, width = 1/6) +
	 scale_color_manual(values = mycolors[1:7]) + 
	 scale_fill_manual(values = mycolors[1:7]) +
	 labs(x = "Days", y = "Crimes") +
	 ggtitle("Central Philly\nCrimes Per Week v. Days") +
	 theme_light(base_size = 15) +
	 theme(legend.position = "none")

p2 = ggplot(data = DF, aes(x = SUMMER, y = CRIMES, fill = SUMMER, color = SUMMER)) +
	 geom_violin(alpha = 1/4, width = 3/4) +
	 geom_jitter(alpha = 1/2, width = 1/5) +
	 scale_color_manual(values = mycolors[8:9]) + 
	 scale_fill_manual(values = mycolors[8:9]) +
	 labs(x = "Summer", y = "Crimes") +
	 ggtitle("Central Philly\nCrimes Per Week v. Summer") +
	 theme_light(base_size = 15) +
	 theme(legend.position = "none")

grid.arrange(p1, p2, ncol = 2)

DF$SUMMER = as.numeric(DF$SUMMER) - 1
DF$DAYS = as.numeric(DF$DAYS)

rm(p1, p2, mycolors)

# lets look at the residuals of a linear model fit using the xreg variables

mod = lm(CRIMES ~ DAYS + SUMMER, data = DF)
mod

# lets look at differenced residuals to find values for d

par(mfrow = c(2,3))

plot(resid(mod), main = "d = 0", type = "b")
lapply(1:5, function(i) plot(diff(resid(mod), differences = i), main = paste0("d = ", i), type ="b"))

d = 1

# lets look at the acf and pacf plots of the differenced residuals to find values of p and q

d = 1
grid.arrange(ggacf(x = diff(resid(mod), differences = d), partial = TRUE, n = 54, basefont = 15, main = paste0("PACF Plot of Central Philly Residuals\nDifference of Order ", d)),
			 ggacf(x = diff(resid(mod), differences = d), n = 54, basefont = 15, main = paste0("ACF Plot of Central Philly Residuals\nDifference of Order ", d)),
			 ncol = 2)

p = 0:2
q. = 0:2

# lets look at seasonaly differenced residuals to find values for D

par(mfrow = c(2,3))

plot(resid(mod), main = "d = 0", type = "b")
lapply(1:5, function(i) plot(diff(resid(mod), lag = 54, differences = i), main = paste0("D = ", i), type ="b"))

D = 1

rm(d, p, q., D, mod)

# compare auto fits

arima.auto1 = auto.arima(ts(DF$CRIMES, frequency = 54), xreg = DF[,c("SUMMER", "DAYS")])
arima.auto1

c(0, 1, 1, 0, 0, 0)

arima.auto2 = auto.arima(ts(DF$CRIMES, frequency = 54), xreg = DF[,c("SUMMER")])
arima.auto2

c(1, 1, 1, 0, 0, 1)

arima.auto3 = auto.arima(ts(DF$CRIMES, frequency = 54), xreg = DF[,c("DAYS")])
arima.auto3

c(0, 1, 1, 0, 0, 0)

autofits = data.frame("Variables" = rbind("SUMMER,DAYS", "SUMMER", "DAYS"),
					  rbind(accuracy(arima.auto1), accuracy(arima.auto2), accuracy(arima.auto3)),
					  "AIC" = rbind(arima.auto1$aic, arima.auto2$aic, arima.auto3$aic),
					  "AICc" = rbind(arima.auto1$aicc, arima.auto2$aicc, arima.auto3$aicc),
					  "BIC" = rbind(arima.auto1$bic, arima.auto2$bic, arima.auto3$bic))

autofits = autofits[,-c(7, 8)]
autofits

rm(arima.auto1, arima.auto2, arima.auto3)

# lets create our DOE based on the auto fits, acf, pacf, and differencing plots

DOE = expand.grid("p" = 0:2, 
				  "d" = 1, 
				  "q" = 0:2,
				  "P" = 0:1,
				  "D" = 1,
				  "Q" = 0:1)

DOE = rbind(c(0, 1, 1, 0, 0, 0), DOE)
DOE = DOE[!duplicated(DOE),]
rownames(DOE) = 1:NROW(DOE)

head(DOE)
NROW(DOE)

# lets compute the actual and fitted values for each arima model in the DOE

arimas = lapply(1:NROW(DOE), function(i)
								 data.frame("actual" = DF$CRIMES[55:565],
											"fit" = as.numeric(fitted.if(Arima.if(ts(DF$CRIMES, frequency = 54), 
																			order = c(DOE$p[i], DOE$d[i], DOE$q[i]), 
																			seasonal = list(order = c(DOE$P[i], DOE$D[i], DOE$Q[i]), period = 54),
																			xreg = DF[,c("SUMMER", "DAYS")],
																			optim.control = list(maxit = 1000)))[55:565])))
	
# filter out invalid DOE scenarios

invalid = sapply(1:NROW(DOE), function(i) all(is.na(arimas[[i]]$fit)))

DOE = DOE[which(invalid == FALSE),]
rownames(DOE) = 1:NROW(DOE)

arimas = arimas[which(invalid == FALSE)]

rm(invalid)
	
# lower values are ideal

params = sapply(1:NROW(DOE), function(i) 
							 DOE$p[i] + DOE$q[i] + DOE$P[i] + DOE$Q[i])

# values greater than 0.05 are ideal 

ttests = sapply(1:NROW(DOE), function(i) 
							 t.test(arimas[[i]]$actual - arimas[[i]]$fit)$p.value)

# values close to zero are ideal 

stats = data.frame(t(sapply(1:NROW(DOE), function(i)
										 data.frame(accuracy(f = arimas[[i]]$fit, 
															 x = arimas[[i]]$actual), 
													row.names = 1))))

# lets add these results to the DOE

DOE$params = params
DOE$t.pval = ttests
DOE = cbind(DOE, stats)

cnames = colnames(DOE)
DOE[,7:13] = sapply(7:13, function(i) DOE[,i] = as.numeric(DOE[,i]))
colnames(DOE) = cnames

DOE

rm(arimas, params, ttests, stats, cnames)

# lets filter out the top models

# lets plot the histograms of these reponses to find filter limits

par(mfrow = c(2, 4))
lapply(7:13, function(i) hist(DOE[,i], main = colnames(DOE)[i]))

DOE

# apply filters

x = DOE[which(DOE$MAPE <= 9 &
			  DOE$RMSE <= 60 & 
			  DOE$params <= 3),]

par(mfrow = c(2, 4))
lapply(7:13, function(i) hist(x[,i], main = colnames(x)[i]))
x

# lets compare the auto fit (scenario 1) with the chosen fit from the DOE (scenario 21)

DOE[c(1, 21),]

# lets build the final models

arimas = lapply(c(1, 21), function(i)
							Arima(ts(DF$CRIMES, frequency = 54), 
									order = c(DOE$p[i], DOE$d[i], DOE$q[i]), 
									seasonal = list(order = c(DOE$P[i], DOE$D[i], DOE$Q[i]), period = 54),
									xreg = DF[,c("SUMMER", "DAYS")],						
									optim.control = list(maxit = 1000)))

# lets plot the model & residuals for the auto fit and chosen fit 

	# auto: scenario 1

plot_crimes_weekly_arimas1 = predplot(obj = arimas[[1]],
							  xlab = "Observation Number", 
							  ylab = "Crimes", 
							  main = "Crimes Per Week: Central Philly\nARIMA(0, 1, 1)X(0, 0, 0)54\nRegression Variables: 26.3 * SUMMER, 84.1 * DAYS")

plot_crimes_weekly_arimas1

plot_doe_resid1 = residplots(actual = DF$CRIMES[55:565],
							 fit = as.numeric(fitted(arimas[[1]])[55:565]),
							 histlabel.y = -2,
							 n = 54,
							 basefont = 15)

grid.arrange(plot_doe_resid1[[1]], 
			 plot_doe_resid1[[2]], 
			 plot_doe_resid1[[3]], 
			 plot_doe_resid1[[4]], 
			 plot_doe_resid1[[5]], 
			 plot_doe_resid1[[6]], 
			 ncol = 3)		  
		  
	# chosen: scenario 21

plot_crimes_weekly_arimas2 = predplot(obj = arimas[[2]],
							  xlab = "Observation Number", 
							  ylab = "Crimes", 
							  main = "Crimes Per Week: Central Philly\nARIMA(0, 1, 1)X(0, 1, 1)54\nRegression Variables: 7.9 * SUMMER, 81.4 * DAYS")

plot_crimes_weekly_arimas2

plot_doe_resid2 = residplots(actual = DF$CRIMES[55:565],
							 fit = as.numeric(fitted(arimas[[2]])[55:565]), 
							 histlabel.y = -2,
							 n = 54, 
							 basefont = 15)

grid.arrange(plot_doe_resid2[[1]], 
			 plot_doe_resid2[[2]], 
			 plot_doe_resid2[[3]], 
			 plot_doe_resid2[[4]], 
			 plot_doe_resid2[[5]], 
			 plot_doe_resid2[[6]], 
			 ncol = 3)		  

# lets store our results into lists

DR = list()
weekly = list()
plots2 = list()
fits = list()
tables2 = list()

fits$crimes_weekly_CP = arimas[[1]]
plots2$plot_crimes_weekly_CP = plot_crimes_weekly_arimas1
plots2$plot_crimes_weekly_residual_CP = plot_doe_resid1
tables2$model_data_CP = DF
tables2$model_doe_CP = DOE

weekly$fits = fits
weekly$plots = plots2
weekly$tables = tables2

DR$weekly = weekly
models$DR = DR

rm(x, DR, weekly, fits, plots2, tables2, arimas, plot_crimes_weekly_arimas1, plot_crimes_weekly_arimas2, 
   plot_doe_resid1, plot_doe_resid2, DOE, DF)

models$DR$weekly$tables$autofits_CP = autofits

rm(autofits)

# ---- forecast weekly crimes -------------------------------------------------------------------------------------------------------------------------------------------

# lets set up our data set and external variables

DF = data.frame(tables$crime_Area1_tables$crimes_weekly[Area1 == "CP"])
DF$SUMMER = as.numeric(DF$SUMMER) - 1
DF$DAYS = as.numeric(DF$DAYS)

# lets look at the residuals of a linear model fit using the xreg variables

mod = lm(CRIMES ~ DAYS + SUMMER, data = DF[1:478,])
mod

# lets look at differenced residuals to find values for d

par(mfrow = c(2,3))

plot(resid(mod), main = "d = 0", type = "b")
lapply(1:5, function(i) plot(diff(resid(mod), differences = i), main = paste0("d = ", i), type ="b"))

d = 1

# lets look at the acf and pacf plots of the differenced residuals to find values of p and q

d = 1
grid.arrange(ggacf(x = diff(resid(mod), differences = d), partial = TRUE, n = 54, basefont = 15, main = paste0("PACF Plot of Central Philly Residuals\nDifference of Order ", d)),
			 ggacf(x = diff(resid(mod), differences = d), n = 54, basefont = 15, main = paste0("ACF Plot of Central Philly Residuals\nDifference of Order ", d)),
			 ncol = 2)

p = 0:3
q. = 0:1

# lets look at seasonaly differenced residuals to find values for D

par(mfrow = c(2,3))

plot(resid(mod), main = "d = 0", type = "b")
lapply(1:5, function(i) plot(diff(resid(mod), lag = 54, differences = i), main = paste0("D = ", i), type ="b"))

D = 1

rm(d, p, q., D, mod)

# auto fit

arima.auto = auto.arima(ts(DF$CRIMES[1:478], frequency = 54), xreg = DF[1:478, c("SUMMER", "DAYS")])
arima.auto

c(0, 1, 1, 0, 0, 0)

rm(arima.auto)

# lets create our DOE based on the auto.arima, acf, pacf, and differencing plots

DOE = expand.grid("p" = 0:3, 
				  "d" = 1, 
				  "q" = 0:1,
				  "P" = 0:1,
				  "D" = 1,
				  "Q" = 0:1)

DOE = rbind(c(0, 1, 1, 0, 0, 0), DOE)
DOE = DOE[!duplicated(DOE),]
rownames(DOE) = 1:NROW(DOE)

head(DOE)
NROW(DOE)

# lets compute the arima models from the DOE

arimas = lapply(1:NROW(DOE), function(i)
								Arima.if(ts(DF$CRIMES[1:478], frequency = 54), 
											order = c(DOE$p[i], DOE$d[i], DOE$q[i]), 
											seasonal = list(order = c(DOE$P[i], DOE$D[i], DOE$Q[i]), period = 54),
											xreg = DF[1:478, c("SUMMER", "DAYS")],
											optim.control = list(maxit = 1000)))
		
# filter out invalid DOE scenarios

invalid = sapply(1:NROW(DOE), function(i) is.na(arimas[i]))

DOE = DOE[which(invalid == FALSE),]
rownames(DOE) = 1:NROW(DOE)

arimas = arimas[which(invalid == FALSE)]

rm(invalid)		

# lets compute the actual and fitted values for the testing data of each arima model in the DOE

arimas = lapply(1:NROW(DOE), function(i)
								data.frame("actual" = DF$CRIMES[479:565],
											"fit" = as.numeric(forecast(arimas[[i]],
																		h = 87,
																		xreg = DF[479:565, c("SUMMER", "DAYS")],
																		level = c(0, 0))$mean)))

# lower values are ideal

params = sapply(1:NROW(DOE), function(i) 
							 DOE$p[i] + DOE$q[i] + DOE$P[i] + DOE$Q[i])

# values greater than 0.05 are ideal 

ttests = sapply(1:NROW(DOE), function(i) 
							 t.test(arimas[[i]]$actual - arimas[[i]]$fit)$p.value)

# values close to zero are ideal 

stats = data.frame(t(sapply(1:NROW(DOE), function(i)
										 data.frame(accuracy(f = arimas[[i]]$fit, 
															 x = arimas[[i]]$actual), 
													row.names = 1))))

# lets add these results to the DOE

DOE$params = params
DOE$t.pval = ttests
DOE = cbind(DOE, stats)

cnames = colnames(DOE)
DOE[,7:13] = sapply(7:13, function(i) DOE[,i] = as.numeric(DOE[,i]))
colnames(DOE) = cnames

DOE

rm(arimas, params, ttests, stats, cnames)

# lets filter out the top models

# lets plot the histograms of these reponses to find filter limits

par(mfrow = c(2, 4))
lapply(7:13, function(i) hist(DOE[,i], main = colnames(DOE)[i]))

DOE

# apply filters

x = DOE[which(DOE$MAPE <= 12 &
			  DOE$t.pval > 0.2),]

par(mfrow = c(2, 4))
lapply(7:13, function(i) hist(x[,i], main = colnames(x)[i]))

x

# lets compare the auto fit (scenario 1) with the chosen fit from the DOE (scenario 17)

DOE[c(1, 17),]

# lets build the final models

arimas = lapply(c(1, 17), function(i)
							Arima(ts(DF$CRIMES[1:478], frequency = 54), 
									order = c(DOE$p[i], DOE$d[i], DOE$q[i]), 
									seasonal = list(order = c(DOE$P[i], DOE$D[i], DOE$Q[i]), period = 54),
									xreg = DF[1:478, c("SUMMER", "DAYS")],
									optim.control = list(maxit = 1000)))

# lets plot the model & residuals for the auto fit and chosen fit 

	# auto: scenario 1

plot_crimes_weekly_arimas1 = predplot(obj = arimas[[1]],
							  n.ahead = 87,
							  xreg = DF[479:565, c("SUMMER", "DAYS")],
							  actuals = DF$CRIMES,
							  xlab = "Observation Number", 
							  ylab = "Crimes", 
							  main = "Crimes Per Week: Central Philly\nARIMA(0, 1, 1)X(0, 0, 0)54\nRegression Variable: 31.5 * SUMMER, 86.6 * DAYS")

plot_crimes_weekly_arimas1

plot_doe_resid1 = residplots(actual = DF$CRIMES[479:565],
							 fit = as.numeric(predict(arimas[[1]], newxreg = DF[479:565, c("SUMMER", "DAYS")], n.ahead = 87, se.fit = FALSE)),
							 histlabel.y = -0.5,
							 n = 54,
							 binwidth = 25,
							 basefont = 15)

grid.arrange(plot_doe_resid1[[1]], 
			 plot_doe_resid1[[2]], 
			 plot_doe_resid1[[3]], 
			 plot_doe_resid1[[4]], 
			 plot_doe_resid1[[5]], 
			 plot_doe_resid1[[6]], 
			 ncol = 3)		  
		  
	# chosen: scenario 17

plot_crimes_weekly_arimas2 = predplot(obj = arimas[[2]],
							  n.ahead = 87,
							  xreg = DF[479:565, c("SUMMER", "DAYS")],
							  actuals = DF$CRIMES,
							  xlab = "Observation Number", 
							  ylab = "Crimes", 
							  main = "Crimes Per Week: Central Philly\nARIMA(1, 1, 0)X(0, 1, 1)54\nRegression Variable: 9.9 * SUMMER, 82.1 * DAYS")

plot_crimes_weekly_arimas2

plot_doe_resid2 = residplots(actual = DF$CRIMES[479:565],
							 fit = as.numeric(predict(arimas[[2]], newxreg = DF[479:565, c("SUMMER", "DAYS")], n.ahead = 87, se.fit = FALSE)),
							 histlabel.y = -0.5,
							 n = 54,
							 binwidth = 30,
							 basefont = 15)

grid.arrange(plot_doe_resid2[[1]], 
			 plot_doe_resid2[[2]], 
			 plot_doe_resid2[[3]], 
			 plot_doe_resid2[[4]], 
			 plot_doe_resid2[[5]], 
			 plot_doe_resid2[[6]], 
			 ncol = 3)		  

# lets store our results into lists

models$DR$weekly$fits$crimes_weekly_CP_forecast = arimas[[2]]
models$DR$weekly$plots$plot_crimes_weekly_CP_forecast = plot_crimes_weekly_arimas2
models$DR$weekly$plots$plot_crimes_weekly_residual_CP_forecast = plot_doe_resid2
models$DR$weekly$tables$model_data_CP_forecast = DF
models$DR$weekly$tables$model_doe_CP_forecast = DOE		  

rm(x, arimas, plot_crimes_weekly_arimas1, plot_crimes_weekly_arimas2, 
   plot_doe_resid1, plot_doe_resid2, DOE, DF)		  

}

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ---- Dynamic Regression Forecasting with Lagged Regression Variables: Central Philly -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

{

# ---- fit weekly crimes -------------------------------------------------------------------------------------------------------------------------------------------

# lets set up our data set and external variables
# we will compare three different lags for each of the two regression variables
	# a lag of 1, 4, and 54 weeks
	
DF = data.frame(tables$crime_Area1_tables$crimes_weekly[Area1 == "CP"])
DF$SUMMER = as.numeric(DF$SUMMER) - 1
DF$DAYS = as.numeric(DF$DAYS)

DF = data.table(DF)

DF[,SUMMER1 := shift(SUMMER, 1)]
DF[,SUMMER4 := shift(SUMMER, 4)]
DF[,SUMMER54 := shift(SUMMER, 54)]

DF[,DAYS1 := shift(DAYS, 1)]
DF[,DAYS4 := shift(DAYS, 4)]
DF[,DAYS54 := shift(DAYS, 54)]

setcolorder(DF, c(1, 2, 3, 6, 7, 8, 4, 9, 10, 11, 5))
DF = data.frame(DF)

head(DF)

# lets create a DOE to test all possible regression lag combinations
# we will choose the best combination based on the R^2-Pred, AIC, and BIC values of linear regression models

DOE = expand.grid("SUMMER" = 4:6,
				  "DAYS" = 8:10)

DOE$r2pred = sapply(1:NROW(DOE), function(i)
					statslm(lm(DF$CRIMES[55:565] ~ DF[55:565,DOE$DAYS[i]] + DF[55:565,DOE$SUMMER[i]]))$Prediction)

DOE = cbind(DOE, t(sapply(1:NROW(DOE), function(i)
							statslm(lm(DF$CRIMES[55:565] ~ DF[55:565,DOE$DAYS[i]] + DF[55:565,DOE$SUMMER[i]]))$Fitness)))

DOE$SUMMER = colnames(DF)[DOE$SUMMER]
DOE$DAYS = colnames(DF)[DOE$DAYS]

DOE

# update the dataset with the best lagged regression variables

DF = DF[,c(1,2,6,10,11)]
head(DF)

rm(DOE)

# lets look at the residuals of a linear model fit using the xreg variables

mod = lm(CRIMES ~ DAYS54 + SUMMER54, data = DF)
mod

# lets look at differenced residuals to find values for d

par(mfrow = c(2,3))

plot(resid(mod), main = "d = 0", type = "b")
lapply(1:5, function(i) plot(diff(resid(mod), differences = i), main = paste0("d = ", i), type ="b"))

d = 1

# lets look at the acf and pacf plots of the differenced residuals to find values of p and q

d = 1
grid.arrange(ggacf(x = diff(resid(mod), differences = d), partial = TRUE, n = 54, basefont = 15, main = paste0("PACF Plot of Central Philly Residuals\nDifference of Order ", d)),
			 ggacf(x = diff(resid(mod), differences = d), n = 54, basefont = 15, main = paste0("ACF Plot of Central Philly Residuals\nDifference of Order ", d)),
			 ncol = 2)

p = 0:4
q. = 0:1

# lets look at seasonaly differenced residuals to find values for D

par(mfrow = c(2,3))

plot(resid(mod), main = "d = 0", type = "b")
lapply(1:5, function(i) plot(diff(resid(mod), lag = 54, differences = i), main = paste0("D = ", i), type ="b"))

D = 1

rm(d, p, q., D, mod)

# lets look at the auto fit

arima.auto = auto.arima(ts(DF$CRIMES, frequency = 54), xreg = DF[,c("SUMMER54", "DAYS54")])
arima.auto

c(0, 1, 1, 1, 0, 0)

rm(arima.auto)

# lets create our DOE based on the auto fit, acf, pacf, and differencing plots

DOE = expand.grid("p" = 0:4, 
				  "d" = 1, 
				  "q" = 0:1,
				  "P" = 1,
				  "D" = 1,
				  "Q" = 0:1)

DOE = rbind(c(0, 1, 1, 1, 0, 0), DOE)
DOE = DOE[!duplicated(DOE),]
rownames(DOE) = 1:NROW(DOE)

head(DOE)
NROW(DOE)

# lets compute the actual and fitted values for each arima model in the DOE

arimas = lapply(1:NROW(DOE), function(i)
								 data.frame("actual" = DF$CRIMES[109:565],
											"fit" = as.numeric(fitted.if(Arima.if(ts(DF$CRIMES, frequency = 54), 
																			order = c(DOE$p[i], DOE$d[i], DOE$q[i]), 
																			seasonal = list(order = c(DOE$P[i], DOE$D[i], DOE$Q[i]), period = 54),
																			xreg = DF[,c("SUMMER54", "DAYS54")],
																			optim.control = list(maxit = 1000)))[109:565])))
	
# filter out invalid DOE scenarios

invalid = sapply(1:NROW(DOE), function(i) all(is.na(arimas[[i]]$fit)))

DOE = DOE[which(invalid == FALSE),]
rownames(DOE) = 1:NROW(DOE)

arimas = arimas[which(invalid == FALSE)]

rm(invalid)
	
# lower values are ideal

params = sapply(1:NROW(DOE), function(i) 
							 DOE$p[i] + DOE$q[i] + DOE$P[i] + DOE$Q[i])

# values greater than 0.05 are ideal 

ttests = sapply(1:NROW(DOE), function(i) 
							 t.test(arimas[[i]]$actual - arimas[[i]]$fit)$p.value)

# values close to zero are ideal 

stats = data.frame(t(sapply(1:NROW(DOE), function(i)
										 data.frame(accuracy(f = arimas[[i]]$fit, 
															 x = arimas[[i]]$actual), 
													row.names = 1))))

# lets add these results to the DOE

DOE$params = params
DOE$t.pval = ttests
DOE = cbind(DOE, stats)

cnames = colnames(DOE)
DOE[,7:13] = sapply(7:13, function(i) DOE[,i] = as.numeric(DOE[,i]))
colnames(DOE) = cnames

DOE

rm(arimas, params, ttests, stats, cnames)

# lets filter out the top models

# lets plot the histograms of these reponses to find filter limits

par(mfrow = c(2, 4))
lapply(7:13, function(i) hist(DOE[,i], main = colnames(DOE)[i]))

DOE

# apply filters

x = DOE[which(DOE$MAPE <= 13 & 
			  DOE$t.pval >= 0.6),]

par(mfrow = c(2, 4))
lapply(7:13, function(i) hist(x[,i], main = colnames(x)[i]))
x

# lets compare the auto fit (scenario 1) with the chosen fit from the DOE (scenario 13)

DOE[c(1, 13),]

# lets build the final models

arimas = lapply(c(1, 13), function(i)
							Arima(ts(DF$CRIMES, frequency = 54), 
									order = c(DOE$p[i], DOE$d[i], DOE$q[i]), 
									seasonal = list(order = c(DOE$P[i], DOE$D[i], DOE$Q[i]), period = 54),
									xreg = DF[,c("SUMMER54", "DAYS54")],						
									optim.control = list(maxit = 1000)))

# lets plot the model & residuals for the auto fit and chosen fit 

	# auto: scenario 1

plot_crimes_weekly_arimas1 = predplot(obj = arimas[[1]],
							  xlab = "Observation Number", 
							  ylab = "Crimes", 
							  main = "Crimes Per Week: Central Philly\nARIMA(0, 1, 1)X(1, 0, 0)54\nRegression Variables: 24.6 * SUMMER54, 3.0 * DAYS54")

plot_crimes_weekly_arimas1

plot_doe_resid1 = residplots(actual = DF$CRIMES[109:565],
							 fit = as.numeric(fitted(arimas[[1]])[109:565]),
							 histlabel.y = -3,
							 n = 54,
							 basefont = 15)

grid.arrange(plot_doe_resid1[[1]], 
			 plot_doe_resid1[[2]], 
			 plot_doe_resid1[[3]], 
			 plot_doe_resid1[[4]], 
			 plot_doe_resid1[[5]], 
			 plot_doe_resid1[[6]], 
			 ncol = 3)		  
		  
	# chosen: scenario 21

plot_crimes_weekly_arimas2 = predplot(obj = arimas[[2]],
							  xlab = "Observation Number", 
							  ylab = "Crimes", 
							  main = "Crimes Per Week: Central Philly\nARIMA(2, 1, 0)X(1, 1, 1)54\nRegression Variables: 9.8 * SUMMER54, -7.3 * DAYS54")

plot_crimes_weekly_arimas2

plot_doe_resid2 = residplots(actual = DF$CRIMES[109:565],
							 fit = as.numeric(fitted(arimas[[2]])[109:565]), 
							 histlabel.y = -2,
							 n = 54, 
							 basefont = 15)

grid.arrange(plot_doe_resid2[[1]], 
			 plot_doe_resid2[[2]], 
			 plot_doe_resid2[[3]], 
			 plot_doe_resid2[[4]], 
			 plot_doe_resid2[[5]], 
			 plot_doe_resid2[[6]], 
			 ncol = 3)		  

# lets store our results into lists

models$DR$weekly$fits$crimes_weekly_CP_lag = arimas[[1]]
models$DR$weekly$plots$plot_crimes_weekly_CP_lag = plot_crimes_weekly_arimas1
models$DR$weekly$plots$plot_crimes_weekly_residual_CP_lag = plot_doe_resid1
models$DR$weekly$tables$model_data_CP_lag = DF
models$DR$weekly$tables$model_doe_CP_lag = DOE		

rm(x, arimas, plot_crimes_weekly_arimas1, plot_crimes_weekly_arimas2, 
   plot_doe_resid1, plot_doe_resid2, DOE, DF)		  

# ---- forecast weekly crimes -------------------------------------------------------------------------------------------------------------------------------------------

# lets set up our data set and external variables

DF = data.frame(tables$crime_Area1_tables$crimes_weekly[Area1 == "CP"])
DF$SUMMER = as.numeric(DF$SUMMER) - 1
DF$DAYS = as.numeric(DF$DAYS)

DF = data.table(DF)
DF[,SUMMER54 := shift(SUMMER, 54)]
DF[,DAYS54 := shift(DAYS, 54)]
setcolorder(DF, c(1:4, 6, 7, 5))

DF = data.frame(DF)
DF = DF[,-c(3:4)]

# lets look at the residuals of a linear model fit using the xreg variables

mod = lm(CRIMES ~ DAYS54 + SUMMER54, data = DF[1:478,])
mod

# lets look at differenced residuals to find values for d

par(mfrow = c(2,3))

plot(resid(mod), main = "d = 0", type = "b")
lapply(1:5, function(i) plot(diff(resid(mod), differences = i), main = paste0("d = ", i), type ="b"))

d = 1

# lets look at the acf and pacf plots of the differenced residuals to find values of p and q

d = 1
grid.arrange(ggacf(x = diff(resid(mod), differences = d), partial = TRUE, n = 54, basefont = 15, main = paste0("PACF Plot of Central Philly Residuals\nDifference of Order ", d)),
			 ggacf(x = diff(resid(mod), differences = d), n = 54, basefont = 15, main = paste0("ACF Plot of Central Philly Residuals\nDifference of Order ", d)),
			 ncol = 2)

p = 0:4
q. = 0:1

# lets look at seasonaly differenced residuals to find values for D

par(mfrow = c(2,3))

plot(resid(mod), main = "d = 0", type = "b")
lapply(1:5, function(i) plot(diff(resid(mod), lag = 54, differences = i), main = paste0("D = ", i), type ="b"))

D = 1

rm(d, p, q., D, mod)

# auto fit

arima.auto = auto.arima(ts(DF$CRIMES[1:478], frequency = 54), xreg = DF[1:478, c("SUMMER54", "DAYS54")])
arima.auto

c(5, 1, 0, 1, 0, 0)

rm(arima.auto)

# lets create our DOE based on the auto.arima, acf, pacf, and differencing plots

DOE = expand.grid("p" = 0:4, 
				  "d" = 1, 
				  "q" = 0:1,
				  "P" = 1,
				  "D" = 1,
				  "Q" = 0:1)

DOE = rbind(c(5, 1, 0, 1, 0, 0), DOE)
DOE = DOE[!duplicated(DOE),]
rownames(DOE) = 1:NROW(DOE)

head(DOE)
NROW(DOE)

# lets compute the arima models from the DOE

arimas = lapply(1:NROW(DOE), function(i)
								Arima.if(ts(DF$CRIMES[1:478], frequency = 54), 
											order = c(DOE$p[i], DOE$d[i], DOE$q[i]), 
											seasonal = list(order = c(DOE$P[i], DOE$D[i], DOE$Q[i]), period = 54),
											xreg = DF[1:478, c("SUMMER54", "DAYS54")],
											optim.control = list(maxit = 1000)))
		
# filter out invalid DOE scenarios

invalid = sapply(1:NROW(DOE), function(i) is.na(arimas[i]))

DOE = DOE[which(invalid == FALSE),]
rownames(DOE) = 1:NROW(DOE)

arimas = arimas[which(invalid == FALSE)]

rm(invalid)		

# lets compute the actual and fitted values for the testing data of each arima model in the DOE

arimas = lapply(1:NROW(DOE), function(i)
								data.frame("actual" = DF$CRIMES[479:565],
											"fit" = as.numeric(forecast(arimas[[i]],
																		h = 87,
																		xreg = DF[479:565, c("SUMMER54", "DAYS54")],
																		level = c(0, 0))$mean)))

# lower values are ideal

params = sapply(1:NROW(DOE), function(i) 
							 DOE$p[i] + DOE$q[i] + DOE$P[i] + DOE$Q[i])

# values greater than 0.05 are ideal 

ttests = sapply(1:NROW(DOE), function(i) 
							 t.test(arimas[[i]]$actual - arimas[[i]]$fit)$p.value)

# values close to zero are ideal 

stats = data.frame(t(sapply(1:NROW(DOE), function(i)
										 data.frame(accuracy(f = arimas[[i]]$fit, 
															 x = arimas[[i]]$actual), 
													row.names = 1))))

# lets add these results to the DOE

DOE$params = params
DOE$t.pval = ttests
DOE = cbind(DOE, stats)

cnames = colnames(DOE)
DOE[,7:13] = sapply(7:13, function(i) DOE[,i] = as.numeric(DOE[,i]))
colnames(DOE) = cnames

DOE

rm(arimas, params, ttests, stats, cnames)

# lets filter out the top models

# lets plot the histograms of these reponses to find filter limits

par(mfrow = c(2, 4))
lapply(7:13, function(i) hist(DOE[,i], main = colnames(DOE)[i]))

DOE

# apply filters

x = DOE[which(DOE$MAPE <= 20),]

par(mfrow = c(2, 4))
lapply(7:13, function(i) hist(x[,i], main = colnames(x)[i]))

x

# lets compare the auto fit (scenario 1) with the chosen fit from the DOE (scenario 15)

DOE[c(1, 15),]

# lets build the final models

arimas = lapply(c(1, 15), function(i)
							Arima(ts(DF$CRIMES[1:478], frequency = 54), 
									order = c(DOE$p[i], DOE$d[i], DOE$q[i]), 
									seasonal = list(order = c(DOE$P[i], DOE$D[i], DOE$Q[i]), period = 54),
									xreg = DF[1:478, c("SUMMER54", "DAYS54")],
									optim.control = list(maxit = 1000)))

# lets plot the model & residuals for the auto fit and chosen fit 

	# auto: scenario 1

plot_crimes_weekly_arimas1 = predplot(obj = arimas[[1]],
							  n.ahead = 87,
							  xreg = DF[479:565, c("SUMMER54", "DAYS54")],
							  actuals = DF$CRIMES,
							  xlab = "Observation Number", 
							  ylab = "Crimes", 
							  main = "Crimes Per Week: Central Philly\nARIMA(5, 1, 0)X(1, 0, 0)54\nRegression Variable: 21.1 * SUMMER54, -5.3 * DAYS54")

plot_crimes_weekly_arimas1

plot_doe_resid1 = residplots(actual = DF$CRIMES[479:565],
							 fit = as.numeric(predict(arimas[[1]], newxreg = DF[479:565, c("SUMMER54", "DAYS54")], n.ahead = 87, se.fit = FALSE)),
							 histlabel.y = -0.5,
							 n = 54,
							 binwidth = 25,
							 basefont = 15)

grid.arrange(plot_doe_resid1[[1]], 
			 plot_doe_resid1[[2]], 
			 plot_doe_resid1[[3]], 
			 plot_doe_resid1[[4]], 
			 plot_doe_resid1[[5]], 
			 plot_doe_resid1[[6]], 
			 ncol = 3)		  
		  
	# chosen: scenario 15

plot_crimes_weekly_arimas2 = predplot(obj = arimas[[2]],
							  n.ahead = 87,
							  xreg = DF[479:565, c("SUMMER54", "DAYS54")],
							  actuals = DF$CRIMES,
							  xlab = "Observation Number", 
							  ylab = "Crimes", 
							  main = "Crimes Per Week: Central Philly\nARIMA(0, 1, 1)X(1, 1, 1)54\nRegression Variable: -1.1 * SUMMER54, -8.0 * DAYS54")

plot_crimes_weekly_arimas2

plot_doe_resid2 = residplots(actual = DF$CRIMES[479:565],
							 fit = as.numeric(predict(arimas[[2]], newxreg = DF[479:565, c("SUMMER54", "DAYS54")], n.ahead = 87, se.fit = FALSE)),
							 histlabel.y = -0.5,
							 n = 54,
							 binwidth = 25,
							 basefont = 15)

grid.arrange(plot_doe_resid2[[1]], 
			 plot_doe_resid2[[2]], 
			 plot_doe_resid2[[3]], 
			 plot_doe_resid2[[4]], 
			 plot_doe_resid2[[5]], 
			 plot_doe_resid2[[6]], 
			 ncol = 3)		  

# lets store our results into lists

models$DR$weekly$fits$crimes_weekly_CP_lag_forecast = arimas[[2]]
models$DR$weekly$plots$plot_crimes_weekly_CP_lag_forecast = plot_crimes_weekly_arimas2
models$DR$weekly$plots$plot_crimes_weekly_residual_CP_lag_forecast = plot_doe_resid2
models$DR$weekly$tables$model_data_CP_lag_forecast = DF
models$DR$weekly$tables$model_doe_CP_lag_forecast = DOE		  

rm(x, arimas, plot_crimes_weekly_arimas1, plot_crimes_weekly_arimas2, 
   plot_doe_resid1, plot_doe_resid2, DOE, DF)		  

}

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ---- Dynamic Regression Forecasting: West Philly -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

{

# ---- fit weekly crimes -------------------------------------------------------------------------------------------------------------------------------------------

# lets set up our data set and external variables

DF = data.frame(tables$crime_Area1_tables$crimes_weekly[Area1 == "WP"])

# lets look at plots of the dataset

mycolors = colorRampPalette(brewer.pal(n = 8, name = "Dark2"))(9)

p1 = ggplot(data = DF, aes(x = DAYS, y = CRIMES, fill = DAYS, color = DAYS)) +
	 geom_violin(alpha = 1/4, width = 5/4) +
	 geom_jitter(alpha = 1/2, width = 1/6) +
	 scale_color_manual(values = mycolors[1:7]) + 
	 scale_fill_manual(values = mycolors[1:7]) +
	 labs(x = "Days", y = "Crimes") +
	 ggtitle("West Philly\nCrimes Per Week v. Days") +
	 theme_light(base_size = 15) +
	 theme(legend.position = "none")

p2 = ggplot(data = DF, aes(x = SUMMER, y = CRIMES, fill = SUMMER, color = SUMMER)) +
	 geom_violin(alpha = 1/4, width = 3/4) +
	 geom_jitter(alpha = 1/2, width = 1/5) +
	 scale_color_manual(values = mycolors[8:9]) + 
	 scale_fill_manual(values = mycolors[8:9]) +
	 labs(x = "Summer", y = "Crimes") +
	 ggtitle("West Philly\nCrimes Per Week v. Summer") +
	 theme_light(base_size = 15) +
	 theme(legend.position = "none")

grid.arrange(p1, p2, ncol = 2)

DF$SUMMER = as.numeric(DF$SUMMER) - 1
DF$DAYS = as.numeric(DF$DAYS)

rm(p1, p2, mycolors)

# lets look at the residuals of a linear model fit using the xreg variables

mod = lm(CRIMES ~ DAYS + SUMMER, data = DF)
mod

# lets look at differenced residuals to find values for d

par(mfrow = c(2,3))

plot(resid(mod), main = "d = 0", type = "b")
lapply(1:5, function(i) plot(diff(resid(mod), differences = i), main = paste0("d = ", i), type ="b"))

d = 1

# lets look at the acf and pacf plots of the differenced residuals to find values of p and q

d = 1
grid.arrange(ggacf(x = diff(resid(mod), differences = d), partial = TRUE, n = 54, basefont = 15, main = paste0("PACF Plot of West Philly Residuals\nDifference of Order ", d)),
			 ggacf(x = diff(resid(mod), differences = d), n = 54, basefont = 15, main = paste0("ACF Plot of West Philly Residuals\nDifference of Order ", d)),
			 ncol = 2)

p = 0:2
q. = 0:1

# lets look at seasonaly differenced residuals to find values for D

par(mfrow = c(2,3))

plot(resid(mod), main = "d = 0", type = "b")
lapply(1:5, function(i) plot(diff(resid(mod), lag = 54, differences = i), main = paste0("D = ", i), type ="b"))

D = 1

rm(d, p, q., D, mod)

# compare auto fits

arima.auto1 = auto.arima(ts(DF$CRIMES, frequency = 54), xreg = DF[,c("SUMMER", "DAYS")])
arima.auto1

c(0, 1, 1, 0, 0, 2)

arima.auto2 = auto.arima(ts(DF$CRIMES, frequency = 54), xreg = DF[,c("SUMMER")])
arima.auto2

c(1, 1, 1, 0, 0, 2)

arima.auto3 = auto.arima(ts(DF$CRIMES, frequency = 54), xreg = DF[,c("DAYS")])
arima.auto3

c(1, 1, 1, 0, 0, 2)

autofits = data.frame("Variables" = rbind("SUMMER,DAYS", "SUMMER", "DAYS"),
					  rbind(accuracy(arima.auto1), accuracy(arima.auto2), accuracy(arima.auto3)),
					  "AIC" = rbind(arima.auto1$aic, arima.auto2$aic, arima.auto3$aic),
					  "AICc" = rbind(arima.auto1$aicc, arima.auto2$aicc, arima.auto3$aicc),
					  "BIC" = rbind(arima.auto1$bic, arima.auto2$bic, arima.auto3$bic))

autofits = autofits[,-c(7, 8)]
autofits

rm(arima.auto1, arima.auto2, arima.auto3)

# lets create our DOE based on the auto fits, acf, pacf, and differencing plots

DOE = expand.grid("p" = 0:2, 
				  "d" = 1, 
				  "q" = 0:1,
				  "P" = 0:1,
				  "D" = 1,
				  "Q" = 0:1)

DOE = rbind(c(0, 1, 1, 0, 0, 2), DOE)
DOE = DOE[!duplicated(DOE),]
rownames(DOE) = 1:NROW(DOE)

head(DOE)
NROW(DOE)

# lets compute the actual and fitted values for each arima model in the DOE

arimas = lapply(1:NROW(DOE), function(i)
								 data.frame("actual" = DF$CRIMES[55:565],
											"fit" = as.numeric(fitted.if(Arima.if(ts(DF$CRIMES, frequency = 54), 
																			order = c(DOE$p[i], DOE$d[i], DOE$q[i]), 
																			seasonal = list(order = c(DOE$P[i], DOE$D[i], DOE$Q[i]), period = 54),
																			xreg = DF[,c("SUMMER", "DAYS")],
																			optim.control = list(maxit = 1000)))[55:565])))

# filter out invalid DOE scenarios

invalid = sapply(1:NROW(DOE), function(i) all(is.na(arimas[[i]]$fit)))

DOE = DOE[which(invalid == FALSE),]
rownames(DOE) = 1:NROW(DOE)

arimas = arimas[which(invalid == FALSE)]

rm(invalid)
	
# lower values are ideal

params = sapply(1:NROW(DOE), function(i) 
							 DOE$p[i] + DOE$q[i] + DOE$P[i] + DOE$Q[i])

# values greater than 0.05 are ideal 

ttests = sapply(1:NROW(DOE), function(i) 
							 t.test(arimas[[i]]$actual - arimas[[i]]$fit)$p.value)

# values close to zero are ideal 

stats = data.frame(t(sapply(1:NROW(DOE), function(i)
										 data.frame(accuracy(f = arimas[[i]]$fit, 
															 x = arimas[[i]]$actual), 
													row.names = 1))))

# lets add these results to the DOE

DOE$params = params
DOE$t.pval = ttests
DOE = cbind(DOE, stats)

cnames = colnames(DOE)
DOE[,7:13] = sapply(7:13, function(i) DOE[,i] = as.numeric(DOE[,i]))
colnames(DOE) = cnames

DOE

rm(arimas, params, ttests, stats, cnames)

# lets filter out the top models

# lets plot the histograms of these reponses to find filter limits

par(mfrow = c(2, 4))
lapply(7:13, function(i) hist(DOE[,i], main = colnames(DOE)[i]))

DOE

# apply filters

x = DOE[which(DOE$MAPE <= 8),]

par(mfrow = c(2, 4))
lapply(7:13, function(i) hist(x[,i], main = colnames(x)[i]))
x

# lets compare the auto fit (scenario 1) with the chosen fit from the DOE (scenario 16)

DOE[c(1, 16),]

# lets build the final models

arimas = lapply(c(1, 16), function(i)
							Arima(ts(DF$CRIMES, frequency = 54), 
									order = c(DOE$p[i], DOE$d[i], DOE$q[i]), 
									seasonal = list(order = c(DOE$P[i], DOE$D[i], DOE$Q[i]), period = 54),
									xreg = DF[,c("SUMMER", "DAYS")],						
									optim.control = list(maxit = 1000)))

# lets plot the model & residuals for the auto fit and chosen fit 

	# auto: scenario 1

plot_crimes_weekly_arimas1 = predplot(obj = arimas[[1]],
							  xlab = "Observation Number", 
							  ylab = "Crimes", 
							  main = "Crimes Per Week: West Philly\nARIMA(0, 1, 1)X(0, 0, 2)54\nRegression Variables: 21.3 * SUMMER, 81.8 * DAYS")

plot_crimes_weekly_arimas1

plot_doe_resid1 = residplots(actual = DF$CRIMES[55:565],
							 fit = as.numeric(fitted(arimas[[1]])[55:565]),
							 histlabel.y = -2,
							 n = 54,
							 basefont = 15)

grid.arrange(plot_doe_resid1[[1]], 
			 plot_doe_resid1[[2]], 
			 plot_doe_resid1[[3]], 
			 plot_doe_resid1[[4]], 
			 plot_doe_resid1[[5]], 
			 plot_doe_resid1[[6]], 
			 ncol = 3)		  
		  
	# chosen: scenario 16

plot_crimes_weekly_arimas2 = predplot(obj = arimas[[2]],
							  xlab = "Observation Number", 
							  ylab = "Crimes", 
							  main = "Crimes Per Week: West Philly\nARIMA(0, 1, 1)X(0, 1, 1)54\nRegression Variables: 3.0 * SUMMER, 80.9 * DAYS")

plot_crimes_weekly_arimas2

plot_doe_resid2 = residplots(actual = DF$CRIMES[55:565],
							 fit = as.numeric(fitted(arimas[[2]])[55:565]), 
							 histlabel.y = -2,
							 n = 54, 
							 basefont = 15)

grid.arrange(plot_doe_resid2[[1]], 
			 plot_doe_resid2[[2]], 
			 plot_doe_resid2[[3]], 
			 plot_doe_resid2[[4]], 
			 plot_doe_resid2[[5]], 
			 plot_doe_resid2[[6]], 
			 ncol = 3)		  

# lets store our results into lists

models$DR$weekly$fits$crimes_weekly_WP = arimas[[2]]
models$DR$weekly$plots$plot_crimes_weekly_WP = plot_crimes_weekly_arimas2
models$DR$weekly$plots$plot_crimes_weekly_residual_WP = plot_doe_resid2
models$DR$weekly$tables$model_data_WP = DF
models$DR$weekly$tables$model_doe_WP = DOE		

rm(x, arimas, plot_crimes_weekly_arimas1, plot_crimes_weekly_arimas2, 
   plot_doe_resid1, plot_doe_resid2, DOE, DF)		  

models$DR$weekly$tables$autofits_WP = autofits

rm(autofits)

# ---- forecast weekly crimes -------------------------------------------------------------------------------------------------------------------------------------------

# lets set up our data set and external variables

DF = data.frame(tables$crime_Area1_tables$crimes_weekly[Area1 == "WP"])
DF$SUMMER = as.numeric(DF$SUMMER) - 1
DF$DAYS = as.numeric(DF$DAYS)

# lets look at the residuals of a linear model fit using the xreg variables

mod = lm(CRIMES ~ DAYS + SUMMER, data = DF[1:478,])
mod

# lets look at differenced residuals to find values for d

par(mfrow = c(2,3))

plot(resid(mod), main = "d = 0", type = "b")
lapply(1:5, function(i) plot(diff(resid(mod), differences = i), main = paste0("d = ", i), type ="b"))

d = 1

# lets look at the acf and pacf plots of the differenced residuals to find values of p and q

d = 1
grid.arrange(ggacf(x = diff(resid(mod), differences = d), partial = TRUE, n = 54, basefont = 15, main = paste0("PACF Plot of West Philly Residuals\nDifference of Order ", d)),
			 ggacf(x = diff(resid(mod), differences = d), n = 54, basefont = 15, main = paste0("ACF Plot of West Philly Residuals\nDifference of Order ", d)),
			 ncol = 2)

p = 0:2
q. = 0:1

# lets look at seasonaly differenced residuals to find values for D

par(mfrow = c(2,3))

plot(resid(mod), main = "d = 0", type = "b")
lapply(1:5, function(i) plot(diff(resid(mod), lag = 54, differences = i), main = paste0("D = ", i), type ="b"))

D = 1

rm(d, p, q., D, mod)

# auto fit

arima.auto = auto.arima(ts(DF$CRIMES[1:478], frequency = 54), xreg = DF[1:478, c("SUMMER", "DAYS")])
arima.auto

c(0, 1, 1, 0, 0, 0)

rm(arima.auto)

# lets create our DOE based on the auto.arima, acf, pacf, and differencing plots

DOE = expand.grid("p" = 0:2, 
				  "d" = 1, 
				  "q" = 0:1,
				  "P" = 0:1,
				  "D" = 1,
				  "Q" = 0:1)

DOE = rbind(c(0, 1, 1, 0, 0, 0), DOE)
DOE = DOE[!duplicated(DOE),]
rownames(DOE) = 1:NROW(DOE)

head(DOE)
NROW(DOE)

# lets compute the arima models from the DOE

arimas = lapply(1:NROW(DOE), function(i)
								Arima.if(ts(DF$CRIMES[1:478], frequency = 54), 
											order = c(DOE$p[i], DOE$d[i], DOE$q[i]), 
											seasonal = list(order = c(DOE$P[i], DOE$D[i], DOE$Q[i]), period = 54),
											xreg = DF[1:478, c("SUMMER", "DAYS")],
											optim.control = list(maxit = 1000)))
		
# filter out invalid DOE scenarios

invalid = sapply(1:NROW(DOE), function(i) is.na(arimas[i]))

DOE = DOE[which(invalid == FALSE),]
rownames(DOE) = 1:NROW(DOE)

arimas = arimas[which(invalid == FALSE)]

rm(invalid)		

# lets compute the actual and fitted values for the testing data of each arima model in the DOE

arimas = lapply(1:NROW(DOE), function(i)
								data.frame("actual" = DF$CRIMES[479:565],
											"fit" = as.numeric(forecast(arimas[[i]],
																		h = 87,
																		xreg = DF[479:565, c("SUMMER", "DAYS")],
																		level = c(0, 0))$mean)))

# lower values are ideal

params = sapply(1:NROW(DOE), function(i) 
							 DOE$p[i] + DOE$q[i] + DOE$P[i] + DOE$Q[i])

# values greater than 0.05 are ideal 

ttests = sapply(1:NROW(DOE), function(i) 
							 t.test(arimas[[i]]$actual - arimas[[i]]$fit)$p.value)

# values close to zero are ideal 

stats = data.frame(t(sapply(1:NROW(DOE), function(i)
										 data.frame(accuracy(f = arimas[[i]]$fit, 
															 x = arimas[[i]]$actual), 
													row.names = 1))))

# lets add these results to the DOE

DOE$params = params
DOE$t.pval = ttests
DOE = cbind(DOE, stats)

cnames = colnames(DOE)
DOE[,7:13] = sapply(7:13, function(i) DOE[,i] = as.numeric(DOE[,i]))
colnames(DOE) = cnames

DOE

rm(arimas, params, ttests, stats, cnames)

# lets filter out the top models

# lets plot the histograms of these reponses to find filter limits

par(mfrow = c(2, 4))
lapply(7:13, function(i) hist(DOE[,i], main = colnames(DOE)[i]))

DOE

# apply filters

x = DOE[which(DOE$MAPE <= 11 &
			  DOE$t.pval > 0.2),]

par(mfrow = c(2, 4))
lapply(7:13, function(i) hist(x[,i], main = colnames(x)[i]))

x

# lets compare the auto fit (scenario 1) with the chosen fit from the DOE (scenario 16)

DOE[c(1, 16),]

# lets build the final models

arimas = lapply(c(1, 16), function(i)
							Arima(ts(DF$CRIMES[1:478], frequency = 54), 
									order = c(DOE$p[i], DOE$d[i], DOE$q[i]), 
									seasonal = list(order = c(DOE$P[i], DOE$D[i], DOE$Q[i]), period = 54),
									xreg = DF[1:478, c("SUMMER", "DAYS")],
									optim.control = list(maxit = 1000)))

# lets plot the model & residuals for the auto fit and chosen fit 

	# auto: scenario 1

plot_crimes_weekly_arimas1 = predplot(obj = arimas[[1]],
							  n.ahead = 87,
							  xreg = DF[479:565, c("SUMMER", "DAYS")],
							  actuals = DF$CRIMES,
							  xlab = "Observation Number", 
							  ylab = "Crimes", 
							  main = "Crimes Per Week: West Philly\nARIMA(0, 1, 1)X(0, 0, 0)54\nRegression Variable: 16.7 * SUMMER, 83.6 * DAYS")

plot_crimes_weekly_arimas1

plot_doe_resid1 = residplots(actual = DF$CRIMES[479:565],
							 fit = as.numeric(predict(arimas[[1]], newxreg = DF[479:565, c("SUMMER", "DAYS")], n.ahead = 87, se.fit = FALSE)),
							 histlabel.y = -0.5,
							 n = 54,
							 binwidth = 25,
							 basefont = 15)

grid.arrange(plot_doe_resid1[[1]], 
			 plot_doe_resid1[[2]], 
			 plot_doe_resid1[[3]], 
			 plot_doe_resid1[[4]], 
			 plot_doe_resid1[[5]], 
			 plot_doe_resid1[[6]], 
			 ncol = 3)		  
		  
	# chosen: scenario 16

plot_crimes_weekly_arimas2 = predplot(obj = arimas[[2]],
							  n.ahead = 87,
							  xreg = DF[479:565, c("SUMMER", "DAYS")],
							  actuals = DF$CRIMES,
							  xlab = "Observation Number", 
							  ylab = "Crimes", 
							  main = "Crimes Per Week: West Philly\nARIMA(0, 1, 1)X(0, 1, 1)54\nRegression Variable: -4.9 * SUMMER, 81.4 * DAYS")

plot_crimes_weekly_arimas2

plot_doe_resid2 = residplots(actual = DF$CRIMES[479:565],
							 fit = as.numeric(predict(arimas[[2]], newxreg = DF[479:565, c("SUMMER", "DAYS")], n.ahead = 87, se.fit = FALSE)),
							 histlabel.y = -0.5,
							 n = 54,
							 binwidth = 30,
							 basefont = 15)

grid.arrange(plot_doe_resid2[[1]], 
			 plot_doe_resid2[[2]], 
			 plot_doe_resid2[[3]], 
			 plot_doe_resid2[[4]], 
			 plot_doe_resid2[[5]], 
			 plot_doe_resid2[[6]], 
			 ncol = 3)		  

# lets store our results into lists

models$DR$weekly$fits$crimes_weekly_WP_forecast = arimas[[2]]
models$DR$weekly$plots$plot_crimes_weekly_WP_forecast = plot_crimes_weekly_arimas2
models$DR$weekly$plots$plot_crimes_weekly_residual_WP_forecast = plot_doe_resid2
models$DR$weekly$tables$model_data_WP_forecast = DF
models$DR$weekly$tables$model_doe_WP_forecast = DOE		  

rm(x, arimas, plot_crimes_weekly_arimas1, plot_crimes_weekly_arimas2, 
   plot_doe_resid1, plot_doe_resid2, DOE, DF)		  

}

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ---- Dynamic Regression Forecasting with Lagged Regression Variables: West Philly -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

{

# ---- fit weekly crimes -------------------------------------------------------------------------------------------------------------------------------------------

# lets set up our data set and external variables
# we will compare three different lags for each of the two regression variables
	# a lag of 1, 4, and 54 weeks
	
DF = data.frame(tables$crime_Area1_tables$crimes_weekly[Area1 == "WP"])
DF$SUMMER = as.numeric(DF$SUMMER) - 1
DF$DAYS = as.numeric(DF$DAYS)

DF = data.table(DF)

DF[,SUMMER1 := shift(SUMMER, 1)]
DF[,SUMMER4 := shift(SUMMER, 4)]
DF[,SUMMER54 := shift(SUMMER, 54)]

DF[,DAYS1 := shift(DAYS, 1)]
DF[,DAYS4 := shift(DAYS, 4)]
DF[,DAYS54 := shift(DAYS, 54)]

setcolorder(DF, c(1, 2, 3, 6, 7, 8, 4, 9, 10, 11, 5))
DF = data.frame(DF)

head(DF)

# lets create a DOE to test all possible regression lag combinations
# we will choose the best combination based on the R^2-Pred, AIC, and BIC values of linear regression models

DOE = expand.grid("SUMMER" = 4:6,
				  "DAYS" = 8:10)

DOE$r2pred = sapply(1:NROW(DOE), function(i)
					statslm(lm(DF$CRIMES[55:565] ~ DF[55:565,DOE$DAYS[i]] + DF[55:565,DOE$SUMMER[i]]))$Prediction)

DOE = cbind(DOE, t(sapply(1:NROW(DOE), function(i)
							statslm(lm(DF$CRIMES[55:565] ~ DF[55:565,DOE$DAYS[i]] + DF[55:565,DOE$SUMMER[i]]))$Fitness)))

DOE$SUMMER = colnames(DF)[DOE$SUMMER]
DOE$DAYS = colnames(DF)[DOE$DAYS]

DOE

# update the dataset with the best lagged regression variables

DF = DF[,c(1,2,4,10,11)]
head(DF)

rm(DOE)

# lets look at the residuals of a linear model fit using the xreg variables

mod = lm(CRIMES ~ DAYS54 + SUMMER1, data = DF)
mod

# lets look at differenced residuals to find values for d

par(mfrow = c(2,3))

plot(resid(mod), main = "d = 0", type = "b")
lapply(1:5, function(i) plot(diff(resid(mod), differences = i), main = paste0("d = ", i), type ="b"))

d = 1

# lets look at the acf and pacf plots of the differenced residuals to find values of p and q

d = 1
grid.arrange(ggacf(x = diff(resid(mod), differences = d), partial = TRUE, n = 54, basefont = 15, main = paste0("PACF Plot of West Philly Residuals\nDifference of Order ", d)),
			 ggacf(x = diff(resid(mod), differences = d), n = 54, basefont = 15, main = paste0("ACF Plot of West Philly Residuals\nDifference of Order ", d)),
			 ncol = 2)

p = 0:4
q. = 0:1

# lets look at seasonaly differenced residuals to find values for D

par(mfrow = c(2,3))

plot(resid(mod), main = "d = 0", type = "b")
lapply(1:5, function(i) plot(diff(resid(mod), lag = 54, differences = i), main = paste0("D = ", i), type ="b"))

D = 1

rm(d, p, q., D, mod)

# lets look at the auto fit

arima.auto = auto.arima(ts(DF$CRIMES, frequency = 54), xreg = DF[,c("SUMMER1", "DAYS54")])
arima.auto

c(5, 1, 0, 1, 0, 0)

rm(arima.auto)

# lets create our DOE based on the auto fit, acf, pacf, and differencing plots

DOE = expand.grid("p" = 0:4, 
				  "d" = 1, 
				  "q" = 0:1,
				  "P" = 1,
				  "D" = 1,
				  "Q" = 0:1)

DOE = rbind(c(5, 1, 0, 1, 0, 0), DOE)
DOE = DOE[!duplicated(DOE),]
rownames(DOE) = 1:NROW(DOE)

head(DOE)
NROW(DOE)

# lets compute the actual and fitted values for each arima model in the DOE

arimas = lapply(1:NROW(DOE), function(i)
								 data.frame("actual" = DF$CRIMES[109:565],
											"fit" = as.numeric(fitted.if(Arima.if(ts(DF$CRIMES, frequency = 54), 
																			order = c(DOE$p[i], DOE$d[i], DOE$q[i]), 
																			seasonal = list(order = c(DOE$P[i], DOE$D[i], DOE$Q[i]), period = 54),
																			xreg = DF[,c("SUMMER1", "DAYS54")],
																			optim.control = list(maxit = 1000)))[109:565])))
	
# filter out invalid DOE scenarios

invalid = sapply(1:NROW(DOE), function(i) all(is.na(arimas[[i]]$fit)))

DOE = DOE[which(invalid == FALSE),]
rownames(DOE) = 1:NROW(DOE)

arimas = arimas[which(invalid == FALSE)]

rm(invalid)
	
# lower values are ideal

params = sapply(1:NROW(DOE), function(i) 
							 DOE$p[i] + DOE$q[i] + DOE$P[i] + DOE$Q[i])

# values greater than 0.05 are ideal 

ttests = sapply(1:NROW(DOE), function(i) 
							 t.test(arimas[[i]]$actual - arimas[[i]]$fit)$p.value)

# values close to zero are ideal 

stats = data.frame(t(sapply(1:NROW(DOE), function(i)
										 data.frame(accuracy(f = arimas[[i]]$fit, 
															 x = arimas[[i]]$actual), 
													row.names = 1))))

# lets add these results to the DOE

DOE$params = params
DOE$t.pval = ttests
DOE = cbind(DOE, stats)

cnames = colnames(DOE)
DOE[,7:13] = sapply(7:13, function(i) DOE[,i] = as.numeric(DOE[,i]))
colnames(DOE) = cnames

DOE

rm(arimas, params, ttests, stats, cnames)

# lets filter out the top models

# lets plot the histograms of these reponses to find filter limits

par(mfrow = c(2, 4))
lapply(7:13, function(i) hist(DOE[,i], main = colnames(DOE)[i]))

DOE

# apply filters

x = DOE[which(DOE$MAPE <= 13 & 
			  DOE$t.pval >= 0.6),]

par(mfrow = c(2, 4))
lapply(7:13, function(i) hist(x[,i], main = colnames(x)[i]))
x

# lets compare the auto fit (scenario 1) with the chosen fit from the DOE (scenario 13)

DOE[c(1, 13),]

# lets build the final models

arimas = lapply(c(1, 13), function(i)
							Arima(ts(DF$CRIMES, frequency = 54), 
									order = c(DOE$p[i], DOE$d[i], DOE$q[i]), 
									seasonal = list(order = c(DOE$P[i], DOE$D[i], DOE$Q[i]), period = 54),
									xreg = DF[,c("SUMMER1", "DAYS54")],						
									optim.control = list(maxit = 1000)))

# lets plot the model & residuals for the auto fit and chosen fit 

	# auto: scenario 1

plot_crimes_weekly_arimas1 = predplot(obj = arimas[[1]],
							  xlab = "Observation Number", 
							  ylab = "Crimes", 
							  main = "Crimes Per Week: West Philly\nARIMA(5, 1, 0)X(1, 0, 0)54\nRegression Variables: -0.4 * SUMMER1, 8.4 * DAYS54")

plot_crimes_weekly_arimas1

plot_doe_resid1 = residplots(actual = DF$CRIMES[109:565],
							 fit = as.numeric(fitted(arimas[[1]])[109:565]),
							 histlabel.y = -3,
							 n = 54,
							 basefont = 15)

grid.arrange(plot_doe_resid1[[1]], 
			 plot_doe_resid1[[2]], 
			 plot_doe_resid1[[3]], 
			 plot_doe_resid1[[4]], 
			 plot_doe_resid1[[5]], 
			 plot_doe_resid1[[6]], 
			 ncol = 3)		  
		  
	# chosen: scenario 13

plot_crimes_weekly_arimas2 = predplot(obj = arimas[[2]],
							  xlab = "Observation Number", 
							  ylab = "Crimes", 
							  main = "Crimes Per Week: West Philly\nARIMA(3, 1, 0)X(1, 1, 1)54\nRegression Variables: -21.8 * SUMMER1, 5.3 * DAYS54")

plot_crimes_weekly_arimas2

plot_doe_resid2 = residplots(actual = DF$CRIMES[109:565],
							 fit = as.numeric(fitted(arimas[[2]])[109:565]), 
							 histlabel.y = -2,
							 n = 54, 
							 basefont = 15)

grid.arrange(plot_doe_resid2[[1]], 
			 plot_doe_resid2[[2]], 
			 plot_doe_resid2[[3]], 
			 plot_doe_resid2[[4]], 
			 plot_doe_resid2[[5]], 
			 plot_doe_resid2[[6]], 
			 ncol = 3)		  

# lets store our results into lists

models$DR$weekly$fits$crimes_weekly_WP_lag = arimas[[2]]
models$DR$weekly$plots$plot_crimes_weekly_WP_lag = plot_crimes_weekly_arimas2
models$DR$weekly$plots$plot_crimes_weekly_residual_WP_lag = plot_doe_resid2
models$DR$weekly$tables$model_data_WP_lag = DF
models$DR$weekly$tables$model_doe_WP_lag = DOE		

rm(x, arimas, plot_crimes_weekly_arimas1, plot_crimes_weekly_arimas2, 
   plot_doe_resid1, plot_doe_resid2, DOE, DF)		  

# ---- forecast weekly crimes -------------------------------------------------------------------------------------------------------------------------------------------

# lets set up our data set and external variables

DF = data.frame(tables$crime_Area1_tables$crimes_weekly[Area1 == "WP"])
DF$SUMMER = as.numeric(DF$SUMMER) - 1
DF$DAYS = as.numeric(DF$DAYS)

DF = data.table(DF)
DF[,SUMMER1 := shift(SUMMER, 1)]
DF[,DAYS54 := shift(DAYS, 54)]
setcolorder(DF, c(1:4, 6, 7, 5))

DF = data.frame(DF)
DF = DF[,-c(3:4)]

# lets look at the residuals of a linear model fit using the xreg variables

mod = lm(CRIMES ~ DAYS54 + SUMMER1, data = DF[1:478,])
mod

# lets look at differenced residuals to find values for d

par(mfrow = c(2,3))

plot(resid(mod), main = "d = 0", type = "b")
lapply(1:5, function(i) plot(diff(resid(mod), differences = i), main = paste0("d = ", i), type ="b"))

d = 1

# lets look at the acf and pacf plots of the differenced residuals to find values of p and q

d = 1
grid.arrange(ggacf(x = diff(resid(mod), differences = d), partial = TRUE, n = 54, basefont = 15, main = paste0("PACF Plot of West Philly Residuals\nDifference of Order ", d)),
			 ggacf(x = diff(resid(mod), differences = d), n = 54, basefont = 15, main = paste0("ACF Plot of West Philly Residuals\nDifference of Order ", d)),
			 ncol = 2)

p = 0:3
q. = 0:1

# lets look at seasonaly differenced residuals to find values for D

par(mfrow = c(2,3))

plot(resid(mod), main = "d = 0", type = "b")
lapply(1:5, function(i) plot(diff(resid(mod), lag = 54, differences = i), main = paste0("D = ", i), type ="b"))

D = 1

rm(d, p, q., D, mod)

# auto fit

arima.auto = auto.arima(ts(DF$CRIMES[1:478], frequency = 54), xreg = DF[1:478, c("SUMMER1", "DAYS54")])
arima.auto

c(5, 1, 0, 1, 0, 0)

rm(arima.auto)

# lets create our DOE based on the auto.arima, acf, pacf, and differencing plots

DOE = expand.grid("p" = 0:3, 
				  "d" = 1, 
				  "q" = 0:1,
				  "P" = 1,
				  "D" = 1,
				  "Q" = 0:1)

DOE = rbind(c(5, 1, 0, 1, 0, 0), DOE)
DOE = DOE[!duplicated(DOE),]
rownames(DOE) = 1:NROW(DOE)

head(DOE)
NROW(DOE)

# lets compute the arima models from the DOE

arimas = lapply(1:NROW(DOE), function(i)
								Arima.if(ts(DF$CRIMES[1:478], frequency = 54), 
											order = c(DOE$p[i], DOE$d[i], DOE$q[i]), 
											seasonal = list(order = c(DOE$P[i], DOE$D[i], DOE$Q[i]), period = 54),
											xreg = DF[1:478, c("SUMMER1", "DAYS54")],
											optim.control = list(maxit = 1000)))
		
# filter out invalid DOE scenarios

invalid = sapply(1:NROW(DOE), function(i) is.na(arimas[i]))

DOE = DOE[which(invalid == FALSE),]
rownames(DOE) = 1:NROW(DOE)

arimas = arimas[which(invalid == FALSE)]

rm(invalid)		

# lets compute the actual and fitted values for the testing data of each arima model in the DOE

arimas = lapply(1:NROW(DOE), function(i)
								data.frame("actual" = DF$CRIMES[479:565],
											"fit" = as.numeric(forecast(arimas[[i]],
																		h = 87,
																		xreg = DF[479:565, c("SUMMER1", "DAYS54")],
																		level = c(0, 0))$mean)))

# lower values are ideal

params = sapply(1:NROW(DOE), function(i) 
							 DOE$p[i] + DOE$q[i] + DOE$P[i] + DOE$Q[i])

# values greater than 0.05 are ideal 

ttests = sapply(1:NROW(DOE), function(i) 
							 t.test(arimas[[i]]$actual - arimas[[i]]$fit)$p.value)

# values close to zero are ideal 

stats = data.frame(t(sapply(1:NROW(DOE), function(i)
										 data.frame(accuracy(f = arimas[[i]]$fit, 
															 x = arimas[[i]]$actual), 
													row.names = 1))))

# lets add these results to the DOE

DOE$params = params
DOE$t.pval = ttests
DOE = cbind(DOE, stats)

cnames = colnames(DOE)
DOE[,7:13] = sapply(7:13, function(i) DOE[,i] = as.numeric(DOE[,i]))
colnames(DOE) = cnames

DOE

rm(arimas, params, ttests, stats, cnames)

# lets filter out the top models

# lets plot the histograms of these reponses to find filter limits

par(mfrow = c(2, 4))
lapply(7:13, function(i) hist(DOE[,i], main = colnames(DOE)[i]))

DOE

# apply filters

x = DOE[which(DOE$MAPE <= 20),]

par(mfrow = c(2, 4))
lapply(7:13, function(i) hist(x[,i], main = colnames(x)[i]))

x

# lets compare the auto fit (scenario 1) with the chosen fit from the DOE (scenario 13)

DOE[c(1, 13),]

# lets build the final models

arimas = lapply(c(1, 13), function(i)
							Arima(ts(DF$CRIMES[1:478], frequency = 54), 
									order = c(DOE$p[i], DOE$d[i], DOE$q[i]), 
									seasonal = list(order = c(DOE$P[i], DOE$D[i], DOE$Q[i]), period = 54),
									xreg = DF[1:478, c("SUMMER1", "DAYS54")],
									optim.control = list(maxit = 1000)))

# lets plot the model & residuals for the auto fit and chosen fit 

	# auto: scenario 1

plot_crimes_weekly_arimas1 = predplot(obj = arimas[[1]],
							  n.ahead = 87,
							  xreg = DF[479:565, c("SUMMER1", "DAYS54")],
							  actuals = DF$CRIMES,
							  xlab = "Observation Number", 
							  ylab = "Crimes", 
							  main = "Crimes Per Week: West Philly\nARIMA(5, 1, 0)X(1, 0, 0)54\nRegression Variable: 3.9 * SUMMER1, 3.8 * DAYS54")

plot_crimes_weekly_arimas1

plot_doe_resid1 = residplots(actual = DF$CRIMES[479:565],
							 fit = as.numeric(predict(arimas[[1]], newxreg = DF[479:565, c("SUMMER1", "DAYS54")], n.ahead = 87, se.fit = FALSE)),
							 histlabel.y = -0.5,
							 n = 54,
							 binwidth = 25,
							 basefont = 15)

grid.arrange(plot_doe_resid1[[1]], 
			 plot_doe_resid1[[2]], 
			 plot_doe_resid1[[3]], 
			 plot_doe_resid1[[4]], 
			 plot_doe_resid1[[5]], 
			 plot_doe_resid1[[6]], 
			 ncol = 3)		  
		  
	# chosen: scenario 13

plot_crimes_weekly_arimas2 = predplot(obj = arimas[[2]],
							  n.ahead = 87,
							  xreg = DF[479:565, c("SUMMER1", "DAYS54")],
							  actuals = DF$CRIMES,
							  xlab = "Observation Number", 
							  ylab = "Crimes", 
							  main = "Crimes Per Week: West Philly\nARIMA(0, 1, 1)X(1, 1, 1)54\nRegression Variable: -18.7 * SUMMER1, -1.4 * DAYS54")

plot_crimes_weekly_arimas2

plot_doe_resid2 = residplots(actual = DF$CRIMES[479:565],
							 fit = as.numeric(predict(arimas[[2]], newxreg = DF[479:565, c("SUMMER1", "DAYS54")], n.ahead = 87, se.fit = FALSE)),
							 histlabel.y = -0.5,
							 n = 54,
							 binwidth = 25,
							 basefont = 15)

grid.arrange(plot_doe_resid2[[1]], 
			 plot_doe_resid2[[2]], 
			 plot_doe_resid2[[3]], 
			 plot_doe_resid2[[4]], 
			 plot_doe_resid2[[5]], 
			 plot_doe_resid2[[6]], 
			 ncol = 3)		  

# lets store our results into lists

models$DR$weekly$fits$crimes_weekly_WP_lag_forecast = arimas[[2]]
models$DR$weekly$plots$plot_crimes_weekly_WP_lag_forecast = plot_crimes_weekly_arimas2
models$DR$weekly$plots$plot_crimes_weekly_residual_WP_lag_forecast = plot_doe_resid2
models$DR$weekly$tables$model_data_WP_lag_forecast = DF
models$DR$weekly$tables$model_doe_WP_lag_forecast = DOE		  

rm(x, arimas, plot_crimes_weekly_arimas1, plot_crimes_weekly_arimas2, 
   plot_doe_resid1, plot_doe_resid2, DOE, DF)		  

}

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ---- Dynamic Regression Forecasting: North Philly -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

{

# ---- fit weekly crimes -------------------------------------------------------------------------------------------------------------------------------------------

# lets set up our data set and external variables

DF = data.frame(tables$crime_Area1_tables$crimes_weekly[Area1 == "NP"])

# lets look at plots of the dataset

mycolors = colorRampPalette(brewer.pal(n = 8, name = "Dark2"))(9)

p1 = ggplot(data = DF, aes(x = DAYS, y = CRIMES, fill = DAYS, color = DAYS)) +
	 geom_violin(alpha = 1/4, width = 5/4) +
	 geom_jitter(alpha = 1/2, width = 1/6) +
	 scale_color_manual(values = mycolors[1:7]) + 
	 scale_fill_manual(values = mycolors[1:7]) +
	 labs(x = "Days", y = "Crimes") +
	 ggtitle("North Philly\nCrimes Per Week v. Days") +
	 theme_light(base_size = 15) +
	 theme(legend.position = "none")

p2 = ggplot(data = DF, aes(x = SUMMER, y = CRIMES, fill = SUMMER, color = SUMMER)) +
	 geom_violin(alpha = 1/4, width = 3/4) +
	 geom_jitter(alpha = 1/2, width = 1/5) +
	 scale_color_manual(values = mycolors[8:9]) + 
	 scale_fill_manual(values = mycolors[8:9]) +
	 labs(x = "Summer", y = "Crimes") +
	 ggtitle("North Philly\nCrimes Per Week v. Summer") +
	 theme_light(base_size = 15) +
	 theme(legend.position = "none")

grid.arrange(p1, p2, ncol = 2)

DF$SUMMER = as.numeric(DF$SUMMER) - 1
DF$DAYS = as.numeric(DF$DAYS)

rm(p1, p2, mycolors)

# lets look at the residuals of a linear model fit using the xreg variables

mod = lm(CRIMES ~ DAYS + SUMMER, data = DF)
mod

# lets look at differenced residuals to find values for d

par(mfrow = c(2,3))

plot(resid(mod), main = "d = 0", type = "b")
lapply(1:5, function(i) plot(diff(resid(mod), differences = i), main = paste0("d = ", i), type ="b"))

d = 1

# lets look at the acf and pacf plots of the differenced residuals to find values of p and q

d = 1
grid.arrange(ggacf(x = diff(resid(mod), differences = d), partial = TRUE, n = 54, basefont = 15, main = paste0("PACF Plot of North Philly Residuals\nDifference of Order ", d)),
			 ggacf(x = diff(resid(mod), differences = d), n = 54, basefont = 15, main = paste0("ACF Plot of North Philly Residuals\nDifference of Order ", d)),
			 ncol = 2)

p = 0:3
q. = 0:1

# lets look at seasonaly differenced residuals to find values for D

par(mfrow = c(2,3))

plot(resid(mod), main = "d = 0", type = "b")
lapply(1:5, function(i) plot(diff(resid(mod), lag = 54, differences = i), main = paste0("D = ", i), type ="b"))

D = 1

rm(d, p, q., D, mod)

# compare auto fits

arima.auto1 = auto.arima(ts(DF$CRIMES, frequency = 54), xreg = DF[,c("SUMMER", "DAYS")])
arima.auto1

c(2, 1, 3, 1, 0, 0)

arima.auto2 = auto.arima(ts(DF$CRIMES, frequency = 54), xreg = DF[,c("SUMMER")])
arima.auto2

c(1, 1, 1, 0, 0, 1)

arima.auto3 = auto.arima(ts(DF$CRIMES, frequency = 54), xreg = DF[,c("DAYS")])
arima.auto3

c(2, 1, 3, 1, 0, 0)

autofits = data.frame("Variables" = rbind("SUMMER,DAYS", "SUMMER", "DAYS"),
					  rbind(accuracy(arima.auto1), accuracy(arima.auto2), accuracy(arima.auto3)),
					  "AIC" = rbind(arima.auto1$aic, arima.auto2$aic, arima.auto3$aic),
					  "AICc" = rbind(arima.auto1$aicc, arima.auto2$aicc, arima.auto3$aicc),
					  "BIC" = rbind(arima.auto1$bic, arima.auto2$bic, arima.auto3$bic))

autofits = autofits[,-c(7, 8)]
autofits

rm(arima.auto1, arima.auto2, arima.auto3)

# lets create our DOE based on the auto fits, acf, pacf, and differencing plots

DOE = expand.grid("p" = 0:3, 
				  "d" = 1, 
				  "q" = 0:1,
				  "P" = 0:1,
				  "D" = 1,
				  "Q" = 0:1)

DOE = rbind(c(2, 1, 3, 1, 0, 0), DOE)
DOE = DOE[!duplicated(DOE),]
rownames(DOE) = 1:NROW(DOE)

head(DOE)
NROW(DOE)

# lets compute the actual and fitted values for each arima model in the DOE

arimas = lapply(1:NROW(DOE), function(i)
								 data.frame("actual" = DF$CRIMES[55:565],
											"fit" = as.numeric(fitted.if(Arima.if(ts(DF$CRIMES, frequency = 54), 
																			order = c(DOE$p[i], DOE$d[i], DOE$q[i]), 
																			seasonal = list(order = c(DOE$P[i], DOE$D[i], DOE$Q[i]), period = 54),
																			xreg = DF[,c("SUMMER", "DAYS")],
																			optim.control = list(maxit = 1000)))[55:565])))
	
# filter out invalid DOE scenarios

invalid = sapply(1:NROW(DOE), function(i) all(is.na(arimas[[i]]$fit)))

DOE = DOE[which(invalid == FALSE),]
rownames(DOE) = 1:NROW(DOE)

arimas = arimas[which(invalid == FALSE)]

rm(invalid)
	
# lower values are ideal

params = sapply(1:NROW(DOE), function(i) 
							 DOE$p[i] + DOE$q[i] + DOE$P[i] + DOE$Q[i])

# values greater than 0.05 are ideal 

ttests = sapply(1:NROW(DOE), function(i) 
							 t.test(arimas[[i]]$actual - arimas[[i]]$fit)$p.value)

# values close to zero are ideal 

stats = data.frame(t(sapply(1:NROW(DOE), function(i)
										 data.frame(accuracy(f = arimas[[i]]$fit, 
															 x = arimas[[i]]$actual), 
													row.names = 1))))

# lets add these results to the DOE

DOE$params = params
DOE$t.pval = ttests
DOE = cbind(DOE, stats)

cnames = colnames(DOE)
DOE[,7:13] = sapply(7:13, function(i) DOE[,i] = as.numeric(DOE[,i]))
colnames(DOE) = cnames

DOE

rm(arimas, params, ttests, stats, cnames)

# lets filter out the top models

# lets plot the histograms of these reponses to find filter limits

par(mfrow = c(2, 4))
lapply(7:13, function(i) hist(DOE[,i], main = colnames(DOE)[i]))

DOE

# apply filters

x = DOE[which(DOE$MAPE <= 7 &
			  DOE$RMSE <= 130 &
			  DOE$params <= 3),]

par(mfrow = c(2, 4))
lapply(7:13, function(i) hist(x[,i], main = colnames(x)[i]))
x

# lets compare the auto fit (scenario 1) with the chosen fit from the DOE (scenario 17)

DOE[c(1, 17),]

# lets build the final models

arimas = lapply(c(1, 17), function(i)
							Arima(ts(DF$CRIMES, frequency = 54), 
									order = c(DOE$p[i], DOE$d[i], DOE$q[i]), 
									seasonal = list(order = c(DOE$P[i], DOE$D[i], DOE$Q[i]), period = 54),
									xreg = DF[,c("SUMMER", "DAYS")],						
									optim.control = list(maxit = 1000)))

# lets plot the model & residuals for the auto fit and chosen fit 

	# auto: scenario 1

plot_crimes_weekly_arimas1 = predplot(obj = arimas[[1]],
							  xlab = "Observation Number", 
							  ylab = "Crimes", 
							  main = "Crimes Per Week: North Philly\nARIMA(2, 1, 3)X(1, 0, 0)54\nRegression Variables: 29.2 * SUMMER, 234.2 * DAYS")

plot_crimes_weekly_arimas1

plot_doe_resid1 = residplots(actual = DF$CRIMES[55:565],
							 fit = as.numeric(fitted(arimas[[1]])[55:565]),
							 histlabel.y = -2,
							 n = 54,
							 basefont = 15)

grid.arrange(plot_doe_resid1[[1]], 
			 plot_doe_resid1[[2]], 
			 plot_doe_resid1[[3]], 
			 plot_doe_resid1[[4]], 
			 plot_doe_resid1[[5]], 
			 plot_doe_resid1[[6]], 
			 ncol = 3)		  
		  
	# chosen: scenario 17

plot_crimes_weekly_arimas2 = predplot(obj = arimas[[2]],
							  xlab = "Observation Number", 
							  ylab = "Crimes", 
							  main = "Crimes Per Week: North Philly\nARIMA(2, 1, 0)X(0, 1, 1)54\nRegression Variables: 7.6 * SUMMER, 234.9 * DAYS")

plot_crimes_weekly_arimas2

plot_doe_resid2 = residplots(actual = DF$CRIMES[55:565],
							 fit = as.numeric(fitted(arimas[[2]])[55:565]), 
							 histlabel.y = -2,
							 n = 54, 
							 basefont = 15)

grid.arrange(plot_doe_resid2[[1]], 
			 plot_doe_resid2[[2]], 
			 plot_doe_resid2[[3]], 
			 plot_doe_resid2[[4]], 
			 plot_doe_resid2[[5]], 
			 plot_doe_resid2[[6]], 
			 ncol = 3)		  

# lets store our results into lists

models$DR$weekly$fits$crimes_weekly_NP = arimas[[2]]
models$DR$weekly$plots$plot_crimes_weekly_NP = plot_crimes_weekly_arimas2
models$DR$weekly$plots$plot_crimes_weekly_residual_NP = plot_doe_resid2
models$DR$weekly$tables$model_data_NP = DF
models$DR$weekly$tables$model_doe_NP = DOE		

rm(x, arimas, plot_crimes_weekly_arimas1, plot_crimes_weekly_arimas2, 
   plot_doe_resid1, plot_doe_resid2, DOE, DF)		  

models$DR$weekly$tables$autofits_NP = autofits

rm(autofits)

# ---- forecast weekly crimes -------------------------------------------------------------------------------------------------------------------------------------------

# lets set up our data set and external variables

DF = data.frame(tables$crime_Area1_tables$crimes_weekly[Area1 == "NP"])
DF$SUMMER = as.numeric(DF$SUMMER) - 1
DF$DAYS = as.numeric(DF$DAYS)

# lets look at the residuals of a linear model fit using the xreg variables

mod = lm(CRIMES ~ DAYS + SUMMER, data = DF[1:478,])
mod

# lets look at differenced residuals to find values for d

par(mfrow = c(2,3))

plot(resid(mod), main = "d = 0", type = "b")
lapply(1:5, function(i) plot(diff(resid(mod), differences = i), main = paste0("d = ", i), type ="b"))

d = 1

# lets look at the acf and pacf plots of the differenced residuals to find values of p and q

d = 1
grid.arrange(ggacf(x = diff(resid(mod), differences = d), partial = TRUE, n = 54, basefont = 15, main = paste0("PACF Plot of North Philly Residuals\nDifference of Order ", d)),
			 ggacf(x = diff(resid(mod), differences = d), n = 54, basefont = 15, main = paste0("ACF Plot of North Philly Residuals\nDifference of Order ", d)),
			 ncol = 2)

p = 0:3
q. = 0:1

# lets look at seasonaly differenced residuals to find values for D

par(mfrow = c(2,3))

plot(resid(mod), main = "d = 0", type = "b")
lapply(1:5, function(i) plot(diff(resid(mod), lag = 54, differences = i), main = paste0("D = ", i), type ="b"))

D = 1

rm(d, p, q., D, mod)

# auto fit

arima.auto = auto.arima(ts(DF$CRIMES[1:478], frequency = 54), xreg = DF[1:478, c("SUMMER", "DAYS")])
arima.auto

c(0, 1, 1, 1, 0, 0)

rm(arima.auto)

# lets create our DOE based on the auto.arima, acf, pacf, and differencing plots

DOE = expand.grid("p" = 0:3, 
				  "d" = 1, 
				  "q" = 0:1,
				  "P" = 0:1,
				  "D" = 1,
				  "Q" = 0:1)

DOE = rbind(c(0, 1, 1, 1, 0, 0), DOE)
DOE = DOE[!duplicated(DOE),]
rownames(DOE) = 1:NROW(DOE)

head(DOE)
NROW(DOE)

# lets compute the arima models from the DOE

arimas = lapply(1:NROW(DOE), function(i)
								Arima.if(ts(DF$CRIMES[1:478], frequency = 54), 
											order = c(DOE$p[i], DOE$d[i], DOE$q[i]), 
											seasonal = list(order = c(DOE$P[i], DOE$D[i], DOE$Q[i]), period = 54),
											xreg = DF[1:478, c("SUMMER", "DAYS")],
											optim.control = list(maxit = 1000)))
		
# filter out invalid DOE scenarios

invalid = sapply(1:NROW(DOE), function(i) is.na(arimas[i]))

DOE = DOE[which(invalid == FALSE),]
rownames(DOE) = 1:NROW(DOE)

arimas = arimas[which(invalid == FALSE)]

rm(invalid)		

# lets compute the actual and fitted values for the testing data of each arima model in the DOE

arimas = lapply(1:NROW(DOE), function(i)
								data.frame("actual" = DF$CRIMES[479:565],
											"fit" = as.numeric(forecast(arimas[[i]],
																		h = 87,
																		xreg = DF[479:565, c("SUMMER", "DAYS")],
																		level = c(0, 0))$mean)))

# lower values are ideal

params = sapply(1:NROW(DOE), function(i) 
							 DOE$p[i] + DOE$q[i] + DOE$P[i] + DOE$Q[i])

# values greater than 0.05 are ideal 

ttests = sapply(1:NROW(DOE), function(i) 
							 t.test(arimas[[i]]$actual - arimas[[i]]$fit)$p.value)

# values close to zero are ideal 

stats = data.frame(t(sapply(1:NROW(DOE), function(i)
										 data.frame(accuracy(f = arimas[[i]]$fit, 
															 x = arimas[[i]]$actual), 
													row.names = 1))))

# lets add these results to the DOE

DOE$params = params
DOE$t.pval = ttests
DOE = cbind(DOE, stats)

cnames = colnames(DOE)
DOE[,7:13] = sapply(7:13, function(i) DOE[,i] = as.numeric(DOE[,i]))
colnames(DOE) = cnames

DOE

rm(arimas, params, ttests, stats, cnames)

# lets filter out the top models

# lets plot the histograms of these reponses to find filter limits

par(mfrow = c(2, 4))
lapply(7:13, function(i) hist(DOE[,i], main = colnames(DOE)[i]))

DOE

# apply filters

x = DOE[which(DOE$MAPE <= 8 &
			  DOE$t.pval >= 0.3),]

par(mfrow = c(2, 4))
lapply(7:13, function(i) hist(x[,i], main = colnames(x)[i]))

x

# lets compare the auto fit (scenario 1) with the chosen fit from the DOE (scenario 20)

DOE[c(1, 20),]

# lets build the final models

arimas = lapply(c(1, 20), function(i)
							Arima(ts(DF$CRIMES[1:478], frequency = 54), 
									order = c(DOE$p[i], DOE$d[i], DOE$q[i]), 
									seasonal = list(order = c(DOE$P[i], DOE$D[i], DOE$Q[i]), period = 54),
									xreg = DF[1:478, c("SUMMER", "DAYS")],
									optim.control = list(maxit = 1000)))

# lets plot the model & residuals for the auto fit and chosen fit 

	# auto: scenario 1

plot_crimes_weekly_arimas1 = predplot(obj = arimas[[1]],
							  n.ahead = 87,
							  xreg = DF[479:565, c("SUMMER", "DAYS")],
							  actuals = DF$CRIMES,
							  xlab = "Observation Number", 
							  ylab = "Crimes", 
							  main = "Crimes Per Week: North Philly\nARIMA(0, 1, 1)X(1, 0, 0)54\nRegression Variable: 33.7 * SUMMER, 240.7 * DAYS")

plot_crimes_weekly_arimas1

plot_doe_resid1 = residplots(actual = DF$CRIMES[479:565],
							 fit = as.numeric(predict(arimas[[1]], newxreg = DF[479:565, c("SUMMER", "DAYS")], n.ahead = 87, se.fit = FALSE)),
							 histlabel.y = -0.5,
							 n = 54,
							 binwidth = 40,
							 basefont = 15)

grid.arrange(plot_doe_resid1[[1]], 
			 plot_doe_resid1[[2]], 
			 plot_doe_resid1[[3]], 
			 plot_doe_resid1[[4]], 
			 plot_doe_resid1[[5]], 
			 plot_doe_resid1[[6]], 
			 ncol = 3)		  
		  
	# chosen: scenario 20

plot_crimes_weekly_arimas2 = predplot(obj = arimas[[2]],
							  n.ahead = 87,
							  xreg = DF[479:565, c("SUMMER", "DAYS")],
							  actuals = DF$CRIMES,
							  xlab = "Observation Number", 
							  ylab = "Crimes", 
							  main = "Crimes Per Week: North Philly\nARIMA(1, 1, 1)X(0, 1, 1)54\nRegression Variable: 0.19 * SUMMER, 232.8 * DAYS")

plot_crimes_weekly_arimas2

plot_doe_resid2 = residplots(actual = DF$CRIMES[479:565],
							 fit = as.numeric(predict(arimas[[2]], newxreg = DF[479:565, c("SUMMER", "DAYS")], n.ahead = 87, se.fit = FALSE)),
							 histlabel.y = -0.5,
							 n = 54,
							 binwidth = 40,
							 basefont = 15)

grid.arrange(plot_doe_resid2[[1]], 
			 plot_doe_resid2[[2]], 
			 plot_doe_resid2[[3]], 
			 plot_doe_resid2[[4]], 
			 plot_doe_resid2[[5]], 
			 plot_doe_resid2[[6]], 
			 ncol = 3)		  

# lets store our results into lists

models$DR$weekly$fits$crimes_weekly_NP_forecast = arimas[[2]]
models$DR$weekly$plots$plot_crimes_weekly_NP_forecast = plot_crimes_weekly_arimas2
models$DR$weekly$plots$plot_crimes_weekly_residual_NP_forecast = plot_doe_resid2
models$DR$weekly$tables$model_data_NP_forecast = DF
models$DR$weekly$tables$model_doe_NP_forecast = DOE		  

rm(x, arimas, plot_crimes_weekly_arimas1, plot_crimes_weekly_arimas2, 
   plot_doe_resid1, plot_doe_resid2, DOE, DF)

}

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ---- Dynamic Regression Forecasting with Lagged Regression Variables: North Philly -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

{

# ---- fit weekly crimes -------------------------------------------------------------------------------------------------------------------------------------------

# lets set up our data set and external variables
# we will compare three different lags for each of the two regression variables
	# a lag of 1, 4, and 54 weeks
	
DF = data.frame(tables$crime_Area1_tables$crimes_weekly[Area1 == "NP"])
DF$SUMMER = as.numeric(DF$SUMMER) - 1
DF$DAYS = as.numeric(DF$DAYS)

DF = data.table(DF)

DF[,SUMMER1 := shift(SUMMER, 1)]
DF[,SUMMER4 := shift(SUMMER, 4)]
DF[,SUMMER54 := shift(SUMMER, 54)]

DF[,DAYS1 := shift(DAYS, 1)]
DF[,DAYS4 := shift(DAYS, 4)]
DF[,DAYS54 := shift(DAYS, 54)]

setcolorder(DF, c(1, 2, 3, 6, 7, 8, 4, 9, 10, 11, 5))
DF = data.frame(DF)

head(DF)

# lets create a DOE to test all possible regression lag combinations
# we will choose the best combination based on the R^2-Pred, AIC, and BIC values of linear regression models

DOE = expand.grid("SUMMER" = 4:6,
				  "DAYS" = 8:10)

DOE$r2pred = sapply(1:NROW(DOE), function(i)
					statslm(lm(DF$CRIMES[55:565] ~ DF[55:565,DOE$DAYS[i]] + DF[55:565,DOE$SUMMER[i]]))$Prediction)

DOE = cbind(DOE, t(sapply(1:NROW(DOE), function(i)
							statslm(lm(DF$CRIMES[55:565] ~ DF[55:565,DOE$DAYS[i]] + DF[55:565,DOE$SUMMER[i]]))$Fitness)))

DOE$SUMMER = colnames(DF)[DOE$SUMMER]
DOE$DAYS = colnames(DF)[DOE$DAYS]

DOE

# update the dataset with the best lagged regression variables

DF = DF[,c(1,2,4,10,11)]
head(DF)

rm(DOE)

# lets look at the residuals of a linear model fit using the xreg variables

mod = lm(CRIMES ~ DAYS54 + SUMMER1, data = DF)
mod

# lets look at differenced residuals to find values for d

par(mfrow = c(2,3))

plot(resid(mod), main = "d = 0", type = "b")
lapply(1:5, function(i) plot(diff(resid(mod), differences = i), main = paste0("d = ", i), type ="b"))

d = 1

# lets look at the acf and pacf plots of the differenced residuals to find values of p and q

d = 1
grid.arrange(ggacf(x = diff(resid(mod), differences = d), partial = TRUE, n = 54, basefont = 15, main = paste0("PACF Plot of North Philly Residuals\nDifference of Order ", d)),
			 ggacf(x = diff(resid(mod), differences = d), n = 54, basefont = 15, main = paste0("ACF Plot of North Philly Residuals\nDifference of Order ", d)),
			 ncol = 2)

p = 0:5
q. = 0:2

# lets look at seasonaly differenced residuals to find values for D

par(mfrow = c(2,3))

plot(resid(mod), main = "d = 0", type = "b")
lapply(1:5, function(i) plot(diff(resid(mod), lag = 54, differences = i), main = paste0("D = ", i), type ="b"))

D = 1

rm(d, p, q., D, mod)

# lets look at the auto fit

arima.auto = auto.arima(ts(DF$CRIMES, frequency = 54), xreg = DF[,c("SUMMER1", "DAYS54")])
arima.auto

c(1, 1, 1, 0, 0, 2)

# lets create our DOE based on the auto fit, acf, pacf, and differencing plots

DOE = expand.grid("p" = 0:4, 
				  "d" = 1, 
				  "q" = 0:1,
				  "P" = 0:1,
				  "D" = 1,
				  "Q" = 1)

DOE = rbind(c(1, 1, 1, 0, 0, 2), DOE)
DOE = DOE[!duplicated(DOE),]
rownames(DOE) = 1:NROW(DOE)

head(DOE)
NROW(DOE)

# lets compute the actual and fitted values for each arima model in the DOE

arimas = lapply(1:NROW(DOE), function(i)
								 data.frame("actual" = DF$CRIMES[109:565],
											"fit" = as.numeric(fitted.if(Arima.if(ts(DF$CRIMES, frequency = 54), 
																			order = c(DOE$p[i], DOE$d[i], DOE$q[i]), 
																			seasonal = list(order = c(DOE$P[i], DOE$D[i], DOE$Q[i]), period = 54),
																			xreg = DF[,c("SUMMER1", "DAYS54")],
																			optim.control = list(maxit = 1000)))[109:565])))
	
# filter out invalid DOE scenarios

invalid = sapply(1:NROW(DOE), function(i) all(is.na(arimas[[i]]$fit)))

DOE = DOE[which(invalid == FALSE),]
rownames(DOE) = 1:NROW(DOE)

arimas = arimas[which(invalid == FALSE)]

rm(invalid)
	
# lower values are ideal

params = sapply(1:NROW(DOE), function(i) 
							 DOE$p[i] + DOE$q[i] + DOE$P[i] + DOE$Q[i])

# values greater than 0.05 are ideal 

ttests = sapply(1:NROW(DOE), function(i) 
							 t.test(arimas[[i]]$actual - arimas[[i]]$fit)$p.value)

# values close to zero are ideal 

stats = data.frame(t(sapply(1:NROW(DOE), function(i)
										 data.frame(accuracy(f = arimas[[i]]$fit, 
															 x = arimas[[i]]$actual), 
													row.names = 1))))

# lets add these results to the DOE

DOE$params = params
DOE$t.pval = ttests
DOE = cbind(DOE, stats)

cnames = colnames(DOE)
DOE[,7:13] = sapply(7:13, function(i) DOE[,i] = as.numeric(DOE[,i]))
colnames(DOE) = cnames

DOE

rm(arimas, params, ttests, stats, cnames)

# lets filter out the top models

# lets plot the histograms of these reponses to find filter limits

par(mfrow = c(2, 4))
lapply(7:13, function(i) hist(DOE[,i], main = colnames(DOE)[i]))

DOE

# apply filters

x = DOE[which(DOE$RMSE <= 220),]

par(mfrow = c(2, 4))
lapply(7:13, function(i) hist(x[,i], main = colnames(x)[i]))
x

# lets compare the auto fit (scenario 1) with the chosen fit from the DOE (scenario 17)

DOE[c(1, 17),]

# lets build the final models

arimas = lapply(c(1, 17), function(i)
							Arima(ts(DF$CRIMES, frequency = 54), 
									order = c(DOE$p[i], DOE$d[i], DOE$q[i]), 
									seasonal = list(order = c(DOE$P[i], DOE$D[i], DOE$Q[i]), period = 54),
									xreg = DF[,c("SUMMER1", "DAYS54")],						
									optim.control = list(maxit = 1000)))

# lets plot the model & residuals for the auto fit and chosen fit 

	# auto: scenario 1

plot_crimes_weekly_arimas1 = predplot(obj = arimas[[1]],
							  xlab = "Observation Number", 
							  ylab = "Crimes", 
							  main = "Crimes Per Week: North Philly\nARIMA(1, 1, 1)X(0, 0, 2)54\nRegression Variables: 28.2 * SUMMER1, 6.5 * DAYS54")

plot_crimes_weekly_arimas1

plot_doe_resid1 = residplots(actual = DF$CRIMES[109:565],
							 fit = as.numeric(fitted(arimas[[1]])[109:565]),
							 histlabel.y = -3,
							 n = 54,
							 basefont = 15)

grid.arrange(plot_doe_resid1[[1]], 
			 plot_doe_resid1[[2]], 
			 plot_doe_resid1[[3]], 
			 plot_doe_resid1[[4]], 
			 plot_doe_resid1[[5]], 
			 plot_doe_resid1[[6]], 
			 ncol = 3)		  
		  
	# chosen: scenario 21

plot_crimes_weekly_arimas2 = predplot(obj = arimas[[2]],
							  xlab = "Observation Number", 
							  ylab = "Crimes", 
							  main = "Crimes Per Week: North Philly\nARIMA(0, 1, 1)X(1, 1, 1)54\nRegression Variables: -12.2 * SUMMER1, 0.9 * DAYS54")

plot_crimes_weekly_arimas2

plot_doe_resid2 = residplots(actual = DF$CRIMES[109:565],
							 fit = as.numeric(fitted(arimas[[2]])[109:565]), 
							 histlabel.y = -2,
							 n = 54, 
							 basefont = 15)

grid.arrange(plot_doe_resid2[[1]], 
			 plot_doe_resid2[[2]], 
			 plot_doe_resid2[[3]], 
			 plot_doe_resid2[[4]], 
			 plot_doe_resid2[[5]], 
			 plot_doe_resid2[[6]], 
			 ncol = 3)		  

# lets store our results into lists

models$DR$weekly$fits$crimes_weekly_NP_lag = arimas[[2]]
models$DR$weekly$plots$plot_crimes_weekly_NP_lag = plot_crimes_weekly_arimas2
models$DR$weekly$plots$plot_crimes_weekly_residual_NP_lag = plot_doe_resid2
models$DR$weekly$tables$model_data_NP_lag = DF
models$DR$weekly$tables$model_doe_NP_lag = DOE		

rm(x, arimas, plot_crimes_weekly_arimas1, plot_crimes_weekly_arimas2, 
   plot_doe_resid1, plot_doe_resid2, DOE, DF)		  

# ---- forecast weekly crimes -------------------------------------------------------------------------------------------------------------------------------------------

# lets set up our data set and external variables

DF = data.frame(tables$crime_Area1_tables$crimes_weekly[Area1 == "NP"])
DF$SUMMER = as.numeric(DF$SUMMER) - 1
DF$DAYS = as.numeric(DF$DAYS)

DF = data.table(DF)
DF[,SUMMER1 := shift(SUMMER, 1)]
DF[,DAYS54 := shift(DAYS, 54)]
setcolorder(DF, c(1:4, 6, 7, 5))

DF = data.frame(DF)
DF = DF[,-c(3:4)]

# lets look at the residuals of a linear model fit using the xreg variables

mod = lm(CRIMES ~ DAYS54 + SUMMER1, data = DF[1:478,])
mod

# lets look at differenced residuals to find values for d

par(mfrow = c(2,3))

plot(resid(mod), main = "d = 0", type = "b")
lapply(1:5, function(i) plot(diff(resid(mod), differences = i), main = paste0("d = ", i), type ="b"))

d = 1

# lets look at the acf and pacf plots of the differenced residuals to find values of p and q

d = 1
grid.arrange(ggacf(x = diff(resid(mod), differences = d), partial = TRUE, n = 54, basefont = 15, main = paste0("PACF Plot of North Philly Residuals\nDifference of Order ", d)),
			 ggacf(x = diff(resid(mod), differences = d), n = 54, basefont = 15, main = paste0("ACF Plot of North Philly Residuals\nDifference of Order ", d)),
			 ncol = 2)

p = 0:3
q. = 0:1

# lets look at seasonaly differenced residuals to find values for D

par(mfrow = c(2,3))

plot(resid(mod), main = "d = 0", type = "b")
lapply(1:5, function(i) plot(diff(resid(mod), lag = 54, differences = i), main = paste0("D = ", i), type ="b"))

D = 1

rm(d, p, q., D, mod)

# auto fit

arima.auto = auto.arima(ts(DF$CRIMES[1:478], frequency = 54), xreg = DF[1:478, c("SUMMER1", "DAYS54")])
arima.auto

c(0, 1, 2, 1, 0, 0)

rm(arima.auto)

# lets create our DOE based on the auto.arima, acf, pacf, and differencing plots

DOE = expand.grid("p" = 0:3, 
				  "d" = 1, 
				  "q" = 0:1,
				  "P" = 1,
				  "D" = 1,
				  "Q" = 0:1)

DOE = rbind(c(0, 1, 2, 1, 0, 0), DOE)
DOE = DOE[!duplicated(DOE),]
rownames(DOE) = 1:NROW(DOE)

head(DOE)
NROW(DOE)

# lets compute the arima models from the DOE

arimas = lapply(1:NROW(DOE), function(i)
								Arima.if(ts(DF$CRIMES[1:478], frequency = 54), 
											order = c(DOE$p[i], DOE$d[i], DOE$q[i]), 
											seasonal = list(order = c(DOE$P[i], DOE$D[i], DOE$Q[i]), period = 54),
											xreg = DF[1:478, c("SUMMER1", "DAYS54")],
											optim.control = list(maxit = 1000)))
		
# filter out invalid DOE scenarios

invalid = sapply(1:NROW(DOE), function(i) is.na(arimas[i]))

DOE = DOE[which(invalid == FALSE),]
rownames(DOE) = 1:NROW(DOE)

arimas = arimas[which(invalid == FALSE)]

rm(invalid)		

# lets compute the actual and fitted values for the testing data of each arima model in the DOE

arimas = lapply(1:NROW(DOE), function(i)
								data.frame("actual" = DF$CRIMES[479:565],
											"fit" = as.numeric(forecast(arimas[[i]],
																		h = 87,
																		xreg = DF[479:565, c("SUMMER1", "DAYS54")],
																		level = c(0, 0))$mean)))

# lower values are ideal

params = sapply(1:NROW(DOE), function(i) 
							 DOE$p[i] + DOE$q[i] + DOE$P[i] + DOE$Q[i])

# values greater than 0.05 are ideal 

ttests = sapply(1:NROW(DOE), function(i) 
							 t.test(arimas[[i]]$actual - arimas[[i]]$fit)$p.value)

# values close to zero are ideal 

stats = data.frame(t(sapply(1:NROW(DOE), function(i)
										 data.frame(accuracy(f = arimas[[i]]$fit, 
															 x = arimas[[i]]$actual), 
													row.names = 1))))

# lets add these results to the DOE

DOE$params = params
DOE$t.pval = ttests
DOE = cbind(DOE, stats)

cnames = colnames(DOE)
DOE[,7:13] = sapply(7:13, function(i) DOE[,i] = as.numeric(DOE[,i]))
colnames(DOE) = cnames

DOE

rm(arimas, params, ttests, stats, cnames)

# lets filter out the top models

# lets plot the histograms of these reponses to find filter limits

par(mfrow = c(2, 4))
lapply(7:13, function(i) hist(DOE[,i], main = colnames(DOE)[i]))

DOE

# apply filters

x = DOE[which(DOE$MAPE <= 20 &
			  DOE$RMSE <= 230),]

par(mfrow = c(2, 4))
lapply(7:13, function(i) hist(x[,i], main = colnames(x)[i]))

x

# lets compare the auto fit (scenario 1) with the chosen fit from the DOE (scenario 13)

DOE[c(1, 13),]

# lets build the final models

arimas = lapply(c(1, 13), function(i)
							Arima(ts(DF$CRIMES[1:478], frequency = 54), 
									order = c(DOE$p[i], DOE$d[i], DOE$q[i]), 
									seasonal = list(order = c(DOE$P[i], DOE$D[i], DOE$Q[i]), period = 54),
									xreg = DF[1:478, c("SUMMER1", "DAYS54")],
									optim.control = list(maxit = 1000)))

# lets plot the model & residuals for the auto fit and chosen fit 

	# auto: scenario 1

plot_crimes_weekly_arimas1 = predplot(obj = arimas[[1]],
							  n.ahead = 87,
							  xreg = DF[479:565, c("SUMMER1", "DAYS54")],
							  actuals = DF$CRIMES,
							  xlab = "Observation Number", 
							  ylab = "Crimes", 
							  main = "Crimes Per Week: North Philly\nARIMA(0, 1, 2)X(1, 0, 0)54\nRegression Variable: 46.2 * SUMMER1, -4.0 * DAYS54")

plot_crimes_weekly_arimas1

plot_doe_resid1 = residplots(actual = DF$CRIMES[479:565],
							 fit = as.numeric(predict(arimas[[1]], newxreg = DF[479:565, c("SUMMER1", "DAYS54")], n.ahead = 87, se.fit = FALSE)),
							 histlabel.y = -0.5,
							 n = 54,
							 binwidth = 50,
							 basefont = 15)

grid.arrange(plot_doe_resid1[[1]], 
			 plot_doe_resid1[[2]], 
			 plot_doe_resid1[[3]], 
			 plot_doe_resid1[[4]], 
			 plot_doe_resid1[[5]], 
			 plot_doe_resid1[[6]], 
			 ncol = 3)		  
		  
	# chosen: scenario 13

plot_crimes_weekly_arimas2 = predplot(obj = arimas[[2]],
							  n.ahead = 87,
							  xreg = DF[479:565, c("SUMMER1", "DAYS54")],
							  actuals = DF$CRIMES,
							  xlab = "Observation Number", 
							  ylab = "Crimes", 
							  main = "Crimes Per Week: North Philly\nARIMA(0, 1, 1)X(1, 1, 1)54\nRegression Variable: -16.3 * SUMMER1, -18.9 * DAYS54")

plot_crimes_weekly_arimas2

plot_doe_resid2 = residplots(actual = DF$CRIMES[479:565],
							 fit = as.numeric(predict(arimas[[2]], newxreg = DF[479:565, c("SUMMER1", "DAYS54")], n.ahead = 87, se.fit = FALSE)),
							 histlabel.y = -0.5,
							 n = 54,
							 binwidth = 50,
							 basefont = 15)

grid.arrange(plot_doe_resid2[[1]], 
			 plot_doe_resid2[[2]], 
			 plot_doe_resid2[[3]], 
			 plot_doe_resid2[[4]], 
			 plot_doe_resid2[[5]], 
			 plot_doe_resid2[[6]], 
			 ncol = 3)		  

# lets store our results into lists

models$DR$weekly$fits$crimes_weekly_NP_lag_forecast = arimas[[2]]
models$DR$weekly$plots$plot_crimes_weekly_NP_lag_forecast = plot_crimes_weekly_arimas2
models$DR$weekly$plots$plot_crimes_weekly_residual_NP_lag_forecast = plot_doe_resid2
models$DR$weekly$tables$model_data_NP_lag_forecast = DF
models$DR$weekly$tables$model_doe_NP_lag_forecast = DOE		  

rm(x, arimas, plot_crimes_weekly_arimas1, plot_crimes_weekly_arimas2, 
   plot_doe_resid1, plot_doe_resid2, DOE, DF)		  

}

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ---- Project Part 3 Working Space -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

{

# ------------------------------------------------

 + scale_x_continuous(limits = c(479, 565))

 + theme_light(18) + theme(legend.position = "none") + scale_x_continuous(limits = c(479, 565))

# ------------------------------------------------

myplot = models$ARIMA$weekly$plots$plot_crimes_weekly_residual_NP

grid.arrange(myplot[[1]], 
			 myplot[[2]], 
			 myplot[[3]], 
			 myplot[[5]], 
			 ncol = 2)	

myplot[[6]]

rm(myplot)

# ------------------------------------------------

hwdoe = models$HW$weekly$tables$model_doe_NP_forecast
arimadoe = models$ARIMA$weekly$tables$model_doe_NP_forecast

hwbest = 7257
arimabest = 31

mydoe = data.frame("Model" = c("HW", "ARIMA"),
				   rbind(hwdoe[hwbest,4:9],
						 arimadoe[arimabest,8:13]))

mydoe

rm(hwdoe, arimadoe, hwbest, arimabest, mydoe)

# ------------------------------------------------

drdoe = models$DR$weekly$tables$model_doe_NP
drlagdoe = models$DR$weekly$tables$model_doe_NP_lag

drbest = 17
drlagbest = 17

mydoe = data.frame("Model" = c("DR", "DR.Lag"),
				   rbind(drdoe[drbest,8:13],
						 drlagdoe[drlagbest,8:13]))

mydoe

rm(drdoe, drlagdoe, drbest, drlagbest, mydoe)

# ------------------------------------------------

hwdoe = models$HW$weekly$tables$model_doe_NP_forecast
arimadoe = models$ARIMA$weekly$tables$model_doe_NP_forecast
drdoe = models$DR$weekly$tables$model_doe_NP_forecast
drlagdoe = models$DR$weekly$tables$model_doe_NP_lag_forecast

hwbest = 7257
arimabest = 31
drbest = 20
drlagbest = 13

mydoe = data.frame("Model" = c("HW", "ARIMA", "DR", "DR.Lag"),
				   rbind(hwdoe[hwbest,4:9],
						 arimadoe[arimabest,8:13],
						 drdoe[drbest,8:13],
						 drlagdoe[drlagbest,8:13]))
						 
mydoe					 

rm(hwdoe, arimadoe, drdoe, drlagdoe, hwbest, arimabest, drbest, drlagbest, mydoe)

}












