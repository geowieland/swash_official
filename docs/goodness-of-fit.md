# Goodness-of-Fit Metrics

## Background

Goodness-of-fit metrics are used to evaluate how well a model's predicted values match the observed values. They quantify different aspects of model performance, including the proportion of explained variance, the magnitude of prediction errors, relative deviations, and the distribution of errors across observations[<a href="#ref-1">1</a>][<a href="#ref-2">2</a>]. The model fit indicators are relevant for the empirical parametrization approaches of the models of infection growth over time, e.g. [logistic](logistic-growth.md) or [exponential growth model](exponential-growth.md). In these cases, we are dealing with [time series](infections-data.md#time-series-data) consisting of $$T$$ time points; the number of observations is therefore equal to $$T$$.

Throughout the following definitions, the notation is:

- $$y_i$$: Observed value of observation $$i$$
- $$\hat{y}_i$$: Expected (predicted) value of observation $$i$$
- $$e_i = y_i - \hat{y}_i$$: Residual (prediction error)
- $$\bar{y}$$: Mean of the observed values
- $$n$$: Number of valid observed–expected pairs
- $$p$$: Percentage scaling factor (typically 100)

The following list and the formulas are based on academic literature from the fields of econometrics and machine learning[<a href="#ref-1">1</a>][<a href="#ref-2">2</a>]:

## Variance-Based Metrics

Variance-based metrics evaluate the extent to which a model explains the variability in the observed data. They compare the residual variation against the total variation of the observations.

### Coefficient of Determination (R-squared, R²)

The coefficient of determination measures the proportion of the total variance in the observed data that is explained by the model. Values closer to 1 indicate a better fit, whereas a value of 0 indicates that the model performs no better than predicting the mean of the observations. Negative values are possible and indicate that the model performs worse than this mean-based reference prediction.

$$R^2 = 1 - \frac{\sum_{i=1}^{n}(y_i-\hat{y}_i)^2} {\sum_{i=1}^{n}(y_i-\bar{y})^2}$$

**Interpretation:** Higher values indicate better model performance. The theoretical maximum is 1, while negative values are unbounded.

## Deviation-Based Metrics

Deviation-based metrics quantify the magnitude of the differences between observed and predicted values. They are expressed in the original unit of measurement (except for squared error metrics) and provide a direct assessment of prediction accuracy.

### Sum of Squared Residuals (SSR)

The sum of squared residuals represents the total squared difference between observed and predicted values. Squaring the residuals gives larger errors disproportionately more weight, making this metric particularly sensitive to outliers. Lower values indicate a better fit.

$$SSR = \sum_{i=1}^{n}(y_i-\hat{y}_i)^2$$

**Interpretation:** Values range from 0 to infinity, with 0 representing a perfect fit.

### Sum of Absolute Residuals (SAR)

The sum of absolute residuals measures the total magnitude of prediction errors without considering their direction. Unlike squared residuals, absolute residuals do not disproportionately penalize large errors. Lower values indicate a better fit.

$$SAR = \sum_{i=1}^{n}|y_i-\hat{y}_i|$$

**Interpretation:** Values range from 0 to infinity, with 0 representing a perfect fit.

### Mean Squared Error (MSE)

The mean squared error calculates the average squared difference between observed and predicted values. Due to the squaring operation, large prediction errors receive a disproportionately high weight. Lower values indicate better model performance.

$$MSE = \frac{1}{n}\sum_{i=1}^{n}(y_i-\hat{y}_i)^2$$

**Interpretation:** Values range from 0 to infinity, with 0 indicating a perfect fit. The metric is expressed in squared units of the observed variable.

### Root Mean Squared Error (RMSE)

The root mean squared error is the square root of the mean squared error and therefore expresses prediction errors in the original unit of measurement. Similar to MSE, it penalizes larger errors more strongly than smaller ones. Lower values indicate better model performance.

$$RMSE = \sqrt{\frac{1}{n}\sum_{i=1}^{n}(y_i-\hat{y}_i)^2}$$

**Interpretation:** Values range from 0 to infinity, with 0 indicating a perfect fit.

### Mean Absolute Error (MAE)

The mean absolute error represents the average absolute difference between observed and predicted values. It provides an intuitive measure of the typical prediction error in the original unit of measurement and is less sensitive to outliers than MSE or RMSE. Lower values indicate better model performance.

$$MAE = \frac{1}{n}\sum_{i=1}^{n}|y_i-\hat{y}_i|$$

**Interpretation:** Values range from 0 to infinity, with 0 indicating a perfect fit.

## Relative and Percentage-Based Metrics

Relative and percentage-based metrics express prediction errors in relation to the magnitude of the observed or predicted values. They facilitate comparisons between datasets or variables with different scales and units.

### Mean Absolute Percentage Error (MAPE)

The mean absolute percentage error measures the average absolute prediction error relative to the observed values, expressed as a percentage. It provides a scale-independent measure of prediction accuracy. However, MAPE is undefined for observations with an observed value of zero and can become disproportionately large when observed values are close to zero.

$$MAPE = \frac{p}{n^*} \sum_{i \in I^*} \left|\frac{y_i-\hat{y}_i}{y_i}\right|$$

where $$I^*$$ is the set of observations for which $$y_i \neq 0$$, and $$n^*$$ is the number of such observations.

**Interpretation:** Lower values indicate better model performance, with 0 representing a perfect fit. The metric is expressed as a percentage when \(p=100\).

Observations with zero observed values must be excluded from the MAPE calculation.

### Absolute Percentage Error (APE)

The absolute percentage error represents the absolute deviation of an individual prediction relative to its observed value. It is calculated for each observation and can be used to examine the distribution of prediction errors across individual data points. Lower values indicate more accurate predictions.

$$APE_i = \left|\frac{y_i-\hat{y}_i}{y_i}\right|\cdot p$$

**Interpretation:** Values are non-negative, with 0% indicating a perfect prediction. APE is undefined for observed values equal to zero and is therefore excluded for these observations in the implementation.

APE threshold statistics describe the proportion of predictions whose absolute percentage error falls below a predefined threshold. They provide an intuitive way to assess how frequently a model achieves a specified level of relative prediction accuracy. Higher percentages indicate better performance.

For a threshold $$t$$, the statistic is calculated as:

$$APE_{<t} = \frac{1}{n} \sum_{i=1}^{n} \mathbb{1}(APE_i<t)\cdot p$$

where $$\mathbb{1}(\cdot)$$ is the indicator function, which equals 1 if the condition is satisfied and 0 otherwise.

**Interpretation:** This calculation allows us to answer the following question: *For what percentage of $$N$$ empirical observations is the deviation in the model at most $$t$$ percent?* Values range from 0% to 100% when $$p=100$$. Higher values indicate that a larger proportion of predictions fall below the specified percentage-error threshold.

Observations with zero observed values must be excluded from the APE calculation.


## References

<a id="ref-1"></a>
[1] Boehmke B, Greenwell B (2020) 
Hands-On Machine Learning with R. 1 ed. 
New York, NY: Taylor & Francis. [10.1201/9780367816377](https://doi.org/10.1201/9780367816377)

<a id="ref-2"></a>
[2] Buturac G (2022)
Measurement of Economic Forecast Accuracy: A Systematic Overview of the Empirical Literature. 
*Journal of Risk and Financial Management* 15(1), 1. [10.3390/jrfm15010001](https://doi.org/10.3390/jrfm15010001)