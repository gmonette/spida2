# Let's randomly generate some nonlinear data:

set.seed(12345) # for reproducibility
x <- runif(200, 0, 10)
x <- sort(x)
Ey <- cos(1.25*(x + 1)) + x/5
y <- Ey + 0.5*rnorm(200)
plot(x, y)
f <- function(x) cos(1.25*(x + 1)) + x/5
curve(f, add=TRUE, lty=2, lwd=2) # population regression
```

The line on the graph represents the "true" population regression curve.

_Piecewise polynomials_ begin by dividing the range of $x$ into non-overlapping intervals at $k$ ordered $x$-values $t_j$ called _knots_: $(-\inf, t_1], (t_1, t_2], \ldots, (t_k, \inf)$. Then a degree-$p$ polynomial is fit by least-squares regression to the data in each interval. The simplest case is an degree-0 polynomial, which generates a _piecewise-constant_ fit by computing the mean $y$-value in each interval. For the example, we set $k = 2$ knots at $t_1 = 10/3$ and $t_2 = 2\times10/3$, which are the 1/3 and 2/3 quantiles of the uniform[0, 10] distribution from which the $x$-values were generated. For real data, we can pick knots using sample $x$-quantiles or (as we discuss in another vignette) theoretically interesting transitional points in the data.
```{r}
plot(x, y, main="Piecewise Constant Fit")
curve(f, add=TRUE, lty=2, lwd=2)
t <- c(10/3, 20/3) # knots
abline(v=t, lty=3)
mtext(c(expression(t[1]), expression(t[2])), side=1, line=0.5, at=t) 
group <- cut(x, breaks=c(-Inf, t, Inf))
means <- tapply(y, group, mean) # within-interval means
x0 <- c(0, t, 10) # knots augmented with boundaries
for (i in 1:3) lines(x0[c(i, i+1)], rep(means[i], 2), lwd=2)
yhat1 <- means[group] # fitted values
mean((Ey - yhat1)^2) # RMSE
```

RMSE is the root-mean-squared error of estimating $E(y)$ for the cases in the data. It provides a rough guide to how well the model recovers the population regression curve, with the caveats that (1) it is an _estimate_ because we have only $n = 200$ sampled cases, and (2) we haven't penalized the RMSE for the number of parameters in the model.

A next step is to generalize the model to a _piecewise-linear_ fit (that is, a piecewise degree-1 polynomial):
```{r}
mods2 <- by(data.frame(x=x, y=y, group=group), group, 
    function(df) lm(y ~ x, data=df)) # within-group linear regressions
plot(x, y, main="Piecewise Linear Fit")
curve(f, add=TRUE, lty=2, lwd=2)
abline(v=t, lty=3)
mtext(c(expression(t[1]), expression(t[2])), side=1, line=0.5, at=t) 
for (j in 1:3){
    x1 <- x0[c(j, j + 1)]
    y1 <- predict(mods2[[j]], list(x=x1))
    lines(x1, y1, lwd=2) # draw within-group regression lines
}
yhat2 <- unlist(lapply(mods2, fitted))
mean((Ey - yhat2)^2) # RMSE
```

We convert the piecewise-linear fit into a _linear regression spline_ by constraining the regression lines on the two sides of each knot to be equal at the knot. We can do this by fitting the linear model
$$
y_i = \beta_0 + \beta_1 x_{i1} + \beta_2 x_{i2} + \beta_3 x_{i3} + \varepsilon_i
$$
where $x_{i1} = x_i$,
$$
x_{i2} = \left\{  {\begin{array}{ll}
   0 & \mathrm{for} \: x_i \le t_1 \\
   x_i - t_1 & \mathrm{for} \: x_i > t_1 \\
  \end{array} } \right.
$$
and
$$
x_{i3} = \left\{  {\begin{array}{ll}
   0 & \mathrm{for} \: x_i \le t_2 \\
   x_i - t_2 & \mathrm{for} \: x_i > t_2 \\
  \end{array} } \right.
$$
Notice that the regression-spline coefficients each have a straightforward interpretation: $\beta_0$ is the expectation of the response $y$ when $x = 0$; $\beta_1$ is the slope of the regression line in the first interval; $\beta_2$ is the _change_ in slope between the first and second interval; and $\beta_3$ is the change in slope between the second and third interval.

Fitting the linear regression spline to our artificial data produces the following result:
```{r}
x1 <- x # define spline regressors
x2 <- (x - t[1]) * (x > t[1])
x3 <- (x - t[2]) * (x > t[2])
mod3 <- lm(y ~ x1 + x2 + x3) # linear regression spline model
plot(x, y, main="Linear Regression Spline")
curve(f, add=TRUE, lty=2, lwd=2)
abline(v=t, lty=3)
mtext(c(expression(t[1]), expression(t[2])), side=1, line=0.5, at=t)
x01 <- x0
x02 <- (x0 - t[1]) * (x0 > t[1])
x03 <- (x0 - t[2]) * (x0 > t[2])
lines(x0, predict(mod3, list(x1=x01, x2=x02, x3=x03)), lwd=2) # draw fitted lines
yhat3 <- fitted(mod3)
mean((Ey - yhat3)^2) # RMSE
```

Because of the 2 constraints that the regression lines join at the 2 knots, the linear regression spline uses 2 fewer parameters than the unconstrained within-interval linear regressions: 4 parameters for the linear regression spline, $3 \times 2 = 6$ parameters for the unconstrained linear regressions.

This approach generalizes readily to higher-degree polyomials, the most common of which is the third-degree or _cubic regression spline_. Consider first  unconstrained within-interval cubic regressions:
```{r}
mods4 <- by(data.frame(x=x, y=y, group=group), group, # within-interval cubic regressions
    function(df) lm(y ~ poly(x, degree=3, raw=TRUE), data=df))
plot(x, y, main="Piecewise Cubic Fit")
curve(f, add=TRUE, lty=2, lwd=2)
abline(v=t, lty=3)
mtext(c(expression(t[1]), expression(t[2])), side=1, line=0.5, at=t)
xx <- seq(0, 10, length=500)
for (j in 1:3){
    xj <- xx[xx >= x0[j] & xx < x0[j + 1]]
    yj <- predict(mods4[[j]], list(x=xj))
    lines(xj, yj, lwd=2)  # draw within-interval cubics
}
yhat4 <- unlist(lapply(mods4, fitted))
mean((Ey - yhat4)^2) # RMSE
```

This model uses $4 \times 3 = 12$ parameters.

As before, we can constrain the cubic polynomials to join at the knots, generating the model
\begin{align*}
y_i = \beta_0 &+ \beta_{11} x_{i1} + \beta_{12} x_{i1}^2 + \beta_{13} x_{i1}^3 \\
              &+ \beta_{21} x_{i2} + \beta_{22} x_{i2}^2 + \beta_{23} x_{i2}^3 \\
              &+ \beta_{31} x_{i3} + \beta_{32} x_{i3}^2 + \beta_{33} x_{i2}^3 + \varepsilon_i
\end{align*}
```{r}
mod5 <- lm(y ~ poly(x1, degree=3, raw=TRUE) + 
               poly(x2, degree=3, raw=TRUE) + 
               poly(x3, degree=3, raw=TRUE))
plot(x, y, main="Cubic Regression Spline with Order-0 Smoothness")
curve(f, add=TRUE, lty=2, lwd=2)
abline(v=t, lty=3)
mtext(c(expression(t[1]), expression(t[2])), side=1, line=0.5, at=t)
x01 <- xx
x02 <- (xx - t[1]) * (xx > t[1])
x03 <- (xx - t[2]) * (xx > t[2])
lines(xx, predict(mod5, list(x1=x01, x2=x02, x3=x03)), lwd=2)
yhat5 <- fitted(mod5)
mean((Ey - yhat5)^2) # RMSE
```

Because this model uses only 1 intercept, it has 2 fewer parameters than the unconstrained within-interval cubic regressions. The cubic spline is continuous but, as close inspection reveals, it isn't smooth at the knots, where the slope of the fitted regression is free to change. We say that it has _order-0_ smoothness. 

If we eliminate the linear terms $\beta_{21} x_{i2}$ and $\beta_{31} x_{i3}$ from the model we force the slope---that is, the first derivative---of the regression function to be equal on both sides of each knot, producing _order-1 smoothness_:
\begin{align*}
y_i = \beta_0 + \beta_{11} x_{i1} &+ \beta_{12} x_{i1}^2 + \beta_{13} x_{i1}^3 \\
              &+ \beta_{22} x_{i2}^2 + \beta_{23} x_{i2}^3 \\
              &+ \beta_{32} x_{i3}^2 + \beta_{33} x_{i2}^3 + \varepsilon_i
\end{align*}

```{r}
mod6 <- lm(y ~ poly(x1, degree=3, raw=TRUE) + I(x2^2) + I(x3^2) + I(x2^3) + I(x3^3))
plot(x, y, main="Cubic Regression Spline with Order-1 Smoothness")
curve(f, add=TRUE, lty=2, lwd=2)
abline(v=t, lty=3)
mtext(c(expression(t[1]), expression(t[2])), side=1, line=0.5, at=t)
lines(xx, predict(mod6, list(x1=x01, x2=x02, x3=x03)), lwd=2)
yhat6 <- fitted(mod6)
mean((Ey - yhat6)^2) # RMSE
```

Finally, removing both the linear and the quadratic terms $\beta_{22} x_{i2}^2$ and $\beta_{32} x_{i3}^2$ forces the slope (first derivative) and curvature (second derivative) to be equal on both sides of each knot, producing _order-2 smoothess_ and, except for parametrization (discussed below), a traditional cubic regression spline:
$$
y_i = \beta_0 + \beta_{11} x_{i1} + \beta_{12} x_{i1}^2 + \beta_{13} x_{i1}^3 + \beta_{23} x_{i2}^3 + \beta_{33} x_{i3}^3 + \varepsilon_i
$$
This model uses $4 + 2 = 7$, or more generally $4 + k$, parameters, considerably fewer than the $4 \times 3 = 12$, or more generally $4 \times (k + 1)$ parameters used by the unconstained piecewise cubic model. 
```{r}
mod7 <- lm(y ~ poly(x1, degree=3, raw=TRUE) + I(x2^3) + I(x3^3))
plot(x, y, main="Cubic Regression Spline with Order-2 Smoothness")
curve(f, add=TRUE, lty=2, lwd=2)
abline(v=t, lty=3)
mtext(c(expression(t[1]), expression(t[2])), side=1, line=0.5, at=t)
lines(xx, predict(mod7, list(x1=x01, x2=x02, x3=x03)), lwd=2)
yhat7 <- fitted(mod7)
mean((Ey - yhat7)^2) # RMSE
```

## 3. Using `gspline()` for Fitting Generalized Regression Spline Models

The `gspline()` (*generalized regression spline*) function in the **gspline** package can fit all of the piecewise polynomial and regression spline models in the preceding section, and has substantial flexibility beyond what we have thus far discussed. We'll focus initially on three arguments to `gspline()` that allow us to reproduce the results of the previous section:

* `knots`: a vector giving the location(s) of the knots; this argument has no default value.

* `degree`: the degree of the regression spline, defaulting to `3` (a cubic spline); `degree` is often a single number but it can also be a vector of length equal to the number of knots plus 1, corresponding to the intervals between the extremes and knots. Thus, for example, if we have 2 knots at `knots = c(3/10, 6/10)` generating 3 intervals, specifying `degree = c(2, 3, 2)` would fit  (possibly constrained) quadratics in the first and last interval and a cubic in the middle interval (a model that `gspline()` can accommodate but for which it isn't straightforward to code regressors).

* `smooth`: the order of smoothness at the knots, specifying the number of derivatives to be matched on both sides of a knot. Thus `smooth = 2` implies that the fitted regression is continuous at a knot and has equal first (slope) and second (curvature) derivatives on both sides of the knot. The value `smooth = 0`, therefore, insures continuity but nothing more, and the special value `smooth = -1` permits discontinuity, as in a piecewise fit. The default for `smooth` is 1 less than `degree`; thus, for example, for a cubic regression spline the default is `smooth = 2`, matching first and second derivatives. If `degree` is a vector, then the default for `smooth` is 1 less than the minimum degree, but `smooth` may also be a vector of length equal to the number of knots, specifying the number of derivatives to be matched at each knot.

`gspline()` returns a function (or, more technically, a _closure_) that has several arguments, the first of which, `x`, is normally the predictor vector for which we wish to construct a regression spline, and the only argument that we'll describe here. The function returned by `gspline()` can also be used to generate hypothesis matrices for testing characteristics of the resulting fitted model, as we explain in another vignette.

Here are the examples from the preceding section, with checks to verify that we obtain identical results (within rounding error); to reduce redundancy, we've graphed the first example, but not the others:

### 3.1 Piecewise Constant Fit (Discontinuous)
```{r}
library(gspline)
sp_pc <- gspline(knots=t, degree=0, smoothness=-1)
sp_pc(0:10)
mod.gsp1 <- lm(y ~ sp_pc(x))
all.equal(as.vector(yhat1), as.vector(fitted(mod.gsp1)))
plot(x, y, main="Piecewise Constant Fit Using gspline()")
curve(f, add=TRUE, lty=2, lwd=2)
abline(v=t, lty=3)
mtext(c(expression(t[1]), expression(t[2])), side=1, line=0.5, at=t) 
lines(seq(0, 10, length=1000), 
      predict(mod.gsp1, newdata=list(x=seq(0, 10, length=1000))),
      lwd=2)
```

### 3.2 Piecewise Linear Fit (Discontinuous)
```{r}
sp_pl <- gspline(knots=t, degree=1, smoothness=-1)
sp_pl(0:10)
mod.gsp2 <- lm(y ~ sp_pl(x))
all.equal(as.vector(yhat2), as.vector(fitted(mod.gsp2)))
```

### 3.3 Linear Regression Spline (With Continuity at the Knots)
```{r}
sp_l <- gspline(knots=t, degree=1, smoothness=0)
sp_l(0:10)
mod.gsp3 <- lm(y ~ sp_l(x))
all.equal(as.vector(yhat3), as.vector(fitted(mod.gsp3)))
```

### 3.4 Piecewise Cubic Fit (Discontinuous)
```{r}
sp_pc <- gspline(knots=t, degree=3, smoothness=-1)
sp_pc(0:10)
mod.gsp4 <- lm(y ~ sp_pc(x))
all.equal(as.vector(yhat4), as.vector(fitted(mod.gsp4)))
```

### 3.5 Cubic Regression Spline with Smoothness 0 (With Continuity Only)

```{r}
sp_c0 <- gspline(knots=t, degree=3, smoothness=0)
sp_c0(0:10)
mod.gsp5 <- lm(y ~ sp_c0(x))
all.equal(as.vector(yhat5), as.vector(fitted(mod.gsp5)))
```


### 3.6 Cubic Regression Spline with Smoothness 1 (Matching Slope)
```{r}
sp_c1 <- gspline(knots=t, degree=3, smoothness=1)
sp_c1(0:10)
mod.gsp6 <- lm(y ~ sp_c1(x))
all.equal(as.vector(yhat6), as.vector(fitted(mod.gsp6)))
```

### 3.7 Cubic Regression Spline with Smoothness 2 (Matching Slope and Curvature)
```{r}
sp_c <- gspline(knots=t, degree=3, smoothness=2)
sp_c(0:10)
mod.gsp7 <- lm(y ~ sp_c(x))
all.equal(as.vector(yhat7), as.vector(fitted(mod.gsp7)))
```

## 4. Fitting Regression Splines With `bs()`

The `bs()` function in the standard R **splines** package, for fitting B-splines, takes several arguments, some of which are redundant (in the sense of being alternatives). As for `gspline()` and the functions that it generates, we discuss here only the arguments used in the examples:

* `x`: the predictor vector for which we wish to construct a regression spline, and the first argument to `bs()`.

* `knots`: a vector giving the location(s) of the knots; if it is not specified (and if the `df` argument, for degrees of freedom and not discussed here, is also not specified), then a global polynomial regression is fit.

* `degree`: the degree of the regression spline, defaulting to `3` (a cubic spline); unlike for `gspline()`, this is a single positive integer > 0. Also unlike `gspline()`, the order of smoothness at the knots isn't controllable and is always `degree` - 1.

Only 2 of the models considered above can be fit as B-splines using `bs()`; for brevity, we graph only the first of these  models but check in both cases that the fitted values are the same as computed by `gspline()`:

### 4.1 Linear Regression Spline (With Continuity at the Knots)
```{r}
library(splines)
mod.bs3 <- lm(y ~ bs(x, knots=t, degree=1))
all.equal(as.vector(fitted(mod.gsp3)), as.vector(fitted(mod.bs3)))
plot(x, y, main="Linear Regression Spline Using bs()")
curve(f, add=TRUE, lty=2, lwd=2)
abline(v=t, lty=3)
mtext(c(expression(t[1]), expression(t[2])), side=1, line=0.5, at=t)
lines(seq(0, 10, length=1000), 
      predict(mod.bs3, newdata=list(x=seq(min(x), max(x), length=1000))),
      lwd=2)
```


### 4.2 Cubic Regression Spline with Smoothness 2 (Matching Slope and Curvature)
```{r}
mod.bs7 <- lm(y ~ bs(x, knots=t, degree=3))
all.equal(as.vector(fitted(mod.gsp7)), as.vector(fitted(mod.bs7)))
```

## 5. How `gspline()` and `bs()` Construct Bases for Fitting Regression Splines

We've demonstrated how (1) directly coding spline regressors, (2) using `gspline()`, and (3) using `bs()` construct equivalent regression-spline bases that produce, within rounding error, identical fits to the data for cases, such as a standard cubic regressions spline, to which all three apply. In this section, we look more closely at the spline bases constructed directly and by `gspline()` and `bs()`.

We return to the cubic regression spline with two knots at $t_1$ and $t_2$ discussed in Section 2:
$$
y_i = \beta_0 + \beta_1 x_{i1} + \beta_2 x_{i2} + \beta_3 x_{i3} + \varepsilon_i
$$
where $x_{i1} = x_i$,
$$
x_{i2} = \left\{  {\begin{array}{ll}
   0 & \mathrm{for} \: x_i \le t_1 \\
   x_i - t_1 & \mathrm{for} \: x_i > t_1 \\
  \end{array} } \right.
$$
and
$$
x_{i3} = \left\{  {\begin{array}{ll}
   0 & \mathrm{for} \: x_i \le t_2 \\
   x_i - t_2 & \mathrm{for} \: x_i > t_2 \\
  \end{array} } \right.
$$

Here are the bases constructed directly, by `gspline()`, and by `bs()`, respectively named `X`, `X.gsp`, and `X.bs`, for a simple $x$ vector with the integers from 0 to 10:
```{r}
x <- 0:10
t <- c(10/3, 20/3)
x1 <- x # define spline regressors
x2 <- (x - t[1]) * (x > t[1])
x3 <- (x - t[2]) * (x > t[2])
X <- cbind(x1, x1^2, x1^3, x2^3, x3^3)
colnames(X) <- c("x1", "x1^2", "x1^3", "x2^3", "x3^3")
round(X, 5)

sp <- gspline(knots=t, degree=3, smoothness=2)
X.gsp <- sp(x)
X.gsp

X.bs <- bs(x, knots=t, degree=3)
X.bs

```
We'll discuss _boundary knots_ in the next section when we take up natural splines. The other attributes associated with the matrix returned by `bs()` are self-explanatory; for further information see `?bs`.

The coding in `X` is entirely straightforward, directly reflecting the definitions of the variables $x_1$, $x_2$, and $x_3$ above. Comparison of `X` and `X.gsp` suggests that the columns of the latter are multiples of those of the former as we can easily verify:
```{r}
X / X.gsp
```
The `NaN`s are produced by dividing 0 by 0 and so are innocuous here. We can produce `X.gsp` from `X` by post-multiplying the latter by the diagonal weight matrix $\mathrm{diag}(1, 1/2, 1/6, 1/6, 1/6)$:
```{r}
round(X %*% diag(c(1, 1/2, 1/6, 1/6, 1/6)), 5)
```
We explain below why `gspline()` uses this scaling of the columns.

Demonstrating the equivalence of the bases produces by `bs()` and `gspline()` is a bit harder because the coefficients for the individual columns of the B-spline basis aren't straightforwardly interpretable (which, by the way, is the essential advantage of using `gspline()`), but an indirect way to do this is to project the columns of `X.bs` onto the column space of `X.gsp`; if the two span the same subspace, then the projection should be equal to `X.bs`. We can easily perform the projection by least-squares regression, extracting the projected columns as the "fitted values":
```{r}
all.equal(fitted(lm(X.bs ~ X.gsp - 1)), X.bs, check.attributes=FALSE)
```

To understand the coding used by `gspline()` let's rewrite the simple cubic-spline model with 2 knots in terms of $x$, omitting the subscript $i$ for cases and taking the expectation $\mu = E(y)$ of the response:
$$
\mu = \left\{  {\begin{array}{ll}
    \beta_0 + \beta_{11}x + \beta_{12}x^2 + \beta_{13}x^3 & \mathrm{for} \: x \le t_1 \\
    \beta_0 + \beta_{11}x + \beta_{12}x^2 + \beta_{13}x^3 + \beta_{23}(x - t_1)^3 & \mathrm{for} \: t_1 < x \le t_2 \\
    \beta_0 + \beta_{11}x + \beta_{12}x^2 + \beta_{13}x^3 + \beta_{23}(x - t_1)^3 + \beta_{33}(x - t_2)^3 & \mathrm{for} \: x > t_2 \\
  \end{array} } \right.
$$
Now differentiate $\mu$ repeatedly with respect to $x$ in each of the 3 intervals generated by the 2 knots:
$$
\begin{array}{ll}
\dfrac{d\mu}{dx} &= \left\{  {\begin{array}{ll}
    \beta_{11} + 2\beta_{12}x + 3\beta_{13}x^2 & \mathrm{for} \: x \le t_1 \\
    \beta_{11} + 2\beta_{12}x + 3\beta_{13}x^2 + 3\beta_{23}(x - t_1)^2 & \mathrm{for} \: t_1 < x \le t_2 \\
    \beta_{11} + 2\beta_{12}x + 3\beta_{13}x^2 + 3\beta_{23}(x - t_1)^2 + 3\beta_{33}(x - t_2)^2 & \mathrm{for} \: x > t_2 \\
  \end{array} } \right.\\
  \\
\dfrac{d^2\mu}{dx^2} &= \left\{  {\begin{array}{ll}
    2\beta_{12} + 6\beta_{13}x & \mathrm{for} \: x \le t_1 \\
    2\beta_{12} + 6\beta_{13}x + 6\beta_{23}(x - t_1) & \mathrm{for} \: t_1 < x \le t_2 \\
    2\beta_{12} + 6\beta_{13}x + 6\beta_{23}(x - t_1) + 6\beta_{33}(x - t_2) & \mathrm{for} \: x > t_2 \\
  \end{array} } \right.\\
  \\
\dfrac{d^3\mu}{dx^3} &= \left\{  {\begin{array}{ll}
    6\beta_{13} & \mathrm{for} \: x \le t_1 \\
    6\beta_{13} + 6\beta_{23} & \mathrm{for} \: t_1 < x \le t_2 \\
    6\beta_{13} + 6\beta_{23} + 6\beta_{33} & \mathrm{for} \: x > t_2 \\
  \end{array} } \right.\\
\end{array}
$$
The change in the third derivative at each knot is consequently 6 times the coefficient of the corresponding cubic regressor:
$$
\begin{array}{ll}
\dfrac{d^3\mu}{dx^3}\biggr\rvert_{x = t_{1+}} - \dfrac{d^3\mu}{dx^3}\biggr\rvert_{x = t_{1}} &= 6\beta_{23}\\
\dfrac{d^3\mu}{dx^3}\biggr\rvert_{x = t_{2+}} - \dfrac{d^3\mu}{dx^3}\biggr\rvert_{x = t_{2}} &= 6\beta_{33}\\
\end{array}
$$
where $\dfrac{d^3\mu}{dx^3}\biggr\rvert_{x = t_{j+}}$ is the third derivative evaluated just to the right of knot $t_j$. 

Moreover, had the spline been a _quadratic_ rather than a cubic regression spline, then the differences in the _second_ derivative at the knots would have been 2 times the corresponding coefficients of the quadratic spline regressors, $2\beta_{22}$ and $2\beta_{23}$; and, finally, had the spline been a _linear_ regression spline, then the differences in the _first_ derivative at the knots would simply have been equal to the corresponding coefficients of the linear spline regressors, $\beta_{12}$ and $\beta_{13}$. 

These observations explain the coding used by `gspline()`, which therefore produces coefficients that directly give differences in the highest-order derivative at the knots (e.g., of the third derivative for a cubic regression spline). This result is reflected in the names of the columns of the model matrix generated by `gspline()`. In the example, `D1|0`, `D2|0`, and `D3|0` are respectively the first, second, and third derivatives ("`D`") at the left boundary (i.e., $x = 0$);  `C3|3.33` is the change ("`C`") in the third derivative (symbolized by "`3`") at the first interior knot ($t_1 = 10/3 \approx 3.33$), and `C3|6.67` is the change in the third derivative at the second interior knot ($t_1 = 20/3 \approx 6.67$).

### 5.1 Numerical Conditioning

Although the spline regression bases constructed directly by `gspline()` and by `bs()` are equivalent in the sense that they span the same subspace, these bases have different numerical properties. One way to see this is by computing the correlations among the regressors in each case:
```{r}
round(cor(X), 3)
round(cor(X.gsp), 3)
round(cor(X.bs), 3)
```
The correlations for the first two bases are the same because, as we just explained, the columns of `X.gsp` are multiples of those of `X`, and some are quite large, while the correlations among the columns of `X.bs` are much smaller.

Even more directly, checking the condition numbers of the three basis matrices, each augmented with a column of 1s for the regression intercept, shows that both `X` and `X.gsp` are much more ill-conditioned than `X.bs`:
```{r}
kappa(cbind(1, X))
kappa(cbind(1, X.gsp))
kappa(cbind(1, X.bs))
```
Notice that the condition numbers of the augmented `X` and `X.gsp` matrices are _not_ identical because of the different scalings of the columns.

These observations suggest that if we're interested in regression splines simply for fitting smooth nonlinear regressions, and if the model is expressible as a B-spline, then there's an advantage in using the B-spline basis in preference to the basis constructed by `gspline()`. As we mentioned, the advantages of `gspline()` are in supporting a wider variety of piecewise polynomial and regression spline models and, crucially, in producing interpretable coefficients and hypothesis tests. We illustrate uses of the latter in another vignette.

The `gspline()` function provides a mechanism for generating numerically stable general regression spline bases while maintaining interpretable coefficients---allowing us to have our cake and eat it too. It does so by scaling the predictor $x$ to the range $-10$ to $+10$; orthonormalizing the resulting interpretable spline basis (as described above); and transforming the spline regressors to orthogonality with the constant regressor (i.e., for the intercept). `gspline()` also saves the information necessary to undo these transformations. The `wald()` function uses this information to produce naturally specified and interpretable estimates and tests. The details are given in another vignette, but we develop simple illustrations here.

`gspline()` needs to know the values of the predictor $x$ in order to compute a numerically stable basis for the spline regressors generated from $x$. It does this through the initial (and optional) argument `x`; if `x` is specified then the resulting spline-generating function returns a scaled and orthonormal basis; if `x` is absent, then `gspline()` behaves as in the preceding examples in this vignette.

Providing a bit more detail, the creation of a stable spline basis is controlled by three arguments:

1. `stable`: if `TRUE` (the default if the `x` argument is specified and inapplicable if it is not), transform the spline basis into a equivalent orthonormal basis after possible rescaling of `x`.

2. `rescale`: if `TRUE` (the default is the value of `stable`), rescale `x` to the range between $-10$ and $10$.

3. `ortho2intercept`: if `TRUE` (the default is the value of `stable`), the spline basis is replaced with a basis for the orthogonal complement to the 1-vector in the space spanned by the 1-vector and the spline basis. 

The defaults for these arguments are chosen to provide a numerically stable basis in the most common circumstances, but the `ortho2intercept` argument in particular must be set carefully: If the model violates the principle of marginality (for example by failing to include the constant regressor), then the resulting basis for the regression spline will not be constructed correctly.

Consider the earlier examples in this section, for which
```{r}
x
t
```
Then
```{r}
sp.stable <- gspline(x, knots=t, degree=3, smoothness=2)
X.gspst <- sp.stable(x)
X.gspst
all.equal(fitted(lm(cbind(1, X.gsp) ~ X.gspst)), 
          cbind(1, X.gsp), 
          check.attributes=FALSE) # span the same space
zapsmall(crossprod(X.gspst)) # orthonormal
round(cor(X.gspst), 3) # uncorrelated
kappa(cbind(1, X.gspst)) # better conditioned
```

The rows and columns of the matrix returned by `sp.stable()` are unnamed, because they no longer have the straightforward interpretations that they had in the 'raw' spline bases considered previously. As mentioned, however, `wald()` can reconstruct the parameters implied by the raw parametrization. To see how this works, let's return to the artificial data generated in Section 2 of the vignette (and regenerated here):
```{r}
set.seed(12345)
x <- runif(200, 0, 10)
x <- sort(x)
Ey <- cos(1.25*(x + 1)) + x/5
y <- Ey + 0.5*rnorm(200)
```

Now let's use `gspline()` to fit two equivalent models, one with raw and the other with stable cubic regression splines:
```{r}
sp.raw <- gspline(knots=c(10/3, 20/3))
sp.st <- gspline(x, knots=c(10/3, 20/3))
mod.raw <- lm(y ~ sp.raw(x))
mod.st <- lm(y ~ sp.st(x))
all.equal(fitted(mod.raw), fitted(mod.st))
kappa(model.matrix(mod.raw))
kappa(model.matrix(mod.st))
```
The fits to the data are equivalent, but the model matrix for the stable fit is *much* better conditioned. 

Of course, the coefficients for the two fits differ---those for the raw fit are directly interpretable, while those for the numerically stable fit are not: 
```{r}
coef(mod.raw)
coef(mod.st)
```
Yet `wald()` can recover the raw coefficients from the numerically stable fit:
```{r}
wald(mod.raw)
wald(mod.st) # identical
```
The first line of output in each case shows the omnibus $F$-test for the intercept and the 5 regression spline coefficients, and the subsequent table gives identical estimates, standard errors, tests, and confidence limits for the individual interpretable coefficients. By default, the `print()` method for `"wald"` objects shows "adjusted" as well as "raw" two-sided $p$-values for the individual coefficient tests. The adjustment is made using the standard R function `p.adjust()` with (again by default) the Holm method. In this example, the raw and adjusted $p$-values are indistinguishable because the very small $p$-values for the spline coefficients are all rounded to $<.00001$.

## Natural Splines

Natural splines are B-splines with additional _boundary knots_, typically placed at the extremes of $x$, and the additional constraints that the fitted regression is linear to the left of the lower boundary knot and to the right of the upper boundary knots. Derivatives up to order 1 less than the degree of the interior spline polynomials also match at the boundaries. Natural splines in effect place 2 additional constraints on the model and therefore are representable using 2 fewer parameters than corresponding B-splines: A cubic natural regression spline with $k$ knots uses $1 + k$ parameters (and regressors), excluding the intecept, 2 fewer than the corresponding B-spline.

Because of the condition that the regression is linear beyond the boundaries, natural splines tend to produce less wild extrapolations than B-splines. Natural splines may also behave more reasonably in the interior near the boundaries, but are less flexible than B-splines. We can compensate for the lower flexibility of natural splines within the range of the data by introducing 1 or 2 additional interior knots, perhaps spending the degrees of freedom more wisely than for a B-spline with the same number of parameters. The `ns()` function in the standard R **splines** package only fits natural _cubic_ regression splines.

To illustrate, we continue the example from the preceding section, with $x = (0, 1, 2, \ldots, 10)$, $k = 2$ interior knots at 10/3 and 20/3, and boundary knots at 0 and 10:
```{r}
x <- 0:10
X.ns <- ns(x, knots=c(10/3, 20/3), Boundary.knots=c(0, 10))
round(X.ns, 5)
```
The `Boundary.knots` argument to `ns()` defaults to the range of `x`, and so we didn't have to specify it explicitly but have done so for clarity. The `knots` argument gives the interior knots.

We can specify an equivalent basis using `gspline()` by specifying linear (i.e., degree-1) polynomials beyond the boundaries along with interior cubics, and by matching derivatives up to order-2:
```{r}
sp_ns <- gspline(knots=c(0, 10/3, 20/3, 10), 
             degree=c(1, 3, 3, 3, 1), 
             smoothness=c(2, 2, 2, 2))
X.gsp.ns <- sp_ns(x)
X.gsp.ns

all.equal(fitted(lm(X.ns ~ X.gsp - 1)), X.ns, check.attributes=FALSE)
```
Notice that `gspline()` doesn't distinguish between interior and boundary knots, and that the distinction is imposed implicitly by the distribution of the $x$-values and the specification of the `degree` and `smoothness` arguments. Because all of the values of `smoothness` are equal to 2, we could have specified `smoothness=2`, but we entered the argument as a 4-element vector to emphasize the number of constraints placed on the derivatives.

## References 
