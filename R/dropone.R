#' Drop each row or cluster from a fitted object
#' 
#' Drop each row or each cluster in a fitted object and 
#' collect the resulting linear estimates. Return 
#' the drop one estimates along with variable from original
#' data that are invariant within clusters, as well as the maximum
#' value of DFBETAS for each cluster dropped.
#' 
#' @param fit a fitted object with an \code{\link{update}} method
#'        and a \code{\link{getData}} method.
#' @param form a formula evaluated in the data frame for 'fit'
#'        that defines clusters to be dropped one at a time. If NULL, the 
#'        default, each row is dropped successively.
#' @param FUN function used to extract coefficients. Default:
#'        \code{\link[spida2]{fixef}} for class 'lme' and 
#'        \code{\link{coef}} otherwise. 
#' @param data  data frame in which to evaluate the refitted object. Default:
#'        'getData(fit)'.
#' @export
dropone <-
function (fit, form = NULL, FUN = if (inherits(fit, "lme")) fixef else coef, 
    data = getData(fit), ...) 
{
    if (is.null(form)) {
        by <- factor(1:nrow(data))
        data$by <- by
        dframe <- data
    }
    else {
        by <- model.frame(form, data)
        by <- do.call(paste, c(by, sep = "/"))
        data$by <- factor(by)
        dframe <- up(data, form)
    }
    values <- dframe$by
    names(values) <- values
    ret <- lapply(values, function(v) {
        ret <- try(update(fit, data = data[by != v, , drop = FALSE], 
            ...))
        if (inherits(ret, "try-error")) 
            NA
        else FUN(ret)
    })
    ret <- do.call(rbind, ret)
    max_dfbetas <- apply(scale(ret), 1, function(x) max(abs(x), na.rm = T))
    colnames(ret) <- paste0("b_", colnames(ret))
    ret <- cbind(ret, dframe)
    ret$max_dfbetas <- max_dfbetas
    ret
}
#' @describeIn dropone parallel computing version of dropone using the parallel package
#' @param mc.cores number of cores to use. Default: \code{parallel::detectCores()} \eqn{\exp^n}
#' @export
parDropone <-
function (fit, form = NULL, FUN = if (inherits(fit, "lme")) fixef else coef, 
    data = getData(fit), mc.cores = detectCores(), ...) 
{
    library(parallel)
    if (is.null(form)) {
        by <- factor(1:nrow(data))
        data$by <- by
        dframe <- data
    }
    else {
        by <- model.frame(form, data)
        by <- do.call(paste, c(by, sep = "/"))
        data$by <- factor(by)
        dframe <- up(data, form)
    }
    values <- dframe$by
    names(values) <- values
    ret <- mclapply(values, function(v) {
        ret <- try(update(fit, data = data[by != v, , drop = FALSE], 
            ...))
        if (inherits(ret, "try-error")) 
            NA
        else FUN(ret)
    }, mc.cores = mc.cores)
    ret <- do.call(rbind, ret)
    max_dfbetas <- apply(scale(ret), 1, function(x) max(abs(x), na.rm = T))
    colnames(ret) <- paste0("b_", colnames(ret))
    ret <- cbind(ret, dframe)
    ret$max_dfbetas <- max_dfbetas
    ret
}
