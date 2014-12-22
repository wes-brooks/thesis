ecdf.w <- function (x, wt) 
{
    ord <- order(x)
    n <- length(x)
    sum.w = sum(wt)

    x <- x[ord]
    wt <- wt[ord]

    if (n < 1) 
        stop("'x' must have 1 or more non-missing values")

    cum <- cbind(x, wt) %>% data.frame %>% group_by(x) %>% summarize(w=sum(wt))
    cum$cum.wt = cumsum(cum$w) / sum(cum$w)
    
    rval <- approxfun(cum$x, cum$cum.wt, method = "constant", yleft = 0, yright = 1, f = 0, ties = "ordered")
    class(rval) <- c("ecdf", "stepfun", class(rval))
    assign("nobs", n, envir = environment(rval))
    assign("weights", w, envir = environment(rval))
    attr(rval, "call") <- sys.call()
    rval
}