## Function 'tracer' constructs a tracer object containing a 'trace' and 
## a 'get' function. The 'trace' function can print object values from its
## calling environment (the parent frame) when evaluated, and it can also 
## store those object values in a local list. The time between 'trace' calls
## can be recorded as well. It is measured by the 'hires_time' function from the 
## bench package. The 'get' function can subsequently access the traced values
## in the list, though the user will typically do so via the S3 methods 'print' 
## or 'summary' below.
## 
## The call of the 'trace' function can be manually inserted into the body of 
## any function, it can be inserted using 'base::trace', or it can be passed
## as an argument to any function with a callback argument. 
##
## Arguments of 'tracer' are:
##
## object: a character vector of names of the objects in the calling 
##         environment of the 'trace' function that are to be traced. Objects 
##         created by the 'expr' argument can be traced. 
## N:      an integer specifying if and how often trace information is printed. 
##         N = 0 means never, and otherwise trace information is printed every 
##         N-th iteration. N = 1 is the default.
## save:   a logical value. Sets if the trace information is saved.
## time:   a logical value. Sets if run time information is saved.
## expr:   an expression that will be evaluated in the calling environment 
##         of the 'trace' function.
## ...:    other arguments passed to `format` for printing.

#install.packages("tracer")
#install.packages("bench")
library(microbenchmark)

library(bench)
tracer <- function(objects = NULL, N = 1, save = TRUE, time = TRUE, expr = NULL, ...) {
  n <- 1
  values_save <- list()
  last_time <- bench::hires_time()
  trace <- function() {
    time_diff <- bench::hires_time() - last_time
    if(is.expression(expr))
      eval(expr, envir = parent.frame())
    if(is.null(objects))
      objects <- ls(parent.frame())
    
    values <- mget(objects, envir = parent.frame(), ifnotfound = list(NA))
    if(N && (n == 1 || n %% N == 0))
      cat("n = ", n, ": ",  paste(names(values), " = ", format(values, ...), 
                                  "; ", sep = ""), "\n", sep = "")
    if(save) {
      if(time)
        values[[".time"]] <- time_diff
      values_save[[n]] <<- values
    }
    n <<- n + 1
    last_time <<- bench::hires_time()
    invisible(NULL)
  }
  get <- function(simplify = FALSE) {
    if(simplify) {
      col_names <- unique(unlist(lapply(values_save, names)))
      values_save <- lapply(
        col_names, 
        function(x) do.call(rbind, unlist(lapply(values_save, function(y) y[x]), 
                                          recursive = FALSE))
      )
      names(values_save) <- col_names
      values_save <- lapply(col_names, function(x) {
        x_val <- values_save[[x]] 
        if(!is.null(ncol(x_val)) && ncol(x_val) == 1) {
          colnames(x_val) <- x
        } else {
          if(is.null(colnames(x_val)))
            colnames(x_val) <- 1:ncol(x_val)
          colnames(x_val) <- paste(x, ".", colnames(x_val), sep = "")
        }
        x_val
      })
      values_save <- do.call(cbind, values_save)
      row.names(values_save) <- 1:nrow(values_save)
    }
    values_save
  }
  structure(list(trace = trace, get = get), class = "tracer")
}


## Methods for subsetting, printing and summarizing tracer objects
'[.tracer' <- function(x, i, j, ..., drop = TRUE) {
  values <- x$get(...)[i] 
  if (drop && length(i) == 1)
    values <- values[[1]]
  values
}
print.tracer <- function(x, ...) print(x$get(...))
summary.tracer <- function(x, ...) {
  x <- suppressWarnings(x$get(simplify = TRUE))
  x[, ".time"] <- c(0, cumsum(x[-1, ".time"]))
  as.data.frame(x)
}

