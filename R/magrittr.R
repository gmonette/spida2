###
### Extra pipes for magrittr
###

#' Pipes for magrittr that passes LHS invisibly
#'
#' Also the code for \code{%>%} copied from magrittr (and/or dplyr)
#' 
#' @param lhs
#' @param rhs
#' @name pipes
#' @examples
#' data(hs)
#' hs %T% dim  %T% names  %>% head 
#' @export
# "%>%" <-
# function (lhs, rhs) 
# {
#   if(!require(magrittr)) stop("Install the magrittr package")
#   parent <- parent.frame()
#   env <- new.env(parent = parent)
#   chain_parts <- split_chain(match.call(), env = env)
#   pipes <- chain_parts[["pipes"]]
#   rhss <- chain_parts[["rhss"]]
#   lhs <- chain_parts[["lhs"]]
#   env[["_function_list"]] <- lapply(1:length(rhss), function(i) wrap_function(rhss[[i]], 
#                                                                               pipes[[i]], parent))
#   env[["_fseq"]] <- `class<-`(eval(quote(function(value) freduce(value, 
#                                                                  `_function_list`)), env, env), c("fseq", "function"))
#   env[["freduce"]] <- freduce
#   if (is_placeholder(lhs)) {
#     env[["_fseq"]]
#   }
#   else {
#     env[["_lhs"]] <- eval(lhs, parent, parent)
#     result <- withVisible(eval(quote(`_fseq`(`_lhs`)), env, 
#                                env))
#     if (is_compound_pipe(pipes[[1L]])) {
#       eval(call("<-", lhs, result[["value"]]), parent, 
#            parent)
#     }
#     else {
#       if (result[["visible"]]) 
#         result[["value"]]
#       else invisible(result[["value"]])
#     }
#   }
# }
# alternative if we're going to use 'magrittr' anyways:
"%|%" <-
function (lhs, rhs) {
    if(!require(magrittr)) stop("Install the magrittr package")
    lhs %>% rhs
}

`%T%` <-
function (lhs, rhs) 
{
  library(magrittr)
  parent <- parent.frame()
  env <- new.env(parent = parent)
  chain_parts <- magrittr:::split_chain(match.call(), env = env)
  pipes <- chain_parts[["pipes"]]
  rhss <- chain_parts[["rhss"]]
  lhs <- chain_parts[["lhs"]]
  disp(length(rhss))
  disp(length(pipes))
  env[["_function_list"]] <- lapply(seq_along(rhss), function(i) magrittr:::wrap_function(rhss[[i]], 
                                                                              pipes[[i]], parent))
  env[["_fseq"]] <- `class<-`(eval(quote(function(value) magrittr:::freduce(value, 
                                                                 `_function_list`)), env, env), c("fseq", "function"))
  env[["freduce"]] <- freduce
  if (magrittr:::is_placeholder(lhs)) {
    env[["_fseq"]]
  }
  else {
    env[["_lhs"]] <- eval(lhs, parent, parent)
    result <- invisible(eval(quote(`_fseq`(`_lhs`)), env, 
                               env))
    if (magrittr:::is_compound_pipe(pipes[[1L]])) {
      invisible(eval(call("<-", lhs, result[["value"]]), parent, 
           parent))
    }
    else {
      if (result[["visible"]]) 
        invisible(result[["value"]])
      else invisible(result[["value"]])
    }
  }
}



`%T%` <- function(lhs,rhs) {
  library(magrittr)
  lhs.n <- deparse(substitute(lhs))
  rhs.n <- deparse(substitute(rhs))
  lhs <- lhs
  p <- (lhs  %>% rhs)
  print(p)
  cat("----------\n")
  invisible(lhs)
}


'%T%' <-
function (lhs, rhs) 
{
  parent <- parent.frame()
  env <- new.env(parent = parent)
  chain_parts <- split_chain(match.call(), env = env)
  pipes <- chain_parts[["pipes"]]
  rhss <- chain_parts[["rhss"]]
  lhs <- chain_parts[["lhs"]]
  env[["_function_list"]] <- lapply(1:length(rhss), function(i) wrap_function(rhss[[i]], 
                                                                              pipes[[i]], parent))
  env[["_fseq"]] <- `class<-`(eval(quote(function(value) freduce(value, 
                                                                 `_function_list`)), env, env), c("fseq", "function"))
  env[["freduce"]] <- freduce
  if (is_placeholder(lhs)) {
    env[["_fseq"]]
  }
  else {
    env[["_lhs"]] <- eval(lhs, parent, parent)
    result <- withVisible(eval(quote(`_fseq`(`_lhs`)), env, 
                               env))
    if (is_compound_pipe(pipes[[1L]])) {
      eval(call("<-", lhs, result[["value"]]), parent, 
           parent)
    }
    else {
      if (result[["visible"]]) 
        result[["value"]]
      else invisible(result[["value"]])
    }
  }
}

