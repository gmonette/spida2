#' Stack data frames vertically
#'
#' Stack data frames vertically on top of each other with columns identified by variable names. Matrices are coerced to data frames.
#'
#' @param \dots data frames or objects that can be coerced to data frames
#' @examples
#' zd1 <- data.frame(x = 1:3, y = 2, z = 3, .unique. = 'a')
#' zd2 <- data.frame(x = 1:4, y = 'y', w = 1, m = matrix(1:8, nrow = 4))
#' zm1 <- cbind(x = 1:3, y = 2, z = 3)
#' zm2 <- cbind(x = 1:4, y = 3, w = 1)
#' Rbind(zm1, zm2, zd1, zd1, zd2)
#' Rbind(Rbind(zm1, zm2), Rbind(zd1, zd2))
#' @export
Rbind <- function(..., verbose = FALSE, id = ".index") {
  x <- list(...)
  x <- lapply(x, as.data.frame)
  nams <- unique(unlist(sapply(x, names)))
  # disp(nams)
  while( any(id %in% nams)) { # keep adding '.' until name does not conflict
    id <- paste0(id,".")
  } 
  if(verbose) print(id)
  x <- lapply(seq_along(x), function(ii) {
    x[[ii]][[id]] <- ii
    x[[ii]]
  })
  fun <- function(x,y) merge(x, y, all = TRUE, sort = FALSE)
  ret <- Reduce(fun,x)
  # ret[[nam]] <- NULL
  ret
}
#' Fix .x and .y variables in a merged data frame
#' 
#' Function to combine consistent .x and .y variable names after a merge
#' 
#' Some pairs will be identical and they should become one variable without the suffix.
#' Some will be compatible in the sense that non-missing values are identical and
#' values missing in one but not the other can be merged
#' Some will be incompatible and they should remain in their .x and .y version.
#' 
#' The idea is that after applying this function, any .x and .y variable pairs
#' are, in fact, distinct variables that need to be looked at before being used
#' in any analyses.
#' 
#' @param data  a data frame that is the result of a merge operation, possibly
#'              with .x and .y variables
#'              
#' Note that factors will merge into a common variable only if they have the same
#' set of values. To merge factors with different sets of values, consider
#' converting them to character strings.
#' @examples
#' da <- data.frame(
#'   id = 1:5, 
#'   g = c(1,2,NA,4,5), 
#'   h = factor(c('a','b','c','d','e')),
#'   hh = c('a','b','c','d','e'),
#'   m = factor(c('M',NA,NA,'F','M'))  # factor with same set of values
#' )
#' db <- data.frame(
#'   id = c(1:3,6), 
#'   g = c(1,2,NA,6), 
#'   h = factor(c('a','b','c','f')),  # factor with different sets of values
#'   hh = c('a','b','c','f'),
#'   k = c(1:3,6),
#'   m = factor(c('M','M',NA,'F'))
#' )
#' da  
#' db  
#' dm <- dplyr::full_join(da, db, by = 'id' )  
#' dm
#' fix_xy(dm)
#' @export              
fix_xy <- function(data) {
  # data is a data.frame, ususally the result of a merge operation such as
  # base::merge or dplyr::left_join, etc., which will create separate .x and
  # .y variable names for variables whose names appear in both source files
  # but are not among the 'by' variables for merging.
  # 
  # Note:
  # - factors will merge only if they have the same levels
  # 
  names.x <- sub('\\.x$','',grepv('\\.x$', names(data)))
  names.y <- sub('\\.y$','',grepv('\\.y$', names(data)))
  names.xy <- intersect(names.x, names.y)
  
  for (nn in names.xy){
    nx <- paste0(nn,'.x')
    ny <- paste0(nn,'.y')
    if(identical(data[[nx]], data[[ny]])){
      data[[nn]] <- data[[nx]]
      data[[nx]] <- NULL
      data[[ny]] <- NULL
    } else {
      both_valid <- (!is.na(data[[nx]])) & (!is.na(data[[ny]]))
      if(identical(data[[nx]][both_valid], data[[ny]][both_valid])) {
        new_var <- data[[nx]]
        new_var[is.na(new_var)] <- data[[ny]][is.na(new_var)] # copy from data[[ny]]
        data[[nn]] <- new_var
        data[[nx]] <- NULL
        data[[ny]] <- NULL
      }
    }
  }
  data
} 
#' 
#' Function to resolve .x and .y variable names after a merge
#' 
#' Some pairs will be identical and they should become one variable without the suffix.
#' Some will be compatible in the sense that non-missing values are identical and
#' values missing in one but not the other can be merged
#' Some will be incompatible and they should remain in their .x and .y version.
#' 
#' The idea is that after applying this function, any .x and .y variable pairs
#' are, in fact, distinct variables that need to be looked at before being used
#' in any analyses.
#' 
fix_xy <- function(data) {
  # data is a data.frame, usually the result of a merge operation such as
  # base::merge or dplyr::left_join, etc., which will create separate .x and
  # .y variable names for variables whose names appear in both source files
  # but are not among the 'by' variables for merging.
  # 
  # Note:
  # - factors will merge only if they have the same levels
  # 
  names.x <- sub('\\.x$','',grepv('\\.x$', names(data)))
  names.y <- sub('\\.y$','',grepv('\\.y$', names(data)))
  names.xy <- intersect(names.x, names.y)
  
  for (nn in names.xy){
    nx <- paste0(nn,'.x')
    ny <- paste0(nn,'.y')
    if(identical(data[[nx]], data[[ny]])){
      data[[nn]] <- data[[nx]]
      data[[nx]] <- NULL
      data[[ny]] <- NULL
    } else {
      both_valid <- (!is.na(data[[nx]])) & (!is.na(data[[ny]]))
      if(identical(data[[nx]][both_valid], data[[ny]][both_valid])) {
        new_var <- data[[nx]]
        new_var[is.na(new_var)] <- data[[ny]][is.na(new_var)] # copy from data[[ny]]
        data[[nn]] <- new_var
        data[[nx]] <- NULL
        data[[ny]] <- NULL
      }
    }
  }
  data
} 

if(FALSE) { # tests
  da <- data.frame(
    id = 1:5, 
    g = c(1,2,NA,4,5), 
    h = factor(c('a','b','c','d','e')),
    hh = c('a','b','c','d','e'),
    m = factor(c('M',NA,NA,'F','M'))  # factor with same set of values
  )
  db <- data.frame(
    id = c(1:3,6), 
    g = c(1,2,NA,6), 
    h = factor(c('a','b','c','f')),  # factor with different sets of values
    hh = c('a','b','c','f'),
    k = c(1:3,6),
    m = factor(c('M','M',NA,'F'))
  )
  da  
  db  
  dm <- dplyr::full_join(da, db, by = 'id' )  
  dm
  fix_xy(dm)
}


