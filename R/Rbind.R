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
#' @export
Rbind <- function(..., verbose = FALSE) {
  nam <- '.unique.'
  x <- list(...)
  x <- lapply(x, as.data.frame)
  nams <- unique(unlist(sapply(x, names)))
  print(nams)
  while( nam %in% nams) nam <- paste0(nam,".")
  if(verbose) print(nam)
  x <- lapply(seq_along(x), function(ii) {
    x[[ii]][[nam]] <- ii
    x[[ii]]
  })
  fun <- function(x,y) merge(x, y, all = TRUE, sort = FALSE)
  ret <- Reduce(fun,x)
  ret[[nam]] <- NULL
  ret
}
