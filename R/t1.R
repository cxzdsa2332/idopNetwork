#' @title calculate log-sum-exp values
#' @param v numeric vector
#' @return log-sum-exp values
#' @examples logsumexp(c(100, 1000, 10000))
#' @export
logsumexp <- function(v) {
  vm = max(v)
  log(sum(exp(v-vm))) + vm
}