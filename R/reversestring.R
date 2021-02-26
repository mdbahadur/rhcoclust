# This is the function for reverse string
#' @export
reversestring <- function (string, n = 1)
  {
    lens <- nchar(string)
    sapply(1:length(string), function(i) {
      #tmp <- substring(string[i], 1, )
      paste0(c(substring(string[i], seq(lens[i] - n + 1, 1,
                                      -n), seq(lens[i], n, -n)), substr(string[i], 1,
                                                                           lens[i]%%n)), collapse = "")
    }, USE.NAMES = F)
  }
