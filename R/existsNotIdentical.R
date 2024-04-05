#' @title existsNotIdentical
#'
#' @description Determines if the input variable name exists and is identical to a given variable. This is used to determine if a saved variable needs to be re-created or not with the additional note that if the variable does not exist, it must be (re-)created anyway.
#'
#' @param variable String of variable name to investigate
#' @param compare Actual variable that the string's variable should be compared with
#'
#' @return If false if the variable exists and is identical, and true if either the variable does not exist, or is different from the comparison
#' @examples
#' x <- "hello"
#' existsNotIdentical("x", "hello") # FALSE
#' existsNotIdentical("x", "hello.world") # TRUE
#' existsNotIdentical("y", "hello") # TRUE
#' @export
#'
#'

existsNotIdentical <- function(variable, compare) {
  if (exists(variable)) { # The variable exists...
    if (identical(get(variable),compare)) { # ...and it is identical
      return(FALSE)
    } else { # ...and it is not identical
      return(TRUE)
    }
  } else { # The variable straight up does not exist
    return(TRUE)
  }
}
