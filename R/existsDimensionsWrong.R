existsDimensionsWrong <- function(variable, compare) {
  if (exists(variable)) { # The variable exists...
    if (all(dim(get(variable))==length(compare))) { # ...and the dimensions are all correct
      return(FALSE)
    } else { # ...and one of the dimensions is incorrect
      return(TRUE)
    }
  } else { # The variable straight up does not exist
    return(TRUE)
  }
}
