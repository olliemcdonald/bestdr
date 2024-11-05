# Function where user defines for the given data frame what the parent and offspring vector
# types are for current and previous times
# as well as the predictor variable (i.e. concentration/logconc)
# renames the concentration vector to x (or whatever the predictor was defined in the model above)

#' attach_data
#' Function to set up the data for import into stan based on the model and data
#' frame defined by the user.
#'
#' @param model bp_model defined by the user
#' @param dat data frame that contains counts for current and previous time and
#'      predictor variables
#' @param count_vars character vector of column names in dat that are the
#'      current time object counts
#' @param prevcount_vars character vector of column names in dat that are
#'      the previous time object counts
#' @param predictor_vars character of variable name in dat that is the predictor
#'      variable for the statistical model
#' @param hierarchical_group character of variable name representing the group
#'      variable if the model has random effects
#'
#' @export
attach_data <- function(model, dat,
                        time_var = "dt",
                        count_vars = NULL,
                        prevcount_vars = NULL,
                        predictor_vars = NULL,
                        hierarchical_group = NULL) {
  
  dat <- as.data.frame(dat)
  
  if(is.null(predictor_vars)){
    data.list <- list(
      N = nrow(dat),
      x_idx = rep(1, nrow(dat)),
      xval = 1,
      nx = 1
    )
  }else{
    x <- get(predictor_vars, dat)
    xvals <- unique(x)
    data.list <- list(
      N = nrow(dat),
      x_idx = match(x, xvals),
      xval = xvals,
      nx = length(xvals)
    )
  }
  if(model$birthdeath){
    if(length(count_vars) > 1) stop("more count variable types than 1")
    addtl.data.list <-  list(count = dat[, count_vars],
                          count_prev = dat[, prevcount_vars],
                          dt = dat[, time_var])
  } else {
    addtl.data.list <- list(ntypes = model$ntypes,
                         nevents = model$nevents,
                         n_dt_unique = length(unique(dat$dt)),
                         count = as.matrix(dat[, count_vars]),
                         count_prev = as.matrix(dat[, prevcount_vars]),
                         e_mat = model$event_matrix,
                         p_vec = model$parent_vec,
                         times_idx = match(dat$dt, unique(dat$dt)),
                         dt = array(unique(dat[, time_var])))
  }
  data.list <- c(data.list,
                 addtl.data.list)
  if (!is.null(model$hierarchical)) {
    # Check that hierarchical group is defined from the data.frame
    groups <- get(hierarchical_group, dat)
    g_idx <- as.numeric(as.factor(groups))
    ng <- length(unique(groups))
    data.list <- c(data.list, list(g_idx = g_idx, ng = ng))
  }
  return(data.list)
}
