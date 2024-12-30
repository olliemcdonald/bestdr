#' new_transition
#' constructor for object of class \code{transition}
#'
#' @param name name for the transition event
#' @param parent the parent type for the transition
#' @param offspring the update to the system after the transition
#' @param model the user-defined bp_model output of process_model
#'
#' @return A new branching process rate object
#' @export
new_transition <- function(name, parent, offspring, model) {
  trans <- list("name" = name, "parent" = parent, "offspring" = offspring, "model" = model)
  class(trans) <- "transition"
  return(trans)
}

#' validate_transition
#' verifies the correctness of an object of class \code{transition}
#' 
#' @param trans_obj a transition object
#'
#' @return The \code{transition} object if it is valid, throws an error otherwise
validate_transition <- function(trans_obj) {
  if (class(trans_obj) != "transition") {
    stop("not a valid estipop transition object!")
  }
  if (!is.numeric(trans_obj$parent) || !is.numeric(trans_obj$offspring)) {
    stop("parent and offpspring must be numeric!")
  }
  if (length(trans_obj$parent) > 1) {
    stop("parent should not be a vector!")
  }
  if (trans_obj$parent <= 0 || any(trans_obj$offspring < 0)) {
    stop("parent must be positive and offspring must be nonnegative!")
  }
  if (trans_obj$parent > length(trans_obj$offspring)) {
    stop("parent index is larger than length of offspring vector!")
  }
  if (class(trans_obj$model) != "formula") {
    stop("transition statistical model is not a formula")
  }
  return(trans_obj)
}

#' transition
#' constructs and validates an object of class \code{transition}
#'
#' @param name a name for the transition event
#' @param parent the parent type for the transition
#' @param offspring the update to the system after the transition
#' @param model the user-defined bp_model output of process_model
#'
#' @return The \code{transition} object if it is valid, throws an error otherwise
#' @export
transition <- function(name = NULL, parent, offspring, model=NULL) {
  return(validate_transition(new_transition(name, parent, offspring, as.formula(model))))
}

#' constructor for object of class \code{bp_model}
#'
#' @param ... a transition list structured as a vector of transitions
#' 
#' @return A new branching process model object
#' 
bp_model <- function(...) {
  pm <- list(transition_list = list(...))
  class(pm) <- "bp_model"
  return(pm)
}

#' validate_process_model
#' verifies the correctness of an object of class \code{bp_model}
#' 
#' @param bp_model a bp_model object containing user-defined transitions
#' @param parameter_constraints vector of bounds on the parameters
#' @param priors vector of characters listing priors for any parameters
#' @param predictor_names character vector for the name of the predictors used in the formulates (ex. x, c, concentration)
#' @param observation_error logical to include error due to observation
#' @param hierarchical vector of parameters that are random effects and have
#'     hyperparameters
#'
#' @return The \code{process_model} object if it is valid, throws an error otherwise
#' 
#' 
validate_process_model <- function(bp_model,
                                   parameter_constraints = NULL,
                                   priors = NULL,
                                   predictor_names = NULL,
                                   observation_error = FALSE,
                                   hierarchical = NULL) {
  if (!is.list(bp_model$transition_list) || length(bp_model$transition_list) == 0) {
    stop("transition_list must be a list with a positive number of elements!")
  }
  if (class(bp_model) != "bp_model") {
    stop("invalid process model!")
  }
  lapply(bp_model$transition_list, function(b) {
    if (class(b) != "transition") {
      stop("invalid transition object!")
    }
  })

  bp_model$ntypes <- length(bp_model$transition_list[[1]]$offspring)
  bp_model$nevents <- length(bp_model$transition_list)
  bp_model$parameter_constraints <- parameter_constraints
  bp_model$priors <- priors
  bp_model$predictor_names <- predictor_names
  bp_model$observation_error <- observation_error
  bp_model$hierarchical <- hierarchical

  # Test that there are a consistent number of offspring in each transition
  lapply(bp_model$transition_list, function(b) {
    if (length(b$offspring) != bp_model$ntypes) {
      stop("transitions have inconsistent number of types!")
    }
  })

  # Check that the parameter names in constraints all exist in the formulas
  params_from_constraints <- sub("\\[.*\\]", "", bp_model$parameter_constraints)
  params_from_formulas <- character()
  for (tr in bp_model$transition_list) {
    params_from_formulas <- c(params_from_formulas, all.vars(tr$model), tr$name)
  }
  if (!all(params_from_constraints %in% params_from_formulas)) {
    stop("some parameters defined in the constraint list aren't part of the model formulas")
  }

  return(bp_model)
}


#' process_model
#' constructs and validates an object of class \code{process_model}
#'
#' @param bp_model a bp_model object containing user-defined transitions
#' @param parameter_constraints vector of bounds on the parameters
#' @param priors vector of characters listing priors for any parameters
#' @param predictor_names character vector for the name of the predictors used in the formulates (ex. x, c, concentration)
#' @param observation_error logical to include error due to observation
#' @param hierarchical vector of parameters that are random effects and have
#'     hyperparameters
#'
#' @return The \code{process_model} object if it is valid, throws an error otherwise
#' @export
process_model <- function(bp_model,
                          parameter_constraints = NULL,
                          priors = NULL,
                          predictor_names = NULL,
                          observation_error = FALSE,
                          hierarchical = NULL) {
  proc_model <- validate_process_model(
    bp_model,
    parameter_constraints,
    priors,
    predictor_names,
    observation_error,
    hierarchical
  )
  proc_model <- parse_model(proc_model)
  return(proc_model)
}
