#' Format the constraints from a processed model
#'
#' @param constraint a constraint on a parameter or data variable
#' 
#'
format_constraint <- function(constraint) {
  constraint <- gsub(" ", "", constraint, fixed = TRUE)
  bounds <- stringi::stri_split(constraint, fixed = ",", simplify = T)

  lowerb <- bounds[1]
  upperb <- bounds[2]
  lowerflag <- nchar(lowerb) > 0
  upperflag <- nchar(upperb) > 0

  if (lowerflag && upperflag) {
    sprintf("<lower=%s,upper=%s>", lowerb, upperb)
  } else if (lowerflag && !upperflag) {
    sprintf("<lower=%s>", lowerb)
  } else if (!lowerflag && upperflag) {
    sprintf("<upper=%s>", upperflag)
  } else {
    ""
  }
}

#' Separate the name of the variable from the constraint
#'
#' @param name_with_constraint a name and constraint text chunk
#'
#'
separate_name_and_constraint <- function(name_with_constraint) {
  name <- gsub("\\[.*\\]", "", name_with_constraint)
  constraint <- stringi::stri_extract(name_with_constraint, regex = "(?<=\\[).*?(?=\\])")
  if (is.na(constraint)) {
    return(list(name = name, constraint = ""))
  }
  constraint <- format_constraint(constraint)
  return(list(name = name, constraint = constraint))
}
 
#' Function to parse parameters and separate name and constraint.
#' Requirement here to include all parameters and any possible constraints
#'
#' @param parameters parameters block from a model
#'
#'
parse_parameters <- function(parameters) {
  param_df <- data.frame(name = character(), param_block = character(), constraint = character())
  for (param in parameters) {
    param <- separate_name_and_constraint(param)
    param$param_block <- paste0("real", param$constraint, " ", param$name, ";")
    param_df <- rbind(param_df, param)
  }
  return(param_df)
}

#' Function to parse formula and extract parts
#'
#' @param formulas formulas from a model
#' @param hierarchical hierarchical parameters from a model
#'
#'
parse_priors <- function(formulas, hierarchical) {
  parsed_priors <- data.frame()

  for (formula in formulas) {
    ## PARAMETER
    # Extracting parameters
    prior_as_formula <- as.formula(formula)
    param_with_constraint <- deparse(prior_as_formula[[2]])
    param <- separate_name_and_constraint(param_with_constraint)

    # PRIORS
    # Extract prior distribution
    full_prior <- deparse(prior_as_formula[[3]])
    prior_dist <- gsub("\\(.*", "", full_prior)

    # Extract prior hyperparameters
    hyperparams_with_constraints <- stringi::stri_extract(full_prior, regex = "(?<=\\().*?(?=\\))")
    hyperparams_with_constraints <- stringi::stri_match_all(hyperparams_with_constraints, regex = "[A-Za-z0-9_\\.]+(?>\\[[^\\[\\])]*,[^\\[\\])]*\\])*")
    hyperparam_df <- data.frame(name = character(), constraint = character())
    for (hp in hyperparams_with_constraints[[1]]) {
      hyperparam_df <- rbind(hyperparam_df, separate_name_and_constraint(hp))
    }

    # CONVERT PARAMETER AND PRIOR SPECIFICATION TO STAN CODE
    # model block <- paramer ~ prior_dist(hyperparam_list)
    # data block <- Hyperparameters
    # parameter block <- parameters
    params <- param$name
    hyperparams <- hyperparam_df$name
    model_block <- paste0(param$name, " ~ ", prior_dist, "(", paste(hyperparam_df$name, collapse = ","), ");")
    data_block <- paste(paste0("real", hyperparam_df$constraint, " ", hyperparam_df$name, ";"), collapse = "\n")
    prior_blocks <- data.frame(
      paramNames = params,
      hyperparamNames = I(list(hyperparams)),
      model_block = model_block,
      data_block = data_block,
      hierarchical = params %in% hierarchical
    )
    parsed_priors <- rbind(parsed_priors, prior_blocks)
  }

  return(parsed_priors)
}

#' Function to parse a single formula formula and extract parts
#'
#' @param formula formula block
#' @param predictor_names names of the predictor values (concentration, etc)
#'
#'
parse_single_formula <- function(formula, predictor_names){
  varNames <- all.vars(formula) # parameter and variable names from formula
  varNames_wConstraints <- as.vector(stringi::stri_match_all(deparse(formula[[3]]), regex = "[A-Za-z0-9_\\.]+(?>\\[[^\\[\\])]*,[^\\[\\])]*\\])", omit_no_match = TRUE)[[1]])

  varNamesData <- predictor_names
  ## for prediction we will need to know those which are in RHS
  form2 <- formula
  form2[[2L]] <- 0
  varNamesRHS <- all.vars(form2)
  rateNamesLHS <- formula[[2L]]
  varNamesRHS <- varNamesRHS[!(varNamesRHS %in% varNamesData)]

  # Substitute in the constrained variables into varNamesRHS
  varNames_noConstraint <- varNamesRHS[!varNamesRHS %in% sub("\\[.*\\]", "", varNames_wConstraints)]
  varNamesRHS <- c(varNames_wConstraints, varNames_noConstraint)

  # Remove constraints from RHS
  updateRHS <- stringi::stri_replace_all(deparse(form2[[3L]]), regex = "(?>\\[[^\\[\\])]*,[^\\[\\])]*\\])+", "")

  # Substitute X for the data variable
  # XXX NEED TO FIX THIS LATER
  update_predictor <- str2lang("xval[xi]")
  updateRHS <- do.call("substitute", list(str2lang(updateRHS), `names<-`(list(update_predictor), predictor_names)))
  return(list(updateRHS=updateRHS,
              rateNamesLHS=deparse(rateNamesLHS),
              varNamesRHS=varNamesRHS))
}

#' Function to parse formula when there is no model
#'
#' @param formula formula block
#' @param rateName string name to identify the rate
#'
#'
parse_single_formula_nomodel <- function(formula, rateName){
  return(list(updateRHS=str2lang(rateName),
              rateNamesLHS=rateName,
              varNamesRHS=rateName))
}

#' Function to parse all formulas from a model
#'
#' @param model formula block
#'
#'
parse_formulas <- function(model) {
  formula_df <- data.frame()
  for (tr in model$transition_list) {
    ## PARAMETER
    # Extracting parameters
    formula <- tr$model
    formula <- as.formula(formula)

    if(length(formula) == 0){
      updateForm <- parse_single_formula_nomodel(formula, tr$name)
    } else{
      updateForm <- parse_single_formula(formula, model$predictor_names)
    }

    # Add hierarchical pieces
    if (!is.null(model$hierarchical)) {
      hierarchicalVars <- sapply(as.list(paste0(model$hierarchical, "[g]")), str2lang)
      names(hierarchicalVars) <- sapply(model$hierarchical, str2lang)
      #if(length(formula) != 0){
        updateForm$updateRHS <- do.call("substitute", list(updateForm$updateRHS, hierarchicalVars))
      #}
      updateForm$hierarchical <- any(model$hierarchical %in% updateForm$varNamesRHS)
    }else{
      updateForm$hierarchical <- FALSE
    }
    updateForm$updateRHS <- deparse(updateForm$updateRHS)

    split_formula <- data.frame(
      rateNames = updateForm$rateNamesLHS,
      RHS = updateForm$updateRHS,
      paramNames = I(list(updateForm$varNamesRHS)),
      hierarchical = updateForm$hierarchical
    )

    formula_df <- rbind(formula_df, split_formula)
  }
  return(formula_df)
}

#' Function to merge all constraints from the blocks of the model definition
#'
#' @param params_from_constraints all parameters pulled from the constraints arguments
#' @param params_from_model all parameters pulled from the model arguments
#'
#'
merge_constraints <- function(params_from_constraints, params_from_model) {
  # Merge parameters from model and constraints
  # Parameters from Constraints will supercede those from the model, but need to merge
  all_params <- union(params_from_constraints$name, params_from_model$name)
  #params_from_model[!params_from_model$name %in% params_from_constraints$name, ]
  all_params <- merge(params_from_model, params_from_constraints,
    by = "name",
    nodups = FALSE,
    all = T,
    suffixes = c(".model", "")
  )
  # Go through parameters and compare constraints and issue warnings when there is an issue
  for (i in seq_len(nrow(all_params))) {
    paramName <- all_params$name[i]
    constraint.model <- all_params$constraint.model[i]
    constraint <- all_params$constraint[i]

    if (is.na(constraint.model)) stop(paste0("parameter ", paramName, " not specified in model"))

    if (is.na(constraint) | (nchar(constraint) == 0 && nchar(constraint.model) > 0)) {
      # If no specified constraint then use the constraint from the model block (if it exists)
      all_params$constraint[i] <- all_params$constraint.model[i]
      all_params$param_block[i] <- all_params$param_block.model[i]
    }
    # Otherwise do nothing and just take those last columns
  }

  return(all_params[, c("name", "constraint", "param_block")])
}

#' Function to build the parants and list of events into a vector and matrix
#'
#' @param model bpmodel structure
#'
#'
build_parents_and_eventmatrix <- function(model) {
  model_deets <- data.frame()
  for (tr in model$transition_list) {
    model_deets <- rbind(model_deets, data.frame(name = tr$name, parent = tr$parent, offspring = I(list(tr$offspring))))
  }
  list(
    name = model_deets$name,
    parent_vec = model_deets$parent,
    e_mat = do.call("rbind", model_deets$offspring)
  )
}

#' Function to parse through the model and get the required pieces to construct
#' the stan code
#'
#' @param model model structure
#'
#'
parse_model <- function(model) {
  # Priors
  priors_blocks <- parse_priors(model$priors, model$hierarchical)

  # Parameters from Constraints
  params_from_constraints <- parse_parameters(model$parameter_constraints)
  # Rate and formula from model
  rate_blocks <- parse_formulas(model)
  params_from_model <- parse_parameters(unlist(rate_blocks$paramNames))
  param_blocks <- merge_constraints(params_from_constraints, params_from_model)

  # Adjust parameters for hierarchical model
  if(!is.null(model$hierarchical)){
    hierarcIdx <- param_blocks$name %in% model$hierarchical
    param_blocks$param_block[hierarcIdx] <-
      paste0("vector", param_blocks$constraint[hierarcIdx], "[ng] ", param_blocks$name[hierarcIdx], ";")
    #paste0("vector<lower=0>", param_blocks$constraint[hierarcIdx], "[ng] ", param_blocks$name[hierarcIdx], ";")
  }


  # Extract model info
  rate_info <- build_parents_and_eventmatrix(model)
  model$event_matrix <- rate_info$e_mat
  model$parent_vec <- rate_info$parent_vec

  rateIndex <- match(rate_blocks$rateNames, rate_info$name)
  # order rate_block section to follow the defined matrix (should already be in order)
  rate_blocks <- rate_blocks[rateIndex, ]

  bd_eventFlag <- all(model$parent_vec == 1) && (identical(model$event_matrix, matrix(c(2,0), ncol=1)) || identical(model$event_matrix, matrix(c(0,2), ncol=1)))
  model$birthdeath <- ifelse(bd_eventFlag, TRUE, FALSE)

  if(model$birthdeath){
    rate_blocks$stan_index <- rep(NA, length = nrow(rate_blocks))
    for(i in 1:nrow(rate_blocks)){
      if(rate_blocks$hierarchical[i] && (rate_blocks$rateNames[i] == rate_blocks$RHS[i])) rate_blocks$stan_index[i] <- "[g]"
      else if(!rate_blocks$hierarchical[i] && (rate_blocks$rateNames[i] == rate_blocks$RHS[i])) rate_blocks$stan_index[i] <- ""
      else if(rate_blocks$hierarchical[i] && (rate_blocks$rateNames[i] != rate_blocks$RHS[i])) rate_blocks$stan_index[i] <- "[g, xi]"
      else if(!rate_blocks$hierarchical[i] && (rate_blocks$rateNames[i] != rate_blocks$RHS[i])) rate_blocks$stan_index[i] <- "[xi]"
    }


    temp <- paste0(rate_blocks$rateNames, rate_blocks$stan_index)
    rate_blocks$stan_block <- paste0(temp, " = ", rate_blocks$RHS, ";")
    rate_blocks$stan_block[rate_blocks$rateNames == rate_blocks$RHS] <- ""
    rate_blocks$bd <- ifelse(model$event_matrix == 2, "b", "d")
  }else{
    rate_blocks$stan_block <- paste0("r_mat[", 1:nrow(rate_blocks), ", p_vec[", 1:nrow(rate_blocks), "]] = ", rate_blocks$RHS, ";")
  }


  model$blocks <- list(
    rate_model = rate_blocks,
    params = param_blocks,
    priors = priors_blocks
  )

  return(model)
}
