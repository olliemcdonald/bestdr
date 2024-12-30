#' Generate Stan code from a provided bp model
#'
#' @param model A bp model that has been processed and parsed
#'
#' @export
#'
create_stan_code <- function(model) {

  if (!model$birthdeath && !is.null(model$hierarchical)) {
    fileName <- system.file("stan_custommodel/ktype_template_me.stan", package="bestdr")
  } else if(!model$birthdeath && is.null(model$hierarchical)) {
    fileName <- system.file("stan_custommodel/ktype_template_nome.stan", package="bestdr")
  } else if(model$birthdeath && !is.null(model$hierarchical)) {
    fileName <- system.file("stan_custommodel/onetype_template_me.stan", package="bestdr")
  } else if(model$birthdeath && is.null(model$hierarchical)) {
    fileName <- system.file("stan_custommodel/onetype_template_nome.stan", package="bestdr")
  } else{
    stop("Model does not fit the templates provided")
  }
  stan_template <- readChar(fileName, file.info(fileName)$size)


  # Fill in blocks from model data
  datablock <- create_stan_data_block(model)
  paramblock <- create_stan_parameter_block(model)

  modelblock <- create_stan_statmodel_block(model)
  priorblock <- create_stan_priors_block(model)
  error_additions <- create_stan_error_block(model)

  # Concatenate the error piece onto the other blocks
  datablock <- paste0(datablock, error_additions$data_block)
  paramblock <- paste0(paramblock, error_additions$parameter_block)
  priorblock <- paste0(priorblock, error_additions$prior_block)
  countmodelblock <- error_additions$count_model_block


  if(model$birthdeath){
    tparamblock <- create_transparameter_block(model)
    bpmodelblock <- create_stan_bdmodel_block(model)
    return(sprintf(stan_template, datablock, paramblock, tparamblock, modelblock, priorblock, bpmodelblock, countmodelblock))
  }else{
    return(sprintf(stan_template, datablock, paramblock, modelblock, priorblock, countmodelblock))
  }

}

#' Create Stan parameter block from parsed text
#'
#' @param model A bp model that has been processed and parsed
#'
#'
create_stan_parameter_block <- function(model) {
  paramtext <- paste(model$blocks$params$param_block, collapse = "\n")
  hierparamtext <- paste(model$blocks$priors$data_block[model$blocks$priors$hierarchical == TRUE], collapse = "\n")
  return(paste(paramtext, hierparamtext, sep = "\n"))
}

#' Create Stan transformed parameter block from parsed text
#'
#' @param model A bp model that has been processed and parsed
#'
#'
create_transparameter_block <- function(model) {
  birthVar <- model$blocks$rate_model$rateNames[model$blocks$rate_model$bd == "b"]
  deathVar <- model$blocks$rate_model$rateNames[model$blocks$rate_model$bd == "d"]

  param_type <- ifelse(model$blocks$rate_model$hierarchical, "matrix<lower=0>[ng,nx]", "vector<lower=0>[nx]")
  b_idx <- param_type[model$blocks$rate_model$bd == "b"]
  d_idx <- param_type[model$blocks$rate_model$bd == "d"]

  paramDecl <- c(paste0(b_idx, " ", birthVar, ";"), paste0(d_idx, " ", deathVar, ";"))
  paramDecl <- ifelse(model$blocks$rate_model$rateNames == model$blocks$rate_model$RHS, "", paramDecl)
  return(paste(paramDecl, collapse = "\n"))
}

#' Create Stan data block from parsed text
#'
#' @param model A bp model that has been processed and parsed
#'
#'
create_stan_data_block <- function(model) {
  return(paste(model$blocks$priors$data_block[model$blocks$priors$hierarchical == FALSE], collapse = "\n"))
}

#' Create Stan priors block from parsed text
#'
#' @param model A bp model that has been processed and parsed
#'
#'
create_stan_priors_block <- function(model) {
  return(paste(model$blocks$priors$model_block, collapse = "\n"))
}

#' Create Stan statistical model block from parsed text
#'
#' @param model A bp model that has been processed and parsed
#'
#'
create_stan_statmodel_block <- function(model) {
  return(paste(model$blocks$rate_model$stan_block, collapse = "\n"))
}

#' Create Stan observed error chunk for model from parsed text
#'
#' @param model A bp model that has been processed and parsed
#'
#'
create_stan_error_block <- function(model) {
  if (model$birthdeath) {
    if (model$observation_error) {
      return(list(
        data_block = "\nreal<lower=0> prior_s_obs_err;\n",
        parameter_block = "\nreal<lower=0> obs_err;\n",
        prior_block = "\nobs_err ~ normal(0, prior_s_obs_err);\n",
        count_model_block = "\nz ~ normal(mu, sigma + obs_err * mu);\n"
      ))
    } else {
      return(list(
        data_block = "",
        parameter_block = "",
        prior_block = "",
        count_model_block = "\nz ~ normal(mu, sigma);\n"
      ))
    }
  } else {
    if (model$observation_error) {
      return(list(
        data_block = "\nreal<lower=0> prior_s_obs_err;\n",
        parameter_block = "\nreal<lower=0> obs_err;\n",
        prior_block = "\nobs_err ~ normal(0, prior_s_obs_err);\n",
        count_model_block = "\nz[n] ~ multi_normal(mu_t, sigma_t + obs_err*diag_matrix(mu_t));\n"
      ))
    } else {
      return(list(
        data_block = "",
        parameter_block = "",
        prior_block = "",
        count_model_block = "\nz[n] ~ multi_normal(mu_t, sigma_t);\n"
      ))
    }
  }
}

#' Create Stan mechanstic model block from parsed text
#'
#' @param model A bp model that has been processed and parsed
#'
#'
create_stan_bdmodel_block <- function(model){

  birthVar <- model$blocks$rate_model$rateNames[model$blocks$rate_model$bd == "b"]
  deathVar <- model$blocks$rate_model$rateNames[model$blocks$rate_model$bd == "d"]

  param_idx <- rep(NA, nrow(model$blocks$rate_model))
  for(i in seq_along(param_idx)){
    if(model$blocks$rate_model$rateNames[i] == model$blocks$rate_model$RHS[i]){
      if(model$blocks$rate_model$hierarchical[i]) param_idx[i] <- "[g_idx[i]]"
      else param_idx[i] <- ""
    } else{
      if(model$blocks$rate_model$hierarchical[i]) param_idx[i] <- "[g_idx[i], x_idx[i]]"
      else param_idx[i] <- "[x_idx[i]]"
    }
  }

  birthVar <- paste0(birthVar, param_idx[model$blocks$rate_model$bd == "b"])
  deathVar <- paste0(deathVar, param_idx[model$blocks$rate_model$bd == "d"])

  model_block <- c(
    paste0("sigma[i] = sqrt(count_prev[i] * (", birthVar, "+", deathVar, ") / (", birthVar, "-", deathVar, ") * (exp(2*(", birthVar, "-", deathVar, ")*dt[i]) - exp((", birthVar, "-", deathVar, ")*dt[i])));"),
    paste0("mu[i] = count_prev[i] * exp((", birthVar, "-", deathVar, ") * dt[i]);")
  )

  return(paste(model_block, collapse = "\n"))
}
