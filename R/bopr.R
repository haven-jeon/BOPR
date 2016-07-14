pkg.env <- new.env()

pkg.env$MAX_SURPRISE <- 5
pkg.env$MIN_SURPRISE <- pkg.env$MAX_SURPRISE * -1
pkg.env$PRIOR_VARIANCE <- 1.0

#' BOPR
#'
#' Bayesian online learning scheme for probit regression (BOPR)
#'
#' @param x a  matrix of predictors.
#' @param y a factor vector with 2 level
#' @param beta scaling parameter
#' @param prior_prob prior of initial parametes
#' @param epsilon parameter to apply dynamics
#' @param formula an optional data frame in which to interpret the variables named in the formula.
#' @param data an optional data frame in which to interpret the variables named in the formula.
#' @param subset optional expression saying that only a subset of the rows of the data should be used in the fit.(currently it's not working.)
#' @param na.action a function which indicates what should happen when the data contain \code{NA}s.
#' @param ... not used
#' @return S3 \code{BOPR} object; a list of consisting of
#'  \item{beta_matrix}{beta matrix with mean and variance}
#'  \item{beta}{scaling parameter}
#'  \item{prior_prob}{prior of initial parametes}
#'  \item{epsilon}{parameter to apply dynamics}
#'  \item{formula}{formula}
#' @author Heewon Jeon \email{madjakarta@@gmail.com}
#' @references Graepel, Thore, et al. "Web-scale bayesian click-through rate prediction for sponsored search advertising in microsoft's bing search engine." Proceedings of the 27th International Conference on Machine Learning (ICML-10). 2010.
#' He, X., Bowers, S., Candela, J. Q., Pan, J., Jin, O., Xu, T.,Herbrich, R. (2014). Practical Lessons from Predicting Clicks on Ads at Facebook. Proceedings of 20th ACM SIGKDD Conference on Knowledge Discovery and Data Mining - ADKDD '14, 1-9.
#' @examples
#' idx  <- sample(1:nrow(credit_approval))
#' first_train_set  <- credit_approval[idx[1:200],]
#' second_train_set  <- credit_approval[idx[201:400],]
#' test_set <- credit_approval[idx[401:690],]
#'
#' bopr_mdl <- BOPR(A16 ~ A1 + A4 + A5 + A7 + A9 + A10 + A12 + A13 , first_train_set)
#' bopr_mdl_up  <- online_leraning(bopr_mdl, second_train_set)
#' pred  <- predict(bopr_mdl_up, test_set)
#' @export
#' @import stats
BOPR <- function(x, ...) UseMethod("BOPR")

#' @rdname BOPR
#' @export
BOPR.default <- function(x, y,
                                beta=0.05,
                                prior_prob=0.5,
                                epsilon=0.05,
                                subset=NULL,
                                ...)
{
  #init bias prior
  num_feat <- ncol(x)
  weight_mat <- matrix(nrow = 2,ncol = num_feat,data = 0, dimnames=list(c("mean", "variance"), colnames(x)))
  bias_mean <- prior_mean(prior_prob, beta, num_feat)
  #added bias mean and variance
  weight_mat[,seq(1,num_feat)] <- c(bias_mean, pkg.env$PRIOR_VARIANCE)
  if(!is.null(subset) & is.numeric(subset)){
    x <- x[subset,]
    y <- y[subset]
  }
  #train
  y_lab <- ifelse(y == 1, 1, -1)

  #active_mean_variance
  for(i in 1:nrow(x)){
    total_mean <- sum(weight_mat[1,])
    total_variance <- sum(weight_mat[2,]) + beta ** 2
    t <- y_lab[i] * total_mean / sqrt(total_variance)
    t <- clip(t, pkg.env$MIN_SURPRISE, pkg.env$MAX_SURPRISE )
    v <- dnorm(t) / pnorm(t)
    w <- v * (v + t)
    col_nm <- colnames(x)
    if(length(col_nm) != ncol(weight_mat))
      stop("num of features is not match.")
    for(j in col_nm){
        if(x[i,j] == 0) next
        mean_delta <- y_lab[i] * weight_mat['variance',j] / sqrt(total_variance) * v
        variance_multiplier <- 1.0 - weight_mat['variance',j] / total_variance * w
        updated <- c(
          mean = as.numeric(weight_mat['mean',j] + mean_delta),
          variance = as.numeric(weight_mat['variance',j] * variance_multiplier))
        #apply dynamics
        prior <- c(mean=0, variance=pkg.env$PRIOR_VARIANCE)
        adjusted_variance <- updated['variance'] * prior['variance'] /
        ((1.0 - epsilon) * prior['variance'] +
           epsilon * updated['variance'])
        adjusted_mean <- adjusted_variance * (
          (1.0 - epsilon) * updated['mean'] / updated['variance'] +
            epsilon * prior['mean'] / prior['variance'])

        weight_mat[,j] <- c(mean=adjusted_mean, variance=adjusted_variance)
    }

  }
  mdl_mat <- list(beta_matrix=weight_mat, beta=beta, prior_prob=prior_prob,epsilon=epsilon)
  class(mdl_mat) <- 'BOPR'
  return(mdl_mat)
}



#' @rdname BOPR
#' @export
BOPR.formula <- function(formula, data, subset=NULL, na.action =na.pass, beta=0.05,
                                prior_prob=0.5,
                                epsilon=0.05, ...)
{
  Call <- match.call()

  indx <- match(c("formula", "data", "subset"),
                names(Call), nomatch = 0L)
  if (indx[1] == 0L) stop("a 'formula' argument is required")
  temp <- Call[c(1L, indx)]      # only keep the arguments we wanted
  temp$na.action <- na.action    # This one has a default
  temp[[1L]] <- quote(stats::model.frame) # change the function called
  m <- eval.parent(temp)

  Terms <- attr(m, "terms")

  y <- model.extract(m, "response")

  m <- m[,-1,drop = FALSE]

  if(any(sapply(m, is.numeric) == TRUE)){
    stop("All variables need to be factor!")
  }

  character_vars <- lapply(m, class) == "character"
  m[, character_vars] <- lapply(m[, character_vars], as.factor)

  mdl_mat <- model.matrix(as.formula(paste("~",paste0(colnames(m), collapse=' + '))),
              m, contrasts.arg = lapply(m,contrasts, contrasts=FALSE))[,-1]

  object <- BOPR.default(x=mdl_mat, y=y, beta=beta, prior_prob=prior_prob,
                      epsilon=epsilon,subset=subset, ...)
  object$formula <- formula
  return(object)
}


.make_training_data <- function(formula, data, subset=NULL, na.action =na.pass)
{
  Call <- match.call()

  indx <- match(c("formula", "data", "subset"),
                names(Call), nomatch = 0L)
  if (indx[1] == 0L) stop("a 'formula' argument is required")
  temp <- Call[c(1L, indx)]      # only keep the arguments we wanted
  temp$na.action <- na.action    # This one has a default
  temp[[1L]] <- quote(stats::model.frame) # change the function called
  m <- eval.parent(temp)

  Terms <- attr(m, "terms")

  y <- model.extract(m, "response")

  m <- m[,-1,drop = FALSE]

  if(any(sapply(m, is.numeric) == TRUE)){
    stop("All variables need to be factor!")
  }

  character_vars <- lapply(m, class) == "character"
  m[, character_vars] <- lapply(m[, character_vars], as.factor)

  mdl_mat <- model.matrix(as.formula(paste("~",paste0(colnames(m), collapse=' + '))),
                          m, contrasts.arg = lapply(m,contrasts, contrasts=FALSE))[,-1]
  return(list(y=y, x=mdl_mat))
}



#' Predict new samples using a BOPR model
#'
#' This function produces predicted values from a BOPR model.
#'
#' @param object \code{BOPR} object
#' @param newdata a data.frame of predictors
#' @param type either "class" for the predicted class or "prob" for model confidence values.
#' @param na.action when using a formula for the original model fit, how should missing values be handled?
#' @param ... other options (not currently used)
#'
#' @return when type = "class", a factor vector is returned. When type = "response", probability is returned.
#' @export
predict.BOPR <- function(object, newdata = NULL, type = "response", na.action = na.pass, ...){
  if(class(object) != 'BOPR'){
    stop('object is not BOPR class!')
  }
  if(any(class(newdata) == 'data.frame')) {
    variables <- labels(terms(object$formula))
    newdata_ <- newdata[,variables]
    character_vars <- lapply(newdata_, class) == "character"
    newdata_[, character_vars] <- lapply(newdata_[, character_vars], as.factor)
    mdl_mat <- model.matrix(as.formula(paste("~",paste0(variables, collapse=' + '))),
              newdata_, contrasts.arg = lapply(newdata_,contrasts, contrasts=FALSE))[,-1]
    if(!setequal(colnames(object$beta_matrix), colnames(mdl_mat))){
      inter_params <- intersect(colnames(object$beta_matrix), colnames(mdl_mat))
      beta_mat     <- object$beta_matrix[,inter_params]
      newdata.mat  <- mdl_mat[,inter_params]
    }else{
      beta_mat <- object$beta_matrix
      newdata.mat <- mdl_mat
    }
  }else if(any(class(newdata) == 'matrix')){
    newdata.mat <- newdata[,colnames(object$beta_matrix)]
  }else{
    stop('newdata must be one of matrix or data.frame!')
  }

  return(pnorm((newdata.mat %*% beta_mat['mean',])/(newdata.mat %*% beta_mat['variance',] + object$beta ** 2)))
}


#' Online learning with stream of samples.
#'
#' @param object \code{BOPR} object
#' @param newdata a data.frame of predictors.
#' @param allow.new allow new levels of current variable
#'
#' @return S3 \code{BOPR} object; a list of consisting of
#'  \item{beta_matrix}{beta matrix with mean and variance}
#'  \item{beta}{scaling parameter}
#'  \item{prior_prob}{prior of initial parametes}
#'  \item{epsilon}{parameter to apply dynamics}
#' @export
online_leraning <- function(object, newdata=NULL, allow.new=TRUE){
  if(!any(class(newdata) == 'data.frame')){
    stop('newdata must be data.frame!')
  }
  res <- .make_training_data(object$formula, data = newdata)
  return(online_leraning.matrix(object, x=res$x, y=res$y, allow.new=allow.new))
}

#' Online learning with stream of samples.
#'
#' @param object \code{BOPR} object
#' @param x a  matrix of predictors.
#' @param y a factor vector with 2 level
#' @param allow.new allow new levels of current variable
#'
#' @return S3 \code{BOPR} object; a list of consisting of
#'  \item{beta_matrix}{beta matrix with mean and variance}
#'  \item{beta}{scaling parameter}
#'  \item{prior_prob}{prior of initial parametes}
#'  \item{epsilon}{parameter to apply dynamics}
#' @export
online_leraning.matrix <- function(object, x = NULL, y=NULL, allow.new=TRUE){
  if(class(object) != 'BOPR'){
    stop('object is not BOPR class!')
  }
  y_lab <- ifelse(y == 1, 1, -1)

  #if there is other column in x
  if(!setequal(colnames(object$beta_matrix), colnames(x))){
    added_col <- setdiff(colnames(x), colnames(object$beta_matrix))
    print(sprintf('will added %s features!',length(added_col)))
    for(nm in added_col){
      bias_mean <- prior_mean(object$prior_prob, object$beta, ncol(object$beta_matrix) + 1)
      object$beta_matrix <- cbind(object$beta_matrix, c(mean=bias_mean, variance=pkg.env$PRIOR_VARIANCE))
      colnames(object$beta_matrix)[ncol(object$beta_matrix)] <- nm
    }
  }


  for(i in 1:nrow(x)){
    total_mean <- sum(object$beta_matrix[1,])
    total_variance <- sum(object$beta_matrix[2,]) + object$beta ** 2
    t <- y_lab[i] * total_mean / sqrt(total_variance)
    t <- clip(t, pkg.env$MIN_SURPRISE, pkg.env$MAX_SURPRISE )
    v <- dnorm(t) / pnorm(t)
    w <- v * (v + t)
    col_nm <- colnames(x)
    #if(length(col_nm) != ncol(object$beta_matrix))
    #  stop("num of features is not match.")
    for(j in col_nm){
      if(x[i,j] == 0) next
      mean_delta <- y_lab[i] * object$beta_matrix['variance',j] / sqrt(total_variance) * v
      variance_multiplier <- 1.0 - object$beta_matrix['variance',j] / total_variance * w
      updated <- c(
        mean = as.numeric(object$beta_matrix['mean',j] + mean_delta),
        variance = as.numeric(object$beta_matrix['variance',j] * variance_multiplier))
      #apply dynamics
      prior <- c(mean=0, variance=1)
      adjusted_variance <- updated['variance'] * prior['variance'] /
        ((1.0 - object$epsilon) * prior['variance'] +
           object$epsilon * updated['variance'])
      adjusted_mean <- adjusted_variance * (
        (1.0 - object$epsilon) * updated['mean'] / updated['variance'] +
          object$epsilon * prior['mean'] / prior['variance'])

      object$beta_matrix[,j] <- c(mean=adjusted_mean, variance=adjusted_variance)
    }

  }
  return(object)
}



