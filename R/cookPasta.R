#' @title Calculate predicted values for fixed and random effects
#' @name cookPasta
#'
#' @description
#'   \code{cookPasta()} creates dataframes from fixed and random effects of 'lme4' or 'glmmTMB' models (e.g., for plotting)
#'
#' @param model lme4 or glmmTMB model object
#' @param predictor (Character) Name of predictor (e.g., "time" or "math_score"), as it is present in the model
#' @param nested.in (Character) Name of the variable your time points or subjects are nested in (e.g.,"school" or "id")
#' @param group (Optional, character) The name of your grouping variable (e.g., "condition" or "gender")
#' @param ci.lvl (Optional, numeric) Set level of confidence (prediction) intervals (default: 0.95). Requires ci.int to be set to TRUE
#' @param ci.int (Optional, boolean) Enable confidence (prediction) intervals, disabled by default
#' @return Returns a list of two dataframes, in which the first element is the fixed effects dataframe and the second element the random effects dataframe

cookPasta <- function(model = NULL,
                      predictor = NULL,
                      nested.in = NULL,
                      group = NULL,
                      ci.int = FALSE,
                      ci.lvl = ci.lvl
) {
  # check model type
  if (inherits(model, "glmmTMB")) {
    if (!requireNamespace("glmmTMB", quietly = TRUE)) {
      stop("Package 'glmmTMB' required for this function to work, please install it.")
    }

    rand.ef <- glmmTMB::ranef(model)[[1]][[1]]
    fix.ef <- as.data.frame(lme4::fixef(model)[[1]])
  } else if (inherits(model, "lmerMod")) {
    if (!requireNamespace("lme4", quietly = TRUE)) {
      stop("Package 'lme4' required for this function to work, please install it.")
    }
    rand.ef <- lme4::ranef(model)[[1]]
    fix.ef <- as.data.frame(lme4::fixef(model))
  }
  else {
    stop("Model needs to be of class 'lmerMod' or 'glmmTMB'.")
  }
  if(length(rand.ef) > 2) stop("More than 2 random effects. Can only plot one random intercept and/or one random slope")
  # check arguments
  if(is.null(predictor)) stop("Please specify the name of your predictor as character corresponding to the name in the passed model.")
  # extract the model data
  model_frame <- model.frame(model)

  # # create dataframes for fixed and random effects
  if (!(predictor %in% names(model_frame))) stop("Predictor does not exist in the model.")
  if (is.double(model_frame[[predictor]])) {
    range_1 <- range(model_frame[[predictor]])
    vector_1 <- range_1[1]:range_1[2]
  } else {
    stop("Predictor is not double/numeric. Must be double/numeric.")
  }
  condition <- is.null(group)
  if (!condition) {
    if (!(group %in% names(model_frame))) stop("Grouping variable does not exist in the model.")
    if (is.double(model_frame[[group]])) {
      range_2 <- range(model_frame[[group]])
      vector_2 <- range_2[1]:range_2[2]
      df_fixef <- expand.grid(vector_1, vector_2)
      pred_names <- c("predictor", "group")
      names(df_fixef) <- pred_names
    } else {
      stop("Grouping variable is not numeric. Must be numeric.")
    }
  } else {
    df_fixef <- as.data.frame(vector_1)
    pred_names <- c("predictor")
    names(df_fixef) <- pred_names
  }

  df_conf <- as.data.frame(ggeffects::ggpredict(model, terms = c(paste(predictor,sep = " ", "[all]"), group), ci_level = ci.lvl))
  names(df_conf)[names(df_conf) == "x"] <- "predictor"
  names(df_conf)[names(df_conf) == "predicted"] <- "pred"
  df_fixef <- merge(df_fixef, df_conf, by = pred_names)

  df <- expand.grid(unique(rownames(rand.ef)), vector_1)
  names(df) <- c(nested.in, "predictor")
  if (!(condition)) {
    df_2 <- model_frame[!duplicated(model_frame[nested.in]), c(nested.in, group)]
    names(df_2) <- c(nested.in, "group")
    df <- merge(df, df_2, by = nested.in)
  } else {
    df_2 <- as.data.frame(model_frame[!duplicated(model_frame[nested.in]), c(nested.in)])
    names(df_2) <- c(nested.in)
  }
  df_ran <- as.data.frame(rand.ef)
  df_ran[[nested.in]] <- rownames(df_ran)
  df_ran <- merge(df, df_ran, by = nested.in)
  names(df_ran)[names(df_ran) == nested.in] <- "id"
  if("(Intercept)" %in% colnames(df_ran)) {
    df_ran[["pred"]] <- (fix.ef["(Intercept)", 1] + df_ran[["(Intercept)"]])
  } else {
    df_ran[["pred"]] <- fix.ef["(Intercept)", 1]
  }

  if(predictor %in% colnames(df_ran)) {
    df_ran[["pred"]] <- df_ran[["pred"]] + (fix.ef[predictor, 1] + df_ran[[predictor]])*df_ran[["predictor"]] +
      ifelse(!condition, ((fix.ef[group, 1]*df_ran[["group"]])), 0)
  } else {
    df_ran[["pred"]] <- df_ran[["pred"]] + (fix.ef[predictor, 1])*df_ran[["predictor"]] +
      ifelse(!condition, ((fix.ef[group, 1]*df_ran[["group"]])), 0)
  }
  return(list(df_fixef, df_ran))
}

