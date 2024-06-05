#' @title Spaghetti-plot fixed and random effects of linear mixed models
#' @name pastaPlot
#'
#' @description
#'   \code{pastaPlot()} plots slopes for both fixed and random effects of linear mixed models from 'lme4' or 'glmmTMB' packages as a single spaghetti plot, optionally between conditions including confidence bands for fixed effects.
#'
#' @param model lme4 or glmmTMB model object
#' @param predictor (Character) Name of predictor (e.g., "time" or "math_score"), as it is present in the model
#' @param nested.in (Character) Name of the variable your time points or subjects are nested in (e.g.,"school" or "id")
#' @param group (Optional, character) The name of your grouping variable (e.g., "condition" or "gender")
#' @param legend.title (Optional, character) Name of legend in plot (e.g., "Condition", or "Gender")
#' @param group.labels (Optional, vector of characters) Names of group labels to be displayed in the plot (e.g., c("Control", "Intervention"))
#' @param xlab (Optional, character) Label of x-axis (predictor) (e.g., "Time (days)")
#' @param ylab (Optional, character) Label of y-axis (dependant variable) (e.g., "GAF")
#' @param font.family (Optional, character) Name of the font family (e.g. "serif")
#' @param colors (Optional, vector of characters) Set color of slopes. Length of vector should correspond to number of values in group variable (e.g., c("#5e9aff", "blue")). If no group variable is specified, pass a single color.
#' @param ci.lvl (Optional, numeric) Set confidence interval (default: 0.95)
#' @param ci.int (Optional, logical) Enable confidence (prediction) intervals, disabled by default
#' @param ci.linetype (Optional, numeric) Set linetype of confidence bands outline (default: 0)
#' @param lwd.fix (Optional, numeric) Line width of fixed effects (default: 1)
#' @param lwd.ran (Optional, numeric) Line width of random effects (default: 0.5)
#' @param xlab.inc (Optional, numeric) Increment the displayed values of your predictor (e.g., xlab_int = 1 changes range of x from 0-29 to 1-30), set to 0 by default
#' @param xlab.int (Optional, numeric) Interval between displayed predictor values on x-axis (e.g., "1"), disabled by default
#' @param ylim (Optional, numeric vector) Limited range of values on y-axis (e.g. c(1,5.5))
#' @param opacity.ci (Optional, numeric) Set opacity of confidence bands in the range of 0 to 1 (default = 0.1)
#' @param opacity.ran (Optional, numeric) Set opacity of random slopes in the range of 0 to 1 (default = 0.4)
#' @param colors.ci (Optional, vector of characters) Set color of confidence bands. Length of vector should correspond to number of values in group variable (e.g., c("#5e9aff", "blue")). If no group variable is specified, pass a single color.
#' @return Returns a ggplot2 plot object to further be modified
#' @export
#'
#' @examples lme4_model <- lme4::lmer(CO2 ~ 1 + time*condition + (1 + time | id),
#' data=ecovia_data, REML = FALSE, control = lme4::lmerControl(optimizer = "bobyqa"))
#' @examples pastaPlot(lme4_model, "time", "id", group = "condition", legend.title = "Condition",
#' group.labels = c("Control", "Intervention"), ci.int = TRUE, xlab = "Time (days)",
#' ylab = "CO2")
#'
#' @examples glmmTMB_model <- glmmTMB::glmmTMB(math_score_y3 ~ 1 + math_score_y1*gender +
#' (1 + math_score_y1 | school), data=jsp_data, REML = FALSE)
#' @examples pastaPlot(glmmTMB_model, "math_score_y1", "school", group = "gender",
#' legend.title = "Gender", group.labels = c("Male", "Female"), ci.int = FALSE,
#' xlab = "Math score (year 1)", ylab = "Math score (year 3)")


pastaPlot <- function(model = NULL,
                      predictor = NULL,
                      nested.in = NULL,
                      group = NULL,
                      legend.title = "Legend",
                      group.labels = NULL,
                      xlab = NULL,
                      ylab = NULL,
                      font.family = NULL,
                      colors = NULL,
                      ci.lvl = .95,
                      ci.int = FALSE,
                      ci.linetype = 0,
                      lwd.fix = 1,
                      lwd.ran = 0.5,
                      xlab.inc = 0,
                      xlab.int = NULL,
                      ylim = NULL,
                      opacity.ci = 0.25,
                      opacity.ran = 0.3,
                      colors.ci = NULL) {

  # create data frames for fixed and random effects
  cookedPasta <- cookPasta(model = model, predictor = predictor, nested.in = nested.in, group = group, ci.int = ci.int, ci.lvl = ci.lvl)
  fix.ef <- cookedPasta[[1]]
  rand.ef <- cookedPasta[[2]]
  group.exists <- !(is.null(group))
  colors.exist <- !(is.null(colors))
  colors.ci.exist <- !(is.null(colors.ci))
  ngroups <- nrow(unique(cookedPasta[[1]]["group"]))

  model_frame <- model.frame(model)
  range_1 <- range(model_frame[[predictor]])
  range_2 <- range(model_frame[[group]])

  if(is.null(ylab)) ylab <- names(model_frame[1])
  if(is.null(xlab)) xlab <- predictor
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' required for this function to work, please install it.")
  }
  if(length(colors) > ngroups) stop("More colors set than values in grouping variable. Check your colors.")
  if(group.exists) {
    p1 <- ggplot2::ggplot() + ggplot2::geom_line(ggplot2::aes(x = rand.ef$predictor, y = rand.ef$pred, group = rand.ef$id, color = factor(rand.ef$group)), lwd = lwd.ran, alpha = opacity.ran)
    if(ci.int == TRUE) {
      p1 <- p1 + ggplot2::geom_ribbon(ggplot2::aes(x = fix.ef$predictor, y = fix.ef$pred, ymin = fix.ef$conf.low, ymax = fix.ef$conf.high, fill = factor(fix.ef$group)), linetype = ci.linetype, alpha = opacity.ci)
    }
    p1 <- p1 + ggplot2::geom_line(ggplot2::aes(x = fix.ef$predictor, y = fix.ef$pred, color = factor(fix.ef$group)), lwd = lwd.fix)
    if(colors.exist){
      if(length(colors) < ngroups) stop("Less colors passed to colors than levels of your grouping variable. Add some color (to your life) or remove the argument to get default colors!")
      p1 <- p1 + ggplot2::scale_color_manual(values = colors, labels = group.labels) +
        ggplot2::labs(color = legend.title)
      if(colors.ci.exist) {
        if(length(colors.ci) < ngroups) stop("Less colors passed to colors.ci than levels of your grouping variable. Add some color (to your life) or remove the argument to get default colors!")
        p1 <- p1 + ggplot2::scale_fill_manual(values = colors.ci, guide = "none")
      } else {
        p1 <- p1 + ggplot2::scale_fill_manual(values = colors, guide = "none")
      }
    } else {
      if(colors.ci.exist) {
        p1 <- p1 + ggplot2::scale_color_discrete(labels = group.labels) +
        ggplot2::scale_fill_manual(values = colors.ci, labels = group.labels, guide = "none") +
        ggplot2::labs(color = legend.title)
      } else {
        p1 <- p1 + ggplot2::scale_color_discrete(labels = group.labels) +
        ggplot2::scale_fill_discrete(guide = "none") +
        ggplot2::labs(color = legend.title)
      }
    }
  } else {
    p1 <- ggplot2::ggplot() + ggplot2::geom_line(ggplot2::aes(x = rand.ef$predictor, y = rand.ef$pred, group = rand.ef$id), lwd = lwd.ran, alpha = opacity.ran)
    if(ci.int == TRUE) {
      p1 <- p1 + ggplot2::geom_ribbon(ggplot2::aes(x = fix.ef$predictor, y = fix.ef$pred, ymin = fix.ef$conf.low, ymax = fix.ef$conf.high), linetype = ci.linetype, alpha = opacity.ci)
    }
    p1 <- p1 + ggplot2::geom_line(ggplot2::aes(x = fix.ef$predictor, y = fix.ef$pred), lwd = lwd.fix)
    if(colors.exist){
      p1 <- ggplot2::ggplot() + ggplot2::geom_line(ggplot2::aes(x = rand.ef$predictor, y = rand.ef$pred, group = rand.ef$id, color = colors), lwd = lwd.ran, alpha = opacity.ran)
      if(ci.int == TRUE) {
        p1 <- p1 + ggplot2::geom_ribbon(ggplot2::aes(x = fix.ef$predictor, y = fix.ef$pred, ymin = fix.ef$conf.low, ymax = fix.ef$conf.high, fill = colors), linetype = ci.linetype, alpha = opacity.ci)
      }
      p1 <- p1 + ggplot2::geom_line(ggplot2::aes(x = fix.ef$predictor, y = fix.ef$pred, color = colors), lwd = lwd.fix)
      p1 <- p1 + ggplot2::scale_color_manual(values = colors, guide = "none")
      if(colors.ci.exist) {
        p1 <- p1 + ggplot2::scale_fill_manual(values = colors.ci, guide = "none")
      } else {
        p1 <- p1 + ggplot2::scale_fill_manual(values = colors, guide = "none")
      }
    }
  }
  p1 <- p1 +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), plot.title.position = "plot") +
    ggplot2::xlab(xlab) +
    ggplot2::ylab(ylab)

  if(!(is.null(xlab.int))){
    p1 <- p1 + ggplot2::scale_x_continuous(breaks = seq(range_1[1], range_1[2], by = xlab.int), labels = seq((range_1[1]+xlab.inc), (range_1[2]+xlab.inc), by = xlab.int), limits = c(range_1[1], range_1[2]), expand = c(0, 0))
  } else {
    p1 <- p1 + ggplot2::scale_x_continuous(expand = c(0,0))
  }

  if(!(is.null(font.family))){
    p1 <- p1 + ggplot2::theme_minimal(base_family = font.family)
  }
  if(!(is.null(ylim))) {
    p1 <- p1 + ggplot2::ylim(ylim)
  }
  return(p1)
}

