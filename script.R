library(metafor)
library(metaviz)
library(readr)
library(ggplot2)
library(ggthemes)
library(patchwork)
library(TOSTER)
library(dplyr)
library(tidyr)
library(devtools)
library(ggrepel)
devtools::install_github("MathiasHarrer/dmetar")
library(dmetar)

ma_pipe_sei <-
  function(dat,
           true_effect,
           rep_upper,
           rep_lower,
           analysis_title,
           plot = TRUE) {
    dat <- read.csv(dat) # Read the file 
    res <- rma(yi,
               sei = sei,
               data = dat,
               method = "DL") # Perform a meta-analysis
    sunset <-  viz_sunset(
      x = dat[, c("yi", "sei")],
      contours = FALSE,
      true_effect = true_effect,
      power_contours = "continuous"
    ) # Create sunset plot
    dat[["power_observed"]] <-
      (1 - stats::pnorm(stats::qnorm(1 - 0.05 / 2) * dat[["sei"]],
                        abs(true_effect), dat[["sei"]])) +
      stats::pnorm(stats::qnorm(0.05 / 2) *
                     dat[["sei"]], abs(true_effect),
                   dat[["sei"]]) # Calculate power for each study
    dat[["power_e33"]] <-
      (1 - stats::pnorm(stats::qnorm(1 - 0.05 / 2) * dat[["sei"]],
                        abs(0.33), dat[["sei"]])) +
      stats::pnorm(stats::qnorm(0.05 / 2) *
                     dat[["sei"]], abs(0.33),
                   dat[["sei"]]) # Calculate power for each study for an effect of 0.33
    dat[["power_e66"]] <-
      (1 - stats::pnorm(stats::qnorm(1 - 0.05 / 2) * dat[["sei"]],
                        abs(0.66), dat[["sei"]])) +
      stats::pnorm(stats::qnorm(0.05 / 2) *
                     dat[["sei"]], abs(0.66),
                   dat[["sei"]]) # Calculate power for each study for an effect of 0.66
    power_median_dat <- data.frame(observed=numeric(1),e33=numeric(1),e66=numeric(1))
    power_median_dat[["observed"]] <- median(dat[["power_observed"]])
    power_median_dat[["e33"]] <- median(dat[["power_e33"]])
    power_median_dat[["e66"]] <- median(dat[["power_e66"]])
    power_median_dat <- as.data.frame(power_median_dat)
    power_median_dat <- mutate(power_median_dat, 
                               analysis = analysis_title) # Calculate median power
    rep_se <-
      ((rep_upper) - (rep_lower)) / (2 * 1.96) # Calculate SE for the summary effect
    sink("/dev/null")
    et <- TOSTmeta(
      ES = true_effect,
      se = rep_se,
      low_eqbound_d = -0.2,
      high_eqbound_d = 0.2,
      plot = TRUE
    ) # Perform an equivalence test
    et <- as.data.frame(et)
    et <- mutate(et, Analysis = analysis_title) # Calculate median power
    sink()
    if (plot == TRUE) {
      effect_size <- et$ES
      et <- as.data.frame(et)
      et_plot <-
        ggplot(et,
               aes(
                 x = "",
                 y = ES,
                 ymin = low_eqbound_d,
                 ymax = high_eqbound_d
               )) +
        geom_hline(aes(yintercept = 0), linetype = 'solid', size = 0.5) +
        geom_linerange(size = 5,
                       colour = "#00AFBB",
                       alpha = 0.5) +
        annotate(
          geom = "point",
          x = "",
          y = effect_size,
          color = "black",
          shape = 18,
          size = 4
        ) +
        theme(axis.text.x = element_text(size = 8)) +
        coord_flip() +
        theme_minimal() +
        theme(strip.text = element_text(
          size = 8,
          face = "bold",
          angle = 90
        )) +
        theme(strip.background = element_rect(colour = "black", fill = "white")) +
        theme(axis.title.y = element_text(face = "bold", size = 12)) +
        theme(axis.title.x = element_text(face = "bold", size = 12)) +
        theme(axis.title.y = element_blank()) +
        ylab(expression("Effect size"))
      et_plot <-
        et_plot + geom_linerange(aes(ymin = UL_CI_ZTEST, ymax = LL_CI_ZTEST), size = 0.5)
      et_plot <-
        et_plot + geom_linerange(aes(ymin = UL_CI_TOST, ymax = LL_CI_TOST),
                                 size = 1.5,
                                 colour = "black")
      et_plot # Creates an equivalence test plot
    }
    
    value <- list(
      res = res,
      dat = dat,
      power_median_dat = power_median_dat,
      sunset = sunset,
      et = et,
      et_plot = et_plot
    ) # Create a list of output objects
    attr(value, "class") <- "ma_pipe_sei"
    value
  }


combine_et <-
  function(dat){
    et_plot_combine <-
      ggplot(dat,
             aes(
               x = reorder(Analysis, -ES),
               y = ES,
               ymin = -0.6,
               ymax = 0.6
             )) +
      geom_rect(xmin = -Inf, ymin = -0.1, xmax = Inf, ymax = 0.1,
                fill = "#00AFBB", alpha = 0.2) +
      geom_rect(xmin = -Inf, ymin = -0.2, xmax = Inf, ymax = 0.2,
                fill = "#00AFBB", alpha = 0.1) +
      geom_rect(xmin = -Inf, ymin = -0.5, xmax = Inf, ymax = 0.5,
                fill = "#00AFBB", alpha = 0.05) +
      geom_point(color = "black",
                 shape = 18,
                 size = 4) +
      theme(axis.text.x = element_text(size = 8)) +
      coord_flip() +
      theme_minimal() +
      theme(strip.text = element_text(
        size = 8,
        face = "bold",
        angle = 90
      )) +
      theme(strip.background = element_rect(colour = "black", fill = "white")) +
      theme(axis.title.y = element_text(face = "bold", size = 12)) +
      theme(axis.title.x = element_text(face = "bold", size = 12)) +
      theme(axis.title.y = element_blank()) +
      ylab(expression("Effect size"))
    et_plot_combine <-
      et_plot_combine + geom_linerange(aes(ymin = UL_CI_ZTEST, ymax = LL_CI_ZTEST), size = 0.5)
    et_plot_combine <-
      et_plot_combine + geom_linerange(aes(ymin = UL_CI_TOST, ymax = LL_CI_TOST),
                                       size = 1.5,
                                       colour = "black")
    et_plot_combine <-
      et_plot_combine + scale_y_continuous(breaks=c(-.5, -.20,-.1, 0, .1, 0.2, 0.5))
    et_plot_combine <-
      et_plot_combine + geom_hline(yintercept=0, linetype="dashed", color = "black")
    et_plot_combine # Creates an equivalence test plot
  }

leppanen_2017_1_results <- ma_pipe_sei(
  "dat_leppanen_2017_1.csv",
  true_effect = 0.09,
  rep_lower = -0.06,
  rep_upper = 0.24,
  analysis_title = "Leppanen_2017_1",
  plot = TRUE
)

leppanen_2017_1_power_median_dat <- leppanen_2017_1_results$power_median_dat
leppanen_2017_1_power_median_dat$ES <- 0.09

leppanen_2017_1_et_dat <- leppanen_2017_1_results$et

# One outlier was removed from the dataset in the original analysis
dat_leppanen_2017_2 <- 
  read.csv("dat_leppanen_2017_2.csv") # Load raw data
dat_leppanen_2017_2 <- 
  dat_leppanen_2017_2[!(dat_leppanen_2017_2$Trial %in% c(16)), ] # Remove outlier
write.csv(dat_leppanen_2017_2, 
          "dat_leppanen_2017_2_outlier_removed.csv") # Create new .csv file

leppanen_2017_2_results <- ma_pipe_sei(
  "dat_leppanen_2017_2_outlier_removed.csv",
  true_effect = 0.18,
  rep_lower = 0.06,
  rep_upper = 0.29,
  analysis_title = "Leppanen_2017_2",
  plot = TRUE
)

leppanen_2017_2_power_median_dat <- leppanen_2017_2_results$power_median_dat
leppanen_2017_2_power_median_dat$ES <- 0.18


leppanen_2017_2_et_dat <- leppanen_2017_2_results$et

leppanen_2017_3_results <- ma_pipe_sei(
  "dat_leppanen_2017_3.csv",
  true_effect = 0.05,
  rep_lower = -0.05,
  rep_upper = 0.15,
  analysis_title = "Leppanen_2017_3",
  plot = TRUE
)

leppanen_2017_3_power_median_dat <- leppanen_2017_3_results$power_median_dat
leppanen_2017_3_power_median_dat$ES <- 0.05

leppanen_2017_3_et_dat <- leppanen_2017_3_results$et

leppanen_2017_4_results <- ma_pipe_sei(
  "dat_leppanen_2017_4.csv",
  true_effect = 0.21,
  rep_lower = 0.07,
  rep_upper = 0.34,
  analysis_title = "Leppanen_2017_4",
  plot = TRUE
)

leppanen_2017_4_power_median_dat <- leppanen_2017_4_results$power_median_dat
leppanen_2017_4_power_median_dat$ES <- 0.21

leppanen_2017_4_et_dat <- leppanen_2017_4_results$et

leppanen_2017_5_results <- ma_pipe_sei(
  "dat_leppanen_2017_5.csv",
  true_effect = 0.18,
  rep_lower = -0.02,
  rep_upper = 0.39,
  analysis_title = "Leppanen_2017_5",
  plot = TRUE
)

leppanen_2017_5_power_median_dat <- leppanen_2017_5_results$power_median_dat
leppanen_2017_5_power_median_dat$ES <- 0.18

leppanen_2017_5_et_dat <- leppanen_2017_5_results$et

leppanen_2017_6_results <- ma_pipe_sei(
  "dat_leppanen_2017_6.csv",
  true_effect = 0.04,
  rep_lower = -0.1,
  rep_upper = 0.17,
  analysis_title = "Leppanen_2017_6",
  plot = TRUE
)

leppanen_2017_6_power_median_dat <- leppanen_2017_6_results$power_median_dat
leppanen_2017_6_power_median_dat$ES <- 0.04


leppanen_2017_6_et_dat <- leppanen_2017_6_results$et


com1 <- rbind(leppanen_2017_1_et_dat,
              leppanen_2017_2_et_dat,
              leppanen_2017_3_et_dat,
              leppanen_2017_4_et_dat,
              leppanen_2017_5_et_dat, 
              leppanen_2017_6_et_dat)
com1 <- combine_et(com1)
com1 <- com1 + ggtitle("Summary effect sizes and equivalence bounds") +
  theme(plot.title = element_text(hjust = 0.5))

power_med <- rbind(leppanen_2017_1_power_median_dat,
                   leppanen_2017_2_power_median_dat,
                   leppanen_2017_3_power_median_dat,
                   leppanen_2017_4_power_median_dat,
                   leppanen_2017_5_power_median_dat,
                   leppanen_2017_6_power_median_dat)

power_med_long <- gather(power_med, effect, power, observed:e66, factor_key=TRUE)

tile <- ggplot(data = power_med_long, aes(x = effect, y = reorder(analysis, -ES),)) + 
  geom_tile(aes(fill = power)) + 
  coord_equal(ratio = 0.8) + 
  scale_fill_gradient(name = "Power", low = "#FDF0FF", high = "#E200FD") + theme_tufte(base_family="Helvetica")

tile <- tile + theme(axis.text.y = element_blank(),
                     axis.ticks.y = element_blank(),
                     axis.title.y = element_blank())

tile <- tile + ggtitle("Median statistical power") +
  xlab("Effect size") +
  theme(plot.title = element_text(hjust = 0.5))

tile <- tile + scale_x_discrete(labels=c("observed" = "Observed \n effect", "e33" = "0.33",
                                         "e66" = "0.66"))
tile
com1 + tile