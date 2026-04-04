#######################################################
##### Web appendix A: plots of PROMs sum scores  ######
#######################################################

library(dplyr)
library(ggplot2)
library(patchwork)

to_num <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

eq_items <- c(
  "Mobility",
  "SelfCare",
  "UsualActivities",
  "PainorDiscomfort",
  "Anxietyordepression"
)

koos_items <- c(
  "Gettingoutofbed",
  "Puttingonsocks",
  "Gettingupfromsitting",
  "Bendingtowardthefloor_pickinguponobjectfromtheground",
  "Twisting_pivotingontheinjuredknee",
  "Kneeling",
  "Squatting"
)

set.seed(1234)

ids_ok <- dat %>%
  distinct(ID, visit_idx) %>%
  filter(visit_idx %in% c(0, 1, 2)) %>%
  count(ID, name = "n_vis") %>%
  filter(n_vis == 3) %>%
  pull(ID)

ids100 <- sample(ids_ok, 100)

dat_plot <- dat %>%
  filter(ID %in% ids100, visit_idx %in% c(0, 1, 2)) %>%
  mutate(
    time_months = case_when(
      visit_idx == 0 ~ 0,
      visit_idx == 1 ~ 6,
      visit_idx == 2 ~ 12,
      TRUE ~ NA_real_
    ),
    eq5d_sum = rowSums(
      across(all_of(eq_items), ~ to_num(.x)),
      na.rm = FALSE
    ),
    koosps_sum = rowSums(
      across(all_of(koos_items), ~ to_num(.x)),
      na.rm = FALSE
    ),
    phs = Perceivedcurrenthealthstatus
  ) %>%
  arrange(ID, time_months)

base_family <- "serif"
base_size <- 24
title_size <- 22
axis_title_size <- 20
axis_text_size <- 18
plot_title_size <- 20

theme_fig <- theme_minimal(base_family = base_family, base_size = base_size) +
  theme(
    plot.title = element_text(size = title_size, hjust = 0),
    axis.title = element_text(size = axis_title_size),
    axis.text = element_text(size = axis_text_size, color = "grey30"),
    panel.grid.minor = element_line(color = "grey88", linewidth = 0.25),
    panel.grid.major = element_line(color = "grey82", linewidth = 0.35),
    plot.margin = margin(6, 12, 6, 6),
    aspect.ratio = 0.55
  )
make_panel <- function(data, y, y_lab, main_lab, y_breaks, y_limits) {
  summ <- data %>%
    group_by(time_months) %>%
    summarise(
      n = sum(!is.na(.data[[y]])),
      mean_y = mean(.data[[y]], na.rm = TRUE),
      sd_y = sd(.data[[y]], na.rm = TRUE),
      se_y = sd_y / sqrt(n),
      tcrit = qt(0.975, df = pmax(n - 1, 1)),
      ci_low = mean_y - tcrit * se_y,
      ci_high = mean_y + tcrit * se_y,
      .groups = "drop"
    )
  
  ggplot(data, aes(x = time_months, y = .data[[y]], group = ID)) +
    geom_line(
      linewidth = 0.25,
      alpha = 0.18,
      color = "black",
      linetype = "solid",
      na.rm = TRUE
    ) +
    geom_line(
      data = summ,
      aes(x = time_months, y = mean_y, group = 1),
      inherit.aes = FALSE,
      linewidth = 1.2,
      color = "#3569ff",
      linetype = "dashed"
    ) +
    geom_point(
      data = summ,
      aes(x = time_months, y = mean_y),
      inherit.aes = FALSE,
      shape = 16,
      size = 2.4,
      color = "#3569ff"
    ) +
    geom_errorbar(
      data = summ,
      aes(x = time_months, ymin = ci_low, ymax = ci_high),
      inherit.aes = FALSE,
      width = 0.35,
      linewidth = 0.5,
      color = "#3569ff"
    ) +
    scale_x_continuous(
      name = "Months",
      breaks = c(0, 6, 12),
      limits = c(-0.5, 12.5)
    ) +
    scale_y_continuous(
      name = y_lab,
      breaks = y_breaks,
      limits = y_limits
    ) +
    labs(
      title = paste0("", main_lab)
    ) +
    theme_fig
}

p_eq5d <- make_panel(
  dat_plot,
  y = "eq5d_sum",
  y_lab = "EQ-5D Sum (5 dimensions)",
  main_lab = "",
  y_breaks = c(0, 2.5, 5, 7.5, 10),
  y_limits = c(0, 10)
)

p_koos <- make_panel(
  dat_plot,
  y = "koosps_sum",
  y_lab = "KOOS-PS Sum (7 items)",
  main_lab = "",
  y_breaks = c(0, 7, 14, 21, 28),
  y_limits = c(0, 28)
)

p_phs <- make_panel(
  dat_plot,
  y = "phs",
  y_lab = "Perceived Health Status",
  main_lab = "",
  y_breaks = c(0, 25, 50, 75, 100),
  y_limits = c(0, 100)
)

p_final <- (p_eq5d / p_koos / p_phs) +
  plot_annotation(
    title = "",
    theme = theme(
      plot.title = element_text(
        family = base_family,
        size = plot_title_size,
        hjust = 0
      )
    )
  )

print(p_final)

tikzDevice::tikz("/home/nc/MEGA/BioPsychometrics/biomlatex/scores.tex", width=8, height=17)
print(p_final)
dev.off()
