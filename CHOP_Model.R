# CHOP-INTEND nonlinear model with 95% delta-method CI

library(tidyverse)
library(minpack.lm)
library(ggplot2)
library(MASS)   # for mvrnorm

# 1. Load data
setwd("~/Downloads")
chop <- read.csv("CHOP_Dataset.csv")

chop <- chop %>%
  mutate(
    Treatment_Sequence = factor(Treatment_Sequence),
    Change_Mean = as.numeric(Change_Mean),
    Followup_Time_Months = as.numeric(Followup_Time_Months)
  )

# 2. Indicator variables
chop <- chop %>%
  mutate(
    I_Z_N = ifelse(Treatment_Sequence == "Z>N", 1, 0),
    I_N_R = ifelse(Treatment_Sequence == "N>R", 1, 0),
    I_N_Z = ifelse(Treatment_Sequence == "N>Z", 1, 0)
  )

# 3. Fixed rate
k_fixed <- 0.25

# 4. Fit nonlinear model
chop_model <- nlsLM(
  Change_Mean ~
    (A_ZN * I_Z_N +
       A_NR * I_N_R +
       A_NZ * I_N_Z) *
    (1 - exp(-k_fixed * Followup_Time_Months)),
  
  data = chop,
  
  start = list(
    A_ZN = 6,
    A_NR = 6,
    A_NZ = 6
  ),
  
  control = nls.lm.control(maxiter = 1000)
)

coef(chop_model)
vc <- vcov(chop_model)

# 5. Prediction grid
time_grid <- seq(
  0,
  max(chop$Followup_Time_Months),
  length.out = 100
)

pred_grid <- expand.grid(
  Followup_Time_Months = time_grid,
  Treatment_Sequence = levels(chop$Treatment_Sequence)
)

pred_grid <- pred_grid %>%
  mutate(
    I_Z_N = ifelse(Treatment_Sequence == "Z>N", 1, 0),
    I_N_R = ifelse(Treatment_Sequence == "N>R", 1, 0),
    I_N_Z = ifelse(Treatment_Sequence == "N>Z", 1, 0)
  )

# 6. Parametric simulation for CI
set.seed(123)
n_sim <- 5000

theta_sim <- mvrnorm(
  n = n_sim,
  mu = coef(chop_model),
  Sigma = vc
)

pred_mat <- matrix(NA, nrow = nrow(pred_grid), ncol = n_sim)

for (i in seq_len(n_sim)) {
  
  A_ZN <- theta_sim[i, "A_ZN"]
  A_NR <- theta_sim[i, "A_NR"]
  A_NZ <- theta_sim[i, "A_NZ"]
  
  pred_mat[, i] <-
    (A_ZN * pred_grid$I_Z_N +
       A_NR * pred_grid$I_N_R +
       A_NZ * pred_grid$I_N_Z) *
    (1 - exp(-k_fixed * pred_grid$Followup_Time_Months))
}

pred_grid$fit <- rowMeans(pred_mat)
pred_grid$lower <- apply(pred_mat, 1, quantile, 0.025)
pred_grid$upper <- apply(pred_mat, 1, quantile, 0.975)

# 7. Plot (ALL sequences on SAME plot, ALL CIs visible)
ggplot(pred_grid, aes(
  x = Followup_Time_Months,
  y = fit,
  color = Treatment_Sequence,
  fill = Treatment_Sequence
)) +
  geom_ribbon(
    aes(ymin = lower, ymax = upper),
    alpha = 0.30,
    color = NA
  ) +
  geom_line(linewidth = 1.3) +
  labs(
    title = "CHOP-INTEND Improvement by Treatment Sequence",
    subtitle = "Asymptotic nonlinear model with 95% parametric confidence intervals",
    x = "Follow-up Time (Months)",
    y = "Change in CHOP-INTEND Score"
  ) +
  theme_minimal(base_size = 14)

ggsave(
  filename = "CHOP_model_95CI.png",
  plot = last_plot(),
  width = 8,
  height = 6,
  dpi = 300
)

# STEP 2: CHOP fixed-time predictions with 95% CI
library(dplyr)
library(MASS)   # for mvrnorm

# 1. Timepoints of interest
timepoints <- c(6, 12, 18, 24)

# 2. Prediction grid
step2_chop <- expand.grid(
  Followup_Time_Months = timepoints,
  Treatment_Sequence  = levels(chop$Treatment_Sequence)
)

step2_chop <- step2_chop %>%
  mutate(
    I_Z_N = ifelse(Treatment_Sequence == "Z>N", 1, 0),
    I_N_R = ifelse(Treatment_Sequence == "N>R", 1, 0),
    I_N_Z = ifelse(Treatment_Sequence == "N>Z", 1, 0)
  )

# 3. Point predictions
step2_chop$fit <-
  (coef(chop_model)["A_ZN"] * step2_chop$I_Z_N +
     coef(chop_model)["A_NR"] * step2_chop$I_N_R +
     coef(chop_model)["A_NZ"] * step2_chop$I_N_Z) *
  (1 - exp(-k_fixed * step2_chop$Followup_Time_Months))

# 4. Delta-method CI (parametric simulation)
set.seed(123)
n_sim <- 5000

theta_sim <- mvrnorm(
  n     = n_sim,
  mu    = coef(chop_model),
  Sigma = vcov(chop_model)
)

pred_mat <- matrix(NA, nrow = nrow(step2_chop), ncol = n_sim)

for (i in seq_len(n_sim)) {
  
  A_ZN <- theta_sim[i, "A_ZN"]
  A_NR <- theta_sim[i, "A_NR"]
  A_NZ <- theta_sim[i, "A_NZ"]
  
  pred_mat[, i] <-
    (A_ZN * step2_chop$I_Z_N +
       A_NR * step2_chop$I_N_R +
       A_NZ * step2_chop$I_N_Z) *
    (1 - exp(-k_fixed * step2_chop$Followup_Time_Months))
}

# 5. 95% CI extraction
step2_chop$lower <- apply(
  pred_mat, 1, quantile, probs = 0.025
)

step2_chop$upper <- apply(
  pred_mat, 1, quantile, probs = 0.975
)

step2_chop

# Find time when LOWER 95% CI crosses +4.0 CHOP-INTEND
threshold <- 4.0

# ---- Z>N ----
cross_time_ZN <- min(
  pred_grid$Followup_Time_Months[
    pred_grid$Treatment_Sequence == "Z>N" &
      pred_grid$lower >= threshold
  ],
  na.rm = TRUE
)

# ---- N>Z ----
cross_time_NZ <- min(
  pred_grid$Followup_Time_Months[
    pred_grid$Treatment_Sequence == "N>Z" &
      pred_grid$lower >= threshold
  ],
  na.rm = TRUE
)

cross_time_ZN
cross_time_NZ

library(tidyverse)
library(minpack.lm)
library(MASS)
library(grid)

# 0) Load + preprocess data
setwd("~/Downloads")
chop <- read.csv("CHOP_Dataset.csv")

chop <- chop %>%
  mutate(
    Treatment_Sequence = factor(Treatment_Sequence),
    Change_Mean = as.numeric(Change_Mean),
    Followup_Time_Months = as.numeric(Followup_Time_Months)
  ) %>%
  filter(!is.na(Change_Mean), !is.na(Followup_Time_Months), !is.na(Treatment_Sequence))

# 1) Indicator variables
chop <- chop %>%
  mutate(
    I_Z_N = ifelse(Treatment_Sequence == "Z>N", 1, 0),
    I_N_R = ifelse(Treatment_Sequence == "N>R", 1, 0),
    I_N_Z = ifelse(Treatment_Sequence == "N>Z", 1, 0)
  )

# 2) Fixed rate
k_fixed <- 0.25

# 3) Fit nonlinear model
chop_model <- nlsLM(
  Change_Mean ~
    (A_ZN * I_Z_N +
       A_NR * I_N_R +
       A_NZ * I_N_Z) *
    (1 - exp(-k_fixed * Followup_Time_Months)),
  data = chop,
  start = list(A_ZN = 6, A_NR = 6, A_NZ = 6),
  control = nls.lm.control(maxiter = 1000)
)

summary(chop_model)

# 4) Prediction grid (0–24 months)
time_grid <- seq(0, 24, length.out = 200)

pred_grid <- expand.grid(
  Followup_Time_Months = time_grid,
  Treatment_Sequence = levels(chop$Treatment_Sequence)
) %>%
  mutate(
    I_Z_N = ifelse(Treatment_Sequence == "Z>N", 1, 0),
    I_N_R = ifelse(Treatment_Sequence == "N>R", 1, 0),
    I_N_Z = ifelse(Treatment_Sequence == "N>Z", 1, 0)
  )

# 5) Parametric simulation for 95% uncertainty band
set.seed(123)
n_sim <- 5000

theta_sim <- mvrnorm(
  n = n_sim,
  mu = coef(chop_model),
  Sigma = vcov(chop_model)
)

pred_mat <- matrix(NA, nrow = nrow(pred_grid), ncol = n_sim)

for (i in seq_len(n_sim)) {
  A_ZN <- theta_sim[i, "A_ZN"]
  A_NR <- theta_sim[i, "A_NR"]
  A_NZ <- theta_sim[i, "A_NZ"]
  
  pred_mat[, i] <-
    (A_ZN * pred_grid$I_Z_N +
       A_NR * pred_grid$I_N_R +
       A_NZ * pred_grid$I_N_Z) *
    (1 - exp(-k_fixed * pred_grid$Followup_Time_Months))
}

pred_grid$fit   <- rowMeans(pred_mat)
pred_grid$lower <- apply(pred_mat, 1, quantile, probs = 0.025)
pred_grid$upper <- apply(pred_mat, 1, quantile, probs = 0.975)

# 6) Conservative time-to-+4 (lower 95% crossing)
threshold_pos <-  4.0
threshold_neg <- -4.0

cross_tbl <- pred_grid %>%
  group_by(Treatment_Sequence) %>%
  summarise(
    cross_time_lower = suppressWarnings(min(
      Followup_Time_Months[lower >= threshold_pos], na.rm = TRUE
    )),
    .groups = "drop"
  ) %>%
  mutate(
    cross_time_lower = ifelse(is.infinite(cross_time_lower), NA, cross_time_lower),
    label = ifelse(is.na(cross_time_lower), NA, paste0(round(cross_time_lower, 1), " mo"))
  )

print(cross_tbl)

y_top <- max(pred_grid$upper, na.rm = TRUE)
cross_tbl <- cross_tbl %>%
  mutate(
    label_y = y_top - 1.2,      
    label_x = cross_time_lower + 0.5
  )

x_right <- max(pred_grid$Followup_Time_Months, na.rm = TRUE)

p2 <- ggplot(pred_grid, aes(
  x = Followup_Time_Months,
  y = fit,
  color = Treatment_Sequence,
  fill = Treatment_Sequence
)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.25, color = NA) +
  geom_line(linewidth = 1.5) +
  
  geom_hline(yintercept = threshold_pos, linetype = "dashed", linewidth = 0.9) +
  geom_hline(yintercept = threshold_neg, linetype = "dashed", linewidth = 0.9) +
  
  geom_vline(
    data = cross_tbl,
    aes(xintercept = cross_time_lower, color = Treatment_Sequence),
    linetype = "dotted",
    linewidth = 1.5,
    show.legend = FALSE,
    na.rm = TRUE
  ) +
  
  geom_text(
    data = cross_tbl,
    aes(x = label_x, y = label_y, label = label, color = Treatment_Sequence),
    angle = 90,
    vjust = -0.2,
    size = 5.5,
    fontface = "bold",
    show.legend = FALSE,
    na.rm = TRUE
  ) +
  
  annotate(
    "text",
    x = x_right - 0.4, y = threshold_pos + 0.65,
    label = "Clinically meaningful (+4)",
    fontface = "bold",
    size = 4,
    hjust = 1
  ) +
  annotate(
    "text",
    x = x_right - 0.4, y = threshold_neg - 0.65,
    label = "Clinically meaningful (-4)",
    fontface = "bold",
    size = 4,
    hjust = 1
  ) +
  
  labs(
    title = "CHOP-INTEND Improvement by Treatment Sequence",
    subtitle = "Shaded band = 95% uncertainty; dotted lines = time to reach +4",
    x = "Months since treatment initiation",
    y = "Predicted change in CHOP-INTEND",
    color = NULL,
    fill  = NULL
  ) +
  
  theme_minimal(base_size = 18) +
  theme(
    plot.title = element_text(face = "bold", size = 24, margin = margin(b = 2)),
    plot.subtitle = element_text(size = 16, margin = margin(b = 2)),
    
    # Compact header: legend in top margin
    legend.position = c(0.60, 1.06),
    legend.direction = "horizontal",
    legend.justification = c(0.5, 1),
    legend.text = element_text(size = 14),
    
    plot.margin = margin(t = 6, r = 10, b = 8, l = 10)
  )

print(p2)

ggsave("ISEF_Fig2_CHOP_TimeToPlus4_FINAL.png", p2, width = 16, height = 6.5, dpi = 300)


library(dplyr)
library(ggplot2)


if (!("w" %in% names(chop))) {
  # Try to infer SD + N names commonly used
  # Change these if needed:
  sd_col <- dplyr::case_when(
    "Change_SD" %in% names(chop) ~ "Change_SD",
    "SD_Change" %in% names(chop) ~ "SD_Change",
    TRUE ~ NA_character_
  )
  n_col <- dplyr::case_when(
    "Sample_Size(N)" %in% names(chop) ~ "Sample_Size(N)",
    "Sample_Size.N." %in% names(chop) ~ "Sample_Size.N.",
    "N" %in% names(chop) ~ "N",
    TRUE ~ NA_character_
  )
  
  if (is.na(sd_col) || is.na(n_col)) {
    stop("CHOP weights not found. Please confirm CHOP SD and N column names so we can compute inverse-variance weights.")
  }
  
  chop <- chop %>%
    mutate(
      SE = .data[[sd_col]] / sqrt(.data[[n_col]]),
      w  = 1 / (SE^2)
    )
}

# Helper: Weighted RMSE
w_rmse <- function(y, yhat, w) {
  sqrt(sum(w * (y - yhat)^2, na.rm = TRUE) / sum(w, na.rm = TRUE))
}

# 1) Fit baseline models (weighted)
lin_fit_chop  <- lm(Change_Mean ~ 0 + Followup_Time_Months, data = chop, weights = w)
quad_fit_chop <- lm(Change_Mean ~ 0 + Followup_Time_Months + I(Followup_Time_Months^2),
                    data = chop, weights = w)

chop <- chop %>%
  mutate(
    pred_lin  = predict(lin_fit_chop,  newdata = chop),
    pred_quad = predict(quad_fit_chop, newdata = chop),
    pred_ml   = predict(chop_model,    newdata = chop)
  )

# 2) Compute weighted RMSE
rmse_lin  <- w_rmse(chop$Change_Mean, chop$pred_lin,  chop$w)
rmse_quad <- w_rmse(chop$Change_Mean, chop$pred_quad, chop$w)
rmse_ml   <- w_rmse(chop$Change_Mean, chop$pred_ml,   chop$w)

rmse_tbl_chop <- data.frame(
  Model = c("Linear baseline", "Quadratic baseline", "Biologically constrained ML"),
  Weighted_RMSE = c(rmse_lin, rmse_quad, rmse_ml)
) %>%
  mutate(
    # baselines first, ML last
    Model = factor(Model, levels = c("Linear baseline", "Quadratic baseline", "Biologically constrained ML")),
    Type  = ifelse(Model == "Biologically constrained ML", "ML", "Baseline"),
    label_val = sprintf("%.3f", Weighted_RMSE)
  )

print(rmse_tbl_chop)

# Save table
write.csv(rmse_tbl_chop, "CHOP_ModelComparison_RMSE_Table.csv", row.names = FALSE)

# 3) “× lower error” callout vs mean baseline
baseline_mean <- mean(c(rmse_lin, rmse_quad))
ratio <- baseline_mean / rmse_ml

# y-limit padding
ymax <- max(rmse_tbl_chop$Weighted_RMSE) * 1.22

# Colors (match your HFMSE style)
col_ml <- "turquoise4"
col_baseline <- "gray80"

# 4) Plot (poster-ready)
p_rmse_chop <- ggplot(rmse_tbl_chop, aes(x = Model, y = Weighted_RMSE, fill = Type)) +
  geom_col(width = 0.72) +
  
  geom_text(
    aes(label = label_val),
    vjust = -0.6,
    size = 6,
    fontface = "bold",
    color = "black"
  ) +
  
  scale_fill_manual(values = c("Baseline" = col_baseline, "ML" = col_ml)) +
  
  annotate(
    "text",
    x = 3,
    y = rmse_ml + 0.15,
    label = paste0("≈ ", round(ratio, 1), "× lower error"),
    fontface = "bold",
    size = 5.6,
    color = col_ml
  ) +
  
  labs(
    title = "CHOP-INTEND: Model vs Baseline Comparisons",
    subtitle = "Validation that the biologically constrained ML model reduces prediction error vs standard baselines",
    x = NULL,
    y = "Weighted RMSE (lower is better)",
    fill = NULL,
    caption = "Weights: inverse-variance (w = 1/SE²). Baselines constrained to Δ(0)=0."
  ) +
  
  coord_cartesian(ylim = c(0, ymax), clip = "off") +
  
  theme_minimal(base_size = 18) +
  theme(
    plot.title = element_text(face = "bold", size = 26, margin = margin(b = 6)),
    plot.subtitle = element_text(size = 15, margin = margin(b = 12)),
    axis.text.y = element_text(size = 18),
    axis.title = element_text(face = "bold", size = 18),
    legend.position = "none",
    panel.grid.minor = element_blank(),
    plot.margin = margin(t = 10, r = 30, b = 10, l = 10),
    plot.caption = element_text(size = 13, hjust = 0)
  )

print(p_rmse_chop)

ggsave(
  filename = "ISEF_Fig6_CHOP_ModelComparison_RMSE.png",
  plot = p_rmse_chop,
  width = 12,
  height = 6,
  dpi = 300
)

