
# 1. Install & load packages
library(minpack.lm)
library(dplyr)
library(ggplot2)
# 2. Load and preprocess data
setwd("~/Downloads")
hfmse_raw <- read.csv("HFMSE_Dataset.csv")

hfmse <- hfmse_raw %>%
  rename(
    Time       = Followup_Time_Months,
    Change     = Change_Mean,
    Baseline   = Baseline_Mean,
    Sequence   = Treatment_Sequence,
    SD_Change  = Change_SD,
    N          = Sample_Size.N.
  ) %>%
  filter(!is.na(Time), !is.na(Change),
         !is.na(SD_Change), !is.na(N)) %>%
  group_by(Sequence) %>%
  filter(n() >= 2) %>%  
  ungroup()


hfmse <- hfmse %>%
  mutate(
    A_N_R = ifelse(Sequence == "N>R", 1, 0),
    A_Z_R = ifelse(Sequence == "Z>R", 1, 0)
  )

# 3. Inverse-variance weights
hfmse <- hfmse %>%
  mutate(
    SE = SD_Change / sqrt(N),
    w  = 1 / (SE^2)
  )

# 4. Fit weighted nonlinear model (Δ = 0 at Time = 0)
asym_model <- nlsLM(
  Change ~ 
    (
      (A_N * A_N_R + A_Z * A_Z_R) +
        b * Baseline
    ) * (1 - exp(-k * Time)),
  
  data    = hfmse,
  weights = w,
  
  start = list(
    A_N = -3,
    A_Z =  5,
    b   =  0.02,
    k   =  0.15
  ),
  
  control = nls.lm.control(maxiter = 1000)
)

summary(asym_model)

# 5. Prediction grid
newdat <- expand.grid(
  Time = seq(0, max(hfmse$Time), length.out = 100),
  Sequence = c("N>R", "Z>R"),
  Baseline = 6.6
)

newdat <- newdat %>%
  mutate(
    A_N_R = ifelse(Sequence == "N>R", 1, 0),
    A_Z_R = ifelse(Sequence == "Z>R", 1, 0)
  )

newdat$fit <- predict(asym_model, newdat)

# 6. Bootstrapping (weights recomputed each time)
set.seed(123)

n_boot <- 1000
boot_preds <- matrix(NA, nrow = nrow(newdat), ncol = n_boot)

for (i in seq_len(n_boot)) {
  
  idx <- sample(seq_len(nrow(hfmse)), replace = TRUE)
  boot_data <- hfmse[idx, ]
  
  boot_data <- boot_data %>%
    mutate(
      SE = SD_Change / sqrt(N),
      w  = 1 / (SE^2)
    )
  
  boot_fit <- try(
    nlsLM(
      Change ~ 
        (
          (A_N * A_N_R + A_Z * A_Z_R) +
            b * Baseline
        ) * (1 - exp(-k * Time)),
      
      data    = boot_data,
      weights = w,
      start   = coef(asym_model),
      control = nls.lm.control(maxiter = 500)
    ),
    silent = TRUE
  )
  
  if (!inherits(boot_fit, "try-error")) {
    boot_preds[, i] <- predict(boot_fit, newdat)
  }
}

# 7. 95% bootstrap prediction intervals
newdat$lower <- apply(
  boot_preds, 1, quantile, probs = 0.025, na.rm = TRUE
)

newdat$upper <- apply(
  boot_preds, 1, quantile, probs = 0.975, na.rm = TRUE
)

# 8. Plot
ggplot(newdat, aes(Time, fit, color = Sequence, fill = Sequence)) +
  geom_ribbon(
    aes(ymin = lower, ymax = upper),
    alpha = 0.25,
    color = NA
  ) +
  geom_line(linewidth = 1.3) +
  labs(
    title = "Nonlinear Prediction of HFMSE Change Over Time",
    subtitle = "Inverse-variance weighted model with correct bootstrap uncertainty",
    x = "Months Since Treatment Initiation",
    y = "Predicted Change in HFMSE"
  ) +
  theme_minimal()

coef(asym_model)

median_baseline <- median(hfmse$Baseline, na.rm = TRUE)
median_baseline

names(coef(asym_model))


# STEP 2: HFMSE predictions + bootstrap CI at fixed times

library(dplyr)

baseline_fixed <- median(hfmse$Baseline, na.rm = TRUE)

timepoints <- c(6, 12, 18, 24)

step2_grid <- expand.grid(
  Time = timepoints,
  Sequence = c("N>R", "Z>R"),
  Baseline = baseline_fixed
)

step2_grid <- step2_grid %>%
  mutate(
    A_N_R = ifelse(Sequence == "N>R", 1, 0),
    A_Z_R = ifelse(Sequence == "Z>R", 1, 0)
  )

step2_grid$fit <- predict(asym_model, step2_grid)

set.seed(123)
n_boot <- 1000

boot_mat <- matrix(NA, nrow = nrow(step2_grid), ncol = n_boot)

for (i in seq_len(n_boot)) {
  
  idx <- sample(seq_len(nrow(hfmse)), replace = TRUE)
  boot_data <- hfmse[idx, ]
  
  boot_data <- boot_data %>%
    mutate(
      SE = SD_Change / sqrt(N),
      w  = 1 / (SE^2)
    )
  
  boot_fit <- try(
    nlsLM(
      Change ~ 
        (
          (A_N * A_N_R + A_Z * A_Z_R) +
            b * Baseline
        ) * (1 - exp(-k * Time)),
      
      data    = boot_data,
      weights = w,
      start   = coef(asym_model),
      control = nls.lm.control(maxiter = 500)
    ),
    silent = TRUE
  )
  
  if (!inherits(boot_fit, "try-error")) {
    boot_mat[, i] <- predict(boot_fit, step2_grid)
  }
}

# 95% prediction intervals
step2_grid$lower <- apply(
  boot_mat, 1, quantile, probs = 0.025, na.rm = TRUE
)

step2_grid$upper <- apply(
  boot_mat, 1, quantile, probs = 0.975, na.rm = TRUE
)

step2_grid

# Find time when LOWER CI crosses +3.0 (Z>R only)

baseline_val <- 15.1

newdat_ZR <- data.frame(
  Time = seq(0, 24, by = 0.01),
  Sequence = "Z>R",
  Baseline = baseline_val
) %>%
  mutate(
    A_N_R = 0,
    A_Z_R = 1
  )

newdat_ZR$fit <- predict(asym_model, newdat_ZR)

# Bootstrap predictions for Z>R

boot_preds_ZR <- matrix(NA, nrow = nrow(newdat_ZR), ncol = n_boot)

for (i in seq_len(n_boot)) {
  
  idx <- sample(seq_len(nrow(hfmse)), replace = TRUE)
  boot_data <- hfmse[idx, ] %>%
    mutate(
      SE = SD_Change / sqrt(N),
      w  = 1 / (SE^2)
    )
  
  boot_fit <- try(
    nlsLM(
      Change ~ 
        (
          (A_N * A_N_R + A_Z * A_Z_R) +
            b * Baseline
        ) * (1 - exp(-k * Time)),
      
      data    = boot_data,
      weights = w,
      start   = coef(asym_model),
      control = nls.lm.control(maxiter = 500)
    ),
    silent = TRUE
  )
  
  if (!inherits(boot_fit, "try-error")) {
    boot_preds_ZR[, i] <- predict(boot_fit, newdat_ZR)
  }
}

# Lower 95% CI
newdat_ZR$lower <- apply(
  boot_preds_ZR, 1, quantile, probs = 0.025, na.rm = TRUE
)

# Time when LOWER CI crosses +3.0
cross_time_ZR <- min(
  newdat_ZR$Time[newdat_ZR$lower >= 3.0],
  na.rm = TRUE
)

cross_time_ZR

# Minimum baseline HFMSE guaranteeing +3.0 improvement
# Z>R sequence only

threshold <- 3.0
time_grid <- seq(0, 12, by = 0.05)

# Plausible baseline range (conservative)
baseline_grid <- seq(0, 40, by = 0.1)

min_baseline <- NA

for (b0 in baseline_grid) {
  
  # Prediction grid for this baseline
  newdat <- data.frame(
    Time = time_grid,
    Sequence = "Z>R",
    Baseline = b0,
    A_N_R = 0,
    A_Z_R = 1
  )
  
  # Bootstrap predictions
  boot_preds <- matrix(NA, nrow = nrow(newdat), ncol = n_boot)
  
  for (i in seq_len(n_boot)) {
    
    idx <- sample(seq_len(nrow(hfmse)), replace = TRUE)
    boot_data <- hfmse[idx, ] %>%
      mutate(
        SE = SD_Change / sqrt(N),
        w  = 1 / (SE^2)
      )
    
    boot_fit <- try(
      nlsLM(
        Change ~ 
          (
            (A_N * A_N_R + A_Z * A_Z_R) +
              b * Baseline
          ) * (1 - exp(-k * Time)),
        
        data    = boot_data,
        weights = w,
        start   = coef(asym_model),
        control = nls.lm.control(maxiter = 300)
      ),
      silent = TRUE
    )
    
    if (!inherits(boot_fit, "try-error")) {
      boot_preds[, i] <- predict(boot_fit, newdat)
    }
  }
  
  # Lower CI across time
  lower_ci <- apply(
    boot_preds, 1, quantile, probs = 0.025, na.rm = TRUE
  )
  
  # Check if lower CI ever exceeds +3.0
  if (any(lower_ci >= threshold, na.rm = TRUE)) {
    min_baseline <- b0
    break
  }
}

min_baseline

# Minimum baseline guaranteeing NO clinically meaningful
# decline (lower CI > -3.0) for N>R — HFMSE

threshold <- -3.0
time_grid <- seq(0, 24, by = 0.05)

# Plausible baseline range
baseline_grid <- seq(0, 40, by = 0.1)

min_baseline_safe <- NA

for (b0 in baseline_grid) {
  
  # Prediction grid for this baseline (N>R only)
  newdat <- data.frame(
    Time = time_grid,
    Sequence = "N>R",
    Baseline = b0,
    A_N_R = 1,
    A_Z_R = 0
  )
  
  boot_preds <- matrix(NA, nrow = nrow(newdat), ncol = n_boot)
  
  for (i in seq_len(n_boot)) {
    
    idx <- sample(seq_len(nrow(hfmse)), replace = TRUE)
    boot_data <- hfmse[idx, ] %>%
      mutate(
        SE = SD_Change / sqrt(N),
        w  = 1 / (SE^2)
      )
    
    boot_fit <- try(
      nlsLM(
        Change ~ 
          (
            (A_N * A_N_R + A_Z * A_Z_R) +
              b * Baseline
          ) * (1 - exp(-k * Time)),
        
        data    = boot_data,
        weights = w,
        start   = coef(asym_model),
        control = nls.lm.control(maxiter = 300)
      ),
      silent = TRUE
    )
    
    if (!inherits(boot_fit, "try-error")) {
      boot_preds[, i] <- predict(boot_fit, newdat)
    }
  }
  
  # Lower CI across time
  lower_ci <- apply(
    boot_preds, 1, quantile, probs = 0.025, na.rm = TRUE
  )
  
  # Check if lower CI ALWAYS stays above -3.0
  if (all(lower_ci > threshold, na.rm = TRUE)) {
    min_baseline_safe <- b0
    break
  }
}

min_baseline_safe


# TIER 1 RESULT 2A — HFMSE baseline dependence

# Baseline grid to evaluate
baseline_grid <- seq(0, 40, by = 0.1)

# Fixed time horizon
T_eval <- 24

# Prediction dataframe
base_pred <- expand.grid(
  Baseline = baseline_grid,
  Sequence = c("N>R", "Z>R")
)

base_pred <- base_pred %>%
  mutate(
    Time = T_eval,
    A_N_R = ifelse(Sequence == "N>R", 1, 0),
    A_Z_R = ifelse(Sequence == "Z>R", 1, 0)
  )

# Bootstrap predictions at 24 months
boot_preds_base <- matrix(NA, nrow = nrow(base_pred), ncol = n_boot)

for (i in seq_len(n_boot)) {
  if (all(is.na(boot_preds[, i]))) next
  
  boot_preds_base[, i] <- boot_preds[, i][
    match(
      paste(base_pred$Time, base_pred$Sequence),
      paste(newdat$Time, newdat$Sequence)
    )
  ]
}

# Compute CI
base_pred$lower <- apply(
  boot_preds_base, 1, quantile, probs = 0.025, na.rm = TRUE
)

# Find minimum baseline for stability
baseline_threshold_hfmse <- base_pred %>%
  filter(lower >= -3) %>%
  group_by(Sequence) %>%
  summarise(
    min_baseline = min(Baseline)
  )

baseline_threshold_hfmse

# ----------------------------
# Additional Plot
# ----------------------------

# Force facet order
newdat1$BaselineGroup <- factor(
  newdat1$BaselineGroup,
  levels = c(lbl_low, lbl_med, lbl_high)
)

# ADDITION: To compite the time-to-+3 for Z→R (per baseline facet)
# Conservative: first time Z→R LOWER 95% band reaches +3  

threshold_pos <- 3.0

time_to_plus3_tbl <- newdat1 %>%
  filter(Sequence == "Z>R") %>%
  group_by(BaselineGroup) %>%
  summarise(
    time_to_plus3 = suppressWarnings(min(Time[lower >= threshold_pos], na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    time_to_plus3 = ifelse(is.infinite(time_to_plus3), NA, time_to_plus3),
    label = ifelse(is.na(time_to_plus3), NA, paste0(round(time_to_plus3, 1), " mo")),
    label_x = time_to_plus3 + 1.3,
    label_y = 11.2
  )

# Plotting again
p1 <- ggplot(newdat1, aes(x = Time, y = fit, color = Sequence, fill = Sequence)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.25, color = NA) +
  geom_line(linewidth = 1.5) +
  geom_hline(yintercept =  3, linetype = "dashed", linewidth = 0.9) +
  geom_hline(yintercept = -3, linetype = "dashed", linewidth = 0.9) +
  
geom_vline(
  data = time_to_plus3_tbl,
  aes(xintercept = time_to_plus3),
  linetype = "dotted",
  linewidth = 1.0,
  inherit.aes = FALSE,
  na.rm = TRUE
) +

geom_text(
  data = time_to_plus3_tbl,
  aes(x = label_x, y = label_y, label = label),
  angle = 90,
  vjust = -0.2,
  size = 5.2,
  fontface = "bold",
  inherit.aes = FALSE,
  color = "turquoise4",
  na.rm = TRUE
) +
  
  facet_wrap(~ BaselineGroup, nrow = 1) +
  coord_cartesian(ylim = c(-4.5, 12)) +
  
  annotate("text", x = 18, y =  3.55, label = "Clinically meaningful (+3)", fontface = "bold", size = 4) +
  annotate("text", x = 18, y = -3.55, label = "Clinically meaningful (-3)", fontface = "bold", size = 4) +
  annotate("text", x = 0.3, y = 11.2, label = "Z→R: Improves", fontface = "bold", size = 4.5, hjust = 0) +
  annotate("text", x = 0.3, y = -1.5, label = "N→R: Stable",  fontface = "bold", size = 4.5, hjust = 0) +
  
  labs(
    title = "HFMSE Change After Sequential Therapy",
    subtitle = "95% uncertainty band; dotted lines = time to reach +3 (Z→R)",
    x = "Months since treatment initiation",
    y = "Predicted change in HFMSE",
    color = NULL,
    fill  = NULL
  ) +
  
  theme_minimal(base_size = 18) +
  
  theme(
    plot.title = element_text(face = "bold", size = 22, margin = margin(b = 2)),
    plot.subtitle = element_text(size = 14, margin = margin(b = 2)),
    plot.margin = margin(t = 4, r = 10, b = 8, l = 10),
    
    legend.position = c(0.52, 1.18),
    legend.direction = "horizontal",
    legend.justification = c(0.5, 1),
    legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0),
    legend.text = element_text(size = 18),
    
    strip.text = element_text(face = "bold", size = 14),
    
    plot.title.position = "plot",
    panel.spacing.x = unit(1.18, "cm")
  )

print(p1)

ggsave("ISEF_Fig1_HFMSE_FINAL_compactHeader.png", p1, width = 16, height = 6.5, dpi = 300)

library(ggplot2)
library(dplyr)
library(grid)

USE_COLOR <- TRUE

col_base_dark  <- "#1f77b4"
col_base_light <- "#9ecae1"

col_seq_dark   <- "#2ca02c"
col_seq_light  <- "#98df8a"

col_black <- "black"

# ----------------------------
# 1) Demo parameters (testing shape; not from data)
# ----------------------------
k_demo   <- 0.20
A_high   <- 7.0
A_low    <- 3.5
b_demo   <- 0.06
B_low    <- 10
B_highB  <- 25

tgrid <- seq(0, 24, by = 0.1)

curve_fun <- function(t, A, k) A * (1 - exp(-k * t))
curve_fun_base <- function(t, A, b, Baseline, k) (A + b * Baseline) * (1 - exp(-k * t))


left_df <- data.frame(
  t = tgrid,
  y = curve_fun(tgrid, A = A_high, k = k_demo)
)

t_half <- log(2) / k_demo
y_half <- curve_fun(t_half, A_high, k_demo)

right_base <- bind_rows(
  data.frame(
    t = tgrid,
    y = curve_fun_base(tgrid, A = A_low, b = b_demo, Baseline = B_low, k = k_demo),
    line = paste0("Lower baseline (", B_low, ")")
  ),
  data.frame(
    t = tgrid,
    y = curve_fun_base(tgrid, A = A_low, b = b_demo, Baseline = B_highB, k = k_demo),
    line = paste0("Higher baseline (", B_highB, ")")
  )
) %>%
  mutate(Panel = "Baseline shifts magnitude (via b·Baseline)")

right_seq <- bind_rows(
  data.frame(
    t = tgrid,
    y = curve_fun(tgrid, A = A_high, k = k_demo),
    line = "Higher asymptote"
  ),
  data.frame(
    t = tgrid,
    y = curve_fun(tgrid, A = A_low, k = k_demo),
    line = "Lower asymptote"
  )
) %>%
  mutate(Panel = "Sequence shifts magnitude (via A)")

right_all <- bind_rows(right_base, right_seq) %>%
  mutate(
    lt = case_when(
      Panel == "Baseline shifts magnitude (via b·Baseline)" & grepl("^Higher baseline", line) ~ "dashed",
      Panel == "Baseline shifts magnitude (via b·Baseline)" ~ "solid",
      Panel == "Sequence shifts magnitude (via A)" & line == "Higher asymptote" ~ "dashed",
      TRUE ~ "solid"
    ),
    
    # ---- CHANGE HERE: different tones within each facet ----
    col = case_when(
      !USE_COLOR ~ col_black,
      
      Panel == "Baseline shifts magnitude (via b·Baseline)" &
        grepl("^Higher baseline", line) ~ col_base_dark,
      Panel == "Baseline shifts magnitude (via b·Baseline)" ~ col_base_light,
      
      Panel == "Sequence shifts magnitude (via A)" &
        line == "Higher asymptote" ~ col_seq_dark,
      TRUE ~ col_seq_light
    )
  )

p_left <- ggplot(left_df, aes(t, y)) +
  geom_line(linewidth = 1.8, color = "black") +
  geom_hline(yintercept = A_high, linetype = "dashed", linewidth = 0.95) +
  geom_vline(xintercept = t_half, linetype = "dotted", linewidth = 0.95) +
  geom_point(aes(x = t_half, y = y_half), size = 2.8) +
  
  annotate("text",
           x = 6.8, y = A_high + 0.45,
           label = "Asymptote (A)",
           fontface = "bold", size = 5, hjust = 0) +
  
  annotate("text",
           x = t_half + 0.45, y = 3.40,
           label = paste0("Half-time ≈ ln(2)/k\n≈ ", round(t_half, 1), " mo"),
           fontface = "bold", size = 4.6, hjust = 0, vjust = 1) +
  
  # Δ(0)=0 moved near the origin (0,0)
  annotate("label",
           x = 0.75, y = 0.2,
           label = "Δ(0) = 0",
           fontface = "bold", size = 4.6, hjust = 0,
           label.size = 0, fill = "white",
           label.padding = unit(0.15, "lines")) +
  annotate("label",
           x = 0.75, y = -0.55,
           label = "Asymptotic recovery",
           size = 4.2, hjust = 0,
           label.size = 0, fill = "white",
           label.padding = unit(0.12, "lines")) +
  
  labs(
    title = "Model shape (biologically constrained)",
    x = "Months since initiation",
    y = "Predicted change"
  ) +
  coord_cartesian(ylim = c(-1.25, A_high + 1.4)) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title  = element_text(face = "bold", size = 18, margin = margin(b = 6)),
    axis.title  = element_text(face = "bold", size = 14),  # x and y same size
    plot.margin = margin(t = 6, r = 10, b = 10, l = 10)
  )

x_lab <- 23.1

label_df <- bind_rows(
  data.frame(
    Panel = "Baseline shifts magnitude (via b·Baseline)",
    x = x_lab,
    y = curve_fun_base(x_lab, A_low, b_demo, B_highB, k_demo) + 0.75,
    lab = paste0("Higher baseline (", B_highB, ")")
  ),
  data.frame(
    Panel = "Baseline shifts magnitude (via b·Baseline)",
    x = x_lab,
    y = curve_fun_base(x_lab, A_low, b_demo, B_low, k_demo) - 0.85,
    lab = paste0("Lower baseline (", B_low, ")")
  ),
  data.frame(
    Panel = "Sequence shifts magnitude (via A)",
    x = x_lab,
    y = curve_fun(x_lab, A_high, k_demo) + 0.80,
    lab = "Higher asymptote"
  ),
  data.frame(
    Panel = "Sequence shifts magnitude (via A)",
    x = x_lab,
    y = curve_fun(x_lab, A_low, k_demo) - 1.15,
    lab = "Lower asymptote"
  )
)

p_right <- ggplot(right_all, aes(t, y)) +
  geom_line(aes(linetype = lt, color = col), linewidth = 1.6) +
  scale_linetype_identity() +
  scale_color_identity() +
  facet_wrap(~Panel, ncol = 1, scales = "free_y") +
  
  geom_text(
    data = data.frame(
      Panel = c("Baseline shifts magnitude (via b·Baseline)",
                "Sequence shifts magnitude (via A)"),
      x = 0.6, y = c(Inf, Inf),
      lab = "Same k (rate)"
    ),
    aes(x = x, y = y, label = lab),
    inherit.aes = FALSE,
    vjust = 1.2, hjust = 0,
    fontface = "bold",
    size = 4.4
  ) +
  
  geom_label(
    data = label_df,
    aes(x = x, y = y, label = lab),
    inherit.aes = FALSE,
    fontface = "bold",
    size = 4.2,
    hjust = 1,
    label.size = 0,
    label.padding = unit(0.28, "lines"),
    fill = "white"
  ) +
  
  labs(
    title = "How the model separates effects:",
    x = "Months since initiation",
    y = "Predicted change"
  ) +
  coord_cartesian(xlim = c(0, 24), clip = "off") +
  theme_minimal(base_size = 16) +
  theme(
    plot.title  = element_text(face = "bold", size = 18, margin = margin(b = 6)),
    axis.title  = element_text(face = "bold", size = 14),  # x and y same size
    strip.text  = element_text(face = "bold", size = 15),
    legend.position = "none",
    panel.spacing   = unit(1.0, "lines"),
    plot.margin = margin(t = 6, r = 40, b = 10, l = 10)
  )

# ----------------------------
g_left  <- ggplotGrob(p_left)
g_right <- ggplotGrob(p_right)

out_file <- "ISEF_Fig3_Model_Schematic_CLEAN_v6.png"
png(out_file, width = 2200, height = 1100, res = 200)

grid.newpage()

lay <- grid.layout(
  nrow = 2, ncol = 2,
  heights = unit(c(0.18, 0.82), "npc"),
  widths  = unit(c(1, 1.08), "null")
)
pushViewport(viewport(layout = lay))

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1:2))
grid.text(
  expression(Delta(t) == ((A[seq]) + b %.% Baseline) * (1 - e^{-k*t})),
  gp = gpar(fontface = "bold", fontsize = 22)
)
popViewport()

pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
grid.draw(g_left)
popViewport()

pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2))
grid.draw(g_right)
popViewport()

dev.off()

message("Saved: ", out_file)

# =========================================================
# External validation
# =========================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(minpack.lm)


# ----------------------------
# 1) Manually inputing external validation data
# ----------------------------
ext_patients <- tibble::tribble(
  ~Patient,    ~Baseline, ~Time, ~ObsChange,
  "Patient 1", 25,         0,     0,
  "Patient 1", 25,        11,     7,
  "Patient 1", 25,        24,     7,
  
  "Patient 2", 19,         0,     0,
  "Patient 2", 19,        13,     2,
  "Patient 2", 19,        24,    11,
  
  "Patient 3", 15,         0,     0,
  "Patient 3", 15,        10,     4,
  "Patient 3", 15,        17,     3
)

# ----------------------------
# 2) Building model prediction grid per patient (Z→R only)
# ----------------------------
tgrid <- seq(0, 24, length.out = 220)

pred_grid <- ext_patients %>%
  distinct(Patient, Baseline) %>%
  tidyr::crossing(Time = tgrid) %>%
  mutate(
    Sequence = "Z>R",
    A_N_R = 0,
    A_Z_R = 1
  )

# Point predictions from your fitted model
pred_grid$fit <- predict(asym_model, newdata = pred_grid)


# 3) Bootstrap 95% uncertainty band for the prediction curves
set.seed(123)
n_boot <- 1000

boot_mat <- matrix(NA_real_, nrow = nrow(pred_grid), ncol = n_boot)

for (i in seq_len(n_boot)) {
  
  idx <- sample(seq_len(nrow(hfmse)), replace = TRUE)
  boot_data <- hfmse[idx, ] %>%
    mutate(
      SE = SD_Change / sqrt(N),
      w  = 1 / (SE^2)
    )
  
  boot_fit <- try(
    nlsLM(
      Change ~ ((A_N * A_N_R + A_Z * A_Z_R) + b * Baseline) * (1 - exp(-k * Time)),
      data    = boot_data,
      weights = w,
      start   = coef(asym_model),
      control = nls.lm.control(maxiter = 500)
    ),
    silent = TRUE
  )
  
  if (!inherits(boot_fit, "try-error")) {
    boot_mat[, i] <- predict(boot_fit, newdata = pred_grid)
  }
}

pred_grid$lower <- apply(boot_mat, 1, quantile, probs = 0.025, na.rm = TRUE)
pred_grid$upper <- apply(boot_mat, 1, quantile, probs = 0.975, na.rm = TRUE)

patient_cols_dark <- c(
  "Patient 1" = "#F8766D",
  "Patient 2" = "#00BA38",
  "Patient 3" = "#619CFF"
)

# 5) Plot
p4 <- ggplot() +
  
  geom_ribbon(
    data = pred_grid,
    aes(x = Time, ymin = lower, ymax = upper, fill = Patient),
    alpha = 0.18,
    color = NA
  ) +
  
  geom_line(
    data = pred_grid,
    aes(x = Time, y = fit, color = Patient),
    linewidth = 1.7,
    linetype = "solid"
  ) +
  
  geom_line(
    data = ext_patients,
    aes(x = Time, y = ObsChange, color = Patient, group = Patient),
    linewidth = 1.5,
    linetype = "dashed",
    alpha = 0.45
  ) +
  
  geom_point(
    data = ext_patients,
    aes(x = Time, y = ObsChange, color = Patient),
    size = 3.4
  ) +
  
  # Clinically meaningful line (+3)
  geom_hline(yintercept = 3, linetype = "dashed", linewidth = 0.9) +
  annotate(
    "text",
    x = 17.5, y = 3.55,
    label = "Clinically meaningful (+3)",
    fontface = "bold",
    size = 4.4,
    hjust = 0
  ) +
  
  labs(
    title = "External Validation (Z→R): Patient Trajectories vs Model",
    subtitle = "Solid (dark) = model prediction (95% band); dashed (light) = external patients",
    x = "Months since treatment initiation",
    y = "Change in HFMSE (ΔHFMSE)",
    color = NULL,
    fill  = NULL
  ) +
  
  scale_color_manual(values = patient_cols_dark) +
  scale_fill_manual(values = patient_cols_dark) +
  
  coord_cartesian(xlim = c(-0.5, 24.5)) +
  
  theme_minimal(base_size = 18) +
  theme(
    plot.title.position = "plot",
    plot.title = element_text(face = "bold", size = 20, margin = margin(b = 3), hjust = 0.5),
    plot.subtitle = element_text(size = 12, margin = margin(b = 8), hjust = 0.5),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.text = element_text(size = 15),
    panel.grid.minor = element_line(linewidth = 0.4),
    plot.margin = margin(t = 10, r = 18, b = 10, l = 12)
  )

print(p4)

# 6) Save (horizontally compressed vs your prior wide version)
ggsave(
  filename = "ISEF_Fig4_HFMSE_ExternalValidation_ZR.png",
  plot = p4,
  width = 8.5,     # narrower to "compress" horizontally
  height = 6.8,
  dpi = 300
)

library(dplyr)
library(ggplot2)

# Helper: Weighted RMSE
w_rmse <- function(y, yhat, w) {
  sqrt(sum(w * (y - yhat)^2, na.rm = TRUE) / sum(w, na.rm = TRUE))
}

# 1) Fit baseline models (same weights w)
lin_fit  <- lm(Change ~ 0 + Time, data = hfmse, weights = w)  # force through origin
quad_fit <- lm(Change ~ 0 + Time + I(Time^2), data = hfmse, weights = w)  # force through origin

# Predictions on observed cohort-time rows
hfmse <- hfmse %>%
  mutate(
    pred_lin  = predict(lin_fit,  newdata = hfmse),
    pred_quad = predict(quad_fit, newdata = hfmse),
    pred_ml   = predict(asym_model, newdata = hfmse)
  )

# 2) Compute weighted RMSE
rmse_lin  <- w_rmse(hfmse$Change, hfmse$pred_lin,  hfmse$w)
rmse_quad <- w_rmse(hfmse$Change, hfmse$pred_quad, hfmse$w)
rmse_ml   <- w_rmse(hfmse$Change, hfmse$pred_ml,   hfmse$w)

rmse_tbl <- data.frame(
  Model = c("Linear baseline", "Quadratic baseline", "Biologically constrained ML"),
  Weighted_RMSE = c(rmse_lin, rmse_quad, rmse_ml)
)

rmse_tbl <- rmse_tbl %>%
  mutate(
    Model = factor(Model, levels = c("Linear baseline", "Quadratic baseline", "Biologically constrained ML")),
    Type  = ifelse(Model == "Biologically constrained ML", "ML", "Baseline")
  )

print(rmse_tbl)

# Save table
write.csv(rmse_tbl, "HFMSE_ModelComparison_RMSE_Table.csv", row.names = FALSE)

baseline_mean <- mean(c(rmse_lin, rmse_quad))
ratio <- baseline_mean / rmse_ml   # ≈ 2.7

rmse_tbl <- rmse_tbl %>%
  mutate(
    label_val = sprintf("%.3f", Weighted_RMSE)
  )
col_ml <- "turquoise4"
col_baseline <- "gray80"

ymax <- max(rmse_tbl$Weighted_RMSE) * 1.22

# 4) Plot
p_rmse <- ggplot(rmse_tbl, aes(x = Model, y = Weighted_RMSE, fill = Type)) +
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
    y = rmse_ml + 0.35,
    label = paste0("≈ ", round(ratio, 1), "× lower error"),
    fontface = "bold",
    size = 5.6,
    color = col_ml
  ) +
  
  labs(
    title = "HFMSE: Model vs Baseline Comparisons",
    subtitle = "Validation that the biologically constrained ML model reduces prediction error vs standard baselines",
    x = NULL,
    y = "Weighted RMSE (lower is better)",
    fill = NULL,
    caption = "Weights: SE = SD/√N; w = 1/SE². All models constrained to Δ(0)=0."
  ) +
  
  coord_cartesian(ylim = c(0, ymax), clip = "off") +
  
  theme_minimal(base_size = 18) +
  theme(
    plot.title = element_text(face = "bold", size = 26, margin = margin(b = 6)),
    plot.subtitle = element_text(size = 15, margin = margin(b = 12)),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 18),
    axis.title = element_text(face = "bold", size = 18),
    legend.position = "none",
    panel.grid.minor = element_blank(),
    plot.margin = margin(t = 10, r = 30, b = 10, l = 10),
    plot.caption = element_text(size = 13, hjust = 0)
  )

print(p_rmse)

ggsave(
  filename = "ISEF_Fig5_HFMSE_ModelComparison_RMSE.png",
  plot = p_rmse,
  width = 12,
  height = 6,
  dpi = 300
)
