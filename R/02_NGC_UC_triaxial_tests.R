###############################################################
#### Combined New Griffith Criterion (NGC) + Ucar Criterion
#### with triaxial test points and Mohr circles
####
#### Purpose:
####   - Generate the combined NGC–UC failure envelope
####   - Add experimental triaxial points in sigma_3 - sigma_1 space
####   - Add Mohr circles in sigma_n - tau space
####
#### Method:
####   - NGC: new rigorous single-integral formulation
####   - UC : Ucar non-linear compressive criterion
####
#### Outputs:
####   - CSV tables in results/tables/
####   - PNG/PDF figures in results/figures/
####
#### Suggested file name:
####   02_plot_NGC_UC_with_triaxial_tests.R
###############################################################

# ============================================================
# 1. Install and load libraries
# ============================================================

paquetes <- c("ggplot2", "dplyr")

instalar_faltantes <- paquetes[
  !sapply(paquetes, requireNamespace, quietly = TRUE)
]

if (length(instalar_faltantes) > 0) {
  install.packages(instalar_faltantes)
}

lapply(paquetes, library, character.only = TRUE)


# ============================================================
# 2. User-defined parameters
# ============================================================

sigma_t_adj <- -9.5    # T0, tensile strength, MPa; negative in tension
sigma_c_adj <- 130     # C0, uniaxial compressive strength, MPa

xi_adj <- sigma_t_adj / sigma_c_adj

n_ngc <- 1000
n_uc  <- 5000

# UC curve generation domain
uc_sigma3_max_factor <- 100
sigma3_uc_max <- uc_sigma3_max_factor * sigma_c_adj

use_ductile_transition_limit <- TRUE
sigma1_sigma3_min_ratio <- 3.4


# ============================================================
# 2.1. Triaxial test data
# ============================================================

s_3_test <- c(-7.9, -7.5, -7.8, -7.8, -10.6, -9.7, -9.8, -6.9, -4.0, -3.0, 0, 2.5, 4.3)

s_1_test <- c(7.5, 15, 30, 60, 70, 80, 80,90, 100, 120, 130, 140, 150)



df_tri <- data.frame(
  id = seq_along(s_3_test),
  sigma_3 = s_3_test,
  sigma_1 = s_1_test
)


# ============================================================
# 2.2. Plot limits based on experimental data
# ============================================================

plot_margin <- 0.10

use_default_plot_limits <- TRUE

# sigma_3 - sigma_1 plot limits
default_xlim_sigma3_sigma1 <- c(
  min(c(sigma_t_adj, df_tri$sigma_3), na.rm = TRUE) * (1 + plot_margin),
  max(c(df_tri$sigma_3, 0.1 * sigma_c_adj), na.rm = TRUE) * (1 + plot_margin)
)

default_ylim_sigma3_sigma1 <- c(
  0,
  max(c(df_tri$sigma_1, sigma_c_adj), na.rm = TRUE) * (1 + plot_margin)
)

# sigma_n - tau plot limits, including Mohr circles
mohr_centres <- 0.5 * (df_tri$sigma_1 + df_tri$sigma_3)
mohr_radii   <- 0.5 * (df_tri$sigma_1 - df_tri$sigma_3)

default_xlim_sigma_tau <- c(
  min(c(sigma_t_adj, df_tri$sigma_3), na.rm = TRUE) * (1 + plot_margin),
  max(c(df_tri$sigma_1, sigma_c_adj), na.rm = TRUE) * (1 + plot_margin)
)

default_ylim_sigma_tau <- c(
  0,
  max(mohr_radii, na.rm = TRUE) * (1 + plot_margin)
)

# Manual limits used only if use_default_plot_limits <- FALSE
manual_xlim_sigma3_sigma1 <- default_xlim_sigma3_sigma1
manual_ylim_sigma3_sigma1 <- default_ylim_sigma3_sigma1

manual_xlim_sigma_tau <- default_xlim_sigma_tau
manual_ylim_sigma_tau <- default_ylim_sigma_tau


# ============================================================
# 2.3. Output folders
# ============================================================

dir.create("results", showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)


# ============================================================
# 3. Basic checks
# ============================================================

stopifnot(sigma_t_adj < 0)
stopifnot(sigma_c_adj > 0)
stopifnot(abs(sigma_t_adj) < sigma_c_adj)
stopifnot(length(s_3_test) == length(s_1_test))


# ============================================================
# 4. Ucar criterion
# ============================================================

k1 <- (-(1 + abs(xi_adj)) +
         sqrt((1 + 7 * abs(xi_adj)) * (1 - abs(xi_adj)))) /
  (2 * abs(xi_adj))

k2 <- (1 - k1 * (-xi_adj)) / sqrt(pmax(1e-12, -xi_adj))

sigma_3 <- seq(0, sigma3_uc_max, length.out = n_uc)
sigma3_bar <- sigma_3 / sigma_c_adj

sigma1_bar <- k1 * (sigma3_bar - xi_adj) +
  k2 * sqrt(pmax(0, sigma3_bar - xi_adj))

tan_alpha2_ucar <- k1 +
  k2 / (2 * sqrt(pmax(1e-12, sigma3_bar - xi_adj)))

alpha_rad <- atan(sqrt(tan_alpha2_ucar))

sigma_n_bar <- sigma1_bar * cos(alpha_rad)^2 +
  sigma3_bar * sin(alpha_rad)^2

tau_bar <- 0.5 * (sigma1_bar - sigma3_bar) * sin(2 * alpha_rad)

curve_uc <- data.frame(
  branch    = "UC",
  sigma_3   = sigma_3,
  sigma_1   = sigma1_bar * sigma_c_adj,
  sigma_n   = sigma_n_bar * sigma_c_adj,
  tau       = tau_bar * sigma_c_adj,
  alpha_deg = alpha_rad * 180 / pi,
  beta_deg  = 2 * alpha_rad * 180 / pi - 90
) %>%
  mutate(
    sigma1_sigma3_ratio = ifelse(
      sigma_3 > 0,
      sigma_1 / sigma_3,
      Inf
    )
  )

if (use_ductile_transition_limit) {
  
  first_ductile_index <- which(
    curve_uc$sigma1_sigma3_ratio <= sigma1_sigma3_min_ratio
  )[1]
  
  if (!is.na(first_ductile_index)) {
    curve_uc <- curve_uc[1:first_ductile_index, ]
  } else {
    warning(
      "The UC curve did not reach the selected sigma_1/sigma_3 ratio. ",
      "Consider increasing uc_sigma3_max_factor."
    )
  }
}


# ============================================================
# 5. NGC parameter calculation
# ============================================================

calcular_beta1 <- function(a, ratio) {
  
  arg <- 0.75 * a * (-ratio)^(-0.25)
  
  if (!is.finite(arg) || arg < 0) {
    return(NA_real_)
  }
  
  beta1 <- 2 * (atan(sqrt(arg)) - pi / 4)
  
  return(beta1)
}


ecuacion_objetivo <- function(a, ratio, xi_adj) {
  
  if (!is.finite(a) || a <= 0) {
    return(1e12)
  }
  
  lambda <- tryCatch(
    (1 - a * (-ratio)^0.75) / (-ratio),
    error = function(e) return(1e12)
  )
  
  if (!is.finite(lambda)) {
    return(1e12)
  }
  
  beta_1 <- calcular_beta1(a, ratio)
  
  if (is.na(beta_1)) {
    return(1e12)
  }
  
  arg2 <- 1 -
    (
      xi_adj /
        (
          (xi_adj / (1 - sin(beta_1))) -
            (0.5 - 1 / exp(1.008))
        )
    )
  
  arg2 <- pmin(1, pmax(-1, arg2))
  beta_2 <- asin(arg2)
  
  tan_term <- tan(pi / 4 + 0.5 * beta_2)
  
  sigma_3_aux <- tryCatch(
    (tan_term^2 / (0.75 * a))^(-4) + ratio,
    error = function(e) return(1e12)
  )
  
  if (!is.finite(sigma_3_aux)) {
    return(1e12)
  }
  
  term1 <- 3 * a * sigma_3_aux
  term2 <- 4 * a * (sigma_3_aux - ratio)
  term3 <- 4 * (sigma_3_aux - ratio)^0.25 * lambda * abs(ratio)
  
  residuo <- term1 + term2 + term3
  
  if (!is.finite(residuo)) {
    return(1e12)
  }
  
  penalizacion <- 0
  
  if (sigma_3_aux > 0) {
    penalizacion <- penalizacion + 1e6 * sigma_3_aux
  }
  
  if (sigma_3_aux < xi_adj) {
    penalizacion <- penalizacion + 1e6 * abs(sigma_3_aux - xi_adj)
  }
  
  return(residuo^2 + penalizacion)
}


calcular_parametros_ngc <- function(xi_adj, sigma_c_adj, sigma_t_adj) {
  
  ratio <- sigma_t_adj / sigma_c_adj
  
  opt <- optim(
    par = 5,
    fn = ecuacion_objetivo,
    method = "L-BFGS-B",
    lower = 2,
    upper = 10,
    ratio = ratio,
    xi_adj = xi_adj
  )
  
  if (opt$convergence != 0) {
    warning("optim did not converge")
  }
  
  a <- opt$par
  lambda <- (1 - a * (-ratio)^0.75) / (-ratio)
  
  beta_1 <- calcular_beta1(a, ratio)
  
  arg2 <- 1 -
    (
      xi_adj /
        (
          (xi_adj / (1 - sin(beta_1))) -
            (0.5 - 1 / exp(1.008))
        )
    )
  
  arg2 <- pmin(1, pmax(-1, arg2))
  beta_2 <- asin(arg2)
  
  tan_term <- tan(pi / 4 + 0.5 * beta_2)
  
  sigma_3_ngc <- (tan_term^2 / (0.75 * a))^(-4) + ratio
  
  list(
    a = a,
    lambda = lambda,
    beta_1_rad = beta_1,
    beta_2_rad = beta_2,
    beta_1_deg = beta_1 * 180 / pi,
    beta_2_deg = beta_2 * 180 / pi,
    sigma_3_ngc_bar = sigma_3_ngc,
    sigma_3_ngc = sigma_3_ngc * sigma_c_adj,
    residual = opt$value,
    convergence = opt$convergence
  )
}


ngc_params <- calcular_parametros_ngc(
  xi_adj = xi_adj,
  sigma_c_adj = sigma_c_adj,
  sigma_t_adj = sigma_t_adj
)


# ============================================================
# 6. Rigorous NGC curve: new formulation
# ============================================================

a <- ngc_params$a
beta_1_NGC <- ngc_params$beta_1_rad

Omega <- (3 * a / 4)^4

Omega_1 <- (0.5 - xi_adj / (1 - sin(beta_1_NGC))) -
  (7 / 6) * Omega *
  tan(pi / 4 - beta_1_NGC / 2)^6

beta_vals <- seq(pi / 2, beta_1_NGC, length.out = n_ngc)

curve_ngc <- data.frame(beta_rad = beta_vals) %>%
  mutate(
    branch = "NGC",
    beta_deg = beta_rad * 180 / pi,
    alpha_deg = 45 + 0.5 * beta_deg,
    
    tan_term = tan(pi / 4 - beta_rad / 2),
    
    sigma_n_pri =
      (
        xi_adj +
          (1 - sin(beta_rad)) *
          (
            Omega_1 +
              (7 * Omega / 6) * tan_term^6
          )
      ) / abs(xi_adj),
    
    tau_pri =
      (
        cos(beta_rad) *
          (
            Omega_1 +
              Omega * tan_term^6 *
              (
                7 / 6 - 1 / (1 + sin(beta_rad))
              )
          )
      ) / abs(xi_adj),
    
    sigma_n = sigma_n_pri * abs(sigma_t_adj),
    tau     = tau_pri * abs(sigma_t_adj)
  )

eps_cos <- 1e-10

curve_ngc <- curve_ngc %>%
  mutate(
    delta_sigma = ifelse(
      abs(cos(beta_rad)) < eps_cos,
      abs(sigma_t_adj),
      2 * tau / cos(beta_rad)
    ),
    sum_sigma = 2 * sigma_n + delta_sigma * sin(beta_rad),
    sigma_1 = 0.5 * (sum_sigma + delta_sigma),
    sigma_3 = 0.5 * (sum_sigma - delta_sigma)
  )

curve_ngc$sigma_3[1] <- sigma_t_adj
curve_ngc$sigma_1[1] <- 0
curve_ngc$sigma_n[1] <- sigma_t_adj
curve_ngc$tau[1]     <- 0

curve_ngc <- curve_ngc %>%
  select(
    branch,
    sigma_3,
    sigma_1,
    sigma_n,
    tau,
    alpha_deg,
    beta_deg,
    beta_rad,
    sigma_n_pri,
    tau_pri
  )


# ============================================================
# 7. Combined curve and export-ready filtered tables
# ============================================================

curve_combined <- bind_rows(
  curve_ngc %>%
    select(branch, sigma_3, sigma_1, sigma_n, tau, alpha_deg, beta_deg),
  curve_uc %>%
    select(branch, sigma_3, sigma_1, sigma_n, tau, alpha_deg, beta_deg)
)

if (use_default_plot_limits) {
  
  xlim_sigma3_sigma1 <- default_xlim_sigma3_sigma1
  ylim_sigma3_sigma1 <- default_ylim_sigma3_sigma1
  
  xlim_sigma_tau <- default_xlim_sigma_tau
  ylim_sigma_tau <- default_ylim_sigma_tau
  
} else {
  
  xlim_sigma3_sigma1 <- manual_xlim_sigma3_sigma1
  ylim_sigma3_sigma1 <- manual_ylim_sigma3_sigma1
  
  xlim_sigma_tau <- manual_xlim_sigma_tau
  ylim_sigma_tau <- manual_ylim_sigma_tau
}

curve_combined_visible <- curve_combined %>%
  filter(
    sigma_3 >= xlim_sigma3_sigma1[1],
    sigma_3 <= xlim_sigma3_sigma1[2],
    sigma_1 >= ylim_sigma3_sigma1[1],
    sigma_1 <= ylim_sigma3_sigma1[2]
  )


# ============================================================
# 8. Mohr circles from triaxial tests
# ============================================================

make_mohr_circles <- function(df, n = 361) {
  
  theta <- seq(0, 2 * pi, length.out = n)
  
  do.call(
    rbind,
    lapply(seq_len(nrow(df)), function(i) {
      
      s1  <- df$sigma_1[i]
      s3  <- df$sigma_3[i]
      cen <- 0.5 * (s1 + s3)
      rad <- 0.5 * (s1 - s3)
      
      data.frame(
        id = df$id[i],
        sigma = cen + rad * cos(theta),
        tau = rad * sin(theta)
      )
    })
  )
}

df_mohr <- make_mohr_circles(df_tri, n = 361)

df_mohr_labels <- df_tri %>%
  mutate(
    sigma = sigma_1,
    tau = 0
  )


# ============================================================
# 9. Export CSV files
# ============================================================

write.csv(
  df_tri,
  "results/tables/triaxial_test_points.csv",
  row.names = FALSE
)

write.csv(
  df_mohr,
  "results/tables/mohr_circles_from_triaxial_tests.csv",
  row.names = FALSE
)

write.csv(
  curve_ngc,
  "results/tables/curve_NGC_rigorous_new.csv",
  row.names = FALSE
)

write.csv(
  curve_uc,
  "results/tables/curve_UC.csv",
  row.names = FALSE
)

write.csv(
  curve_combined,
  "results/tables/curve_combined_NGC_UC.csv",
  row.names = FALSE
)

write.csv(
  curve_combined_visible,
  "results/tables/curve_combined_NGC_UC_visible_range.csv",
  row.names = FALSE
)

parameters_table <- data.frame(
  sigma_t_adj = sigma_t_adj,
  sigma_c_adj = sigma_c_adj,
  xi_adj = xi_adj,
  k1 = k1,
  k2 = k2,
  a = ngc_params$a,
  lambda = ngc_params$lambda,
  beta_1_deg = ngc_params$beta_1_deg,
  beta_2_deg = ngc_params$beta_2_deg,
  sigma_3_ngc_bar = ngc_params$sigma_3_ngc_bar,
  sigma_3_ngc = ngc_params$sigma_3_ngc,
  residual = ngc_params$residual,
  convergence = ngc_params$convergence,
  Omega = Omega,
  Omega_1 = Omega_1,
  uc_sigma3_max_factor = uc_sigma3_max_factor,
  sigma3_uc_max = sigma3_uc_max,
  use_ductile_transition_limit = use_ductile_transition_limit,
  sigma1_sigma3_min_ratio = sigma1_sigma3_min_ratio,
  uc_final_sigma_3 = max(curve_uc$sigma_3, na.rm = TRUE),
  uc_final_sigma_1 = max(curve_uc$sigma_1, na.rm = TRUE),
  uc_final_sigma1_sigma3_ratio = tail(curve_uc$sigma1_sigma3_ratio, 1)
)

write.csv(
  parameters_table,
  "results/tables/parameters_NGC_UC_with_triaxial_tests.csv",
  row.names = FALSE
)


# ============================================================
# 10. Plot sigma_n - tau with Mohr circles
# ============================================================

plot_sigma_tau_mohr <- ggplot() +
  geom_line(
    data = curve_ngc,
    aes(x = sigma_n, y = tau, color = branch),
    linewidth = 1.2
  ) +
  geom_line(
    data = curve_uc,
    aes(x = sigma_n, y = tau, color = branch),
    linewidth = 1.2,
    alpha = 0.8
  ) +
  geom_path(
    data = df_mohr,
    aes(x = sigma, y = tau, group = id),
    linewidth = 0.5,
    alpha = 0.75
  ) +
  geom_text(
    data = df_mohr_labels,
    aes(x = sigma, y = tau, label = id),
    inherit.aes = FALSE,
    nudge_y = 0.03 * diff(ylim_sigma_tau),
    size = 3,
    check_overlap = TRUE
  ) +
  scale_color_manual(
    values = c("NGC" = "blue", "UC" = "orange")
  ) +
  coord_fixed(
    ratio = 1,
    xlim = xlim_sigma_tau,
    ylim = ylim_sigma_tau,
    expand = FALSE
  ) +
  labs(
    title = "Combined NGC + UC failure envelope with Mohr circles",
    x = expression(sigma[n]~"(MPa)"),
    y = expression(tau~"(MPa)"),
    color = "Criterion"
  ) +
  theme_gray()

print(plot_sigma_tau_mohr)

ggsave(
  "results/figures/combined_NGC_UC_with_Mohr_circles_sigma_tau.png",
  plot = plot_sigma_tau_mohr,
  width = 8,
  height = 6,
  dpi = 300
)

ggsave(
  "results/figures/combined_NGC_UC_with_Mohr_circles_sigma_tau.pdf",
  plot = plot_sigma_tau_mohr,
  width = 8,
  height = 6
)


# ============================================================
# 11. Plot sigma_3 - sigma_1 with triaxial points
# ============================================================

plot_sigma3_sigma1_tests <- ggplot() +
  geom_line(
    data = curve_ngc,
    aes(x = sigma_3, y = sigma_1, color = branch),
    linewidth = 1.2,
    alpha = 0.85
  ) +
  geom_line(
    data = curve_uc,
    aes(x = sigma_3, y = sigma_1, color = branch),
    linewidth = 1.2,
    alpha = 0.85
  ) +
  geom_point(
    data = df_tri,
    aes(x = sigma_3, y = sigma_1),
    inherit.aes = FALSE,
    size = 2.5,
    color = "red"
  ) +
  scale_color_manual(
    values = c("NGC" = "blue", "UC" = "orange")
  ) +
  coord_fixed(
    ratio = 1,
    xlim = xlim_sigma3_sigma1,
    ylim = ylim_sigma3_sigma1,
    expand = FALSE
  ) +
  labs(
    title = "Combined NGC + UC failure curve with triaxial test points",
    x = expression(sigma[3]~"(MPa)"),
    y = expression(sigma[1]~"(MPa)"),
    color = "Criterion"
  ) +
  theme_gray()

print(plot_sigma3_sigma1_tests)

ggsave(
  "results/figures/combined_NGC_UC_with_triaxial_points_sigma3_sigma1.png",
  plot = plot_sigma3_sigma1_tests,
  width = 8,
  height = 6,
  dpi = 300
)

ggsave(
  "results/figures/combined_NGC_UC_with_triaxial_points_sigma3_sigma1.pdf",
  plot = plot_sigma3_sigma1_tests,
  width = 8,
  height = 6
)


# ============================================================
# 12. Console checks
# ============================================================

cat("\n========================================\n")
cat("Combined NGC + UC with triaxial tests\n")
cat("========================================\n")
cat("sigma_t_adj =", sigma_t_adj, "MPa\n")
cat("sigma_c_adj =", sigma_c_adj, "MPa\n")
cat("xi_adj      =", xi_adj, "\n\n")

cat("Triaxial tests:\n")
cat("Number of tests =", nrow(df_tri), "\n")
cat("sigma_3 range   =", range(df_tri$sigma_3), "MPa\n")
cat("sigma_1 range   =", range(df_tri$sigma_1), "MPa\n\n")

cat("NGC parameters:\n")
cat("a           =", ngc_params$a, "\n")
cat("lambda      =", ngc_params$lambda, "\n")
cat("beta_1      =", ngc_params$beta_1_deg, "deg\n")
cat("beta_2      =", ngc_params$beta_2_deg, "deg\n")
cat("sigma_3_NGC =", ngc_params$sigma_3_ngc, "MPa\n")
cat("residual    =", ngc_params$residual, "\n\n")

cat("Plot limits:\n")
cat("sigma_n - tau xlim:", xlim_sigma_tau, "\n")
cat("sigma_n - tau ylim:", ylim_sigma_tau, "\n")
cat("sigma_3 - sigma_1 xlim:", xlim_sigma3_sigma1, "\n")
cat("sigma_3 - sigma_1 ylim:", ylim_sigma3_sigma1, "\n\n")

cat("Rows exported:\n")
cat("NGC curve             :", nrow(curve_ngc), "\n")
cat("UC curve              :", nrow(curve_uc), "\n")
cat("Combined curve        :", nrow(curve_combined), "\n")
cat("Visible combined curve:", nrow(curve_combined_visible), "\n")
cat("Mohr circle rows      :", nrow(df_mohr), "\n\n")

cat("Files written to:\n")
cat("results/tables/\n")
cat("results/figures/\n")
cat("========================================\n")