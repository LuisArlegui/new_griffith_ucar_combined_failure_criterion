###############################################################
#### Combined New Griffith Criterion (NGC) + Ucar Criterion
#### Official stable implementation
####
#### Outputs:
####   - CSV files with NGC, UC and combined curves
####   - PNG/PDF plots in sigma_n - tau and sigma_3 - sigma_1 spaces
####
#### Version:
####   v12_stable_candidate
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
# 2. User-defined parameters to EDIT
# ============================================================

sigma_t_adj <- -20     # T0, tensile strength, MPa; negative in tension
sigma_c_adj <-160     # C0, uniaxial compressive strength, MPa

xi_adj <- sigma_t_adj / sigma_c_adj

n_ngc <- 3000    # edit to densify the plots or tables
n_uc  <- 5000   # edit to densify the plots or tables

# UC curve generation domain
# A dense UC curve is first generated up to:
# sigma3_uc_max = uc_sigma3_max_factor * sigma_c_adj
# Then, if requested, it is truncated when sigma_1 / sigma_3
# reaches the selected ductile-transition threshold.

uc_sigma3_max_factor <- 100
sigma3_uc_max <- uc_sigma3_max_factor * sigma_c_adj

use_ductile_transition_limit <- TRUE
sigma1_sigma3_min_ratio <- 3.4 # Mogi (1966) threshold to ductile failure (Schwartz, 1964 puts it at 4)

# ============================================================
# 2.1. User-defined plot limits
# ============================================================

use_default_plot_limits <- TRUE

# Manual limits used only if use_default_plot_limits <- FALSE

manual_xlim_sigma_tau <- c(1.1 * sigma_t_adj, 1.1 * sigma_c_adj)
manual_ylim_sigma_tau <- c(0, 2*sigma_c_adj)

manual_xlim_sigma3_sigma1 <- c(2 * sigma_t_adj, 0.6 * sigma_c_adj)
manual_ylim_sigma3_sigma1 <- c(0, 3 * sigma_c_adj)

# ============================================================
# 2.2. Output folders
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

# Truncate UC curve at the selected ductile-transition threshold
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


# Recover sigma_3 and sigma_1 from sigma_n, tau and beta

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

# Force exact tensile endpoint
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
# 7. Combined NGC + UC curve
# ============================================================

curve_combined <- bind_rows(
  curve_ngc %>%
    select(branch, sigma_3, sigma_1, sigma_n, tau, alpha_deg, beta_deg),
  curve_uc %>%
    select(branch, sigma_3, sigma_1, sigma_n, tau, alpha_deg, beta_deg)
)


# ============================================================
# 8. Plot limits
# ============================================================

if (use_default_plot_limits) {
  
  # sigma_n - tau plot
  xlim_sigma_tau <- c(1.1 * sigma_t_adj, sigma_c_adj)
  
  tau_at_sigma_c <- curve_combined %>%
    filter(abs(sigma_n - sigma_c_adj) == min(abs(sigma_n - sigma_c_adj))) %>%
    summarise(tau_ref = max(tau, na.rm = TRUE)) %>%
    pull(tau_ref)
  
  ylim_sigma_tau <- c(0, 1.1 * tau_at_sigma_c)
  
  # sigma_3 - sigma_1 plot
  xlim_sigma3_sigma1 <- c(2 * sigma_t_adj, 0.6*sigma_c_adj)
  ylim_sigma3_sigma1 <- c(0, 3 * sigma_c_adj)
  
} else {
  
  xlim_sigma_tau <- manual_xlim_sigma_tau
  ylim_sigma_tau <- manual_ylim_sigma_tau
  
  xlim_sigma3_sigma1 <- manual_xlim_sigma3_sigma1
  ylim_sigma3_sigma1 <- manual_ylim_sigma3_sigma1
}


# ============================================================
# 9. Export CSV files
# ============================================================

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
  "results/tables/parameters_NGC_UC.csv",
  row.names = FALSE
)


# ============================================================
# 10. Plot in sigma_n - tau space
# ============================================================

plot_sigma_tau <- ggplot() +
  geom_line(
    data = curve_ngc,
    aes(x = sigma_n, y = tau, color = branch),
    linewidth = 1.2
  ) +
  geom_line(
    data = curve_uc,
    aes(x = sigma_n, y = tau, color = branch),
    linewidth = 1.2
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
    title = "Combined NGC + UC failure envelope",
    x = expression(sigma[n]~"(MPa)"),
    y = expression(tau~"(MPa)"),
    color = "Criterion"
  ) +
  theme_gray()

print(plot_sigma_tau)

ggsave(
  filename = "results/figures/combined_NGC_UC_sigma_tau.png",
  plot = plot_sigma_tau,
  width = 7,
  height = 6,
  dpi = 300
)

ggsave(
  filename = "results/figures/combined_NGC_UC_sigma_tau.pdf",
  plot = plot_sigma_tau,
  width = 7,
  height = 6
)


# ============================================================
# 11. Plot in sigma_3 - sigma_1 space
# ============================================================

plot_sigma3_sigma1 <- ggplot() +
  geom_line(
    data = curve_ngc,
    aes(x = sigma_3, y = sigma_1, color = branch),
    linewidth = 1.2
  ) +
  geom_line(
    data = curve_uc,
    aes(x = sigma_3, y = sigma_1, color = branch),
    linewidth = 1.2
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
    title = "Combined NGC + UC failure envelope",
    x = expression(sigma[3]~"(MPa)"),
    y = expression(sigma[1]~"(MPa)"),
    color = "Criterion"
  ) +
  theme_gray()

print(plot_sigma3_sigma1)

ggsave(
  filename = "results/figures/combined_NGC_UC_sigma3_sigma1.png",
  plot = plot_sigma3_sigma1,
  width = 7,
  height = 6,
  dpi = 300
)

ggsave(
  filename = "results/figures/combined_NGC_UC_sigma3_sigma1.pdf",
  plot = plot_sigma3_sigma1,
  width = 7,
  height = 6
)


# ============================================================
# 12. Console checks
# ============================================================

cat("\n==============================\n")
cat("Combined NGC + UC calculation\n")
cat("==============================\n")
cat("sigma_t_adj =", sigma_t_adj, "MPa\n")
cat("sigma_c_adj =", sigma_c_adj, "MPa\n")
cat("xi_adj      =", xi_adj, "\n\n")

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
cat("NGC curve      :", nrow(curve_ngc), "\n")
cat("UC curve       :", nrow(curve_uc), "\n")
cat("Combined curve :", nrow(curve_combined), "\n\n")

cat("Files written to:\n")
cat("data/processed/\n")
cat("figures/\n")
cat("==============================\n")