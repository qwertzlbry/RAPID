simulate_microdata <- function(n = 1000, seed = NULL, dep = 1,
                               age_mean = 45, age_sd = 15,
                               age_min = 18, age_max = 90) {
  if (!is.null(seed)) set.seed(seed)

  signal_weight <- sqrt(dep / (1 + dep))
  noise_weight  <- sqrt(1 / (1 + dep))

  rtruncnorm <- function(n, mean, sd, lo, hi) {
    x <- rnorm(n, mean, sd)
    while (any(out <- x < lo | x > hi)) {
      x[out] <- rnorm(sum(out), mean, sd)
    }
    x
  }

  age <- rtruncnorm(n, mean = age_mean, sd = age_sd, lo = age_min, hi = age_max)
  age_z <- scale(age)[, 1]
  latent_SES <- rnorm(n)

  edu_latent <- signal_weight * (0.8 * latent_SES - 0.4 * age_z) + noise_weight * rnorm(n)
  education_num <- cut(edu_latent, breaks = c(-Inf, -0.3, 0.7, Inf), labels = c(0,1,2)) |> as.integer() - 1
  education <- factor(c("low", "medium", "high")[education_num + 1], levels = c("low", "medium", "high"))

  lin_income <- signal_weight * (0.5 * latent_SES + 0.3 * age_z + 0.25 * education_num) + noise_weight * rnorm(n)
  income <- exp(10 + lin_income)

  lin_health <- signal_weight * (0.6 * latent_SES - 0.5 * age_z + 0.2 * education_num + 0.2 * scale(log(income))[,1]) + noise_weight * rnorm(n)
  health_score <- 100 / (1 + exp(-lin_health))

  health_z <- scale(health_score)[, 1]
  log_income_z <- scale(log(income))[, 1]

  eta_diab <- -1.5 + dep * (0.8 * age_z - 0.3 * log_income_z - 0.2 * education_num)
  eta_hyp <- -1.3 + dep * (1.0 * age_z - 0.2 * log_income_z - 0.1 * education_num)
  eta_healthy <- 0

  exp_eta <- cbind(healthy = exp(eta_healthy), diabetic = exp(eta_diab), hypertensive = exp(eta_hyp))
  prob_mat <- exp_eta / rowSums(exp_eta)

  draw_multinom <- function(p, levels) {
    idx <- apply(p, 1, function(prob) sample(seq_along(prob), size = 1, prob = prob))
    factor(levels[idx], levels = levels)
  }
  disease_status <- draw_multinom(prob_mat, c("healthy", "diabetic", "hypertensive"))

  eta_gender <- signal_weight * (0.3 * latent_SES - 0.2 * age_z + 0.2 * education_num)
  prob_male <- 1 / (1 + exp(-eta_gender))
  gender <- factor(rbinom(n, 1, prob_male), labels = c("female", "male"))

  data.frame(gender = gender, age = age, education = education,
             income = income, health_score = health_score, disease_status = disease_status)
}
