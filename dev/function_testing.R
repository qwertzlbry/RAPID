library(tidyverse)
library(progressr)
library(synthpop)
devtools::load_all()
# edge cases -------------------------------------------------------------------
ng <- function(gi, bi){
  out <- (gi-bi)/(1-bi)
  round(out,3)
}

# Create sequences
gi_vals <- seq(0, 1, by = 0.1)
bi_vals <- seq(0.1, 1, by = 0.1)

# Expand grid of all combinations
grid <- expand.grid(gi = gi_vals, bi = bi_vals)

grid$ng <- mapply(ng, grid$gi, grid$bi)

# Make matrix (rows = gi, columns = bi)
ng_matrix <- with(grid, tapply(ng, list(gi, bi), mean))

grid$at_risk <- grid$ng >= 0.3


ggplot(grid, aes(x = bi, y = gi, fill = ng)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(ng, 2)), size = 2.5) +
  scale_fill_gradientn(
    colors = c("darkgreen", "green", "yellow", "orange", "red", "darkred"),
    values = scales::rescale(c(-9, -2, 0, 0.3, 0.6, 1)),
    limits = c(-9, 1)
  ) +
  geom_contour(aes(z = ng), breaks = 0.3, color = "black", linewidth = 0.5) +
  labs(
    title = "Normalized Gain: Edge Case Analysis",
    subtitle = expression(
      "normalized_gain = " * frac(g[i] - b[i], 1 - b[i]) *
        " | Black line = " * tau * " = 0.3 | Green = safe, Red = at risk"
    ),
    x = expression("Baseline (" * b[i] * ")"),
    y = expression("Predicted Probability (" * g[i] * ")"),
    fill = "Normalized\nGain"
  ) +
  theme_minimal(base_size = 14) +
  coord_fixed()
