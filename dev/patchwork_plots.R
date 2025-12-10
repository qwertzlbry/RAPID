library(patchwork)





p <- ggplot(sim1_long, aes(x = kappa, y = value, color = metric, group = metric)) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = value - sd_value, ymax = value + sd_value, fill = metric),
              alpha = 0.2, color = NA) +
  scale_color_manual(values = c("Confidence Rate" = "steelblue",
                                "Accuracy" = "firebrick")) +
  scale_fill_manual(values = c("Confidence Rate" = "steelblue",
                               "Accuracy" = "firebrick")) +
  labs(
    # title = "Attribute Disclosure Risk and Accuracy vs Dependency Strength",
    # subtitle = sprintf("Mean ± 1 SD over repeated simulations (n = %d per dep)", n_sim),
    x = "κ (Dependency Strength)",
    y = "RAPID (Proportion at Risk)"
  ) +
  scale_x_continuous(breaks = seq(0, 100, 10)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank()
  )


p1 <- ggplot(sim2_summary, aes(x = tau, y = RAPID_mean, color = kappa_label, group = kappa_label)) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  scale_color_manual(
    values = c("κ = 0" = "gray50",
               "κ = 5" = "steelblue",
               "κ = 10" = "dodgerblue",
               "κ = 20" = "firebrick",
               "κ = 50" = "darkred")
  ) +
  labs(
    #title = "RAPID vs Confidence Threshold for Different Dependency Strengths",
    #subtitle = sprintf("Mean over %d simulations per (κ, τ) combination", n_sim),
    x = "Confidence Threshold (τ)",
    y = "",#"RAPID (Proportion at Risk)",
    color = "Dependency\nStrength",
    #caption = "Synthesizer: CART | Attacker: Random Forest"
  ) +
  scale_x_continuous(breaks = seq(0, 1, 0.1)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom",
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

p1 <- p1 + theme(
  legend.position = "bottom",
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank()
)

combined_plot <- p + p1 +
  plot_annotation(tag_levels = 'a', tag_suffix = ')') +
  plot_layout(widths = c(1, 1))

combined_plot
ggsave("sim1_and_sim2_combined.pdf", combined_plot, width = 10, height = 6)
