library(tidyverse)
library(terra)
library(ozmaps)
library(paletteer)
library(patchwork)

# Shortcut plotting functions
prettyplot <- function() {
  theme_classic() +
    theme(
      text = element_text(family = "serif", size = 12, colour = "black"),
      legend.position = "none",
      axis.title.y = element_text(vjust = 1.5),
      axis.title.x = element_text(vjust = 1.5),
      legend.title = element_blank()
    )
}

rm_y_axis <- function() {
  theme(
    axis.title.y = element_blank(), 
    axis.text.y = element_blank(), 
    axis.ticks.y = element_blank(),
    axis.line = element_blank()
  )
}

rm_x_axis <- function() {
  theme(
    axis.title.x = element_blank(), 
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(),
    axis.line = element_blank()
  )
}

plot_sim_results <- function(df, CS = NA) {  
  # --- helper: base theme ---
  base_theme <- function(d) {
    ggplot(d, aes(days, value)) +
      theme_classic(base_size = 10) +
      theme(
        legend.position  = "bottom", 
        legend.title     = element_blank(),
        strip.background = element_blank(),
        strip.text       = element_blank()
      ) +
      labs(x = "Day")
  }
  
  # --- single-variable panels ---
  p_weight <- base_theme(filter(df, output_var == "weight")) +
    geom_line() + labs(title = "Individiaul weight (g)", y = NULL) 
  
  if (!is.na(CS)) {
    p_weight <- p_weight +
      geom_hline(yintercept = CS, linewidth = 0.75, colour = "red", linetype = "dashed") +
      annotate("text", x = min(df$days), y = CS, label = sprintf("CS: %d", CS), hjust = -0.1, vjust = -0.5, colour = "red")
  }
  
  p_dw <- base_theme(filter(df, output_var == "dw")) +
    geom_line() + labs(title = "Daily weight change (g/d)", y = NULL)
  
  p_sgr <- base_theme(filter(df, output_var == "SGR")) +
    geom_line() + labs(title = "SGR (%)", y = NULL)
  
  # --- energy: overlaid lines ---
  p_energy <- base_theme(
    df %>% 
      filter(output_var %in% c("E_assim", "E_somat", "metab"))
    ) +
    geom_line(aes(colour = output_var)) +
    labs(title = "Energy (J/g/d)", y = NULL)
  
  # --- environment: facetted, free y ---
  p_env <- base_theme(
      df %>% 
        filter(output_var %in% c("water_temp", "T_response", "rel_feeding"))
    ) +
    geom_line(aes(colour = output_var)) +
    facet_wrap(~output_var, scales = "free_y", ncol = 1) +
    labs(title = "Temperature & feeding", y = NULL)
  
  # --- food / ingestion: overlaid ---
  p_food <- base_theme(
      df %>% 
        filter(output_var %in% c("food_prov", "food_enc", "ing_pot", "ing_act", "weight")) %>% 
        pivot_wider(names_from = output_var, values_from = value) %>% 
        mutate(
          food_prov = food_prov/weight, 
          food_enc = food_enc/weight, 
          ing_pot = ing_pot/weight, 
          ing_act = ing_act/weight
        ) %>% 
        select(-weight) %>% 
        pivot_longer(cols = -days, names_to = "output_var", values_to = "value")
    ) +
    geom_line(aes(colour = output_var)) +
    labs(title = "Feeding & ingestion (g/g/d)", y = NULL)
  
  # --- O2 & NH4: overlaid ---
  p_o2_nh4 <- base_theme(filter(df, output_var %in% c("O2", "NH4"))) +
    geom_line(aes(colour = output_var)) +
    labs(title = "O2 & NH4 (g/g/d)", y = NULL)
  
  # --- excretion / uneaten: totals black, components coloured ---
  d_excr <- df %>% 
    filter(!str_detect(output_var, "carbon|nitrogen")) %>%
    filter(str_detect(output_var, "excr|uneat")) %>%
    mutate(
      type = str_split_i(output_var, "_", 1) %>% as.factor(),
      out = str_split_i(output_var, "_", 2) %>% as.factor()
    )
  
  p_excr <- base_theme(d_excr) +
    geom_line(aes(colour = type, linetype = out)) +
    scale_colour_manual(values = c("total" = "black", setNames(scales::hue_pal()(3), c("P", "L", "C")))) +
    labs(title = "Excretion & uneaten feed (g/g/d)", y = NULL)
  
  # --- elemental budgets: overlaid ---
  df_elem <- df %>% 
    filter(str_detect(output_var, "carbon|nitrogen")) %>%
    mutate(
      type = case_when(str_detect(output_var, "carbon") ~ "carbon", T ~ "nitrogen") %>% as.factor(),
      out = case_when(str_detect(output_var, "total") ~ "total", str_detect(output_var, "excr") ~ "excr", T ~ "uneat") %>% as.factor()
    )

  p_elem <- base_theme(df_elem) +
    geom_line(aes(colour = type, linetype = out)) +
    scale_linetype_manual(values = c("total" = "solid", "excr" = "dashed", "uneat" = "dotted")) +
    labs(title = "Carbon & nitrogen (g/g/d)", y = NULL)
  
  # --- assemble with patchwork ---
  (p_weight | p_dw | p_sgr) /
    (p_energy | p_env | p_food) /
    (p_o2_nh4 | p_excr | p_elem)
}
