library(tidyverse)
library(Polychrome)
library(igraph)
library(colorspace)

here::here("src", "dirs.R") %>% source()

assigned_farms <- file.path(prepdata_path, "assigned_farms.qs") %>% 
  qs2::qs_read() %>% 
  sf::st_drop_geometry()

# which model_names co-occur within a country
af <- assigned_farms %>% 
  mutate(model_name = case_when(model_name == "japanese_seabass" ~ "asian_seabass",  T ~ model_name))
all_models <- af %>% 
  pull(model_name) %>% unique()
co <- af %>% distinct(country, model_name)
n  <- length(all_models)
adj <- matrix(FALSE, n, n, dimnames = list(all_models, all_models))
for (grp in split(co$model_name, co$country)) {
  i <- match(grp, all_models)
  adj[i, i] <- TRUE
}
diag(adj) <- FALSE

# greedy colouring, most-constrained model_names first (Welsh–Powell-ish)
color_id <- setNames(integer(n), all_models)
for (i in order(rowSums(adj), decreasing = TRUE)) {
  used <- color_id[which(adj[i, ])]
  color_id[i] <- min(setdiff(seq_len(n), used[used > 0]))
}

k <- max(color_id)
base <- createPalette(k, c("#3B82F6", "#EF4444", "#10B981"))
base <- desaturate(darken(base, amount = 0.1), amount = 0.2)
models_pal <- setNames(base[color_id], all_models) # still a fixed named vector


# Golden pompano not being gold is throwing me
col1 <- models_pal["golden_pompano"]
col2 <- models_pal["grouper"]
models_pal["grouper"] <- col1
models_pal["golden_pompano"] <- col2


# edge list: every within-country species pair
edges <- co %>%
  group_by(country) %>% 
  filter(n() >= 2) %>% 
  reframe(as.data.frame(t(combn(sort(unique(model_name)), 2)))) %>%
  select(from = V1, to = V2) %>%
  distinct()

g <- graph_from_data_frame(edges, directed = FALSE, vertices = data.frame(name = all_models))

V(g)$color <- models_pal[V(g)$name]   # the exact colours your plots use
V(g)$size  <- 18

# plot(
#   g,
#   layout            = layout_with_fr,   # force-directed; clusters cliques together
#   vertex.label.color = "black",
#   vertex.frame.color = "grey30",
#   edge.color        = "grey80"
# )

# # Focal node
# gf <- g
# focal <- "golden_pompano"
# inc <- incident(gf, V(gf)[focal])   # edge IDs touching that node
# E(gf)$width <- 1
# E(gf)$width[inc] <- 3
# E(gf)$color <- "grey85"
# E(gf)$color[inc] <- "grey20"

# plot(
#   gf,
#   layout             = layout_with_fr,
#   vertex.label.color = "black",
#   vertex.frame.color = "grey30"
# )
