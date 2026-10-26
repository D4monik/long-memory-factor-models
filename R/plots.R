# Figures for the climate application.

suppressPackageStartupMessages({
  library(ggplot2)
  library(sf)
})
stopifnot(requireNamespace("rnaturalearthdata", quietly = TRUE))  # map data used by rnaturalearth

DEG_C <- "°C"
BLUE <- "#0070B2"; RED <- "#D55E00"; GRAY <- "#525252"
YEAR_BREAKS <- as.Date(
  paste0(seq(1960, 2020, by = 20), "-01-01")
)

save_pdf <- function(plot, file, width, height) {
  ggsave(file.path(FIG_DIR, file), plot, width = width, height = height,
         units = "in", device = grDevices::cairo_pdf)
}

smooth_ma <- function(x, window = 24) {
  as.numeric(stats::filter(x, rep(1 / window, window), sides = 2))
}

theme_paper <- function() {
  theme_minimal(base_size = 24) + theme(panel.grid.minor = element_blank())
}

# Figure 1: station series (thin) and their cross-sectional mean (bold)
plot_panel <- function(Y, dates) {
  df <- data.frame(date = rep(dates, ncol(Y)), station = rep(colnames(Y), each = nrow(Y)),
                   value = as.vector(Y))
  ggplot(df, aes(date, value, group = station)) +
    geom_line(alpha = 0.05, linewidth = 0.15, colour = "grey35") +
    geom_line(data = data.frame(date = dates, value = rowMeans(Y)), aes(date, value),
              inherit.aes = FALSE, linewidth = 1.1, colour = "black") +
    scale_x_date(breaks = YEAR_BREAKS, date_labels = "%Y") +
    labs(x = NULL, y = paste0("Temperature (", DEG_C, ")")) +
    theme_paper()
}

# Figures 2(a) and 8(a): selected number of factors and its standard deviation
plot_abc <- function(path, s_max) {
  ratio <- if (max(path$sd) > 0) s_max / max(path$sd) else 1
  ggplot(path, aes(x = c)) +
    geom_line(aes(y = s_hat), colour = "black", linewidth = 1) +
    geom_line(aes(y = sd * ratio), colour = BLUE, linewidth = 1) +
    scale_y_continuous(name = expression(hat(s)(c)), breaks = seq(0, s_max, by = 2), limits = c(0, s_max),
                       sec.axis = sec_axis(~ . / ratio, name = "Standard deviation")) +
    scale_x_continuous(name = "Penalty multiplier c") +
    theme_paper()
}

# Figure 2(b): leading eigenvalues; the first one is drawn above an axis break if it
# exceeds twice the second.
plot_eigenvalues <- function(eigenvalues, n_show = 10) {
  df <- data.frame(index = seq_len(n_show), eigenvalue = eigenvalues[seq_len(n_show)])
  df$y <- df$eigenvalue
  p <- ggplot(df, aes(index, y))
  second_max <- max(df$eigenvalue[-1])
  if (df$eigenvalue[1] > 2 * second_max) {
    break_at <- 1.15 * second_max
    df$y[1] <- 1.45 * second_max
    breaks <- pretty(c(0, second_max), n = 5)
    p <- ggplot(df, aes(index, y)) +
      geom_hline(yintercept = break_at, linetype = "dashed", linewidth = 0.4) +
      annotate("text", x = n_show, y = break_at, label = "axis break",
               hjust = 1, vjust = -0.4, size = 5.5) +
      scale_y_continuous(breaks = breaks[breaks <= 1.10 * second_max])
  }
  p + geom_line(linewidth = 0.6) + geom_point(size = 2) +
    scale_x_continuous(breaks = df$index) +
    labs(x = "Eigenvalue index", y = "Eigenvalue") +
    theme_paper()
}

# Figure 3: factors relative to their baseline means, in standard deviations
plot_factors <- function(Z, dates) {
  Z <- as.matrix(Z)
  
  df <- data.frame(date = dates, value = as.numeric(Z[, 1]))
  df$smoothed <- smooth_ma(df$value)
  
  ggplot(df, aes(date)) +
    geom_line(aes(y = value), colour = "grey60", linewidth = 0.35) +
    geom_line(aes(y = smoothed), linewidth = 1.1, na.rm = TRUE) +
    scale_x_date(breaks = YEAR_BREAKS, date_labels = "%Y") +
    labs(x = NULL, y = "Standard deviations") +
    theme_paper()
}

# Figure 4(a): sample-average common component relative to its baseline mean
plot_common_signal <- function(x, dates) {
  df <- data.frame(date = dates, value = x, smoothed = smooth_ma(x))
  ggplot(df, aes(date)) +
    geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey50") +
    geom_line(aes(y = value), colour = "grey60", linewidth = 0.35) +
    geom_line(aes(y = smoothed), linewidth = 1.15, na.rm = TRUE) +
    scale_x_date(breaks = YEAR_BREAKS, date_labels = "%Y") +
    labs(x = NULL, y = paste0("Common component (", DEG_C, ")")) +
    theme_paper()
}

# Figure 4(b): contribution of each factor to the sample-average common component
plot_factor_contributions <- function(M, dates) {
  df <- data.frame(date = rep(dates, ncol(M)), factor = rep(colnames(M), each = nrow(M)),
                   value = as.vector(M))
  df$smoothed <- ave(df$value, df$factor, FUN = smooth_ma)
  ggplot(df, aes(date)) +
    geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey50") +
    geom_line(aes(y = value), colour = "grey60", linewidth = 0.35) +
    geom_line(aes(y = smoothed), linewidth = 1.1, na.rm = TRUE) +
    facet_wrap(~ factor, ncol = 2, scales = "free_y") +
    scale_x_date(breaks = YEAR_BREAKS, date_labels = "%Y") +
    labs(x = NULL, y = paste0("Factor contribution (", DEG_C, ")")) +
    theme_paper() +
    theme(strip.text = element_text(face = "bold"))
}

# Figure 8(b): levels vs. cumulated differences
plot_levels_vs_differences <- function(levels, differences, dates) {
  lab <- c("Levels, two factors", "Differences, one factor")
  df <- data.frame(date = rep(dates, 2), value = c(levels, differences),
                   model = factor(rep(lab, each = length(dates)), levels = lab))
  df$ma <- ave(df$value, df$model, FUN = smooth_ma)
  ggplot(df, aes(date, colour = model, linetype = model)) +
    geom_hline(yintercept = 0, colour = GRAY, linewidth = 0.3) +
    geom_line(aes(y = value), linewidth = 0.3, alpha = 0.35) +
    geom_line(aes(y = ma), linewidth = 1.0, na.rm = TRUE) +
    scale_colour_manual(values = c(BLUE, RED), name = NULL) +
    scale_linetype_manual(values = c("solid", "dashed"), name = NULL) +
    scale_x_date(breaks = YEAR_BREAKS, date_labels = "%Y") +
    labs(x = NULL, y = paste0(DEG_C, " relative to ", BASELINE[1], "–", BASELINE[2])) +
    theme_minimal(base_size = 18) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank())
}

# Figure 6(b): regional means of station warming
plot_regions <- function(reg) {
  reg$region <- reorder(reg$region, reg$warming)
  ggplot(reg, aes(warming, region)) +
    geom_col(width = 0.65) +
    geom_text(aes(label = paste0("n=", n)), hjust = -0.1, size = 5.5) +
    scale_x_continuous(limits = c(0, 1.30 * max(reg$warming)), expand = expansion(mult = c(0, 0.02))) +
    labs(x = paste0("Mean common change (", DEG_C, ")"), y = NULL) +
    theme_paper()
}

# Figure 7: likelihood-ratio statistics for the memory parameter
plot_memory <- function(prof) {
  ggplot(prof, aes(d, LR)) +
    geom_hline(yintercept = qchisq(0.95, 1), colour = GRAY, linewidth = 0.4) +
    geom_vline(xintercept = 1, colour = GRAY, linewidth = 0.3, linetype = "dotted") +
    geom_line(linewidth = 0.8, colour = BLUE) +
    facet_wrap(~ factor) +
    coord_cartesian(xlim = c(0, 2.5), ylim = c(0, 12)) +
    labs(x = expression("memory parameter " * d), y = "likelihood-ratio statistic") +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(), strip.text = element_text(face = "bold"))
}

# ---- Maps (Figures 5, 6(a) and 9) ----

base_europe_map <- function() {
  bbox <- sf::st_as_sfc(sf::st_bbox(c(xmin = -25, xmax = 45, ymin = 30, ymax = 72),
                                    crs = sf::st_crs(4326)))
  countries <- suppressWarnings(sf::st_crop(
    rnaturalearth::ne_countries(scale = "medium", returnclass = "sf"), bbox))
  ggplot() +
    geom_sf(data = countries, fill = "grey96", colour = "grey65", linewidth = 0.25) +
    coord_sf(xlim = c(EUROPE_BBOX$lon_min, EUROPE_BBOX$lon_max),
             ylim = c(EUROPE_BBOX$lat_min, EUROPE_BBOX$lat_max), expand = FALSE) +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 13) +
    theme(panel.grid.major = element_line(colour = "grey90", linewidth = 0.2),
          panel.grid.minor = element_blank(), axis.text = element_text(size = 18))
}

as_points <- function(df) sf::st_as_sf(df, coords = c("lon", "lat"), crs = 4326, remove = FALSE)

size_breaks <- function(x, n = 4) {
  br <- pretty(range(x), n = n)
  br <- br[br >= min(x) & br <= max(x)]
  if (length(br) == 0L) range(x) else unique(br)
}

size_guide <- function() {
  guide_legend(order = 2, title.position = "top",
               override.aes = list(fill = "grey80", colour = "grey35", alpha = 1))
}

# Figure 5: sign (colour) and absolute size of the loadings
plot_loading_map <- function(df, value) {
  df$sign <- ifelse(df[[value]] >= 0, "positive", "negative")
  df$abs_value <- abs(df[[value]])
  base_europe_map() +
    geom_sf(data = as_points(df), aes(fill = sign, size = abs_value),
            shape = 21, colour = "white", stroke = 0.25, alpha = 0.95) +
    scale_fill_manual(values = c(positive = "#B2182B", negative = "#0070B2")) +
    scale_size_continuous(range = c(1.2, 5.5), breaks = size_breaks(df$abs_value)) +
    theme(legend.position = "none")
}

# Figure 6(a): station-specific common warming
plot_warming_map <- function(df, value) {
  df$value <- df[[value]]
  df$abs_value <- abs(df$value)
  v_max <- max(df$value)
  base_europe_map() +
    geom_sf(data = as_points(df), aes(fill = value, size = abs_value),
            shape = 21, colour = "grey25", stroke = 0.25, alpha = 0.95) +
    scale_fill_gradientn(colours = c("#FDD835", "#F9844A", "#D7191C"),
                         values = scales::rescale(c(1, mean(c(1, v_max)), v_max)),
                         limits = range(df$value), oob = scales::squish,
                         name = paste0("Change (", DEG_C, ")")) +
    scale_size_continuous(name = paste0("Absolute change (", DEG_C, ")"), range = c(1.4, 6.0),
                          breaks = size_breaks(df$abs_value)) +
    guides(fill = guide_colorbar(order = 1, title.position = "top", barheight = grid::unit(4, "cm")),
           size = size_guide()) +
    theme(
      legend.title = element_text(size = 18),
      legend.text = element_text(size = 16)
    )
}

# Figure 9: station warming after differencing minus warming in levels
plot_difference_map <- function(df) {
  df$abs_diff <- abs(df$diff)
  m <- max(df$abs_diff)
  base_europe_map() +
    geom_sf(data = as_points(df), aes(fill = diff, size = abs_diff),
            shape = 21, colour = "grey25", stroke = 0.25, alpha = 0.95) +
    scale_fill_gradient2(low = BLUE, mid = "grey97", high = RED, midpoint = 0, limits = c(-m, m),
                         name = paste0("Differenced minus levels (", DEG_C, ")")) +
    scale_size_continuous(name = paste0("Absolute value (", DEG_C, ")"), range = c(1.0, 5.0),
                          breaks = size_breaks(df$abs_diff)) +
    guides(fill = guide_colorbar(order = 1, title.position = "top", barheight = grid::unit(4, "cm")),
           size = size_guide()) +
    theme(
      axis.text = element_text(size = 13)
    )
}
