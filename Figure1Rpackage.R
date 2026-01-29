
library(ggplot2)
library(grid)


df <- data.frame(
  Infection_pattern = factor(
    c("Single-type infection",
      "Double-type co-infection",
      "Triple-type co-infection",
      "≥Quadruple-type co-infection"),
    levels = c("Single-type infection",
               "Double-type co-infection",
               "Triple-type co-infection",
               "≥Quadruple-type co-infection")
  ),
  Positive_samples = c(722, 115, 28, 4),
  Positive_rate    = c(83.08, 13.23, 3.22, 0.46)
)

step_left  <- 50
step_right <- 5

left_raw   <- ceiling(max(df$Positive_samples) / step_left)  * step_left
right_raw  <- ceiling(max(df$Positive_rate)    / step_right) * step_right

left_limit  <- left_raw  + 2 * step_left
right_limit <- right_raw + 2 * step_right

scale_factor <- left_limit / right_limit

SCALE <- 0.75   

pdf_file <- "Figure_dual_axis_WIDE.pdf"
pdf_w    <- 14 * SCALE
pdf_h    <- 6  * SCALE


BASE_SIZE       <- 18  * SCALE
AXIS_TEXT_SIZE  <- 25  * SCALE
AXIS_TITLE_SIZE <- 25  * SCALE
LEGEND_TEXT_SZ  <- 25  * SCALE

PCT_LABEL_SIZE  <- 8 * SCALE
PCT_NUDGE_X     <- 0.12
PCT_NUDGE_Y     <- 12

YL_VJUST  <- 1.0
YR_VJUST  <- 1.0
YL_MARGIN <- 40     
YR_MARGIN <- 40     


BAR_COL   <- "#F4B5BD"
LINE_COL  <- "#AD002A"
PCT_COL   <- "black"   


p <- ggplot(df, aes(x = Infection_pattern)) +
  geom_col(aes(y = Positive_samples, fill = "Number of positive samples"),
           width = 0.6) +
  geom_line(aes(y = Positive_rate * scale_factor, color = "Positive rate (%)"),
            group = 1, linewidth = 1.1) +
  geom_point(aes(y = Positive_rate * scale_factor, color = "Positive rate (%)"),
             size = 3.0) +
  geom_text(
    aes(y = Positive_rate * scale_factor,
        label = sprintf("%.2f%%", Positive_rate)),
    color = PCT_COL,
    nudge_x = PCT_NUDGE_X,
    nudge_y = PCT_NUDGE_Y,
    size  = PCT_LABEL_SIZE,
    show.legend = FALSE
  ) +
  scale_y_continuous(
    name   = "Number of positive samples",
    limits = c(0, left_limit),
    breaks = seq(0, left_limit, by = step_left),
    expand = c(0, 0),
    sec.axis = sec_axis(
      ~ . / scale_factor,
      name   = "Positive rate (%)",
      breaks = seq(0, right_limit, by = step_right)
    )
  ) +
  scale_fill_manual(
    name   = NULL,
    values = c("Number of positive samples" = BAR_COL)
  ) +
  scale_color_manual(
    name   = NULL,
    values = c("Positive rate (%)" = LINE_COL)
  ) +
  labs(x = NULL) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = BASE_SIZE) +
  theme(
    axis.text.x       = element_text(size = AXIS_TEXT_SIZE, angle = 0, vjust = 0.5),
    axis.text.y       = element_text(size = AXIS_TEXT_SIZE),
    axis.text.y.right = element_text(size = AXIS_TEXT_SIZE),
    
    
    axis.title.y = element_text(
      size = AXIS_TITLE_SIZE,
      vjust = YL_VJUST,
      margin = margin(r = YL_MARGIN)
    ),
    axis.title.y.right = element_text(
      size = AXIS_TITLE_SIZE,
      angle = 90,
      vjust = YR_VJUST,
      margin = margin(l = YR_MARGIN)
    ),
    
    legend.position      = c(0.03, 1),
    legend.justification = c("left", "top"),
    legend.background    = element_rect(fill = "transparent", color = NA),
    legend.text          = element_text(size = LEGEND_TEXT_SZ),
    
    plot.margin          = margin(t = 8, r = 60, b = 12, l = 18)
  )


print(p)


ggsave(pdf_file, plot = p, device = cairo_pdf,
       width = pdf_w, height = pdf_h, units = "in")
