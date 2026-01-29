library(ggplot2)
library(dplyr)
library(tibble)

dumbbell_df <- tibble(
  Age       = factor(c("35–39", "40–44", "45–49", "50–54", "55–59", "60–64"),
                     levels = c("35–39", "40–44", "45–49", "50–54", "55–59", "60–64"),
                     ordered = TRUE),
  HPV_Rate  = c(58.8, 58.0, 63.5, 63.3,  72.9,  77.4),
  STI_Rate  = c(48.5, 52.5, 57.2, 46.2,  39.4,  34.0)
) %>%
  mutate(
    P_label_upper = c("Ref", "ns", "ns", "ns", "***", "***"),
    P_label_lower = c("Ref", "ns", "ns", "ns", "*", "**")
  ) %>%
  mutate(
    upper_y = pmax(HPV_Rate, STI_Rate) + 1,
    lower_y = pmin(HPV_Rate, STI_Rate) - 1
  )

y_max <- max(dumbbell_df$upper_y) + 1
y_min <- min(0, min(dumbbell_df$lower_y) - 1)


PCT_NUDGE_X <- 0.15   
PCT_SIZE    <- 6      
PCT_COL     <- "black"


Y_TITLE_VJUST <- 0.5   
Y_TITLE_RMARG <- 25    


p <- ggplot(dumbbell_df, aes(x = Age)) +
 
  geom_segment(
    aes(x = Age, xend = Age, y = HPV_Rate, yend = STI_Rate),
    colour = "black",
    linewidth = 0.7
  ) +
 
  geom_point(aes(y = HPV_Rate, colour = "HPV"), size = 4) +
  geom_point(aes(y = STI_Rate, colour = "STIs"), size = 4) +

  geom_text(
    aes(y = HPV_Rate, label = sprintf("%.1f%%", HPV_Rate)),
    nudge_x = PCT_NUDGE_X, hjust = 0,
    size = PCT_SIZE, color = PCT_COL
  ) +
  geom_text(
    aes(y = STI_Rate, label = sprintf("%.1f%%", STI_Rate)),
    nudge_x = PCT_NUDGE_X, hjust = 0,
    size = PCT_SIZE, color = PCT_COL
  ) +
 
  geom_text(
    aes(y = upper_y, label = P_label_upper),
    fontface = "bold", size = 8, vjust = -0.75
  ) +
  geom_text(
    aes(y = lower_y, label = P_label_lower),
    fontface = "bold", size = 8, vjust = 1.5
  ) +
 
  scale_colour_manual(
    name   = NULL,
    values = c("HPV" = "#E18727", "STIs" = "#20854E"),
    breaks = c("HPV", "STIs"),
    labels = c("HR-HPV positive",
               "C. trachomatis/U. urealyticum/N. gonorrhoeae positive"),
    guide = guide_legend(keyheight = unit(2, "lines"))
  ) +
  scale_y_continuous(
    name   = "Prevalence Rate (%)",
    limits = c(y_min, y_max),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(x = NULL) +
  theme_classic(base_size = 12) +
  theme(
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(
      size = 20,
      angle = 90,
      vjust = Y_TITLE_VJUST,
      margin = margin(r = 35)  
    ),
    axis.text    = element_text(size = 20),
    legend.position      = c(0.97, 0.2),
    legend.justification = c(1, 1),
    legend.text          = element_text(size = 20),
    legend.background    = element_rect(fill = "white", colour = NA)
  )

p
