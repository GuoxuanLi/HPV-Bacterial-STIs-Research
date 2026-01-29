library(ggplot2)
library(dplyr)
library(tibble)
library(patchwork)
library(grid)
library(cowplot)


axis_min <- 0.01
axis_max <- 10


y_offset <- 0.18


x_or_col <- axis_max * 3       
x_p_col  <- axis_max * 20      
right_margin_pt <- 300        


cap_h <- 0.30


na_or_hjust <- 0.5
na_or_vjust <- 0.5
na_or_dx    <- 0     
na_or_dy    <- 0     
na_or_size  <- 6


na_p_hjust  <- 0.5
na_p_vjust  <- 0.5
na_p_dx     <- 0
na_p_dy     <- 0
na_p_size   <- 6


age_note_text <- "Adjusted for age"
age_note_x    <- 0.04
age_note_y    <- 0.85
age_note_size <- 18
age_note_face <- "bold"

add_age_note <- function(p,
                         x = age_note_x, y = age_note_y,
                         size = age_note_size, face = age_note_face,
                         text = age_note_text) {
  cowplot::ggdraw(p) +
    cowplot::draw_label(
      text,
      x = x, y = y,
      hjust = 0, vjust = 1,
      size = size,
      fontface = face
    )
}


prep_forest_data <- function(df, group_levels, y_label) {
  df %>%
    mutate(
      Group = factor(Group, levels = group_levels, ordered = TRUE),
      Pathogen = factor(
        Pathogen,
        levels = c("CT (Chlamydia trachomatis)",
                   "UU (Ureaplasma urealyticum)")
      ),
      
    
      y_base = as.numeric(Group),
      y_plot = y_base + case_when(
        Pathogen == "CT (Chlamydia trachomatis)"  ~  y_offset,
        Pathogen == "UU (Ureaplasma urealyticum)" ~ -y_offset,
        TRUE ~ 0
      ),
      
     
      OR_plot    = if_else(is.na(OR), NA_real_, pmin(pmax(OR, axis_min), axis_max)),
      Lower_plot = if_else(is.na(Lower), NA_real_, pmax(Lower, axis_min)),
      Upper_plot = if_else(is.na(Upper), NA_real_, pmin(Upper, axis_max)),
      
      
      OR_label = if_else(is.na(OR), "NA", sprintf("%.2f (%.2f–%.2f)", OR, Lower, Upper)),
      p_num = case_when(
        is.na(P_val)  ~ "NA",
        P_val < 0.001 ~ "<0.001",
        TRUE          ~ sprintf("%.3f", P_val)
      ),
      p_star = case_when(
        is.na(P_val)  ~ "",
        P_val < 0.001 ~ "***",
        P_val < 0.01  ~ "**",
        P_val < 0.05  ~ "*",
        TRUE          ~ ""
      ),
      P_label  = trimws(paste(p_num, p_star)),
      fontface = case_when(is.na(P_val) ~ "plain", P_val < 0.05 ~ "bold", TRUE ~ "plain"),
      Y_label  = y_label
    )
}


plot_forest_panel <- function(df_prepped) {
  glev     <- levels(df_prepped$Group)
  n_group  <- length(glev)
  y_header <- n_group + 0.7
  y_lab    <- unique(df_prepped$Y_label)
  
  df_draw  <- df_prepped %>% filter(!is.na(OR))
  data_na_or <- df_prepped %>% filter(is.na(OR))
  data_na_p  <- df_prepped %>% filter(is.na(P_val))
  
  ggplot(df_prepped, aes(x = OR_plot, y = y_plot, colour = Pathogen)) +
    
    
    geom_errorbarh(
      data = df_draw,
      aes(xmin = Lower_plot, xmax = Upper_plot),
      height    = 0.15,
      linewidth = 0.7
    ) +
    geom_point(
      data = df_draw,
      size  = 3,
      shape = 15
    ) +
    
   
    geom_vline(
      xintercept = 1,
      linetype   = "dashed",
      colour     = "grey40",
      linewidth  = 0.6
    ) +
    

  geom_text(
    data = df_prepped %>% filter(!is.na(OR)),
    aes(x = x_or_col, label = OR_label),
    hjust       = 0,
    size        =5,
    colour      = "black",
    show.legend = FALSE
  ) +

    geom_text(
      data = data_na_or,
      aes(x = x_or_col + na_or_dx, y = y_base + na_or_dy, label = OR_label),
      hjust       = -2.2,
      vjust       = -0.4,
      size        = 5,
      colour      = "black",
      show.legend = FALSE,
      inherit.aes = FALSE
    ) +
    

  geom_text(
    data = df_prepped %>% filter(!is.na(P_val)),
    aes(x = x_p_col, label = P_label, fontface = fontface),
    hjust       = 0,
    size        = 5,
    colour      = "black",
    show.legend = FALSE
  ) +
 
    geom_text(
      data = data_na_p,
      aes(x = x_p_col + na_p_dx, y = y_base + na_p_dy, label = P_label, fontface = fontface),
      hjust       = -0.5,
      vjust       = -0.4,
      size        = 5,
      colour      = "black",
      show.legend = FALSE,
      inherit.aes = FALSE
    ) +
    
 
    annotate("text", x = x_or_col, y = y_header, label = "OR (95% CI)",
             hjust = 0, vjust = 0, size = 6, fontface = "bold") +
    annotate("text", x = x_p_col,  y = y_header, label = "P value",
             hjust = 0.15, vjust = 0, size = 6, fontface = "bold") +
    
 
    scale_colour_manual(
      name   = NULL,
      values = c(
        "CT (Chlamydia trachomatis)"  = "#20854E",
        "UU (Ureaplasma urealyticum)" = "#E18727"
      ),
      breaks = c("CT (Chlamydia trachomatis)",
                 "UU (Ureaplasma urealyticum)"),
      labels = c(
        "CT (Chlamydia trachomatis)"  = "C. trachomatis",
        "UU (Ureaplasma urealyticum)" = "U. urealyticum"
      )
    ) +
    guides(colour = guide_legend(override.aes = list(size = 5))) +
    
 
    scale_x_continuous(
      trans  = "log10",
      breaks = c(0.01, 0.05, 0.1, 0.5, 1, 2, 5, 10),
      labels = c("0.01", "0.05", "0.1", "0.5", "1", "2", "5", "10"),
      name   = NULL
    ) +
    
   
    scale_y_continuous(
      breaks = seq_len(n_group),
      labels = glev,
      expand = expansion(mult = c(0.02, 0.25))
    ) +
    
    labs(y = y_lab) +
    theme_classic(base_size = 12) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text        = element_text(size = 17),
      plot.margin      = margin(t = 10, r = right_margin_pt, b = 5, l = 5)
    ) +
    coord_cartesian(xlim = c(axis_min, axis_max), clip = "off")
}

cyto_raw <- tribble(
  ~Pathogen, ~Group,                      ~OR,      ~Lower,   ~Upper,    ~P_val,
  "CT (Chlamydia trachomatis)",  "HSIL",   1.89,     0.102,   10.065,    0.547,
  "CT (Chlamydia trachomatis)",  "LSIL",   0.668,    0.037,    3.286,    0.696,
  "CT (Chlamydia trachomatis)",  "ASC-US", 3.037,    1.184,    6.853,    0.012,
  "CT (Chlamydia trachomatis)",  "ASC-H",  NA_real_, NA_real_, NA_real_, NA_real_,
  "CT (Chlamydia trachomatis)",  "Overall abnormal cytology", 1.961, 0.86, 4.063, 0.086,
  
  "UU (Ureaplasma urealyticum)", "HSIL",   0.82,     0.258,    2.541,    0.729,
  "UU (Ureaplasma urealyticum)", "LSIL",   1.562,    0.775,    3.307,    0.224,
  "UU (Ureaplasma urealyticum)", "ASC-US", 1.119,    0.652,    1.934,    0.685,
  "UU (Ureaplasma urealyticum)", "ASC-H",  1.045,    0.187,    5.856,    0.958,
  "UU (Ureaplasma urealyticum)", "Overall abnormal cytology", 1.177, 0.787, 1.769, 0.429
)

cyto_levels <- c("Overall abnormal cytology", "HSIL", "ASC-H", "LSIL", "ASC-US")
cyto_df <- prep_forest_data(cyto_raw, cyto_levels, "Cytology group")

p_cyto <- plot_forest_panel(cyto_df) +
  theme(
    legend.position      = c(0.03, 1),
    legend.justification = c(0, 1),
    legend.text          = element_text(size = 16),
    legend.key.size      = grid::unit(0.6, "cm"),
    axis.title.y         = element_text(size = 17, face = "bold")
  )
p_cyto <- add_age_note(p_cyto)


histo_raw <- tribble(
  ~Pathogen, ~Group,                       ~OR,      ~Lower,   ~Upper,    ~P_val,
  "CT (Chlamydia trachomatis)",  "CIN1",    1.947,    0.715,    4.503,    0.149,
  "CT (Chlamydia trachomatis)",  "CIN2+",   0.539,    0.03,     2.609,    0.547,
  "CT (Chlamydia trachomatis)",  "SCC",     NA_real_, NA_real_, NA_real_, NA_real_,
  "CT (Chlamydia trachomatis)",  "Overall abnormal histology", 1.306, 0.521, 2.851, 0.532,
  
  "UU (Ureaplasma urealyticum)", "CIN1",    1.705,    1.028,    2.888,    0.042,
  "UU (Ureaplasma urealyticum)", "CIN2+",   1.071,    0.563,    2.063,    0.835,
  "UU (Ureaplasma urealyticum)", "SCC",     0.855,    0.206,    3.339,    0.820,
  "UU (Ureaplasma urealyticum)", "Overall abnormal histology", 1.372, 0.925, 2.052, 0.119
)

histo_levels <- c("Overall abnormal histology", "SCC", "CIN2+", "CIN1")
histo_df <- prep_forest_data(histo_raw, histo_levels, "Pathology group")

p_histo <- plot_forest_panel(histo_df) +
  theme(
    legend.position      = c(0.03, 1),
    legend.justification = c(0, 1),
    legend.text          = element_text(size = 16),
    legend.key.size      = grid::unit(0.6, "cm"),
    axis.title.y         = element_text(size = 17, face = "bold")
  )
p_histo <- add_age_note(p_histo)


p_all <- (p_cyto / p_histo) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 18, face = "bold"))

p_all
