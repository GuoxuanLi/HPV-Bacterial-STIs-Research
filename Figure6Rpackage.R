

library(readxl)
library(dplyr)
library(tibble)
library(writexl)
library(ggplot2)
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
na_or_size  <- 5.5

na_p_hjust  <- 0.5
na_p_vjust  <- 0.5
na_p_dx     <- 0
na_p_dy     <- 0
na_p_size   <- 5.5


unadj_text <- "Unadjusted"
unadj_x    <- 0.05   
unadj_y    <- 0.88   
unadj_size <- 18
unadj_face <- "bold"

add_unadj_note <- function(p,
                           x = unadj_x, y = unadj_y,
                           size = unadj_size, face = unadj_face,
                           text = unadj_text) {
  cowplot::ggdraw(p) +
    cowplot::draw_label(
      text,
      x = x, y = y,
      hjust = 0, vjust = 1,
      size = size,
      fontface = face
    )
}


to01 <- function(x){
  x2 <- trimws(as.character(x))
  ifelse(x2 %in% c("阳性","positive","POS","1"), 1L,
         ifelse(x2 %in% c("阴性","negative","NEG","0"), 0L, NA_integer_))
}


fit_one_logit <- function(dat, var_name, pathogen_label, group_label){

  
  if (length(unique(dat$结局)) < 2 ||
      length(unique(dat[[var_name]])) < 2) {
    return(tibble(
      Pathogen = pathogen_label,
      Group    = group_label,
      OR       = NA_real_,
      Lower    = NA_real_,
      Upper    = NA_real_,
      P_val    = NA_real_
    ))
  }
  
  dat2 <- dat %>%
    mutate(暴露 = .data[[var_name]])
  
  fit <- glm(结局 ~ 暴露,
             family = binomial(link = "logit"),
             data   = dat2)
  
  sm <- summary(fit)$coefficients
  

  ci_vec <- tryCatch(
    {
      ci <- suppressMessages(confint(fit, parm = "暴露"))
      c(ci[1], ci[2])
    },
    error = function(e) c(NA_real_, NA_real_)
  )
  
  OR     <- unname(exp(coef(fit)["暴露"]))
  CI_low <- exp(ci_vec[1])
  CI_up  <- exp(ci_vec[2])
  p_val  <- sm["暴露", "Pr(>|z|)"]
  
  tibble(
    Pathogen = pathogen_label,
    Group    = group_label,
    OR       = OR,
    Lower    = CI_low,
    Upper    = CI_up,
    P_val    = p_val
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
      

      valid = !(is.na(OR) | is.na(Lower) | is.na(Upper) |
                  !is.finite(OR) | !is.finite(Lower) | !is.finite(Upper) |
                  Lower <= 0),   
      
      y_base = as.numeric(Group),
      y_plot = y_base + case_when(
        Pathogen == "CT (Chlamydia trachomatis)"  ~  y_offset,
        Pathogen == "UU (Ureaplasma urealyticum)" ~ -y_offset,
        TRUE ~ 0
      ),
      
 
      OR_plot    = if_else(valid, pmin(pmax(OR, axis_min), axis_max), NA_real_),
      Lower_plot = if_else(valid, pmax(Lower, axis_min),                NA_real_),
      Upper_plot = if_else(valid, pmin(Upper, axis_max),                NA_real_),
      
 
      OR_label = if_else(
        valid,
        sprintf("%.2f (%.2f–%.2f)", OR, Lower, Upper),
        "NA"
      ),
      
      p_num = case_when(
        !valid        ~ "NA",          
        is.na(P_val)  ~ "NA",
        P_val < 0.001 ~ "<0.001",
        TRUE          ~ sprintf("%.3f", P_val)
      ),
      p_star = case_when(
        !valid        ~ "",
        is.na(P_val)  ~ "",
        P_val < 0.001 ~ "***",
        P_val < 0.01  ~ "**",
        P_val < 0.05  ~ "*",
        TRUE          ~ ""
      ),
      
      P_label  = trimws(paste(p_num, p_star)),
      fontface = case_when(
        !valid       ~ "plain",
        is.na(P_val) ~ "plain",
        P_val < 0.05 ~ "bold",
        TRUE         ~ "plain"
      ),
      
      Y_label = y_label
    )
}



plot_forest_panel <- function(df_prepped) {
  
  glev     <- levels(df_prepped$Group)
  n_group  <- length(glev)
  y_header <- n_group + 0.7
  y_lab    <- unique(df_prepped$Y_label)
  
  df_draw    <- df_prepped %>% filter(!is.na(OR))
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
      size        = 5,
      colour      = "black",
      show.legend = FALSE
    ) +
 
    geom_text(
      data = data_na_or,
      aes(x = x_or_col + na_or_dx, y = y_base + na_or_dy, label = OR_label),
      hjust       = -2.3,
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
      aes(x = x_p_col + na_p_dx, y = y_base + na_p_dy,
          label = P_label, fontface = fontface),
      hjust       = -0.4,
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



file_in <- "XX"
df <- read_excel(file_in, sheet = 1)



normal_cyto <- "未见上皮内病变细胞和恶性细胞"
exclude_cyto <- c("含16、18高危型，不做细胞学检查", "HPV阴性未做")


cyto_map <- tribble(
  ~Group,                      ~Cyto_label,
  "HSIL",                      "高度鳞状上皮内病变（HSIL）",
  "ASC-H",                     "不典型鳞状上皮细胞-不除外高度鳞状上皮内病变（ASC-H）",
  "LSIL",                      "低度鳞状上皮内病变（LSIL）",
  "ASC-US",                    "未明确意义的不典型鳞状上皮细胞（ASC-US）"
)

df_cyto <- df %>%
  mutate(
    HPV阳性 = if_else(`HPV亚型` == "HPV阴性" | is.na(`HPV亚型`), 0L, 1L),
    CT阳性  = to01(CT),
    UU阳性  = to01(UU)
  ) %>%
  filter(
    HPV阳性 == 1,
    !is.na(CT阳性),
    !is.na(UU阳性),
    !is.na(细胞学),
    !细胞学 %in% exclude_cyto
  )


pathogen_list <- tribble(
  ~Pathogen,                         ~var_name,
  "CT (Chlamydia trachomatis)",     "CT阳性",
  "UU (Ureaplasma urealyticum)",    "UU阳性"
)


run_cyto_overall <- function(pathogen, var_name){
  dat <- df_cyto %>%
    mutate(
      结局 = if_else(细胞学 == normal_cyto, 0L, 1L)
    )
  fit_one_logit(dat, var_name, pathogen, "Overall abnormal cytology")
}


run_cyto_each <- function(pathogen, var_name){
  purrr::map_dfr(
    1:nrow(cyto_map),
    function(i){
      grp_eng <- cyto_map$Group[i]
      grp_ch  <- cyto_map$Cyto_label[i]
      
      dat <- df_cyto %>%
        filter(细胞学 %in% c(normal_cyto, grp_ch)) %>%
        mutate(
          结局 = if_else(细胞学 == grp_ch, 1L, 0L)
        )
      
      fit_one_logit(dat, var_name, pathogen, grp_eng)
    }
  )
}

cyto_raw <- purrr::map_dfr(
  1:nrow(pathogen_list),
  function(i){
    pathogen <- pathogen_list$Pathogen[i]
    var_name <- pathogen_list$var_name[i]
    
    bind_rows(
      run_cyto_overall(pathogen, var_name),
      run_cyto_each(pathogen, var_name)
    )
  }
)


df_histo <- df %>%
  mutate(
    HPV阳性 = if_else(`HPV亚型` == "HPV阴性" | is.na(`HPV亚型`), 0L, 1L),
    CT阳性  = to01(CT),
    UU阳性  = to01(UU)
  ) %>%
  filter(
    HPV阳性 == 1,
    !is.na(CT阳性),
    !is.na(UU阳性),
    !is.na(病理学),
    病理学 != "HPV阴性未做"
  )



run_histo_overall <- function(pathogen, var_name){
  dat <- df_histo %>%
    filter(病理学 %in% c("CIN0", "CIN1", "CIN2+", "SCC")) %>%
    mutate(
      结局 = if_else(病理学 == "CIN0", 0L, 1L)
    )
  fit_one_logit(dat, var_name, pathogen, "Overall abnormal histology")
}

run_histo_each <- function(pathogen, var_name){
  histo_groups <- c("CIN2+", "CIN1", "SCC")
  
  purrr::map_dfr(
    histo_groups,
    function(g){
      dat <- df_histo %>%
        filter(病理学 %in% c("CIN0", g)) %>%
        mutate(
          结局 = if_else(病理学 == g, 1L, 0L)
        )
      
      fit_one_logit(dat, var_name, pathogen, g)
    }
  )
}

histo_raw <- purrr::map_dfr(
  1:nrow(pathogen_list),
  function(i){
    pathogen <- pathogen_list$Pathogen[i]
    var_name <- pathogen_list$var_name[i]
    
    bind_rows(
      run_histo_overall(pathogen, var_name),
      run_histo_each(pathogen, var_name)
    )
  }
)



cyto_levels  <- c("Overall abnormal cytology", "HSIL", "ASC-H", "LSIL", "ASC-US")
cyto_df      <- prep_forest_data(cyto_raw, cyto_levels, "Cytology group")
p_cyto <- plot_forest_panel(cyto_df) +
  theme(
    legend.position      = c(0.03, 1),
    legend.justification = c(0, 1),
    legend.text          = element_text(size = 16),
    legend.key.size      = grid::unit(0.6, "cm"),
    axis.title.y         = element_text(size = 17, face = "bold")
  )
p_cyto <- add_unadj_note(p_cyto)

histo_levels <- c("Overall abnormal histology", "SCC", "CIN2+", "CIN1")
histo_df     <- prep_forest_data(histo_raw, histo_levels, "Pathology group")
p_histo <- plot_forest_panel(histo_df) +
  theme(
    legend.position      = c(0.03, 1),
    legend.justification = c(0, 1),
    legend.text          = element_text(size = 16),
    legend.key.size      = grid::unit(0.6, "cm"),
    axis.title.y         = element_text(size = 17, face = "bold")
  )
p_histo <- add_unadj_note(p_histo)

p_all <- (p_cyto / p_histo) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(size = 18, face = "bold"),
    plot.tag.position = c(0.01, 0.99)
  )

p_all
