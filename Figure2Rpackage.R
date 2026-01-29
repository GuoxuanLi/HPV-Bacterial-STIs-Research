library(ggplot2)
library(dplyr)
library(tibble)
library(patchwork)
library(readxl)
library(stringr)

file_path <- "XX"
sheet_id  <- 1

col_hpv <- "HPV亚型"
col_ct  <- "CT"
col_uu  <- "UU"
col_ng  <- "NG"
col_age <- "年龄"

HPV_NEG_PATTERN  <- "阴性|未做"
PATH_POS_PATTERN <- "阳性"

EVENTS_FISHER_CUTOFF <- 5


parse_age_num <- function(x){
  if(is.numeric(x)) return(as.numeric(x))
  x <- str_squish(as.character(x))
  x[x==""] <- NA
  suppressWarnings(as.numeric(str_extract(x, "\\d+\\.?\\d*")))
}

df <- read_excel(file_path, sheet = sheet_id)

dat <- df %>%
  mutate(
    HPV_raw = str_squish(as.character(.data[[col_hpv]])),
    CT_raw  = str_squish(as.character(.data[[col_ct]])),
    UU_raw  = str_squish(as.character(.data[[col_uu]])),
    NG_raw  = str_squish(as.character(.data[[col_ng]])),
    Age_raw = .data[[col_age]],
    Age     = parse_age_num(Age_raw),
    
    HPV_pos = ifelse(is.na(HPV_raw) | HPV_raw=="" | str_detect(HPV_raw, HPV_NEG_PATTERN), 0, 1),
    
    CT_pos  = ifelse(!is.na(CT_raw) & str_detect(CT_raw, PATH_POS_PATTERN), 1, 0),
    UU_pos  = ifelse(!is.na(UU_raw) & str_detect(UU_raw, PATH_POS_PATTERN), 1, 0),
    NG_pos  = ifelse(!is.na(NG_raw) & str_detect(NG_raw, PATH_POS_PATTERN), 1, 0)
  )

cat("==== 数据检查 ====\n")
cat("总人数 n =", nrow(dat), "\n")
cat("HPV+ =", sum(dat$HPV_pos==1, na.rm=TRUE), "；HPV- =", sum(dat$HPV_pos==0, na.rm=TRUE), "\n")
cat("Age缺失 =", sum(is.na(dat$Age)), "\n\n")

rate_pct <- function(x) round(mean(x == 1, na.rm = TRUE) * 100, 2)

df_group <- tibble(
  Pathogen = factor(rep(c("CT", "UU", "NG"), each = 2), levels = c("CT","UU","NG")),
  Group    = rep(c("HPV Positive", "HPV Negative"), times = 3),
  Rate     = c(
    rate_pct(dat %>% filter(HPV_pos==1) %>% pull(CT_pos)),
    rate_pct(dat %>% filter(HPV_pos==0) %>% pull(CT_pos)),
    rate_pct(dat %>% filter(HPV_pos==1) %>% pull(UU_pos)),
    rate_pct(dat %>% filter(HPV_pos==0) %>% pull(UU_pos)),
    rate_pct(dat %>% filter(HPV_pos==1) %>% pull(NG_pos)),
    rate_pct(dat %>% filter(HPV_pos==0) %>% pull(NG_pos))
  )
)

df_overall <- tibble(
  Pathogen = factor(c("CT","UU","NG"), levels = c("CT","UU","NG")),
  Group    = "Overall",
  Rate     = c(rate_pct(dat$CT_pos), rate_pct(dat$UU_pos), rate_pct(dat$NG_pos))
)

df_plot <- bind_rows(df_overall, df_group) %>%
  mutate(Group = factor(Group, levels = c("Overall", "HPV Positive", "HPV Negative")))


get_p_age_adjusted <- function(y_col){
  d0 <- dat %>%
    transmute(
      Y = .data[[y_col]],
      HPV_pos = factor(HPV_pos, levels = c(0,1)),  
      Age = Age
    ) %>%
    filter(!is.na(Y), !is.na(HPV_pos), !is.na(Age))
  
  if(nrow(d0) == 0) return(list(p=NA_real_, method="Age-adjusted (no data)", events=0))
  if(length(unique(d0$Y)) < 2) return(list(p=NA_real_, method="Age-adjusted (no variation)", events=sum(d0$Y==1)))
  
  events <- sum(d0$Y == 1, na.rm = TRUE)
  
  
  if(requireNamespace("brglm2", quietly = TRUE)){
    fit <- tryCatch(
      glm(Y ~ HPV_pos + Age, data=d0, family=binomial("logit"), method="brglmFit"),
      error = function(e) NULL
    )
    if(!is.null(fit)){
      sm <- summary(fit)$coefficients
      rn <- grep("^HPV_pos", rownames(sm), value = TRUE)[1]  
      cn <- grep("^Pr\\(", colnames(sm), value = TRUE)[1]
      if(!is.na(rn) && !is.na(cn)){
        return(list(p=as.numeric(sm[rn, cn]), method="Age-adjusted (brglmFit)", events=events))
      }
      return(list(p=NA_real_, method="Age-adjusted (term missing)", events=events))
    }
    return(list(p=NA_real_, method="Age-adjusted (fit failed)", events=events))
  }
  
 
  fit2 <- tryCatch(glm(Y ~ HPV_pos + Age, data=d0, family=binomial("logit")), error=function(e) NULL)
  if(is.null(fit2)) return(list(p=NA_real_, method="Age-adjusted glm (fit failed)", events=events))
  sm2 <- summary(fit2)$coefficients
  rn2 <- grep("^HPV_pos", rownames(sm2), value = TRUE)[1]
  cn2 <- grep("^Pr\\(", colnames(sm2), value = TRUE)[1]
  if(is.na(rn2) || is.na(cn2)) return(list(p=NA_real_, method="Age-adjusted glm (term missing)", events=events))
  list(p=as.numeric(sm2[rn2, cn2]), method="Age-adjusted (glm)", events=events)
}

get_p_chisq <- function(y_col){
  d1 <- dat %>%
    transmute(Y = .data[[y_col]], HPV_pos = HPV_pos) %>%
    filter(!is.na(Y), !is.na(HPV_pos))
  
  if(nrow(d1) == 0) return(list(p=NA_real_, method="Pearson χ² (no data)", events=0))
  
  # 固定成 2x2（避免缺列/缺行）
  tab <- matrix(0, nrow = 2, ncol = 2,
                dimnames = list(HPV_pos = c("0","1"), Y = c("0","1")))
  tt <- table(factor(d1$HPV_pos, levels = c(0,1)),
              factor(d1$Y,       levels = c(0,1)))
  tab[rownames(tt), colnames(tt)] <- tt
  
  events <- sum(d1$Y == 1, na.rm = TRUE)
  
  p <- tryCatch(
    suppressWarnings(chisq.test(tab, correct = FALSE)$p.value),  
    error = function(e) NA_real_
  )
  list(p = as.numeric(p), method = "Pearson χ²", events = events)
}


get_p_age_adjusted <- get_p_age_adjusted
get_p_fisher <- get_p_fisher
get_p_with_rule <- get_p_with_rule

tmp <- tibble(
  Pathogen = c("CT","UU","NG"),
  res = list(
    get_p_chisq("CT_pos"),      
    get_p_chisq("UU_pos"),      
    get_p_with_rule("NG_pos")   
  )
) %>%
  mutate(
    events = sapply(res, `[[`, "events"),
    method = sapply(res, `[[`, "method"),
    p_value = sapply(res, `[[`, "p"),
    Significance = case_when(
      is.na(p_value)   ~ "NA",
      p_value <= 0.001 ~ "***",
      p_value <= 0.01  ~ "**",
      p_value <= 0.05  ~ "*",
      TRUE             ~ "ns"
    )
  )

p_df <- tmp %>% select(Pathogen, events, method, p_value, Significance)
print(df_plot)
print(p_df)


p_df <- tmp %>% select(Pathogen, events, method, p_value, Significance)
print(df_plot)
print(p_df)


groupwidth <- 0.9
pd <- position_dodge2(width = groupwidth, padding = 0)

SHOW_PERCENT <- TRUE
LABEL_FMT <- function(x) sprintf("%.2f%%", x)

BAR_LABEL_SIZE <- 7
BAR_LABEL_NUDGE_UU   <- 1.2
BAR_LABEL_NUDGE_CTNG <- 0.18

BRACKET_LWD <- 1.3
STAR_SIZE   <- 8

LEGEND_KEY_SIZE <- 0.9
LEGEND_BOX_PAD  <- 0.3


group_levels <- levels(df_plot$Group)
n_group      <- length(group_levels)
offsets      <- -groupwidth/2 + ((1:n_group) - 0.5) * groupwidth / n_group
names(offsets) <- group_levels

make_bracket_df <- function(df_one_pathogen, p_df,
                            bottom_offset = 1.0,
                            height        = 1.0,
                            star_offset   = 0.8) {
  
  df_one <- df_one_pathogen
  patho  <- as.character(unique(df_one$Pathogen))
  center <- as.numeric(df_one$Pathogen[1])
  
  x_hpvp   <- center + offsets["HPV Positive"]
  x_hpvneg <- center + offsets["HPV Negative"]
  
  y_hpvp <- df_one %>% filter(Group == "HPV Positive") %>% pull(Rate)
  y_bottom <- y_hpvp + bottom_offset
  y_top    <- y_bottom + height
  y_star   <- y_top + star_offset
  x_mid <- (x_hpvp + x_hpvneg) / 2
  
  sig_lab <- p_df %>% filter(Pathogen == patho) %>% pull(Significance)
  
  tibble(Pathogen = patho,
         x_left   = x_hpvp,
         x_right  = x_hpvneg,
         x_mid    = x_mid,
         y_bottom = y_bottom,
         y_top    = y_top,
         y_star   = y_star,
         label    = sig_lab)
}


df_uu <- df_plot %>%
  filter(Pathogen == "UU") %>%
  mutate(Pathogen = factor(Pathogen))

bracket_uu <- make_bracket_df(
  df_one_pathogen = df_uu,
  p_df            = p_df,
  bottom_offset   = 3,
  height          = 1.5,
  star_offset     = 0.8
)

p_uu <- ggplot(df_uu, aes(x = Pathogen, y = Rate, fill = Group)) +
  geom_col(position = pd, width = groupwidth) +
  {if (SHOW_PERCENT) geom_text(
    aes(label = LABEL_FMT(Rate)),
    position = pd,
    vjust = -0.7,
    nudge_y = BAR_LABEL_NUDGE_UU,
    size = BAR_LABEL_SIZE
  )} +
  geom_segment(data = bracket_uu,
               aes(x = x_left, xend = x_left, y = y_bottom, yend = y_top),
               inherit.aes = FALSE, linewidth = BRACKET_LWD) +
  geom_segment(data = bracket_uu,
               aes(x = x_right, xend = x_right, y = y_bottom, yend = y_top),
               inherit.aes = FALSE, linewidth = BRACKET_LWD) +
  geom_segment(data = bracket_uu,
               aes(x = x_left, xend = x_right, y = y_top, yend = y_top),
               inherit.aes = FALSE, linewidth = BRACKET_LWD) +
  geom_text(data = bracket_uu,
            aes(x = x_mid, y = y_star, label = label),
            inherit.aes = FALSE,
            size = STAR_SIZE,
            fontface = "bold") +
  scale_fill_manual(
    name   = "Group",
    values = c("Overall"      = "#F4B5BD",
               "HPV Positive" = "#AD002A",
               "HPV Negative" = "#4A6F8A")
  ) +
  scale_y_continuous(
    name   = "Infection Rate (%)",
    limits = c(0, 60),
    expand = expansion(mult = c(0, 0.08))
  ) +
  labs(x = NULL) +
  scale_x_discrete(labels = "U. urealyticum") +
  theme_classic(base_size = 12) +
  theme(
    axis.title.y = element_text(size = 20),
    axis.text    = element_text(size = 20),
    axis.title.x = element_blank(),
    legend.position = "none"
  )


df_ctng <- df_plot %>%
  filter(Pathogen %in% c("CT", "NG")) %>%
  mutate(Pathogen = factor(Pathogen, levels = c("CT","NG")))

bracket_ct <- make_bracket_df(
  df_one_pathogen = df_ctng %>% filter(Pathogen == "CT"),
  p_df            = p_df,
  bottom_offset   = 0.4,
  height          = 0.15,
  star_offset     = 0.10
)

bracket_ng <- make_bracket_df(
  df_one_pathogen = df_ctng %>% filter(Pathogen == "NG"),
  p_df            = p_df,
  bottom_offset   = 0.4,
  height          = 0.15,
  star_offset     = 0.15
)

brackets_ctng <- bind_rows(bracket_ct, bracket_ng)

p_ctng <- ggplot(df_ctng, aes(x = Pathogen, y = Rate, fill = Group)) +
  geom_col(position = pd, width = groupwidth) +
  {if (SHOW_PERCENT) geom_text(
    aes(label = LABEL_FMT(Rate)),
    position = pd,
    vjust = -1,
    nudge_y = BAR_LABEL_NUDGE_CTNG,
    size = BAR_LABEL_SIZE
  )} +
  geom_segment(data = brackets_ctng,
               aes(x = x_left, xend = x_left, y = y_bottom, yend = y_top),
               inherit.aes = FALSE, linewidth = BRACKET_LWD) +
  geom_segment(data = brackets_ctng,
               aes(x = x_right, xend = x_right, y = y_bottom, yend = y_top),
               inherit.aes = FALSE, linewidth = BRACKET_LWD) +
  geom_segment(data = brackets_ctng,
               aes(x = x_left, xend = x_right, y = y_top, yend = y_top),
               inherit.aes = FALSE, linewidth = BRACKET_LWD) +
  geom_text(data = brackets_ctng,
            aes(x = x_mid, y = y_star, label = label),
            inherit.aes = FALSE,
            size = STAR_SIZE,
            fontface = "bold") +
  scale_fill_manual(
    name   = NULL,
    values = c("Overall"      = "#F4B5BD",
               "HPV Positive" = "#AD002A",
               "HPV Negative" = "#4A6F8A")
  ) +
  scale_y_continuous(
    name   = "Infection Rate (%)",
    limits = c(0, 6),
    expand = expansion(mult = c(0, 0.10))
  ) +
  labs(x = NULL) +
  scale_x_discrete(labels = c("CT" = "C. trachomatis",
                              "NG" = "N. gonorrhoeae")) +
  theme_classic(base_size = 12) +
  theme(
    axis.title.y = element_text(size = 20),
    axis.text    = element_text(size = 20),
    axis.title.x = element_blank(),
    legend.position      = c(0.98, 0.55),
    legend.justification = c(1, 1),
    legend.background    = element_rect(fill = "white", colour = NA),
    legend.text          = element_text(size = 20),
    legend.key.size      = unit(LEGEND_KEY_SIZE, "cm"),
    legend.box.margin    = margin(LEGEND_BOX_PAD, LEGEND_BOX_PAD, LEGEND_BOX_PAD, LEGEND_BOX_PAD, "cm")
  )

# =======================================
# 8) 组合图
# =======================================
(p_uu | p_ctng) + plot_layout(widths = c(1, 2))
p_df %>%
  mutate(
    p_value = ifelse(is.na(p_value), NA, signif(p_value, 3))
  ) %>%
  print(n = Inf)

