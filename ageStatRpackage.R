pkgs <- c("readxl","dplyr","stringr","tidyr","tibble","purrr","openxlsx")
ins <- pkgs[!pkgs %in% installed.packages()[,"Package"]]
if(length(ins)>0) install.packages(ins)
invisible(lapply(pkgs, library, character.only = TRUE))


file_in   <- "XX"
sheet_id  <- 1                     
col_age   <- "年龄"                
col_hpv   <- "HPV亚型"             
file_out  <- "年龄_35-64_HPV阳性vs阴性_统计结果.xlsx"

df0 <- read_excel(file_in, sheet = sheet_id)
names(df0) <- str_squish(names(df0))

if(!all(c(col_age, col_hpv) %in% names(df0))){
  stop("找不到列名：请检查 col_age / col_hpv 是否与Excel第一行完全一致。当前列名有：\n",
       paste(names(df0), collapse = " | "))
}

df <- df0 %>%
  transmute(
    年龄_raw = .data[[col_age]],
    HPV亚型_raw = .data[[col_hpv]]
  ) %>%
  mutate(
    年龄 = suppressWarnings(as.numeric(as.character(年龄_raw))),
    HPV亚型 = str_squish(as.character(HPV亚型_raw)),
    HPV组 = case_when(
      is.na(HPV亚型) ~ NA_character_,
      HPV亚型 == "" ~ NA_character_,
      HPV亚型 == "HPV阴性" ~ "HPV阴性",
      TRUE ~ "HPV阳性"   
    )
  ) %>%
  filter(!is.na(年龄), 年龄 >= 35, 年龄 <= 64, !is.na(HPV组)) %>%
  mutate(HPV组 = factor(HPV组, levels = c("HPV阳性","HPV阴性")))


test_continuous_2group <- function(x, g){
  dat <- tibble(x=x, g=g) %>% filter(!is.na(x), !is.na(g))
  g_levels <- levels(droplevels(dat$g))
  if(length(g_levels) != 2) stop("分组不是2组，无法检验")
  
  x1 <- dat$x[dat$g == g_levels[1]]
  x2 <- dat$x[dat$g == g_levels[2]]
  n1 <- length(x1); n2 <- length(x2)
  
  sh1 <- if(n1 >= 3 && n1 <= 5000) tryCatch(shapiro.test(x1)$p.value, error=function(e) NA_real_) else NA_real_
  sh2 <- if(n2 >= 3 && n2 <= 5000) tryCatch(shapiro.test(x2)$p.value, error=function(e) NA_real_) else NA_real_
  
  use_wilcox <- ( (n1 < 30 || n2 < 30) && ( (!is.na(sh1) && sh1 < 0.05) || (!is.na(sh2) && sh2 < 0.05) ) )
  
  if(use_wilcox){
    p <- wilcox.test(x ~ g, data = dat)$p.value
    method <- "Wilcoxon秩和检验（非参数）"
  } else {
    p <- t.test(x ~ g, data = dat, var.equal = FALSE)$p.value
    method <- "Welch独立样本t检验（方差不齐稳健）"
  }
  list(p_value = p, method = method)
}


test_2x2 <- function(m){
  exp <- tryCatch(chisq.test(m, correct = FALSE)$expected, error=function(e) NULL)
  if(is.null(exp) || any(exp < 5)){
    out <- fisher.test(m)
    list(p_value = out$p.value, method = "Fisher精确检验（2×2）")
  } else {
    out <- chisq.test(m, correct = FALSE)
    list(p_value = out$p.value, method = "Pearson χ²检验（2×2）")
  }
}


test_rx2 <- function(tab){
  exp <- tryCatch(chisq.test(tab, correct = FALSE)$expected, error=function(e) NULL)
  if(is.null(exp) || any(exp < 5)){
    p <- fisher.test(tab)$p.value
    method <- "Fisher精确检验（R×2，期望数不足）"
  } else {
    p <- chisq.test(tab, correct = FALSE)$p.value
    method <- "Pearson χ²检验（R×2）"
  }
  list(p_value=p, method=method)
}


overall_test <- test_continuous_2group(df$年龄, df$HPV组)

overall_sum <- df %>%
  group_by(HPV组) %>%
  summarise(
    n = n(),
    mean = mean(年龄, na.rm = TRUE),
    sd   = sd(年龄, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(stat = sprintf("%.2f ± %.2f", mean, sd))

pos <- overall_sum %>% filter(HPV组 == "HPV阳性")
neg <- overall_sum %>% filter(HPV组 == "HPV阴性")

overall_out <- tibble(
  指标 = "年龄（35-64岁）",
  HPV阳性 = paste0("n=", pos$n, "; ", pos$stat),
  HPV阴性 = paste0("n=", neg$n, "; ", neg$stat),
  P值 = overall_test$p_value,
  统计学方法 = overall_test$method
)

breaks <- c(35,40,45,50,55,60,65)
labels <- c("35-39","40-44","45-49","50-54","55-59","60-64")

df2 <- df %>%
  mutate(
    年龄层 = cut(年龄, breaks = breaks, right = FALSE, include.lowest = TRUE, labels = labels),
    年龄层 = as.character(年龄层)
  ) %>%
  filter(!is.na(年龄层))


tab_rx2 <- table(df2$年龄层, df2$HPV组)
dist_test <- test_rx2(tab_rx2)


denoms <- df2 %>% count(HPV组, name = "组内总数")

age_strata <- df2 %>%
  count(年龄层, HPV组, name = "n") %>%
  left_join(denoms, by="HPV组") %>%
  mutate(
    pct = ifelse(组内总数>0, 100*n/组内总数, NA_real_),
    `n(%)` = sprintf("%d (%.2f%%)", n, pct)
  ) %>%
  select(年龄层, HPV组, `n(%)`)

age_wide <- age_strata %>%
  pivot_wider(names_from = HPV组, values_from = `n(%)`) %>%
  arrange(factor(年龄层, levels = labels))


get_stratum_p <- function(stratum){
  in_stratum <- df2$年龄层 == stratum
  a <- sum(df2$HPV组=="HPV阳性" & in_stratum)
  b <- sum(df2$HPV组=="HPV阴性" & in_stratum)
  c <- sum(df2$HPV组=="HPV阳性" & !in_stratum)
  d <- sum(df2$HPV组=="HPV阴性" & !in_stratum)
  m <- matrix(c(a,b,c,d), nrow=2, byrow=TRUE)
  out <- test_2x2(m)
  tibble(年龄层=stratum, P值=out$p_value, 分层P值方法=out$method)
}

p_each <- map_dfr(labels, get_stratum_p)

age_table_out <- age_wide %>%
  left_join(p_each, by="年龄层") %>%
  mutate(
    分层总体分布P值 = ifelse(年龄层==labels[1], dist_test$p_value, NA_real_),
    分层总体方法 = ifelse(年龄层==labels[1], dist_test$method, NA_character_)
  )

wb <- createWorkbook()

addWorksheet(wb, "Overall_age")
writeData(wb, "Overall_age", overall_out)

addWorksheet(wb, "Age_5yr_strata")
writeData(wb, "Age_5yr_strata", age_table_out)

addWorksheet(wb, "Notes_methods")
notes <- tibble(
  项目 = c("总体年龄比较", "5岁分层总体分布", "每个年龄层单独P值"),
  统计学方法 = c(
    "Welch独立样本t检验；若任一组样本量<30且任一组Shapiro正态性检验P<0.05，则改用Wilcoxon秩和检验",
    "年龄层(6层) × HPV组(2组)：优先Pearson χ²检验；若任一格期望数<5，则用Fisher精确检验",
    "2×2：该年龄层 vs 其他年龄层 × HPV阳性/阴性；优先Pearson χ²；期望数<5用Fisher"
  )
)
writeData(wb, "Notes_methods", notes)

saveWorkbook(wb, file_out, overwrite = TRUE)
message("已输出：", normalizePath(file_out))
