pkgs <- c("readxl","dplyr","stringr","tibble")
ins <- pkgs[!pkgs %in% installed.packages()[,"Package"]]
if(length(ins)>0) install.packages(ins)
invisible(lapply(pkgs, library, character.only = TRUE))


file_path <- "C:/Users/黎先生/Desktop/2025.12.19 初稿（公卫版）/Raw data-20251130.xlsx"
sheet_id  <- 1

df <- readxl::read_excel(file_path, sheet = sheet_id)


age_col <- if ("年龄" %in% names(df)) "年龄" else if ("Age" %in% names(df)) "Age" else stop("找不到年龄列：请命名为 年龄 或 Age")

hpv_col <- dplyr::case_when(
  "HPV亚型" %in% names(df) ~ "HPV亚型",
  "HPV类型" %in% names(df) ~ "HPV类型",
  "HPV_Subtypes" %in% names(df) ~ "HPV_Subtypes",
  TRUE ~ NA_character_
)
if (is.na(hpv_col)) stop("找不到HPV列：请命名为 HPV亚型 / HPV类型 / HPV_Subtypes")

need_cols <- c("CT","UU","NG")
miss <- setdiff(need_cols, names(df))
if (length(miss) > 0) stop(paste("缺少列：", paste(miss, collapse = ", ")))


is_pos <- function(x){
  x <- str_trim(as.character(x))
  str_detect(x, "阳") & !str_detect(x, "阴")
}

df2 <- df %>%
  mutate(
    Age = as.numeric(.data[[age_col]]),
    HPV_raw = str_trim(as.character(.data[[hpv_col]])),
    
   
    HPV_Pos = ifelse(!is.na(HPV_raw) & HPV_raw != "" & !str_detect(HPV_raw, "阴性"), 1L, 0L),
    
    CT_Pos = ifelse(is_pos(CT), 1L, 0L),
    UU_Pos = ifelse(is_pos(UU), 1L, 0L),
    NG_Pos = ifelse(is_pos(NG), 1L, 0L),
    STIs_Pos = ifelse((CT_Pos + UU_Pos + NG_Pos) >= 1, 1L, 0L),
    
    AgeGroup = cut(
      Age,
      breaks = c(35, 40, 45, 50, 55, 60, 65),
      right = FALSE,
      labels = c("35–39","40–44","45–49","50–54","55–59","60–64")
    )
  ) %>%
  filter(!is.na(AgeGroup)) %>%
  mutate(AgeGroup = factor(AgeGroup, levels = c("35–39","40–44","45–49","50–54","55–59","60–64")))


auto_test_2x2 <- function(tab_2x2){
  cs <- suppressWarnings(chisq.test(tab_2x2, correct = FALSE))
  if (min(cs$expected) < 5) {
    ft <- fisher.test(tab_2x2)
    list(p = ft$p.value, method = "Fisher’s exact test")
  } else {
    list(p = cs$p.value, method = "Pearson’s chi-square test")
  }
}


p_vs_ref <- function(data, outcome_col, ref_level = "35–39"){
  lvls <- levels(data$AgeGroup)
  out <- tibble(
    AgeGroup = lvls,
    p = NA_real_,
    method = NA_character_
  )
  out$p[out$AgeGroup == ref_level] <- NA_real_
  out$method[out$AgeGroup == ref_level] <- NA_character_
  
  for (g in lvls) {
    if (g == ref_level) next
    sub <- data %>% filter(AgeGroup %in% c(ref_level, g)) %>%
      mutate(Group2 = ifelse(AgeGroup == g, g, ref_level))
    tab <- table(sub$Group2, sub[[outcome_col]]) 
    
    if (!all(c("0","1") %in% colnames(tab))) {
     
      tab2 <- matrix(0, nrow = 2, ncol = 2,
                     dimnames = list(rownames(tab), c("0","1")))
      tab2[rownames(tab), colnames(tab)] <- tab
      tab <- tab2
    }
    tt <- auto_test_2x2(tab)
    out$p[out$AgeGroup == g] <- tt$p
    out$method[out$AgeGroup == g] <- tt$method
  }
  out
}


summ <- df2 %>%
  group_by(AgeGroup) %>%
  summarise(
    n = n(),
    HPV_pos_n  = sum(HPV_Pos == 1),
    HPV_pos_pct = 100 * HPV_pos_n / n,
    STIs_pos_n = sum(STIs_Pos == 1),
    STIs_pos_pct = 100 * STIs_pos_n / n,
    .groups = "drop"
  ) %>%
  mutate(
    `HPV阳性`  = sprintf("%d (%.2f%%)", HPV_pos_n,  HPV_pos_pct),
    `STIs阳性` = sprintf("%d (%.2f%%)", STIs_pos_n, STIs_pos_pct)
  ) %>%
  select(`年龄组` = AgeGroup, `例数` = n, `HPV阳性`, `STIs阳性`)

p_hpv  <- p_vs_ref(df2, "HPV_Pos",  ref_level = "35–39") %>%
  transmute(`年龄组` = AgeGroup, `P值(HPV vs 35–39)` = ifelse(is.na(p), "", format.pval(p, digits = 3, eps = 0.001)),
            `方法(HPV)` = ifelse(is.na(method), "", method))

p_stis <- p_vs_ref(df2, "STIs_Pos", ref_level = "35–39") %>%
  transmute(`年龄组` = AgeGroup, `P值(STIs vs 35–39)` = ifelse(is.na(p), "", format.pval(p, digits = 3, eps = 0.001)),
            `方法(STIs)` = ifelse(is.na(method), "", method))

table6_ref <- summ %>%
  left_join(p_hpv,  by = "年龄组") %>%
  left_join(p_stis, by = "年龄组")

print(table6_ref)


write.csv(table6_ref, "Table6_vs35_39_Pvalues.csv", row.names = FALSE, fileEncoding = "UTF-8")


cat("\n统计学方法：各年龄组与35–39岁组比较，分类资料采用2×2列联表检验；")
cat("当任一格期望频数<5时采用Fisher精确检验，否则采用Pearson卡方检验。\n")
