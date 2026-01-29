library(readxl)
library(dplyr)
library(stringr)
library(openxlsx)


file_path <- "XX"
sheet_id  <- 1
out_xlsx  <- "HPV_亚型感染数量_分布_合并表.xlsx"


df <- read_excel(file_path, sheet = sheet_id, col_names = TRUE)


fix_name <- function(x){
  x <- as.character(x)
  x <- str_replace_all(x, "[\u00A0\u3000]", " ")  
  x <- str_replace_all(x, "[\r\n\t]", " ")
  x <- str_squish(x)
  x
}
names(df) <- fix_name(names(df))

if(!"HPV亚型" %in% names(df)){
  stop("找不到列名：HPV亚型。当前列名为：\n", paste(names(df), collapse = " | "))
}


dat <- df %>%
  mutate(
    hpv_raw = as.character(.data[["HPV亚型"]]),
    
   
    hpv_raw = str_replace_all(hpv_raw, "[，、;；\\s]+", ","),
    
   
    hpv_raw = str_replace_all(hpv_raw, ",{2,}", ","),
    hpv_raw = str_replace_all(hpv_raw, "^,|,$", ""),
    
    
    hpv_raw = ifelse(
      is.na(hpv_raw) | hpv_raw %in% c("", "无", "阴性", "HPV阴性", "HPV 阴性", "NEG", "Negative"),
      NA, hpv_raw
    ),
    
    
    n_types = ifelse(is.na(hpv_raw), 0L, str_count(hpv_raw, ",") + 1L),
    
   
    n_group = case_when(
      n_types == 0 ~ NA_character_,
      n_types == 1 ~ "单一感染",
      n_types == 2 ~ "2种混合感染",
      n_types == 3 ~ "3种混合感染",
      n_types >= 4 ~ "4种及以上混合感染"
    )
  )


pos  <- dat %>% filter(!is.na(n_group))
Npos <- nrow(pos)


add_binom_ci <- function(tab, total_n){
  
  get_ci_str <- function(x, n){
    ci <- binom.test(as.integer(x), as.integer(n))$conf.int
    paste0(round(ci[1]*100, 2), "%–", round(ci[2]*100, 2), "%")
  }
  
  tab %>%
    mutate(
      例数 = as.integer(例数),
      百分率 = round(例数 / total_n * 100, 2),
      `95%CI` = vapply(例数, get_ci_str, character(1), n = total_n)
    )
}




tab_group <- pos %>%
  count(n_group, name = "例数") %>%
  mutate(n_group = factor(n_group, levels = c("单一感染","2种混合感染","3种混合感染","4种及以上混合感染"))) %>%
  arrange(n_group) %>%
  transmute(分类 = "分层（1/2/3/≥4）", 项目 = as.character(n_group), 例数) %>%
  add_binom_ci(total_n = Npos)


tab_exact <- pos %>%
  count(n_types, name = "例数") %>%
  arrange(n_types) %>%
  transmute(分类 = "精细分布（1,2,3,4...）", 项目 = as.character(n_types), 例数) %>%
  add_binom_ci(total_n = Npos)


tab_all <- bind_rows(tab_group, tab_exact)


wb <- createWorkbook()
addWorksheet(wb, "HPV感染数量分布")

writeData(wb, "HPV感染数量分布", tab_all, startRow = 1, startCol = 1)
addStyle(wb, "HPV感染数量分布", createStyle(textDecoration = "bold"),
         rows = 1, cols = 1:ncol(tab_all), gridExpand = TRUE)

note_row <- nrow(tab_all) + 3
writeData(wb, "HPV感染数量分布", paste0("分母：HPV阳性样本数 N = ", Npos),
          startRow = note_row, startCol = 1)
writeData(wb, "HPV感染数量分布",
          "统计学方法：采用描述性统计计算构成比（百分率）；百分率的95%置信区间采用二项分布精确法（Clopper–Pearson）计算（R函数 binom.test）。",
          startRow = note_row + 1, startCol = 1)

saveWorkbook(wb, out_xlsx, overwrite = TRUE)
message("✅ HPV列固定为：HPV亚型；HPV阳性 N = ", Npos, "；已导出：", out_xlsx)
