library(readxl)
library(dplyr)
library(tidyr)
library(openxlsx)


file_path <-  "XX"
sheet_id  <- 1                                                   

df_raw <- read_excel(file_path, sheet = sheet_id)
names(df_raw) <- trimws(names(df_raw))   


df2 <- df_raw %>%
  mutate(

    HPV亚型 = trimws(as.character(HPV亚型)),
    HPV分组 = case_when(
      HPV亚型 == "HPV阴性"                      ~ "HPV阴性",
      !is.na(HPV亚型) & HPV亚型 != "" & HPV亚型 != "HPV阴性" ~ "HPV阳性",
      TRUE                                      ~ NA_character_
    ),
    HPV分组 = factor(HPV分组, levels = c("HPV阴性","HPV阳性")),
    
    
    文化程度 = trimws(as.character(文化程度)),
    教育分组 = case_when(
      文化程度 %in% c("文盲","小学","初中") ~ "低于高中学历",
      文化程度 %in% c("中专","高中")       ~ "高中学历",
      文化程度 %in% c("大专","大学本科","本科","研究生","研究生及以上") ~ "大学学历",
      TRUE ~ NA_character_
    ),
    教育分组 = factor(教育分组,
                  levels = c("低于高中学历","高中学历","大学学历"))
  )


df_use <- df2 %>%
  filter(!is.na(教育分组), !is.na(HPV分组))


edu3_cnt <- df_use %>%
  count(教育分组, HPV分组) %>%            
  group_by(HPV分组) %>%
  mutate(百分比 = 100 * n / sum(n)) %>%  
  ungroup()

edu3_tab <- edu3_cnt %>%
  mutate(
    列名 = paste0(HPV分组, " 例数(%)"),
    值   = sprintf("%d (%.1f%%)", n, 百分比)
  ) %>%
  select(教育分组, 列名, 值) %>%
  pivot_wider(names_from = 列名, values_from = 值) %>%
  arrange(教育分组)


tbl_chi <- table(df_use$教育分组, df_use$HPV分组)

chi <- suppressWarnings(chisq.test(tbl_chi))

if (any(chi$expected < 5)) {
  if (all(dim(tbl_chi) == c(2, 2))) {
    test_res   <- fisher.test(tbl_chi)
    method_txt <- "Fisher确切概率法"
  } else {
    test_res   <- fisher.test(tbl_chi, simulate.p.value = TRUE, B = 1e5)
    method_txt <- "Fisher确切概率法（Monte Carlo模拟）"
  }
} else {
  test_res   <- chi
  method_txt <- "卡方检验"
}

p_val      <- test_res$p.value
p_val_char <- ifelse(p_val < 0.001, "<0.001", sprintf("%.3f", p_val))

edu3_tab_final <- edu3_tab %>%
  mutate(
    P值              = ifelse(row_number() == 1, p_val_char, ""),
    `统计学方法(HPV)` = ifelse(row_number() == 1, method_txt, "")
  )


wb <- createWorkbook()
addWorksheet(wb, "学历三线表")
writeData(wb, sheet = "学历三线表", x = edu3_tab_final)
saveWorkbook(wb, "学历三档_HPV分组_三线表.xlsx", overwrite = TRUE)
