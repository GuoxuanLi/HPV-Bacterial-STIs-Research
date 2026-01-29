
pkgs <- c("readxl","dplyr","stringr","tidyr","tibble","openxlsx")
ins <- pkgs[!pkgs %in% installed.packages()[,"Package"]]
if(length(ins)>0) install.packages(ins)
invisible(lapply(pkgs, library, character.only = TRUE))

file_in   <-  "XX"
sheet_id  <- 1
col_hpv   <- "HPV亚型"   
col_eth   <- "民族"
file_out  <- "民族_HPV阳性vs阴性_人数百分比_总体P值.xlsx"


df0 <- read_excel(file_in, sheet = sheet_id)
names(df0) <- str_squish(names(df0))

if(!all(c(col_hpv, col_eth) %in% names(df0))){
  stop("找不到列名：请检查 col_hpv / col_eth 是否与Excel第一行完全一致。当前列名有：\n",
       paste(names(df0), collapse = " | "))
}

df <- df0 %>%
  transmute(
    HPV_raw = .data[[col_hpv]],
    民族_raw = .data[[col_eth]]
  ) %>%
  mutate(
    HPV亚型 = str_squish(as.character(HPV_raw)),
    民族    = str_squish(as.character(民族_raw)),
    HPV组 = case_when(
      is.na(HPV亚型)        ~ NA_character_,
      HPV亚型 == ""         ~ NA_character_,
      HPV亚型 == "HPV阴性"  ~ "HPV阴性",
      TRUE                  ~ "HPV阳性"
    )
  ) %>%
  filter(!is.na(HPV组), !is.na(民族), 民族 != "") %>%
  mutate(HPV组 = factor(HPV组, levels = c("HPV阳性","HPV阴性")))


df_minority <- df %>% filter(民族 != "汉族")

if(nrow(df_minority) == 0){
  stop("数据中没有少数民族，无法计算少数民族的组间P值。")
}

tab_min <- table(df_minority$民族, df_minority$HPV组)

exp_min <- tryCatch(chisq.test(tab_min, correct = FALSE)$expected,
                    error = function(e) NULL)

if(is.null(exp_min) || any(exp_min < 5)){
  p_all     <- fisher.test(tab_min)$p.value
  method_all <- "Fisher精确检验（少数民族×HPV组）"
} else {
  chi_res   <- chisq.test(tab_min, correct = FALSE)
  p_all     <- chi_res$p.value
  method_all <- "Pearson χ²检验（少数民族×HPV组）"
}

p_all_char <- ifelse(p_all < 0.001, "<0.001", sprintf("%.3f", p_all))


denoms <- df %>%
  count(HPV组, name = "组内总数")   

eth_tbl <- df %>%
  count(民族, HPV组, name = "n") %>%
  left_join(denoms, by = "HPV组") %>%
  mutate(
    pct   = 100 * n / 组内总数,
    `n(%)` = sprintf("%d (%.2f%%)", n, pct)
  ) %>%
  select(民族, HPV组, `n(%)`) %>%
  pivot_wider(
    names_from  = HPV组,
    values_from = `n(%)`,
    values_fill = "0 (0.00%)"
  ) %>%
  arrange(民族) %>%
  mutate(
    总体P值     = p_all_char,
    统计学方法 = method_all
  )


eth_tbl$总体P值[-1]     <- NA
eth_tbl$统计学方法[-1] <- NA


wb <- createWorkbook()
addWorksheet(wb, "Ethnicity")
writeData(wb, "Ethnicity", eth_tbl)
saveWorkbook(wb, file_out, overwrite = TRUE)

message("已输出：", normalizePath(file_out))
