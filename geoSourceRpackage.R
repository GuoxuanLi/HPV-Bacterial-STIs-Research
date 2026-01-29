pkgs <- c("readxl","dplyr","stringr","tidyr","tibble","openxlsx")
ins <- pkgs[!pkgs %in% installed.packages()[,"Package"]]
if(length(ins)>0) install.packages(ins)
invisible(lapply(pkgs, library, character.only = TRUE))

file_in   <-  "XX"
sheet_id  <- 1
col_hpv   <- "HPV亚型"      
col_city  <- "送检市县"     
file_out  <- "送检市县_HPV阳性vs阴性_人数百分比_总体P值.xlsx"

df0 <- read_excel(file_in, sheet = sheet_id)
names(df0) <- str_squish(names(df0))

if(!all(c(col_hpv, col_city) %in% names(df0))){
  stop("找不到列名：请检查 col_hpv / col_city 是否与Excel第一行完全一致。当前列名有：\n",
       paste(names(df0), collapse = " | "))
}

df <- df0 %>%
  transmute(
    HPV_raw  = .data[[col_hpv]],
    City_raw = .data[[col_city]]
  ) %>%
  mutate(
    HPV亚型   = str_squish(as.character(HPV_raw)),
    送检市县  = str_squish(as.character(City_raw)),
    HPV组 = case_when(
      is.na(HPV亚型)        ~ NA_character_,
      HPV亚型 == ""         ~ NA_character_,
      HPV亚型 == "HPV阴性"  ~ "HPV阴性",
      TRUE                  ~ "HPV阳性"
    )
  ) %>%
  filter(!is.na(HPV组), !is.na(送检市县), 送检市县 != "") %>%
  mutate(
    HPV组 = factor(HPV组, levels = c("HPV阳性","HPV阴性"))
  )

tab <- table(df$送检市县, df$HPV组)

chi_res <- suppressWarnings(chisq.test(tab, correct = FALSE))
exp_tab <- chi_res$expected

if(any(exp_tab < 5)){
  
  if(all(dim(tab) == c(2, 2))){
    test_res   <- fisher.test(tab)
    method_all <- "Fisher确切概率法（2×2）"
  } else {
   
    test_res   <- fisher.test(tab, simulate.p.value = TRUE, B = 1e5)
    method_all <- "Fisher确切概率法（Monte Carlo 模拟，送检市县×HPV组）"
  }
} else {
  test_res   <- chi_res
  method_all <- "Pearson χ²检验（送检市县×HPV组）"
}

p_all      <- test_res$p.value
p_all_char <- ifelse(p_all < 0.001, "<0.001", sprintf("%.3f", p_all))

denoms <- df %>%
  count(HPV组, name = "组内总数")

city_tbl <- df %>%
  count(送检市县, HPV组, name = "n") %>%
  left_join(denoms, by = "HPV组") %>%
  mutate(
    pct   = 100 * n / 组内总数,
    `n(%)` = sprintf("%d (%.2f%%)", n, pct)
  ) %>%
  select(送检市县, HPV组, `n(%)`) %>%
  pivot_wider(
    names_from  = HPV组,
    values_from = `n(%)`,
    values_fill = "0 (0.00%)"
  ) %>%
  arrange(送检市县) %>%
  mutate(
    总体P值     = p_all_char,
    统计学方法 = method_all
  )


city_tbl$总体P值[-1]     <- NA
city_tbl$统计学方法[-1] <- NA


wb <- createWorkbook()
addWorksheet(wb, "送检市县")
writeData(wb, "送检市县", city_tbl)
saveWorkbook(wb, file_out, overwrite = TRUE)

message("已输出：", normalizePath(file_out))
