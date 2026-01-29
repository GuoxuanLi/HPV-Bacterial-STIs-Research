pkgs <- c("readxl","dplyr","stringr","tidyr","tibble","openxlsx")
ins <- pkgs[!pkgs %in% installed.packages()[,"Package"]]
if(length(ins)>0) install.packages(ins)
invisible(lapply(pkgs, library, character.only = TRUE))

file_in   <- "XX"
sheet_id  <- 1
col_hpv   <- "HPV亚型"           
col_vac   <- "HPV疫苗接种史"     
file_out  <- "HPV疫苗接种史_HPV阳性vs阴性_人数百分比_组间P值.xlsx"


df0 <- read_excel(file_in, sheet = sheet_id)
names(df0) <- str_squish(names(df0))

if(!all(c(col_hpv, col_vac) %in% names(df0))){
  stop("找不到列名：请检查 col_hpv / col_vac 是否与Excel第一行完全一致。当前列名有：\n",
       paste(names(df0), collapse = " | "))
}

df <- df0 %>%
  transmute(
    HPV_raw = .data[[col_hpv]],
    Vac_raw = .data[[col_vac]]
  ) %>%
  mutate(
    
    HPV亚型 = str_squish(as.character(HPV_raw)),
    HPV组 = case_when(
      is.na(HPV亚型)        ~ NA_character_,
      HPV亚型 == ""         ~ NA_character_,
      HPV亚型 == "HPV阴性"  ~ "HPV阴性",
      TRUE                  ~ "HPV阳性"
    ),
    
    
    Vac_str = str_squish(as.character(Vac_raw)),
    疫苗组 = case_when(
      Vac_str == "是"                              ~ "接种",
      Vac_str %in% c("否","未接种","无","从未接种") ~ "未接种",
      TRUE                                         ~ NA_character_
    )
  ) %>%
  filter(!is.na(HPV组), !is.na(疫苗组)) %>%    
  mutate(
    HPV组  = factor(HPV组,  levels = c("HPV阳性","HPV阴性")),
    疫苗组 = factor(疫苗组, levels = c("未接种","接种"))
  )


tab <- table(df$疫苗组, df$HPV组)

chi <- suppressWarnings(chisq.test(tab, correct = FALSE))
exp <- chi$expected

if(any(exp < 5)){
  test_res   <- fisher.test(tab)
  method_all <- "Fisher确切概率法（2×2）"
} else {
  test_res   <- chi
  method_all <- "Pearson χ²检验（2×2）"
}

p_all      <- test_res$p.value
p_all_char <- ifelse(p_all < 0.001, "<0.001", sprintf("%.3f", p_all))


denoms <- df %>%
  count(HPV组, name = "组内总数")

vac_tbl <- df %>%
  count(疫苗组, HPV组, name = "n") %>%
  left_join(denoms, by = "HPV组") %>%
  mutate(
    pct   = 100 * n / 组内总数,
    `n(%)` = sprintf("%d (%.2f%%)", n, pct)
  ) %>%
  select(疫苗组, HPV组, `n(%)`) %>%
  pivot_wider(
    names_from  = HPV组,
    values_from = `n(%)`,
    values_fill = "0 (0.00%)"
  ) %>%
  arrange(疫苗组) %>%
  mutate(
    组间P值     = p_all_char,
    统计学方法 = method_all
  )


vac_tbl$组间P值[-1]     <- NA
vac_tbl$统计学方法[-1] <- NA


wb <- createWorkbook()
addWorksheet(wb, "HPV疫苗接种史")
writeData(wb, "HPV疫苗接种史", vac_tbl)
saveWorkbook(wb, file_out, overwrite = TRUE)

message("已输出：", normalizePath(file_out))
