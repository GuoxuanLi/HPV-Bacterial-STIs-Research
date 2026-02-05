

pkgs <- c("readxl","dplyr","stringr","tibble","purrr","broom","writexl")
ins <- pkgs[!pkgs %in% installed.packages()[,"Package"]]
if(length(ins)>0) install.packages(ins)
invisible(lapply(pkgs, library, character.only = TRUE))

file_path <- "XX"
sheet_id  <- 1

col_hpv   <- "HPV亚型"
col_path  <- "病理学"
col_ct    <- "CT"             
col_uu    <- "UU"             
col_age   <- "年龄"
col_edu   <- "文化程度"
col_city  <- "送检市县"
col_vax   <- "HPV疫苗接种史"


squish_chr <- function(x){
  x <- as.character(x)
  x <- stringr::str_squish(x)
  x[x==""] <- NA
  x
}


to_posneg <- function(x){
  x <- squish_chr(x)
  out <- dplyr::case_when(
    is.na(x) ~ NA_character_,
    stringr::str_detect(x, "阳性|\\+") ~ "阳性",
    stringr::str_detect(x, "阴性|\\-|未检出|未检测|未见") ~ "阴性",
    TRUE ~ NA_character_
  )
  factor(out, levels = c("阴性","阳性"))
}


to_num <- function(x){
  if (is.numeric(x)) return(as.numeric(x))
  x <- squish_chr(x)
  suppressWarnings(as.numeric(stringr::str_extract(x, "\\d+\\.?\\d*")))
}

make_hpv_risk_group <- function(hpv_str){
  s <- squish_chr(hpv_str)

  has_1618 <- stringr::str_detect(s, "(?<!\\d)(16|18)(?!\\d)")
  has_5258 <- stringr::str_detect(s, "(?<!\\d)(52|58)(?!\\d)")
  
  grp <- dplyr::case_when(
    is.na(s) ~ NA_character_,
    has_1618 ~ "16/18型",
    has_5258 ~ "52/58型",
    TRUE     ~ "其他高危"
  )
  factor(grp, levels = c("其他高危","16/18型","52/58型"))
}

df0 <- readxl::read_excel(file_path, sheet = sheet_id)
names(df0) <- stringr::str_squish(names(df0))

need_cols <- c(col_hpv,col_path,col_ct,col_uu,col_age,col_edu,col_city,col_vax)
miss <- setdiff(need_cols, names(df0))
if(length(miss)>0){
  stop("缺少列名：", paste(miss, collapse = "、"),
       "\n请检查 Excel 第一行列名是否与代码一致。")
}

dat <- df0 %>%
  mutate(
    HPV_raw  = squish_chr(.data[[col_hpv]]),
    PATH_raw = squish_chr(.data[[col_path]]),
    CT       = to_posneg(.data[[col_ct]]),
    UU       = to_posneg(.data[[col_uu]]),
    年龄     = to_num(.data[[col_age]]),
    文化程度 = factor(squish_chr(.data[[col_edu]])),
    送检市县 = factor(squish_chr(.data[[col_city]])),
    HPV疫苗接种史 = factor(squish_chr(.data[[col_vax]])),
    HPV_Risk_Group = make_hpv_risk_group(HPV_raw)
  ) %>%

  filter(!is.na(HPV_raw), HPV_raw != "", !stringr::str_detect(HPV_raw, "阴性")) %>%
 
  filter(is.na(PATH_raw) | PATH_raw != "HPV阴性未做")



fit_one <- function(data, outcome_group, case_level, control_level = "CIN0"){
  d <- data %>%
    mutate(
      Y = dplyr::case_when(
        PATH_raw %in% case_level      ~ 1L,
        PATH_raw == control_level     ~ 0L,
        TRUE                          ~ NA_integer_
      )
    ) %>%
    filter(!is.na(Y)) %>%
 
    filter(!is.na(CT), !is.na(UU), !is.na(年龄),
           !is.na(文化程度), !is.na(送检市县), !is.na(HPV疫苗接种史), !is.na(HPV_Risk_Group))
  
  if(nrow(d) == 0) {
    return(tibble(
      Outcome_Group = outcome_group,
      Term = character(),
      OR = numeric(),
      CI_Lower = numeric(),
      CI_Upper = numeric(),
      P.value = numeric()
    ))
  }
  
  fml <- as.formula("Y ~ CT + UU + 年龄 + 文化程度 + 送检市县 + HPV疫苗接种史 + HPV_Risk_Group")
  m <- glm(fml, data = d, family = binomial())
  

  ci_mat <- tryCatch(
    suppressMessages(confint(m)),
    error = function(e) suppressMessages(confint.default(m))
  )
  ci_df <- as.data.frame(ci_mat) %>%
    tibble::rownames_to_column("term") %>%
    rename(conf.low = `2.5 %`, conf.high = `97.5 %`)
  
  out <- broom::tidy(m) %>%
    left_join(ci_df, by = c("term")) %>%
    mutate(
      Outcome_Group = outcome_group,
      OR       = exp(estimate),
      CI_Lower = exp(conf.low),
      CI_Upper = exp(conf.high),
      P.value  = p.value
    ) %>%

    filter(stringr::str_detect(term, paste0("^", col_ct)) |
             stringr::str_detect(term, paste0("^", col_uu))) %>%
    transmute(
      Outcome_Group,
      Term = term,
      OR, CI_Lower, CI_Upper, P.value
    )
  
  out
}


model_defs <- list(
  list(outcome_group = "CIN1",   case = c("CIN1")),
  list(outcome_group = "CIN2+",  case = c("CIN2+")),
  list(outcome_group = "SCC",    case = c("SCC")),
  list(outcome_group = "Overall",case = c("CIN1","CIN2+","SCC"))
)

res_all <- purrr::map_dfr(
  model_defs,
  ~fit_one(dat, outcome_group = .x$outcome_group, case_level = .x$case, control_level = "CIN0")
)


res_all_fmt <- res_all %>%
  mutate(
    OR       = round(OR, 3),
    CI_Lower = round(CI_Lower, 3),
    CI_Upper = round(CI_Upper, 3),
    P.value  = signif(P.value, 3)
  )


print(res_all_fmt)


out_xlsx <- "Subgroup_Logistic_Pathology_CT_UU.xlsx"
writexl::write_xlsx(list("CT_UU_Subgroup" = res_all_fmt), path = out_xlsx)
message("已导出：", out_xlsx)
