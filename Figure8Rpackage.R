pkgs <- c("readxl","dplyr","stringr","broom","pROC","ragg")
ins <- pkgs[!pkgs %in% installed.packages()[,"Package"]]
if(length(ins) > 0) install.packages(ins)
invisible(lapply(pkgs, library, character.only = TRUE))


file_path <- "XX"
sheet_id  <- 1

col_hpv  <- "HPV亚型"
col_ct   <- "CT"
col_uu   <- "UU"
col_path <- "病理学"      

label_cin0 <- "CIN0"
label_cin1 <- "CIN1"

out_base <- "Figure_ROC_4models_CIN1"


roc_cols  <- c(A="#4DBBD5", B="#E64B35", C="#3C5488", D="#00A087")
base_family <- "sans" 

CE_XAXIS <- 1      
CE_LAB  <- 1       
CE_MAIN <- 1       
CE_AUC  <- 1       
CE_PANEL<- 1       


AUC_X_FRAC <- 0.35   
AUC_Y_FRAC <- 0.05  


AUC_ADJ_LABEL <- "Unadjusted"
AUC_ADJ_DX <- -10    
AUC_ADJ_DY <- 3    


PANEL_X_FRAC <- -0.13  
PANEL_Y_FRAC <-  1.09  


df <- readxl::read_excel(file_path, sheet = sheet_id)

dat <- df %>%
  dplyr::mutate(
    HPV_raw  = stringr::str_squish(as.character(.data[[col_hpv]])),
    CT_raw   = stringr::str_squish(as.character(.data[[col_ct]])),
    UU_raw   = stringr::str_squish(as.character(.data[[col_uu]])),
    PATH_raw = stringr::str_squish(as.character(.data[[col_path]])),
    
    HPV_pos  = ifelse(is.na(HPV_raw) | HPV_raw=="" | stringr::str_detect(HPV_raw, "阴性|未做"), 0, 1),
    HPV52_58 = ifelse(!is.na(HPV_raw) & stringr::str_detect(HPV_raw, "HPV\\s*52|HPV\\s*58"), 1, 0),
    
    CT_pos   = ifelse(!is.na(CT_raw) & stringr::str_detect(CT_raw, "阳性"), 1, 0),
    UU_pos   = ifelse(!is.na(UU_raw) & stringr::str_detect(UU_raw, "阳性"), 1, 0)
  ) %>%
  dplyr::filter(HPV_pos == 1) %>%
  dplyr::filter(PATH_raw %in% c(label_cin0, label_cin1)) %>%
  dplyr::mutate(Y_cin1 = ifelse(PATH_raw == label_cin1, 1, 0)) %>%
  dplyr::filter(!is.na(Y_cin1), !is.na(CT_pos), !is.na(UU_pos), !is.na(HPV52_58))

cat("==== 分析集检查 ====\n")
cat("n =", nrow(dat), "\n")
cat("CIN1 =", sum(dat$Y_cin1==1), "；CIN0 =", sum(dat$Y_cin1==0), "\n\n")
if(length(unique(dat$Y_cin1)) < 2) stop("Y_cin1 只有一个类别（全0或全1），无法做ROC。")


fit_one <- function(df, formula, model_name){
  fit <- glm(formula, data = df, family = binomial())
  prob <- predict(fit, type = "response")
  
  roc_obj <- pROC::roc(df$Y_cin1, prob, quiet = TRUE, direction = "<")
  auc_obj <- pROC::auc(roc_obj)
  ci_obj  <- pROC::ci.auc(roc_obj)
  
  best <- pROC::coords(roc_obj, x="best", best.method="youden",
                       ret=c("threshold","sensitivity","specificity"),
                       transpose=FALSE)
  
  or_tbl <- broom::tidy(fit, conf.int = TRUE, exponentiate = TRUE) %>%
    dplyr::filter(term != "(Intercept)") %>%
    dplyr::mutate(Model = model_name) %>%
    dplyr::transmute(Model, Term = term, OR = estimate, L95 = conf.low, U95 = conf.high, P = p.value)
  
  auc_tbl <- data.frame(
    Model = model_name,
    N = nrow(df),
    CIN1_n = sum(df$Y_cin1==1),
    CIN0_n = sum(df$Y_cin1==0),
    AUC = as.numeric(auc_obj),
    AUC_L95 = as.numeric(ci_obj[1]),
    AUC_U95 = as.numeric(ci_obj[3]),
    Best_threshold = as.numeric(best["threshold"]),
    Sensitivity = as.numeric(best["sensitivity"]),
    Specificity = as.numeric(best["specificity"])
  )
  
  list(roc=roc_obj, auc=auc_obj, ci=ci_obj, or=or_tbl, auc_tbl=auc_tbl)
}

mA <- fit_one(dat, Y_cin1 ~ UU_pos,            "A: CIN1 ~ UU")
mB <- fit_one(dat, Y_cin1 ~ UU_pos + HPV52_58, "B: CIN1 ~ UU + HPV52/58")
mC <- fit_one(dat, Y_cin1 ~ CT_pos,            "C: CIN1 ~ CT")
mD <- fit_one(dat, Y_cin1 ~ CT_pos + HPV52_58, "D: CIN1 ~ CT + HPV52/58")

OR_all  <- dplyr::bind_rows(mA$or, mB$or, mC$or, mD$or)
AUC_all <- dplyr::bind_rows(mA$auc_tbl, mB$auc_tbl, mC$auc_tbl, mD$auc_tbl)


test_AB <- pROC::roc.test(mA$roc, mB$roc, method = "delong", paired = TRUE)
test_CD <- pROC::roc.test(mC$roc, mD$roc, method = "delong", paired = TRUE)

auc_diff_tbl <- data.frame(
  Comparison = c("A vs B (UU vs UU+HPV52/58)", "C vs D (CT vs CT+HPV52/58)"),
  Z = c(as.numeric(test_AB$statistic), as.numeric(test_CD$statistic)),
  P = c(as.numeric(test_AB$p.value),   as.numeric(test_CD$p.value))
)

write.csv(auc_diff_tbl, paste0(out_base, "_AUC_differences_DeLong.csv"), row.names = FALSE)
write.csv(OR_all,       paste0(out_base, "_OR.csv"),  row.names = FALSE)
write.csv(AUC_all,      paste0(out_base, "_AUC.csv"), row.names = FALSE)


plot_roc_panel <- function(roc_obj, auc_obj, ci_obj,
                           panel="A", main_txt="", roc_col="#4DBBD5"){
  
  plot(roc_obj,
       legacy.axes = TRUE,
       lwd = 2.8,
       col = roc_col,
       main = main_txt,
       xlab = "1 - Specificity",
       ylab = "Sensitivity",
       cex.axis = CE_XAXIS,
       cex.lab  = CE_LAB,
       cex.main = CE_MAIN,
       grid = FALSE)
  
  usr <- par("usr")
  x1 <- usr[1]; x2 <- usr[2]
  y1 <- usr[3]; y2 <- usr[4]
  w  <- x2 - x1
  h  <- y2 - y1
  

  auc_txt <- sprintf("AUC = %.3f (95%% CI %.3f–%.3f)",
                     as.numeric(auc_obj), as.numeric(ci_obj[1]), as.numeric(ci_obj[3]))
  
  auc_x <- x1 + AUC_X_FRAC * w
  auc_y <- y1 + AUC_Y_FRAC * h
  

  line_h <- strheight("A", cex = CE_AUC, units = "user")
  

  text(auc_x + AUC_ADJ_DX * line_h,
       auc_y + AUC_ADJ_DY * line_h,
       labels = AUC_ADJ_LABEL,
       adj = c(0, 0), cex = CE_AUC, font = 1, xpd = NA)
  

  text(auc_x, auc_y, labels = auc_txt,
       adj = c(0, 0), cex = CE_AUC, font = 1, xpd = NA)
  

  panel_x <- x1 + PANEL_X_FRAC * w
  panel_y <- y1 + PANEL_Y_FRAC * h
  
  text(panel_x, panel_y, labels = panel,
       adj = c(0, 1), cex = CE_PANEL, font = 2, xpd = NA)
}


show_4panel <- function(){
  while(!is.null(dev.list())) dev.off()
  
  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)
  
  par(mfrow = c(2,2),
      mar = c(4.2, 4.2, 2.4, 1.1),
      oma = c(0.5, 0.5, 0.5, 0.5),
      family = base_family)
  
  plot_roc_panel(mC$roc, mC$auc, mC$ci, panel="A",
                 main_txt="C. trachomatis predicting CIN1 among HR-HPV positive women",
                 roc_col=roc_cols["C"])
  
  plot_roc_panel(mD$roc, mD$auc, mD$ci, panel="B",
                 main_txt="C. trachomatis + HR-HPV52/58 predicting CIN1 among HR-HPV positive women",
                 roc_col=roc_cols["D"])
  
  plot_roc_panel(mA$roc, mA$auc, mA$ci, panel="C",
                 main_txt="U. urealyticum predicting CIN1 among HR-HPV positive women",
                 roc_col=roc_cols["A"])
  
  plot_roc_panel(mB$roc, mB$auc, mB$ci, panel="D",
                 main_txt="U. urealyticum + HR-HPV52/58 predicting CIN1 among HR-HPV-positive women",
                 roc_col=roc_cols["B"])
}


draw_4panel <- function(){
  

  ragg::agg_tiff(paste0(out_base, ".tiff"),
                 width = 3300, height = 2800, units = "px",
                 res = 300, compression = "lzw")
  par(mfrow = c(2,2),
      mar = c(4.2, 4.2, 2.4, 1.1),
      oma = c(0.5, 0.5, 0.5, 0.5),
      family = base_family)
  
  plot_roc_panel(mC$roc, mC$auc, mC$ci, panel="A",
                 main_txt="C. trachomatis predicting CIN1 among HR-HPV positive women",
                 roc_col=roc_cols["C"])
  plot_roc_panel(mD$roc, mD$auc, mD$ci, panel="B",
                 main_txt="C. trachomatis + HR-HPV52/58 predicting CIN1 among HR-HPV positive women",
                 roc_col=roc_cols["D"])
  plot_roc_panel(mA$roc, mA$auc, mA$ci, panel="C",
                 main_txt="U. urealyticum predicting CIN1 among HR-HPV positive women",
                 roc_col=roc_cols["A"])
  plot_roc_panel(mB$roc, mB$auc, mB$ci, panel="D",
                 main_txt="U. urealyticum + HR-HPV52/58 predicting CIN1 among HR-HPV positive women",
                 roc_col=roc_cols["B"])
  dev.off()
  

  ragg::agg_png(paste0(out_base, ".png"),
                width = 3300, height = 2800, units = "px", res = 200)
  par(mfrow = c(2,2),
      mar = c(4.2, 4.2, 2.4, 1.1),
      oma = c(0.5, 0.5, 0.5, 0.5),
      family = base_family)
  
  plot_roc_panel(mC$roc, mC$auc, mC$ci, panel="A",
                 main_txt="C. trachomatis predicting CIN1 among HR-HPV positive women",
                 roc_col=roc_cols["C"])
  plot_roc_panel(mD$roc, mD$auc, mD$ci, panel="B",
                 main_txt="C. trachomatis + HR-HPV52/58 predicting CIN1 among HR-HPV positive women",
                 roc_col=roc_cols["D"])
  plot_roc_panel(mA$roc, mA$auc, mA$ci, panel="C",
                 main_txt="U. urealyticum predicting CIN1 among HR-HPV positive women",
                 roc_col=roc_cols["A"])
  plot_roc_panel(mB$roc, mB$auc, mB$ci, panel="D",
                 main_txt="U. urealyticum + HR-HPV52/58 predicting CIN1 among HR-HPV positive women",
                 roc_col=roc_cols["B"])
  dev.off()
  
 
  pdf(paste0(out_base, ".pdf"), width = 10, height = 8, family = base_family)
  par(mfrow = c(2,2),
      mar = c(4.2, 4.2, 2.4, 1.1),
      oma = c(0.5, 0.5, 0.5, 0.5),
      family = base_family)
  
  plot_roc_panel(mC$roc, mC$auc, mC$ci, panel="A",
                 main_txt="C. trachomatis predicting CIN1 among HR-HPV positive women",
                 roc_col=roc_cols["C"])
  plot_roc_panel(mD$roc, mD$auc, mD$ci, panel="B",
                 main_txt="C. trachomatis + HR-HPV52/58 predicting CIN1 among HR-HPV positive women",
                 roc_col=roc_cols["D"])
  plot_roc_panel(mA$roc, mA$auc, mA$ci, panel="C",
                 main_txt="U. urealyticum predicting CIN1 among HR-HPV positive women",
                 roc_col=roc_cols["A"])
  plot_roc_panel(mB$roc, mB$auc, mB$ci, panel="D",
                 main_txt="U. urealyticum + HR-HPV52/58 predicting CIN1 among HR-HPV positive women",
                 roc_col=roc_cols["B"])
  dev.off()
  
  cat("\n==== 已输出 ====\n")
  cat(file.path(getwd(), paste0(out_base, ".tiff")), "\n")
  cat(file.path(getwd(), paste0(out_base, ".png")),  "\n")
  cat(file.path(getwd(), paste0(out_base, ".pdf")),  "\n")
  cat(file.path(getwd(), paste0(out_base, "_OR.csv")),  "\n")
  cat(file.path(getwd(), paste0(out_base, "_AUC.csv")), "\n")
  cat(file.path(getwd(), paste0(out_base, "_AUC_differences_DeLong.csv")), "\n")
}


show_4panel()
draw_4panel()
