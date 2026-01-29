library(pwr)


n_case    <- 869   
n_control <- 473   


p1_ct <- 0.047      
p2_ct <- 0.004      


p1_uu <- 0.529     
p2_uu <- 0.326      


p1_ng <- 0.002      
p2_ng <- 0.000      



run_analysis <- function(p1, p2, label) {
  cat(paste0("\n================== 分析指标：", label, " ==================\n"))
  

  if(abs(p1 - p2) < 0.001) {
    cat("  [提示] 两组发生率几乎相同，无法计算效能 (需要无限样本)。\n")
    return(NULL)
  }
  

  h_val <- ES.h(p1, p2)
  cat(sprintf("  两组发生率: Case=%.1f%% vs Control=%.1f%%\n", p1*100, p2*100))
  cat(sprintf("  效应量 (Cohen's h): %.3f\n", h_val))

  apriori <- power.prop.test(p1 = p1, p2 = p2, power = 0.80, sig.level = 0.05)
  n_theory_total <- ceiling(apriori$n) * 2
  
  cat("\n  【A Priori - 理论需求】\n")
  cat(sprintf("  如果要达到 80%% 的效能 (Alpha=0.05):\n"))
  cat(sprintf("  >>> 理论最少总样本量: %d 例 (每组 %d 例)\n", 
              n_theory_total, ceiling(apriori$n)))

  posthoc <- pwr.2p2n.test(h = h_val, n1 = n_case, n2 = n_control, sig.level = 0.05)
  power_val <- posthoc$power
  
  cat("\n  【Post-hoc - 实际效能】\n")
  cat(sprintf("  基于当前样本量 (Case=%d, Control=%d):\n", n_case, n_control))
  if(power_val > 0.999) {
    cat("  >>> 实际统计效能 (Power): > 0.999 (极高)\n")
  } else {
    cat(sprintf("  >>> 实际统计效能 (Power): %.5f\n", power_val))
  }

  cat("\n  【结论验证】\n")
  total_n <- n_case + n_control
  if(total_n > n_theory_total) {
    cat("  √ 实际样本量 > 理论需求量，抽样充足，设计合理。\n")
  } else {
    cat("  ! 注意：实际样本量未达到理论 80% 效能需求。\n")
  }
}



run_analysis(p1_ct, p2_ct, "CT (沙眼衣原体)")


run_analysis(p1_uu, p2_uu, "UU (解脲支原体)")


run_analysis(p1_ng, p2_ng, "NG (淋球菌)")