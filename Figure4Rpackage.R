pkgs <- c("tidyverse","readxl","circlize")
ins <- pkgs[!pkgs %in% installed.packages()[,"Package"]]
if(length(ins) > 0) install.packages(ins)
invisible(lapply(pkgs, library, character.only = TRUE))

excel_file <- "XX"
sheet_id   <- 1

out_prefix_main <- "Fig4_cooccurrence_count"
out_prefix_sens <- "FigSx_cooccurrence_jaccard"

label_cex <- 2
jaccard_display_scale <- 100
jaccard_min_cutoff <- 0.00

df_raw <- read_excel(excel_file, sheet = sheet_id)


names(df_raw) <- stringr::str_squish(names(df_raw))

pick_idx <- function(nms, patterns){
  hit <- which(stringr::str_detect(nms, stringr::regex(paste(patterns, collapse="|"), ignore_case = TRUE)))
  if(length(hit) == 0){
    stop(
      paste0("❌ 未识别到列：", paste(patterns, collapse="/"),
             "\n当前列名：", paste(nms, collapse="，")),
      call. = FALSE
    )
  }
  if(length(hit) > 1) message("⚠️ 多个候选列：", paste(nms[hit], collapse="；"), "；默认使用：", nms[hit[1]])
  hit[1]
}

uu_idx  <- pick_idx(names(df_raw), c("^UU$", "UU", "解脲", "Ureaplasma"))
hpv_idx <- pick_idx(names(df_raw), c("^HPV亚型$", "^HPV_Subtypes$", "HPV.*(亚型|分型|subtype|genotype)", "(亚型|分型|subtype|genotype).*HPV", "^HPV$"))
ct_idx  <- pick_idx(names(df_raw), c("^CT$", "CT", "衣原体", "Chlamydia"))
ng_idx  <- pick_idx(names(df_raw), c("^NG$", "NG", "淋球菌", "Gonorr"))

target_idx <- c(uu_idx, hpv_idx, ct_idx, ng_idx)
target_nm  <- c("UU","HPV_Subtypes","CT","NG")


for(nm in target_nm){
  pos <- which(names(df_raw) == nm)
  if(length(pos) > 0){
    pos_other <- setdiff(pos, target_idx)
    if(length(pos_other) > 0){
      names(df_raw)[pos_other] <- paste0(nm, "_old")
    }
  }
}


names(df_raw)[uu_idx]  <- "UU"
names(df_raw)[hpv_idx] <- "HPV_Subtypes"
names(df_raw)[ct_idx]  <- "CT"
names(df_raw)[ng_idx]  <- "NG"


df_raw <- df_raw %>% mutate(ID = row_number())

df_hpv <- df_raw %>%
  select(ID, HPV_Subtypes) %>%
  mutate(HPV_Subtypes = as.character(HPV_Subtypes)) %>%
  separate_rows(HPV_Subtypes, sep = "[,，、\\s]+") %>%
  mutate(Pathogen = str_trim(HPV_Subtypes)) %>%
  filter(!is.na(Pathogen), Pathogen != "", !str_detect(Pathogen, "阴性")) %>%
  mutate(Pathogen = str_replace(Pathogen, "^HPV", "HR-HPV")) %>%
  select(ID, Pathogen)

df_uu <- df_raw %>% filter(str_detect(UU, "阳|\\+")) %>% mutate(Pathogen = "UU") %>% select(ID, Pathogen)
df_ct <- df_raw %>% filter(str_detect(CT, "阳|\\+")) %>% mutate(Pathogen = "CT") %>% select(ID, Pathogen)
df_ng <- df_raw %>% filter(str_detect(NG, "阳|\\+")) %>% mutate(Pathogen = "NG") %>% select(ID, Pathogen)

df_all <- bind_rows(df_hpv, df_uu, df_ct, df_ng) %>%
  filter(!is.na(ID), !is.na(Pathogen), Pathogen != "") %>%
  distinct(ID, Pathogen)

if(nrow(df_all) == 0) stop("数据提取为空：请检查 Excel 内容/阳性标记/分隔符。")

tab_id_path <- table(df_all$ID, df_all$Pathogen)
M <- (tab_id_path > 0) * 1

A_count <- t(M) %*% M
diag(A_count) <- 0

n_i <- colSums(M)

A_jaccard <- matrix(0,
                    nrow = ncol(M), ncol = ncol(M),
                    dimnames = list(colnames(M), colnames(M))
)

for(i in seq_len(ncol(M))){
  for(j in seq_len(ncol(M))){
    inter <- A_count[i, j]
    uni   <- n_i[i] + n_i[j] - inter
    A_jaccard[i, j] <- ifelse(uni > 0, inter / uni, 0)
  }
}
diag(A_jaccard) <- 0

if(jaccard_min_cutoff > 0){
  A_jaccard[A_jaccard < jaccard_min_cutoff] <- 0
}

all_pathogens <- colnames(M)

grid_col <- setNames(rep("#4DBBD5", length(all_pathogens)), all_pathogens)
if ("UU" %in% all_pathogens) grid_col["UU"] <- "#E64B35"
if ("CT" %in% all_pathogens) grid_col["CT"] <- "#00A087"
if ("NG" %in% all_pathogens) grid_col["NG"] <- "#3C5488"

label_map <- c(
  "UU" = "U. urealyticum",
  "CT" = "C. trachomatis",
  "NG" = "N. gonorrhoeae"
)

dist_hpv <- 1.00
dist_special_map <- c(CT = 1.42, UU = 1.48, NG = 1.54)
x_shift_map      <- c(CT = 0.00, UU = 0.00, NG = 0.00)

draw_chord <- function(mat, title_text = NULL){
  
  circos.clear()
  circos.par(
    start.degree = 90,
    gap.degree = 1.5,
    track.margin = c(0.01, 0.01),
    cell.padding = c(0, 0, 0, 0),
    points.overflow.warning = FALSE
  )
  
  chordDiagram(
    mat,
    grid.col = grid_col,
    transparency = 0.4,
    annotationTrack = "grid",
    annotationTrackHeight = 0.06,
    preAllocateTracks = list(track.height = 0.30)
  )
  
  circos.track(track.index = 1, panel.fun = function(x, y) {
    
    xlim <- get.cell.meta.data("xlim")
    ylim <- get.cell.meta.data("ylim")
    sector.name <- get.cell.meta.data("sector.index")
    
    dist    <- dist_hpv
    x_shift <- 0
    
    if(sector.name %in% names(dist_special_map)) dist <- dist_special_map[sector.name]
    if(sector.name %in% names(x_shift_map))      x_shift <- x_shift_map[sector.name]
    
    final_y <- ylim[1] + dist
    x_pos   <- mean(xlim) + x_shift * diff(xlim)
    
    label_to_show <- ifelse(sector.name %in% names(label_map), label_map[sector.name], sector.name)
    
    circos.text(
      x_pos, final_y, label_to_show,
      facing = "reverse.clockwise",
      niceFacing = TRUE,
      adj = c(0, 0.5),
      col = "black",
      font = 2,
      cex = label_cex
    )
    
  }, bg.border = NA)
  
  if(!is.null(title_text)){
    mtext(title_text, side = 3, line = -1.5, cex = 1.0, font = 2)
  }
  
  circos.clear()
}

while (!is.null(dev.list())) dev.off()

tiff(paste0(out_prefix_main, ".tiff"),
     width = 12, height = 12, units = "in", res = 300, compression = "lzw")
draw_chord(A_count, title_text = "Co-occurrence network (absolute counts)")
dev.off()

pdf(paste0(out_prefix_main, ".pdf"),
    width = 12, height = 12, useDingbats = FALSE)
draw_chord(A_count, title_text = "Co-occurrence network (absolute counts)")
dev.off()

A_jaccard_display <- A_jaccard * jaccard_display_scale

tiff(paste0(out_prefix_sens, ".tiff"),
     width = 12, height = 12, units = "in", res = 300, compression = "lzw")
draw_chord(A_jaccard_display)
dev.off()

pdf(paste0(out_prefix_sens, ".pdf"),
    width = 12, height = 12, useDingbats = FALSE)
draw_chord(A_jaccard_display)
dev.off()

message("🎉 已生成：",
        out_prefix_main, ".tiff / .pdf  +  ",
        out_prefix_sens, ".tiff / .pdf")
