
library(tidyverse)
library(readxl)
library(ComplexHeatmap)
library(grid)

excel_file <- "XX"
sheet_id   <- 1
df_raw <- read_excel(excel_file, sheet = sheet_id)


names(df_raw) <- names(df_raw) %>%
  stringr::str_replace_all("\ufeff", "") %>%                                  
  stringr::str_replace_all("[\u00A0\u2000-\u200B\u202F\u205F\u3000]", " ") %>%  
  stringr::str_replace_all("[\r\n\t]", " ") %>%                                
  stringr::str_squish()

pick_col <- function(nms, patterns) {
  hit <- which(stringr::str_detect(nms, stringr::regex(paste(patterns, collapse="|"), ignore_case = TRUE)))
  if (length(hit) == 0) {
    stop(
      paste0("❌ 未识别到列。当前列名：", paste(nms, collapse = "，"),
             "\n匹配关键字：", paste(patterns, collapse = "/")),
      call. = FALSE
    )
  }
  if (length(hit) > 1) message("⚠️ 匹配到多个候选列：", paste(nms[hit], collapse="；"), "；默认使用：", nms[hit[1]])
  nms[hit[1]]
}

uu_col  <- pick_col(names(df_raw), c("^UU$", "UU", "解脲", "Ureaplasma"))
hpv_col <- pick_col(names(df_raw), c("^HPV亚型$", "^HPV_Subtypes$", "HPV.*(亚型|分型|subtype|genotype)", "(亚型|分型|subtype|genotype).*HPV", "^HPV$"))

names(df_raw)[names(df_raw) == uu_col]  <- "UU"
names(df_raw)[names(df_raw) == hpv_col] <- "HPV_Subtypes"

message("✅ UU列识别为：", uu_col, " -> 重命名为 UU")
message("✅ HPV列识别为：", hpv_col, " -> 重命名为 HPV_Subtypes")
stopifnot(all(c("UU","HPV_Subtypes") %in% names(df_raw)))


plot_data <- df_raw %>%
  mutate(ID = row_number()) %>%
  mutate(UU = as.character(UU), HPV_Subtypes = as.character(HPV_Subtypes)) %>%
  filter(!is.na(HPV_Subtypes)) %>%
  mutate(HPV_Subtypes = str_squish(HPV_Subtypes)) %>%
  filter(HPV_Subtypes != "") %>%
  filter(!str_detect(HPV_Subtypes, "阴性")) %>%  
  mutate(
    Combined_Infection = ifelse(str_detect(UU, "阳") | str_detect(UU, "\\+"),
                                paste0(HPV_Subtypes, ", UU"),
                                HPV_Subtypes)
  ) %>%
  separate_rows(Combined_Infection, sep = "[,，、\\s]+") %>%
  mutate(Pathogen = str_trim(Combined_Infection)) %>%
  filter(Pathogen != "") %>%
  filter(!str_detect(Pathogen, "阴性")) %>%
  filter(!str_detect(Pathogen, "\\.\\.")) %>%
  select(ID, Pathogen) %>%
  mutate(Value = 1) %>%
  distinct() %>%
  pivot_wider(names_from = Pathogen, values_from = Value, values_fill = 0) %>%
  select(-ID) %>%
  as.data.frame()


colnames(plot_data) <- ifelse(grepl("^HPV", colnames(plot_data)),
                              sub("^HPV", "HR-HPV", colnames(plot_data)),
                              colnames(plot_data))

colnames(plot_data) <- ifelse(colnames(plot_data) == "UU",
                              "U. urealyticum",
                              colnames(plot_data))


m_full <- make_comb_mat(as.matrix(plot_data))
keep_n <- min(40, length(comb_size(m_full)))
m <- m_full[order(comb_size(m_full), decreasing = TRUE)[1:keep_n]]


bright_palette <- c(
  "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99",
  "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A",
  "#FFFF99", "#B15928", "#8DD3C7", "#FFFFB3", "#BEBADA",
  "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9"
)

final_set_order <- order(set_size(m), decreasing = TRUE)
ordered_names <- set_name(m)[final_set_order]

n_colors <- length(ordered_names)
bg_base_colors <- rep(bright_palette, length.out = n_colors)
row_bg_colors <- paste0(bg_base_colors, "33")
names(row_bg_colors) <- ordered_names

all_sets <- set_name(m)
right_bar_colors <- rep("#2C3E50", length(all_sets))
highlight_list <- c("U. urealyticum", "HR-HPV52", "HR-HPV58")
right_bar_colors[all_sets %in% highlight_list] <- "#C0392B"

n_top_bars <- length(comb_size(m))
top_bar_colors <- rep("#2C3E50", n_top_bars)
top_bar_colors[1:min(5, n_top_bars)] <- "#C0392B"


LEFT_SETNAME_FONTSIZE  <- 20
TOP_AXIS_FONTSIZE      <- 20
TOP_NUM_FONTSIZE       <- 20
RIGHT_AXIS_FONTSIZE    <- 20
RIGHT_NUM_FONTSIZE     <- 20
INTERSECTION_TITLE_FS  <- 20


right_anno <- upset_right_annotation(
  m,
  gp = gpar(fill = right_bar_colors),
  add_numbers = FALSE,
  bar_width = 0.5,
  annotation_name_side = "bottom",
  annotation_name_gp = gpar(fontsize = 20, fontface = "bold", col = "black"),
  axis_param = list(side = "bottom",
                    gp = gpar(fontsize = RIGHT_AXIS_FONTSIZE, fontface = "bold"))
)

ht <- UpSet(
  m,
  name = "main_matrix",
  set_order = final_set_order,
  row_names_gp = gpar(col = "black", fontsize = LEFT_SETNAME_FONTSIZE, fontface = "bold"),
  right_annotation = right_anno,
  top_annotation = upset_top_annotation(
    m,
    add_numbers = FALSE,
    show_annotation_name = FALSE,
    gp = gpar(fill = top_bar_colors),
    height = unit(7, "cm"),
    axis_param = list(gp = gpar(fontsize = TOP_AXIS_FONTSIZE, fontface = "bold"))
  ),
  pt_size = unit(3.5, "mm"),
  lwd = 2,
  comb_order = order(comb_size(m), decreasing = TRUE)
)


draw_and_decorate <- function(ht_obj, m_obj, ordered_names, row_bg_colors) {
  
  draw(ht_obj, padding = unit(c(2, -15, 10, 12), "mm"))
  
  decorate_heatmap_body("main_matrix", {
    n_rows <- length(ordered_names)
    for (i in seq_len(n_rows)) {
      current_name <- ordered_names[i]
      bg_col <- row_bg_colors[current_name]
      grid.rect(
        x = 0.5, y = 1 - (i - 0.5) / n_rows,
        width = 1, height = 1 / n_rows,
        gp = gpar(fill = bg_col, col = NA)
      )
    }
  })
  
  decorate_annotation("intersection_size", {
    
    grid.text("Intersection Size",
              x = unit(-25, "mm"),
              rot = 90,
              just = "bottom",
              gp = gpar(fontsize = INTERSECTION_TITLE_FS, fontface = "bold"))
    
    vals <- comb_size(m_obj)
    TOP_NUM_OFFSET_MM <- 1.2
    
    grid.text(
      label = vals,
      x = seq_len(length(vals)),
      y = unit(vals, "native") + unit(TOP_NUM_OFFSET_MM, "mm"),
      default.units = "native",
      just = c("left", "bottom"),
      rot = 45,
      gp = gpar(fontsize = TOP_NUM_FONTSIZE, col = "black", fontface = "bold")
    )
  })
  
  decorate_annotation("set_size", {
    vals_set <- set_size(m_obj)[ordered_names]
    n <- length(vals_set)
    y_pos <- rev(seq_len(n))
    RIGHT_NUM_OFFSET_MM <- 1.2
    
    grid.text(
      label = vals_set,
      x = unit(vals_set, "native") + unit(RIGHT_NUM_OFFSET_MM, "mm"),
      y = y_pos,
      default.units = "native",
      just = c("left", "center"),
      gp = gpar(fontsize = RIGHT_NUM_FONTSIZE, fontface = "bold", col = "black")
    )
  })
}


draw_and_decorate(ht, m, ordered_names, row_bg_colors)
