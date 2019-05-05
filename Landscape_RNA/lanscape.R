library(grid)
library(ComplexHeatmap)



alter_fun_list = list(
  background = function(x, y, w, h) {grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#CCCCCC", col = NA))},
  Nonsense=function(x, y, w, h) {grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#0071C0", col = NA))},
  up=function(x, y, w, h) {grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "red", col = NA))},
  down= function(x, y, w, h) {grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "green", col = NA))},
  In_Frame =function(x, y, w, h) {grid.rect(x, y, w*0.95, h*0.5, gp = gpar(fill = "#92D14F", col = NA))},
  Frame_Shift=function(x, y, w, h) {grid.rect(x, y, w*0.95, h*0.5, gp = gpar(fill = "#00B64E", col = NA))})
  col = c('Nonsense'='#0071C0','up'='red','down'='green','In_Frame'='#92D14F','Frame_Shift'='#00B64E')

mat_1<-read.table('lanscape_sum_result.xls',sep='\t',header=T)
mat_1[is.na(mat_1)] = ""
rownames(mat_1) = mat_1[, 1]
mat_1 = mat_1[, -1]
mat_1 = as.matrix(mat_1)

pdf("DMR_lanscape.pdf")

P1=oncoPrint(mat_1, get_type = function(x) strsplit(x, ",")[[1]],
              alter_fun = alter_fun_list , col= col,
              row_order = NULL,
              pct_gp = gpar(col="white",fontsize = 0.01), row_names_gp = gpar(fontsize = 12),row_names_side = "left",
              column_title = "",column_title_gp=gpar(fontsize=10),show_row_barplot =FALSE ,
              show_column_names = TRUE,show_heatmap_legend=T,
              column_names_gp=gpar(fontsize = 12),
              heatmap_legend_param = list(title = "Alternations",
              at = c("up","down"),
              labels = c("up","down"))
              )

draw(P1)

dev.off()

