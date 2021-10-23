# Set wd to the current directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

library(ComplexHeatmap) # For Heatmap function
library(circlize) # for generating colors using colorremp2
library(gplots) # for heatmap.2 which has better row names
library(dplyr)
library(data.table) # For %like% function

# file = "nfur.D1m_v_nfur.PreD.Y.txt"
# col.prefix = "D1m_PreD"
ReadLipidFile = function(file, col.prefix, columns) {
  
  df = read.table(file=file, head = T, sep = "\t")
  df$FattyAcid = paste0(df$class, df$FattyAcid) 
  head(df)
  colnames(df)
  df = df[, columns]
  head(df)
  colnames(df) = c("FattyAcid", paste0(col.prefix, ".padj"), paste0(col.prefix, ".FC"))
  head(df)
  return(df)
}

# Read lipid DE files ------------------------------------------------------------------------------------------------------------------
# nfur.D1m_v_nfur.PreD.Y
D1m_PreD = ReadLipidFile("nfur.D1m_v_nfur.PreD.Y.txt", "D1m_PreD", c(3, 13, 14))
head(D1m_PreD)

# nfur.D1m_v_NonD
D1m_NonD = ReadLipidFile("nfur.D1m_v_NonD.txt", "D1m_NonD", c(3, 12, 13)) # For samples with NonD exclude one column
head(D1m_NonD)

# nfur.D1m_v_nfur.Dexit
D1m_Dexit = ReadLipidFile("nfur.D1m_v_nfur.Dexit.txt", "D1m_Dexit", c(3, 13, 14)) # For samples with NonD exclude one column
head(D1m_Dexit)

# nfur.D1m_v_nfur.D6d
D1m_D6d = ReadLipidFile("nfur.D1m_v_nfur.D6d.txt", "D1m_D6d", c(3, 13, 14)) # For samples with NonD exclude one column
head(D1m_D6d)

# nfur.D6D_v_nfur.PreD.Y
D6d_PreD = ReadLipidFile("nfur.D6D_v_nfur.PreD.Y.txt", "D6d_PreD", c(3, 13, 14))
head(D6d_PreD)

# nfur.D6d_v_NonD
D6d_NonD = ReadLipidFile("nfur.D6d_v_NonD.txt", "D6d_NonD", c(3, 12, 13)) # For samples with NonD exclude one column
head(D6d_NonD)

# nfur.D6d_v_nfur.Dexit
D6d_Dexit = ReadLipidFile("nfur.D6d_v_nfur.Dexit.txt", "D6d_Dexit", c(3, 13, 14)) # For samples with NonD exclude one column
head(D6d_Dexit)

# NonD_v_nfur.PreD.Y
NonD_PreD = ReadLipidFile("NonD_v_nfur.PreD.Y.txt", "NonD_PreD", c(3, 12, 13)) # For samples with NonD exclude one column
head(NonD_PreD)

# nonD_v_nfur.Dexit
NonD_Dexit = ReadLipidFile("nonD_v_nfur.Dexit.txt", "NonD_Dexit", c(3, 12, 13)) # For samples with NonD exclude one column
head(NonD_Dexit)

# nfur.PreD.Y_v_nfur.Dexit
PreDY_Dexit = ReadLipidFile("nfur.PreD.Y_v_nfur.Dexit.txt", "PreDY_Dexit", c(3, 13, 14))
head(PreDY_Dexit)

# nfur.PreD.Y_v_Ast.PreD.Y
nfur_ast = ReadLipidFile("nfur.PreD.Y_v_Ast.PreD.Y.txt", "Nfur_Ast", c(3, 13, 14))
head(nfur_ast)

# nfur.PreD.Y_v_Nfur.PreD.O
nfur_PreD_YO = ReadLipidFile("nfur.PreD.Y_v_nfur.PreD.O.txt", "NfurY_NfurO", c(3, 13, 14))
head(nfur_PreD_YO)

merged <- Reduce(function(x, y) merge(x, y, by="FattyAcid", all = TRUE), 
                 list(D1m_PreD, D1m_NonD, D1m_Dexit, D1m_D6d, 
                      D6d_PreD, D6d_NonD, D6d_Dexit, 
                      NonD_PreD, NonD_Dexit, 
                      PreDY_Dexit, nfur_PreD_YO,
                      nfur_ast))
rownames(merged) = merged$FattyAcid
head(merged)

# take mean of the FC
merged$D6d.mean.FC = rowMeans(merged[,c('D6d_PreD.FC', 'D6d_NonD.FC')])
merged$D1m.mean.FC = rowMeans(merged[,c('D1m_PreD.FC', 'D1m_NonD.FC')])

write.csv(merged, "nfur.ast.merged_All_2.csv")

filtered_all = merged[((merged$D6d_PreD.padj <= 0.05 |
                   merged$D6d_NonD.padj <= 0.05 |
                   merged$D1m_PreD.padj <= 0.05 |
                   merged$D1m_NonD.padj <= 0.05) &
                   (merged$NonD_PreD.padj > 0.05)) ,]
dim(filtered_all)
filtered_all <- filtered_all[!is.infinite(rowSums(filtered_all[, 2:25])),] # Remove non infinited rows created by dvision by 0
dim(filtered_all)
colnames(filtered_all)

write.csv(filtered_all, "nfur.ast.merged_filtered_2.csv")

# Heatmap with raw values - ALL ---------------------------------------
colnames(filtered_all)

All_DE_heatmap = Heatmap(as.matrix(filtered_all[, c(3, 5, 9, 11, 13, 17)]), # Nfur                         
                         heatmap_legend_param = list(title = "Fold change"),                 
                         cluster_rows = T, 
                         cluster_columns = T,
                         show_row_names = T,
                         row_names_gp = gpar(fontsize = 6),
                         row_names_max_width = unit(30, "cm")
                         )
All_DE_heatmap

Ast_DE_heatmap = Heatmap(as.matrix(filtered_all[, c(25, 25)]), # Ast plotted in the same order as in above
        heatmap_legend_param = list(title = "Fold change"),
        cluster_rows = F, 
        show_row_names = T,
        row_order = unlist(row_order(All_DE_heatmap)),
        row_names_gp = gpar(fontsize = 6),
        row_names_max_width = unit(30, "cm")
)
All_combined = All_DE_heatmap + Ast_DE_heatmap

pdf (file = "Nfur_All_Combined_FC.pdf", width = 4.5, height = 25)
draw(All_combined)
dev.off()

# colnames(filtered_all)
# 
# # With average values ALL ---------------------------------------------------------------------
# Heatmap(as.matrix(filtered_all[, c(11, 14, 15)]),
#         cluster_rows = T,
#         cluster_columns = T,
#         show_row_names = T,
#         row_names_gp = gpar(fontsize = 6),
#         row_names_max_width = unit(30, "cm")
# )

# -------------------------------------------------------------------------------------
# # Divide them into two parts: TG and non-TG
head(filtered_all)
colnames(filtered_all)
filtered_TG = filtered_all[filtered_all$FattyAcid %like% "TG",]
filtered_TG

TG_heatmap = Heatmap(as.matrix(filtered_TG[, c(3, 5, 9, 11, 13, 17)]), # Nfur                     
                     heatmap_legend_param = list(title = "Fold change"),
                     cluster_rows = T,
                     cluster_columns = T,
                     show_row_names = T,
                     row_names_gp = gpar(fontsize = 6),
                     row_names_max_width = unit(30, "cm")
)
TG_heatmap

TG_heatmap_Ast = Heatmap(as.matrix(filtered_TG[, c(25, 25)]), # Ast plotted in the same order as in above
        heatmap_legend_param = list(title = "Fold change"),
        colorRamp2(c(-10, 0, 10), c("blue", "white", "red")),
        cluster_rows = F,
        show_row_names = T,
        row_order = unlist(row_order(TG_heatmap)),
        row_names_gp = gpar(fontsize = 6),
        row_names_max_width = unit(30, "cm")
)
TG_heatmap_Ast

TG_combined = TG_heatmap + TG_heatmap_Ast

pdf (file = "Nfur-Ast_TG_Combined_SM.pdf", width = 4.5, height = 10)
draw(TG_combined)
dev.off()
