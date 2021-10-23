# Set wd to the current directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(dplyr) # for data manipulation
library(ggfortify) # for autoplot function to plot pca
library(ggplot2)
library(ggrepel)
library(gridExtra)  

# Protein concentrations
filtered_lipids = read.table("nfur.ast.merged_filtered_2.csv", header = T, sep = ",")
head(filtered_lipids)

# Stripe the fatty acid and get the class
filtered_lipids$Class = gsub("\\(.*\\)", "", filtered_lipids$FattyAcid)
head(filtered_lipids)

df = filtered_lipids %>% group_by(Class) %>% tally()
colnames(df) = c("Class", "Counts")
df$Percentage = round(((df$Counts*100)/(sum(df$Counts))), digit = 1)
df

bar = ggplot(data=df, aes(x=reorder(Class, -Counts), y = Counts)) +
  geom_bar(stat="identity") + 
  theme_bw() + 
  theme(text = element_text(size=20)) +
  geom_text(aes(label = paste0(Percentage, "%")), vjust = -1)
bar

df <- df %>% 
  arrange(desc(Class)) %>%
  mutate(prop = Counts / sum(df$Counts) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )

pie = ggplot(data=df, aes(x="", y=prop, fill=Class)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + 
  theme(legend.position="none") +
  geom_label_repel(data = df, aes(y = ypos, 
                   label = paste0(Class, " (", Percentage, "%)")), color = "black", size=6, nudge_x = 1, 
                   show.legend = FALSE) +
  scale_fill_brewer(palette="Set3")
pie

pdf(file = "LipidClass_change_contribution.pdf", height = 5, width = 15)
grid.arrange(bar, pie, ncol = 2) 
dev.off()
