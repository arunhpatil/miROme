setwd("D:/Halushka_lab/Arun/microRNAome_manuscript/Figures")
data <- read.table("forScatterPlot2.txt", header = T, sep="\t")
#x = log10(data$Filtered_miRNA_Reads)
x = data$Filtered_miRNA_Reads
y = data$Unique_miRNAs
plot(x, y, pch = 19, col = factor(data$groups))
group =  data$groups
legend("topleft",
       legend = levels(factor(group)),
       pch = 19,
       col = factor(levels(factor(group))))


library(ggplot2)
ggplot2.scatterplot(data=data, xName='Filtered miRNA reads',yName='Unique miRNAs', 
                    groupName='Class', size=3,
                    backgroundColor="white",
                    groupColors=c('#999999','#E69F00', '#56B4E9'))  

df=as.data.frame(data)
ggplot(df,aes(x=log10(Filtered_miRNA_Reads),y=Unique_miRNAs, col=groups, shape = unique(groups)))+geom_point()

ggplot(df,aes(x=log10(Filtered_miRNA_Reads),y=Unique_miRNAs)) + geom_point(aes(shape=groups,color=groups)) + scale_shape_manual(name="shape", values=c(0,15,1, 16, 2, 17, 10, 11, 12,13,14, 25, 23)) + scale_color_manual(name="colour", values=c("red","blue","yellow", "green", "purple", "black","Pink", "Maroon", "Brown", "Lavender", "darkgreen", "orange", "skyblue"))
ggplot(df,aes(x=log10(Filtered_miRNA_Reads),y=Unique_miRNAs)) + geom_point(aes(shape=groups,color=groups), size=2) + scale_shape_manual(name="shape", values=c(0,15,1, 16, 2, 17, 10, 11, 12,13,14, 25, 23)) + scale_color_manual(name="colour", values=c("red","blue","yellow", "green", "purple", "black","Pink", "Maroon", "Brown", "Lavender", "#469990", "orange", "#aaffc3"))
ggplot(df,aes(x=log10(Filtered_miRNA_Reads),y=Unique_miRNAs)) + geom_point(aes(color=groups), size=1.5) + scale_color_manual(name="colour", values=c("red","blue","yellow", "green", "purple", "black","Pink", "Maroon", "Brown", "Lavender", "#469990", "orange", "#aaffc3"))
ggplot(df,aes(x=Filtered_miRNA_Reads,y=Unique_miRNAs, col=groups))+geom_point() + scale_x_continuous(trans='log10')

ggplot(df,aes(x=Filtered_miRNA_Reads,y=Unique_miRNAs)) + geom_point(aes(shape=groups,color=groups), size=2) + scale_shape_manual(name="shape", values=c(0,15,1, 16, 2, 17, 10, 11, 12,13,14, 25, 23)) + scale_color_manual(name="colour", values=c("red","blue","yellow", "green", "purple", "black","Pink", "Maroon", "Brown", "Lavender", "#469990", "orange", "#aaffc3")) + scale_x_continuous(trans='log10')
ggplot(df,aes(x=Filtered_miRNA_Reads,y=Unique_miRNAs)) + geom_point(aes(color=groups), size=1.5) + scale_color_manual(name="colour", values=c("red","blue","yellow", "green", "purple", "black","Pink", "Maroon", "Brown", "Lavender", "#469990", "orange", "#aaffc3")) + scale_x_continuous(trans='log10')


ggplot(df,aes(x=Filtered_miRNA_Reads,y=Unique_miRNAs)) + geom_point(aes(color=groups), size=1.5) + scale_color_manual(name="colour", values=c("#F8766D", "#E18A00", "#BE9C00", "#8CAB00", "#24B700", "#00BE70", "#00C1AB", "#00BBDA", "#00ACFC", "#8B93FF", "#D575FE", "#F962DD", "#FF65AC")) + scale_x_continuous(trans='log10')
ggplot(df,aes(x=Filtered_miRNA_Reads,y=Unique_miRNAs)) + geom_point(aes(color=groups), size=1.5) + scale_color_manual(name="colour", values=c("#F8766D", "#E18A00", "#BE9C00", "#8CAB00", "#24B700", "#00BE70", "#00C1AB", "#00BBDA", "#00ACFC", "#8B93FF", "#D575FE", "#F962DD", "#FF65AC")) + scale_x_continuous(trans='log10') + theme_bw()

# command used for scatter plot
ggplot(df,aes(x=Filtered_miRNA_Reads,y=Unique_miRNAs)) + geom_point(aes(color=groups), size=1.5) + scale_color_manual(name="colour", values=c("#F8766D", "#E18A00", "#BE9C00", "#8CAB00", "#24B700", "#00BE70", "#00C1AB", "#00BBDA", "#00ACFC", "#8B93FF", "#D575FE", "#F962DD", "#FF65AC")) + scale_x_continuous(trans='log10') + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

library(scales)
show_col(hue_pal()(13))


## BOX PLOT
# ggplot(df, aes(y=groups, x=Unique_miRNAs)) + 
#   geom_boxplot(notch = TRUE)+
#   geom_jitter(position=position_jitter(0.2))

nsx <- list(counts = as.numeric(df$Unique_miRNAs), group = as.factor(df$groups))
nsx <- as.tibble(nsx)
q<- ggplot(nsx, aes(group, counts, fill = group)) + geom_boxplot() 
q + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#q + theme(axis.text.x = element_text(size=6, angle = 90, vjust = 0.5, hjust=1)) + coord_flip()
#q + theme(axis.text.x = element_text(size=6, angle = 90, vjust = 0.5, hjust=1))  + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
q + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),axis.text.x = element_text(size=12, angle = 90, vjust = 0.5, hjust=1) )

