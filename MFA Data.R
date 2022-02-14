library(readxl)
SNPs <- read_excel("~/Path/To/File", 
                                          sheet = "MFA")


install.packages(c("FactoMineR", "factoextra"))
library("FactoMineR")
library("factoextra")
library(tidyverse)
SNPs2 = SNPs %>% remove_rownames %>% column_to_rownames(var="Strains") %>% as.data.frame()

res.mfa <- MFA(SNPs2, 
               group = c(6, 4), 
               type = c("s", "s"),
               name.group = c( "SNPs","Linear"
                              ),
               graph = FALSE)



fviz_mfa_var(res.mfa, "quanti.var", palette = "jco", 
             col.var.sup = "violet", repel = TRUE)

fviz_mfa_var(res.mfa, "quanti.var", palette = "jco", 
             col.var.sup = "violet", repel = TRUE,
             geom = c("point", "text"), legend = "bottom")

fviz_mfa_var(res.mfa, "quanti.var", col.var = "contrib", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             col.var.sup = "violet", repel = TRUE,
             geom = c("point", "text"))


fviz_mfa_ind(res.mfa, col.ind = "cos2", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE)


fviz_mfa_ind(res.mfa, 
             habillage = "Label", # color by groups 
             palette = c("#00AFBB", "#E7B800", "#FC4E07"),
             addEllipses = TRUE, ellipse.type = "confidence", 
             repel = TRUE # Avoid text overlapping
) 

res.pca.1 <- prcomp(SNPs2, scale = T)

fviz_pca_ind(res.pca.1,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)


fviz_pca_var(res.pca.1,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE,     # Avoid text overlapping
             geom = c("point", "text"))



SNPsPCA <- read_excel("~/Path/To/File", 
                   sheet = "PCA")

SNPsPCA2 = SNPsPCA %>% remove_rownames %>% column_to_rownames(var="Strains") %>% as.data.frame()

res.pca <- prcomp(SNPsPCA2, scale = F)

fviz_pca_ind(res.pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

fviz_pca_var(res.pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)


