require(tidyverse)
require(dplyr)
require(reshape)

celltypist_lung_preds <- read.csv('Datasets/LungMap/CellTypist_annotations_from_LungMap.csv')
# celltypist <- celltypist_lung_preds$predicted_labels %>% unique() # 58
celltypist <- celltypist_lung_preds$predicted_labels


azimuth_lung_preds <- read.csv('Datasets/LungMap/azimuth_pred.tsv', sep='\t')
# azimuth <- azimuth_lung_preds$predicted.ann_finest_level %>% unique() # 45
azimuth <- azimuth_lung_preds$predicted.ann_finest_level

setdiff(celltypist, azimuth)


results <- matrix(data=NA, nrow=length(azimuth), ncol=length(celltypist))
for (i in seq_along(azimuth)) {
  for (j in seq_along(celltypist)) {
    print(paste0(azimuth[i], ' & ', celltypist[j], ' --> ', ifelse(azimuth[i]==celltypist[j], 1, 0)))
    results[i,j] <- ifelse(azimuth[i]==celltypist[j], 1, 0)
  }
}


res <- azimuth==celltypist
barplot(table(res), main='Matching predictions for LungMap dataset', sub='CellTypist vs Azimuth', col=c('red','green'))



colnames(results) <- celltypist
results <- results %>%
  as.data.frame() %>%
  mutate(rnames=factor(azimuth, levels=(azimuth)[order(azimuth)]))



results %>%
  select(order(colnames(results))) %>%
  melt() %>%
  ggplot(aes(x=rnames, y=variable, fill=as.factor(value))) +
  geom_tile(color='white', lwd=.1, linetype=1) +
  scale_fill_manual(name='CT-label match', values=c('grey','darkgreen')) +
  labs(
    title='LungMap predictions from Azimuth and CellTypist', 
    x=paste0('Azimuth (',length(azimuth), ' CTs)'), 
    y=paste0('CellTypist (',length(celltypist), ' CTs)')
    ) +
  coord_fixed() +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))




# Azimuth doesnt have Mesothelium, CellTypist has Mesothelium.











































########################## VENN DIAGRAM FOR COMPARING PREDICTION LABELS FROM 2 TOOLS
require(VennDiagram)
require(jsonlite)
require(tidyverse)



azimuth_preds <- read.csv('Datasets/LungMap/azimuth_pred.tsv', sep='\t')
celltypist_preds <- read.csv('Datasets/LungMap/CellTypist_annotations_from_LungMap.csv')



celltypist_preds <- celltypist_preds %>%
  select(majority_voting) %>%
  unique() %>%
  as.vector()
celltypist_preds <- tolower(celltypist_preds$majority_voting)


azimuth_preds <- azimuth_preds %>%
  select(predicted.ann_finest_level) %>%
  unique() %>%
  as.vector()
azimuth_preds <- tolower(azimuth_preds$predicted.ann_finest_level)


popv_preds <- c('type ii pneumocyte','neutrophil','cd4-positive alpha-beta t cell','cd8-positive alpha-beta t cell','nk cell','bronchial vessel endothelial cell','smooth muscle cell','adventitial cell','macrophage','respiratory mucous cell','basophil','classical monocyte','endothelial cell of artery','lung ciliated cell','dendritic cell','cd4-positive, alpha-beta t cell','respiratory goblet cell','basal cell','plasma cell','b cell','cd8-positive, alpha-beta t cell','serous cell of epithelium of bronchus','non-classical monocyte','capillary endothelial cell','capillary aerocyte','type i pneumocyte','vein endothelial cell','alveolar fibroblast','endothelial cell of lymphatic vessel','bronchial smooth muscle cell','plasmacytoid dendritic cell','club cell','lung microvascular endothelial cell','vascular associated smooth muscle cell','myofibroblast cell','pericyte cell','mesothelial cell','fibroblast','intermediate monocyte','pulmonary ionocyte')



# TESTING OUT ONCLASS STANDARDIZATIONS
celltypist_preds <- c('fibroblast of the aortic adventitia', 'fibroblast of areolar connective tissue', 'alveolar macrophage', 'cerebellar stellate cell', 'omentum preadipocyte', 'omentum preadipocyte', 'type i pneumocyte', 'type ii pneumocyte', 'columnar neuron', 'b cell', 'tip cell', 'cd4-positive type i nk t cell', 'r8 photoreceptor cell', 'classical monocyte', 'cap cell', 'non-branched duct epithelial cell', 'conventional dendritic cell', 'plasmacytoid dendritic cell', 'omentum preadipocyte', 'flight muscle cell', 'arcuate artery cell', 'cap cell', 'pulmonary artery endothelial cell', 'arcuate artery cell', 'fibroblast of villous mesenchyme', 'lung goblet cell', 'nasal mucosa goblet cell', 'conjunctiva goblet cell', 'omentum preadipocyte', 'ionocyte', 'lymph gland plasmatocyte', 'lymphangioblast', 'lymph gland plasmatocyte', 'mast cell', 'perirenal preadipocyte', 'migratory trunk neural crest cell', 'lamellocyte', 'polymodal neuron', 'unimodal nocireceptor', 'cystoblast', 'neuroepidermoblast', 'natural killer cell', 'non-branched duct epithelial cell', 'fibroblast of the aortic adventitia', 'perirenal preadipocyte', 'plasma cell', 'plasmatocyte', 'stuff accumulating cell', 'urethra cell', 'mucus secreting cell', 'fibroblast of dermis', 'obsolete pstab/alc', 'smooth muscle cell of bladder', 'fibroblast of the aortic adventitia', 'omentum preadipocyte', 'keratin accumulating cell', 'type ii pinealocyte', 'lamellocyte')
azimuth_preds <- c('fibroblast of the aortic adventitia', 'fibroblast of areolar connective tissue', 'alveolar macrophage', 'omentum preadipocyte', 'type i pneumocyte', 'type ii pneumocyte', 'columnar neuron', 'b cell', 'tip cell', 'cd4-positive type i nk t cell', 'r8 photoreceptor cell', 'classical monocyte', 'cap cell', 'non-branched duct epithelial cell', 'plasmacytoid dendritic cell', 'omentum preadipocyte', 'flight muscle cell', 'arcuate artery cell', 'cap cell', 'pulmonary artery endothelial cell', 'arcuate artery cell', 'nasal mucosa goblet cell', 'omentum preadipocyte', 'ionocyte', 'lymph gland plasmatocyte', 'lymphangioblast', 'mast cell', 'lamellocyte', 'polymodal neuron', 'unimodal nocireceptor', 'cystoblast', 'natural killer cell', 'non-branched duct epithelial cell', 'fibroblast of the aortic adventitia', 'perirenal preadipocyte', 'plasma cell', 'plasmatocyte', 'stuff accumulating cell', 'mucus secreting cell', 'smooth muscle cell of bladder', 'omentum preadipocyte', 'keratin accumulating cell', 'type ii pinealocyte', 'lamellocyte')
popv_preds <- c('adventitial cell', 'fibroblast of areolar connective tissue', 'b cell', 'basal cell', 'basophil', 'bronchial smooth muscle cell', 'peritubular capillary endothelial cell', 'plasmatocyte', 'capillary endothelial cell', 'cd4-negative cd8-negative gamma-delta intraepithelial t cell', 'cd4-negative cd8-negative gamma-delta intraepithelial t cell', 'cd8-positive, alpha-beta t cell', 'cd8-positive, alpha-beta t cell', 'classical monocyte', 'club cell', 'dendritic cell', 'endothelial cell of artery', 'endothelial cell of lymphatic vessel', 'fibroblast', 'intermediate monocyte', 'lung ciliated cell', 'lung microvascular endothelial cell', 'macrophage', 'mesothelial cell', 'myofibroblast cell', 'neutrophil', 'natural killer cell', 'non-branched duct epithelial cell', 'pericyte', 'plasma cell', 'plasmacytoid dendritic cell', 'fibroblast of pulmonary artery', 'respiratory goblet cell', 'lung secretory cell', 'serous cell of epithelium of bronchus', 'smooth muscle cell', 'type i pneumocyte', 'type ii pneumocyte', 'vascular associated smooth muscle cell', 'vein endothelial cell')



# v <- venn.diagram(
#   x=list(celltypist_ref_organs, asctb_organs),
#   category.names=c('CellTypist','PopV'),
#   fill=c('yellow', 'light green'),
#   alpha=c(0.5, 0.5),
#   cat.cex=1.5, cex=1,
#   filename=NULL,
#   disable.logging=T,
#   main=paste0('CellTypist predictions vs PopV')
# )


# # in celltypist_ref_organs only
# v[[5]]$label  <- paste(v[[5]]$label, '\n', paste0(sort(setdiff(asctb_organs, celltypist_ref_organs)), collapse="\n"))
# # in asctb_organs only
# v[[6]]$label <- paste(v[[6]]$label, '\n', paste0(sort(setdiff(celltypist_ref_organs, asctb_organs)), collapse="\n"))
# # intesection
# v[[7]]$label <- paste(v[[7]]$label, '\n', paste0(sort(intersect(celltypist_ref_organs, asctb_organs)), collapse="\n"))




v <- venn.diagram(
  x=list(azimuth_preds, celltypist_preds, popv_preds),
  category.names=c('Azimuth.finest','CellTypist','PopV'),
  fill=c('light blue', 'yellow', 'light green'),
  alpha=c(.5, .5, .5),
  cat.cex=1.5, cex=2,
  filename=NULL,
  disable.logging=T,
  main=paste0('Predictions for Azimuth.finest vs CellTypist vs PopV')
)




update_labels_3_way_intersection <- function(v, A, B, C){
  A_intersect_B_intersect_C <- intersect(A, intersect(B, C))
  B_intersect_C <- intersect(B, C)
  A_intersect_B <- intersect(B, A)
  C_intersect_A <- intersect(A, C)
  
  
  label7 <- setdiff(B, union(B_intersect_C, union(A_intersect_B, A_intersect_B_intersect_C)))
  v[[7]]$label <- paste(v[[7]]$label, '\n', paste0(sort(label7), collapse="\n"))
  
  label8 <- setdiff(B_intersect_C, A_intersect_B_intersect_C)
  v[[8]]$label <- paste(v[[8]]$label, '\n', paste0(sort(label8), collapse="\n"))
  
  label9 <- setdiff(C, union(B_intersect_C, setdiff(C_intersect_A, A_intersect_B_intersect_C)))
  v[[9]]$label <- paste(v[[9]]$label, '\n', paste0(sort(label9), collapse="\n"))
  
  label10 <- setdiff(A_intersect_B, A_intersect_B_intersect_C)
  v[[10]]$label <- paste(v[[10]]$label, '\n', paste0(sort(label10), collapse="\n"))
  
  label11 <- A_intersect_B_intersect_C
  v[[11]]$label <- paste(v[[11]]$label, '\n', paste0(sort(label11), collapse="\n"))
  
  label12 <- setdiff(C_intersect_A, A_intersect_B_intersect_C)
  v[[12]]$label <- paste(v[[12]]$label, '\n', paste0(sort(label12), collapse="\n"))
  
  label13 <- setdiff(A, union(C_intersect_A, setdiff(A_intersect_B, A_intersect_B_intersect_C)))
  v[[13]]$label <- paste(v[[13]]$label, '\n', paste0(sort(label13), collapse="\n"))
  
  return (v)
}





# Azimuth     : A
# CellTypist  : B
# ASCTB/TS    : C
# v <- update_labels_3_way_intersection(v, B=celltypist_preds, A=azimuth_preds, C=popv_preds)
v <- update_labels_3_way_intersection(v, B=c(), A=c(), C=c())


# plot  
grid.newpage()
grid.draw(v)