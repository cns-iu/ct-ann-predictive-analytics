azimuth_organs <- c('PBMC','Brain','Pancreas','Fetal Development','Kidney','Bone Marrow','Lung','Tonsil','Adipose','Heart')

tabula_sapiens_organs <- c('Bladder', 'Blood', 'Bone Marrow', 'Eye', 'Fat', 'Heart', 'Kidney', 'Colon (Large Intestine)', 'Liver', 'Lung', 'Lymph Nodes', 'Mammary Gland', 'Muscle', 'Pancreas', 'Prostate', 'Salivary Gland', 'Skin', 'Intestine (Small Intestine)', 'Spleen', 'Thymus', 'Tongue', 'Trachea', 'Uterus', 'Vasculature')

celltypist_ref_organs <- c('Upper airway','Skin', 'Tonsil','Eye', 'Placenta','Oesophaegus', 'Trachea','Omentum adipose tissue', 'Lung-draining Lymph node','Kidney','Mesenteric Lymph node','Lung','Intestine (Small Intestine)','Decidua','Colon (Large Intestine)','Liver','Spleen','Bone Marrow','Blood','Thymus')

# asctb_organs <- c('Eye','Fallopian Tube','Knee', 'Liver','Lymph Vasculature','Ovary','Pancreas', 'Peripheral Nervous System','Prostate', 'Small Intestine','Intestine','Ureter', 'Urinary Bladder','Uterus', 'Blood','Blood Vasculature','Bone Marrow','Brain','Heart','Kidney', 'Large Intestine','Lung', 'Lymph Node','Skin', 'Spleen','Thymus', 'Bone Marrow Blood','Brain','Heart','Kidney', 'Large Intestine','Lung', 'Lymph Node','Skin', 'Spleen','Thymus', 'Vasculature')


# install.packages("jsonlite")
# install.packages("VennDiagram")

require(VennDiagram)
require(jsonlite)
require(tidyverse)



asctb_organs_rows <- fromJSON("https://hubmapconsortium.github.io/ccf-asct-reporter/assets/sheet-config.json")



ASCTB_VERSION <- 'v1.2'
res <- as.data.frame(asctb_organs_rows[0, c('name','display')])

for (i in 1:nrow(asctb_organs_rows)){
  x <- asctb_organs_rows[i, ]
  version_list <- as.data.frame(x$version)
  
  tryCatch({
    if (grep(pattern=ASCTB_VERSION, x=version_list$hraVersion)){
        res <- as.data.frame(rbind(res, c(x$name, x$display)))
        }
  },
  error=function(e){
    cat(paste('\nDidn\'t find either `hraVersion` as a field for `',x$display,'` or something went wrong...'))
    print(e)
  })
}

colnames(res) <- c('name', 'display')
res <- res %>%
  filter(!grepl('All', display)) %>%
  select(display)


res[res$display=='Small Intestine', 'display'] <- 'Intestine (Small Intestine)'
res[res$display=='Large Intestine', 'display'] <- 'Colon (Large Intestine)'

asctb_organs <- res$display


# original asctb vs celltypist_ref_organs only
# v <- venn.diagram(
#   x=list(celltypist_ref_organs, asctb_organs),
#   category.names=c('CellTypist','ASCTB'),
#   fill=c('yellow', 'light green'),
#   alpha=c(0.5, 0.5),
#   cat.cex=1.5, cex=1,
#   filename=NULL,
#   disable.logging=T,
#   main=paste0('CellTypist reference-organs vs ASCTB ', ASCTB_VERSION, ' organs', '\nOrgan-name in brackets is the ASCTB name.')
# )

# v[[5]]$label  <- paste(v[[5]]$label, '\n', paste0(setdiff(asctb_organs, celltypist_ref_organs), collapse="\n"))
# # in asctb_organs only
# v[[6]]$label <- paste(v[[6]]$label, '\n', paste0(setdiff(celltypist_ref_organs, asctb_organs)  , collapse="\n"))
# # intesection
# v[[7]]$label <- paste(v[[7]]$label, '\n', paste0(intersect(celltypist_ref_organs, asctb_organs), collapse="\n"))




v <- venn.diagram(
  x=list(celltypist_ref_organs, tabula_sapiens_organs, azimuth_organs),
  category.names=c('CellTypist','Tabula Sapiens', 'Azimuth'),
  fill=c('yellow', 'green', 'light blue'),
  alpha=c(0.5, 0.5, 0.5),
  cat.cex=1.5, 
  cex=1,
  filename=NULL,
  disable.logging=T,
  main=paste0('Reference-organ comparison for CellTypist vs Azimuth vs Tabula Sapiens', '\nOrgan-name in brackets is a manually assigned common name.')
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
v <- update_labels_3_way_intersection(v, A=azimuth_organs, B=celltypist_ref_organs, C=tabula_sapiens_organs)


# plot  
grid.newpage()
grid.draw(v)
