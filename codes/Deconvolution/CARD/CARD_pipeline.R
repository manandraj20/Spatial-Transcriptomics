library(anndata)
libaray(CARD)
import os
import sys
scrna_path = sys.argv[1]//genes*cells
spatial_path = sys.argv[2] //genes*locations
celltype_key = sys.argv[3]
locations_path=sys.argv[4]
output_path = sys.argv[5]
st_anndata<- read_h5ad(scrna_path)
sc_anndata<- read_h5ad(spatial_path)
st_adata<- as.matrix(st_anndata)
sc_adata<-as.matrix(sc_anndata)
cell_type<- read.csv(celltype_key)
cellID<- factor(sc_adata@Dimnames[[2]])
sc_annmeta<- data.frame(cellID = cellID, cell_type= cell_type)
sc_annmeta$sampleInfo=sc_anndata$var$sample
rownames(sc_annmeta)<-sc_adata@Dimnames[[2]]
rownames(spatial_location_anndata)<-st_adata@Dimnames[[2]]
CARD_obj = createCARDObject(
  sc_count = sc_adata,
  sc_meta = sc_annmeta,
  spatial_count = st_adata,
  spatial_location = spatial_location_anndata,
  ct.varname = "cell_type",
  ct.select = unique(sc_annmeta$cell_type),
  sample.varname = "sampleInfo",
  minCountGene = 100,
  minCountSpot = 5)
CARD_obj = CARD_deconvolution(CARD_object = CARD_obj)
write.csv(CARD_obj@Proportion_CARD,file=output_path+"/CARD_result.csv")
p <- CARD.visualize.prop(
  proportion = CARD_obj@Proportion_CARD,        
  spatial_location = CARD_obj@spatial_location, 
  ct.visualize = cell_type,                
  colors = c("lightblue","lightyellow","red"), 
  NumCols = 4)   
print(p)