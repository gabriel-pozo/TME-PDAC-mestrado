#Pre-process Peng

metadata <- read.csv("D:/Scanpy/peng_to_seurat/metadata.csv")

rownames(metadata) <- metadata$X

metadata$X <- NULL

counts <- readMM("D:/Scanpy/peng_to_seurat/matrix.mtx")
genes <- read_tsv("D:/Scanpy/peng_to_seurat/features.tsv", col_names = FALSE)
cell_ids <- read_tsv("D:/Scanpy/peng_to_seurat/barcodes.tsv", col_names = FALSE)

rownames(counts) <- genes$X1
colnames(counts) <- cell_ids$X1

peng <- CreateSeuratObject(counts = counts, meta.data = metadata)

View(peng@meta.data)

peng[["percent.mt"]] <- PercentageFeatureSet(peng, pattern = "^MT-")
peng[["percent.hb"]] <- PercentageFeatureSet(peng, pattern = "^HBA|^HBB")
peng[["percent.rp"]] <- PercentageFeatureSet(peng, pattern = "^RPS|^RPL")
peng$author <- "Peng_2019"
metadata <- peng@meta.data
metadata <- metadata %>% select(-c(total_counts, total_counts_mt, n_genes_by_counts, pct_counts_mt, orig.ident))
peng@meta.data <- metadata

nFeature_lower <- 500
nFeature_upper <- 10000
nCount_lower <- 1000
nCount_upper <- 100000
percent.mt_lower <- 0
percent.mt_upper <- 30
percent.hb_lower <- 0
percent.hb_upper <- 5

#Filter low quality cells
peng <- subset(peng, subset = nFeature_RNA > nFeature_lower & nFeature_RNA < nFeature_upper & nCount_RNA > nCount_lower & nCount_RNA < nCount_upper & percent.mt < percent.mt_upper & percent.hb < percent.hb_upper)

peng$barcode <- colnames(peng)
write.csv(peng@meta.data, file='D:Scanpy/Peng_2019/metadata.csv', quote=F, row.names=F)

counts_matrix <- GetAssayData(peng, assay='RNA', slot='counts')
writeMM(counts_matrix, file='D:Scanpy/Peng_2019/matrix.mtx')

write.table(
  data.frame('gene'=rownames(counts_matrix)),file='D:Scanpy/Peng_2019/genes.csv',
  quote=F,row.names=F,col.names=F
)

#Pre-process Chen

T25 <- Read10X(data.dir = "C:/Users/Gabriel/Downloads/GSE212966_RAW/PDAC1/")
T25 <- CreateSeuratObject(counts = T25, project = 'T25', min.cells = 3, min.features = 200)

T26 <- Read10X(data.dir = "C:/Users/Gabriel/Downloads/GSE212966_RAW/PDAC2/")
T26 <- CreateSeuratObject(counts = T26, project = 'T26', min.cells = 3, min.features = 200)

T27 <- Read10X(data.dir = "C:/Users/Gabriel/Downloads/GSE212966_RAW/PDAC3/")
T27 <- CreateSeuratObject(counts = T27, project = 'T27', min.cells = 3, min.features = 200)

T28 <- Read10X(data.dir = "C:/Users/Gabriel/Downloads/GSE212966_RAW/PDAC4/")
T28 <- CreateSeuratObject(counts = T28, project = 'T28', min.cells = 3, min.features = 200)

T29 <- Read10X(data.dir = "C:/Users/Gabriel/Downloads/GSE212966_RAW/T35/")
T29 <- CreateSeuratObject(counts = T29, project = 'T29', min.cells = 3, min.features = 200)

T30 <- Read10X(data.dir = "C:/Users/Gabriel/Downloads/GSE212966_RAW/PDAC6/")
T30 <- CreateSeuratObject(counts = T30, project = 'T30', min.cells = 3, min.features = 200)

ADJ4 <- Read10X(data.dir = "C:/Users/Gabriel/Downloads/GSE212966_RAW/ADJ1/")
ADJ4 <- CreateSeuratObject(counts = ADJ4, project = 'ADJ4', min.cells = 3, min.features = 200)

ADJ5 <- Read10X(data.dir = "C:/Users/Gabriel/Downloads/GSE212966_RAW/ADJ2/")
ADJ5 <- CreateSeuratObject(counts = ADJ5, project = 'ADJ5', min.cells = 3, min.features = 200)


ADJ6 <- Read10X(data.dir = "C:/Users/Gabriel/Downloads/GSE212966_RAW/ADJ6/")
ADJ6 <- CreateSeuratObject(counts = ADJ6, project = 'ADJ6', min.cells = 3, min.features = 200)

T25$Condition <- "T"
T25$Patient <- "T25"
T26$Condition <- "T"
T26$Patient <- "T26"
T27$Condition <- "T"
T27$Patient <- "T27"
T28$Condition <- "T"
T28$Patient <- "T28"
T29$Condition <- "T"
T29$Patient <- "T29"
T30$Condition <- "T"
T30$Patient <- "T30"
ADJ4$Condition <- "ADJ"
ADJ4$Patient <- "ADJ4"
ADJ5$Condition <- "ADJ"
ADJ5$Patient <- "ADJ5"
ADJ6$Condition <- "ADJ"
ADJ6$Patient <- "ADJ6"

obj <- merge(T25, y = c(T26,T27,T28,T29,T30,ADJ4,ADJ5,ADJ6), add.cell.ids = c("T25","T26","T27","T28","T29","T30","ADJ4","ADJ5","ADJ6"),
             project = "GSE212966")

#Mito/Hemoglobin/Ribosomal ratio
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
obj[["percent.hb"]] <- PercentageFeatureSet(obj, pattern = "^HBA|^HBB")
obj[["percent.rp"]] <- PercentageFeatureSet(obj, pattern = "^RPS|^RPL")
obj$author <- "Chen_2022"

View(obj@meta.data)

obj@meta.data$orig.ident <- NULL

#Filter low quality cells
obj <- subset(obj, subset = nFeature_RNA > nFeature_lower & nFeature_RNA < nFeature_upper & nCount_RNA > nCount_lower & nCount_RNA < nCount_upper & percent.mt < percent.mt_upper & percent.hb < percent.hb_upper)

obj$barcode <- colnames(obj)
write.csv(obj@meta.data, file='D:Scanpy/Chen_2022/metadata.csv', quote=F, row.names=F)

counts_matrix <- GetAssayData(obj, assay='RNA', slot='counts')
writeMM(counts_matrix, file='D:Scanpy/Chen_2022/matrix.mtx')

write.table(
  data.frame('gene'=rownames(counts_matrix)),file='D:Scanpy/Chen_2022/genes.csv',
  quote=F,row.names=F,col.names=F
)


#Pre-process Steele

ADJ1 <- Read10X(data.dir = "D:scRNA - Human PDAC/AdjNorm_TISSUE_1")
ADJ1 <- CreateSeuratObject(counts = ADJ1, project = 'ADJ1', min.cells = 3, min.features = 200)

ADJ2 <- Read10X(data.dir = "D:scRNA - Human PDAC/AdjNorm_TISSUE_2")
ADJ2 <- CreateSeuratObject(counts = ADJ2, project = 'ADJ2', min.cells = 3, min.features = 200)

ADJ3 <- Read10X(data.dir = "D:scRNA - Human PDAC/AdjNorm_TISSUE_3")
ADJ3 <- CreateSeuratObject(counts = ADJ3, project = 'ADJ3', min.cells = 3, min.features = 200)

T31 <- Read10X(data.dir = "D:scRNA - Human PDAC/PDAC_1")
T31 <- CreateSeuratObject(counts = T31, project = 'T31', min.cells = 3, min.features = 200)

T32 <- Read10X(data.dir = "D:scRNA - Human PDAC/PDAC_2")
T32 <- CreateSeuratObject(counts = T32, project = 'T32', min.cells = 3, min.features = 200)

T33 <- Read10X(data.dir = "D:scRNA - Human PDAC/PDAC_3")
T33 <- CreateSeuratObject(counts = T33, project = 'T33', min.cells = 3, min.features = 200)

T34 <- Read10X(data.dir = "D:scRNA - Human PDAC/PDAC_4")
T34 <- CreateSeuratObject(counts = T34, project = 'T34', min.cells = 3, min.features = 200)

T35 <- Read10X(data.dir = "D:scRNA - Human PDAC/PDAC_5")
T35 <- CreateSeuratObject(counts = T35, project = 'T35', min.cells = 3, min.features = 200)

T36 <- Read10X(data.dir = "D:scRNA - Human PDAC/PDAC_6")
T36 <- CreateSeuratObject(counts = T36, project = 'T36', min.cells = 3, min.features = 200)

T37 <- Read10X(data.dir = "D:scRNA - Human PDAC/PDAC_7")
T37 <- CreateSeuratObject(counts = T37, project = 'T37', min.cells = 3, min.features = 200)

T38 <- Read10X(data.dir = "D:scRNA - Human PDAC/PDAC_8")
T38 <- CreateSeuratObject(counts = T38, project = 'T38', min.cells = 3, min.features = 200)

T39 <- Read10X(data.dir = "D:scRNA - Human PDAC/PDAC_9")
T39 <- CreateSeuratObject(counts = T39, project = 'T39', min.cells = 3, min.features = 200)

T40 <- Read10X(data.dir = "D:scRNA - Human PDAC/PDAC_10")
T40 <- CreateSeuratObject(counts = T40, project = 'T40', min.cells = 3, min.features = 200)

T41 <- Read10X(data.dir = "D:scRNA - Human PDAC/PDAC_11A")
T41 <- CreateSeuratObject(counts = T41, project = 'T41', min.cells = 3, min.features = 200)

T42 <- Read10X(data.dir = "D:scRNA - Human PDAC/PDAC_11B")
T42 <- CreateSeuratObject(counts = T42, project = 'T42', min.cells = 3, min.features = 200)

T43 <- Read10X(data.dir = "D:scRNA - Human PDAC/PDAC_12")
T43 <- CreateSeuratObject(counts = T43, project = 'T43', min.cells = 3, min.features = 200)

T44 <- Read10X(data.dir = "D:scRNA - Human PDAC/PDAC_13")
T44 <- CreateSeuratObject(counts = T44, project = 'T44', min.cells = 3, min.features = 200)

T45 <- Read10X_h5(filename = "D:scRNA - Human PDAC/PDAC_14/filtered_feature_bc_matrix.h5")
T45 <- CreateSeuratObject(counts = T45, project = "T45", min.cells = 3, min.features = 200)

T46 <- Read10X(data.dir = "D:scRNA - Human PDAC/PDAC_15")
T46 <- CreateSeuratObject(counts = T46, project = 'T46', min.cells = 3, min.features = 200)

T47 <- Read10X(data.dir = "D:scRNA - Human PDAC/PDAC_16")
T47 <- CreateSeuratObject(counts = T47, project = 'T47', min.cells = 3, min.features = 200)

T31$Condition <- "T"
T31$Patient <- "T31"
T32$Condition <- "T"
T32$Patient <- "T32"
T33$Condition <- "T"
T33$Patient <- "T33"
T34$Condition <- "T"
T34$Patient <- "T34"
T35$Condition <- "T"
T35$Patient <- "T35"
T36$Condition <- "T"
T36$Patient <- "T36"
T37$Condition <- "T"
T37$Patient <- "T37"
T38$Condition <- "T"
T38$Patient <- "T38"
T39$Condition <- "T"
T39$Patient <- "T39"
T40$Condition <- "T"
T40$Patient <- "T40"
T41$Condition <- "T"
T41$Patient <- "T41"
T42$Condition <- "T"
T42$Patient <- "T42"
T43$Condition <- "T"
T43$Patient <- "T43"
T44$Condition <- "T"
T44$Patient <- "T44"
T45$Condition <- "T"
T45$Patient <- "T45"
T46$Condition <- "T"
T46$Patient <- "T46"
T47$Condition <- "T"
T47$Patient <- "T47"

ADJ1$Condition <- "ADJ"
ADJ1$Patient <- "ADJ1"
ADJ2$Condition <- "ADJ"
ADJ2$Patient <- "ADJ2"
ADJ3$Condition <- "ADJ"
ADJ3$Patient <- "ADJ3"


obj <- merge(T31, y = c(T32,T33,T34,T35,T36,T37,T38,T39,T40,T41,T42,T43,T44,T45,T46,T47,ADJ1,ADJ2,ADJ3), add.cell.ids = c("T31","T32","T33","T34","T35","T36","T37","T38","T39",
                                                                                                                          "T40","T41","T42","T43","T44","T45","T46","T47","ADJ1",
                                                                                                                          "ADJ2","ADJ3"),
             project = "Steele")

#Mito/Hemoglobin/Ribosomal ratio
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
obj[["percent.hb"]] <- PercentageFeatureSet(obj, pattern = "^HBA|^HBB")
obj[["percent.rp"]] <- PercentageFeatureSet(obj, pattern = "^RPS|^RPL")
obj$author <- "Steele_2020"

View(obj@meta.data)

obj@meta.data$orig.ident <- NULL

#Filter low quality cells
obj <- subset(obj, subset = nFeature_RNA > nFeature_lower & nFeature_RNA < nFeature_upper & nCount_RNA > nCount_lower & nCount_RNA < nCount_upper & percent.mt < percent.mt_upper & percent.hb < percent.hb_upper)

obj$barcode <- colnames(obj)
write.csv(obj@meta.data, file='D:Scanpy/Steele_2020/metadata.csv', quote=F, row.names=F)

counts_matrix <- GetAssayData(obj, assay='RNA', slot='counts')
writeMM(counts_matrix, file='D:Scanpy/Steele_2020/matrix.mtx')

write.table(
  data.frame('gene'=rownames(counts_matrix)),file='D:Scanpy/Steele_2020/genes.csv',
  quote=F,row.names=F,col.names=F
)


#Pre-process Lin and Lee

T48 <- Read10X(data.dir = "D:scRNA - Human PDAC/PDAC_17")
T48 <- CreateSeuratObject(counts = T48, project = 'T48', min.cells = 3, min.features = 200)

T49 <- Read10X(data.dir = "D:scRNA - Human PDAC/PDAC_18")
T49 <- CreateSeuratObject(counts = T49, project = 'T49', min.cells = 3, min.features = 200)

T50 <- Read10X(data.dir = "D:scRNA - Human PDAC/PDAC_19")
T50 <- CreateSeuratObject(counts = T50, project = 'T50', min.cells = 3, min.features = 200)

T51 <- Read10X(data.dir = "D:scRNA - Human PDAC/PDAC_20")
T51 <- CreateSeuratObject(counts = T51, project = 'T51', min.cells = 3, min.features = 200)

T52 <- Read10X(data.dir = "D:scRNA - Human PDAC/PDAC_21")
T52 <- CreateSeuratObject(counts = T52, project = 'T52', min.cells = 3, min.features = 200)

T53 <- Read10X(data.dir = "D:scRNA - Human PDAC/PDAC_22")
T53 <- CreateSeuratObject(counts = T53, project = 'T53', min.cells = 3, min.features = 200)

T54 <- Read10X(data.dir = "D:scRNA - Human PDAC/PDAC_23")
T54 <- CreateSeuratObject(counts = T54, project = 'T54', min.cells = 3, min.features = 200)

T55 <- Read10X(data.dir = "D:scRNA - Human PDAC/PDAC_24")
T55 <- CreateSeuratObject(counts = T55, project = 'T55', min.cells = 3, min.features = 200)

T56 <- Read10X(data.dir = "D:scRNA - Human PDAC/PDAC_25")
T56 <- CreateSeuratObject(counts = T56, project = 'T56', min.cells = 3, min.features = 200)

T57 <- Read10X(data.dir = "D:scRNA - Human PDAC/PDAC_26")
T57 <- CreateSeuratObject(counts = T57, project = 'T57', min.cells = 3, min.features = 200)

MET1 <- Read10X(data.dir = "D:scRNA - Human PDAC/MET_01")
MET1 <- CreateSeuratObject(counts = MET1, project = 'MET1', min.cells = 3, min.features = 200)

MET3 <- Read10X(data.dir = "D:scRNA - Human PDAC/MET_03")
MET3 <- CreateSeuratObject(counts = MET3, project = 'MET3', min.cells = 3, min.features = 200)

MET4 <- Read10X(data.dir = "D:scRNA - Human PDAC/MET_04/")
MET4 <- CreateSeuratObject(counts = MET4, project = 'MET4', min.cells = 3, min.features = 200)

MET5 <- Read10X(data.dir = "D:scRNA - Human PDAC/MET_05/")
MET5 <- CreateSeuratObject(counts = MET5, project = 'MET5', min.cells = 3, min.features = 200)

MET6 <- Read10X(data.dir = "D:scRNA - Human PDAC/MET_06/")
MET6 <- CreateSeuratObject(counts = MET6, project = 'MET6', min.cells = 3, min.features = 200)

MET7 <- Read10X(data.dir = "D:scRNA - Human PDAC/GSE156405/Liver_MET")
MET7 <- CreateSeuratObject(counts = MET7, project = 'MET7', min.cells = 3, min.features = 200)

MET8 <- Read10X(data.dir = "D:scRNA - Human PDAC/GSE156405/Lung_MET")
MET8 <- CreateSeuratObject(counts = MET8, project = 'MET8', min.cells = 3, min.features = 200)

MET9 <- Read10X(data.dir = "D:scRNA - Human PDAC/GSE156405/Vaginal_MET")
MET9 <- CreateSeuratObject(counts = MET9, project = 'MET9', min.cells = 3, min.features = 200)

T58 <- Read10X(data.dir = "D:scRNA - Human PDAC/GSE156405/Biopsy_P1")
T58 <- CreateSeuratObject(counts = T58, project = 'T58', min.cells = 3, min.features = 200)

T59 <- Read10X(data.dir = "D:scRNA - Human PDAC/GSE156405/Biopsy_P2")
T59 <- CreateSeuratObject(counts = T59, project = 'T59', min.cells = 3, min.features = 200)

T60 <- Read10X(data.dir = "D:scRNA - Human PDAC/GSE156405/Biopsy_P3")
T60 <- CreateSeuratObject(counts = T60, project = 'T60', min.cells = 3, min.features = 200)

T61 <- Read10X(data.dir = "D:scRNA - Human PDAC/GSE156405/Biopsy_P4")
T61 <- CreateSeuratObject(counts = T61, project = 'T61', min.cells = 3, min.features = 200)

T62 <- Read10X(data.dir = "D:scRNA - Human PDAC/GSE156405/Biopsy_P5")
T62 <- CreateSeuratObject(counts = T62, project = 'T62', min.cells = 3, min.features = 200)


counts <- readMM("D:scRNA - Human PDAC/GSE156405/Liver_MET/matrix.mtx")
genes <- read_tsv("D:scRNA - Human PDAC/GSE156405/Liver_MET/features.tsv", col_names = FALSE)
cell_ids <- read_tsv("D:scRNA - Human PDAC/GSE156405/Liver_MET/barcodes.tsv", col_names = FALSE)

rownames(counts) <- genes$X2
colnames(counts) <- cell_ids$X1

MET7 <- CreateSeuratObject(counts = counts, project = "MET7", min.cells = 3, min.features = 200)

T48$Condition <- "T"
T48$Patient <- "T48"
T49$Condition <- "T"
T49$Patient <- "T49"
T50$Condition <- "T"
T50$Patient <- "T50"
T51$Condition <- "T"
T51$Patient <- "T51"
T52$Condition <- "T"
T52$Patient <- "T52"
T53$Condition <- "T"
T53$Patient <- "T53"
T54$Condition <- "T"
T54$Patient <- "T54"
T55$Condition <- "T"
T55$Patient <- "T55"
T56$Condition <- "T"
T56$Patient <- "T56"
T57$Condition <- "T"
T57$Patient <- "T57"


MET1$Condition <- "MET"
MET1$Patient <- "MET1"
MET3$Condition <- "MET"
MET3$Patient <- "MET3"
MET4$Condition <- "MET"
MET4$Patient <- "MET4"
MET5$Condition <- "MET"
MET5$Patient <- "MET5"
MET6$Condition <- "MET"
MET6$Patient <- "MET6"



obj <- merge(T48, y = c(T49,T50,T51,T52,T53,T54,T55,T56,T57,MET1,MET3,MET4,MET5,MET6), add.cell.ids = c("T48","T49","T50","T51","T52","T53","T54","T55","T56",
                                                                                                                          "T57","MET1","MET3","MET4","MET5","MET6"),
             project = "Lin")

#Mito/Hemoglobin/Ribosomal ratio
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
obj[["percent.hb"]] <- PercentageFeatureSet(obj, pattern = "^HBA|^HBB")
obj[["percent.rp"]] <- PercentageFeatureSet(obj, pattern = "^RPS|^RPL")
obj$author <- "Lin_2020"

View(obj@meta.data)

obj@meta.data$orig.ident <- NULL

#Filter low quality cells
obj <- subset(obj, subset = nFeature_RNA > nFeature_lower & nFeature_RNA < nFeature_upper & nCount_RNA > nCount_lower & nCount_RNA < nCount_upper & percent.mt < percent.mt_upper & percent.hb < percent.hb_upper)

obj$barcode <- colnames(obj)
write.csv(obj@meta.data, file='D:Scanpy/Lin_2020/metadata.csv', quote=F, row.names=F)

counts_matrix <- GetAssayData(obj, assay='RNA', slot='counts')
writeMM(counts_matrix, file='D:Scanpy/Lin_2020/matrix.mtx')

write.table(
  data.frame('gene'=rownames(counts_matrix)),file='D:Scanpy/Lin_2020/genes.csv',
  quote=F,row.names=F,col.names=F
)



T58$Condition <- "T"
T58$Patient <- "T58"
T59$Condition <- "T"
T59$Patient <- "T59"
T60$Condition <- "T"
T60$Patient <- "T60"
T61$Condition <- "T"
T61$Patient <- "T61"
T62$Condition <- "T"
T62$Patient <- "T62"


MET7$Condition <- "MET"
MET7$Patient <- "MET7"
MET8$Condition <- "MET"
MET8$Patient <- "MET8"
MET9$Condition <- "MET"
MET9$Patient <- "MET9"




obj <- merge(T58, y = c(T59,T60,T61,T62,MET7,MET8,MET9), add.cell.ids = c("T58","T59","T60","T61","T62", "MET7","MET8","MET9"),
             project = "Lee")

#Mito/Hemoglobin/Ribosomal ratio
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
obj[["percent.hb"]] <- PercentageFeatureSet(obj, pattern = "^HBA|^HBB")
obj[["percent.rp"]] <- PercentageFeatureSet(obj, pattern = "^RPS|^RPL")
obj$author <- "Lee_2021"

View(obj@meta.data)

obj@meta.data$orig.ident <- NULL

#Filter low quality cells
obj <- subset(obj, subset = nFeature_RNA > nFeature_lower & nFeature_RNA < nFeature_upper & nCount_RNA > nCount_lower & nCount_RNA < nCount_upper & percent.mt < percent.mt_upper & percent.hb < percent.hb_upper)

obj$barcode <- colnames(obj)
write.csv(obj@meta.data, file='D:Scanpy/Lee_2021/metadata.csv', quote=F, row.names=F)

counts_matrix <- GetAssayData(obj, assay='RNA', slot='counts')
writeMM(counts_matrix, file='D:Scanpy/Lee_2021/matrix.mtx')

write.table(
  data.frame('gene'=rownames(counts_matrix)),file='D:Scanpy/Lee_2021/genes.csv',
  quote=F,row.names=F,col.names=F
)