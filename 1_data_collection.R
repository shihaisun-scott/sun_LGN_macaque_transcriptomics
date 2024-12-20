# collect and save data and metadata clustering ids
# make sure to set the working directory to this file's folder

rm(list = ls())

require(Matrix)

# make sure the following files are downloaded and saved in /data
# macaque_LGN_2021_exon-matrix.csv
# macaque_LGN_2021_intron-matrix.csv
# data/macaque_LGN_2021_metadata.csv

# load and save transcripts
exon<- read.csv("data/macaque_LGN_2021_exon-matrix.csv", header = TRUE, row.names = 1)  # Adjust file name and path as necessary
intron <- read.csv("data/macaque_LGN_2021_intron-matrix.csv", header = TRUE, row.names = 1)  # Adjust file name and path as necessary
mac_dat=exon+intron

mac_dat=sweep(mac_dat,2,colSums(mac_dat),"/")*10^6
mac_dat2=Matrix(as.matrix(mac_dat))
rownames(mac_dat2)=rownames(mac_dat)
colnames(mac_dat2)=colnames(mac_dat)
mac_dat2@x=mac_dat2@x/rep.int(Matrix::colSums(mac_dat2),diff(mac_dat2@p))*10^6
datlist=mac_dat2

# remove unnecessary genes
keepgen = setdiff(rownames(datlist), grep("^LOC|^MT-|^RPL|^RPS|^SNAR|^CT45|^MIR|^MTND|^TRNA", rownames(datlist), val = TRUE))
datlist=datlist[keepgen,]

metadata <- read.csv("data/macaque_LGN_2021_metadata.csv")
metadata$sample_name = gsub("-", ".", metadata$sample_name) # R script changes - to .
metadata$cluster_label[metadata$cluster_label == "Pulv"] <- "_Pulv"# change pulvinar labeling to _pul for easier identification down the track
metadata$cluster_label[grepl("^G", metadata$cluster_label)] <- "GABA"# combined GABAs
rownames(metadata)=metadata$sample_name

# remove low-quality cells
idx = metadata$class_label != 'Low Quality'
keepcells1 = colnames(datlist[,idx])
metadata <- metadata[keepcells1,]
datlist <- datlist[,keepcells1]

metalist=data.frame(metadata)


save(datlist,metalist,file="data/data.rda")



# load and save precluster information (clustering from Bakken et al.)
precluster <- data.frame(sample_name = metadata$sample_name, 
                    cluster_id = metadata$cluster_id, 
                    cluster_label = metadata$cluster_label, 
                    cluster_color = metadata$cluster_color, 
                    class_label = metadata$class_label, 
                    donor = metadata$donor, 
                    roi = metadata$roi, 
                    gender = metadata$gender, 
                    age = metadata$age, 
                    species = metadata$species_label)
rownames(precluster) <- precluster$sample_name


save(precluster,file="data/precluster_info.rda")

