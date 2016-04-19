LINCS <- read.table('reports/LINCS2000/result_UPTAG_summly_n7503.txt',header=T,sep="\t", stringsAsFactors = F, 
                    quote="", comment.char = "", na.strings = "")

LINCS <- LINCS[LINCS$pert_type=='trt_cp',]

unique(unlist(strsplit(LINCS$cell_id,split="\\|")))


length(unique(LINCS$pert_iname))

table(abs(LINCS$mean_rankpt_4) > 95)
table(LINCS$mean_rankpt_4 > 95)
table(LINCS$mean_rankpt_4 < -95)
