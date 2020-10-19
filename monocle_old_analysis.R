OPC_d77=read.table("OPC_d77_dge.txt.gz", header = T, row.names=1)
OPC_d89=read.table("OPC_d89_dge.txt.gz", header = T, row.names=1)
OPC_d104=read.table("OPC_d104_dge.txt.gz", header = T, row.names=1)
OPC_O4=read.table("O4_dge.txt.gz", header =T, row.names=1)

d77_cells=data.frame(row.names=colnames(OPC_d77), rep.int(77, 1500))

d89_cells=data.frame(row.names=colnames(OPC_d89), rep.int(89, 2500))
O4_cells=data.frame(row.names=colnames(OPC_O4), rep.int(89, 1500))
d104_cells=data.frame(row.names=colnames(OPC_d104), rep.int(104, 1250))
fd_d77=data.frame(row.names=rownames(OPC_d77), x=rownames(OPC_d77))
d77_cells=new("AnnotatedDataFrame", data=d77_cells)

fd_d77=new("AnnotatedDataFrame", data=fd_d77)
OPC_d77_cds=newCellDataSet(as.matrix(OPC_d77), phenoData=d77_cells, featureData=fd_d77)
fd_d89=data.frame(row.names=rownames(OPC_d89), x=rownames(OPC_d89))
fd_O4=data.frame(row.names=rownames(OPC_O4), x=rownames(OPC_O4))
fd_d104=data.frame(row.names=rownames(OPC_d104), x=rownames(OPC_d104))
d89_cells=new("AnnotatedDataFrame", data=d89_cells)
fd_d89=new("AnnotatedDataFrame", data=fd_d89)
fd_O4=new("AnnotatedDataFrame", data=fd_O4)
O4_cells=new("AnnotatedDataFrame", data=O4_cells)
d104_cells=new("AnnotatedDataFrame", data=d104_cells)
fd_d104=new("AnnotatedDataFrame", data=fd_d104)
OPC_d89_cds=newCellDataSet(as.matrix(OPC_d89), phenoData=d89_cells, featureData=fd_d89)
OPC_O4_cds=newCellDataSet(as.matrix(OPC_O4), phenoData=O4_cells, featureData=fd_O4)
OPC_d104_cds=newCellDataSet(as.matrix(OPC_d104), phenoData=d104_cells, featureData=fd_d104)
OPCs=read.table("DS_monocle_merge.csv")

pd=data.frame(row.names=colnames(OPCs), rep(NA, 5186))
pd[1:1104,][is.na(pd[1:1104,])]=89
pd[1105:1954,][is.na(pd[1105:1954,])]=77
pd[1954:4172,][is.na(pd[1954:4172,])]=89
pd[4173:5186,][is.na(pd[4173:5186,])]=104

fd=data.frame(row.names=rownames(OPCs), x=rownames(OPCs))
colnames(fd)="gene"
pd=new("AnnotatedDataFrame", data=pd)
fd=new("AnnotatedDataFrame", data=fd)
OPCs_cds=newCellDataSet(as.matrix(OPCs), phenoData=pd, featureData=fd, expressionFamily=negbinomial.size())
OPCs_cds=estimateSizeFactors(OPCs_cds)
OPCs_cds=estimateDispersions(OPCs_cds)

OPCs_cds=detectGenes(OPCs_cds, min_expr=0.1)
expressed_genes=row.names(subset(fData(OPCs_cds), num_cells_expressed >=3))
valid_cells=row.names(subset(pData(OPCs_cds), num_genes_expressed > 250))

OPCs_cds_all_cells=OPCs_cds
OPCs_cds=OPCs_cds[,valid_cells]

pData(OPCs_cds)$total_transcripts=Matrix::colSums(exprs(OPCs_cds))
upper_bound=10^(mean(log10(pData(OPCs_cds)$total_transcripts)) + 2*sd(log10(pData(OPCs_cds)$total_transcripts)))
lower_bound=10^(mean(log10(pData(OPCs_cds)$total_transcripts)) - 2*sd(log10(pData(OPCs_cds)$total_transcripts)))

OPCs_cds_transcript_filt=OPCs_cds[,pData(OPCs_cds)$total_transcripts < upper_bound & pData(OPCs_cds)$total_transcripts > lower_bound]
OPCs_cds_pre_transcript_filt=OPCs_cds
OPCs_cds_pre_transcript_filt
OPCs_cds=OPCs_cds_transcript_filt
OPCs_cds

# dataset is filtered based on gene counts and transcript counts -> go on to analysis
diff_test_res=differentialGeneTest(OPCs_cds[expressed_genes,], fullModelFormulaStr="~age")
ordering_genes=row.names(subset(diff_test_res, qval < 0.01))
OPCs_cds=setOrderingFilter(OPCs_cds, ordering_genes)
plot_ordering_genes(OPCs_cds)
OPCs_cds=reduceDimension(OPCs_cds, max_components=2, method="DDRTree")
OPCs_cds=orderCells(OPCs_cds)
plot_cell_trajectory(OPCs_cds, color_by="age")
plot_cell_trajectory(OPCs_cds, color_by="State")
genes=c("GFAP", "OLIG1", "OLIG2", "SOX10", "PDGFRA", "MBP", "PLP1", "MIR219A2")
plot_genes_in_pseudotime(OPCs_cds[genes,], color_by="State")
# lower number states are oligodendrocytes, purpley color states are astrocytes
# left branch: oligos, upper branch: astros, right branch: OPCs
# pseudotime out of order -> reverse so that left most point as pseudotime=0
# adapt function from tutorial to set pseudotime=0 in state with most d77 cells
State0=function(cds){
   if (length(unique(pData(cds)$State)) > 1){
     d77_counts <- table(pData(cds)$State, pData(cds)$age)[,"77"]
     return(as.numeric(names(d77_counts)[which
           (d77_counts == max(d77_counts))]))
   } else {
    return (1)
   }
 }
OPCs_cds=orderCells(OPCs_cds, root_state=State0(OPCs_cds))
plot_cell_trajectory(OPCs_cds, color_by="Pseudotime")
node23=BEAM(OPCs_cds, branch_point=23)
node23=node23[order(node23$qval),]
node23=node23[,c("gene", "pval", "qval")]
plot_genes_branched_heatmap(OPCs_cds[row.names(head(node23, 50)),], branch_point=23, show_rownames=T)
pseudotime_diff_res_test=differentialGeneTest(OPCs_cds, fullModelFormulaStr="~sm.ns(Pseudotime)")
# take top 350 genes of node23 (all significant with correction) and cross reference with human TF list
# input into branch heatmap
TFs=scan("TFs_human.txt", what="character")
node23_TFs=intersect(row.names(head(node23, 350)), TFs)
node23_TFs

plot_genes_branched_pseudotime(OPCs_cds[c("TSC22D4", "SOX10", "E2F3", "ZEB2", "RBPJ"),], branch_point=23)
plot_genes_branched_pseudotime(OPCs_cds[c("TCF7L2", "SOX6", "SOX11", "TSC22D1"),], branch_point=23)
plot_genes_branched_pseudotime(OPCs_cds[c("HES1", "ZFP36L1", "YBX1"),], branch_point=23, color_by="State")
plot_genes_branched_pseudotime(OPCs_cds[c("HES1", "ZFP36L1", "YBX1"),], branch_point=23, color_by="State")
plot_cell_trajectory(OPCs_cds, color_by="State")
plot_cell_trajectory(OPCs_cds, color_by="State")
plot_cell_trajectory(OPCs_cds, color_by="age")

# 5/20/18
# get sig lncRNAs and miRs at node 23 for heatmap (set lower sig threshold to include more)
# get list of lncRNAs online and merge with those named "LINC..."
lncRNAs=scan("lncRNA_human.txt", what="character")
lncRNAs_OPCs=union(lncRNAs, grep("LINC", row.names(OPCs), value=T))
node_23_lncRNAs=intersect(node23$gene, lncRNAs_OPCs)
node23[node_23_lncRNAs,]
plot_genes_branched_heatmap(OPCs_cds[row.names(node_23_sig_lncRNAs),], branch_point=23, show_rownames=T)
miRs=grep("^MIR", row.names(OPCs), value=T)
node_23_miRs=intersect(node23$gene, miRs)
node23[node_23_miRs,]
jpeg("rplot_node_23_miRs.jpg", res=100)
plot_genes_branched_heatmap(OPCs_cds[row.names(node_23_miRs),], branch_point=23, show_rownames=T, num_clusters=2)

## subset out OPC and O4 datasets, use previous filtering thresholds
OPC_only=OPCs[,grep(pattern="^OPC", x=colnames(OPCs), value=T)]
O4=OPCs[,grep(pattern="^O4", x=colnames(OPCs), value=T)]
pd_OPCs=data.frame(row.names=colnames(OPC_only), rep(NA, ncol(OPC_only)))
pd_O4=data.frame(row.names=colnames(O4), rep("89", ncol(O4)))
colnames(pd_O4)="age"
pd_OPCs[1:850,][is.na(pd_OPCs[1:850,])]="77"
pd_OPCs[851:3068,][is.na(pd_OPCs[851:3068,])]="89"
pd_OPCs[3069:4082,][is.na(pd_OPCs[3069:4082,])]="104"
pd_OPCs=new("AnnotatedDataFrame", data=pd_OPCs)
pd_O4=new("AnnotatedDataFrame", data=pd_O4)
OPC_only_cds=newCellDataSet(as.matrix(OPC_only), phenoData=pd_OPCs, featureData=fd, expressionFamily=negbinomial.size())
O4_cds=newCellDataSet(as.matrix(O4), phenoData=pd_O4, featureData=fd, expressionFamily=negbinomial.size())
OPC_only_cds=estimateSizeFactors(OPC_only_cds)
OPC_only_cds=estimateDispersions(OPC_only_cds)
O4_cds=estimateSizeFactors(O4_cds)
O4_cds=estimateDispersions(O4_cds)
OPC_only_cds=detectGenes(OPC_only_cds, min_expr=0.1)
O4_cds=detectGenes(O4_cds, min_expr=0.1)
expressed_genes_OPC=rownames(subset(fData(OPC_only_cds), num_cells_expressed>=3))
expressed_genes_O4=rownames(subset(fData(O4_cds), num_cells_expressed>=3))
length(expressed_genes_O4)
valid_cells_OPC=row.names(subset(pData(OPC_only_cds), num_genes_expressed>250))
valid_cells_O4=row.names(subset(pData(O4_cds), num_genes_expressed>250))
OPC_only_cds=OPC_only_cds[,valid_cells_OPC]
O4_cds=O4_cds[,valid_cells_O4]
pData(OPC_only_cds)$total_transcripts=Matrix::colSums(exprs(OPC_only_cds))
pData(O4_cds)$total_transcripts=Matrix::colSums(exprs(O4_cds))
CellDataSet (storageMode: environment)
OPC_only_cds=OPC_only_cds[,pData(OPC_only_cds)$total_transcripts > lower_bound & pData(OPC_only_cds)$total_transcripts < upper_bound]
OPC_only_cds
O4_cds=O4_cds[,pData(O4_cds)$total_transcripts > lower_bound & pData(O4_cds)$total_transcripts < upper_bound]
diff_test_res_OPC=differentialGeneTest(OPC_only_cds[expressed_genes_OPC,], fullModelFormulaStr="~age")
ordering_genes_OPC=row.names(subset(diff_test_res_OPC, qval<0.01))
length(ordering_genes_OPC)
OPC_only_cds=setOrderingFilter(OPC_only_cds, ordering_genes_OPC)
O4_cds=setOrderingFilter(O4_cds, ordering_genes_OPC)
OPC_only_cds=reduceDimension(OPC_only_cds, max_components=2, method="DDRTree")
O4_cds=reduceDimension(O4_cds, max_components=2, method="DDRTree")
OPC_only_cds=orderCells(OPC_only_cds)
O4_cds=orderCells(O4_cds)
plot_cell_trajectory(OPC_only_cds)
plot_cell_trajectory(OPC_only_cds)
plot_cell_trajectory(OPC_only_cds, color_by="age")
plot_cell_trajectory(OPC_only_cds, color_by="age")
plot_cell_trajectory(OPC_only_cds, color_by="Pseudotime")
plot_cell_trajectory(OPC_only_cds, color_by="Pseudotime")
plot_cell_trajectory(O4_cds)
plot_cell_trajectory(O4_cds)
plot_cell_trajectory(O4_cds, color_by="Pseudotime")
plot_cell_trajectory(O4_cds)
# set state 5 as root state
State0=function(cds){
  if (length(unique(pData(cds)$State)) > 1){
    state5_counts <- table(pData(cds)$State, pData(cds)$State)[,"5"]
    return(as.numeric(names(state5_counts)[which
          (state5_counts == max(state5_counts))]))
  } else {
    return (1)
  }
}
O4_cds=orderCells(O4_cds, root_state=State0(O4_cds))
plot_cell_trajectory(O4_cds, color_by="Pseudotime")
plot_cell_trajectory(O4_cds)
plot_cell_trajectory(O4_cds, color_by="Pseudotime")
plot_cell_trajectory(OPC_only_cds)
plot_genes_branched_pseudotime(OPC_only_cds[c(astrocyte_genes, oligo_genes),], branch_point=3, color_by=)
plot_genes_branched_pseudotime(OPC_only_cds[c(astrocyte_genes, oligo_genes),], branch_point=3, color_by=)
plot_genes_branched_pseudotime(O4_cds[c(astrocyte_genes, oligo_genes),], branch_point=2, color_by="State"> plot_cell_trajectory(O4_cds)
plot_genes_branched_pseudotime(O4_cds[c(astrocyte_genes, oligo_genes),], branch_point=2, color_by="State"> plot_cell_trajectory(O4_cds)
plot_genes_branched_pseudotime(O4_cds[c(astrocyte_genes, oligo_genes),], branch_point=2, color_by="State"> dev.off()
plot_cell_trajectory(OPC_only_cds)
BEAM_OPC_only_3=BEAM(OPC_only_cds, branch_point=3)
plot_cell_trajectory(O4_cds)
BEAM_O4_2=BEAM(O4_cds, branch_point=2)
BEAM_OPC_only_3=BEAM_OPC_only_3[c("gene_short_name", "pval", "qval")]
BEAM_O4_2=BEAM_O4_2[,c("gene_short_name", "pval", "qval")]
node23_cor=node23[order(node23$gene),]
cor_fullvsOPC=cor.test(node23_cor$pval, BEAM_OPC_only_3$pval, method="spearman")
cor_fullvsO4=cor.test(node23_cor$pval, BEAM_O4_2$pval, method="spearman")
cor_OPCvsO4=cor.test(BEAM_OPC_only_3$pval, BEAM_O4_2$pval, method="spearman")
cor_fullvsOPC
node23_cor=subset(node23_cor, node23_cor$pval<0.05)
BEAM_OPC_cor=subset(BEAM_OPC_only_3, BEAM_OPC_only_3$pval<0.05)
BEAM_O4_2=subset(BEAM_O4_2, BEAM_O4_2$pval<0.05)
dim(node23_cor)
dim(BEAM_OPC_cor)
dim(BEAM_O4_2)
genes_BEAM_cor=intersect(node23_cor$gene, BEAM_OPC_cor$gene_short_name)
genes_BEAM_cor=intersect(genes_BEAM_cor, BEAM_O4_2$gene_short_name)
length(genes_BEAM_cor)
node23_cor=node23_cor[genes_BEAM_cor,]
BEAM_OPC_cor=BEAM_OPC_cor[genes_BEAM_cor,]
BEAM_O4_2=BEAM_O4_2[genes_BEAM_cor,]
cor_fullvsOPC=cor.test(node23_cor$pval, BEAM_OPC_cor$pval, method="spearman")
cor_fullvsO4=cor.test(node23_cor$pval, BEAM_O4_2$pval, method="spearman")
cor_OPCvsO4=cor.test(BEAM_OPC_cor$pval, BEAM_O4_2$pval, method="spearman")
cor_fullvsOPC
cor_fullvsO4
cor_OPCvsO4

## filter out even more genes based on pval
BEAM_OPC_cor=subset(BEAM_OPC_only_3, BEAM_OPC_only_3$pval<0.0005)
node23_cor=subset(node23, node23$pval<0.0005)
node23_cor=subset(node23_cor, node23_cor$pval<0.0005)
BEAM_O4_cor=subset(BEAM_O4_2, BEAM_O4_2$pval<0.0005)
genes_BEAM_cor=intersect(node23_cor$gene, BEAM_OPC_cor$gene_short_name)
genes_BEAM_cor=intersect(genes_BEAM_cor, BEAM_O4_cor$gene_short_name)
length(genes_BEAM_cor)
node23_cor=node23_cor[genes_BEAM_cor,]
BEAM_OPC_cor=BEAM_OPC_cor[genes_BEAM_cor,]
BEAM_O4_cor=BEAM_O4_cor[genes_BEAM_cor,]
cor_fullvsOPC=cor.test(node23_cor$pval, BEAM_OPC_cor$pval, method="spearman")
cor_fullvsO4=cor.test(node23_cor$pval, BEAM_O4_cor$pval, method="spearman")\
cor_OPCvsO4=cor.test(BEAM_OPC_cor$pval, BEAM_O4_cor$pval, method="spearman")
cor_fullvsOPC
cor_fullvsO4
cor_OPCvsO4

node3_OPC_TFs=BEAM_OPC_only_3[TFs,]
node3_OPC_TFs=node3_OPC_TFs[order(node3_OPC_TFs$pval),]
node3_OPC_top20_TFs=head(node3_OPC_TFs, 20)
node3_OPC_top20_TFs=row.names(node3_OPC_top20_TFs)
node3_OPC_top20_TFs

# 8/14/18
BEAM_OPC_node1=BEAM(OPC_only_cds, branch_point=1)
BEAM_OPC_node2=BEAM(OPC_only_cds, branch_point=2)
BEAM_OPC_2v7=BEAM(OPC_only_cds, branch_states=c(2,7))

# 10/2/18
# look at all astrocyte states (2, 5, 7)
# run BEAM on 2v5, 2v7, 5v7 for differences
# for similarities:
OPC_state=pData(OPC_only_cds)[,c(1,6)]
OPC_state=as.matrix(OPC_state)
OPC_state=OPC_state[,-1] 
OPC_state=as.data.frame(OPC_state)
colnames(OPC_state)="state"
OPC_state2=subset(OPC_state, OPC_state$state==2)

OPC_state5=subset(OPC_state, OPC_state$state==5)
OPC_state7=subset(OPC_state, OPC_state$state==7)
OPC_state2=rownames(OPC_state2)
OPC_state5=rownames(OPC_state5)
OPC_state7=rownames(OPC_state7)
state2_GE=OPCs[,OPC_state2]
state5_GE=OPCs[,OPC_state5]
state7_GE=OPCs[,OPC_state7]
dim(state5_GE)
dim(state7_GE)
state2_total_GE=data.frame(rownames=row.names(state2_GE), data=rowSums(state2_GE)/ncol(state2_GE))
state2_total_GE=as.matrix(state2_total_GE)
dim(state2_total_GE)

state2_total_GE=state2_total_GE[,2]
state2_total_GE=as.data.frame(state2_total_GE)


state5_total_GE=data.frame(rownames=row.names(state5_GE), data=rowSums(state5_GE)/ncol(state5_GE))
state5_total_GE=as.matrix(state5_total_GE)
state5_total_GE=state5_total_GE[,2]
state5_total_GE=as.data.frame(state5_total_GE)

state7_total_GE=data.frame(rownames=row.names(state7_GE), data=rowSums(state7_GE)/ncol(state7_GE))
state7_total_GE=as.matrix(state7_total_GE)
state7_total_GE=state7_total_GE[,2]
state7_total_GE=as.data.frame(state7_total_GE)

state2_total_GE=as.matrix(state2_total_GE)
state5_total_GE=as.matrix(state5_total_GE)
state7_total_GE=as.matrix(state7_total_GE)
OPC_state6=subset(OPC_state, OPC_state$state==6)
OPC_state6=rownames(OPC_state6)
state6_GE=OPCs[, OPC_state6]
state6_total_GE=data.frame(rownames=row.names(state6_GE), data=rowSums(state6_GE)/ncol(state6_GE))
state6_total_GE=as.matrix(state6_total_GE)
state6_total_GE=state6_total_GE[,2]
state6_total_GE=as.matrix(state6_total_GE)

