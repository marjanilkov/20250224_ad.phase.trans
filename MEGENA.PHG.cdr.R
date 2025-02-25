# We will first build MEGENA networks for each CDR level (0, 0.5, 1, 2, 3, 4, 5)
# We will then try to see what differences appear between them as well as
# differences with the extablished AD subtypes by Ryan Neff et al.

rm(list = ls())

setwd("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/MSSM/2025/20250224_ad.phase.trans")

require(MEGENA)


#### Computation and module parameters ####
### These should be the same for all analyses
n.cores <- 7; # number of cores to be used for parallelization (need to MATCH those requested in your LSF).
doPar <- TRUE; # TRUE = perform parallel PFN -> MCA. FALSE = no parallelization
method = "pearson" # Currently “pearson” (Pearson’s correlation) and “spearman” (Spearman’s correlation available) for correlation.
FDR <- 0.05; # FDR threshold to identify significant interactions.
n.cor.perm = 10; # Number of permutations in correlation screening.
n.hub.perm = 100; # Number of permutations to calculate hub significance in MHA
mod.pval = 0.05; # module p-value threshold in MCA
hub.pval = 0.05; # hub p-value threshold in MHA
min.size = 10 # minimum module size
max.size = NULL # maximum module size
### annotation to be done on the downstream
annot.table=NULL
id.col = 1
symbol.col= 2

source("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/MSSM/sinaiRepo/marjanRfunctions.R")
source("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/MSSM/sinaiRepo/R-tomfunctions_static.R")

## Path to filtered RNA-seq data

# Generate Correlation Matrix
# Pearson Correlation Coefficients (PCCs) based on 100 permutations
print("Generating correlation matrix")

datExpr = readRDS("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/WGCNA_from_scratch/anew3/20220829_data/BM_36.RNA.CDR_age_sex_adjusted.std.RDS")
meta = readRDS("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/WGCNA_from_scratch/anew3/20220829_data/metadata.RDS")

cdr = 0

# use only cases with dementia and not MCI + dementia
meta_cdr = meta[meta$CDR == cdr,]

# create a common ID to be able to properly knit together these DFs
commonID = Reduce(intersect,list(meta_cdr$SynapseBrainID, colnames(datExpr)))#,rownames(deSex),rownames(deDxSex)))
# put everything together in a list
meta_cdr = meta_cdr[commonID,]
datExpr = datExpr[,commonID]
# for the downstream analysis we need the rows to be patients
#datExpr = as.data.frame(t(datExpr))
#######################################################################################
################################################################################
################################################################################
# Z-transform the data
# zmatrix = apply(datExpr, 2, Ztransform, 4) # 1 means the apply function is done on rows and 2 would mean columns
# rownames(zmatrix) = rownames(datExpr)
# datExpr = zmatrix
# dim(datExpr)
# rm(zmatrix)
# gc()

# datExpr = as.data.frame(t(datExpr))

ijw <- calculate.correlation(datExpr, 
                             doPerm = n.cor.perm, 
                             doPar = TRUE,
                             num.cores = n.cores,
                             output.corTable = FALSE,
                             output.permFDR = FALSE)


# Register multiple cores if needed: note that set.parallel.backend() is deprecated.
run.par = doPar & (getDoParWorkers() == 1)
if (run.par)
{
  cl <- parallel::makeCluster(n.cores)
  registerDoParallel(cl)
  # check how many workers are there
  cat(paste("number of cores to use:",getDoParWorkers(),"\n",sep = ""))
}

### Major Heavy Memory Calculation (will take a lot of time)
print("Computing PFN")
# Compute PFN
el <- calculate.PFN(ijw[,1:3],doPar = doPar,num.cores = n.cores)
# graph dataframe
g <- graph.data.frame(el,directed = FALSE)

## perform MCA clustering.
print("Performing MCA clustering")
MEGENA.output <- do.MEGENA(g,
                           mod.pval = mod.pval,hub.pval = hub.pval,remove.unsig = TRUE,
                           min.size = 10,max.size = vcount(g)/2,
                           doPar = doPar,num.cores = n.cores,n.perm = n.hub.perm,
                           save.output = FALSE)

# unregister cores as these are not needed anymore.
if (getDoParWorkers() > 1)
{
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}

summary.output <- MEGENA.ModuleSummary(MEGENA.output,
                                       mod.pvalue = mod.pval,hub.pvalue = hub.pval,
                                       min.size = min.size,max.size = vcount(g)/2,
                                       annot.table = annot.table,id.col = id.col,symbol.col = symbol.col,
                                       output.sig = TRUE)

module.output <- module_convert_to_table(MEGENA.output,mod.pval = 0.05,
                                         hub.pval = 0.05,min.size = 10,max.size=vcount(g)/2)

print(paste0("Saving data to file:MEGENA.Results.RData"))
save(summary.output,MEGENA.output,module.output,g,el,file=paste0("MEGENA.PHG.cdr",cdr,".RData"))

