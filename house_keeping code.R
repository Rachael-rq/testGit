# development codes
library(devtools)
tgfb_expr <- read.delim("/Users/ruqianlyu/Desktop/WEHI_Research Project/Packagedevelopment/All_exprs_combat.txt",header = TRUE)
rownames(tgfb_expr) <- tgfb_expr$EntrezID
tgfb_expr <- tgfb_expr[2:ncol(tgfb_expr)]
tgfb_gs <- read.delim("/Users/ruqianlyu/Desktop/WEHI_Research Project/Packagedevelopment/TGF_EMT_signature.txt", header = TRUE)
tgfb_gs_up <- tgfb_gs[tgfb_gs$upDown=='up',1]
length(tgfb_gs_up)
tgfb_gs_dn <- tgfb_gs[tgfb_gs$upDown=='down',1]
length(tgfb_gs_dn)

# subset the TGFb expression matrix and the TGFb-EMT signature's up-set and
# down-set for the examples

toy_expr <- tgfb_expr[1:20,1:2]
toy_up <- sample(rownames(toy_expr),5)
toy_dn <- sample(setdiff(rownames(toy_expr),toy_up),5)
toy_dn <- GSEABase::GeneSet(toy_dn)
toy_up <- GSEABase::GeneSet(toy_up)

use_data(toy_dn, overwrite = TRUE)
use_data(toy_up, overwrite = TRUE)
use_data(toy_expr, overwrite = TRUE)
ranked <- rankGenes(toy_expr)

scoredf <- simpleScore(ranked, upSet = toy_up, downSet = toy_dn)
n_up = length(GSEABase::geneIds(toy_up))
n_down = length(GSEABase::geneIds(toy_dn))
#no_cores <- detectCores() - 1
#This permutation function can be run using parallel scripts, refer to the
#vignette for examples
#cl <- makeCluster(no_cores)
#registerDoParallel(cl)
permuteResult = permute_null(n_up = n_up, n_down = n_down, ranked, B =10,
seed = 1)
pvals <- get_pval(permuteResult,scoredf)
plot_null(permuteResult,scoredf,pvals,sampleLabel = colnames(ranked)[2:5])

# Data set for vignettes

tgfb_expr_20 <- tgfb_expr[,1:20]

document()
