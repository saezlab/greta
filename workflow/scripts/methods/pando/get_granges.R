# Parse args
args <- commandArgs(trailingOnly = F)
path_hg <- args[6]
path_mm <- args[7]

library(EnsDb.Hsapiens.v86)
gene.ranges_hg <- Signac::GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
write.csv(gene.ranges_hg, path_hg, row.names=FALSE)

library(EnsDb.Mmusculus.v79)
gene.ranges_mm <- Signac::GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
write.csv(gene.ranges_mm, path_mm, row.names=FALSE)
