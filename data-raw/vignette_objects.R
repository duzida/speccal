## 预先计算 vignette 中耗时的两步（网格拟合、100 次 Freedman-Lane 置换），存为 inst/extdata/vignette_objects.rds
## 在包根目录运行：OMP_NUM_THREADS=1 Rscript data-raw/vignette_objects.R [cores]
suppressMessages({devtools::load_all(".", quiet = TRUE); library(survey)})
args <- commandArgs(trailingOnly = TRUE); CORES <- if (length(args)) as.integer(args[1]) else 1L
options(survey.lonely.psu = "adjust")
des <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WT_TOTAL, data = nhanes_perio, nest = TRUE)
source("vignettes/axes.R")                         # 与 vignette 共用的轴定义
g   <- spec_fit(ax, nhanes_perio, design = des, exposure = c("WBC", "NLR", "AGR"), cores = CORES)
nul <- spec_null(g, scheme = "fl", B = 100, stratum = "SDMVSTRA", covariates = ax$covariates$full, seed = 1, cores = CORES)
strip <- function(x) { attr(x, "refit") <- NULL; attr(x, "context") <- NULL; x }   # 去掉闭包与数据副本，保持文件小
saveRDS(list(g = strip(g), nul = nul, made = Sys.time(), speccal = as.character(utils::packageVersion("speccal"))),
        "inst/extdata/vignette_objects.rds", compress = "xz")
cat("saved", round(file.size("inst/extdata/vignette_objects.rds") / 1e6, 2), "MB\n")
