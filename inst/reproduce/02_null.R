## 验收 2：用 speccal 的 spec_null(scheme="fl") 以论文原种子重跑若干次置换，与存档分片
##         4_结果/4.2_置换与零分布/FL零分布完整网格/fl_XXX.csv 逐设定核对（不重跑全部 300 次）
## 用法：OMP_NUM_THREADS=1 Rscript inst/reproduce/02_null.R [reps, 如 1:20] [cores]
suppressMessages({library(survey); library(dplyr); devtools::load_all(Sys.getenv("SPECCAL_DIR", "."), quiet = TRUE)})
args <- commandArgs(trailingOnly = TRUE)
REPS <- if (length(args) >= 1) eval(parse(text = args[1])) else 1:20
CORES <- if (length(args) >= 2) as.integer(args[2]) else 1L
B <- Sys.getenv("NH_BASE", "/Users/hp/hp_data/nhanes")
source(file.path(B, "3_分析脚本/spec_core.R")); d <- prep(); options(survey.lonely.psu = "adjust")

## 与 01_grid.R 相同的轴；零分布只在 插补+加权 子网格上算（论文口径）
ax <- spec_axes(
  outcome    = c(cdc_any = "y_cdc_any", cdc_modsev = "y_cdc_modsev", cdc_severe = "y_cdc_severe", cdc_allsite = "y_allsites",
                 efp1 = "y_efp1", efp2 = "y_efp2", efp3 = "y_efp3", efp4 = "y_efp4",
                 st4_pairs = "y_st4_pairs", st4_t20 = "y_st4_t20", st4_jaw10 = "y_st4_jaw10"),
  covariates = list(crude = character(0), demo = COV[["+人口学"]], socio = COV[["+行为社会"]], full = COV[["全模型"]]),
  coding     = c("per_sd", "log2", "q4_vs_q1"),
  missing    = list(imputed = TRUE, complete = quote(complete_case == 1)),
  sample     = list(all = TRUE, wbc_normal = quote(WBC >= 3.5 & WBC <= 11), age40 = quote(Age >= 40)),
  weighting  = c("design", "unweighted"))
des <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WT_TOTAL, data = d, nest = TRUE)
## 观察网格只需子网格（零分布用 refit，不依赖观察网格的内容）
g <- spec_fit(restrict_axes(ax, list(missing = "imputed", weight = "design")), d, des, exposure = EXPO, cores = CORES, quiet = TRUE)
t0 <- Sys.time()
nul <- spec_null(g, scheme = "fl", reps = REPS, stratum = "SDMVSTRA", covariates = COV[["全模型"]],
                 seed = 500000, cores = CORES)     # fl_null_full.R: set.seed(500000 + b)
cat("置换", length(REPS), "次，用时", round(difftime(Sys.time(), t0, units = "mins"), 1), "分钟\n"); print(nul)

MAP <- list(
  outcome = c(cdc_any = "CDC/AAP 任何", cdc_modsev = "CDC/AAP 中重度", cdc_severe = "CDC/AAP 重度", cdc_allsite = "CDC 全位点变体",
              efp1 = "EFP 1mm", efp2 = "EFP 2mm", efp3 = "EFP 3mm", efp4 = "EFP 4mm",
              st4_pairs = "EFP IV期(对颌对<10)", st4_t20 = "EFP IV期(余留<20)", st4_jaw10 = "EFP IV期(单颌<10)"),
  covset = c(crude = "粗模型", demo = "+人口学", socio = "+行为社会", full = "全模型"),
  coding = c(per_sd = "每 SD", log2 = "每翻倍(log2)", q4_vs_q1 = "Q4 vs Q1"),
  sample = c(all = "全样本", wbc_normal = "WBC 3.5-11", age40 = "年龄>=40"))
k <- nul$perm %>% transmute(rep, exposure, outcome = MAP$outcome[outcome], covset = MAP$covset[covset],
                            scale = MAP$coding[coding], sample = MAP$sample[sample], or = exp(estimate), p, rr, p_rr, sep_flag)
ref <- bind_rows(lapply(REPS, function(b) read.csv(file.path(B, "4_结果/4.2_置换与零分布/FL零分布完整网格", sprintf("fl_%03d.csv", b)), fileEncoding = "UTF-8")))
m <- inner_join(ref, k, by = c("rep", "exposure", "outcome", "covset", "scale", "sample"), suffix = c("_ref", "_pkg"))
cat("参照行", nrow(ref), "| 包输出行", nrow(k), "| 匹配", nrow(m), "\n"); stopifnot(nrow(m) == nrow(ref), nrow(k) == nrow(ref))
rel <- sapply(c("or", "p", "rr", "p_rr"), function(v) { a <- m[[paste0(v, "_ref")]]; b <- m[[paste0(v, "_pkg")]]
  ok <- is.finite(a) & is.finite(b); max(abs(a[ok] - b[ok]) / pmax(abs(a[ok]), 1e-12)) })
print(signif(rel, 3)); cat("sep_flag 一致:", mean(m$sep_flag_ref == m$sep_flag_pkg), "\n")
if (all(rel < 1e-6) && all(m$sep_flag_ref == m$sep_flag_pkg)) cat("\n验收通过：spec_null 以原种子重现存档的 FL 零分布分片\n") else stop("与存档不一致")
