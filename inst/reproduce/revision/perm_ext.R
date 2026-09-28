## 第 3 步：置换扩展（C3 / R1-M2）
## 与主分析同一暴露残差零分布（全协变量线性残差化、层内置换、种子 500000+b），b = 301..2300，
## 只拟合合理设定（+行为社会 / 全模型 × 9 套非饱和定义 × 3 编码 × 3 样本；插补；加权）= 每暴露 162 个。
## 这 2,000 次与用于判定合理性的前 300 次互相独立，因此可在其上做不受选择偏倚影响的检验与多重校正。
BASE <- Sys.getenv("NH_BASE"); source(file.path(BASE,"3_分析脚本/spec_core2.R"))
suppressMessages(library(parallel))
B0 <- as.integer(Sys.getenv("PX_FROM","301")); B1 <- as.integer(Sys.getenv("PX_TO","2300")); NC <- as.integer(Sys.getenv("PX_CORES","45"))
d0 <- prep(); DIR <- file.path(BASE,"4_结果/step3/置换扩展"); dir.create(DIR, recursive=TRUE, showWarnings=FALSE)
COVF <- COV[["全模型"]]; Xc <- model.matrix(as.formula(paste("~", paste(COVF, collapse="+"))), data=d0)
L <- sapply(EXPO, function(ex) log(d0[[ex]])); fitL <- Xc %*% qr.coef(qr(Xc), L); resL <- L - fitL
st <- split(seq_len(nrow(d0)), d0$SDMVSTRA)
COV <<- COV[c("+行为社会","全模型")]                                   # 残差化用全集之后再限定拟合的协变量集
OK9 <- OUTC[setdiff(names(OUTC), c("EFP 1mm","EFP 2mm"))]
one <- function(b){
  fp <- file.path(DIR, sprintf("px_%04d.csv", b)); if (file.exists(fp)) return(invisible(NULL))
  set.seed(500000+b); idx <- seq_len(nrow(d0)); for (s in st) idx[s] <- if (length(s)>1) sample(s) else s
  d <- d0; Lp <- fitL + resL[idx,,drop=FALSE]; for (j in seq_along(EXPO)) d[[EXPO[j]]] <- exp(Lp[,j])
  R <- run_grid2(d, miss="MICE 插补", wt="加权", outc=OK9) %>% mutate(rep=b, logor=signif(log(or),5), p=signif(p,5)) %>%
       select(rep, exposure, outcome, covset, scale, sample, logor, p)
  write.csv(R, fp, row.names=FALSE); invisible(NULL)
}
if (nzchar(Sys.getenv("PX_TEST"))) { t1 <- Sys.time(); one(B0); print(Sys.time()-t1); quit() }
t0 <- Sys.time(); cat("置换扩展 |", B0, "-", B1, "| 核", NC, "\n", flush=TRUE)
ORD <- if (nzchar(Sys.getenv("PX_REV"))) B1:B0 else B0:B1
invisible(mclapply(ORD, function(b) try(one(b), silent=TRUE), mc.cores=NC, mc.preschedule=FALSE))
cat("完成", length(list.files(DIR)), "次，用时", round(difftime(Sys.time(),t0,units="mins"),1), "分钟\n")
