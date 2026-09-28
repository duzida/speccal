## B1 敏感性：换残差化模型重建零分布（R1-M1 / R2-M1 / R3-M2）
## RMODEL=socio  : 暴露只对 +行为社会 9 个协变量残差化
## RMODEL=g2core : 只对候选因果图 G2 的 9 个核心混杂（不含 6 个争议协变量）残差化
## RMODEL=spline : 全协变量，两个连续变量用自然样条（年龄 5 df 且与性别交互、HbA1c 4 df；PIR、BMI 在数据中已是分类变量）
## 其余与主分析 fl_null_full.R 完全一致：层内置换残差、插补+加权子网格 3,960 个设定
BASE <- Sys.getenv("NH_BASE"); source(file.path(BASE,"3_分析脚本/spec_core2.R"))
suppressMessages({library(parallel); library(splines)})
RM <- Sys.getenv("RMODEL"); B <- as.integer(Sys.getenv("FL_B","300")); NC <- as.integer(Sys.getenv("FL_CORES","22"))
d0 <- prep(); DIR <- file.path(BASE, "4_结果/残差化敏感性", RM); dir.create(DIR, recursive=TRUE, showWarnings=FALSE)
CORE <- c("Age","Gender","Race","Educational_level","PIR","Marriage_ststus","Smoking_status","Drinking_status","BMI")
rhs <- switch(RM,
  socio  = paste(COV[["+行为社会"]], collapse="+"),
  g2core = paste(CORE, collapse="+"),
  spline = paste(c("ns(Age,5)*Gender","ns(HbA1c,4)",
                   setdiff(COV[["全模型"]], c("Age","Gender","HbA1c"))), collapse="+"))
Xc <- model.matrix(as.formula(paste("~", rhs)), data=d0)
L <- sapply(EXPO, function(ex) log(d0[[ex]])); fitL <- Xc %*% qr.coef(qr(Xc), L); resL <- L - fitL
cat("残差化模型", RM, "| 列数", ncol(Xc), "| R2:", paste(round(1 - colSums(resL^2)/colSums(scale(L,scale=FALSE)^2),3), collapse=" "), "\n", flush=TRUE)
st <- split(seq_len(nrow(d0)), d0$SDMVSTRA)
off <- c(socio=1e6, g2core=2e6, spline=3e6)[[RM]]
one <- function(b){
  fp <- file.path(DIR, sprintf("fl_%03d.csv", b)); if (file.exists(fp)) return(invisible(NULL))
  set.seed(off + b); idx <- seq_len(nrow(d0)); for (s in st) idx[s] <- if (length(s)>1) sample(s) else s
  d <- d0; Lp <- fitL + resL[idx,,drop=FALSE]; for (j in seq_along(EXPO)) d[[EXPO[j]]] <- exp(Lp[,j])
  R <- run_grid2(d, miss="MICE 插补", wt="加权") %>% mutate(rep=b) %>%
       select(rep, exposure, outcome, covset, scale, sample, or, p, rr, p_rr, sep_flag)
  write.csv(R, fp, row.names=FALSE); invisible(NULL)
}
t0 <- Sys.time(); cat("B =", B, "| 核", NC, "\n", flush=TRUE)
ORD <- if (nzchar(Sys.getenv("REV"))) B:1 else 1:B
invisible(mclapply(ORD, function(b) try(one(b), silent=TRUE), mc.cores=NC, mc.preschedule=FALSE))
cat("完成", length(list.files(DIR)), "次，用时", round(difftime(Sys.time(),t0,units="mins"),1), "分钟\n")
