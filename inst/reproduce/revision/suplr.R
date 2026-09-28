## 第 3 步：拐点非线性的 sup-LR 检验（C5，Davies 问题）
## 统计量 = 线性模型偏差 − 41 个候选拐点中最小的两段式模型偏差（与 bp_raowu.R 同一网格、归一化权重）。
## 零分布：参数自助法，从加权线性 logistic 模型的拟合概率生成结局，在每个重复中重新搜索全部候选
## （Hansen 1996 的思路）。与朴素的 1 自由度卡方 P 值对照。重复内不含整群结构（见正文局限）。
BASE <- Sys.getenv("NH_BASE"); source(file.path(BASE,"3_分析脚本/spec_core.R"))
suppressMessages(library(parallel))
B <- as.integer(Sys.getenv("LR_B","999")); NC <- as.integer(Sys.getenv("LR_CORES","10"))
d0 <- prep(); COVF <- COV[["全模型"]]; EXS <- c("SII","dNLR","AGR"); NCAND <- 41
DIR <- file.path(BASE,"4_结果/step3/supLR"); dir.create(DIR, recursive=TRUE, showWarnings=FALSE)
Xc <- model.matrix(as.formula(paste("~", paste(COVF, collapse="+"))), data=d0)
w <- d0$WT_TOTAL/mean(d0$WT_TOTAL)
cand <- lapply(EXS, function(e) as.numeric(quantile(log2(d0[[e]]), seq(.05,.95,length.out=NCAND)))); names(cand) <- EXS
supstat <- function(y, ex){ v <- log2(d0[[ex]]); X0 <- cbind(Xc, v)
  m0 <- glm.fit(X0, y, weights=w, family=quasibinomial(), control=glm.control(maxit=100))
  dv <- sapply(cand[[ex]], function(c) { m <- try(glm.fit(cbind(X0, pmax(v-c,0)), y, weights=w, family=quasibinomial(), control=glm.control(maxit=100)), silent=TRUE)
    if (inherits(m,"try-error") || !m$converged) NA else m$deviance })
  list(stat=m0$deviance - min(dv, na.rm=TRUE), p0=m0$fitted.values) }
obs <- lapply(EXS, function(e) supstat(d0$y_cdc_any, e)); names(obs) <- EXS
## 加权 logistic 偏差需按离散度换算为似然比尺度：用线性模型的 Pearson 离散度
disp <- sapply(EXS, function(e){ v <- log2(d0[[e]]); m0 <- glm.fit(cbind(Xc,v), d0$y_cdc_any, weights=w, family=quasibinomial())
  sum(m0$weights*m0$residuals^2)/m0$df.residual })
one <- function(b){ fp <- file.path(DIR, sprintf("lr_%04d.csv", b)); if (file.exists(fp)) return(invisible(NULL))
  set.seed(700000+b)
  out <- do.call(rbind, lapply(EXS, function(e){ ys <- rbinom(nrow(d0), 1, obs[[e]]$p0); data.frame(b=b, exposure=e, stat=supstat(ys, e)$stat) }))
  write.csv(out, fp, row.names=FALSE); invisible(NULL) }
write.csv(data.frame(exposure=EXS, stat=sapply(obs, `[[`, "stat"), disp=disp, p_naive=pchisq(sapply(obs, `[[`, "stat")/disp, 1, lower.tail=FALSE)),
          file.path(DIR,"观察值.csv"), row.names=FALSE)
if (nzchar(Sys.getenv("LR_TEST"))) { t1 <- Sys.time(); one(1); print(Sys.time()-t1); print(read.csv(file.path(DIR,"观察值.csv"))); quit() }
t0 <- Sys.time(); cat("sup-LR | B =", B, "| 核", NC, "\n", flush=TRUE)
invisible(mclapply(1:B, function(b) try(one(b), silent=TRUE), mc.cores=NC, mc.preschedule=FALSE))
cat("完成", length(list.files(DIR, pattern="^lr_")), "次，用时", round(difftime(Sys.time(),t0,units="mins"),1), "分钟\n")
