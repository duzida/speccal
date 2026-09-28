## 第 3 步：显著设定比例的 Rao–Wu 自助区间（C1）
## 每次重复用 Rao–Wu 再标定权重（subbootstrap：每层抽 n_h-1 个 PSU）重估 162 个合理设定的对数比值比；
## 显著性按原样本的设计基标准误和残差自由度判定（重复内无法再算设计基标准误：多数层只剩 1 个 PSU），
## 因此区间反映估计值的抽样变异，而把标准误视为固定。b = 0 为原权重，用于核对与主分析逐一相同。
BASE <- Sys.getenv("NH_BASE"); source(file.path(BASE,"3_分析脚本/spec_core2.R"))
suppressMessages(library(parallel))
B <- as.integer(Sys.getenv("SB_B","1000")); NC <- as.integer(Sys.getenv("SB_CORES","25"))
d0 <- prep(); DIR <- file.path(BASE,"4_结果/step3/比例自助"); dir.create(DIR, recursive=TRUE, showWarnings=FALSE)
des <- svydesign(id=~SDMVPSU, strata=~SDMVSTRA, weights=~WT_TOTAL, data=d0, nest=TRUE)
set.seed(20260928); W <- weights(as.svrepdesign(des, type="subbootstrap", replicates=B), type="analysis")
CVS <- c("+行为社会","全模型"); OK9 <- OUTC[setdiff(names(OUTC), c("EFP 1mm","EFP 2mm"))]
XC <- lapply(COV[CVS], function(v) model.matrix(as.formula(paste("~", paste(v, collapse="+"))), data=d0)[,-1,drop=FALSE])
## 规格清单与每个规格的设计矩阵行、原始标准误
specs <- list()
for (ex in EXPO){ v <- d0[[ex]]
  for (sc in SCALE){
    x <- if (sc=="每 SD") as.numeric(scale(v)) else if (sc=="每翻倍(log2)") log2(v) else {
      q <- cut(v, quantile(v,0:4/4,na.rm=TRUE), include.lowest=TRUE, labels=paste0("Q",1:4)); ifelse(q=="Q4",1, ifelse(q=="Q1",0,NA)) }
    for (sp in SAMP){ sk <- switch(sp, "全样本"=rep(TRUE,nrow(d0)), "WBC 3.5-11"= d0$WBC>=3.5 & d0$WBC<=11, "年龄>=40"= d0$Age>=40)
      for (ou in names(OK9)){ y <- d0[[OK9[[ou]]]]; keep <- sk & !is.na(x) & !is.na(y)
        for (cv in CVS) specs[[length(specs)+1]] <- list(exposure=ex, scale=sc, sample=sp, outcome=ou, covset=cv,
                                                          keep=which(keep), x=x, y=y) }}}}
cat("规格数", length(specs), "\n", flush=TRUE)
fit1 <- function(s, w){ i <- s$keep; X <- cbind(1, s$x[i], XC[[s$covset]][i,,drop=FALSE]); ww <- w[i]; ww <- ww/mean(ww[ww>0])
  f <- suppressWarnings(glm.fit(X, s$y[i], weights=ww, family=quasibinomial(), control=glm.control(maxit=50))); f$coefficients[2] }
## 原样本：svyglm 得设计基标准误与残差自由度（与 run_grid2 同）
O <- run_grid2(d0, miss="MICE 插补", wt="加权", outc=OK9) %>% filter(covset %in% CVS)
O <- O %>% mutate(logor=log(or), se=(log(hi)-log(or))/1.96)
key <- function(a) paste(a$exposure, a$scale, a$sample, a$outcome, a$covset, sep="|")
O <- O[match(sapply(specs, key), key(O)),]; stopifnot(!anyNA(O$logor))
e0 <- sapply(specs, fit1, w=d0$WT_TOTAL)
cat("b=0 与 svyglm 点估计最大差", signif(max(abs(e0-O$logor)),3), "\n", flush=TRUE)
## 残差自由度：由原 P 值与 t 统计量反解（逐规格；与 svyglm 的 df.residual 一致）
tcrit <- mapply(function(z, p) { if (p <= 1e-12 || p >= 0.9999) return(NA)
  df <- try(uniroot(function(k) 2*pt(-abs(z), k) - p, c(1, 1e4))$root, silent=TRUE); if (inherits(df,"try-error")) NA else qt(0.975, df) },
  O$logor/O$se, O$p)
tcrit[is.na(tcrit)] <- median(tcrit, na.rm=TRUE)
write.csv(data.frame(O[,c("exposure","scale","sample","outcome","covset","logor","se","p")], tcrit=tcrit, e0=e0), file.path(DIR,"原样本.csv"), row.names=FALSE)
one <- function(b){ fp <- file.path(DIR, sprintf("sb_%04d.csv", b)); if (file.exists(fp)) return(invisible(NULL))
  e <- sapply(specs, fit1, w=W[,b]); write.csv(data.frame(b=b, est=signif(e,6)), fp, row.names=FALSE); invisible(NULL) }
if (nzchar(Sys.getenv("SB_TEST"))) { t1 <- Sys.time(); one(1); print(Sys.time()-t1); quit() }
t0 <- Sys.time(); cat("比例自助 | B =", B, "| 核", NC, "\n", flush=TRUE)
invisible(mclapply(1:B, function(b) try(one(b), silent=TRUE), mc.cores=NC, mc.preschedule=FALSE))
cat("完成", length(list.files(DIR, pattern="^sb_")), "次，用时", round(difftime(Sys.time(),t0,units="mins"),1), "分钟\n")
