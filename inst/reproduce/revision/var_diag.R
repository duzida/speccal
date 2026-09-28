## 诊断：全模型在零分布下假阳性率 6.7% 而非 5% 的来源（R2-M1）
## 零分布构造与主分析相同（全协变量线性残差化、层内置换）；只拟合全模型、9 套非饱和定义、log2、全样本 = 每次 90 个设定
## 同一数据上比较 5 种推断：线性化 SE + t(残差自由度)【主分析】、线性化 SE + t(设计自由度)、线性化 SE + 正态、
##   JKn 刀切重复权重 SE + t(设计自由度)、未加权 glm 模型 SE
BASE <- Sys.getenv("NH_BASE"); source(file.path(BASE,"3_分析脚本/spec_core2.R"))
suppressMessages(library(parallel))
B <- as.integer(Sys.getenv("VD_B","100")); NC <- as.integer(Sys.getenv("VD_CORES","8"))
d0 <- prep(); DIR <- file.path(BASE,"4_结果/方差诊断"); dir.create(DIR, recursive=TRUE, showWarnings=FALSE)
COVF <- COV[["全模型"]]; Xc <- model.matrix(as.formula(paste("~", paste(COVF, collapse="+"))), data=d0)
L <- sapply(EXPO, function(ex) log(d0[[ex]])); fitL <- Xc %*% qr.coef(qr(Xc), L); resL <- L - fitL
st <- split(seq_len(nrow(d0)), d0$SDMVSTRA)
OK9 <- setdiff(names(OUTC), c("EFP 1mm","EFP 2mm"))
one <- function(b){
  fp <- file.path(DIR, sprintf("vd_%03d.csv", b)); if (file.exists(fp)) return(invisible(NULL))
  set.seed(500000+b)                       # 与主分析 FL 零分布同一种子 → 前 100 次与主分析逐次对应
  idx <- seq_len(nrow(d0)); for (s in st) idx[s] <- if (length(s)>1) sample(s) else s
  d <- d0; Lp <- fitL + resL[idx,,drop=FALSE]; for (j in seq_along(EXPO)) d[[EXPO[j]]] <- exp(Lp[,j])
  des <- svydesign(id=~SDMVPSU, strata=~SDMVSTRA, weights=~WT_TOTAL, data=d, nest=TRUE)
  rep <- as.svrepdesign(des, type="JKn")
  out <- list()
  for (ex in EXPO){ x <- log2(d[[ex]]); des2 <- update(des, ..x=x); rep2 <- update(rep, ..x=x)
    for (ou in OK9){ y <- OUTC[[ou]]
      f <- as.formula(paste(y, "~ ..x +", paste(COVF, collapse="+")))
      m  <- svyglm(f, design=des2, family=quasibinomial()); s <- summary(m)$coefficients["..x",]
      dfr <- m$df.residual; dfd <- degf(des2)
      mr <- svyglm(f, design=rep2, family=quasibinomial()); sr <- summary(mr)$coefficients["..x",]
      dd <- d; dd$..x <- x; mu <- glm(f, data=dd, family=binomial()); su <- summary(mu)$coefficients["..x",]
      z <- s[1]/s[2]
      out[[length(out)+1]] <- data.frame(rep=b, exposure=ex, outcome=ou, est=s[1], se_lin=s[2], se_jk=sr[2],
        df_resid=dfr, df_design=dfd,
        p_lin_tres=2*pt(-abs(z), dfr), p_lin_tdes=2*pt(-abs(z), dfd), p_lin_norm=2*pnorm(-abs(z)),
        p_jk_tdes=2*pt(-abs(sr[1]/sr[2]), degf(rep2)), p_unw=su[4], row.names=NULL)
    }}
  write.csv(do.call(rbind,out), fp, row.names=FALSE); invisible(NULL)
}
t0 <- Sys.time(); cat("方差诊断 | B =", B, "| 核", NC, "\n", flush=TRUE)
invisible(mclapply(1:B, function(b) try(one(b), silent=FALSE), mc.cores=NC, mc.preschedule=FALSE))
cat("完成", length(list.files(DIR)), "次，用时", round(difftime(Sys.time(),t0,units="mins"),1), "分钟\n")
