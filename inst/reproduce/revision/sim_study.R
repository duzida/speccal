## 模拟研究（R3-M2 / B1）：已知真值下，零分布判据与置换检验的表现
## 以 NHANES 分析样本的真实协变量为底（n=10,108 或随机抽 3,000），合成暴露与结局：
##   x   = Zγ + e（γ 取自真实 log 中性粒细胞对全协变量的回归，e 为正态，方差同真实残差）
##   y_k = Bernoulli(plogis(α_k + Zθ + β·x_sd))，θ 取自真实 CDC/AAP 牙周炎全模型；α_k 使患病率约 42% 与 8%
## 情景：measured（混杂全部已测且线性）| nonlinear（年龄对 x 和 y 均有二次效应，分析模型为线性）|
##       unmeasured（未测混杂 U 同时影响 x 与 y）| collider（全模型含一个受 x 与 y 共同影响的变量）
## 每次重复：4 协变量集 × 2 编码（每 SD、Q4 vs Q1）× 2 结局 = 16 个设定；
##   零分布 = x 对（该情景下全模型的）协变量线性残差化后置换残差，B=99；
##   假阳性率 > 10% 的协变量集判为不合理；统计量 = 合理设定中在众数方向显著的个数；P =(1+#≥)/(B+1)
BASE <- Sys.getenv("NH_BASE"); source(file.path(BASE,"3_分析脚本/spec_core.R"))
suppressMessages(library(parallel))
NC <- as.integer(Sys.getenv("SIM_CORES","12")); NPERM <- as.integer(Sys.getenv("SIM_B","99"))
d0 <- prep(); DIR <- file.path(BASE,"4_结果/模拟研究v2"); dir.create(DIR, recursive=TRUE, showWarnings=FALSE)
mm <- function(v, d) if (length(v)) model.matrix(as.formula(paste("~", paste(v, collapse="+"))), data=d)[,-1,drop=FALSE] else NULL
CS <- list(crude=character(0), demo=COV[["+人口学"]], socio=COV[["+行为社会"]], full=COV[["全模型"]])
Zall <- lapply(CS, mm, d=d0)
Zf <- cbind(1, Zall$full)
xr <- log(d0$ANC); gam <- qr.coef(qr(Zf), xr); sig <- sd(xr - Zf %*% gam)
th <- coef(glm(d0$y_cdc_any ~ Zall$full, family=binomial()))[-1]
agez <- as.numeric(scale(d0$Age))
fitp <- function(X, y){ f <- suppressWarnings(glm.fit(X, y, family=binomial())); if (!f$converged) return(c(NA,NA))
  R <- f$qr; V <- chol2inv(R$qr[1:R$rank,1:R$rank,drop=FALSE]); j <- which(R$pivot==2)
  b <- f$coefficients[2]; se <- sqrt(V[j,j]); c(b, 2*pnorm(-abs(b/se))) }
codes <- function(x){ q <- quantile(x, c(.25,.5,.75))
  list(sd = cbind(as.numeric(scale(x)), TRUE), q41 = cbind(as.numeric(x>q[3]), x<=q[1] | x>q[3])) }
grid_fit <- function(x, Y, Zs){          # 返回 16 行：cov, code, out, b, p
  cx <- codes(x); res <- list()
  for (cv in names(Zs)) for (cd in names(cx)) for (k in 1:2){
    keep <- cx[[cd]][,2]==1; X <- cbind(1, cx[[cd]][keep,1], Zs[[cv]][keep,,drop=FALSE])
    r <- fitp(X, Y[keep,k]); res[[length(res)+1]] <- data.frame(cov=cv, code=cd, out=k, b=r[1], p=r[2]) }
  do.call(rbind, res) }
one_rep <- function(sc, n, OR, r){
  set.seed(1e5*match(sc, SCN) + 1e3*(n==3000) + 1e4*round(10*OR) + r)
  i <- if (n < nrow(d0)) sample(nrow(d0), n) else seq_len(nrow(d0))
  Zs <- lapply(Zall, function(Z) if (is.null(Z)) matrix(0, length(i), 0) else Z[i,,drop=FALSE])
  eta_x <- drop(Zf[i,] %*% gam); lin <- drop(Zs$full %*% th); a2 <- agez[i]^2 - 1
  U <- rnorm(length(i))
  x <- eta_x + rnorm(length(i), 0, sig)
  if (sc=="nonlinear") { x <- x + 0.4*sig*a2; lin <- lin + 0.5*a2 }
  if (sc=="unmeasured"){ x <- x + 0.5*sig*U;  lin <- lin + 0.5*U }
  bx <- log(OR) * (x - mean(x))/sd(x)
  al <- sapply(c(.42,.08), function(p) uniroot(function(a) mean(plogis(a + lin + bx)) - p, c(-15,15))$root)
  ystar <- lin + bx
  Y <- sapply(al, function(a) rbinom(length(i), 1, plogis(a + ystar)))
  if (sc=="collider"){                     # 结局由潜变量阈值产生；对撞变量 C 同时受 x 和结局潜变量（含其随机项）影响
    eps <- rlogis(length(i)); Y <- sapply(al, function(a) as.integer(a + ystar + eps > 0))
    C <- as.numeric(scale(x)) + as.numeric(scale(ystar + eps)) + rnorm(length(i))
    Zs$full <- cbind(Zs$full, C) }
  Zr <- cbind(1, Zs$full); fx <- drop(Zr %*% qr.coef(qr(Zr), x)); ex <- x - fx
  obs <- grid_fit(x, Y, Zs)
  nul <- lapply(1:NPERM, function(b) grid_fit(fx + sample(ex), Y, Zs))
  fpr <- sapply(names(CS), function(cv) mean(unlist(lapply(nul, function(g) g$p[g$cov==cv] < 0.05)), na.rm=TRUE))
  okc <- names(CS)[fpr <= 0.10]
  dom <- sign(sum(sign(obs$b[obs$cov %in% okc]), na.rm=TRUE)); if (dom==0) dom <- 1   # 众数方向取自合理设定
  stat <- function(g, sel) sum(g$cov %in% sel & g$p < 0.05 & sign(g$b)==dom, na.rm=TRUE)
  To <- stat(obs, okc); Tn <- sapply(nul, stat, sel=okc)
  Ta <- stat(obs, names(CS)); Tna <- sapply(nul, stat, sel=names(CS))
  data.frame(scenario=sc, n=length(i), OR=OR, rep=r, fpr_crude=fpr[1], fpr_demo=fpr[2], fpr_socio=fpr[3], fpr_full=fpr[4],
             defensible=paste(okc, collapse="+"), n_def=sum(obs$cov %in% okc), T_def=To,
             P_def=(1+sum(Tn>=To))/(NPERM+1), T_all=Ta, P_all=(1+sum(Tna>=Ta))/(NPERM+1),
             share_full=mean(obs$p[obs$cov=="full"]<0.05 & sign(obs$b[obs$cov=="full"])==dom, na.rm=TRUE))
}
SCN <- c("measured","nonlinear","unmeasured","collider")
COND <- rbind(expand.grid(scenario=SCN, n=10108, OR=1, R=300, stringsAsFactors=FALSE),
              data.frame(scenario="measured", n=3000, OR=1, R=300),
              data.frame(scenario="measured", n=c(10108,10108,3000), OR=c(1.1,1.2,1.2), R=150))
if (nzchar(Sys.getenv("SIM_ONLY"))) COND <- subset(COND, scenario==Sys.getenv("SIM_ONLY") & OR==1)
jobs <- do.call(rbind, lapply(seq_len(nrow(COND)), function(k) data.frame(COND[k, 1:3], r=1:COND$R[k])))
if (nzchar(Sys.getenv("SIM_TEST"))) { t1 <- Sys.time(); for (sc in SCN) print(one_rep(sc, 10108, 1, 1)); print(Sys.time()-t1); quit() }
cat("模拟研究 |", nrow(COND), "个条件 |", nrow(jobs), "次重复 | 核", NC, "\n", flush=TRUE)
t0 <- Sys.time()
chunks <- split(seq_len(nrow(jobs)), ceiling(seq_len(nrow(jobs))/20))
invisible(mclapply(seq_along(chunks), function(ci){
  fp <- file.path(DIR, sprintf("sim_%s%04d.csv", Sys.getenv("SIM_ONLY"), ci)); if (file.exists(fp)) return(NULL)
  out <- lapply(chunks[[ci]], function(j) try(with(jobs[j,], one_rep(scenario, n, OR, r)), silent=TRUE))
  out <- do.call(rbind, Filter(is.data.frame, out))
  write.csv(out, fp, row.names=FALSE); NULL }, mc.cores=NC, mc.preschedule=FALSE))
cat("完成，用时", round(difftime(Sys.time(),t0,units="mins"),1), "分钟\n")
