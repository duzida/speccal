## 第 3 步结果汇总：扩展置换（独立于筛选的 2,000 次）+ 多重校正；Rao–Wu 比例区间；sup-LR
suppressMessages({library(dplyr); library(tidyr)})
B0 <- Sys.getenv("NH_BASE", "."); R3 <- file.path(B0, "4_结果/4.6_第3步推断")
sat <- c("EFP 1mm","EFP 2mm"); adj <- c("+行为社会","全模型")
O <- read.csv(file.path(B0,"4_结果/4.1_规格网格/规格曲线_v2.csv")) %>% filter(missing=="MICE 插补", weight=="加权", covset %in% adj, !outcome %in% sat)
DOM <- read.csv(file.path(B0,"4_结果/4.2_置换与零分布/合理子集_FL推断_300.csv")) %>% select(exposure, dom)
To <- O %>% left_join(DOM, by="exposure") %>% group_by(exposure) %>% summarise(n=n(), T=sum(p<0.05 & sign(log(or))==dom), .groups="drop")
stopifnot(all(To$n==162))
## 1. 扩展置换
P <- bind_rows(lapply(list.files(file.path(R3,"置换扩展"), full.names=TRUE), read.csv))
Tn <- P %>% left_join(DOM, by="exposure") %>% group_by(rep, exposure) %>% summarise(T=sum(p<0.05 & sign(logor)==dom), n=n(), .groups="drop")
stopifnot(all(Tn$n==162)); nb <- length(unique(Tn$rep))
W <- Tn %>% select(rep, exposure, T) %>% pivot_wider(names_from=exposure, values_from=T)
ex <- To$exposure; obsT <- setNames(To$T, ex)
praw <- sapply(ex, function(e) (1+sum(W[[e]]>=obsT[e]))/(nb+1))
## 单步 minP（Westfall–Young）：每次置换内各暴露的置换 P，取最小值的分布
pstar <- sapply(ex, function(e){ x <- W[[e]]; sapply(x, function(v) mean(x>=v)) })
minp <- apply(pstar[, setdiff(ex, c("WBC","ANC","AGR")), drop=FALSE], 1, min)
comp <- setdiff(ex, c("WBC","ANC","AGR"))
p_minP <- sapply(comp, function(e) (1+sum(minp<=praw[e]))/(nb+1))
p_holm <- setNames(p.adjust(praw[comp], "holm"), comp)
old <- read.csv(file.path(B0,"4_结果/4.2_置换与零分布/合理子集_FL推断_300.csv")) %>% select(exposure, P300=FL_P)
res1 <- tibble(exposure=ex, T_obs=obsT, share=round(100*obsT/162,1), null_med=sapply(ex, function(e) median(W[[e]])),
               null95=sapply(ex, function(e) quantile(W[[e]],.95)), P_2000=signif(praw,3),
               P_holm7=signif(p_holm[ex],3), P_minP7=signif(p_minP[ex],3)) %>% left_join(old, by="exposure")
cat("独立置换次数:", nb, "\n"); print(as.data.frame(res1))
write.csv(res1, file.path(R3,"扩展置换_检验.csv"), row.names=FALSE)
## 2. Rao–Wu 比例区间
S0 <- read.csv(file.path(R3,"比例自助/原样本.csv")) %>% left_join(DOM, by="exposure")
fs <- list.files(file.path(R3,"比例自助"), pattern="^sb_", full.names=TRUE)
E <- sapply(fs, function(f) read.csv(f)$est)
sig <- (abs(E/S0$se) > S0$tcrit) & (sign(E)==S0$dom)
sh <- apply(sig, 2, function(s) tapply(s, S0$exposure, mean))
obs_sh <- tapply(S0$p<0.05 & sign(S0$logor)==S0$dom, S0$exposure, mean)
ne <- read.csv(file.path(B0,"4_结果/4.2_置换与零分布/有效独立设定数.csv"))
res2 <- tibble(exposure=rownames(sh), obs=round(100*obs_sh[rownames(sh)],1), boot_mean=round(100*rowMeans(sh),1),
               lo=round(100*apply(sh,1,quantile,.025),1), hi=round(100*apply(sh,1,quantile,.975),1), sd=round(100*apply(sh,1,sd),1))
cat("\nRao–Wu 自助次数:", ncol(E), "\n"); print(as.data.frame(res2)); cat("有效设定数文件列:", names(ne), "\n")
## 两两差：原始计数/AGR 与复合指数之差的区间
pair <- function(a,b) quantile(sh[a,]-sh[b,], c(.025,.975))
ref <- c("WBC","AGR","ANC"); cps <- c("SII","dNLR","NLR")
dd <- expand.grid(a=ref, b=cps, stringsAsFactors=FALSE); dd$lo <- NA; dd$hi <- NA
for (i in seq_len(nrow(dd))) { q <- pair(dd$a[i], dd$b[i]); dd$lo[i] <- round(100*q[1],1); dd$hi[i] <- round(100*q[2],1) }
print(dd)
write.csv(res2, file.path(R3,"比例_RaoWu区间.csv"), row.names=FALSE); write.csv(dd, file.path(R3,"比例差_RaoWu区间.csv"), row.names=FALSE)
## 3. sup-LR
L0 <- read.csv(file.path(R3,"supLR/观察值.csv")); LR <- bind_rows(lapply(list.files(file.path(R3,"supLR"), pattern="^lr_", full.names=TRUE), read.csv))
res3 <- L0 %>% rowwise() %>% mutate(B=sum(LR$exposure==exposure), P_supLR=(1+sum(LR$stat[LR$exposure==exposure]>=stat))/(B+1),
            null95=quantile(LR$stat[LR$exposure==exposure], .95)) %>% ungroup()
cat("\n"); print(as.data.frame(res3)); write.csv(res3, file.path(R3,"supLR_检验.csv"), row.names=FALSE)
