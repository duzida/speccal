## 汇总一个零分布（一组 fl_*.csv）：各层级假阳性率、合理集合、各暴露置换 P
## 用法：Rscript summarise_null.R <零分布目录> <标签>
suppressMessages({library(dplyr); library(tidyr)})
B0 <- Sys.getenv("NH_BASE", ".")
a <- commandArgs(TRUE); DIR <- a[1]; LAB <- a[2]
N <- bind_rows(lapply(list.files(DIR, pattern="^fl_.*csv$", full.names=TRUE), read.csv))
nb <- length(unique(N$rep))
O <- read.csv(file.path(B0,"4_结果/4.1_规格网格/规格曲线_v2.csv")) %>% filter(missing=="MICE 插补", weight=="加权")
adj <- c("+行为社会","全模型"); sat <- c("EFP 1mm","EFP 2mm")
ci <- function(v) sprintf("%.1f (%.1f–%.1f)", 100*mean(v), 100*(mean(v)-1.96*sd(v)/sqrt(length(v))), 100*(mean(v)+1.96*sd(v)/sqrt(length(v))))
per_rep <- function(df, g) df %>% group_by(rep, .data[[g]]) %>% summarise(f=mean(p<0.05), .groups="drop")
fc <- per_rep(N, "covset") %>% group_by(covset) %>% summarise(FPR=100*mean(f), txt=ci(f), .groups="drop")
fo <- per_rep(filter(N, covset %in% adj), "outcome") %>% group_by(outcome) %>% summarise(FPR=100*mean(f), txt=ci(f), .groups="drop")
cat(sprintf("\n==== %s | %d 次置换 ====\n", LAB, nb)); print(as.data.frame(fc), row.names=FALSE)
print(as.data.frame(arrange(fo, desc(FPR))), row.names=FALSE)
bad_c <- fc$covset[fc$FPR > 10]; bad_o <- fo$outcome[fo$FPR > 10]
cat("10% 阈值剔除：协变量集 [", paste(bad_c, collapse=", "), "]；定义 [", paste(bad_o, collapse=", "), "]\n")
## 置换检验：主分析合理集合（固定 162 个）与本零分布自己导出的合理集合
dom <- O %>% group_by(exposure) %>% summarise(dom=ifelse(sum(p<0.05 & or>1) >= sum(p<0.05 & or<1), 1, -1), .groups="drop")
stat <- function(df, keep) df %>% filter(keep(covset, outcome)) %>% left_join(dom, by="exposure") %>%
  group_by(exposure) %>% summarise(T=sum(p<0.05 & sign(log(or))==dom), n=n(), .groups="drop")
k_main <- function(c,o) c %in% adj & !o %in% sat
k_own  <- function(c,o) !c %in% bad_c & !o %in% bad_o
test <- function(keep, nm){
  To <- stat(O, keep); Tn <- N %>% group_by(rep) %>% group_modify(~stat(.x, keep)) %>% ungroup()
  Tn %>% left_join(To %>% select(exposure, To=T, n_obs=n), by="exposure") %>% group_by(exposure) %>%
    summarise(n=first(n_obs), obs=100*first(To)/first(n_obs), null_med=100*median(T)/first(n_obs),
              null95=100*quantile(T,.95)/first(n_obs), P=(1+sum(T>=To))/(n()+1), .groups="drop") %>% mutate(set=nm)
}
R <- bind_rows(test(k_main, "main162"), test(k_own, "own"))
print(as.data.frame(R %>% mutate(across(c(obs,null_med,null95), ~round(.x,1)), P=round(P,3))), row.names=FALSE)
dir.create(file.path(B0,"4_结果/4.2_置换与零分布/残差化敏感性汇总"), showWarnings=FALSE)
write.csv(bind_rows(fc %>% mutate(level=covset, axis="covset"), fo %>% mutate(level=outcome, axis="outcome")) %>%
  select(axis, level, FPR, txt) %>% mutate(null=LAB, B=nb), file.path(B0,"4_结果/4.2_置换与零分布/残差化敏感性汇总", paste0("FPR_",LAB,".csv")), row.names=FALSE)
write.csv(R %>% mutate(null=LAB, B=nb), file.path(B0,"4_结果/4.2_置换与零分布/残差化敏感性汇总", paste0("test_",LAB,".csv")), row.names=FALSE)
