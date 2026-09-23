## 测试用的合成复杂抽样数据与轴定义（testthat 自动加载 helper-*.R）
make_toy <- function(n = 1500, seed = 1) {
  set.seed(seed)
  d <- data.frame(strat = rep(1:10, each = n / 10), psu = rep(1:2, n / 2), w = stats::runif(n, .5, 2),
                  Age = round(stats::runif(n, 30, 80)), Gender = factor(sample(c("F", "M"), n, TRUE)),
                  x1 = stats::rlnorm(n), x2 = stats::rlnorm(n))
  d$y_any <- stats::rbinom(n, 1, stats::plogis(-1 + .4 * log(d$x1) + .01 * d$Age))
  d$y_sev <- stats::rbinom(n, 1, stats::plogis(-2 + .3 * log(d$x1)))
  d
}
toy_axes <- function() spec_axes(
  outcome = c(any = "y_any", severe = "y_sev"),
  covariates = list(crude = character(0), demo = c("Age", "Gender")),
  coding = c("per_sd", "log2", "q4_vs_q1"),
  sample = list(all = TRUE, age40 = quote(Age >= 40)))

