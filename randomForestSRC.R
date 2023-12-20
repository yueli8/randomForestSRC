#https://www.jingege.wang/2023/03/31/%e6%9c%ba%e5%99%a8%e5%ad%a6%e4%b9%a0%e4%b9%8b%e9
#%9a%8f%e6%9c%ba%e6%a3%ae%e6%9e%97%e7%94%9f%e5%ad%98%e5%88%86%e6%9e%90%ef%bc%88
#randomforestsrc%ef%bc%89/
#if (!require(randomForestSRC)) install.packages("randomForestSRC")
#if (!require(survival)) install.packages("survival")
library(randomForestSRC)
library(survival)
data(veteran, package = "randomForestSRC")
head(veteran)
##   trt celltype time status karno diagtime age prior
## 1   1        1   72      1    60        7  69     0
## 2   1        1  411      1    70        5  64    10
## 3   1        1  228      1    60        3  38     0
## 4   1        1  126      1    60        9  63    10
## 5   1        1  118      1    70       11  65    10
## 6   1        1   10      1    20        5  49     0

data(pbc, package = "randomForestSRC")
head(pbc)
##   days status treatment   age sex ascites hepatom spiders edema bili chol
## 1  400      1         1 21464   1       1       1       1   1.0 14.5  261
## 2 4500      0         1 20617   1       0       1       1   0.0  1.1  302
## 3 1012      1         1 25594   0       0       0       0   0.5  1.4  176
## 4 1925      1         1 19994   1       0       1       1   0.5  1.8  244
## 5 1504      0         2 13918   1       0       1       1   0.0  3.4  279
## 6 2503      1         2 24201   1       0       1       0   0.0  0.8  248
##   albumin copper    alk   sgot trig platelet prothrombin stage
## 1    2.60    156 1718.0 137.95  172      190        12.2     4
## 2    4.14     54 7394.8 113.52   88      221        10.6     3
## 3    3.48    210  516.0  96.10   55      151        12.0     4
## 4    2.54     64 6121.8  60.63   92      183        10.3     4
## 5    3.53    143  671.0 113.15   72      136        10.9     3
## 6    3.98     50  944.0  93.00   63       NA        11.0     3

data(wihs, package = "randomForestSRC")
head(wihs)
##   time status ageatfda idu black cd4nadir
## 1 0.02      2       48   0     1     6.95
## 2 0.02      2       35   1     1     2.51
## 3 0.02      2       28   0     1     0.18
## 4 0.02      2       46   1     0     4.65
## 5 0.02      2       31   0     1     0.08
## 6 0.02      2       45   1     1     2.05

#example1
v.obj <- rfsrc(Surv(time, status) ~ ., data = veteran, ntree = 100, nsplit = 10,
               na.action = "na.impute", tree.err = TRUE, importance = TRUE, block.size = 1)
###设置树的个数为3：'
## plot tree number 3
plot(get.tree(v.obj, 3))
###绘制训练森林的结果
## plot results of trained forest
plot(v.obj)
## plot survival curves for first 10 individuals -- direct way
matplot(v.obj$time.interest, 100 * t(v.obj$survival.oob[1:10, ]), xlab = "Time",
        ylab = "Survival", type = "l", lty = 1)
##使用函数plot.survival绘制前10个个体的生存曲线
## plot survival curves for first 10 individuals using function 'plot.survival'
plot.survival(v.obj, subset = 1:10)
## fast nodesize optimization for veteran data optimal nodesize in survival is
## larger than other families see the function 'tune' for more examples
tune.nodesize(Surv(time, status) ~ ., veteran)
tune.nodesize(Surv(time, status) ~ ., veteran)

#example2
vd <- veteran
vd$celltype = factor(vd$celltype)
vd$diagtime = factor(vd$diagtime)
vd.obj <- rfsrc(Surv(time, status) ~ ., vd, ntree = 100, nodesize = 5)
plot(get.tree(vd.obj, 3))

#原发性胆汁性肝硬化(PBC) example1
pbc.obj <- rfsrc(Surv(days, status) ~ ., pbc)
print(pbc.obj)

#example2
pbc.obj2 <- rfsrc(Surv(days, status) ~ ., pbc, nsplit = 10, na.action = "na.impute")
## same as above but iterate the missing data algorithm
pbc.obj3 <- rfsrc(Surv(days, status) ~ ., pbc, na.action = "na.impute", nimpute = 3)
## fast way to impute data (no inference is done) see impute for more details
pbc.imp <- impute(Surv(days, status) ~ ., pbc, splitrule = "random")

#3. 比较RF-SRC和Cox回归
#比较RF-SRC和Cox回归，说明性能的c指数和Brier评分措施，假设加载了“pec”和“survival”
require("survival")
require("pec")
require("prodlim")

## prediction function required for pec
predictSurvProb.rfsrc <- function(object, newdata, times, ...) {
  ptemp <- predict(object, newdata = newdata, ...)$survival
  pos <- sindex(jump.times = object$time.interest, eval.times = times)
  p <- cbind(1, ptemp)[, pos + 1]
  if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
    stop("Prediction failed")
  p
}

## data, formula specifications
data(pbc, package = "randomForestSRC")
pbc.na <- na.omit(pbc)  ##remove NA's
surv.f <- as.formula(Surv(days, status) ~ .)
pec.f <- as.formula(Hist(days, status) ~ 1)

## run cox/rfsrc models for illustration we use a small number of trees
cox.obj <- coxph(surv.f, data = pbc.na, x = TRUE)
rfsrc.obj <- rfsrc(surv.f, pbc.na, ntree = 150)

## compute bootstrap cross-validation estimate of expected Brier score see
## Mogensen, Ishwaran and Gerds (2012) Journal of Statistical Software
set.seed(17743)
prederror.pbc <- pec(list(cox.obj, rfsrc.obj), data = pbc.na, formula = pec.f, splitMethod = "bootcv", B = 50)
print(prederror.pbc)
plot(prederror.pbc)

## compute out-of-bag C-index for cox regression and compare to rfsrc
rfsrc.obj <- rfsrc(surv.f, pbc.na)
cat("out-of-bag Cox Analysis ...", "\n")
## out-of-bag Cox Analysis ...  
cox.err <- sapply(1:100, function(b) {
  if (b%%10 == 0)
    cat("cox bootstrap:", b, "\n")
  train <- sample(1:nrow(pbc.na), nrow(pbc.na), replace = TRUE)
  cox.obj <- tryCatch({
    coxph(surv.f, pbc.na[train, ])
  }, error = function(ex) {
    NULL
  })
  if (!is.null(cox.obj)) {
    get.cindex(pbc.na$days[-train], pbc.na$status[-train], predict(cox.obj, pbc.na[-train, ]))
  } else NA
}) 

#4.妇女跨机构艾滋病毒研究 生存竞争风险比例模型
wihs.obj <- rfsrc(Surv(time, status) ~ ., wihs, nsplit = 3, ntree = 100)
plot.competing.risk(wihs.obj)
cif <- wihs.obj$cif.oob
Time <- wihs.obj$time.interest
idu <- wihs$idu
cif.haart <- cbind(apply(cif[, , 1][idu == 0, ], 2, mean), apply(cif[, , 1][idu ==
                                                                              1, ], 2, mean))
cif.aids <- cbind(apply(cif[, , 2][idu == 0, ], 2, mean), apply(cif[, , 2][idu ==
                                                                             1, ], 2, mean))
matplot(Time, cbind(cif.haart, cif.aids), type = "l", lty = c(1, 2, 1, 2), col = c(4,
                                                                                   4, 2, 2), lwd = 3, ylab = "Cumulative Incidence")
legend("bottomright", legend = c("HAART (Non-IDU)", "HAART (IDU)", "AIDS (Non-IDU)",
                                 "AIDS (IDU)"), lty = c(1, 2, 1, 2), col = c(4, 4, 2, 2), lwd = 3, cex = 0.6)
