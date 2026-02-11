library(devtools)
clean_dll()
document()
load_all()
library(survival)
data(lung, package = "survival")

# 2) 清洗并构造 Surv 对象
# lung$status: 1=censored, 2=dead -> Surv 需要 0/1 的 event
df <- within(lung, {
  event <- as.integer(status == 2)
})

# 3) 选一些常用协变量（你也可以换成更多）
vars <- c("age", "sex", "ph.ecog", "wt.loss")
df2 <- df[, c("time", "event", vars)]
df2 <- na.omit(df2)

y <- Surv(time = df2$time, event = df2$event)

# 4) glmnet 需要矩阵 X（用 model.matrix 自动处理因子/截距）
X <- model.matrix(~ age + sex + ph.ecog + wt.loss, data = df2)[, -1, drop = FALSE]

# 5) 交叉验证选择 lambda（默认 10-fold）
set.seed(123)
sunifit = sulnet2D(X,y,method = "septhresh", family = "cox",eps = 1e-07, maxit = 1000000,alpha = NULL,nlambda = 100)
coxfit = sulnet2D(X,y,method = "suni_2", family = "cox",eps = 1e-07, maxit = 1000000)
sunifit$beta

test = sulnet2D(X,y[,1],method = "septhresh")
