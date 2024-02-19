## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, message = FALSE)
library(WLogit)
library(tibble)
library(ggplot2)
set.seed(123456)

## ----install pacakge, eval=FALSE----------------------------------------------
#  install.packages("WLogit", repos = "http://cran.us.r-project.org")

## ----generate Sigma-----------------------------------------------------------
p <- 500 # number of variables 
d <- 10 # number of actives
n <- 100 # number of samples
actives <- c(1:d)
nonacts <- c(1:p)[-actives]
Sigma <- matrix(0, p, p)
Sigma[actives, actives] <- 0.3
Sigma[-actives, actives] <- 0.5
Sigma[actives, -actives] <- 0.5
Sigma[-actives, -actives] <- 0.7
diag(Sigma) <- rep(1,p)

## ----X, eval=FALSE------------------------------------------------------------
#  X <- MASS::mvrnorm(n = n, mu=rep(0,p), Sigma, tol = 1e-6, empirical = FALSE)
#  beta <- rep(0,p)
#  beta[actives] <- 1
#  pr <- CalculPx(X,beta=beta)
#  y <- rbinom(n,1,pr)

## ---- echo=FALSE--------------------------------------------------------------
data(X)
data(y)
data(beta)

## ----load WLogit, eval=FALSE--------------------------------------------------
#  library(WLogit)

## ----WLogit model, eval=FALSE-------------------------------------------------
mod <- WhiteningLogit(X = X, y = y)

## ---- echo=FALSE--------------------------------------------------------------
data(test)
mod <- test

## ----beta---------------------------------------------------------------------
beta_min <- mod$beta.min
head(beta_min)

## ----variable selection,fig.width=4,fig.height=3------------------------------
beta_min <- mod$beta.min
df_beta <- data.frame(beta_est=beta_min, Status = ifelse(beta==0, "non-active", "active"))
df_plot <- df_beta[which(beta_min!=0), ]
df_plot$index <- which(beta_min!=0)
ggplot2::ggplot(data=df_plot, mapping=aes(y=beta_est, x=index, color=Status))+geom_point()+
  theme_bw()+ylab("Estimated coefficients")+xlab("Indices of selected variables")

a <- round(CalculPx(X+1, beta_min))
table(a,y)

## ----lasso--------------------------------------------------------------------
library(glmnet)
cvfit = cv.glmnet(X, y, family = "binomial", type.measure = "class", intercept=FALSE)

## ----res lasso----------------------------------------------------------------
beta_lasso <- coef(cvfit, s = "lambda.min")
head(beta_lasso)

## ----lasso selection,fig.width=4,fig.height=3---------------------------------
beta_lasso <- as.vector(beta_lasso)[-1]
df_beta <- data.frame(beta_est=beta_lasso, Status = ifelse(beta==0, "non-active", "active"))
df_plot <- df_beta[which(beta_lasso!=0), ]
df_plot$index <- which(beta_lasso!=0)
ggplot2::ggplot(data=df_plot, mapping=aes(y=beta_est, x=index, color=Status))+geom_point()+
  theme_bw()+ylab("Estimated coefficients by glmnet")+xlab("Indices of selected variables")
