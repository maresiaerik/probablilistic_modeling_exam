library("igraph")
library("gRbase")
library("mgm")
library("RBGL")
library("gRim")
library("tidyr")
library("ggplot2")
library("factoextra")
library("dplyr")
library("corrplot")
library("ggfortify")
library(glasso)
library(caret)
install.packages("moments")
library(moments)
library("Rfast")
install.packages("snpStats")
library("snpStats")
library(mgcViz)
library(stringr)
library(CVglasso)
library(SIN)
library(hash)
library(mlbench)
remotes::install_github("vqv/ggbiplot", force = TRUE)




data = read.csv("./data/who.csv", sep=";")
data = data[data$year == 2016,]
names(data)



# categorical, or NA columns
columns_to_drop = c("year", "country","country_code", "doctors", "region", "une_poverty", "une_edu_spend", "une_literacy", "une_school", "gne_capita", "gni_capita", "hospitals", "life_exp60", "une_life", "adult_mortality")
data = data[, !(names(data) %in% columns_to_drop)]

summary(data)
data = drop_na(data)

median(data$life_expect)
mean(data$life_expect)
boxplot(data$life_expect, main="life_expect")
plot(density(data$life_expect))


corrplot(cor(data))

df1 = data

correlationMatrix <- cor(data[,2:length(names(data))])
cor1 = correlationMatrix
# summarize the correlation matrix
corrplot(correlationMatrix)

highlyCorrelated <- findCorrelation(correlationMatrix, cutoff=0.9, names = T)
cols_to_drop = !(names(data) %in% highlyCorrelated)
data = data[cols_to_drop]
cor2 = cor(data[,2:length(names(data))])

par(mfrow=c(1,2))
corrplot(cor1)
mtext("All variables", at=8, line=-3, cex=1.2)
corrplot(cor2)
mtext("Variables after removing highly correlated", at=8, line=-3, cex=1.2)

variables = names(data)
var_df = data.frame(variables)
var_df$Type = "continuous"
colnames(var_df) = c("variable", "continuous")



dim(data)

S.data = cov.wt(data, method = "ML", cor=T)
table(S.data$cov)
corrplot(S.data$cor)

PC.data = cov2cor(S.data$cov)
round(100*PC.data)

mod.sat = cmod(~.^., data)
mod.sat

plot(mod.sat)
round(100*PC.data)



range = seq(0.01, 0.8, by=0.04)
range
allowed = c(0.01, 0.05, 0.17, 0.25)
loglikes = list()

par(mfcol=c(2,2))

for (rho in range) {
  res.lasso = glasso(PC.data, rho=rho)
  AM <- res.lasso$wi!=0 # estimate K matrix
  
  diag(AM) <-FALSE
  g.lasso <- as(AM,"graphNEL")
  nodes(g.lasso) <- names(data)
  glasso.body <- cmod(edgeList(g.lasso),data=data)
  plot(as(glasso.body,"igraph"))
  glasso.data = cmod(edgeList(g.lasso), data=data)
  rho
  if (rho %in% allowed) {
    s = str_interp("rho=${rho}")
    plot(as(glasso.data,"igraph"), main=s)
  }
  key = as.character(rho)
  loglikes[[key]] = -2*glasso.body$fitinfo$logL
}
names(unlist(loglikes))
unlist(loglikes)
plot(unlist(loglikes), xlab=names(loglikes),  type = "l", xaxt = "n")

max(unlist(loglikes))

l
lasso.model = glasso(PC.data, rho=0.17)
AM <- lasso.model$wi!=0 # estimate K matrix

diag(AM) <-FALSE
g.lasso <- as(AM,"graphNEL")
nodes(g.lasso) <- names(data)
lasso.graph <- cmod(edgeList(g.lasso),data=data)
plot(as(lasso.graph,"igraph"))
lasso.graph



par(mfcol=c(2,2))
plot(as(gl01, "igraph"), main="rho = 0.05")
plot(as(gl0, "igraph"), main="rho = 0.1")
plot(as(gl1, "igraph"), main="rho = 0.17")
plot(as(gl2, "igraph"), main="rho = 0.25")



psin.body = sinUG(S.data$cov, n=nrow(data))
plotUGpvalues(psin.body)
gsin.body03 = as(getgraph(psin.body, 0.3), "graphNEL")
gsin.body02 = as(getgraph(psin.body, 0.2), "graphNEL")
par(mfcol=c(1,2))
plot(as(gsin.body03, "igraph"), main="alpha = 0.3")
plot(as(gsin.body02, "igraph"), main="alpha = 0.2")

p.sin.mod = cmod(edgeList(gsin.body03), data=data)
p.sin.mod
key = as.character(alpha)
loglikes[[key]] = -2*p.sin.mod$fitinfo$logL
plot(as(p.sin.mod, "igraph"))
loglikes
names(unlist(loglikes))
plot(unlist(loglikes), xlab=names(loglikes),  type = "l", xaxt = "n")

max(unlist(loglikes))





psin = as(p.sin.mod, "graphNEL")
gl = as(lasso.graph, "graphNEL")

commonedges = graph::intersection(psin, gl)
plot(as(commonedges, "igraph"))
common.mod = cmod(edgeList(commonedges), data=data)
common.mod




## Check single hypotheses

response = "life_expect"


ind.y = names(data)[1]
ind.x = "age1.4mort"
excluded = c(response, ind.x)

verteces_to_exclude = !(names(data) %in% excluded)
verteces = as.list(names(data))
rest_of_v = verteces[verteces_to_exclude]

independence_check = paste("~", paste(excluded, collapse="+"))

v.formula = as.formula(paste(independence_check, "+", paste(rest_of_v, collapse = "+")))

a =ciTest_mvn(list(cov=S.data $cov, n.obs=nrow(data)), set=v.formula, statistic = "F")


p.sin.mod$fitinfo$dev
lasso.graph $fitinfo$dev

lrt <- p.sin.mod$fitinfo$dev- lasso.graph $fitinfo$dev
dfdiff <- p.sin.mod$fitinfo$dimension[4]-lasso.graph$fitinfo$dimension[4]
names(dfdiff) <- NULL

round(pchisq(28.21435, 6, lower.tail = F), 16) 
