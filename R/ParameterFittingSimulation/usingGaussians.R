library(fitness)
library(ks)
library(reshape2)
library(dplyr)
library(ggplot2)

# For the entire set (it doesn't matter, so we choose something)
len <- 100
x <- seq(-1,1,length.out=len)

createDe <- function(x, means, sds) {
  rawde <- rnorm(length(x), mean=means, sd=sds)
  de <- predict(kde(rawde), x=x)
  return(de)
}

nSimPars <- 20
simPars <- expand.grid(
  means=0,#seq(-1,1,length.out=nSimPars),
  sds=seq(0.5,2,length.out=nSimPars)
)
nSolvePars <- 20
solvePars <- expand.grid(
  eta=10^(seq(-0.5,5,length.out=nSolvePars)),
  mu=seq(0.01,2,length.out=nSolvePars)
)

des <- as.data.frame(matrix(NA, nrow=0, ncol=4))
colnames(des) <- c("x","de","means","sds")
ws <- as.data.frame(matrix(NA, nrow=0, ncol=6))
colnames(ws) <- c("x","w","means","sds","eta","mu")
simulations <- nrow(simPars)
solves <- nrow(solvePars)
for (j in 1:simulations) {
  means <- simPars[j,"means"]
  sds <- simPars[j,"sds"]
  de <- createDe(
    x,
    means,
    sds
  )
  # des <- rbind(
  #   des,
  #   data.frame(
  #     x=x,
  #     de=de/mean(de),
  #     means=means,
  #     sds=sds
  #   )
  # )
  W <- as.data.frame(matrix(NA, nrow=0, ncol=3))
  colnames(W) <- c("w","eta","mu")
  for (k in 1:solves) {
    eta <- solvePars[k,"eta"]
    mu <- solvePars[k,"mu"]
    w <- landscape(x, de, eta=eta, mu=mu, sigma=2, kde=F, verbose=F, scale="mean")$w
    realw <- exp(0.5*(1/(mu^2-sds^2)+1/sds^2)*(x-means)^2)
    W <- rbind(
      W,
      data.frame(w=w, realw=realw/mean(realw), eta=eta, mu=mu)
    )
  }
  ws <- rbind(
    ws,
    cbind(
      data.frame(
        x=x,
        means=means,
        sds=sds
      ),
      W
    )
  )
}

# des <- des %>%
#   mutate(
#     meanssd=paste0("mean=",means,";sd=",sds)
#   )
ws <- ws %>%
  mutate(
    meanssd=paste0("mean=",means,";sd=",sds),
    etamu=paste0("eta=",round(eta,2),";mu=",round(mu,2))
  )

# mws <- ws %>%
#   melt(id=c("x","meanssd","etamu"))
# mdes <- des %>%
#   melt(id=c("x","meansd"))


ggplot(ws, aes(x=x,color=as.factor(eta)))+
  geom_line(aes(y=w))+
  geom_line(aes(y=realw),linetype="dashed")+
  facet_grid(round(mu,2)~round(sds,2),scales="free")+
  theme_classic()+theme(legend.position="none")


sws <- ws %>%
  group_by(
    paste0(mu,";",eta,";",means,";",sds)
  ) %>%
  summarize(
    mu=mu,
    eta=eta,
    means=means,
    sds=sds,
    error=sum(w-realw)^2
  )

ggplot(sws, aes(x=eta, y=error))+
  geom_line()+
  facet_grid(round(sds,2)~round(mu,2))+
  theme_classic()







