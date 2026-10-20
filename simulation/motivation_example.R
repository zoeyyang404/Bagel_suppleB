library(bagelR)
pdf(file='simulation/step_2000_motivation.pdf',width=5,height=4)
mean_feature <- function(time, changepoint) {
  if (changepoint == 0) return(matrix(1, 1, 1))
  M <- matrix(0, 2, 1); M[1, 1] <- 1
  if (time > changepoint) M[2, 1] <- 1
  M
}
make_mean_prior <- function(pre_change_mean, post_change_mean, variance) {
  force(pre_change_mean); force(post_change_mean); force(variance)
  function(time) {
    if (time == 1)
      return(list(mu = matrix(pre_change_mean, 1, 1), sigma = matrix(variance, 1, 1)))
    list(mu    = matrix(c(pre_change_mean, post_change_mean - pre_change_mean), 2, 1),
         sigma = variance * matrix(c(1, -1, -1, 2), 2, 2, byrow = TRUE))
  }
}
mean_transform <- function(changepoint) matrix(c(1, 1), 1, 2, byrow = TRUE)
prior <- make_mean_prior(0, 0, 10)
prune_recover <- function(fit) {
    TT=3000
    recovered_weights <- numeric(TT)
    recovered_weights[taus(fit)+1] <- active_weights(fit)
    return(recovered_weights)
  }
set.seed(12345)
x=c(rnorm(3000))
x[1501:3000]=x[1501:3000]+0.2

res11=bagel_uv_approximate(x[1:1600],p=1,p0 =0.5,prior=prior,nu=5,iota=5,
        transform=mean_transform,feature_vector=mean_feature,threshold=Inf,max_particles=50000)
res21=bagel_uv_approximate(x[1:1600],p=1,p0 =0.5,prior=prior,nu=5,iota=5,
        transform=mean_transform,feature_vector=mean_feature,threshold=Inf,max_particles=50)
res31=bagel_uv_approximate(x[1:1600],p=1,p0 =0.5,prior=prior,nu=5,iota=5,
        transform=mean_transform,feature_vector=mean_feature,threshold=Inf,max_particles=50,pruning='LP-P')
res41=bagel_uv_approximate(x[1:1600],p=1,p0 =0.5,prior=prior,nu=5,iota=5,
        transform=mean_transform,feature_vector=mean_feature,threshold=Inf,max_particles=50,pruning='LP-M')
res_cdf21=weights(res31)
res_cdf41=weights(res41)
res12=bagel_uv_approximate(x[1:2000],p=1,p0 =0.5,prior=prior,nu=5,iota=5,
        transform=mean_transform,feature_vector=mean_feature,threshold=Inf,max_particles=50000)
res22=bagel_uv_approximate(x[1:2000],p=1,p0 =0.5,prior=prior,nu=5,iota=5,
        transform=mean_transform,feature_vector=mean_feature,threshold=Inf,max_particles=50)
res32=bagel_uv_approximate(x[1:2000],p=1,p0 =0.5,prior=prior,nu=5,iota=5,
        transform=mean_transform,feature_vector=mean_feature,threshold=Inf,max_particles=50,pruning='LP-P')
res42=bagel_uv_approximate(x[1:2000],p=1,p0 =0.5,prior=prior,nu=5,iota=5,
        transform=mean_transform,feature_vector=mean_feature,threshold=Inf,max_particles=50,pruning='LP-M')    
res_cdf22=weights(res32)
res_cdf42=weights(res42)
loc <- function(res) taus(res) + 1                 # 1-based candidate locations
map_loc <- which.max(weights(res12)[-1]) + 1       # exact MAP location at step 2000

post<-weights(res12)[-1]/sum(weights(res12)[-1])
ord <- order(post, decreasing = TRUE)
cut_idx <- which(cumsum(post[ord]) >= 0.9)[1]
selected <- ord[seq_len(cut_idx)]

par(mar = c(4, 4, 1, 1))
plot(x,col='gray30',ylab='observations',xlab='time step',cex.axis=1.4,cex.lab=1.4)
abline(v=1500,lwd=5,lty=2,col="red")
points(1:1500,rep(0,1500),type='l',lwd=5,col="red")
points(1501:3000,rep(0.25,1500),type='l',lwd=5,col="red")

plot(loc(res11)[-1],weights(res11)[-1],type='l',xlab='changepoint candidates',ylab='posterior probability',
cex.axis=1.4,cex.lab=1.4,lwd=4,col='#D6CFC4',xaxt = 'n',xlim=c(1,1600))  # Turn off automatic x-axis
axis(side = 1, at =c(1,500,1000,1500),cex.axis=1.4,cex.lab=1.4)
points(loc(res31)[-1],active_weights(res11)[loc(res31)[-1]],col='black',pch=19)
#points(map_loc,weights(res11)[map_loc],col='red',pch=8,cex=2)

plot(loc(res12)[-1],weights(res12)[-1],type='l',xlab='changepoint candidates',ylab='posterior probability',
cex.axis=1.4,cex.lab=1.4,lwd=4,col='#D6CFC4',xaxt='n')
axis(side=1,at=c(1,500,1000,1500,2000),cex.axis=1.4,cex.lab=1.4)
points(loc(res32)[-1],active_weights(res12)[loc(res32)[-1]],col='black',pch=19)

plot(cumsum(weights(res11)[-1]/sum(weights(res11)[-1])),type='l',lwd=6,ylab='CDF',xlab='changepoint candidates',col='#D6CFC4',cex.axis=1.4,cex.lab=1.4,xaxt='n')
axis(side = 1, at =c(1,500,1000,1500),cex.axis=1.4,cex.lab=1.4)
points(cumsum(weights(res21)[-1]/sum(weights(res21)[-1])),col='#466CA6',type='l',lty=2,lwd=4)
points(cumsum(prune_recover(res31)[-1]/sum(prune_recover(res31)[-1])),col='black',type='l',lty=2,lwd=4)
points(cumsum(weights(res41)[-1]/sum(weights(res41)[-1])),col='orange',type='l',lty=2,lwd=4)
plot(cumsum(weights(res12)[-1]/sum(weights(res12)[-1])),type='l',lwd=6,ylab='CDF',xlab='changepoint candidates',col='#D6CFC4',cex.axis=1.4,cex.lab=1.4,xaxt='n')
axis(side=1,at=c(1,500,1000,1500,2000),cex.axis=1.4,cex.lab=1.4)
points(cumsum(weights(res22)[-1]/sum(weights(res22)[-1])),col='#466CA6',type='l',lty=2,lwd=4)
points(cumsum(weights(res42)[-1]/sum(weights(res42)[-1])),col='orange',type='l',lty=2,lwd=4)
points(cumsum(prune_recover(res32)[-1]/sum(prune_recover(res32)[-1])),col='black',type='l',lty=2,lwd=4)
dev.off()
