rm(list=ls())
library(inline)
library(Rcpp)
library(RcppArmadillo)
library(dplyr)
library(tidyverse)
library(ggpubr)
library(gridExtra)
library(grid)
library(lattice)
library(plotly)
library(paletteer)
theme_set(theme_pubr())


# setwd("/Users/harmanjaggi/Documents/Research/dispersal/migration code")

sigsq <- c(0, 0.4, 0.8, 1.2)
# eps <- c(1e-10,1e-5,5e-5,1e-4,2e-4,4e-4,6e-4,8e-4,1e-3,2e-3,4e-3,6e-3,8e-3,.01)
# sigsq <- seq(0.1, 0.7, length=4)
sigma <- sqrt(sigsq)

# eps <- c(0,1e-5,5e-5,1e-4,2e-4,4e-4,6e-4,8e-4,1e-3,2e-3,4e-3,6e-3,8e-3,.01)

eps1 <- seq(0, 0.005, length=10)
eps1 <- round(eps1, 4)
eps2 <- c(1e-10, 1e-7, 2e-7, 3e-7, 4e-7, 5e-7,1e-6, 2e-6, 3e-6, 4e-6, 5e-6, 1e-6, 2e-5, 3e-5, 4e-5, 5e-5, 0.0001, 0.0002, 0.0003, 0.0004, 0.0005)
eps3 <- c(eps1, eps2)
eps <- eps3[order(eps3)]
length(eps)

u_bar_func_n <- function(u_bar, At){

  #cv <- seq(0, 1.5, 0.2) #coefficient of variation

  tbl <- data.frame() #final output table

  #sigma[x]

  # for(x in 1:length(sigma)){
  for(x in 1:length(sigsq)){

    site <- 2

    sims <- 100

    amat <- sapply(1:sims, function(y){
      set.seed(10001)
      rand1 <- rnorm(time, mean=0, sd = 1)
      rand2 <- rnorm(time, mean=0, sd = 1)

      z1 <- (rand1-mean(rand1))/sd(rand1) #standard normal
      mean(z1)
      sd(z1)
      z2 <- (rand2-mean(rand2))/sd(rand2)

      X1 <- 1.2  # mean growth rate at optimal site

      Xt1 <- X1+sigma[x]*z1
      mean(Xt1)
      #length(Xt1)
      #sd(Xt1)


      # Xt2 <- (X1-u_bar) + sigma[x]*z2

      Xt2 <- (X1-u_bar) + (sigma[x])*z2

      mean(Xt2)
      tau <- var(Xt2 - Xt1) #tau_bar and tau are same for gaussian
      rho <- u_bar/tau # less than 1/4 for beta <1

      # At <- matrix(c(-1,0.5,0.6,-1), byrow=T, nrow=site, ncol=site) #At is static currently
      #At <- matrix(c(-3,0.5,0.6,-3), byrow=T, nrow=site, ncol=site) #At is static currently
      D.list <- list()

      lambda <- c()
      lambda[1] <- 1
      sum_yt <- c()
      y1 <- c(0.5, 0.5)
      N1 <- c(100, 100)
      sum_yt[1] <- sum(y1)
      a <- c()
      a_n <- c()

      for(j in 1:length(eps))
      {
        yt_mat <- matrix(0, ncol=time, nrow=site)
        N_mat <- matrix(0, ncol=time, nrow=site)

        yt_mat[,1] <- y1


        for(t in 1:(time))
        {
          Bt <- matrix(c(exp(Xt1[t]), 0,0, exp(Xt2[t])), byrow=T, nrow=site)
          # At <- matrix(c(-1, 1, 1, -1), nrow=2)
          I <- matrix(c(1,0,0,1), nrow(At), ncol(At))
          Dt <-  (I+eps[j]*At) %*% Bt
          #print(Dt)
          D.list[[(j-1)*time + t]] <- Dt
        }

      }

      a_mat <- matrix()
      cvar1 <- c()
      cvar2 <- c()
      xxx <- data.frame()

      for(j in 1:length(eps))
      {
        #proportion at both sites
        yt_mat <- matrix(0, ncol=time, nrow=site)
        #    N_mat <- matrix(0, ncol=time, nrow=site)

        yt_mat[,1] <- y1
        N_mat[,1] <- N1

        for(t in 2:time)
        {
          yt_mat[,t] <- (D.list[[((j-1)*time+t-1)]] %*% yt_mat[,t-1])
          sum_yt[t] <- sum(yt_mat[,t])
          yt_mat[,t] <- (yt_mat[,t])/sum_yt[t]


          N_mat[,t] <- (D.list[[((j-1)*time+t-1)]] %*% N_mat[,t-1])

        }

        a[j] <- sum(log(sum_yt))/(time-1)
        cvar1[j] <-cov(yt_mat[1,], Xt1+yt_mat[1,])
        cvar2[j] <- cov(yt_mat[2,], Xt2)
        #                xxx[j,] <- c(a[j], cvar1[j], cvar2[j])
        #length(lambda)
      }
      a #xxx #return a for each of the 100 simulations
    }, simplify = F)

    a_df <- as.data.frame(amat) #create a data drame for each of these a
    junk <- a_df %>% cbind(eps) #bind epsilon row to a_df
    colnames(junk) <- c(1:sims, "epsilon")
    #View(junk)


    a_cal <- junk %>%
      gather(key = sim, value = a, -epsilon) %>%
      group_by(epsilon) %>%
      summarise(a_mean = mean(a), a_check=sum(a)/sims)

    tbl1 <- data.frame(a=round(a_cal$a_mean, 4),
                       eps=a_cal$epsilon,
                       sigma=sigma[x],
                       sigsq=sigsq[x])

    a0 <- tbl1$a[tbl1$eps==0]
    tbl1$a0 <- a0
    tbl1$a_diff <- tbl1$a - tbl1$a0


    #View(a_cal)
    # tbl1 <- data.frame(a=round(a_cal$a_mean, 4),
    #                    eps=a_cal$epsilon, sigma=sigma[x])
    tbl <- rbind(tbl, tbl1)
    #print(c[x])

  }
  return(tbl)
}

# values for difference in expected grwoth rates
ubar <- c(0, 0.15, 0.25, 0.45)
ubar <- c(0.1)
site
time
#plot_list = list()

# eps1 <- seq(0, 0.005, length=10)
# eps1 <- round(eps1, 4)
# eps2 <- c(1e-7, 2e-7, 3e-7, 4e-7, 5e-7,1e-6, 2e-6, 3e-6, 4e-6, 5e-6, 1e-6, 2e-5, 3e-5, 4e-5, 5e-5, 0.0001, 0.0002, 0.0003, 0.0004, 0.0005)
# eps3 <- c(eps1, eps2)
# eps <- eps3[order(eps3)]

# dispersal matrix for two sites
A <- matrix(c(-1, 1, 1, -1), nrow=site)

ubar_df <- data.frame()

for(i in 1:length(ubar))
{
  ubar_df1 <- u_bar_func_n(ubar[i], A)
  ubar_df1$ubar_col <- ubar[i]
  ubar_df <- rbind(ubar_df, ubar_df1)
  print(i)
}

View(ubar_df)
mycol <-c("#00AFBB", "#E69F00", "#D55E00")

# figure 1
#pdf("migplots2.pdf")
plot_list <- list()
for (i in 1:length(ubar)){
  print(i)
  pp <- ggplot(filter(ubar_df, ubar_col==ubar[i]),
               aes(x=eps, y=a, color=sigsq, group=sigsq)) +
    geom_line(linetype="twodash", size=1.1)+
    geom_point(size=2.5, alpha=0.8)+
    labs(title=paste0("\U016B = ", ubar[i]),
         y="stochastic growth rate (a)",
         x=paste0("migration rate", " ", "(m)"))+
    # theme(axis.text.x = element_text(angle = 90))+
    geom_line(linetype="twodash", size=1.1)+
    #    scale_color_gradient(low = "#E7B800", high = "#00AFBB")+
    scale_color_gradientn(colors = mycol, name=expression(sigma^2),
                          #guide = "legend",
                          guide = guide_legend(reverse = TRUE),
                          breaks=c(sigsq),
                          labels=c(sigsq))+
    #labels=c("0.0", "0.2","0.4","0.6","0.8","1.0","1.2","1.4"))+
    # scale_color_manual(name=sigma, values=cc)+
    theme_bw()+
    theme(plot.title = element_text(size=22, face="bold"),
          axis.text.x = element_text(size=16, angle = 90),
          axis.text.y = element_text(size=16),
          axis.title.x = element_text(size=21),
          axis.title.y = element_text(size=21),
          legend.title = element_text(size=24),
          legend.text = element_text(size=16))+
    scale_y_continuous(limits = c(min(ubar_df$a), max(ubar_df$a)))
  plot_list[[i]] <- pp
  print(plot_list[[i]])
}

ggsave("sigsq_mig_twosite_manuscript_slow.png",
       plot = pp,
       # width = width,
       # height = height,
       dpi=400)

# figure 2

plot_list <- list()
for (i in 1:length(ubar)){
  print(i)
  pp <- ggplot(filter(ubar_df, ubar_col==ubar[i]),
               aes(x=log(eps), y=log(a), color=sigsq, group=sigsq)) +
    geom_line(linetype="twodash", size=1.1)+
    geom_point(size=2.5, alpha=0.8)+
    labs(title=paste0("\U016B = ", ubar[i]),
         y="stochastic growth rate (a)",
         x=paste0("migration rate", " ", "(m)"))+
    # theme(axis.text.x = element_text(angle = 90))+
    geom_line(linetype="twodash", size=1.1)+
    #    scale_color_gradient(low = "#E7B800", high = "#00AFBB")+
    scale_color_gradientn(colors = mycol, name=expression(sigma^2),
                          #guide = "legend",
                          guide = guide_legend(reverse = TRUE),
                          breaks=c(sigsq),
                          labels=c(sigsq))+
    #labels=c("0.0", "0.2","0.4","0.6","0.8","1.0","1.2","1.4"))+
    # scale_color_manual(name=sigma, values=cc)+
    theme_bw()+
    theme(plot.title = element_text(size=22, face="bold"),
          axis.text.x = element_text(size=16, angle = 90),
          axis.text.y = element_text(size=16),
          axis.title.x = element_text(size=21),
          axis.title.y = element_text(size=21),
          legend.title = element_text(size=24),
          legend.text = element_text(size=16))+
    xlim(c(-20,-14))
  plot_list[[i]] <- pp
  print(plot_list[[i]])
}


# Figure 10 update!!! kappa plot for multiple sites

sourceCpp('check_kappa_manuscript.cpp')


eps <- c(0,1e-5,5e-5,1e-4,2e-4,4e-4,6e-4,8e-4,1e-3,2e-3,4e-3,6e-3,8e-3,.01)
site <- 4
sigma <- matrix(sqrt(c(0, 0.4,.8,1.2)),nrow = site-1,ncol=4,byrow = TRUE)
# Growth at site 0 is 0 by definition, sigma=0.
time <- 1000
sims <- 1000
kappa <- 2:4



# Function to calculate the stochastic growth rate
# and arrange them in a data frame for subsetting etc later
# Input is At= list of migration matrices
#       and mu = array whose rows correspond to sites, and columns to different choices of mu.
mig_func <- function(At, mu){
  # Check format of migration matrices
  if (!all(unlist(lapply(At,dim))==site)){
    stop("Not all matrices have correct format.")
  }
  if (length(kappa)!=length(At)){
    stop("Wrong number of kappa labels.")
  }
  # Set up identically structured arrays to align parameters and outcomes.
  eps_array <- array(0,c(length(eps),dim(sigma)[2], length(At), dim(mu)[2]))
  sigma_array <- eps_array
  kappa_array <- eps_array
  mu_array <- eps_array
  a_mean <- eps_array #... and to store the final calculations
  a <- numeric(time) # this will be used to store the growth rates for each simulation
  for(i in 1:length(eps)){
    eps_array[i,,,] <- eps[i]
  }
  for(i in 1:dim(sigma)[2]){
    sigma_array[,i,,] <- i
  }
  for(i in 1:length(At)){
    kappa_array[,,i,] <- kappa[i]
  }
  for(i in 1:dim(mu)[2]){
    mu_array[,,,i] <- i # Which mu will be an index pointing to columns of mu
  }
  y0 <- matrix(c(1,rep(0,site-1)),ncol=1)  # Standard initial population vector.
  # Input At as list.
  # Want to use the same Z in every simulation, to make them comparable.
  # mu is a matrix (site-1) x # different choices of mu
  mu <- rbind(0,-mu) # site 0 has relative growth rate 0.
  sigma <- rbind(0,sigma)
  I <- diag(1,nrow=site,ncol=site) #identity matrix

  Z <- array(rnorm((site)*time*sims), c(site,time,sims) ) # These are the growth rates that are going to be scaled appropriately, exponentiated, and placed on the diagonal
  for (i in 1:length(eps)){
    for (which_sigma in 1:dim(sigma)[2]) {
      for (which_mu in 1:dim(mu)[2]){
        Xt <- exp(Z * sigma[,which_sigma]+ mu[,which_mu] ) # Growth rates for given mu and sigma
        for (k in 1:length(At)) {
          migr <- I + eps[i]*At[[k]]
          for(j in 1:sims){
            yt <- lvecprod(c(Xt[,,j]),y0,migr)
            a[j] <- yt/time
          }
          a_mean[i,which_sigma,k,which_mu] <- mean(a)
        }
      }
    }
  }
  data.frame(a=c(a_mean), eps=c(eps_array), sigma = c(sigma_array), kappa = c(kappa_array), mu=c(mu_array) )
}

# ubar_dummy <- c(0.0, 0.1)

#create penalty matrices
#ubar <- c( 0.05, 0.25, 0.50, 0.70)

# New matrices with loops having same rate. matrices corresponding to three kappas

# kappa=2, each site connected to best site
At5 <- matrix(c(-1, 1, 0, 0,
                1.2, -1, .1, .1,
                1, 0, -1, 0,
                1, 0, 0, -1), byrow=T, nrow=site)

# kappa=3, site 1 connected to best site in a cycle of 3
At6 <- matrix(c(-1, 1, 0, 0,
                0, -1, 1, 0,
                1, 0, -1.1, .1,
                1, 0, 0, -1), byrow=T, nrow=site)

# kappa=4, site 1 connected to best site in cycle of 4
At7 <- matrix(c(-1, 1, 0, 0,
                0, -1, 1, 0,
                0, 0, -1, 1,
                1, 0, 0, -1), byrow=T, nrow=site)

# matrices corresponding to three kappas


#At <- list(At2, At3, At4)
# eps1 <- seq(0, 0.005, length=10)
# eps1 <- round(eps1, 4)
# eps2 <- c(1e-10, 1e-7, 2e-7, 3e-7, 4e-7, 5e-7,1e-6, 2e-6, 3e-6, 4e-6, 5e-6, 1e-6, 2e-5, 3e-5, 4e-5, 5e-5, 0.0001, 0.0002, 0.0003, 0.0004, 0.0005)
# eps3 <- c(eps1, eps2)
# eps <- eps3[order(eps3)]
# length(eps)

At <- list(At5,At6,At7)
# At <- list(At2,At3,At4)
mu <- matrix(c(.05,.1,.15,.1,.15,.2,.35,.4,.45),nrow=3)

# kappa <- 2
# mu <- matrix(c(.05,.1,.2),nrow=1)
# At <- list(matrix(c(-1,1,1,-1),2,2))

#### Note that the sites are equidistant here in terms of growth rate difference!
#### can also try out later when they differ not so symmetrically

a_df <- mig_func(At,mu)


# now see what happens for each value of sigma in each case
# loop over sigma and stack plots

rename <- c(`1` = "\U03C3^2 = 0.0",
            `2` = "\U03C3^2 = 0.4",
            `3`=  "\U03C3^2 = 0.8",
            `4` = "\U03C3^2 = 01.2")

my.label_bquote <- function (expr1 = (sigma == .(2)))
{
  quoted1<- substitute(expr1)
  function(variable, value) {
    value <- as.character(value)
    if(variable == 'sigma')
      lapply(value, function(x)
        eval(substitute(bquote(expr1, list(x = x)),list(expr1 = quoted1))))
  }
}

sigmavar <- c(0.0, 0.4, 0.8, 1.2)

a_df6 <- a_df6 %>% mutate(sigmaval = case_when(
  sigma ==1 ~ "0.0",
  sigma ==2 ~ "0.4",
  sigma ==3 ~ "0.8",
  sigma ==4 ~ "1.2",
))

View(a_df6)
mu_low <- ggplot(subset(a_df6,mu==1),
               aes(x=eps, y=a, group=factor(kappa), color=factor(kappa))) +
  geom_line(linetype="twodash", linewidth=0.8)+
  geom_point(size=1.4) +
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  # scale_color_brewer(palette="Dark2")+
  labs(y="stochastic growth rate (a)",
       # x= paste0("migration rate", " ", "(m)"),
       color="kappa",
       title = "low \U016B: sites are similar")+
  facet_wrap(~sigmaval, ncol=7, labeller = label_bquote(sigma^2==.(sigmaval)))+
  theme_bw()+
  theme(plot.title = element_text(size=18, face="bold"),
        axis.text.x = element_text(size=14, angle = 90),
        axis.text.y = element_text(size=14),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_text(size=24),
        legend.text = element_text(size=16))+
  ylim(min(a_df$a), max(a_df$a))
mu_low

mu_high <- ggplot(subset(a_df6,mu==3),
                aes(x=eps, y=a, group=factor(kappa), color=factor(kappa))) +
  geom_line(linetype="twodash", linewidth=0.8)+
  geom_point(size=1.4) +
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  labs(y="stochastic growth rate (a)",
       x= paste0("migration rate", " ", "(m)"),
       color="kappa",
       title = "high \U016B: sites are dissimilar")+
  facet_wrap(~sigmaval, ncol=7, labeller = label_bquote(sigma^2==.(sigmaval)))+
  theme_bw()+
  theme(plot.title = element_text(size=18, face="bold"),
        axis.text.x = element_text(size=14, angle = 90),
        axis.text.y = element_text(size=14),
        axis.title.x = element_text(size=21),
        axis.title.y = element_blank(),
        legend.title = element_text(size=24),
        legend.text = element_text(size=16))+ #,axis.title.y = element_blank())+
   ylim(min(a_df$a), max(a_df$a))
mu_high

kappa_mu_sym_manuscript1 <- ggarrange(mu_low,  mu_high, nrow=2, ncol=1,
                                      common.legend = T,
                                      legend="right")

kappa_mu_sym_manuscript <- annotate_figure(kappa_mu_sym_manuscript1,
                                           left = text_grob("stochastic growth rate (a)", rot = 90, size=20),)

kappa_mu_sym_manuscript
ggsave("kappa_mu_sym_manuscript.png", kappa_mu_sym_manuscript)

eps<- c(1e-10,1e-8,1e-5)

# site <- 2
# kappa <- 2
sigma <- matrix(sqrt(c(0, 0.4,.8,1.2)),nrow = site-1,ncol=4,byrow = TRUE)
# Growth at site 0 is 0 by definition, sigma=0.
time <- 3000
sims <- 3000
eps <- c(1e-10,1e-9,1e-8,1e-4,1e-2,5e-2)
mu <- matrix(c(.05,.1,.2),nrow=1,byrow=TRUE)
At <- list(matrix(c(-1,1,1,-1),2,2))

a_df5 <- mig_func(At,mu)

site <- 4
# # kappa=2, each site connected to best site
# At2 <- matrix(c(-1, 0.4, 0.3, 0.3,
#                 1, -1, 0, 0,
#                 1, 0, -1, 0,
#                 1, 0, 0, -1), byrow=T, nrow=site)
#
# # kappa=3, site 1 connected to best site in a cycle of 3
# At3 <- matrix(c(-1, 1, 0, 0,
#                 0, -1, 0.5, 0.5,
#                 1, 0, -1, 0,
#                 1, 0, 0, -1), byrow=T, nrow=site)
#
# # kappa=4, site 1 connected to best site in cycle of 4
# At4 <- matrix(c(-1, 1, 0, 0,
#                 0, -1, 1, 0,
#                 0, 0, -1, 1,
#                 1, 0, 0, -1), byrow=T, nrow=site)
#
# At <- list(At2,At3,At4)

mu <- matrix(c(.05,.1,.15,.2,.4,rep(.5,10)),nrow=3,ncol=5,byrow=TRUE) # Sites 3 and 4 get big penalty of .5
sigma <- matrix(sqrt(c(0, 0.4,.8,1.2)),nrow = 3,ncol=4,byrow = TRUE)
# Growth at site 0 is 0 by definition, sigma=0.
time <- 2000
sims <- 3000
# eps <- c(1e-10,1e-9,1e-8,1e-4,1e-2,5e-2)
kappa <- c(2,3,4)
a_df6 <- mig_func(At,mu)

# save(a_df5,a_df6,file='mig_sim3')

# Function takes the output of mig_func, and a row number,
# and computes the slope between that row and the next.
#.  It compares that to the upper and lower predictions of the theorem.
slope_comp <- function(adf , rows, kappa){
  adf <- filter(adf, kappa==kappa)
  if (length(rows) == 1){rows <- c(rows,rows+1)}
  s <- log(adf$a[rows[1]]/adf$a[rows[2]])/log(adf$eps[rows[1]]/adf$eps[rows[2]])
  rho <-mu[1,adf$mu[rows[1]]]/sigma[1,adf$sigma[rows[1]]]^2
  upper <- 2*rho*adf$kappa[rows[1]]
  lower <- 2*rho*adf$kappa[rows[1]]/(1+2*rho)
  c(s,upper,lower)
}


s1 <- c()
u1 <- c()
l1 <- c()
x <- c()
for (i in 10:nrow(a_df6)) {
s1[i] <-  slope_comp(a_df, i, 2)[1]
u1[i] <-  slope_comp(a_df, i, 2)[2]
l1[i] <-  slope_comp(a_df, i, 2)[3]
x[i]     <- i
}

slope_dat <- tibble(s1=s1, u1=u1, l1=l1, x=x)

ggplot(slope_dat)+
  geom_point(data=slope_dat, mapping = aes(x=x, y=s1), color="red", size=0.5)+
# geom_line(data=slope_dat, mapping = aes(x=x, y=s1), color="red")
  geom_point(data=slope_dat, mapping = aes(x=x, y=l1), color="darkgreen",  size=0.5)+
geom_point(data=slope_dat, mapping = aes(x=x, y=u1), color="darkblue",  size=0.5)

geom_point(data=slope_dat, mapping = aes(x=x, y=s1))
geom_point(data=slope_dat, mapping = aes(x=x, y=s1))


# facet_wrap(~ubar_col, labeller = as_labeller(ubar_names))


# Figure 11

boxtitle_names <- c(`1`=  "low \U016B low \U03C3^2",
                    `2` = "low \U016B high \U03C3^2",
                    `3` = "high \U016B low \U03C3^2",
                    `4` = "high \U016B high \U03C3^2")
box1<- subset(a_df6, mu==1 & sigmaval==0.4)
box1$group <- "1"
box2<- subset(a_df6, mu==1 & sigmaval==1.2)
box2$group <- "2"
box3<- subset(a_df6, mu==3 & sigmaval==0.4)
box3$group <- "3"
box4<- subset(a_df6, mu==3 & sigmaval==1.2)
box4$group <- "4"

boxdf <- rbind(box1, box2, box3, box4)
View(boxdf)
boxplot_mu_sig_manuscript <- ggplot(boxdf,
                                 aes(x=kappa, y=a, group=kappa)) +
  geom_boxplot(width=0.5,lwd=1, color="#0072B2") +
  geom_jitter(width=0)+
  facet_wrap(~group, labeller = as_labeller(boxtitle_names))+
  labs(x="kappa: length of path",
       y="stochastic growth rate for range of migration rates (m)")+
  # title= "4-site",
  # subtitle= "diff in mean growth rate= (0.01, 0.05, 0.10)", size=0.8)+
  theme_bw()+
  theme(
    plot.title = element_text(size=20, face="bold"),
    plot.subtitle = element_text(size=20),
    axis.text.x = element_text(angle = 90, size=12),
    axis.text.y = element_text(size=12),
    legend.title = element_text( size = 18),
    axis.title.y = element_text(size = 18),
    axis.title.x=element_text(size=18),
    strip.text = element_text(size = 20))+
  scale_x_continuous(breaks=seq(2,4,1))

boxplot_mu_sig_manuscript

ggsave("boxplot_mu_sig_manuscript.png", boxplot_mu_sig_manuscript)

# Figure slope positive negative

sub <- filter(ubar_df, ubar_col==0.1) %>%
  mutate(sign=sign(a_diff)) %>%
  mutate(sign_col = case_when(
    sign == 0 ~ "1",
    sign == 1 ~ "1",
    sign ==-1 ~ "2"
  ))
names(ubar_df)

View(sub)

ss <- ggplot(sub,
             aes(x=factor((sigsq)),
                 y=a_diff, color = (sign_col))) +
  geom_hline(yintercept=0, size=1, color="red", alpha=0.7)+
  geom_point(size=2)+
  geom_line()+
  facet_wrap(~ubar_col)+
  labs(
    y= "Increase in stochastic growth rate: a(m)-a(0)",
    x= expression(sigma^2))+
  scale_color_manual(values = c("#0072B2", "#E69F00"))+    # theme(axis.text.x = element_text(angle = 90))+
  # scale_color_gradient(low = "#E7B800", high = "#00AFBB")+
  # scale_color_gradient(colors = c("#56B4E9"))+
  #                       name="migration rate:\U025B",
  #                       #guide = "legend",
  #                       guide = guide_legend(reverse = TRUE))+
  #labels=c("0.0", "0.2","0.4","0.6","0.8","1.0","1.2","1.4"))+
  # scale_color_manual(name=sigma, values=cc)+
  theme_minimal()+
  theme(axis.text.x = element_text(size=18, angle = 90),
        axis.text.y = element_text(size=18),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        legend.position="none",
        plot.title = element_text(size=18, face="bold"),
        strip.text.x = element_text(size = 20, face = "bold", color="white"),
        strip.text.y = element_text(size = 20, face = "bold" , color="white"))+
  ggtitle("\U016B = 0.1")
# add horizontal line

ss

ggsave("slope_neg_pos_mig_twosite_manuscript.png",
       plot = ss,
       # width = width,
       # height = height,
       dpi=400)

# Figure density


