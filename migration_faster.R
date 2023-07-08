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
getwd()
# theme_set(theme_pubr())
sourceCpp('check_kappa.cpp')

# sourceCpp('test.cpp')

sigsq <- seq(0.1, 1.5, length=5)
sigsq <- round(sigsq,1)
sigma <- sqrt(sigsq)
sigsq <- c(0, 0.2, 0.4, 0.6, 0.8, 1.2)
sigsq <- c(0, 0.4, 0.8, 1.2)

# sigsq <- seq(0.1, 0.7, length=4)

sigma <- sqrt(sigsq)

# eps <- c(0,1e-5,5e-5,1e-4,2e-4,4e-4,6e-4,8e-4,1e-3,2e-3,4e-3,6e-3,8e-3,.01)

eps1 <- seq(0, 0.005, length=10)
eps1 <- round(eps1, 4)
eps2 <- c(1e-7, 2e-7, 3e-7, 4e-7, 5e-7,1e-6, 2e-6, 3e-6, 4e-6, 5e-6, 1e-6, 2e-5, 3e-5, 4e-5, 5e-5, 0.0001, 0.0002, 0.0003, 0.0004, 0.0005)
eps3 <- c(eps1, eps2)
eps <- eps3[order(eps3)]
length(eps)

eps <- round(seq(1e-7, 5e-3, length=30),7)

eps
site <- 2
site <- 4
time <- 1000
sims <- 1000
kappa <- c(2)

# Function to calculate the stochastic growth rate 
# and arrange them in a data frame for subsetting etc later
# Input is At= list of migration matrices
#       and mu = array whose rows correspond to sites, and columns to different choices of mu.
mig_func <- function(At, mu){
  # Set up identically structured arrays to align parameters and outcomes.
  eps_array <- array(0, c(length(eps),
                          length(sigsq), 
                          length(At), 
                          dim(mu)[2]))
  sigsq_array <- eps_array
  kappa_array <- eps_array
  mu_array <- eps_array
  a_mean <- eps_array #... and to store the final calculations
  a <- numeric(time) # this will be used to store the growth rates for each simulation
  for(i in 1:length(eps)){
    eps_array[i,,,] <- eps[i]
  }
  for(i in 1:length(sigsq)){
    sigsq_array[,i,,] <- sigsq[i]
  }
  for(i in 1:length(At)){
    kappa_array[,,i,] <- kappa[i] 
  }
  for(i in 1:dim(mu)[2]){
    mu_array[,,,i] <- i # Which mu will be an index pointing to columns of mu
  }
  result <- numeric(sims) # Temporary storage here.
  # Input At as list.
  # Want to use the same X in every simulation, to make them comparable.
  # mu is a matrix (site-1) x # different choices of mu
  
  mu <- rbind(0,-mu) # site 0 has relative growth rate 0.
  
  tbl <- data.frame() #final output table
  
  Bt<- matrix(0,2,2) # Holder for population growth matrix
  I <- matrix(0,2,2)
  diag(I) <- 1
  Z <- array(rnorm(site*time*sims), c(site,time,sims) )
  # loop over sigma
  
  for (i in 1:length(eps)){
    for (s in 1:length(sigsq)) {
      for (k in 1:length(At)) {
        A <- At[[k]]
        for (which_mu in 1:dim(mu)[2]){
          for(j in 1:sims){
            Xt <- c(exp(Z[,,j] * sqrt(sigsq[s])+ mu[,which_mu] )) # Growth rates for given mu and sigma
            y0 <- matrix(1/site, nrow=site, ncol=1)  # proportion at four sites in the beginning
            migr <- diag(1,nrow = site, ncol = site) + eps[i]*A;
            yt <- lvecprod(Xt,y0,migr)
            a[j] <- log(sum(yt))/time
          }
          a_mean[i,s,k,which_mu] <- mean(a)
        }
      }
    }
  }
  data.frame(a=c(a_mean), 
             eps=c(eps_array), 
             sigsq = c(sigsq_array), 
             kappa = c(kappa_array), 
             mu=c(mu_array) )
}

At1 <- -diag(1, nrow=2)
At1[1,2] <- 1
At1[2,1] <- 1
At1
At <- list(At1)
# At <- list(At2,At3,At4)

# mu <- matrix(c(.05,.1,.15,.1,.15,.2,.35,.4,.45),nrow=3)

mu <- matrix(c(.1),nrow=1)

mu1 <- matrix(c(0),nrow=1)
mu2 <- matrix(c(0.15),nrow=1)
mu3 <- matrix(c(0.25),nrow=1)
mu4 <- matrix(c(0.45),nrow=1)

# mu <- matrix(c(0, 0.15, 0.25, 0.45),nrow=4)

#### Note that the sites are equidistant here in terms of growth rate difference!
#### can also try out later when they differ not so symmetrically

# a_df <- mig_func(At,mu)
sigsq
site

a_df1 <- mig_func(At,mu)

a_df_test1 <- mig_func(At,mu1)
a_df_test1$mu <- mu1[1,1]
a_df_test2 <- mig_func(At,mu2)
a_df_test2$mu <- mu2[1,1]
a_df_test3 <- mig_func(At,mu3)
a_df_test3$mu <- mu3[1,1]
a_df_test4 <- mig_func(At,mu4)
a_df_test4$mu <- mu4[1,1]

a_df <- rbind(a_df_test1, 
      a_df_test2,
      a_df_test3,
      a_df_test4)

# now see what happens for each value of sigma in each case
# loop over sigma and stack plots

# mycol <-c("#00AFBB", "#E69F00", "#D55E00")

mycol <-c("#00AFBB",  "#E69F00", "#D55E00", "#0072B2") #new
#pdf("migplots2.pdf")
plot_list <- list()
for (i in 1:dim(mu)[2]){
  xx <- ggplot(a_df1, 
               aes(x=eps, y=a, 
                   color=sigsq, 
                   group=sigsq)) + 
  
    labs(title=paste0("\U016B = ", mu[i]),
         y="stochastic growth rate (a)",
         x=paste0("migration rate", " ", "(m)"))+
    # theme(axis.text.x = element_text(angle = 90))+
    geom_line(linetype="twodash", size=1.4)+
    geom_point(size=2.4, alpha=0.8)+
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
          legend.text = element_text(size=16)
          )
  # scale_y_continuous(limits = c(min(ubar_df$a), max(ubar_df$a)))
  plot_list[[i]] <- xx
  print(plot_list[[i]])
}


ggsave("sigsq_mig_twosite_manuscript.png",
       plot = xx,
       # width = width, 
       # height = height, 
       dpi=400)

ggsave("sigsq_mig_twosite_plotnew1.png",
       plot = xx,
       # width = width, 
       # height = height, 
       dpi=400)

# log-log plot

plot_list <- list()
for (i in 1:dim(mu)[2]){
  ll <- ggplot(a_df1, 
               aes(x=log(eps), 
                   y=log(a), 
                   color=sigsq, 
                   group=sigsq)) + 
    
    labs(title=paste0("\U016B = ", mu[i]),
         y="log(stochastic growth rate): log(a)",
         x=paste0("log(migration rate):", " ", "log(\U025B)"))+
    # theme(axis.text.x = element_text(angle = 90))+
    geom_line(linetype="twodash", size=1.4)+
    geom_point(size=2.4, alpha=0.8)+
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
          legend.text = element_text(size=16)
    )
  # scale_y_continuous(limits = c(min(ubar_df$a), max(ubar_df$a)))
  plot_list[[i]] <- ll
  print(plot_list[[i]])
}

ggsave("log_sigsq_mig_twosite_plotnew1.png",
       plot = ll,
       # width = width, 
       # height = height, 
       dpi=400)
# multiple mu plot

mu_names <- c(`0` = paste0("\U016B = ", "0.00"),
                  `0.15` = paste0("\U016B = ", "0.15"),
                  `0.25` = paste0("\U016B = ", "0.25"),
                  `0.45` = paste0("\U016B = ", "0.45")
  )

pp <- ggplot(a_df, 
             aes(x=eps, y=a, 
                 color=sigsq, 
                 group=sigsq)) + 
  labs(y="stochastic growth rate (a)",
       x=paste0("migration rate", " ", "(\U025B)"))+
  # theme(axis.text.x = element_text(angle = 90))+
  geom_line(linetype="twodash", size=1)+
  geom_point(size=1.8)+
  facet_wrap(~mu, labeller = as_labeller(mu_names))+
  scale_color_gradientn(colors = mycol, name=expression(sigma^2), 
                        #guide = "legend", 
                        guide = guide_legend(reverse = TRUE),
                        breaks=c(sigsq), 
                        labels=c(sigsq))+
  theme_bw()+
  ylim(c(min(a_df$a)-0.01, max(a_df$a)+0.02))+
  theme(plot.title = element_text(size=22, face="bold"),
        axis.text.x = element_text(size=16, angle = 90),
        axis.text.y = element_text(size=16),
        axis.title.x = element_text(size=21),
        axis.title.y = element_text(size=21),
        legend.title = element_text(size=24),
        legend.text = element_text(size=16),
        # The new stuff
        strip.text = element_text(size = 17, face="bold"))


ggsave("ubar_sigsq_mig_twosite_plotnew1.png",
       plot = pp,
       # width = width, 
       # height = height, 
       dpi=400)

setwd("/Users/harmanjaggi/Documents/Research/dispersal/migration code")
# log plot

xx <- ggplot(a_df_test, 
             aes(x=log(eps), y=log(a), 
                 color=sigsq, 
                 group=sigsq)) + 
  
  labs(title=paste0("\U016B = ", mu[i]),
       y="stochastic growth rate (a)",
       x=paste0("migration rate", " ", "(\U025B)"))+
  # theme(axis.text.x = element_text(angle = 90))+
  geom_line(linetype="twodash", size=1.4)+
  geom_point(size=2.4, alpha=0.8)+
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
        legend.text = element_text(size=16)
  )

xx
# different mubar

mu0 <- ggplot(a_df_test, 
             aes(x=eps, y=a, 
                 color=sigsq, 
                 group=sigsq)) + 
  
  labs(title=paste0("\U016B = ", mu[1,1]),
       y="stochastic growth rate (a)",
       x=paste0("migration rate", " ", "(\U025B)"))+
  # theme(axis.text.x = element_text(angle = 90))+
  geom_line(linetype="twodash", size=1.4)+
  geom_point(size=2.4, alpha=0.8)+
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
        legend.text = element_text(size=16)
  )


mu0

mu1 <- ggplot(a_df_test, 
              aes(x=eps, y=a, 
                  color=sigsq, 
                  group=sigsq)) + 
  
  labs(title=paste0("\U016B = ", mu[1,1]),
       y="stochastic growth rate (a)",
       x=paste0("migration rate", " ", "(\U025B)"))+
  # theme(axis.text.x = element_text(angle = 90))+
  geom_line(linetype="twodash", size=1.4)+
  geom_point(size=2.4, alpha=0.8)+
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
        legend.text = element_text(size=16)
  )


mu1

log_test <- filter(a_df_test, sigsq==c(0.3,0.8,1.0))
plot_list <- list()
for (i in 1:dim(mu)[2]){
  tt <- ggplot(log_test),
               aes(x=log(eps), y=log(a), 
                   color=sigsq, 
                   group=sigsq)) + 
    labs(title=paste0("\U016B = ", mu[i]),
         y="stochastic growth rate (a)",
         x=paste0("migration rate", " ", "(\U025B)"))+
    # theme(axis.text.x = element_text(angle = 90))+
    geom_line(linetype="twodash", size=1.4)+
    geom_point(size=2.4, alpha=0.8)+
    # #    scale_color_gradient(low = "#E7B800", high = "#00AFBB")+
    # scale_color_gradientn(colors = mycol, name=expression(sigma^2), 
    #                       #guide = "legend", 
    #                       guide = guide_legend(reverse = TRUE))+
    #                       breaks=c(sigsq),
    #                       labels=c(sigsq))+
    #labels=c("0.0", "0.2","0.4","0.6","0.8","1.0","1.2","1.4"))+
    # scale_color_manual(name=sigma, values=cc)+
    theme_bw()+
    theme(plot.title = element_text(size=22, face="bold"),
          axis.text.x = element_text(size=16, angle = 90),
          axis.text.y = element_text(size=16),
          axis.title.x = element_text(size=21),
          axis.title.y = element_text(size=21),
          legend.title = element_text(size=24),
          legend.text = element_text(size=16)
    )
  # scale_y_continuous(limits = c(min(ubar_df$a), max(ubar_df$a)))
  plot_list[[i]] <- tt
  print(plot_list[[i]])
}

(0.0001)^(0.1/0.4)
(0.0001)^(0.1/0.6)
(0.0001)^(0.1/0.8)
(0.0001)^(0.1)
(0.0001)^(0.1/1.2)
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


mu_low<-ggplot(subset(a_df,mu==1), 
               aes(x=eps, y=a, group=factor(kappa), color=factor(kappa))) +
  geom_line(linetype="twodash", linewidth=0.8)+
  geom_point(size=1.4) +  
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  # scale_color_brewer(palette="Dark2")+
  labs(y="stochastic growth rate (a)",
       x= paste0("migration rate", " ", "(\U025B)"),
       color="kappa",
       title = "\U016B low: sites are similar")+
  facet_wrap(~sigma, ncol=7)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90))+#,legend.position = "none")+
  ylim(min(a_df$a), max(a_df$a))+
  xlim(c(0, 0.003))+
  theme(plot.title = element_text(size=20, face="bold"),
        plot.subtitle = element_text(size=20),
        axis.text.x = element_text(angle = 90, size=12),
        axis.text.y = element_text(size=12),
        legend.title = element_text( size = 18),
        axis.title.y = element_blank(),
        axis.title.x=element_blank())
mu_low

mu_high<-ggplot(subset(a_df,mu==3), 
                aes(x=eps, y=a, group=factor(kappa), color=factor(kappa))) +
  geom_line(linetype="twodash", linewidth=0.8)+
  geom_point(size=1.4) +
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  labs(y="stochastic growth rate (a)",
       x= paste0("migration rate", " ", "(\U025B)"),
       color="kappa",
       title = "\U016B high: sites are dissimilar")+
  facet_wrap(~sigma, ncol=7)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90))+ #,axis.title.y = element_blank())+
  xlim(c(0, 0.003))+
  ylim(min(a_df$a), max(a_df$a))+
  theme(plot.title = element_text(size=20, face="bold"),
        legend.title = element_text( size = 19, face="bold"),
        axis.text.x = element_text(angle = 90, size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size=17))
mu_high

kappa_mu_plot_draft1 <- ggarrange(mu_low,  mu_high, nrow=2, ncol=1,
                                  common.legend = T,
                                  legend="right")
kappa_mu_plot_draft <- annotate_figure(kappa_mu_plot_draft1,
                                       left = text_grob("stochastic growth rate (a)", rot = 90, size=16),
                                       
)
kappa_mu_plot_draft
ggsave("kappa_mu_draft_plot.png", kappa_mu_sym)



box1<- subset(a_df, mu==1 & sigsq==0.4)
box1$group <- "1" 
box2<- subset(a_df, mu==1 & sigsq==1.2)
box2$group <- "2" 
box3<- subset(a_df, mu==3 & sigsq==0.4)
box3$group <- "3" 
box4<- subset(a_df, mu==3 & sigsq==1.2)
box4$group <- "4" 
sigsq
boxdf <- rbind(box1, box2, box3, box4)

# library(ggtext)
# "low \U016B low \U03C3^2"
boxtitle_names <- c(`1`=  "low \U016B low \U03C3^2",
                    `2` = "low \U016B high \U03C3^2",
                    `3` = "high \U016B low \U03C3^2",
                    `4` = "high \U016B high \U03C3^2")

# facet_wrap(~ubar_col, labeller = as_labeller(ubar_names))


boxplot_mu_sig_plotnew <- ggplot(boxdf, 
                                 aes(x=kappa, y=a, group=kappa)) +
  geom_boxplot(width=0.5,lwd=1, color="#0072B2") +
  geom_jitter(width=0)+
  facet_wrap(~group, labeller = as_labeller(boxtitle_names))+
  labs(x="kappa: length of path", 
       y="stochastic growth rate for range of epsilon")+
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


View(a_df)


# updated David code June 26th 

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

# setwd("/Users/harmanjaggi/Documents/Research/dispersal/migration code")

# theme_set(theme_pubr())
sourceCpp('check_kappa.cpp')

eps <- c(0,1e-5,5e-5,1e-4,2e-4,4e-4,6e-4,8e-4,1e-3,2e-3,4e-3,6e-3,8e-3,.01)
site <- 2
sigma <- matrix(sqrt(c(0, 0.4,.8,1.2)),nrow = site-1,ncol=4,byrow = TRUE)
# Growth at site 0 is 0 by definition, sigma=0.
time <- 1000
sims <- 1000
kappa <- 2:4

# Function to calculate the stochastic growth rate 
# and arrange them in a data frame for subsetting etc later
# Input is At= list of migration matrices
# and mu = array whose rows correspond to sites, 
# and columns to different choices of mu.
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
At <- list(At5,At6,At7)
At <- list(At2,At3,At4)
mu <- matrix(c(.05,.1,.15,.1,.15,.2,.35,.4,.45),nrow=3)

kappa <- 2
mu <- matrix(c(.05,.1,.2),nrow=1)
At <- list(matrix(c(-1,1,1,-1),2,2))

#### Note that the sites are equidistant here in terms of growth rate difference!
#### can also try out later when they differ not so symmetrically

#a_df <- mig_func(At,mu)


# now see what happens for each value of sigma in each case
# loop over sigma and stack plots

mu_low<-ggplot(subset(a_df,mu==1), 
               aes(x=eps, y=a, group=factor(kappa), color=factor(kappa))) +
  geom_line(linetype="twodash", linewidth=0.8)+
  geom_point(size=1.4) +  
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  # scale_color_brewer(palette="Dark2")+
  labs(y="stochastic growth rate (a)",
       x= paste0("migration rate", " ", "(\U025B)"),
       color="kappa",
       title = "low: sites are similar")+
  facet_wrap(~sigma, ncol=7)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90))+#,legend.position = "none")+
  ylim(min(a_df$a), max(a_df$a))
mu_low

mu_high<-ggplot(subset(a_df,mu==3), 
                aes(x=eps, y=a, group=factor(kappa), color=factor(kappa))) +
  geom_line(linetype="twodash", linewidth=0.8)+
  geom_point(size=1.4) +
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  labs(y="stochastic growth rate (a)",
       x= paste0("migration rate", " ", "(\U025B)"),
       color="kappa",
       title = "high: sites are dissimilar")+
  facet_wrap(~sigma, ncol=7)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90))+ #,axis.title.y = element_blank())+
  ylim(min(a_df$a), max(a_df$a))

#kappa_mu_sym <- ggarrange(mu_low,  mu_high, nrow=2, ncol=1)  
#kappa_mu_sym

#ggsave("kappa_mu_sym5.png", kappa_mu_sym)

eps<- c(1e-10,1e-8,1e-5)

site <- 2
kappa <- 2
sigma <- matrix(sqrt(c(0, 0.4,.8,1.2)),nrow = site-1,ncol=4,byrow = TRUE)
# Growth at site 0 is 0 by definition, sigma=0.
time <- 3000
sims <- 3000
eps <- c(1e-10,1e-9,1e-8,1e-4,1e-2,5e-2)
mu <- matrix(c(.05,.1,.2),nrow=1,byrow=TRUE)
At <- list(matrix(c(-1,1,1,-1),2,2))

a_df5 <- mig_func(At,mu)

site <- 4
# kappa=2, each site connected to best site
At2 <- matrix(c(-1, 0.4, 0.3, 0.3, 
                1, -1, 0, 0, 
                1, 0, -1, 0,
                1, 0, 0, -1), byrow=T, nrow=site)

# kappa=3, site 1 connected to best site in a cycle of 3
At3 <- matrix(c(-1, 1, 0, 0, 
                0, -1, 0.5, 0.5, 
                1, 0, -1, 0,
                1, 0, 0, -1), byrow=T, nrow=site)

# kappa=4, site 1 connected to best site in cycle of 4
At4 <- matrix(c(-1, 1, 0, 0, 
                0, -1, 1, 0, 
                0, 0, -1, 1,
                1, 0, 0, -1), byrow=T, nrow=site)

At <- list(At2,At3,At4)

mu <- matrix(c(.05,.1,.15,.2,.4,rep(.5,10)),nrow=3,ncol=5,byrow=TRUE) # Sites 3 and 4 get big penalty of .5
sigma <- matrix(sqrt(c(0, 0.4,.8,1.2)),nrow = 3,ncol=4,byrow = TRUE)
# Growth at site 0 is 0 by definition, sigma=0.
time <- 2000
sims <- 3000
eps <- c(1e-10,1e-9,1e-8,1e-4,1e-2,5e-2)
kappa <- c(2,3,4)
a_df6 <- mig_func(At,mu)

save(a_df5,a_df6,file='mig_sim3')

# Function takes the output of mig_func, and a row number, and computes the slope between that row and the next.
#.  It compares that to the upper and lower predictions of the theorem.
slope_comp <- function(adf , rows){
  if (length(rows) == 1){rows <- c(rows,rows+1)}
  s <- log(adf$a[rows[1]]/adf$a[rows[2]])/log(adf$eps[rows[1]]/adf$eps[rows[2]])
  rho <-mu[1,adf$mu[rows[1]]]/sigma[1,adf$sigma[rows[1]]]^2
  upper <- 2*rho*adf$kappa[rows[1]]
  lower <- 2*rho*adf$kappa[rows[1]]/(1+2*rho)
  c(s,upper,lower)
}


#### David github R history

}
for(i in 1:length(At)){
  kappa_array[,,i,] <- kappa[i]
}
for(i in 1:dim(mu)[2]){
  mu_array[,,,i] <- i # Which mu will be an index pointing to columns of mu
}
result <- numeric(sims) # Temporary storage here.
# Input At as list.
# Want to use the same X in every simulation, to make them comparable.
# mu is a matrix (site-1) x # different choices of mu
mu <- rbind(0,-mu) # site 0 has relative growth rate 0.
tbl <- data.frame() #final output table
Bt<- matrix(0,4,4) # Holder for population growth matrix
I <- matrix(0,4,4)
diag(I) <- 1
Z <- array(rnorm(site*time*sims), c(site,time,sims) )
# loop over sigma
for (i in 1:length(eps)){
  for (s in 1:length(sigma)) {
    for (k in 1:length(At)) {
      A <- At[[k]]
      for (which_mu in 1:dim(mu)[2]){
        for(j in 1:sims){
          Xt <- c(exp(Z[,,j] * sigma[s]+ mu[,which_mu] )) # Growth rates for given mu and sigma
          y0 <- matrix(1/site,nrow=site,ncol=1)  # proportion at four sites in the beginning
          migr <- diag(1,nrow = site,ncol = site)+ eps[i]*A;
          yt <- lvecprod(Xt,y0,migr)
          a[j] <- log(sum(yt))/time
        }
        a_mean[i,s,k,which_mu] <- mean(a)
      }
    }
  }
}
data.frame(a=c(a_mean), eps=c(eps_array), sigma = c(sigma_array), kappa = c(A_array), mu=c(mu_array) )
}
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
# kappa=2, each site connected to best site
At2 <- matrix(c(-1, 0.4, 0.3, 0.3,
                1, -1, 0, 0,
                1, 0, -1, 0,
                1, 0, 0, -1), byrow=T, nrow=site)
# kappa=3, site 1 connected to best site in a cycle of 3
At3 <- matrix(c(-1, 1, 0, 0,
                0, -1, 0.5, 0.5,
                1, 0, -1, 0,
                1, 0, 0, -1), byrow=T, nrow=site)
# kappa=4, site 1 connected to best site in cycle of 4
At4 <- matrix(c(-1, 1, 0, 0,
                0, -1, 1, 0,
                0, 0, -1, 1,
                1, 0, 0, -1), byrow=T, nrow=site)
#At <- list(At2, At3, At4)
At <- list(At5,At6,At7)
Ad <- matrix(c(-1, 0, 0, -1), nrow=2)
mu <- matrix(c(.05,.1,.15,.1,.15,.2,.35,.4,.45),nrow=3)
a_df <- mig_func(At,mu)
# Function to calculate the stochastic growth rate
# and arrange them in a data frame for subsetting etc later
# Input is At= list of migration matrices
#       and mu = array whose rows correspond to sites, and columns to different choices of mu.
mig_func <- function(At, mu){
  # Set up identically structured arrays to align parameters and outcomes.
  eps_array <- array(0,c(length(eps),length(sigma), length(At), dim(mu)[2]))
  sigma_array <- eps_array
  kappa_array <- eps_array
  mu_array <- eps_array
  a_mean <- eps_array #... and to store the final calculations
  a <- numeric(time) # this will be used to store the growth rates for each simulation
  for(i in 1:length(eps)){
    eps_array[i,,,] <- eps[i]
  }
  for(i in 1:length(sigma)){
    sigma_array[,i,,] <- sigma[i]
  }
  for(i in 1:length(At)){
    kappa_array[,,i,] <- kappa[i]
  }
  for(i in 1:dim(mu)[2]){
    mu_array[,,,i] <- i # Which mu will be an index pointing to columns of mu
  }
  result <- numeric(sims) # Temporary storage here.
  # Input At as list.
  # Want to use the same X in every simulation, to make them comparable.
  # mu is a matrix (site-1) x # different choices of mu
  mu <- rbind(0,-mu) # site 0 has relative growth rate 0.
  tbl <- data.frame() #final output table
  Bt<- matrix(0,4,4) # Holder for population growth matrix
  I <- matrix(0,4,4)
  diag(I) <- 1
  Z <- array(rnorm(site*time*sims), c(site,time,sims) )
  # loop over sigma
  for (i in 1:length(eps)){
    for (s in 1:length(sigma)) {
      for (k in 1:length(At)) {
        A <- At[[k]]
        for (which_mu in 1:dim(mu)[2]){
          for(j in 1:sims){
            Xt <- c(exp(Z[,,j] * sigma[s]+ mu[,which_mu] )) # Growth rates for given mu and sigma
            y0 <- matrix(1/site,nrow=site,ncol=1)  # proportion at four sites in the beginning
            migr <- diag(1,nrow = site,ncol = site)+ eps[i]*A;
            yt <- lvecprod(Xt,y0,migr)
            a[j] <- log(sum(yt))/time
          }
          a_mean[i,s,k,which_mu] <- mean(a)
        }
      }
    }
  }
  data.frame(a=c(a_mean), eps=c(eps_array), sigma = c(sigma_array), kappa = c(kappa_array), mu=c(mu_array) )
}
a_df <- mig_func(At,mu)
dim(a_df)

mu_low<-ggplot(subset(a_df,mu==1),
               aes(x=eps, y=a, group=factor(kappa), color=factor(kappa))) +
  geom_line(linetype="twodash", linewidth=0.8)+
  geom_point(size=1.4) +
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  # scale_color_brewer(palette="Dark2")+
  labs(y="stochastic growth rate (a)",
       x= paste0("migration rate", " ", "(\U025B)"),
       color="kappa",
       title = "low: sites are similar")+
  facet_wrap(~sigma, ncol=4)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "none")+
  ylim(min(a_df$a), max(a_df$a))
mu_low
time <- 1000
system.time(a_df <- mig_func(At,mu))
mu_high<-ggplot(subset(a_df,mu==3),
                aes(x=eps, y=a, group=factor(kappa), color=factor(kappa))) +
  geom_line(linetype="twodash", size=0.8)+
  geom_point(size=1.4) +
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  labs(y="stochastic growth rate (a)",
       x= paste0("migration rate", " ", "(\U025B)"),
       color="kappa",
       title = "high: sites are dissimilar")+
  facet_wrap(~sigma, ncol=4)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90),
        axis.title.y = element_blank())+
  ylim(min(a_df$a), max(a_df$a))
mu_high<-ggplot(subset(a_df,mu==3),
                aes(x=eps, y=a, group=factor(kappa), color=factor(kappa))) +
  geom_line(linetype="twodash", linewidth=0.8)+
  geom_point(size=1.4) +
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  labs(y="stochastic growth rate (a)",
       x= paste0("migration rate", " ", "(\U025B)"),
       color="kappa",
       title = "high: sites are dissimilar")+
  facet_wrap(~sigma, ncol=4)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90),
        axis.title.y = element_blank())+
  ylim(min(a_df$a), max(a_df$a))
kappa_mu_sym <- ggarrange(mu_low,  mu_high, nrow=1, ncol=2)
ggsave("kappa_mu_sym2.png", kappa_mu_sym)
a_df <- mig_func(At,mu)
mu_low<-ggplot(subset(a_df,mu==1),
               aes(x=eps, y=a, group=factor(kappa), color=factor(kappa))) +
  geom_line(linetype="twodash", linewidth=0.8)+
  geom_point(size=1.4) +
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  # scale_color_brewer(palette="Dark2")+
  labs(y="stochastic growth rate (a)",
       x= paste0("migration rate", " ", "(\U025B)"),
       color="kappa",
       title = "low: sites are similar")+
  facet_wrap(~sigma, ncol=4)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "none")+
  ylim(min(a_df$a), max(a_df$a))
mu_low
mu_high<-ggplot(subset(a_df,mu==3),
                aes(x=eps, y=a, group=factor(kappa), color=factor(kappa))) +
  geom_line(linetype="twodash", linewidth=0.8)+
  geom_point(size=1.4) +
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  labs(y="stochastic growth rate (a)",
       x= paste0("migration rate", " ", "(\U025B)"),
       color="kappa",
       title = "high: sites are dissimilar")+
  facet_wrap(~sigma, ncol=4)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90),
        axis.title.y = element_blank())+
  ylim(min(a_df$a), max(a_df$a))
kappa_mu_sym <- ggarrange(mu_low,  mu_high, nrow=1, ncol=2)
ggsave("kappa_mu_sym2.png", kappa_mu_sym)
mu_low<-ggplot(subset(a_df,mu==1),
               aes(x=eps, y=a, group=factor(kappa), color=factor(kappa))) +
  geom_line(linetype="twodash", linewidth=0.8)+
  geom_point(size=1.4) +
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  # scale_color_brewer(palette="Dark2")+
  labs(y="stochastic growth rate (a)",
       x= paste0("migration rate", " ", "(\U025B)"),
       color="kappa",
       title = "low: sites are similar")+
  facet_wrap(~sigma, ncol=4)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "none")+
  ylim(min(a_df$a), max(a_df$a))
mu_high<-ggplot(subset(a_df,mu==3),
                aes(x=eps, y=a, group=factor(kappa), color=factor(kappa))) +
  geom_line(linetype="twodash", linewidth=0.8)+
  geom_point(size=1.4) +
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  labs(y="stochastic growth rate (a)",
       x= paste0("migration rate", " ", "(\U025B)"),
       color="kappa",
       title = "high: sites are dissimilar")+
  facet_wrap(~sigma, ncol=4)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90),
        axis.title.y = element_blank())+
  ylim(min(a_df$a), max(a_df$a))
kappa_mu_sym <- ggarrange(mu_low,  mu_high, nrow=1, ncol=2)
ggsave("kappa_mu_sym2.png", kappa_mu_sym)
eps <- c(0,1e-5,5e-5,1e-4,2e-4,4e-4,6e-4,8e-4,1e-3)
a_df <- mig_func(At,mu)
mu_low<-ggplot(subset(a_df,mu==1),
               aes(x=eps, y=a, group=factor(kappa), color=factor(kappa))) +
  geom_line(linetype="twodash", linewidth=0.8)+
  geom_point(size=1.4) +
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  # scale_color_brewer(palette="Dark2")+
  labs(y="stochastic growth rate (a)",
       x= paste0("migration rate", " ", "(\U025B)"),
       color="kappa",
       title = "low: sites are similar")+
  facet_wrap(~sigma, ncol=4)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "none")+
  ylim(min(a_df$a), max(a_df$a))
mu_low
mu_high<-ggplot(subset(a_df,mu==3),
                aes(x=eps, y=a, group=factor(kappa), color=factor(kappa))) +
  geom_line(linetype="twodash", linewidth=0.8)+
  geom_point(size=1.4) +
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  labs(y="stochastic growth rate (a)",
       x= paste0("migration rate", " ", "(\U025B)"),
       color="kappa",
       title = "high: sites are dissimilar")+
  facet_wrap(~sigma, ncol=4)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90),
        axis.title.y = element_blank())+
  ylim(min(a_df$a), max(a_df$a))
kappa_mu_sym <- ggarrange(mu_low,  mu_high, nrow=1, ncol=2)
ggsave("kappa_mu_sym2.png", kappa_mu_sym)
mu_low<-ggplot(subset(a_df,mu==1),
               aes(x=eps, y=a, group=factor(kappa), color=factor(kappa))) +
  geom_line(linetype="twodash", linewidth=0.8)+
  geom_point(size=1.4) +
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  # scale_color_brewer(palette="Dark2")+
  labs(y="stochastic growth rate (a)",
       x= paste0("migration rate", " ", "(\U025B)"),
       color="kappa",
       title = "low: sites are similar")+
  facet_wrap(~sigma, ncol=4)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "none")+
  ylim(min(a_df$a), max(a_df$a))
mu_high<-ggplot(subset(a_df,mu==3),
                aes(x=eps, y=a, group=factor(kappa), color=factor(kappa))) +
  geom_line(linetype="twodash", linewidth=0.8)+
  geom_point(size=1.4) +
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  labs(y="stochastic growth rate (a)",
       x= paste0("migration rate", " ", "(\U025B)"),
       color="kappa",
       title = "high: sites are dissimilar")+
  facet_wrap(~sigma, ncol=4)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90),
        axis.title.y = element_blank())+
  ylim(min(a_df$a), max(a_df$a))
kappa_mu_sym <- ggarrange(mu_low,  mu_high, nrow=1, ncol=2)
ggsave("kappa_mu_sym2.png", kappa_mu_sym)
kappa_mu_sym <- ggarrange(mu_low,  mu_high, nrow=1, ncol=2)
ggsave("kappa_mu_sym2.png", kappa_mu_sym)
At <- list(At2,At3,At4)
a_df <- mig_func(At,mu)
mu_low<-ggplot(subset(a_df,mu==1),
               aes(x=eps, y=a, group=factor(kappa), color=factor(kappa))) +
  geom_line(linetype="twodash", linewidth=0.8)+
  geom_point(size=1.4) +
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  # scale_color_brewer(palette="Dark2")+
  labs(y="stochastic growth rate (a)",
       x= paste0("migration rate", " ", "(\U025B)"),
       color="kappa",
       title = "low: sites are similar")+
  facet_wrap(~sigma, ncol=4)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "none")+
  ylim(min(a_df$a), max(a_df$a))
mu_high<-ggplot(subset(a_df,mu==3),
                aes(x=eps, y=a, group=factor(kappa), color=factor(kappa))) +
  geom_line(linetype="twodash", linewidth=0.8)+
  geom_point(size=1.4) +
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  labs(y="stochastic growth rate (a)",
       x= paste0("migration rate", " ", "(\U025B)"),
       color="kappa",
       title = "high: sites are dissimilar")+
  facet_wrap(~sigma, ncol=4)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90),
        axis.title.y = element_blank())+
  ylim(min(a_df$a), max(a_df$a))
kappa_mu_sym <- ggarrange(mu_low,  mu_high, nrow=1, ncol=2)
ggsave("kappa_mu_sym2.png", kappa_mu_sym)
eps <- c(0,1e-5,5e-5,1e-4,2e-4,4e-4,6e-4,8e-4,1e-3,1e-2)
a_df <- mig_func(At,mu)
mu_low<-ggplot(subset(a_df,mu==1),
               aes(x=eps, y=a, group=factor(kappa), color=factor(kappa))) +
  geom_line(linetype="twodash", linewidth=0.8)+
  geom_point(size=1.4) +
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  # scale_color_brewer(palette="Dark2")+
  labs(y="stochastic growth rate (a)",
       x= paste0("migration rate", " ", "(\U025B)"),
       color="kappa",
       title = "low: sites are similar")+
  facet_wrap(~sigma, ncol=4)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "none")+
  ylim(min(a_df$a), max(a_df$a))
mu_high<-ggplot(subset(a_df,mu==3),
                aes(x=eps, y=a, group=factor(kappa), color=factor(kappa))) +
  geom_line(linetype="twodash", linewidth=0.8)+
  geom_point(size=1.4) +
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  labs(y="stochastic growth rate (a)",
       x= paste0("migration rate", " ", "(\U025B)"),
       color="kappa",
       title = "high: sites are dissimilar")+
  facet_wrap(~sigma, ncol=4)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90),
        axis.title.y = element_blank())+
  ylim(min(a_df$a), max(a_df$a))
kappa_mu_sym <- ggarrange(mu_low,  mu_high, nrow=1, ncol=2)
ggsave("kappa_mu_sym3.png", kappa_mu_sym)
time=100
a_df <- mig_func(At,mu)
mu_low<-ggplot(subset(a_df,mu==1),
               aes(x=eps, y=a, group=factor(kappa), color=factor(kappa))) +
  geom_line(linetype="twodash", linewidth=0.8)+
  geom_point(size=1.4) +
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  # scale_color_brewer(palette="Dark2")+
  labs(y="stochastic growth rate (a)",
       x= paste0("migration rate", " ", "(\U025B)"),
       color="kappa",
       title = "low: sites are similar")+
  facet_wrap(~sigma, ncol=4)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "none")+
  ylim(min(a_df$a), max(a_df$a))
mu_high<-ggplot(subset(a_df,mu==3),
                aes(x=eps, y=a, group=factor(kappa), color=factor(kappa))) +
  geom_line(linetype="twodash", linewidth=0.8)+
  geom_point(size=1.4) +
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  labs(y="stochastic growth rate (a)",
       x= paste0("migration rate", " ", "(\U025B)"),
       color="kappa",
       title = "high: sites are dissimilar")+
  facet_wrap(~sigma, ncol=4)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90),
        axis.title.y = element_blank())+
  ylim(min(a_df$a), max(a_df$a))
kappa_mu_sym <- ggarrange(mu_low,  mu_high, nrow=1, ncol=2)
ggsave("kappa_mu_sym4.png", kappa_mu_sym)
At[[1]]
eps <- c(0,1e-5,5e-5,1e-4,2e-4,4e-4,6e-4,8e-4,1e-3,2e-3,4e-3,6e-3,8e-3,.01)
sigma <- seq(0, 1.4, 0.2)
sigma <- seq(0, 0.4,.6,.8,1,1.2,1.6)
sigma <- c(0, 0.4,.6,.8,1,1.2,1.6)
a_df <- mig_func(At,mu)
time
time <- 1000
a_df <- mig_func(At,mu)
mu_low<-ggplot(subset(a_df,mu==1),
               aes(x=eps, y=a, group=factor(kappa), color=factor(kappa))) +
  geom_line(linetype="twodash", linewidth=0.8)+
  geom_point(size=1.4) +
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  # scale_color_brewer(palette="Dark2")+
  labs(y="stochastic growth rate (a)",
       x= paste0("migration rate", " ", "(\U025B)"),
       color="kappa",
       title = "low: sites are similar")+
  facet_wrap(~sigma, ncol=7)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "none")+
  ylim(min(a_df$a), max(a_df$a))
mu_low
mu_high<-ggplot(subset(a_df,mu==3),
                aes(x=eps, y=a, group=factor(kappa), color=factor(kappa))) +
  geom_line(linetype="twodash", linewidth=0.8)+
  geom_point(size=1.4) +
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  labs(y="stochastic growth rate (a)",
       x= paste0("migration rate", " ", "(\U025B)"),
       color="kappa",
       title = "high: sites are dissimilar")+
  facet_wrap(~sigma, ncol=7)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90),
        axis.title.y = element_blank())+
  ylim(min(a_df$a), max(a_df$a))
mu_high
sigma <- c(0, 0.4,.8,1.2)
eps <- c(0,1e-5,5e-5,1e-4,2e-4,4e-4,6e-4,8e-4,1e-3,2e-3,4e-3,6e-3,8e-3,.01)
site <- 4
time <- 1000
At <- list(At2,At3,At4)
a_df <- mig_func(At,mu)
mu_low<-ggplot(subset(a_df,mu==1),
               aes(x=eps, y=a, group=factor(kappa), color=factor(kappa))) +
  geom_line(linetype="twodash", linewidth=0.8)+
  geom_point(size=1.4) +
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  # scale_color_brewer(palette="Dark2")+
  labs(y="stochastic growth rate (a)",
       x= paste0("migration rate", " ", "(\U025B)"),
       color="kappa",
       title = "low: sites are similar")+
  facet_wrap(~sigma, ncol=7)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "none")+
  ylim(min(a_df$a), max(a_df$a))
mu_high<-ggplot(subset(a_df,mu==3),
                aes(x=eps, y=a, group=factor(kappa), color=factor(kappa))) +
  geom_line(linetype="twodash", linewidth=0.8)+
  geom_point(size=1.4) +
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  labs(y="stochastic growth rate (a)",
       x= paste0("migration rate", " ", "(\U025B)"),
       color="kappa",
       title = "high: sites are dissimilar")+
  facet_wrap(~sigma, ncol=7)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90),
        axis.title.y = element_blank())+
  ylim(min(a_df$a), max(a_df$a))
kappa_mu_sym <- ggarrange(mu_low,  mu_high, nrow=2, ncol=1)
ggsave("kappa_mu_sym5.png", kappa_mu_sym)
mu_low<-ggplot(subset(a_df,mu==1),
               aes(x=eps, y=a, group=factor(kappa), color=factor(kappa))) +
  geom_line(linetype="twodash", linewidth=0.8)+
  geom_point(size=1.4) +
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  # scale_color_brewer(palette="Dark2")+
  labs(y="stochastic growth rate (a)",
       x= paste0("migration rate", " ", "(\U025B)"),
       color="kappa",
       title = "low: sites are similar")+
  facet_wrap(~sigma, ncol=7)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90))+#,legend.position = "none")+
  ylim(min(a_df$a), max(a_df$a))
mu_high<-ggplot(subset(a_df,mu==3),
                aes(x=eps, y=a, group=factor(kappa), color=factor(kappa))) +
  geom_line(linetype="twodash", linewidth=0.8)+
  geom_point(size=1.4) +
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  labs(y="stochastic growth rate (a)",
       x= paste0("migration rate", " ", "(\U025B)"),
       color="kappa",
       title = "high: sites are dissimilar")+
  facet_wrap(~sigma, ncol=7)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90))+ #,axis.title.y = element_blank())+
  ylim(min(a_df$a), max(a_df$a))
kappa_mu_sym <- ggarrange(mu_low,  mu_high, nrow=2, ncol=1)
ggsave("kappa_mu_sym5.png", kappa_mu_sym)

