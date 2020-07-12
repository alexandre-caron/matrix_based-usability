
# Configuration -----------------------------------------------------------

# Librairies
library(tidyverse)
library(rstan)
library(bridgesampling)

# Functions
source(file = "scripts/estimates.R", encoding = "UTF8")
model <- stan_model('scripts/draw_mu_s2.stan')


## Parameters
n_sim <- 2000
n_boot <- 1000
m_seq <- c(20,50,100)
n_seq <- c(15,20,30,40,50)
mu_seq <- logit(c(0.1,0.2))
sigma_seq <- c(0.5,1,2)
screen_limit <- 50


# Generating the discovery matrices --------------------------------------- NOT RUN

# fn_sim <- function(m, mu, sigma, n){
#   pl = logitinv(rnorm(n = m, mean = mu, sd = sigma))
#   d = sapply(pl, function(x) rbinom(n = n, size = 1, prob = x))
#   return(d)
# }
# 
# combinations <- expand.grid(n=n_seq, sigma=sigma_seq, mu=mu_seq, m=m_seq)
# 
# set.seed(1984)
# for(i in 1:n_sim){
#   saveRDS(
#     object = apply(combinations, 1, function(x){fn_sim(m=x["m"], mu=x["mu"], sigma=x["sigma"], n=x["n"])}),
#     file = paste0("outputs/discovery_matrices/sim_", formatC(i, width=5, format="d", flag="0"),".rds")
#     )
# }


# Computing estimates ----------------------------------------------------- NOT RUN


# for(i in 1:n_sim){
#   res <- lapply(
#     X = readRDS(paste0("outputs/discovery_matrices/sim_", formatC(i, width=5, format="d", flag="0"),".rds")),
#     FUN = fn_estimates
#   )
#   saveRDS(object = res, file = paste0("outputs/estimates/res_", formatC(i, width=5, format="d", flag="0"),".rds"))
# }


# Merging grid ------------------------------------------------------------ NOT RUN

# res <- as_tibble(do.call(rbind, lapply(
#   X = 1:n_sim,
#   FUN = function(x){cbind(sim = x, do.call(rbind, readRDS(file = paste0("outputs/estimates/res_", formatC(x, width=5, format="d", flag="0"),".rds"))))}
# ))) %>% bind_cols(expand.grid(n=n_seq, sigma=sigma_seq, mu=mu_seq, m=m_seq, sim2=1:n_sim),.) %>% as_tibble() %>%
#   mutate(m_lnbzt = ifelse(m_lnbzt > m + screen_limit, NA, m_lnbzt))
# 
# 
# saveRDS(object = res %>% select(sim,n,m,sigma,mu,m_hat=m_naive,p_hat=p), file = "outputs/estimates/res_naive.rds")
# saveRDS(object = res %>% select(sim,n,m,sigma,mu,m_hat=m_gt,p_hat=p_gt), file = "outputs/estimates/res_gt.rds")
# saveRDS(object = res %>% select(sim,n,m,sigma,mu,m_hat=m_dd,p_hat=p_dd), file = "outputs/estimates/res_dd.rds")
# saveRDS(object = res %>% select(sim,n,m,sigma,mu,m_hat=m_lnbzt,mu_hat=mu_lnbzt,sigma_hat=sigma_lnbzt), file = "outputs/estimates/res_lnbzt.rds")


# Computing confidence intervals ------------------------------------------

# Naive

# res <- readRDS(file = "outputs/estimates/res_naive.rds") %>% 
#   mutate(`2.5%` = as.numeric(NA), `97.5%` = as.numeric(NA)) %>% as.data.frame()
# for(i in 1:nrow(res)){
#   print(i)
#   set.seed(i)
#   res[i, c("2.5%","97.5%")] <- quantile(boot_naive(m=round(res$m_hat[i]), p=res$p_hat[i], n=res$n[i] ,n_boot = n_boot),probs = c(0.025,0.975), na.rm = TRUE)
# }
# saveRDS(object = res, file = "outputs/confidence_interval/res_naive.rds")

# GT

# res <- readRDS(file = "outputs/estimates/res_gt.rds") %>% 
#   mutate(`2.5%` = as.numeric(NA), `97.5%` = as.numeric(NA)) %>% as.data.frame()
# for(i in 1:nrow(res)){
#   print(i)
#   set.seed(i)
#   res[i, c("2.5%","97.5%")] <- quantile(boot_gt(m=round(res$m_hat[i]), p=res$p_hat[i], n=res$n[i] ,n_boot = n_boot),probs = c(0.025,0.975), na.rm = TRUE)
# }
# saveRDS(object = res, file = "outputs/confidence_interval/res_gt.rds")

# DD

# res <- readRDS(file = "outputs/estimates/res_dd.rds") %>% 
#   mutate(`2.5%` = as.numeric(NA), `97.5%` = as.numeric(NA)) %>% as.data.frame()
# for(i in 1:nrow(res)){
#   print(i)
#   set.seed(i)
#   res[i, c("2.5%","97.5%")] <- quantile(boot_dd(m=round(res$m_hat[i]), p=res$p_hat[i], n=res$n[i] ,n_boot = n_boot),probs = c(0.025,0.975), na.rm = TRUE)
# }
# saveRDS(object = res, file = "outputs/confidence_interval/res_dd.rds")

# LNBzt

# res <- readRDS(file = "outputs/estimates/res_lnbzt.rds") %>% 
#   mutate(`2.5%` = as.numeric(NA), `97.5%` = as.numeric(NA)) %>% as.data.frame()
# for(i in 1:nrow(res)){
#   print(i)
#   set.seed(i)
#   res[i, c("2.5%","97.5%")] <- quantile(boot_lnbzt(m=round(res$m_hat[i]), mu=res$mu_hat[i], sigma=res$sigma_hat[i], n=res$n[i], n_boot = n_boot),probs = c(0.025,0.975), na.rm = TRUE)
# }
# saveRDS(object = res, file = "outputs/confidence_interval/res_lnbzt.rds")

# matrix-based 

# fn_estimate_bayes <- function(d, screen_limit=50){
#   n <- nrow(d)
#   ms <- colSums(d)[colSums(d) > 0]
#   j <- sum(ms>0)
#   bayes <- tryCatch(expr = heterogeneous_bayes_ci(j=j, n=n, ms=ms, M=ncol(d)+screen_limit), error = function(e){list(m=as.integer(NA), posterior_m = NA)})
#   return(bayes)
# }
# 
# for(i in 1:n_sim){
#   res <- lapply(
#     X = readRDS(paste0("outputs/discovery_matrices/sim_", formatC(i, width=5, format="d", flag="0"),".rds"))[1:2],
#     FUN = fn_estimate_bayes
#   )
#   saveRDS(object = res, file = paste0("outputs/confidence_interval/res_ci_bayes", formatC(i, width=5, format="d", flag="0"),".rds"))
# }


# Results -----------------------------------------------------------------

rmse <- function(m, m_hat){
  sqrt(mean((m-m_hat)^2,na.rm=T))
}

recoding <- tibble(
  methods_old = c("m_naive", "m_gt", "m_dd", "m_lnbzt", "m_bayes"),
  methods_new = c("naive", "Good Turing", "double deflation", "LNBzt", "matrix-based")
)

# Merging grid 
res <- as_tibble(do.call(rbind, lapply(
  X = 1:n_sim,
  FUN = function(x){cbind(sim = x, do.call(rbind, readRDS(file = paste0("outputs/estimates/res_", formatC(x, width=5, format="d", flag="0"),".rds"))))}
))) %>% bind_cols(expand.grid(n=n_seq, sigma=sigma_seq, mu=mu_seq, m=m_seq, sim2=1:n_sim),.) %>% as_tibble() %>%
  mutate(m_lnbzt = ifelse(m_lnbzt > 500, NA, m_lnbzt)) %>%
  select(n,m,mu,sigma,starts_with("m_")) %>% 
  gather(key="method", value="m_hat", -m, -n, -mu, -sigma) %>% 
  mutate(method = recode(method, !!!setNames(object = recoding$methods_new, nm = recoding$methods_old)))

# Format results 

results <- res %>% 
  group_by(n,m,mu,sigma, method) %>% 
  summarize(
    mean = mean(m_hat, na.rm=T),
    rmse = rmse(m,m_hat),
    is.na = mean(is.na(m_hat)),
    inf_025 = quantile(m_hat, 0.025, na.rm=T),
    sup_975 = quantile(m_hat, 0.975, na.rm=T)) %>% 
  ungroup() 

borne_sup <- 50
palette <- c("#4daf4a","#984ea3","#377eb8","#e41a1c","#ffb300")

# Plot Mean & 95%CI

results %>% 
  mutate(mean = (mean-m)/m*100) %>% 
  mutate(rmse = rmse/m*100) %>% 
  mutate(inf_025 = (inf_025-m)/m*100) %>% 
  mutate(sup_975 = (sup_975-m)/m*100) %>%
  mutate(arrow_quantile = as.numeric(ifelse(sup_975>borne_sup,borne_sup+2,NA))) %>%
  mutate(sup_975 = ifelse(sup_975>borne_sup,borne_sup+2,sup_975)) %>%
  mutate(m = ordered(m,label = c("m = 20","m = 50","m = 100"))) %>%
  mutate(p = paste0("mu = logit(",round(logitinv(mu),2),")\n", "sigma = ",round(sigma,2))) %>%
  mutate(p = ordered(p)) %>%
  ggplot() +
  geom_hline(yintercept = 0, linetype = 2, color = "black", size = 1) +
  geom_pointrange(mapping = aes(x = n, y = mean, ymin = inf_025,ymax = sup_975,group = method, colour = method, shape = method), position = position_dodge(width=6), size = .5) + 
  geom_point(aes(y= arrow_quantile, x = n, colour = method), show.legend = F, size = 2,shape = 24, position = position_dodge(width=6)) +
  geom_line(mapping = aes(x = n, y = mean, group = method, colour = method, shape = method), position = position_dodge(width=6), size = .8) + 
  scale_shape_manual(values = c(5,4,15,16,3)) + scale_color_manual(values = palette) + ylim(c(-borne_sup-5,borne_sup+5)) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  scale_x_continuous(breaks = seq(10,50,10),labels = function(x) paste0("n=", x)) +
  theme(legend.title = element_blank(), axis.text.x  = element_text(angle=90, vjust = .5),legend.position="top", legend.box = "horizontal") +
  ylab("Mean error and 95% fluctuation interval for the prediction of m (%)") + xlab("Sample size") +
  facet_grid(m~p)

# Plot RMSE

results %>% 
  filter(!(method %in% "LNBzt" & mu<(-2) & m<100)) %>% 
  mutate(mean = (mean-m)/m*100) %>% 
  mutate(rmse = rmse/m*100) %>% 
  mutate(inf_025 = (inf_025-m)/m*100) %>% 
  mutate(sup_975 = (sup_975-m)/m*100) %>%
  mutate(arrow_quantile = as.numeric(ifelse(sup_975>borne_sup,borne_sup+2,NA))) %>%
  mutate(sup_975 = ifelse(sup_975>borne_sup,borne_sup+2,sup_975)) %>%
  mutate(m = ordered(m,label = c("m = 20","m = 50","m = 100"))) %>%
  mutate(p = paste0("mu = logit(",round(logitinv(mu),2),")\n", "sigma = ",round(sigma,2))) %>%
  mutate(p = ordered(p)) %>%
  ggplot() +
  geom_hline(yintercept = 0, linetype = 2, color = "black", size = 1) +
  geom_point(mapping = aes(x = n, y = rmse,group = method, colour = method, shape = method), position = position_dodge(width=6)) + 
  geom_line(mapping = aes(x = n, y = rmse, group = method, colour = method, shape = method), position = position_dodge(width=6), size = .8) + 
  scale_shape_manual(values = c(5,4,15,16,3)) + scale_color_manual(values = palette) + ylim(c(-borne_sup-5,borne_sup+5)) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  scale_x_continuous(breaks = seq(10,50,10),labels = function(x) paste0("n=", x)) +
  theme(legend.title = element_blank(), axis.text.x  = element_text(angle=90, vjust = .5),legend.position="top", legend.box = "horizontal") +
  ylab("Root mean square error for the prediction of m (%)") + xlab("Sample size") +
  ylim(c(0,37)) +
  facet_grid(m~p)  
