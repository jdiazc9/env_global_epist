rm(list = ls())

source('ecoFunctions.R')
library(RColorBrewer)
library(ggarchery)
library(ggExtra)
library(tidyr)

# Plasmodium falciparum carrying subsets of 4 mutations: N51I, C59R, S108N, I164L

file <- '../data/growth_rates_pyr.txt'
#file <- '../data/growth_rates_cyc.txt' # FIXME: uncomment for cycloguanil 
data <- read.table(file)
data[, 2] <- sprintf("%04d", data[, 2])
doses <- c(0, 10^(-2:6))
colnames(data) <- c('genot_id', 'genot_bin', paste('dose_', 0:(length(doses) - 1), sep = ''))

muts <- c('N51I', 'C59R', 'S108N', 'I164L') # in the order they appear in the original data frame
mut_colors <- setNames(c('#e76d57', '#acd152', '#e7da5c', '#7a9dd2'), # nice colors to identify mutations
                       c("C59R", "I164L", "N51I", "S108N"))

# nicer labels for drug dose scale
mylabels <- c(0,
              expression(10^-2),
              expression(10^-1),
              expression(10^0),
              expression(10^1),
              expression(10^2),
              expression(10^3))

if (grepl('cyc', file)) { # label for saving plots ('pyr' or 'cyc' for pyrimethamine or cycloguanil respectively)
  saveplot <- 'cyc'
} else if (grepl('pyr', file)) {
  saveplot <- 'pyr'
}

genot_matrix <- do.call(rbind,
                        lapply(data$genot_bin,
                               FUN = function(genot) {
                                 strsplit(genot, split = '')[[1]]
                               }))
colnames(genot_matrix) <- muts
genot_matrix <- as.data.frame(genot_matrix)

# get global epistasis patterns
ge_data <- data.frame(dose = numeric(0),
                      background = character(0),
                      knock_in = character(0),
                      background_f = numeric(0),
                      d_f = numeric(0))

for (d in 1:length(doses)) {
  
  df <- cbind(genot_matrix, fun = data[, 2+d])
  ge_df <- makeGEdata(matrix2string(df))
  
  ge_data <- rbind(ge_data,
                   cbind(dose = doses[d], ge_df))
  
}

ge_data$dose <- setNames(c('0', paste('1e', -2:6, sep = '')),
                         as.character(doses))[as.character(ge_data$dose)]

# filter out highest doses
ge_data <- ge_data[!(ge_data$dose %in% c('1e4', '1e5', '1e6')), ]




### EPISTASIS COORDINATES

# get slopes, R2, and variances ("coordinates" in epistasis map)
epist_coords <- data.frame(dose = character(0),
                           mut = character(0),
                           slope = numeric(0),
                           R2 = numeric(0),
                           var_df = numeric(0),
                           var_fb = numeric(0))

for (d in unique(ge_data$dose)) {
  for (m in unique(ge_data$knock_in)) {
    
    mylm <- lm(data = ge_data[ge_data$dose == d & ge_data$knock_in == m, ],
               formula = d_f ~ background_f)
    b <- mylm$coefficients[2]
    R2 <- summary(mylm)$r.squared
    var_df <- var(ge_data$d_f[ge_data$dose == d & ge_data$knock_in == m])
    var_fb <- var(ge_data$background_f[ge_data$dose == d & ge_data$knock_in == m])
    epist_coords <- rbind(epist_coords,
                          data.frame(dose = d,
                                     mut = m,
                                     slope = b,
                                     R2 = R2,
                                     var_df = var_df,
                                     var_fb = var_fb))
    
  }
}

epist_coords$var_rel <- epist_coords$var_df / epist_coords$var_fb



### GET CONTRIBUTIONS: DF_J, EPS_IJ

getContributions <- function(mut_i, mut_j, save.plots = F) {
  
  # fitness effect
  df <- ge_data[ge_data$knock_in == mut_j, ]
  df <- df[!grepl(mut_i, df$background), ] # keep only backgrounds of i, B(i)
  df_mean <- aggregate(formula = d_f~dose,
                       data = df[, c('dose', 'd_f')],
                       FUN = mean) # means for each dose
  df_wt <- df[df$background == '', ] # fitness effect of mutation in wild-type
  
  ggplot(df, aes(x = dose, y = d_f, group = 1, color = background)) +
    geom_line(data = cbind(df_mean, background = 'mean')) +
    geom_line(data = df_wt) +
    geom_point(alpha = 0.25,
               cex = 2,
               shape = 16) +
    annotate("text", x = 7.25, y = df_wt$d_f[df_wt$dose == '1e3'],
             label = 'On wild-type',
             hjust = 0, color = 'firebrick1') +
    annotate("text", x = 7.25, y = df_mean$d_f[df_mean$dose == '1e3'],
             label = expression(paste('Average on backgrounds of ', mut_i, ',  ',
                                      symbol("\341")*delta*italic(F)[italic(j)]*symbol("\361")[italic(B)(italic(i))])),
             hjust = 0, color = 'black') +
    scale_x_discrete(name = expression(paste('Drug dose (', mu, 'M)', sep = '')),
                     labels = mylabels) +
    scale_y_continuous(name = paste('Fitness effect\nof mut. ', mut_j, sep = ''),
                       breaks = pretty_breaks(n = 2)) +
    scale_color_manual(values = setNames(c('firebrick1', rep('black', 4)),
                                         c(df$background[df$dose == '0'], 'mean'))) +
    theme_bw() +
    theme(aspect.ratio = 0.6,
          panel.grid = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(hjust = 0,
                                    size = 16),
          strip.text.y = element_text(angle = 0),
          axis.title = element_text(size = 18),
          axis.text = element_text(size = 16),
          legend.position = 'none')
  
  if(save.plots) {
    ggsave(file = paste('../plots/dF_', saveplot, '.pdf', sep = ''),
           device = 'pdf',
           width = 200,
           height = 80,
           units = 'mm')
  }
  
  # interaction
  eps <- data.frame(dose = character(0),
                    background = character(0),
                    f_bg = numeric(0),
                    f_bg_i = numeric(0),
                    f_bg_j = numeric(0),
                    f_bg_ij = numeric(0),
                    eps_ij = numeric(0))
  
  bgs <- unique(ge_data$background[!grepl(mut_i, ge_data$background) & !grepl(mut_j, ge_data$background)])
  
  for (d in unique(ge_data$dose)) {
    
    ge_data_i <- ge_data[ge_data$dose == d, ]
    
    f_bg <- unique(ge_data_i[ge_data_i$background %in% bgs, c('background', 'background_f')])
    colnames(f_bg) <- c('background', 'f')
    f_bg <- f_bg[order(f_bg$background), ]
    
    f_bg_i <- ge_data_i[ge_data_i$background %in% bgs & ge_data_i$knock_in == mut_j, ]
    f_bg_i$f <- f_bg_i$background_f + f_bg_i$d_f
    f_bg_i <- f_bg_i[, c('background', 'f')]
    f_bg_i <- f_bg_i[order(f_bg_i$background), ]
    
    f_bg_j <- ge_data_i[ge_data_i$background %in% bgs & ge_data_i$knock_in == mut_i, ]
    f_bg_j$f <- f_bg_j$background_f + f_bg_j$d_f
    f_bg_j <- f_bg_j[, c('background', 'f')]
    f_bg_j <- f_bg_j[order(f_bg_j$background), ]
    
    f_bg_ij <- ge_data_i[(grepl(mut_i, ge_data_i$background) & ge_data_i$knock_in == mut_j) | (grepl(mut_j, ge_data_i$background) & ge_data_i$knock_in == mut_i), ]
    f_bg_ij$background <- sapply(f_bg_ij$background,
                                 FUN = function(x) {
                                   x <- strsplit(x, split = ',')[[1]]
                                   x <- x[!(x %in% c(mut_i, mut_j))]
                                   x <- paste(x, collapse = ',')
                                   return(x)
                                 })
    f_bg_ij$f <- f_bg_ij$background_f + f_bg_ij$d_f
    f_bg_ij <- aggregate(formula = f~background,
                         data = f_bg_ij[, c('background', 'f')],
                         FUN = mean)
    f_bg_ij <- f_bg_ij[order(f_bg_ij$background), ]
    
    eps <- rbind(eps, data.frame(dose = d,
                                 background = f_bg$background,
                                 f_bg = f_bg$f,
                                 f_bg_i = f_bg_i$f,
                                 f_bg_j = f_bg_j$f,
                                 f_bg_ij = f_bg_ij$f,
                                 eps_ij = f_bg_ij$f - f_bg_i$f - f_bg_j$f + f_bg$f))
    
  }
  
  eps_mean <- aggregate(formula = eps_ij~dose,
                        data = eps[, c('dose', 'eps_ij')],
                        FUN = mean)
  eps_wt <- eps[eps$background == '', ]
  
  ggplot(eps, aes(x = dose, y = eps_ij, group = 1, color = background)) +
    geom_line(data = cbind(eps_mean, background = 'mean')) +
    geom_line(data = eps_wt) +
    geom_point(alpha = 0.25,
               cex = 2,
               shape = 16) +
    annotate("text", x = 7.25, y = eps_wt$eps_ij[eps_wt$dose == '1e3'],
             label = 'On wild-type',
             hjust = 0, color = 'firebrick1') +
    annotate("text", x = 7.25, y = eps_mean$eps_ij[eps_mean$dose == '1e3'],
             label = expression(paste('Average across all backgrounds,  ',
                                      symbol("\341")*italic(epsilon)[italic(ij)]*symbol("\361"))),
             hjust = 0, color = 'black') +
    scale_x_discrete(name = expression(paste('Drug dose (', mu, 'M)', sep = '')),
                     labels = mylabels) +
    scale_y_continuous(name = paste('Epistasis between\n', mut_i, ' and ', mut_j, sep = ''),
                       breaks = pretty_breaks(n = 2)) +
    scale_color_manual(values = setNames(c('firebrick1', rep('black', 4)),
                                         c(eps$background[eps$dose == '0'], 'mean'))) +
    theme_bw() +
    theme(aspect.ratio = 0.6,
          panel.grid = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(hjust = 0,
                                    size = 16),
          strip.text.y = element_text(angle = 0),
          axis.title = element_text(size = 18),
          axis.text = element_text(size = 16),
          legend.position = 'none')
  
  if (save.plots) {
    ggsave(file = paste('../plots/eps-ij_', saveplot, '.pdf', sep = ''),
           device = 'pdf',
           width = 200,
           height = 80,
           units = 'mm')
  }
  
  # data frame to return
  return_df <- merge(df[, c('dose', 'background', 'd_f')],
                     eps[, c('dose', 'background', 'eps_ij')],
                     by = c('dose', 'background'),
                     all = T)
  return_df <- data.frame(dose = return_df$dose,
                          focal_mut = mut_i,
                          mut_j = mut_j,
                          background = return_df$background,
                          dF_j = return_df$d_f,
                          eps_ij = return_df$eps_ij)
  
  return(return_df)
  
}

# get all contributions at all doses
contrib <- data.frame(dose = character(0),
                      focal_mut = character(0),
                      mut_j = character(0),
                      background = character(0),
                      dF_j = numeric(0),
                      eps_ij = numeric(0))

for(mut_i in muts) {
  for (mut_j in muts[muts != mut_i]) {
    contrib <- rbind(contrib,
                     getContributions(mut_i, mut_j))
  }
}

# get beta_ij and w_ij from dF_j and eps_ij
sum_dF_j <- aggregate(formula = dF_j ~ focal_mut + dose,
                      data = contrib,
                      FUN = function(x) sum(x^2))

colnames(sum_dF_j)[3] <- 'sum_dF_j_squared'
contrib <- merge(contrib, sum_dF_j, by = c('dose', 'focal_mut'), all = T)

contrib$w_ij <- contrib$dF_j^2 / contrib$sum_dF_j_squared
contrib$beta_ij <- contrib$eps_ij / contrib$dF_j

contrib$w_ij_times_beta_ij <- contrib$eps_ij * contrib$dF_j / contrib$sum_dF_j_squared
contrib$w_ij_times_beta_ij_squared <- contrib$eps_ij^2 / contrib$sum_dF_j_squared

# get (weighted) mean, mean of squares, and CV of beta_ij
contrib_avg <- aggregate(formula = cbind(w_ij_times_beta_ij, w_ij_times_beta_ij_squared) ~ dose + focal_mut,
                         data = contrib,
                         FUN = sum)
contrib_avg$CVw_beta_ij <- sqrt(contrib_avg$w_ij_times_beta_ij_squared - contrib_avg$w_ij_times_beta_ij^2) / contrib_avg$w_ij_times_beta_ij
colnames(contrib_avg)[3:4] <- c('meanw_beta_ij', 'meanw_beta_ij_squared')

# compare with empirical values for slope, var df / var fb, and R^2
epist_coords <- merge(epist_coords, contrib_avg, by.x = c('dose', 'mut'), by.y = c('dose', 'focal_mut'))
epist_coords$dose_rank <- setNames(1:7, unique(epist_coords$dose))[epist_coords$dose]




### FIG 4B-D: EPISTASIS "COORDINATES" VS EFFECTIVE INTERACTIONS FOR ALL MUTATIONS & DOSES

# slope
xy.limits <- range(c(epist_coords$meanw_beta_ij, epist_coords$slope))
ggplot(epist_coords) +
  geom_abline(slope = 1, intercept = 0, color = 'gray') +
  geom_point(aes(x = meanw_beta_ij, y = slope, color = mut, fill = mut, alpha = dose_rank),
             shape = 16,
             cex = 3) +
  geom_point(aes(x = meanw_beta_ij, y = slope, color = mut, fill = mut),
             alpha = 1,
             shape = 1,
             cex = 3) +
  scale_x_continuous(name = expression(symbol('\341')~italic(beta)~symbol('\361')[italic(omega)]),
                     limits = xy.limits,
                     breaks = pretty_breaks(n = 3)) +
  scale_y_continuous(name = 'Empirical slope',
                     limits = xy.limits,
                     breaks = pretty_breaks(n = 3)) +
  scale_color_manual(name = 'Mutation',
                     values = mut_colors) +
  scale_alpha(range = c(0, 1),
              labels = mylabels,
              name = expression(paste('Drug dose (', mu, 'M)', sep = ''))) +
  theme_bw() +
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 16, vjust = 0),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14, hjust = 0),
        panel.border = element_blank(),
        axis.line = element_line(color = 'black'))

# ggsave(file = paste('../plots/slope_vs_beta_', saveplot, '.pdf', sep = ''),
#        device = 'pdf',
#        width = 120,
#        height = 120,
#        units = 'mm')

# relative variance
xy.limits <- range(c(epist_coords$meanw_beta_ij_squared, epist_coords$var_rel))
ggplot(epist_coords) +
  geom_abline(slope = 1, intercept = 0, color = 'gray') +
  geom_point(aes(x = meanw_beta_ij_squared, y = var_rel, color = mut, fill = mut, alpha = dose_rank),
             shape = 16,
             cex = 3) +
  geom_point(aes(x = meanw_beta_ij_squared, y = var_rel, color = mut, fill = mut),
             alpha = 1,
             shape = 1,
             cex = 3) +
  scale_x_log10(name = expression(symbol('\341')~italic(beta)^2~symbol('\361')[italic(omega)]),
                breaks = 10^(-2:2),
                labels = c(expression(10^-2),
                           expression(10^-1),
                           expression(10^0),
                           expression(10^1),
                           expression(10^2)),
                limits = xy.limits) +
  scale_y_log10(name = expression(paste('Empirical var ', Delta, italic(f), ' / var ', italic(f), '(', italic(B), ')',
                                        sep = '')),
                breaks = 10^(-2:2),
                labels = c(expression(10^-2),
                           expression(10^-1),
                           expression(10^0),
                           expression(10^1),
                           expression(10^2)),
                limits = xy.limits) +
  scale_color_manual(name = 'Mutation',
                     values = mut_colors) +
  scale_alpha(range = c(0, 1),
              labels = mylabels,
              name = expression(paste('Drug dose (', mu, 'M)', sep = ''))) +
  theme_bw() +
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 16, vjust = 0),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14, hjust = 0),
        panel.border = element_blank(),
        axis.line = element_line(color = 'black'))

# ggsave(file = paste('../plots/var_vs_beta_', saveplot, '.pdf', sep = ''),
#        device = 'pdf',
#        width = 125,
#        height = 125,
#        units = 'mm')

# R squared
ggplot(epist_coords) +
  geom_function(fun = function(x) 1 / (1 + x^2), color = 'gray') +
  geom_point(aes(x = abs(CVw_beta_ij), y = R2, color = mut, fill = mut, alpha = dose_rank),
             shape = 16,
             cex = 3) +
  geom_point(aes(x = abs(CVw_beta_ij), y = R2, color = mut, fill = mut),
             alpha = 1,
             shape = 1,
             cex = 3) +
  scale_x_log10(name = expression(paste(CV[italic(omega)], ' (', italic(beta), ')',
                                        sep = '')),
                breaks = c(1, 10, 100),
                labels = c(expression(10^0),
                           expression(10^1),
                           expression(10^2))) +
  scale_y_continuous(name = expression(paste('Empirical ', italic(R)^2,
                                        sep = '')),
                breaks = c(0, 0.3, 0.6)) +
  scale_color_manual(name = 'Mutation',
                     values = mut_colors) +
  scale_alpha(range = c(0, 1),
              labels = mylabels,
              name = expression(paste('Drug dose (', mu, 'M)', sep = ''))) +
  theme_bw() +
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 16, vjust = 0),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14, hjust = 0),
        panel.border = element_blank(),
        axis.line = element_line(color = 'black'))

# ggsave(file = paste('../plots/R2_vs_beta_', saveplot, '.pdf', sep = ''),
#        device = 'pdf',
#        width = 125,
#        height = 125,
#        units = 'mm')




