rm(list = ls())

source('ecoFunctions.R')
library(RColorBrewer)
library(ggExtra)
library(ggtern)
library(ggarchery)

# customization
mut_colors <- setNames(c('#e76d57', '#acd152', '#e7da5c', '#7a9dd2'),
                       c("C59R", "I164L", "N51I", "S108N"))

# load data
file <- '../data/growth_rates_pyr.txt'
#file <- '../data/growth_rates_cyc.txt' # FIXME: uncomment for cycloguanil 
data <- read.table(file)
data[, 2] <- sprintf("%04d", data[, 2])
doses <- c(0, 10^(-2:6))
colnames(data) <- c('genot_id', 'genot_bin', paste('dose_', 0:(length(doses) - 1), sep = ''))

muts <- c('N51I', 'C59R', 'S108N', 'I164L')

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

# filter out highest doses
ge_data$dose <- setNames(c('0', paste('1e', -2:6, sep = '')),
                         as.character(doses))[as.character(ge_data$dose)]
ge_data <- ge_data[!(ge_data$dose %in% c('1e4', '1e5', '1e6')), ]

# plot global epistasis pattern for mutation C59R with no drug
myplot <-
  ggplot(ge_data[ge_data$knock_in == 'C59R' & ge_data$dose == 0, ],
         aes(x = background_f, y = d_f)) +
    geom_abline(slope = 0, intercept = 0, color = '#d1d3d4') +
    geom_point(color = 'black',
               cex = 3) +
    geom_smooth(method = 'lm',
                formula = y~x,
                fullrange = T,
                se = F,
                color = '#1e9fd3') +
    scale_x_continuous(name = expression(paste('Fitness of genetic background, ', italic(f)[italic(B)], sep = '')),
                       breaks = pretty_breaks(n = 3),
                       expand = rep(0.05, 2)) +
    scale_y_continuous(name = expression(paste('Fitness effect\nof focal mutation, ', Delta, italic(f), sep = '')),
                       breaks = pretty_breaks(n = 3),
                       expand = rep(0.05, 2)) +
    theme_bw() +
    theme(aspect.ratio = 0.6,
          panel.grid = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(size = 16, vjust = 0),
          axis.title = element_text(size = 18),
          axis.text = element_text(size = 16),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 14, hjust = 1)) +
    guides(color = guide_colorbar(ticks.colour = NA))
myplot_distrib <- ggMarginal(myplot, type = 'density',
                             xparams = list(bw = 0.075),
                             yparams = list(bw = 0.075))

print(myplot_distrib)
ggsave(myplot_distrib,
       file = '../plots/box1_1.pdf',
       device = 'pdf',
       width = 120,
       height = 100,
       units = 'mm')

# For a given mutation, we want to illustrate:
#   a) how epistatic that mutation is -> sqrt( var(df_i)/var(f_B(i)) ) = gamma_i
#   b) how much of that epistasis is global -> R^2 of linear regression between f_B(i) and df_i
#   c) what is the shape of global epistasis -> slope (b) of linear regression between f_B(i) and df_i
#
# These 3 magnitudes are not independent. For a focal mutation i, we have:
#   R_i^2 = b_i^2 / gamma_i^2
#
# In the case of mut. C59R in the absence of drug:

lmod_i <- lm(formula = d_f ~ background_f,
             data = ge_data[ge_data$knock_in == 'C59R' & ge_data$dose == 0, ])
R_i <- sign(lmod_i$coefficients[2])*sqrt(summary(lmod_i)$r.squared)
b_i <- as.numeric(lmod_i$coefficients[2])
gamma_i <- sqrt( var(ge_data$d_f[ge_data$knock_in == 'C59R' & ge_data$dose == 0]) / var(ge_data$background_f[ge_data$knock_in == 'C59R' & ge_data$dose == 0]) )

x <- 0.5 + log10(gamma_i^2)
y <- 0.5 + log10(R_i^2)
z <- -log10(b_i^2)

# example ternary plot
ggtern(data.frame(x = x, y = y, z = z),
       aes(x, y, z)) +
  geom_point()

# get slopes, gammas, and r squared
slopes <- data.frame(dose = character(0),
                     mut = character(0),
                     slope = numeric(0),
                     pval = numeric(0),
                     R2 = numeric(0),
                     gamma = numeric(0))

for (d in unique(ge_data$dose)) {
  for (m in unique(ge_data$knock_in)) {
    
    fit <- lm(data = ge_data[ge_data$dose == d & ge_data$knock_in == m, ],
              formula = d_f ~ background_f)
    b <- fit$coefficients[2]
    p <- pf(summary(fit)$fstatistic[1],
            summary(fit)$fstatistic[2],
            summary(fit)$fstatistic[3],
            lower.tail=F)
    slopes <- rbind(slopes,
                    data.frame(dose = d,
                               mut = m,
                               slope = b,
                               pval = p,
                               R2 = summary(fit)$r.squared,
                               gamma = sd(ge_data$d_f[ge_data$dose == d & ge_data$knock_in == m]) / sd(ge_data$background_f[ge_data$dose == d & ge_data$knock_in == m])))
    
  }
}

# log-transform coordinates for nernary plot
slopes$x <- 0.5 + log10(slopes$gamma^2)
slopes$y <- 0.5 + log10(slopes$R2)
slopes$z <- -log10(slopes$slope^2)

# add average across all mutations
slopes_avg <- aggregate(formula = . ~ dose,
                        data = slopes[, -c(2)],
                        FUN = mean)
slopes_avg <- cbind(dose = slopes_avg$dose,
                    mut = 'Avg',
                    slopes_avg[, -c(1)])
slopes <- rbind(slopes, slopes_avg)

slopes$dose <- factor(slopes$dose, levels = unique(slopes$dose))
slopes$mut <- factor(slopes$mut, levels = c(names(mut_colors), 'Avg'))

# make segments
plot_this_segments <- data.frame(mut = character(0),
                                 dose = character(0),
                                 x = numeric(0),
                                 xend = numeric(0),
                                 y = numeric(0),
                                 yend = numeric(0),
                                 z = numeric(0),
                                 zend = numeric(0))

for (m in unique(slopes$mut)) {
  
  for (d in 1:(length(unique(slopes$dose)) - 1)) {
    
    plot_this_segments <- rbind(plot_this_segments,
                                data.frame(mut = m,
                                           dose = d,
                                           x = slopes$x[slopes$mut == m & slopes$dose == unique(slopes$dose)[d]],
                                           xend = slopes$x[slopes$mut == m & slopes$dose == unique(slopes$dose)[d + 1]],
                                           y = slopes$y[slopes$mut == m & slopes$dose == unique(slopes$dose)[d]],
                                           yend = slopes$y[slopes$mut == m & slopes$dose == unique(slopes$dose)[d + 1]],
                                           z = slopes$z[slopes$mut == m & slopes$dose == unique(slopes$dose)[d]],
                                           zend = slopes$z[slopes$mut == m & slopes$dose == unique(slopes$dose)[d + 1]]))
    
  }
  
}
plot_this_segments$mut <- factor(plot_this_segments$mut, levels = c(names(mut_colors), 'Avg'))

# ternary plot limits need to be manually adjusted for best visualization
xbreaks <- 0.5 + log10(10^seq(-2, 3, by = 1))
xlabs <- paste('1e', -2:3, sep = '')

ybreaks <- 0.5 + log10(10^seq(-2, 0, by = 1))
ylabs <- paste('1e', -2:0, sep = '')

zbreaks <- -log10(10^seq(-4, 0, by = 1))
zlabs <- paste('1e', -4:0, sep = '')

ggtern(slopes, aes(x, y, z, color = mut, label = dose)) +
  geom_text(size = 2.5) +
  scale_T_continuous(name = expression(italic(R)^2),
                     limits = c(-2.25, 3.50),
                     breaks = ybreaks,
                     labels = ylabs) +
  scale_L_continuous(name = expression(gamma^2),
                     limits = c(-1.75, 4.00),
                     breaks = xbreaks,
                     labels = xlabs) +
  scale_R_continuous(name = expression(italic(b)^2),
                     limits = c(-0.75, 5.00),
                     breaks = zbreaks,
                     labels = zlabs) +
  scale_color_manual(name = 'Mutation',
                     values = c(mut_colors, Avg = 'black')) +
  theme_bw() +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14, hjust = 1))

myplot <-
  ggtern(slopes, aes(x, y, z, color = mut)) +
  scale_T_continuous(name = expression(italic(R)^2),
                     limits = c(-2.25, 3.50),
                     breaks = ybreaks,
                     labels = ylabs) +
  scale_L_continuous(name = expression(gamma^2),
                     limits = c(-1.75, 4.00),
                     breaks = xbreaks,
                     labels = xlabs) +
  scale_R_continuous(name = expression(italic(b)^2),
                     limits = c(-0.75, 5.00),
                     breaks = zbreaks,
                     labels = zlabs) +
  scale_color_manual(name = 'Mutation',
                     values = c(mut_colors, Avg = 'black')) +
  theme_bw() +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14, hjust = 1))

for (d in 1:(length(unique(slopes$dose)) - 1)) myplot <- myplot + geom_segment(data = plot_this_segments,
                                                                               mapping = aes(x = x, xend = xend, y = y, yend = yend, z = z, zend = zend, color = mut),
                                                                               arrow = arrow(angle = 25,
                                                                                             length = unit(0.30, 'cm'),
                                                                                             type = 'closed')) # this gives an error but it seems unavoidable, ggtern and geom_segment do not seem to like each other (the output plot is still correct though)
print(myplot)
ggsave(file = '../plots/box1_tern.pdf',
       device = 'pdf',
       width = 150,
       height = 150,
       units = 'mm')






# plot heatmap
gamma <- seq(0, 1, length.out = 200)
R2 <- seq(0, 1, length.out = 200)
hmap <- do.call(rbind,
               lapply(gamma,
                      FUN = function(g) data.frame(gamma = g,
                                                   R2 = R2,
                                                   b = g*sqrt(R2))))

ggplot(hmap, aes(x = gamma, y = R2, fill = b)) +
  geom_tile() +
  geom_point(data = data.frame(gamma = gamma_i, R2 = R_i^2, b = b_i),
             aes(x = gamma, y = R2),
             cex = 3,
             shape = 1) +
  scale_fill_gradient(name = expression(paste('Slope, ', italic(b[i]))),
                      low = 'white', high = '#1e9fd3',
                      breaks = c(0, 0.5, 1),
                      labels = c('0', '0.5', '1'),
                      limits = c(0, 1)) +
  scale_x_continuous(name = expression(paste('Magnitude of epistatic effects, ', italic(gamma), sep = '')),
                     expand = c(0,0),
                     breaks = seq(0, 1, by = 0.25)) +
  scale_y_continuous(name = expression(paste('Magnitude of\nglobal epistasis, ',~italic(R)^2, sep = '')),
                     expand = c(0,0),
                     breaks = seq(0, 1, by = 0.25)) +
  theme_bw() +
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        strip.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.position = 'top')

ggsave(file = '../plots/box1_hmap.pdf',
       device = 'pdf',
       width = 100,
       height = 100,
       units = 'mm')

# generate pseudo-data to illustrate limit cases

  
