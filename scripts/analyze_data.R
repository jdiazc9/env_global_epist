rm(list = ls())

source('ecoFunctions.R')
library(RColorBrewer)
library(ggarchery)
library(ggExtra)
library(tidyr)
library(lmtest)

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


tst <- ge_data[ge_data$dose == '1e1' & ge_data$knock_in == 'I164L', ]
tst$background <- sapply(tst$background,
                         FUN = function(x) {
                           out <- strsplit(x, split = ',')[[1]]
                           out <- substr(out, 1, 1)
                           out <- paste(out, collapse = '')
                           return(out)
                         })
tst$background[1] <- 'Wild-type'
tst$background <- factor(tst$background,
                         levels = tst$background)

ggplot(tst, aes(x = background, y = d_f, fill = background)) +
  geom_bar(stat = 'identity',
           width = 0.85) +
  scale_fill_manual(values = c('red', rep('gray', 6))) +
  scale_x_discrete(name = 'Genetic background') +
  scale_y_continuous(name = expression(Delta*italic(F))) +
  geom_hline(yintercept = 0) +
  theme_bw() +
  theme(aspect.ratio = 0.5,
        panel.grid = element_blank(),
        legend.position = 'none',
        axis.title = element_text(size = 18),
        axis.text.x = element_text(size = 16,
                                   angle = 45,
                                   hjust = 1),
        axis.text.y = element_text(size = 16),
        plot.title = element_text(size = 18)) +
  ggtitle('Focal mut.: I164L, dose: 10 uM')

ggplot(tst, aes(x = background_f, y = d_f)) +
  geom_hline(yintercept = 0,
             color = 'gray') +
  geom_smooth(method = 'lm',
              formula = y ~ x,
              se = F,
              color = 'firebrick1') +
  geom_point(cex = 3) +
  scale_x_continuous(name = expression(paste('Genetic background fitness, ', italic(F)[B], sep = ''))) +
  scale_y_continuous(name = expression(Delta*italic(F))) +
  theme_bw() +
  theme(aspect.ratio = 0.6,
        panel.grid = element_blank(),
        legend.position = 'none',
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        plot.title = element_text(size = 18)) +
  ggtitle('Focal mut.: I164L, dose: 10 uM')

ggsave(file = paste('../plots/for_Alvaro_2.pdf', sep = ''),
       device = 'pdf',
       width = 120,
       height = 100,
       units = 'mm')






### FIG 1C: FEE OF FOCAL MUTATION (C59R) WITH NO DRUG

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
    scale_x_continuous(name = expression(paste('Fitness of genetic background, ', italic(f), '(', italic(B), ')', sep = '')),
                       breaks = pretty_breaks(n = 3),
                       expand = rep(0.05, 2)) +
    scale_y_continuous(name = expression(paste('Fitness effect\nof mut. C59R, ', Delta, italic(f), sep = '')),
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
       file = paste('../plots/FEEs_C59R_dose0_', saveplot, '.pdf', sep = ''),
       device = 'pdf',
       width = 100,
       height = 100,
       units = 'mm')

# print R squared
mylm <- lm(formula = d_f~background_f,
           data = ge_data[ge_data$knock_in == 'C59R' & ge_data$dose == 0, ])
print(summary(mylm)$r.squared)



### FIG 1D: EPISTASIS "MAP"

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

# plot (only for dose = 0 in this panel)
ggplot(epist_coords[epist_coords$dose == 0, ],
       aes(x = var_rel, y = R2, color = mut, label = mut)) +
  geom_point(cex = 3) +
  geom_text(hjust = -0.5) +
  geom_blank(data = epist_coords) +
  scale_color_manual(values = mut_colors) +
  scale_x_log10(name = expression(paste('var ', Delta, italic(f), ' / ', 'var ', italic(f), '(', italic(B), ')', sep = '')),
                breaks = 10^c(-2, 0, 2),
                labels = c(expression(10^-2),
                           expression(10^0),
                           expression(10^2))) +
  scale_y_continuous(name = expression(italic(R)^2),
                     breaks = c(0, 0.3, 0.6, 0.9)) +
  theme_bw() +
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 16, vjust = 0),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14, hjust = 1),
        legend.position = 'none')

ggsave(file = paste('../plots/epistMap_dose0_', saveplot, '.pdf', sep = ''),
       device = 'pdf',
       width = 75,
       height = 75,
       units = 'mm')




### INSETS IF FIG 1D: ILLUSTRATION OF LIMITS IN EPISTASIS MAP

# auxiliary function: sample two variables x and y with a given correlation
# from: https://stats.stackexchange.com/questions/15011/generate-a-random-variable-with-a-defined-correlation-to-an-existing-variables
complement <- function(y, rho, x) {
  if (missing(x)) x <- rnorm(length(y)) # Optional: supply a default if `x` is not given
  y.perp <- residuals(lm(x ~ y))
  rho * sd(y.perp) * y + y.perp * sd(y) * sqrt(1 - rho^2)
}

set.seed(1)
df <- data.frame(case = character(0),
                 fb = numeric(0),
                 df = numeric(0))

# top-left corner: low epistasis but fitness effect highly correlated with background fitness
fb <- rnorm(8, mean = 1, sd = 0.2)
df <- rbind(df,
            data.frame(case = 'A',
                       fb = fb,
                       df = 0.8 + 0.4*fb + rnorm(8, 0, 0.001)))

# top-right corner: large epistasis & fitness effect largely correlates with background fitness
# negative slope
fb <- rnorm(8, mean = 1, sd = 0.2)
df <- rbind(df,
            data.frame(case = 'B.1',
                       fb = fb,
                       df = 3.2 - 2*fb + rnorm(8, 0, 0.001)))

# positive slope
fb <- rnorm(8, mean = 1, sd = 0.2)
df <- rbind(df,
            data.frame(case = 'B.2',
                       fb = fb,
                       df = -0.5 + 2*fb + rnorm(8, 0, 0.001)))

# bottom-left corner: (quasi-)additive case; low epistasis, same fitness effect across all backgrounds
fb <- rnorm(8, mean = 1, sd = 0.2)
df <- rbind(df,
            data.frame(case = 'C',
                       fb = fb,
                       df = 1 + rnorm(8, 0, 0.001)))

# bottom-right corner: idiosyncratic case; epistasis is high but no correlation between fb and df
fb <- rnorm(8, mean = 1, sd = 0.2)
df <- rbind(df,
            data.frame(case = 'D',
                       fb = fb,
                       df = 1 + 4*complement(fb, 0)))

ggplot(df, aes(x = fb, y = df)) +
  geom_smooth(method = 'lm',
              color = 'gray',
              formula = y~x,
              se = F,
              fullrange = T) +
  geom_point(cex = 3) +
  facet_wrap(~case,
             nrow = 1) +
  scale_x_continuous(name = expression(paste(italic(f), '(', italic(B), ')', sep = ''))) +
  scale_y_continuous(name = expression(paste(Delta, italic(f), sep = ''))) +
  theme_bw() +
  theme(aspect.ratio = 0.6,
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 16, vjust = 0),
        axis.title = element_text(size = 18),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14, hjust = 1),
        legend.position = 'none',
        panel.border = element_blank()) +
  annotate("segment", x=-Inf, xend=-Inf, y=Inf, yend=-Inf) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)

ggsave(file = paste('../plots/epistMap_insets.pdf', sep = ''),
       device = 'pdf',
       width = 150,
       height = 75,
       units = 'mm')



### FIG 2B: HOW DOES EPISTASIS CHANGE AS DRUG DOSE INCREASES FOR MUT. C59R?

# nicer labels for drug dose scale
mylabels <- c(0,
              expression(10^-2),
              expression(10^-1),
              expression(10^0),
              expression(10^1),
              expression(10^2),
              expression(10^3))

ge_data$dose <- factor(ge_data$dose, levels = c('0', paste('1e', -2:3, sep = '')))
ge_data <- merge(ge_data, epist_coords[, c('dose', 'mut', 'slope')], by.x = c('dose', 'knock_in'), by.y = c('dose', 'mut'), all = T)

ggplot(ge_data[ge_data$knock_in == 'C59R', ],
       aes(x = background_f, y = d_f, color = slope)) +
  geom_abline(slope = 0, intercept = 0, color = '#d1d3d4') +
  geom_point(color = 'black') +
  geom_smooth(method = 'lm',
              formula = y~x,
              fullrange = T,
              se = F) +
  facet_wrap(~dose, nrow = 1, scales = 'free_x',
             labeller = label_bquote(.(setNames(mylabels[1:7],
                                                levels(ge_data$dose))[dose]))) +
  scale_x_continuous(name = expression(paste('Fitness of genetic background, ', italic(f), '(', italic(B), ')', sep = '')),
                     breaks = pretty_breaks(n = 2),
                     expand = rep(0.1, 2)) +
  scale_y_continuous(name = expression(paste('Fitness effect\nof mut. C59R, ', Delta, italic(f), sep = '')),
                     breaks = pretty_breaks(n = 3),
                     expand = rep(0.1, 2)) +
  scale_color_gradient2(low = '#1e9fd3',
                        high = '#ef3a37',
                        mid = 'white',
                        midpoint = 0,
                        breaks = c(-0.4, 0, 0.4)) +
  theme_bw() +
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 16, vjust = 0),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14, hjust = 1)) +
  guides(color = guide_colorbar(ticks.colour = NA))

ggsave(file = paste('../plots/FEEs_C59R_', saveplot, '.pdf', sep = ''),
       device = 'pdf',
       width = 250,
       height = 150,
       units = 'mm')



### FIG 2C: EPISTASIS MAP VS DOSE CHANGE

# segments
map_segments <- data.frame(dose = numeric(0),
                           mut = character(0),
                           xstart = numeric(0),
                           xend = numeric(0),
                           ystart = numeric(0),
                           yend = numeric(0))
for (d in 1:(length(unique(epist_coords$dose)) - 1)) {
  for (m in unique(epist_coords$mut)) {
    
    map_segments <- rbind(map_segments,
                          data.frame(dose = unique(epist_coords$dose)[d],
                                     mut = m,
                                     xstart = epist_coords$var_rel[epist_coords$dose == unique(epist_coords$dose)[d] & epist_coords$mut == m],
                                     xend = epist_coords$var_rel[epist_coords$dose == unique(epist_coords$dose)[d + 1] & epist_coords$mut == m],
                                     ystart = epist_coords$R2[epist_coords$dose == unique(epist_coords$dose)[d] & epist_coords$mut == m],
                                     yend = epist_coords$R2[epist_coords$dose == unique(epist_coords$dose)[d + 1] & epist_coords$mut == m]))
    
  }
}

# plot map
ggplot(map_segments, # this ggplot call spits a warning but the plot seems correct
       aes(x = xstart, y = ystart, xend = xend, yend = yend, color = mut)) +
  geom_arrowsegment(arrow_positions = 1,
                    arrows = arrow(angle = 25,
                                   type = 'open',
                                   length = unit(0.25, "cm")),
                    linewidth = 0.5) +
  geom_point(data = epist_coords[epist_coords$dose == '0', ],
             aes(x = var_rel, y = R2, xend = 0, yend = 0, color = mut),
             cex = 3) +
  scale_color_manual(values = mut_colors) +
  scale_x_log10(name = expression(paste('var ', Delta, italic(f), ' / ', 'var ', italic(f), '(', italic(B), ')', sep = '')),
                breaks = 10^c(-2, 0, 2),
                labels = c(expression(10^-2),
                           expression(10^0),
                           expression(10^2))) +
  scale_y_continuous(name = expression(italic(R)^2),
                     breaks = c(0, 0.3, 0.6)) +
  theme_bw() +
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 16, vjust = 0),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14, hjust = 1),
        legend.position = 'none')

ggsave(file = paste('../plots/epistMap_alldoses_', saveplot, '.pdf', sep = ''),
       device = 'pdf',
       width = 75,
       height = 75,
       units = 'mm')



### FIG 2D-F: VAR_DF, R2 AND SLOPE VS. DOSE

epist_coords$dose_rank <- setNames(1:7, unique(epist_coords$dose))[epist_coords$dose]

ggplot(epist_coords,
       aes(x = dose_rank, y = var_rel, color = mut)) +
  geom_line() +
  scale_color_manual(values = mut_colors) +
  scale_x_continuous(name = expression(paste('Drug dose (', mu, 'M)', sep = '')),
                     breaks = c(1, 3, 5, 7),
                     labels = mylabels[c(1, 3, 5, 7)]) +
  scale_y_log10(name = expression(paste('var ', Delta, italic(f), ' / ', 'var ', italic(f), '(', italic(B), ')', sep = '')),
                breaks = 10^c(-2, 0, 2),
                labels = c(expression(10^-2),
                           expression(10^0),
                           expression(10^2))) +
  theme_bw() +
  theme(aspect.ratio = 0.6,
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 16, vjust = 0),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14, hjust = 1),
        legend.position = 'none',
        panel.border = element_blank(),
        axis.line = element_line(color = 'black'))

ggsave(file = paste('../plots/var_vs_doses_', saveplot, '.pdf', sep = ''),
       device = 'pdf',
       width = 75,
       height = 75,
       units = 'mm')

ggplot(epist_coords,
       aes(x = dose_rank, y = R2, color = mut)) +
  geom_line() +
  scale_color_manual(values = mut_colors) +
  scale_x_continuous(name = expression(paste('Drug dose (', mu, 'M)', sep = '')),
                     breaks = c(1, 3, 5, 7),
                     labels = mylabels[c(1, 3, 5, 7)]) +
  scale_y_continuous(name = expression(italic(R)^2),
                     breaks = c(0, 0.3, 0.6)) +
  theme_bw() +
  theme(aspect.ratio = 0.6,
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 16, vjust = 0),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14, hjust = 1),
        legend.position = 'none',
        panel.border = element_blank(),
        axis.line = element_line(color = 'black'))

ggsave(file = paste('../plots/R2_vs_doses_', saveplot, '.pdf', sep = ''),
       device = 'pdf',
       width = 72,
       height = 72,
       units = 'mm')

ggplot(epist_coords,
       aes(x = dose_rank, y = slope, color = mut)) +
  geom_line() +
  scale_color_manual(values = mut_colors) +
  scale_x_continuous(name = expression(paste('Drug dose (', mu, 'M)', sep = '')),
                     breaks = c(1, 3, 5, 7),
                     labels = mylabels[c(1, 3, 5, 7)]) +
  scale_y_continuous(name = 'Slope',
                     breaks = pretty_breaks(n = 3)) +
  theme_bw() +
  theme(aspect.ratio = 0.6,
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 16, vjust = 0),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14, hjust = 1),
        legend.position = 'none',
        panel.border = element_blank(),
        axis.line = element_line(color = 'black'))

ggsave(file = paste('../plots/slope_vs_doses_', saveplot, '.pdf', sep = ''),
       device = 'pdf',
       width = 71,
       height = 71,
       units = 'mm')
  
  



### FIGURE FOR BOX 1: ILLUSTRATE MEANING OF EPSILON_IJ AND DELTA F_J

fitness_colors <- c('white', '#e9672d')

df <- data[, c('genot_bin', 'dose_0')]
colnames(df) <- c('genot_bin', 'fitness')

# reorder mutations for cleaner visualization
df$genot <- sapply(df$genot,
                  FUN = function(G) {
                    
                    G <- strsplit(G, split = '')[[1]]
                    which_muts <- which(as.numeric(G) == 1)
                    G <- paste(sort(substr(muts, 1, 1)[which_muts]), collapse = '')
                    
                    return(G)
                    
                  })
df$genot_bin <- sapply(df$genot_bin,
                       FUN = function(G) {
                         
                         G <- strsplit(G, split = '')[[1]]
                         G <- paste(G[c(2, 4, 1, 3)], collapse = '')
                         
                         return(G)
                         
                       })
df <- rbind(df, data.frame(genot_bin = '0101',
                           fitness = NA,
                           genot = 'IS')) # no data for this allele

df <- df[order(match(df$genot_bin, c('0000',
                                     '0001', '0100', '0010', '1000',  
                                     '0011', '0101', '1001', '0110', '1010', '1100',  
                                     '0111', '1011', '1101', '1110',
                                     '1111'))), ]

# edges of fitness graph
edges <- do.call(rbind,
                 lapply(df$genot[1:(nrow(df) - 1)],
                        FUN = function(Gstart) {
                          
                          Gend <- df$genot[nchar(df$genot) == (1 + nchar(Gstart))]
                          if (Gstart != '') {
                            is_descendent <- sapply(Gend, FUN = function(G) all(strsplit(Gstart, split = '')[[1]] %in% strsplit(G, split = '')[[1]]))
                            Gend <- Gend[is_descendent]
                          }
                          
                          return(data.frame(Gstart = Gstart,
                                            Gend = Gend))
                          
                        }))

# positioning of nodes/edges
df$ypos <- 4 - nchar(df$genot)
df$xpos <- c(0,
             0:3 - 1.5,
             0:5 - 2.5,
             0:3 - 1.5,
             0)

df$genot[df$genot_bin == '0000'] <- 'w.t.'
edges$Gstart[edges$Gstart == ''] <- 'w.t.'

edges$xstart <- setNames(df$xpos, df$genot)[edges$Gstart]
edges$xend <- setNames(df$xpos, df$genot)[edges$Gend]
edges$ystart <- setNames(df$ypos, df$genot)[edges$Gstart]
edges$yend <- setNames(df$ypos, df$genot)[edges$Gend]

#plot
ggplot() +
  geom_segment(data = edges,
               aes(x = xstart, y = ystart, xend = xend, yend = yend),
               color = 'gray') +
  # geom_arrowsegment(data = edges,
  #                   aes(x = xstart, y = ystart, xend = xend, yend = yend),
  #                   arrows = arrow(angle = 25,
  #                                  type = 'open',
  #                                  length = unit(0.25, "cm")),
  #                   arrow_positions = 0.8) +
  geom_point(data = df,
             aes(x = xpos, y = ypos, fill = fitness),
             color = 'black',
             shape = 21,
             cex = 15) +
  geom_text(data = df,
            aes(x = xpos, y = ypos, label = genot),
            color = 'black') +
  scale_fill_gradient(low = fitness_colors[1], high = fitness_colors[2],
                      breaks = c(1, 1.2, 1.4)) +
  theme_void() +
  theme(aspect.ratio = 0.75,
        panel.grid = element_blank(),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14, hjust = 1)) +
  guides(fill = guide_colorbar(ticks.colour = NA))

ggsave(file = paste('../plots/fitnessGraph_dose0_', saveplot, '.pdf', sep = ''),
       device = 'pdf',
       width = 100,
       height = 70,
       units = 'mm')

# 3D plot
df <- df[df$genot %in% c('N', 'IN', 'CN', 'CIN'), c('genot', 'fitness')]
df <- rbind(df, data.frame(genot = 'CIN_additive',
                           fitness = df$fitness[1] + (df$fitness[2] - df$fitness[1]) + (df$fitness[3] - df$fitness[1])))

df$mut_C <- as.numeric(grepl('C', df$genot))
df$mut_I <- as.numeric(grepl('I', df$genot))

col.pal <- colorRampPalette(c(fitness_colors[1], fitness_colors[2]))
mycolors <- col.pal(100)

persp(c(0, 1), c(0, 1), matrix(c(df$fitness[1:3], NA), nrow = 2, byrow = F),
      phi = 10,
      theta = 55,
      xlab = '',
      ylab = '',
      zlab = '',
      zlim = c(min(df$fitness), max(df$fitness)),
      col = fitness_colors[2]
)
par(new=TRUE)
persp(c(0, 1), c(0, 1), matrix(c(NA, df$fitness[2:4]), nrow = 2, byrow = F),
      phi = 10,
      theta = 55,
      xlab = '',
      ylab = '',
      zlab = '',
      zlim = c(min(df$fitness), max(df$fitness)),
      col = fitness_colors[2]
)
par(new=TRUE)
persp(c(0, 1), c(0, 1), matrix(c(NA, df$fitness[2:3], df$fitness[5]), nrow = 2, byrow = F),
      phi = 10,
      theta = 55,
      xlab = '',
      ylab = '',
      zlab = '',
      zlim = c(min(df$fitness), max(df$fitness)),
      col = NA
)




### GET CONTRIBUTIONS: DF_J, EPS_IJ

getContributions <- function(mut_i, mut_j, save.plots = F) {
  
  # fitness effect
  df <- ge_data[ge_data$knock_in == mut_j, ]
  df <- df[!grepl(mut_i, df$background), ] # keep only backgrounds of i, B(i)
  df_mean <- aggregate(d_f~dose,
                       data = df[, c('dose', 'd_f')],
                       FUN = mean) # means for each dose
  df_wt <- df[df$background == '', ] # fitness effect of mutation in wild-type
  
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
    f_bg_ij <- aggregate(f~background,
                         data = f_bg_ij[, c('background', 'f')],
                         FUN = mean)
    f_bg_ij <- f_bg_ij[order(f_bg_ij$background), ]
    
    f <- merge(f_bg, f_bg_i, by = 'background', suffixes = c('_B', '_Bi'))
    f <- merge(f, f_bg_j, by = 'background')
    colnames(f)[4] <- 'f_Bj'
    f <- merge(f, f_bg_ij, by = 'background')
    colnames(f)[5] <- 'f_Bij'
    
    eps <- rbind(eps, data.frame(dose = d,
                                 background = f$background,
                                 f_bg = f$f_B,
                                 f_bg_i = f$f_Bi,
                                 f_bg_j = f$f_Bj,
                                 f_bg_ij = f$f_Bij,
                                 eps_ij = f$f_Bij - f$f_Bi - f$f_Bj + f$f_B))
    
  }
  
  eps_mean <- aggregate(eps_ij~dose,
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

# average across backgrounds
contrib_bgavg <- aggregate(cbind(dF_j, eps_ij) ~ dose + focal_mut + mut_j,
                           data = contrib,
                           FUN = mean)

# get beta_ij and w_ij from dF_j and eps_ij
sum_dF_j <- aggregate(dF_j ~ focal_mut + dose,
                      data = contrib_bgavg,
                      FUN = function(x) sum(x^2))
colnames(sum_dF_j)[3] <- 'sum_dF_j_squared'
contrib_bgavg <- merge(contrib_bgavg, sum_dF_j, by = c('dose', 'focal_mut'), all = T)

contrib_bgavg$w_ij <- contrib_bgavg$dF_j^2 / contrib_bgavg$sum_dF_j_squared
contrib_bgavg$beta_ij <- contrib_bgavg$eps_ij / contrib_bgavg$dF_j

contrib_bgavg$w_ij_times_beta_ij <- contrib_bgavg$w_ij * contrib_bgavg$beta_ij
contrib_bgavg$w_ij_times_beta_ij_squared <- contrib_bgavg$w_ij * contrib_bgavg$beta_ij^2

# get (weighted) mean, mean of squares, and CV of beta_ij
contrib_avg <- aggregate(cbind(w_ij_times_beta_ij, w_ij_times_beta_ij_squared) ~ dose + focal_mut,
                         data = contrib_bgavg,
                         FUN = sum)
contrib_avg$CVw_beta_ij <- sqrt(contrib_avg$w_ij_times_beta_ij_squared - contrib_avg$w_ij_times_beta_ij^2) / contrib_avg$w_ij_times_beta_ij
colnames(contrib_avg)[3:4] <- c('meanw_beta_ij', 'meanw_beta_ij_squared')

# compare with empirical values for slope, var df / var fb, and R^2
epist_coords <- merge(epist_coords, contrib_avg, by.x = c('dose', 'mut'), by.y = c('dose', 'focal_mut'))
epist_coords$dose_rank <- as.numeric(epist_coords$dose_rank)




### FIG 3A-E: EPISTASIS "COORDINATES" VS EFFECTIVE INTERACTIONS FOR FOCAL MUTATION (C59R)

contrib$dose_rank <- setNames(1:7, unique(contrib$dose))[contrib$dose]
contrib_bgavg$dose_rank <- setNames(1:7, unique(contrib_bgavg$dose))[contrib_bgavg$dose]

# dF_j vs. dose
ggplot(contrib[contrib$focal_mut == 'C59R', ],
       aes(x = dose_rank, y = dF_j, color = mut_j)) +
  geom_point() +
  geom_line(data = contrib_bgavg[contrib_bgavg$focal_mut == 'C59R', ]) +
  # geom_smooth(aes(fill = mut_j),
  #             method = 'loess',
  #             formula = y~x,
  #             alpha = 0.25) +
  # facet_wrap(~focal_mut,
  #            ncol = 1,
  #            scales = 'free_y',
  #            labeller = label_bquote(Mut.~italic(i)~'='~.(levels(contrib$focal_mut)[focal_mut]))) +
  scale_x_continuous(name = expression(paste('Drug dose (', mu, 'M)', sep = '')),
                     breaks = c(1, 3, 5, 7),
                     labels = mylabels[c(1,3,5,7)]) +
  scale_y_continuous(name = expression(symbol("\341")*Delta*italic(f)[italic(j)]*symbol("\361")[italic(B)(italic(i))]),
                     breaks = pretty_breaks(n = 2)) +
  scale_color_manual(name = expression(paste('Mut. ', italic(j), sep = '')),
                     values = mut_colors) +
  scale_fill_manual(name = expression(paste('Mut. ', italic(j), sep = '')),
                    values = mut_colors) +
  theme_bw() +
  theme(aspect.ratio = 0.6,
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0,
                                  size = 16),
        strip.text.y = element_text(angle = 0),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14))

ggsave(file = paste('../plots/df_', saveplot, '.pdf', sep = ''),
       device = 'pdf',
       width = 120,
       height = 120,
       units = 'mm')

# eps_ij vs. dose
ggplot(contrib[contrib$focal_mut == 'C59R', ],
       aes(x = dose_rank, y = eps_ij, color = mut_j)) +
  geom_point() +
  geom_line(data = contrib_bgavg[contrib_bgavg$focal_mut == 'C59R', ]) +
  # geom_smooth(aes(fill = mut_j),
  #             method = 'loess',
  #             formula = y~x,
  #             alpha = 0.25) +
  # facet_wrap(~focal_mut,
  #            ncol = 1,
  #            scales = 'free_y',
  #            labeller = label_bquote(Mut.~italic(i)~'='~.(levels(contrib$focal_mut)[focal_mut]))) +
  scale_x_continuous(name = expression(paste('Drug dose (', mu, 'M)', sep = '')),
                     breaks = c(1, 3, 5, 7),
                     labels = mylabels[c(1,3,5,7)]) +
  scale_y_continuous(name = expression(symbol("\341")*epsilon[italic(ij)]*symbol("\361")),
                     breaks = pretty_breaks(n = 2)) +
  scale_color_manual(name = expression(paste('Mut. ', italic(j), sep = '')),
                     values = mut_colors) +
  scale_fill_manual(name = expression(paste('Mut. ', italic(j), sep = '')),
                    values = mut_colors) +
  theme_bw() +
  theme(aspect.ratio = 0.6,
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0,
                                  size = 16),
        strip.text.y = element_text(angle = 0),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14))

ggsave(file = paste('../plots/eps_', saveplot, '.pdf', sep = ''),
       device = 'pdf',
       width = 120,
       height = 120,
       units = 'mm')

# distributions of beta vs. dose
contrib_avg$dose_rank <- setNames(c(1:7),
                                  unique(contrib_avg$dose))[contrib_avg$dose]
contrib_avg$beta_ij <- contrib_avg$meanw_beta_ij
contrib_avg$w_ij <- 1

ylimits <- c(-3, 3)
ggplot(contrib_bgavg[contrib_bgavg$focal_mut == 'C59R' & contrib_bgavg$beta_ij > ylimits[1] & contrib_bgavg$beta_ij < ylimits[2], ],
       aes(x = dose_rank, y = beta_ij, color = mut_j, group = dose_rank, cex = w_ij)) +
  geom_abline(slope = 0, intercept = 0, color = 'gray') +
  # geom_violin(trim = F,
  #             bw = 0.5,
  #             width = 0.5,
  #             color = NA,
  #             fill = 'lightgray',
  #             size = 0) +
  geom_point() +
  geom_point(data = contrib_avg[contrib_avg$focal_mut == 'C59R', ],
             aes(x = dose_rank, y = meanw_beta_ij),
             shape = 4,
             color = 'black') +
  geom_errorbar(data = contrib_avg[contrib_avg$focal_mut == 'C59R', ],
                aes(x = dose_rank, y = meanw_beta_ij, ymin = meanw_beta_ij - meanw_beta_ij*CVw_beta_ij, ymax = meanw_beta_ij + meanw_beta_ij*CVw_beta_ij),
                color = 'black',
                size = 1,
                width = 0) +
  scale_y_continuous(limits = ylimits,
                     name = expression(italic(beta)[italic(ij)]),
                     breaks = pretty_breaks(n = 3)) +
  scale_x_continuous(name = expression(paste('Drug dose (', mu, 'M)', sep = '')),
                     breaks = c(1, 3, 5, 7),
                     labels = mylabels[c(1,3,5,7)],
                     expand = c(0.1, 0.1)) +
  scale_color_manual(name = expression(paste('Mut. ', italic(j), sep = '')),
                     values = mut_colors) +
  theme_bw() +
  theme(aspect.ratio = 0.6,
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0,
                                  size = 16),
        strip.text.y = element_text(angle = 0),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14)) # this ggplot call spits a warning just because it's slightly trimming the edges of the violin plots (since ylimits are specified manually), it is not a problem

ggsave(file = paste('../plots/beta_distrib_vs_dose_', saveplot, '.pdf', sep = ''),
       device = 'pdf',
       width = 140,
       height = 120,
       units = 'mm')

# statistics of distributions: weighted mean, weighted second moment, and weighted CV
epist_coords$invCV <- 1 / (1 + epist_coords$CVw_beta_ij^2)

plot_this <- epist_coords[, c('mut', 'dose_rank', 'slope', 'var_rel', 'R2', 'meanw_beta_ij', 'meanw_beta_ij_squared', 'invCV')]
plot_this <- gather(plot_this, stat, value, slope:invCV)
plot_this$is_empirical <- plot_this$stat %in% c('slope', 'var_rel', 'R2')
plot_this$stat <- gsub('meanw_beta_ij_squared', 'var_rel', plot_this$stat)
plot_this$stat <- gsub('meanw_beta_ij', 'slope', plot_this$stat)
plot_this$stat <- gsub('invCV', 'R2', plot_this$stat)
plot_this$stat <- factor(plot_this$stat, levels = c('slope', 'var_rel', 'R2'))

plot_this$value[plot_this$stat == 'var_rel'] <- log10(plot_this$value[plot_this$stat == 'var_rel'])

ggplot(plot_this,
       aes(x = dose_rank, y = value, group = is_empirical, color = is_empirical)) +
  geom_line(data = plot_this[plot_this$mut == 'C59R', ]) +
  geom_blank(data = plot_this[plot_this$is_empirical == T, ]) + # to homogenize vertical scales with fig 2D-F
  facet_wrap(~stat,
             ncol = 1,
             scales = 'free_y') +
  scale_y_continuous(name = 'Value',
                     breaks = pretty_breaks(n = 3)) +
  scale_x_continuous(name = expression(paste('Drug dose (', mu, 'M)', sep = '')),
                     breaks = c(1, 3, 5, 7),
                     labels = mylabels[c(1,3,5,7)]) +
  scale_color_manual(values = c('gray', as.character(mut_colors['C59R']))) +
  theme_bw() +
  theme(aspect.ratio = 0.3,
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0,
                                  size = 16),
        strip.text.y = element_text(angle = 0),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        panel.border = element_blank(),
        legend.position = 'none') +
  annotate("segment", x=-Inf, xend=-Inf, y=Inf, yend=-Inf, size=0.5) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=0.5)

ggsave(file = paste('../plots/beta_stats_vs_dose_', saveplot, '.pdf', sep = ''),
       device = 'pdf',
       width = 100,
       height = 200,
       units = 'mm')



### FIG 4A: EPISTASIS MAP, EXPECTED VS OBSERVED

# expected segments
map_segments_expected <- data.frame(dose = numeric(0),
                                    mut = character(0),
                                    xstart = numeric(0),
                                    xend = numeric(0),
                                    ystart = numeric(0),
                                    yend = numeric(0))
for (d in 1:(length(unique(epist_coords$dose)) - 1)) {
  for (m in unique(epist_coords$mut)) {
    
    map_segments_expected <- rbind(map_segments_expected,
                                   data.frame(dose = unique(epist_coords$dose)[d],
                                              mut = m,
                                              xstart = epist_coords$meanw_beta_ij_squared[epist_coords$dose == unique(epist_coords$dose)[d] & epist_coords$mut == m],
                                              xend = epist_coords$meanw_beta_ij_squared[epist_coords$dose == unique(epist_coords$dose)[d + 1] & epist_coords$mut == m],
                                              ystart = epist_coords$invCV[epist_coords$dose == unique(epist_coords$dose)[d] & epist_coords$mut == m],
                                              yend = epist_coords$invCV[epist_coords$dose == unique(epist_coords$dose)[d + 1] & epist_coords$mut == m]))
    
  }
}

plot_this <- rbind(cbind(map_segments, is_empirical = T),
                   cbind(map_segments_expected, is_empirical = F))

# plot map
ggplot(plot_this,
       aes(x = xstart, y = ystart, xend = xend, yend = yend, color = mut, linetype = is_empirical)) +
  geom_arrowsegment(arrow_positions = 1,
                    arrows = arrow(angle = 25,
                                   type = 'open',
                                   length = unit(0.25, "cm")),
                    linewidth = 0.5) +
  scale_color_manual(values = mut_colors) +
  scale_linetype_manual(values = c('dashed', 'solid')) +
  scale_x_log10(name = expression(paste('var ', Delta, italic(f), ' / ', 'var ', italic(f), '(', italic(B), ')', sep = '')),
                breaks = 10^c(-2, 0, 2),
                labels = c(expression(10^-2),
                           expression(10^0),
                           expression(10^2))) +
  scale_y_continuous(name = expression(italic(R)^2),
                     breaks = c(0, 0.3, 0.6)) +
  theme_bw() +
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 16, vjust = 0),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14, hjust = 1),
        legend.position = 'none')

ggsave(file = paste('../plots/epistMap_alldoses_theory_vs_empirical_', saveplot, '.pdf', sep = ''),
       device = 'pdf',
       width = 100,
       height = 100,
       units = 'mm')



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

ggsave(file = paste('../plots/slope_vs_beta_', saveplot, '.pdf', sep = ''),
       device = 'pdf',
       width = 120,
       height = 120,
       units = 'mm')

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

ggsave(file = paste('../plots/var_vs_beta_', saveplot, '.pdf', sep = ''),
       device = 'pdf',
       width = 125,
       height = 125,
       units = 'mm')

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

ggsave(file = paste('../plots/R2_vs_beta_', saveplot, '.pdf', sep = ''),
       device = 'pdf',
       width = 125,
       height = 125,
       units = 'mm')



### LINEARITY TESTS

ltests <- data.frame(dose = character(0),
                     focal_mut = character(0),
                     p_non_normal = numeric(0),
                     p_heteroscedasticity = numeric(0),
                     p_rescorr = numeric(0),
                     p_cubic = numeric(0))

for (d in unique(ge_data$dose)) {
  for (m in unique(ge_data$knock_in)) {
    
    df <- ge_data[ge_data$dose == d & ge_data$knock_in == m, ]
    
    linear_model <- lm(d_f ~ background_f, df)
    
    # are residuals normally distributed?
    p_non_normal <- shapiro.test(linear_model$residuals)$p.value
    
    # are residuals homogeneously spread?
    p_heteroscedasticity <- bptest(linear_model)$p.value
    
    # are residuals uncorrelated?
    res_model <- lm(linear_model$residuals[-nrow(df)] ~ linear_model$residuals[-1])
    p_rescorr <- summary(res_model)$coefficients[2,4]
    
    # is a cubic model significantly better than the linear model?
    cubic_model <- lm(d_f ~ poly(background_f, 3), df)
    anova(linear_model, cubic_model)
    p_cubic <- anova(linear_model, cubic_model)$`Pr(>F)`[2]
    
    ltests <- rbind(ltests,
                    data.frame(dose = d,
                               focal_mut = m,
                               p_non_normal = p_non_normal,
                               p_heteroscedasticity = p_heteroscedasticity,
                               p_rescorr = p_rescorr,
                               p_cubic = p_cubic))
    
  }
}

linear_model_ok <- (sum(ltests$p_non_normal > 0.05 & ltests$p_heteroscedasticity > 0.05 & ltests$p_rescorr > 0.05))/nrow(ltests)
cubic_model_better <- sum(ltests$p_cubic < 0.05)/nrow(ltests)







### SUPPLEMENTARY FIGURES

# all global epistasis patterns

plot_this <- NULL # some tweaking to homogenize y scales across rows
for (focal_mut in unique(ge_data$knock_in)) {
  rowlimits <- c(min(ge_data$d_f[ge_data$knock_in == focal_mut]),
                 max(ge_data$d_f[ge_data$knock_in == focal_mut]))
  mean_x <- aggregate(background_f ~ dose,
                      data = ge_data[ge_data$knock_in == focal_mut, ],
                      FUN = mean)
  plot_this <- rbind(plot_this,
                     data.frame(dose = mean_x$dose,
                                knock_in = focal_mut,
                                background_f = mean_x$background_f,
                                d_f = rep(rowlimits, length(unique(ge_data$dose)))))
}
plot_this$slope <- 0

ggplot(ge_data, aes(x = background_f, y = d_f, color = slope)) +
  geom_abline(slope = 0, intercept = 0, color = '#d1d3d4') +
  geom_point(color = 'black') +
  geom_blank(data = plot_this) +
  facet_wrap(knock_in ~ dose,
             nrow = length(mut_colors),
             scales = 'free') +
  geom_smooth(method = 'lm',
              se = F,
              formula = y~x,
              fullrange = T) +
  scale_x_continuous(name = 'F (genetic background)',
                     breaks = pretty_breaks(n = 2)) +
  scale_y_continuous(name = expression(paste(Delta, 'F', sep = '')),
                     breaks = pretty_breaks(n = 2)) +
  scale_color_gradient2(name = 'slope',
                        low = '#1e9fd3',
                        high = '#ef3a37',
                        mid = 'white',
                        midpoint = 0,
                        breaks = c(-0.4, 0, 0.4),
                        limits = c(-0.8, 0.8)) +
  theme_bw() +
  theme(aspect.ratio = 0.6,
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        strip.text.y = element_text(angle = 0),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14, hjust = 1)) +
  guides(color = guide_colorbar(ticks.colour = NA))

ggsave(file = paste('../plots/FEEs_', saveplot, '.pdf', sep = ''),
       device = 'pdf',
       width = 380,
       height = 200,
       units = 'mm')

# bootstrap random landscapes

n_boot <- 500
random_landscapes <- data.frame(id = character(0),
                                dose = character(0),
                                knock_in = character(0),
                                slope = numeric(0),
                                R2 = numeric(0),
                                var_rel = numeric(0))

for (d in 1:length(unique(ge_data$dose))) {
  
  df.rnd <- cbind(genot_matrix, fun = data[, 2 + d])
  
  for (n in 1:n_boot) {
    df.rnd$fun <- sample(df.rnd$fun)
    ge_data.rnd.i <- makeGEdata(matrix2string(df.rnd))
    
    for (m in unique(ge_data.rnd.i$knock_in)) {
      lmod <- lm(d_f ~ background_f, data = ge_data.rnd.i[ge_data.rnd.i$knock_in == m, ])
      random_landscapes <- rbind(random_landscapes,
                                 data.frame(id = paste('rnd', d, n, sep = '.'),
                                            dose = unique(ge_data$dose)[d],
                                            knock_in = m,
                                            slope = as.numeric(lmod$coefficients[2]),
                                            R2 = summary(lmod)$r.squared,
                                            var_rel = var(ge_data.rnd.i$d_f[ge_data.rnd.i$knock_in == m])/var(ge_data.rnd.i$background_f[ge_data.rnd.i$knock_in == m])))
    }
    
  }
  
}

ggplot(random_landscapes,
       aes(x = log10(var_rel), y = R2)) +
  stat_density_2d(aes(fill = after_stat(density)),
                  geom = 'raster',
                  contour = FALSE) +
  geom_arrowsegment(data = map_segments,
                    aes(x = log10(xstart), xend = log10(xend), y = ystart, yend = yend, color = mut),
                    arrow_positions = 1,
                    arrows = arrow(angle = 25,
                                   type = 'open',
                                   length = unit(0.25, "cm")),
                    linewidth = 0.5) +
  geom_blank(data = epist_coords) +
  scale_x_continuous(name = expression(paste('var ', Delta, italic(f), ' / ', 'var ', italic(f), '(', italic(B), ')', sep = '')),
                     breaks = c(-2, 0, 2),
                     labels = c(expression(10^-2),
                                expression(10^0),
                                expression(10^2)),
                     expand = c(0,0)) +
  scale_y_continuous(name = expression(italic(R)^2),
                     breaks = c(0, 0.3, 0.6, 0.9),
                     expand = c(0,0)) +
  scale_fill_gradient2(low = 'white', high = 'grey50') +
  scale_color_manual(values = mut_colors) +
  theme_bw() +
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 16, vjust = 0),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14, hjust = 1),
        legend.position = 'none')

ggsave(file = paste('../plots/epistMap_randomLandscapes_', saveplot, '.pdf', sep = ''),
       device = 'pdf',
       width = 125,
       height = 125,
       units = 'mm')













