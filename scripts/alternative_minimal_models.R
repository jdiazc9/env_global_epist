source('ecoFunctions.R')
library(RColorBrewer)

# Plasmodium falciparum carrying subsets of 4 mutations: N51I, C59R, S108N, I164L (denoted as N, C, S and I)
muts <- c('C', 'I', 'N', 'S')
mut_colors <- setNames(c('#e76d57', '#acd152', '#e7da5c', '#7a9dd2'),
                       muts)

# fitness is additive (delta_i, i = C, I, N or S) + pairwise interactions: C-I, C-N and C-S

# fitness effects scale with environmental variable lambda (L) that goes from -1 (10^-2 uM of drug) to 2 (10^2 uM of drug)
# sigmoid <- function(x) 1/(1+exp(-x))

model_index <- 2 # choose between alternative models 1 to 3 (0 means original model in main text)
if (model_index == 0) {
  
  # original model
  delta <- setNames(c(function(L) 0.00,
                      function(L) -0.1,
                      function(L) -1 + 10^L/(1 + 10^L),
                      function(L) 10^L/(1 + 10^L)),
                    muts)
  
  eps <- setNames(c(function(L) 0.00,
                    function(L) 0.00,
                    function(L) 1 - 10^L/(1 + 10^L),
                    function(L) 10^L/(1 + 10^L)),
                  muts)


} else if (model_index == 1) {
  
  # alternative model 1
  delta <- setNames(c(function(L) 0.00,
                      function(L) -0.1,
                      function(L) -1,
                      function(L) 1),
                    muts)
  
  eps <- setNames(c(function(L) 0.00,
                    function(L) 0.05,
                    function(L) -0.05,
                    function(L) -1 + 2*10^L/(1 + 10^L)),
                  muts)
  
} else if (model_index == 2) {
  
  delta <- setNames(c(function(L) 0.00,
                      function(L) 2.5 - 2*10^L/(1 + 10^L),
                      function(L) -1.5 + 1.5*10^L/(1 + 10^L),
                      function(L) -2 - 0.5*10^L/(1 + 10^L)),
                    muts)
  
  eps <- setNames(c(function(L) 0.00,
                    function(L) -2 + 4*10^L/(1 + 10^L),
                    function(L) -0.5 + 10^L/(1 + 10^L),
                    function(L) 1 - 2*10^L/(1 + 10^L)),
                  muts)
  
}

# initialize array of environmental parameter values and output data frames
L <- seq(-2, 2, length.out = 100)

delta.df <- data.frame(mut = character(0),
                       L = numeric(0),
                       delta = numeric(0))
for (m in muts) {
  delta.df <- rbind(delta.df,
                    data.frame(mut = m,
                               L = L,
                               delta = delta[[m]](L)))
}

eps.df <- data.frame(mut = character(0),
                     L = numeric(0),
                     eps = numeric(0))
for (m in muts) {
  eps.df <- rbind(eps.df,
                  data.frame(mut = m,
                             L = L,
                             eps = eps[[m]](L)))
}

# weighted product
contrib <- merge(delta.df, eps.df)
contrib$term <- contrib$delta*contrib$eps

contrib_sum <- aggregate(formula = term ~ L,
                         data = contrib,
                         FUN = sum)

weights <- aggregate(formula = delta ~ L,
                     data = delta.df[delta.df$mut != 'C', ],
                     FUN = function(x) sum(x^2))
colnames(weights)[2] <- 'weight'

contrib_sum$term <- contrib_sum$term / weights$weight
contrib <- merge(contrib, weights, all = T)
contrib$term <- contrib$term / contrib$weight

# nicer labels for plots
mylabels <- c(expression(10^-2),
              expression(10^-1),
              expression(10^0),
              expression(10^1),
              expression(10^2))

# plot everything in one plot

plot_this_delta <- cbind(variable = 'delta',
                         delta.df[delta.df$mut != 'C', ])
colnames(plot_this_delta) <- c('variable', 'mut', 'L', 'value')

plot_this_eps <- cbind(variable = 'eps',
                       eps.df[eps.df$mut != 'C', ])
colnames(plot_this_eps) <- c('variable', 'mut', 'L', 'value')

plot_this_contrib <- contrib[contrib$mut != 'C', c('mut', 'L', 'term')]
plot_this_contrib <- rbind(plot_this_contrib,
                           data.frame(mut = 'Sum',
                                      L = contrib_sum$L,
                                      term = contrib_sum$term))
plot_this_contrib <- cbind(variable = 'contrib',
                           plot_this_contrib)
colnames(plot_this_contrib) <- c('variable', 'mut', 'L', 'value')

plot_this <- rbind(plot_this_delta, plot_this_eps, plot_this_contrib)
plot_this$variable <- factor(plot_this$variable, levels = c('delta', 'eps', 'contrib'))

ggplot(plot_this, aes(x = L, y = value, color = mut, group = mut, linetype = mut)) +
  geom_line() +
  facet_wrap(~variable,
             nrow = 1,
             scales = 'free_y',
             strip.position = 'left',
             labeller = labeller(variable = setNames(c('delta_i',
                                                       'epsilon_Ci',
                                                       'Effective\ninteraction'),
                                                     levels(plot_this$variable)),
                                 type = label_parsed)) +
  scale_color_manual(name = expression(paste('Mut. ', italic(i), sep = '')),
                     values = c(mut_colors, Sum = 'black')) +
  scale_linetype_manual(values = setNames(c(rep('solid', 4), 'dashed'),
                        c(muts, 'Sum'))) +
  scale_x_continuous(name = expression(paste('Drug dose (', mu, 'M)', sep = '')),
                     breaks = c(-1, 1),
                     labels = mylabels[c(2,4)]) +
  scale_y_continuous(name = 'Value',
                     breaks = pretty_breaks(n = 3)) +
  theme_bw() +
  theme(aspect.ratio = 0.6,
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 18),
        strip.placement = 'outside',
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        axis.title.y = element_blank(),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14))

ggsave(file = paste('../plots/minModel_params_', model_index, '.pdf', sep = ''),
       device = 'pdf',
       width = 250,
       height = 150,
       units = 'mm')
  

# function to generate fitness of a genotype by name

getFitness <- function(genot, dose = -1) {
  
  if (nchar(genot) == 0) {
    fitness <- 1
  } else {
    
    genot <- strsplit(genot, split = ',')[[1]]
    
    # additive component
    fitness <- 1 + sum(sapply(genot, FUN = function(m) delta[[m]](dose)))
    
    # interactions with mut. C (if present, and if there is at least one other mutation)
    if ('C' %in% genot & length(genot) > 1) {
      fitness <- fitness + sum(sapply(genot[genot != 'C'], FUN = function(m) eps[[m]](dose)))
    }
    
  }
  
  return(fitness)
  
}

# generate names of all genotypes
genots <- c('', paste(muts, collapse = ','))
for (nmut in 1:3) {
  
  comb <- combn(4, nmut)
  
  for (i in 1:ncol(comb)) genots <- c(genots,
                                      paste(muts[comb[, i]], collapse = ','))
  
}

# generate fitness of every genotype in a concentration gradient of drug
df <- data.frame(dose = numeric(0),
                 genot = character(0),
                 fitness = numeric(0))
for (d in L) {
  
  df <- rbind(df,
              data.frame(dose = d,
                         genot = genots,
                         fitness = sapply(genots, getFitness, dose = d)))
  
}

# make global epistasis data
gedf <- data.frame(dose = numeric(0),
                   background = character(0),
                   knock_in = character(0),
                   background_f = numeric(0),
                   d_f = numeric(0))
for (d in L) {
  gedf <- rbind(gedf,
                cbind(dose = d,
                      makeGEdata(df[df$dose == d, c('genot', 'fitness')])))
}

# get slopes
slope <- NULL
for (d in L) {
  mylm <- lm(formula = d_f~background_f,
             data = gedf[gedf$knock_in == 'C' & gedf$dose == d, ])
  slope <- c(slope,
             summary(mylm)$coefficients[2])
}

gedf <- merge(gedf,
              data.frame(dose = L, slope = slope))

# plot for dose -1, 0, 1

plot_these_doses <- c(which.min(abs(L - (-1))),
                      which.min(abs(L - (0))),
                      which.min(abs(L - (1))))

plot_this <- gedf[gedf$dose %in% L[plot_these_doses] & gedf$knock_in == 'C', ]
plot_this$dose <- as.character(plot_this$dose)
plot_this$dose <- setNames(paste('1e', -1:1, sep = ''),
                           unique(plot_this$dose))[plot_this$dose]
plot_this$dose <- factor(plot_this$dose,
                         levels = paste('1e', -1:1, sep = ''))

plot_this$which_mut <- sapply(plot_this$background, FUN = function(x) grepl('S', x) + grepl('N', x))
plot_this$which_mut[plot_this$which_mut == 1 & grepl('S', plot_this$background)] <- 'S'
plot_this$which_mut[plot_this$which_mut == 1 & grepl('N', plot_this$background)] <- 'N'
plot_this$which_mut <- as.character(plot_this$which_mut)
plot_this$which_mut <- factor(plot_this$which_mut,
                              levels = c('2', 'N', 'S', '0'))

shape_labels <- sapply(unique(plot_this$background),
                       FUN = function(genot) {
                         genot <- strsplit(genot, split = ',')[[1]]
                         return(paste(as.numeric(muts %in% genot), collapse = ''))
                       })

ggplot(plot_this,
       aes(x = background_f, y = d_f, color = slope, shape = background, group = dose)) +
  geom_abline(slope = 0, intercept = 0, color = '#d1d3d4') +
  geom_point(color = 'black',
             cex = 3) +
  geom_smooth(method = 'lm',
              formula = y~x,
              fullrange = T,
              se = F) +
  facet_wrap(~dose,
             nrow = 1,
             scales = 'free_x',
             labeller = label_bquote(.(setNames(mylabels[2:4],
                                                levels(plot_this$dose))[dose]))) +
  scale_x_continuous(name = expression(paste('Fitness of genetic background, ', italic(f), '(', italic(B), ')', sep = '')),
                     breaks = pretty_breaks(n = 2),
                     expand = rep(0.1, 2)) +
  scale_y_continuous(name = expression(paste('Fitness effect\nof mut. C, ', Delta, italic(f), sep = '')),
                     breaks = pretty_breaks(n = 3),
                     expand = rep(0.1, 2)) +
  scale_color_gradient2(low = '#1e9fd3',
                        high = '#ef3a37',
                        mid = 'white',
                        midpoint = 0,
                        breaks = pretty_breaks(n = 3)) +
  scale_shape_manual(name = 'Background genotype',
                     values = c(1, 3, 16, 0, 4, 15, 2, 17),
                     labels = shape_labels) +
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

ggsave(file = paste('../plots/minModel_FEEs_', model_index, '.pdf', sep = ''),
       device = 'pdf',
       width = 250,
       height = 150,
       units = 'mm')


