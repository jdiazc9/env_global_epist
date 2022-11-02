source('ecoFunctions.R')
library(RColorBrewer)

# Plasmodium falciparum carrying subsets of 4 mutations: N51I, C59R, S108N, I164L (denoted as N, C, S and I)
muts <- c('C', 'I', 'N', 'S')
mut_colors <- setNames(c('#e76d57', '#acd152', '#e7da5c', '#7a9dd2'),
                       muts)

# fitness is additive (delta_i, i = C, I, N or S) + two pairwise interactions: C-N and C-S (eps_CN and eps_CS)

# fitness effects scale with environmental variable lambda (L) that goes from -1 (10^-1 uM of drug) to 3 (10^3 uM of drug)
# sigmoid <- function(x) 1/(1+exp(-x))
# 
# delta <- setNames(c(function(L) 0.00,
#                     function(L) -0.1,
#                     function(L) -1.05 + sigmoid(1.5*L),
#                     function(L) 0.05 + sigmoid(1.5*L)),
#                   muts)
# 
# eps <- setNames(c(function(L) 0.00,
#                   function(L) 0.00,
#                   function(L) 1.05 - sigmoid(1.5*L),
#                   function(L) 0.05 + sigmoid(1.5*L)),
#                 muts)

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

# nicer labels for plots
mylabels <- c(expression(10^-2),
              expression(10^-1),
              expression(10^0),
              expression(10^1),
              expression(10^2))

# plot deltas
ggplot(delta.df[delta.df$mut != 'C', ], aes(x = L, y = delta, color = mut)) +
  geom_line() +
  scale_color_manual(name = expression(paste('Mut. ', italic(i), sep = '')),
                     values = mut_colors[names(mut_colors) != 'C']) +
  scale_x_continuous(name = expression(paste('Drug dose (', mu, 'M)', sep = '')),
                     breaks = c(-2, 0, 2),
                     labels = mylabels[c(1,3,5)]) +
  scale_y_continuous(name = expression(italic(delta)[italic(i)]),
                     breaks = c(-0.8, 0, 0.8)) +
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

ggsave(file = '../plots/minModel_delta.pdf',
       device = 'pdf',
       width = 100,
       height = 80,
       units = 'mm')

# plot epsilons
ggplot(eps.df[eps.df$mut != 'C', ], aes(x = L, y = eps, color = mut)) +
  geom_line() +
  scale_color_manual(name = expression(paste('Mut. ', italic(i), sep = '')),
                     values = mut_colors[names(mut_colors) != 'C']) +
  scale_x_continuous(name = expression(paste('Drug dose (', mu, 'M)', sep = '')),
                     breaks = c(-2, 0, 2),
                     labels = mylabels[c(1,3,5)]) +
  scale_y_continuous(name = expression(italic(epsilon)[C*italic(i)]),
                     breaks = c(-0.8, 0, 0.8),
                     labels = c('-0.8', '-0.0', '0.8')) +
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

ggsave(file = '../plots/minModel_eps.pdf',
       device = 'pdf',
       width = 100,
       height = 80,
       units = 'mm')

# plot weighted product

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

ggplot(contrib[contrib$mut != 'C', ], aes(x = L, y = term, color = mut)) +
  geom_line() +
  geom_line(data = contrib_sum,
            linetype = 'dashed',
            color = 'black') +
  scale_color_manual(name = expression(paste('Mut. ', italic(i), sep = '')),
                     values = mut_colors[names(mut_colors) != 'C']) +
  scale_x_continuous(name = expression(paste('Drug dose (', mu, 'M)', sep = '')),
                     breaks = c(-2, 0, 2),
                     labels = mylabels[c(1,3,5)]) +
  scale_y_continuous(name = 'Effective interaction with mut. C',
                     breaks = c(-0.8, 0, 0.8),
                     labels = c('-0.8', '-0.0', '0.8')) +
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

ggsave(file = '../plots/minModel_eff-inter.pdf',
       device = 'pdf',
       width = 100,
       height = 80,
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

# plot for dose -1, 0, 1, 2, 3

plot_these_doses <- c(which.min(abs(L - (-2))),
                      which.min(abs(L - (-1))),
                      which.min(abs(L - (0))),
                      which.min(abs(L - (1))),
                      which.min(abs(L - (2))))

plot_this <- gedf[gedf$dose %in% L[plot_these_doses] & gedf$knock_in == 'C', ]
plot_this$dose <- as.character(plot_this$dose)
plot_this$dose <- setNames(paste('1e', -1:3, sep = ''),
                           unique(plot_this$dose))[plot_this$dose]
plot_this$dose <- factor(plot_this$dose,
                         levels = paste('1e', -1:3, sep = ''))

plot_this$which_mut <- sapply(plot_this$background, FUN = function(x) grepl('S', x) + grepl('N', x))
plot_this$which_mut[plot_this$which_mut == 1 & grepl('S', plot_this$background)] <- 'S'
plot_this$which_mut[plot_this$which_mut == 1 & grepl('N', plot_this$background)] <- 'N'
plot_this$which_mut <- as.character(plot_this$which_mut)
plot_this$which_mut <- factor(plot_this$which_mut,
                              levels = c('2', 'N', 'S', '0'))

ggplot(plot_this,
       aes(x = background_f, y = d_f, color = slope, shape = which_mut, group = dose)) +
  geom_abline(slope = 0, intercept = 0, color = '#d1d3d4') +
  geom_point(color = 'black',
             cex = 3) +
  geom_smooth(method = 'lm',
              formula = y~x,
              fullrange = T,
              se = F) +
  facet_wrap(~dose, nrow = 1, scales = 'free_x',
             labeller = label_bquote(.(setNames(mylabels,
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
  scale_shape_manual(name = 'Non-focal mutations\npresent in background',
                     values = c(19, 15, 17, 1),
                     labels = c('Both N and S',
                                'N but not S',
                                'S but not N',
                                'None of the two')) +
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

ggsave(file = '../plots/minModel_FEEs.pdf',
       device = 'pdf',
       width = 320,
       height = 150,
       units = 'mm')


