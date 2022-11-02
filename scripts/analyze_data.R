source('ecoFunctions.R')
library(RColorBrewer)

# Plasmodium falciparum carrying subsets of 4 mutations: N51I, C59R, S108N, I164L

file <- '../data/growth_rates_pyr.txt'
file <- '../data/growth_rates_cyc.txt' # FIXME: uncomment for cycloguanil 
data <- read.table(file)
data[, 2] <- sprintf("%04d", data[, 2])
doses <- c(0, 10^(-2:6))
colnames(data) <- c('genot_id', 'genot_bin', paste('dose_', 0:(length(doses) - 1), sep = ''))

genot_matrix <- do.call(rbind,
                        lapply(data$genot_bin,
                               FUN = function(genot) {
                                 strsplit(genot, split = '')[[1]]
                               }))
colnames(genot_matrix) <- c('N51I', 'C59R', 'S108N', 'I164L')
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

### MAKE PLOTS
if (grepl('cyc', file)) {
  saveplot <- 'cyc'
} else if (grepl('pyr', file)) {
  saveplot <- 'pyr'
}

# nicer labels for plots
mylabels <- c(0,
              expression(10^-2),
              expression(10^-1),
              expression(10^0),
              expression(10^1),
              expression(10^2),
              expression(10^3),
              expression(10^4),
              expression(10^5),
              expression(10^6))

# nicer colors to identify mutations
mut_colors <- setNames(c('#e76d57', '#acd152', '#e7da5c', '#7a9dd2'),
                       c("C59R", "I164L", "N51I", "S108N"))

# some tweaking to homogenize y scales across rows
plot_this <- NULL
for (focal_mut in unique(ge_data$knock_in)) {
  rowlimits <- c(min(ge_data$d_f[ge_data$knock_in == focal_mut]),
                 max(ge_data$d_f[ge_data$knock_in == focal_mut]))
  mean_x <- aggregate(formula = background_f ~ dose,
                      data = ge_data[ge_data$knock_in == focal_mut, ],
                      FUN = mean)
  plot_this <- rbind(plot_this,
                     data.frame(dose = mean_x$dose,
                                knock_in = focal_mut,
                                background_f = mean_x$background_f,
                                d_f = rep(rowlimits, length(unique(ge_data$dose)))))
}

# plot FEEs
myplot <- 
ggplot(ge_data, aes(x = background_f, y = d_f, color = knock_in, group = dose)) +
  geom_abline(slope = 0, intercept = 0, color = '#d1d3d4') +
  geom_point(color = 'black') +
  geom_point(data = plot_this,
             color = 'white') +
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
                     #breaks = c(-1, 0, 1), labels = c('-1.0', '0.0', '1.0')) +
  scale_color_manual(name = 'Mutation', values = mut_colors) +
  theme_bw() +
  theme(aspect.ratio = 0.6,
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        strip.text.y = element_text(angle = 0),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16))

print(myplot)
ggsave(myplot,
       file = paste('../plots/FEEs_', saveplot, '.pdf', sep = ''),
       device = 'pdf',
       width = 380,
       height = 200,
       units = 'mm')

# get slopes
slopes <- data.frame(dose = character(0),
                     mut = character(0),
                     slope = numeric(0),
                     pval = numeric(0))

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
                               pval = p))
    
  }
}

slopes$dose_rank <- setNames(1:length(unique(slopes$dose)),
                             unique(slopes$dose))[slopes$dose]

myplot <- 
  ggplot(slopes, aes(x = dose_rank, y = slope, color = mut)) +
    geom_line(size = 0.75) +
    scale_x_continuous(name = expression(paste('Concentration (', mu, 'M)', sep = '')),
                       breaks = unique(slopes$dose_rank),
                       labels = unique(slopes$dose)) +
    scale_y_continuous(name = 'Global epistasis slope') +
    scale_color_brewer(name = 'Mutation', palette = 'Set2') +
    theme_bw() +
    theme(aspect.ratio = 0.6,
          panel.grid = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(hjust = 0),
          axis.title = element_text(size = 18),
          axis.text = element_text(size = 16),
          legend.title = element_text(size = 18),
          legend.text = element_text(size = 16))

print(myplot)
ggsave(myplot,
       file = paste('../plots/rn_', saveplot, '.pdf', sep = ''),
       device = 'pdf',
       width = 150,
       height = 150,
       units = 'mm')

### FEE OF FOCAL MUTATION (C59R) WITH NO DRUG

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

ggsave(file = paste('../plots/FEEs_C59R_dose0_', saveplot, '.pdf', sep = ''),
       device = 'pdf',
       width = 100,
       height = 100,
       units = 'mm')

### FEES AS A FUNCTION OF DRUG DOSE FOR FOCAL MUTATION (C59R)

ge_data$dose <- factor(ge_data$dose, levels = c('0', paste('1e', -2:3, sep = '')))
ge_data <- merge(ge_data, slopes, by.x = c('dose', 'knock_in'), by.y = c('dose', 'mut'), all = T)

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
  theme(aspect.ratio = 0.6,
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
       width = 350,
       height = 150,
       units = 'mm')

# same plot (only for doses 1e-2 to 1e3), dot shapes indicate mutations presence/absence

plot_this <- ge_data[ge_data$knock_in == 'C59R' & ge_data$dose %in% c('1e-2', '1e-1', '1e0', '1e1', '1e2'), ]
plot_this$which_mut <- grepl('N51I', plot_this$background) + grepl('S108N', plot_this$background)
plot_this$which_mut[plot_this$which_mut == 1 & grepl('N51I', plot_this$background)] <- 'N'
plot_this$which_mut[plot_this$which_mut == 1 & grepl('S108N', plot_this$background)] <- 'S'
plot_this$which_mut <- as.character(plot_this$which_mut)
plot_this$which_mut <- factor(plot_this$which_mut, levels = c('2', 'N', 'S', '0'))

ggplot(plot_this,
       aes(x = background_f, y = d_f, color = slope, group = dose, shape = which_mut)) +
  geom_abline(slope = 0, intercept = 0, color = '#d1d3d4') +
  geom_point(color = 'black',
             cex = 3) +
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

ggsave(file = paste('../plots/FEEs_C59R_which-mut_', saveplot, '.pdf', sep = ''),
       device = 'pdf',
       width = 320,
       height = 150,
       units = 'mm')

### FITNESS GRAPH DEMONSTRATING EPISTASIS

# tunable parameters
n_mut <- 4
mycolors <- c('#939598', '#d68f28', '#415ba9', '#a96cad')

# make genotype names
genots <- lapply(0:n_mut, FUN = function(i) t(combn(n_mut, i)))
genots <- lapply(genots, FUN = function(x) sapply(1:nrow(x),
                                                  FUN = function(i) paste(x[i, ], collapse = ',')))
genots <- unlist(genots)

# make edges of fitness graph
nMut <- function(genot) sapply(genot,
                               FUN = function(genot_i) length(strsplit(genot_i, split = ',')[[1]]))

isDescendant <- function(this_genot, of_this_genot) {
  
  this_genot <- strsplit(this_genot, split = ',')[[1]]
  of_this_genot <- strsplit(of_this_genot, split = ',')[[1]]
  
  return(all(of_this_genot %in% this_genot))
  
}

edges <- data.frame(source = character(0),
                    target = character(0),
                    source.nmut = numeric(0),
                    target.nmut = numeric(0))

for(s in genots) {
  
  t <- genots[sapply(genots,
                     isDescendant,
                     of_this_genot = s) & nMut(genots) == nMut(s)+1]
  if(length(t)) {
    edges <- rbind(edges,
                   data.frame(source = s,
                              target = as.character(t),
                              source.nmut = as.numeric(nMut(s)),
                              target.nmut = as.numeric(nMut(s)) + 1))
  }
  
}

edges <- cbind(edge_id = paste('edge_', 1:nrow(edges), sep = ''),
               edges)

# plot landscape (wrapper function)
plotGraph <- function(landscape, save.plot = F) {
  
  df <- cbind(edges,
              source.f = setNames(landscape$f, landscape$genot)[edges$source],
              target.f = setNames(landscape$f, landscape$genot)[edges$target])
  df$source.f[is.na(df$source.f)] <- landscape$f[landscape$genot == '']
  
  if ('color' %in% colnames(landscape)) {
    df <- merge(df, landscape[, c('genot', 'color')], by.x = 'target', by.y = 'genot')
  } else {
    df$color <- 'A'
  }
  df <- df[, c('edge_id', 'source', 'target', 'source.nmut', 'target.nmut', 'source.f', 'target.f', 'color')]
  
  dfx <- gather(df[, c(1, 4, 5)], position, nmut, source.nmut:target.nmut)
  dfx$position <- setNames(c('source', 'target'), c('source.nmut', 'target.nmut'))[dfx$position]
  
  dfy <- gather(df[, c(1, 6, 7)], position, f, source.f:target.f)
  dfy$position <- setNames(c('source', 'target'), c('source.f', 'target.f'))[dfy$position]
  
  dfxy <- merge(dfx, dfy, by = c('edge_id', 'position'))
  
  df <- merge(dfxy, df[, c('edge_id', 'color')], by = 'edge_id')
  
  dy <- min(c(max(landscape$f) - landscape$f[1], landscape$f[1] - min(landscape$f)))
  dy <- round(dy/0.1)*0.1
  ybreaks <- seq(landscape$f[1] - 10*dy, landscape$f[1] + 10*dy, by = dy)
  
  myplot <-
    ggplot(df, aes(x = nmut, y = f, group = edge_id, color = color)) +
    geom_line() +
    scale_x_continuous(name = '# mutations',
                       breaks = 0:n_mut,
                       labels = as.character(0:n_mut)) +
    scale_y_continuous(name = 'Fitness',
                       breaks = pretty_breaks(n = 3),
                       expand = c(0.05, 0.05)) +
    scale_color_manual(values = setNames(mycolors, LETTERS[1:length(mycolors)])) +
    theme_bw() +
    theme(aspect.ratio = 0.6,
          panel.grid = element_blank(),
          panel.border = element_blank(),
          legend.position = 'none',
          axis.title = element_text(size = 18),
          axis.text = element_text(size = 16)) +
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf,size=0.5) +
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf,size=0.5)
  
  if (save.plot != F) {
    ggsave(myplot,
           file = paste('../plots/', save.plot, '.pdf', sep = ''),
           dpi = 600,
           width = 100,
           height = 80,
           units = 'mm')
  }
  
  return(myplot)
  
}

# full fitnesss graph

df <- data.frame(genot = orderName(paste(ge_data$background[ge_data$dose == '0'],
                                         ge_data$knock_in[ge_data$dose == '0'],
                                         sep = ',')),
                 f = ge_data$background_f[ge_data$dose == '0'] + ge_data$d_f[ge_data$dose == '0'])
df <- aggregate(formula = f~genot,
                data = df,
                FUN = mean)
for (i in 1:nrow(df)) if (substr(df$genot[i], 1, 1) == ',') df$genot[i] <- substr(df$genot[i], 2, nchar(df$genot[i]))
df <- rbind(df, data.frame(genot = '', f = data$dose_0[data$genot_bin == '0000']))

landscape <- df
landscape$genot <- gsub('C59R', '1', landscape$genot)
landscape$genot <- gsub('I164L', '2', landscape$genot)
landscape$genot <- gsub('N51I', '3', landscape$genot)
landscape$genot <- gsub('S108N', '4', landscape$genot)

plotGraph(landscape)

ggsave(file = paste('../plots/fitnessGraph-full_', saveplot, '.pdf', sep = ''),
       device = 'pdf',
       width = 100,
       height = 100,
       units = 'mm')

# inset

df <- df[df$genot %in% c('N51I', 'I164L,N51I', 'C59R,N51I', 'C59R,I164L,N51I'), ]
df$nmut <- sapply(df$genot, FUN = function(x) length(strsplit(x, split = ',')[[1]]))

ggplot(df, aes(x = nmut, y = f)) +
  geom_point(cex = 3) +
  scale_y_continuous(name = 'Fitness',
                     limits = c(min(df$f), max(c(df$f[df$nmut == max(df$nmut)],
                                                 sum(df$f[df$nmut == (max(df$nmut) - 1)]) - df$f[df$nmut == min(df$nmut)]))),
                     expand = rep(0.05, 2),
                     breaks = pretty_breaks(n = 3)) +
  scale_x_continuous(expand = rep(0.1, 2)) +
  theme_bw() +
  theme(aspect.ratio = 0.6,
        panel.grid = element_blank(),
        panel.border = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 16, vjust = 0),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14, hjust = 1),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size=0.5)

ggsave(file = paste('../plots/fitnessGraph_', saveplot, '.pdf', sep = ''),
       device = 'pdf',
       width = 100,
       height = 100,
       units = 'mm')


### GET CONTRIBUTIONS TO SLOPE OF MUT. i

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

muts <- unique(ge_data$knock_in)

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

# plot dF_j as a function of drug dose
contrib$focal_mut <- factor(contrib$focal_mut)
contrib$dose <- factor(contrib$dose, levels = c('0', paste('1e', -2:3, sep = '')))

plot_these_focal_muts <- 'C59R'
plot_these_focal_muts <- muts # FIXME: uncomment for plot with all mutations

ggplot(contrib[contrib$focal_mut %in% plot_these_focal_muts, ],
       aes(x = dose, y = dF_j, color = mut_j, group = mut_j)) +
  geom_point() +
  geom_smooth(aes(fill = mut_j),
              method = 'loess',
              formula = y~x,
              alpha = 0.25) +
  facet_wrap(~focal_mut,
             ncol = 1,
             scales = 'free_y',
             labeller = label_bquote(Mut.~italic(i)~'='~.(levels(contrib$focal_mut)[focal_mut]))) +
  scale_x_discrete(name = expression(paste('Drug dose (', mu, 'M)', sep = '')),
                   breaks = c('0', '1e-1', '1e1', '1e3'),
                   labels = mylabels[c(1,3,5,7)]) +
  scale_y_continuous(name = expression(symbol("\341")*delta*italic(f)[italic(j)]*symbol("\361")[italic(B)(italic(i))]),
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
       width = 150,
       #height = 120,
       height = 200, # FIXME: uncomment if all mutations are being plotted
       units = 'mm')

# plot eps_ij as a function of drug dose
ggplot(contrib[contrib$focal_mut %in% plot_these_focal_muts, ],
       aes(x = dose, y = eps_ij, color = mut_j, group = mut_j)) +
  geom_point() +
  geom_smooth(aes(fill = mut_j),
              method = 'loess',
              formula = y~x,
              alpha = 0.25) +
  facet_wrap(~focal_mut,
             ncol = 1,
             scales = 'free_y',
             labeller = label_bquote(Mut.~italic(i)~'='~.(levels(contrib$focal_mut)[focal_mut]))) +
  scale_x_discrete(name = expression(paste('Drug dose (', mu, 'M)', sep = '')),
                   breaks = c('0', '1e-1', '1e1', '1e3'),
                   labels = mylabels[c(1,3,5,7)]) +
  scale_y_continuous(name = expression(symbol("\341")*italic(epsilon)[italic(ij)]*symbol("\361")),
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
       width = 150,
       #height = 120,
       height = 200, # FIXME: uncomment if all mutations are being plotted
       units = 'mm')

# plot total contributions
contrib_avg <- aggregate(formula = cbind(dF_j, eps_ij) ~ dose + focal_mut + mut_j,
                         data = contrib,
                         FUN = mean)
contrib_avg$term <- contrib_avg$dF_j*contrib_avg$eps_ij

contrib_avg$focal_mut <- factor(contrib_avg$focal_mut)
contrib_avg$dose <- factor(contrib_avg$dose, levels = c('0', paste('1e', -2:3, sep = '')))

contrib_sum <- aggregate(formula = term ~ dose + focal_mut,
                         data = contrib_avg,
                         FUN = sum)
plot_this <- rbind(contrib_avg,
                   data.frame(dose = contrib_sum$dose,
                              focal_mut = contrib_sum$focal_mut,
                              mut_j = 'Sum',
                              dF_j = NA,
                              eps_ij = NA,
                              term = contrib_sum$term))

ggplot(plot_this[plot_this$focal_mut %in% plot_these_focal_muts, ],
       aes(x = dose, y = term, color = mut_j, group = mut_j, linetype = mut_j)) +
  geom_line(size = 1) +
  geom_point() +
  facet_wrap(~focal_mut,
             ncol = 1,
             scales = 'free_y',
             labeller = label_bquote(Mut.~italic(i)~'='~.(levels(contrib$focal_mut)[focal_mut]))) +
  scale_x_discrete(name = expression(paste('Drug dose (', mu, 'M)', sep = '')),
                   breaks = c('0', '1e-1', '1e1', '1e3'),
                   labels = mylabels[c(1,3,5,7)]) +
  scale_y_continuous(name = expression(symbol("\341")*delta*italic(F)[italic(j)]*symbol("\361")[italic(B)(italic(i))]%*%symbol("\341")*italic(epsilon)[italic(ij)]*symbol("\361")),
                     breaks = pretty_breaks(n = 2)) +
  scale_color_manual(name = expression(paste('Mut. ', italic(j), sep = '')),
                     values = c(mut_colors, Sum = 'black')) +
  scale_fill_manual(name = expression(paste('Mut. ', italic(j), sep = '')),
                    values = c(mut_colors, Sum = 'black')) +
  scale_linetype_manual(values = c(rep('solid', length(unique(plot_this$mut_j)) - 1),
                                   'dashed')) +
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

ggsave(file = paste('../plots/terms_', saveplot, '.pdf', sep = ''),
       device = 'pdf',
       width = 150,
       #height = 120,
       height = 200, # FIXME: uncomment if all mutations are being plotted
       units = 'mm')

# plot weighted contributions

weights <- aggregate(formula = dF_j ~ dose + focal_mut + mut_j,
                     data = contrib_avg,
                     FUN = function(x) sum(x^2))
weights <- aggregate(formula = dF_j ~ dose + focal_mut,
                     data = weights,
                     FUN = sum)
colnames(weights)[3] <- 'sum_dF_j_squared'

contrib_sum <- merge(contrib_sum, weights, by = c('dose', 'focal_mut'), all = T)
contrib_sum$weighted_term <- contrib_sum$term / contrib_sum$sum_dF_j_squared

contrib_avg <- merge(contrib_avg, weights, by = c('dose', 'focal_mut'), all = T)
contrib_avg$weighted_term <- contrib_avg$term / contrib_avg$sum_dF_j_squared

plot_this <- rbind(contrib_avg,
                   data.frame(dose = contrib_sum$dose,
                              focal_mut = contrib_sum$focal_mut,
                              mut_j = 'Sum',
                              dF_j = NA,
                              eps_ij = NA,
                              term = contrib_sum$term,
                              sum_dF_j_squared = contrib_sum$sum_dF_j_squared,
                              weighted_term = contrib_sum$weighted_term))

ggplot(plot_this[plot_this$focal_mut %in% plot_these_focal_muts, ],
       aes(x = dose, y = weighted_term, color = mut_j, group = mut_j, linetype = mut_j)) +
  geom_abline(slope = 0, intercept = 0, color = '#d1d3d4') +
  geom_line(size = 1) +
  geom_point() +
  facet_wrap(~focal_mut,
             ncol = 1,
             scales = 'free_y',
             labeller = label_bquote(Mut.~italic(i)~'='~.(levels(contrib$focal_mut)[focal_mut]))) +
  scale_x_discrete(name = expression(paste('Drug dose (', mu, 'M)', sep = '')),
                   breaks = c('0', '1e-1', '1e1', '1e3'),
                   labels = mylabels[c(1,3,5,7)]) +
  scale_y_continuous(name = expression(omega[italic(i)]%*%symbol("\341")*delta*italic(F)[italic(j)]*symbol("\361")[italic(B)(italic(i))]%*%symbol("\341")*italic(epsilon)[italic(ij)]*symbol("\361")),
                     breaks = pretty_breaks(n = 3)) +
  scale_color_manual(name = expression(paste('Mut. ', italic(j), sep = '')),
                     values = c(mut_colors, Sum = 'black')) +
  scale_fill_manual(name = expression(paste('Mut. ', italic(j), sep = '')),
                    values = c(mut_colors, Sum = 'black')) +
  scale_linetype_manual(values = c(rep('solid', length(unique(plot_this$mut_j)) - 1),
                                   'dashed')) +
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

ggsave(file = paste('../plots/weighted_terms_', saveplot, '.pdf', sep = ''),
       device = 'pdf',
       width = 150,
       #height = 120,
       height = 200, # FIXME: uncomment if all mutations are being plotted
       units = 'mm')

# plot weights

ggplot(weights, aes(x = dose, y = 1/sum_dF_j_squared, group = focal_mut, color = focal_mut)) +
  geom_line() +
  geom_point() +
  scale_x_discrete(name = expression(paste('Drug dose (', mu, 'M)', sep = '')),
                   breaks = c('0', '1e-1', '1e1', '1e3'),
                   labels = mylabels[c(1,3,5,7)]) +
  scale_y_continuous(name = expression(paste('Weight, ', italic(omega), sep = '')),
                     trans = 'log10') +
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
        legend.text = element_text(size = 14))

ggsave(file = paste('../plots/weights_', saveplot, '.pdf', sep = ''),
       device = 'pdf',
       width = 150,
       #height = 120,
       height = 200, # FIXME: uncomment if all mutations are being plotted
       units = 'mm')
