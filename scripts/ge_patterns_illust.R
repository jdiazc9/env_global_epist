rm(list = ls())
set.seed(0)
source('ecoFunctions.R')

L <- 8 # number of mutations

# mutations and genotype names
muts <- paste('m', 1:L, sep = '')
genots <- orderName(unlist(sapply(1:length(muts),
                                  FUN = function(i) sapply(combn(muts,
                                                                 i,
                                                                 simplify = FALSE),
                                                           FUN = paste, collapse = ','))))
genots <- c('', genots)

lambda_i <- setNames(rnorm(L, mean = 0.3, sd = 0.3), muts)


### CASE 1: fitness is exponential transformation on latent additive variable
f <- function(lambda) 1 - exp(-lambda)

sumLambda <- data.frame(genotype = genots,
                        fitness = sapply(genots,
                                        FUN = function(g) {
                                          
                                          g <- strsplit(g, split = ',')[[1]]
                                          lambda <- sum(lambda_i[g])
                                          return(lambda)
                                          
                                        }))

landscape <- data.frame(genotype = genots,
                        fitness = sapply(genots,
                                        FUN = function(g) {
                                          
                                          g <- strsplit(g, split = ',')[[1]]
                                          lambda <- sum(lambda_i[g])
                                          return(f(lambda))
                                          
                                        }))
gedf <- makeGEdata(landscape)

ggplot(gedf[gedf$knock_in == 'm1', ], aes(x = background_f, y = d_f)) +
  geom_abline(slope = 0, intercept = 0, color = '#d1d3d4') +
  geom_point(color = 'black',
             cex = 3) +
  geom_smooth(method = 'lm',
              formula = y~x,
              fullrange = T,
              se = F,
              color = '#1e9fd3') +
  scale_x_continuous(name = expression(paste(italic(f), '(', italic(B), ')', sep = '')),
                     breaks = pretty_breaks(n = 3),
                     expand = rep(0.05, 2)) +
  scale_y_continuous(name = expression(paste(Delta, italic(f), sep = '')),
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
        legend.text = element_text(size = 14, hjust = 1))

ggsave(file = '../plots/illustrations/case1.pdf',
       device = 'pdf',
       width = 100,
       height = 70,
       units = 'mm')

ggplot() +
  geom_function(fun = f, color = 'orange') +
  scale_x_continuous(name = expression(italic(lambda)),
                     breaks = pretty_breaks(n = 3),
                     limits = c(0, 5),
                     expand = rep(0.05, 2)) +
  scale_y_continuous(name = expression(paste('Fitness, ', italic(f), '(', italic(lambda), ')', sep = '')),
                     breaks = pretty_breaks(n = 3),
                     expand = rep(0.05, 2),
                     limits = c(0, 1)) +
  theme_bw() +
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 16, vjust = 0),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14, hjust = 1))

ggsave(file = '../plots/illustrations/case1_lambda.pdf',
       device = 'pdf',
       width = 100,
       height = 70,
       units = 'mm')



# CASE 2: fitness is non-exponential transformation on latent additive variable
f <- function(lambda) 1 - ((lambda-5)/5)^2
f <- function(lambda) lambda^2 / (1 + lambda^2)

landscape <- data.frame(genotype = genots,
                        fitness = sapply(genots,
                                        FUN = function(g) {
                                          
                                          g <- strsplit(g, split = ',')[[1]]
                                          lambda <- sum(lambda_i[g])
                                          return(f(lambda))
                                          
                                        }))
gedf <- makeGEdata(landscape)

analytic_f <- data.frame(lB = seq(min(sumLambda$fitness[!grepl('m1', landscape$genotype)]),
                                  max(sumLambda$fitness[!grepl('m1', landscape$genotype)]),
                                  length.out = 200))
analytic_f$lBi <- analytic_f$lB + lambda_i['m1']
analytic_f$fB <- f(analytic_f$lB)
analytic_f$fBi <- f(analytic_f$lBi)
analytic_f$df <- analytic_f$fBi - analytic_f$fB

ggplot(gedf[gedf$knock_in == 'm1', ], aes(x = background_f, y = d_f)) +
  geom_abline(slope = 0, intercept = 0, color = '#d1d3d4') +
  geom_point(color = 'black',
             cex = 3) +
  geom_path(data = analytic_f,
            aes(x = fB, y = df),
            color = 'gray',
            linewidth = 0.5) +
  geom_smooth(method = 'lm',
              formula = y~x,
              fullrange = T,
              se = F,
              color = '#1e9fd3') +
  scale_x_continuous(name = expression(paste(italic(f), '(', italic(B), ')', sep = '')),
                     breaks = pretty_breaks(n = 3),
                     expand = rep(0.05, 2)) +
  scale_y_continuous(name = expression(paste(Delta, italic(f), sep = '')),
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
        legend.text = element_text(size = 14, hjust = 1))

ggsave(file = '../plots/illustrations/case2.pdf',
       device = 'pdf',
       width = 100,
       height = 70,
       units = 'mm')

ggplot() +
  geom_function(fun = f, color = 'orange') +
  scale_x_continuous(name = expression(italic(lambda)),
                     breaks = pretty_breaks(n = 3),
                     limits = c(0, 5),
                     expand = rep(0.05, 2)) +
  scale_y_continuous(name = expression(paste('Fitness, ', italic(f), '(', italic(lambda), ')', sep = '')),
                     breaks = pretty_breaks(n = 3),
                     expand = rep(0.05, 2),
                     limits = c(0, 1)) +
  theme_bw() +
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 16, vjust = 0),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14, hjust = 1))

ggsave(file = '../plots/illustrations/case2_lambda.pdf',
       device = 'pdf',
       width = 100,
       height = 70,
       units = 'mm')


# CASE 3: same as case 2 but with noise
landscape$fitness <- landscape$fitness + rnorm(nrow(landscape), mean = 0, sd = 0.05)
gedf <- makeGEdata(landscape)

analytic_f <- data.frame(lB = seq(min(sumLambda$fitness[!grepl('m1', landscape$genotype)]),
                                  max(sumLambda$fitness[!grepl('m1', landscape$genotype)]),
                                  length.out = 200))
analytic_f$lBi <- analytic_f$lB + lambda_i['m1']
analytic_f$fB <- f(analytic_f$lB)
analytic_f$fBi <- f(analytic_f$lBi)
analytic_f$df <- analytic_f$fBi - analytic_f$fB

ggplot(gedf[gedf$knock_in == 'm1', ], aes(x = background_f, y = d_f)) +
  geom_abline(slope = 0, intercept = 0, color = '#d1d3d4') +
  geom_point(color = 'black',
             cex = 3) +
  geom_path(data = analytic_f,
            aes(x = fB, y = df),
            color = 'gray',
            linewidth = 0.5) +
  geom_smooth(method = 'lm',
              formula = y~x,
              fullrange = T,
              se = F,
              color = '#1e9fd3') +
  scale_x_continuous(name = expression(paste(italic(f), '(', italic(B), ')', sep = '')),
                     breaks = pretty_breaks(n = 3),
                     expand = rep(0.05, 2)) +
  scale_y_continuous(name = expression(paste(Delta, italic(f), sep = '')),
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
        legend.text = element_text(size = 14, hjust = 1))

ggsave(file = '../plots/illustrations/case3.pdf',
       device = 'pdf',
       width = 100,
       height = 70,
       units = 'mm')








