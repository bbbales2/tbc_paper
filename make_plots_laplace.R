library(tidyverse)
library(GGally)
library(Hmisc)
library(loo)

read_laplace_csv = function(file) {
  lines = readLines(file)
  starts_with_pound = lapply(lines, function(line) {
    return(substr(line, 1, 1) == "#")
  }) %>%
    unlist
  
  skip = which(starts_with_pound) %>% tail(1)
  n_max = which(starts_with_pound) %>% tail(1) - which(starts_with_pound) %>% tail(2) %>% head(1) - 2
  
  list(df_opt = read_csv(file, comment = "#", n_max = n_max, col_types = cols(.default = col_double())),
       df_laplace = read_csv(file, comment = "#", skip = skip, col_types = cols(.default = col_double())) %>%
         filter(log_p != 0.0))
}

problems = c("bl1", "bl1b", "bl1c", "bl2a")
basenames = c("output.1.csv", "output.2.csv", "output.3.csv", "output.4.csv")

opt_df = lapply(problems, function(problem) {
  filenames = paste0(problem, "_with_priors/", basenames)
  
  lapply(1:length(filenames), function(i) {
    file = filenames[i]
    
    tryCatch({
      df_opt = read_laplace_csv(file)$df_laplace
      
      steps = nrow(df_opt)
      
      return(df_opt %>%
               mutate(chain = i,
                      iteration = row_number() - steps,
                      problem = problem))
    }, error = function(e) {
      print(e)
      print(file)
      return(NULL)
    })
  }) %>% bind_rows
}) %>% bind_rows

opt_df %>%
  filter(lp__ > -10000) %>%
  select(iteration, chain, problem, lp__, c11_coat, c12_coat, c13_coat, c33_coat, c44_coat, c66_coat) %>%
  gather(which, value, -iteration, -chain, -problem) %>%
  ggplot(aes()) +
  geom_point(aes(iteration, value, color = factor(chain)), size = 0.25) +
  facet_grid(which ~ problem, scale = "free")

#df1 = opt_df %>% mutate(log_ratios = log_p - log_g)
#E_loo(df1 %>% pull(c11_coat), loo::psis(df1 %>% pull(log_ratios)),
#      type = "quantile", probs = 0.1, log_ratios = df1 %>% pull(log_ratios))$value

opt_df %>%
  filter(log_p > -500) %>%
  filter(problem == "bl1b") %>%
  mutate(log_ratios = log_p - log_g) %>%
  group_by(chain, problem) %>%
  select(chain, problem, log_ratios, c11_coat, c12_coat, c13_coat, c33_coat, c44_coat, c66_coat) %>%
  gather(which, value, -chain, -problem, -log_ratios) %>%
  group_by(which, chain, problem) %>%
  summarise(m = E_loo(value, loo::psis(log_ratios),
                      type = "mean", log_ratios = log_ratios)$value,
            ql = E_loo(value, loo::psis(log_ratios),
                       type = "quantile", probs = 0.1, log_ratios = log_ratios)$value,
            qh = E_loo(value, loo::psis(log_ratios),
                       type = "quantile", probs = 0.9, log_ratios = log_ratios)$value) %>%
  ggplot(aes()) +
  geom_errorbar(aes(problem, ymin = ql, ymax = qh, group = chain, color = problem), width = 0,
                position = position_dodge(width = 0.5)) +
  geom_point(aes(problem, m, group = chain, color = problem), position=position_dodge(width=0.5)) +
  facet_grid(~ which) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("80% intervals + means + multiple chains\n only keeping samples with lp > -500")

opt_df %>%
  filter(log_p > -500) %>%
  select(chain, problem, log_p, log_g) %>%
  mutate(log_ratio = log_p - log_g) %>%
  group_by(chain, problem) %>%
  summarise(pareto_k = loo::psis(log_ratio)$diagnostics$pareto_k,
            n_eff = loo::psis(log_ratio)$diagnostics$n_eff,
            draws = n()) %>%
  filter(problem == "bl1b")

