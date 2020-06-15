library(tidyverse)
library(GGally)

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
constraint_names = c("both_pos", "c13_pos", "c12_pos", "neither_pos")
basenames = c("output.1.csv", "output.2.csv", "output.3.csv", "output.4.csv")

error_files = NULL
opt_df = lapply(problems, function(problem) {
  lapply(constraint_names, function(constraint) {
    filenames = paste0(problem, "/", constraint, "/", basenames)
  
    lapply(1:length(filenames), function(i) {
      file = filenames[i]
      
      tryCatch({
        df_opt = read_laplace_csv(file)$df_opt
        
        steps = nrow(df_opt)
        
        return(df_opt %>%
          mutate(chain = i,
                 iteration = row_number() - steps,
                 constraint = constraint,
                 problem = problem))
      }, error = function(e) {
        print(e)
        print(file)
        error_files <<- c(error_files, file)
        return(NULL)
      })
    }) %>% bind_rows
  }) %>% bind_rows
}) %>% bind_rows

opt_df %>%
  filter(lp__ > -10000) %>%
  select(iteration, constraint, chain, problem, lp__, c11_coat, c12_coat, c13_coat, c33_coat, c44_coat, c66_coat) %>%
  gather(which, value, -iteration, -chain, -problem, -constraint) %>%
  ggplot(aes()) +
  geom_point(aes(iteration, value, color = constraint), size = 0.05) +
  facet_grid(which ~ problem, scale = "free")

opt_df %>%
  filter(lp__ > -500) %>%
  filter(iteration == 0) %>%
  select(iteration, chain, problem, c11_coat, c12_coat, c13_coat, c33_coat, c44_coat, c66_coat) %>%
  gather(which, value, -iteration, -chain, -problem) %>%
  ggplot(aes()) +
  geom_point(aes(problem, value)) +
  facet_grid(~ which) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
