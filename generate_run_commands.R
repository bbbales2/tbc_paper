library(rstan)

data_files = c("bl1.dat", "bl1b.dat", "bl1c.dat", "bl2a.dat")

priors = list(prior_sd_c11_coat = 1.0,
              prior_sd_c12_coat = 1.0,
              prior_sd_c13_coat = 1.0,
              prior_sd_c33_coat = 1.0,
              prior_sd_c44_coat = 1.0,
              prior_sd_sigma = 0.01)

reps = 4

constraint_names = c("both_pos", "c13_pos", "c12_pos", "neither_pos")
constraints = list(list(c12_positive = 1, c13_positive = 1),
                   list(c12_positive = 0, c13_positive = 1),
                   list(c12_positive = 1, c13_positive = 0),
                   list(c12_positive = 0, c13_positive = 0))

fit_cmds = NULL

for(file in data_files) {
  for(c in 1:length(constraints)) {
    data = paste0(tools::file_path_sans_ext(file), "/", constraint_names[c])
    base_stan_data = read_rdump(file)
    
    stan_data = c(base_stan_data, priors, constraints[[c]])
    
    data_file = paste0(data, "/data.dat")
    
    if(!dir.exists(data)) {
      dir.create(data, recursive = TRUE)
    }
    
    stan_rdump(names(stan_data), data_file, env = list2env(stan_data))
    
    for(r in 1:reps) {
      c66_coat = -1
      C = matrix(0, nrow = 6, ncol = 6)
      while(c66_coat < 0 || any(eigen(C)$values <= 0)) {
          inits = list(c11_coat = runif(1, 0, priors$prior_sd_c11_coat),
                       c12_coat = runif(1, 0, priors$prior_sd_c12_coat),
                       c13_coat = runif(1, 0, priors$prior_sd_c13_coat),
                       c33_coat = runif(1, 0, priors$prior_sd_c33_coat),
                       c44_coat = runif(1, 0, priors$prior_sd_c44_coat),
                       sigma_z = runif(1, 0, 1))
          c66_coat = (inits$c11_coat - inits$c12_coat) / 2.0
  
          C = matrix(0, nrow = 6, ncol = 6)
  
          C[1, 1] = inits$c11_coat
          C[2, 2] = inits$c11_coat
          C[3, 3] = inits$c33_coat
          C[4, 4] = inits$c44_coat
          C[5, 5] = inits$c44_coat
          C[6, 6] = c66_coat
          C[1, 2] = inits$c12_coat
          C[1, 3] = inits$c13_coat
          C[2, 3] = inits$c13_coat
          C[3, 2] = inits$c13_coat
          C[2, 1] = inits$c12_coat
          C[3, 1] = inits$c13_coat
      }
      
      init_file = paste0(data, "/init.", r, ".dat")
      output_file = paste0(data, "/output.", r, ".csv")
      log_file = paste0(data, "/log.", r, ".txt")
      
      stan_rdump(names(inits), init_file, env = list2env(inits))
      
      fit_cmds = c(fit_cmds, paste0("./hex-fit optimize save_iterations=1 laplace_draws=200 init=", init_file, " data file=", data_file, " output file=", output_file, " >& ", log_file))
    }
  }
}
writeLines(fit_cmds, paste0("run_fits.sh"))
