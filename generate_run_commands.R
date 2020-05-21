library(rstan)

data_files = c("bl1.dat", "bl1b.dat", "bl1c.dat", "bl2a.dat")

priors = list(prior_sd_c11_coat = 1.0,
              prior_sd_c12_coat = 1.0,
              prior_sd_c13_coat = 1.0,
              prior_sd_c33_coat = 1.0,
              prior_sd_c44_coat = 1.0,
              prior_sd_sigma = 0.01)

reps = 4

data_files_with_priors = NULL

for(file in data_files) {
  data_with_priors = c(rstan::read_rdump(file), priors)
  data_file_with_prior = paste0(tools::file_path_sans_ext(file), "_with_priors.dat")
  stan_rdump(names(data_with_priors), data_file_with_prior, env = list2env(data_with_priors))
  data_files_with_priors = c(data_files_with_priors, data_file_with_prior)
}

for(file in data_files_with_priors) {
  fit_cmds = NULL

  data = tools::file_path_sans_ext(file)
  
  if(!dir.exists(data)) {
    dir.create(data)
  }

  for(r in 1:reps) {
    inits = list(c11_coat = runif(0, priors$prior_sd_c11_coat),
                 c11_coat = runif(0, priors$prior_sd_c12_coat),
                 c11_coat = runif(0, priors$prior_sd_c13_coat),
                 c11_coat = runif(0, priors$prior_sd_c33_coat),
                 c11_coat = runif(0, priors$prior_sd_c44_coat),
                 sigma_z = runif(0, 1))
    
    init_file = paste0(data, "/init.", r, ".dat")
    output_file = paste0(data, "/output.", r, ".csv")
    log_file = paste0(data, "/log.", r, ".txt")
    
    stan_rdump(names(inits), init_file, env = list2env(inits))
    
    fit_cmds = c(fit_cmds, paste0("./hex-fit optimize save_iterations=1 laplace_draws=200 init=", init_file, " data file=file output file=", output_file, " >& ", log_file))
  }
  writeLines(fit_cmds, paste0("run_", data, "_fits.sh"))
}
