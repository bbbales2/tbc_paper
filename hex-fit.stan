functions {
  int bilayer_lookup_size(int IN, int JN, int KN);
  vector bilayer_init(int IN, int JN, int layer_index, real X, real Y, vector Zs, real density_sub, real density_coat);
  vector bilayer_rus(int N, vector lookup, matrix C1, matrix C2);
  matrix build_C_cubic(real c11, real c12, real c44) {
    matrix[6, 6] C;
    
    for (i in 1:6)
      for (j in 1:6)
        C[i, j] = 0.0;
      
      C[1, 1] = c11;
      C[2, 2] = c11;
      C[3, 3] = c11;
      C[4, 4] = c44;
      C[5, 5] = c44;
      C[6, 6] = c44;
      C[1, 2] = c12;
      C[1, 3] = c12;
      C[2, 3] = c12;
      C[3, 2] = c12;
      C[2, 1] = c12;
      C[3, 1] = c12;
      
      return C;
  }
  matrix build_C_hex(real c11, real c12, real c13, real c33, real c44) {
    real c66;
    matrix[6, 6] C;
    
    c66 = (c11-c12)/2;
    
    for (i in 1:6)
      for (j in 1:6)
        C[i, j] = 0.0;
      
      C[1, 1] = c11;
      C[2, 2] = c11;
      C[3, 3] = c33;
      C[4, 4] = c44;
      C[5, 5] = c44;
      C[6, 6] = c66;
      C[1, 2] = c12;
      C[1, 3] = c13;
      C[2, 3] = c13;
      C[3, 2] = c13;
      C[2, 1] = c12;
      C[3, 1] = c13;
      
      return C;
  }
  matrix build_C_ortho(real c11, real c22, real c33, real c44, real c55, real c66, real c12, real c13, real c23) {
    matrix[6, 6] C;
    for (i in 1:6)
      for (j in 1:6)
        C[i, j] = 0.0;
      
      C[1, 1] = c11;
      C[2, 2] = c22;
      C[3, 3] = c33;
      C[4, 4] = c44;
      C[5, 5] = c55;
      C[6, 6] = c66;
      C[1, 2] = c12;
      C[1, 3] = c13;
      C[2, 3] = c23;
      C[3, 2] = c23;
      C[2, 1] = c12;
      C[3, 1] = c13;
      
      return C;
  }
  //vector mech_init(int P, real X, real Y, real Z, real density);
  matrix mech_rotate(matrix C, vector q);
  //vector mech_rus(int N, vector lookup, matrix C);
  vector cu2qu(vector cu);
}

data {
  int<lower = 1> PX;
  int<lower = 1> PY;
  int<lower = 0> PZ_sub;
  int<lower = 0> PZ_coat;
  int<lower = 1> N; // Number of resonance modes
  
  real<lower = 0.0> X;
  real<lower = 0.0> Y;
  real<lower = 0.0> Z_sub;
  real<lower = 0.0> Z_coat;
  
  real<lower = 0.0> density_sub;
  real<lower = 0.0> density_coat;
  
  real<lower = 0.0> c11_sub;
  real<lower = 0.0> c12_sub;
  real<lower = 0.0> c44_sub;
  
  vector<lower = -1.0725146985555127, upper = 1.0725146985555127>[3] cu;
  
  real measfq[N];

  real<lower = 0.0> prior_sd_c11_coat;
  real<lower = 0.0> prior_sd_c12_coat;
  real<lower = 0.0> prior_sd_c13_coat;
  real<lower = 0.0> prior_sd_c33_coat;
  real<lower = 0.0> prior_sd_c44_coat;
  real<lower = 0.0> prior_sd_sigma;
}

transformed data {
  vector[bilayer_lookup_size(PX, PY, PZ_sub + PZ_coat + 3)] lookup;
  int layer_index;
  {
    vector[PZ_sub + PZ_coat + 3] zs;
    zs[1] = 0.0;
    
    for(i in 1:PZ_sub) {
      zs[i + 1] = Z_sub * i / (PZ_sub + 1.0);
    }
    zs[1 + PZ_sub + 1] = Z_sub;
    
    for(i in 1:PZ_coat) {
      zs[i + 1 + PZ_sub + 1] = Z_sub + Z_coat * i / (PZ_coat + 1.0);
    }
    zs[1 + PZ_sub + 1 + PZ_coat + 1] = Z_coat + Z_sub;
    
    lookup = bilayer_init(PX, PY, 1 + PZ_sub + 1, X, Y, zs, density_sub, density_coat);
  }
}

// See: http://solidmechanics.org/text/Chapter3_2/Chapter3_2.htm   for definitions
parameters {
  real<lower = 0.0> c11_coat;
  real c12_coat;
  real c13_coat;
  real<lower = 0.0> c33_coat;
  real<lower = 0.0> c44_coat;
  real<lower = 0.0> sigma_z;
}

transformed parameters {
  real sigma = sigma_z * prior_sd_sigma;
  real<lower = 0.0> c66_coat = (c11_coat - c12_coat) / 2.0;
  real ct1_coat = c11_coat - (c13_coat*c13_coat/c33_coat);
  
  matrix[6, 6] Csub = mech_rotate(build_C_cubic(c11_sub, c12_sub, c44_sub), cu2qu(cu));
}

model {
  real c11_sub_ = c11_sub;
  vector[N] calcfq =
    bilayer_rus(N, lookup,
		mech_rotate(build_C_cubic(c11_sub_, c12_sub, c44_sub), cu2qu(cu)),
		build_C_hex(c11_coat, c12_coat, c13_coat, c33_coat, c44_coat));
  
  c11_coat ~ normal(0.0, prior_sd_c11_coat);
  c12_coat ~ normal(0.0, prior_sd_c12_coat);
  c13_coat ~ normal(0.0, prior_sd_c13_coat);
  c33_coat ~ normal(0.0, prior_sd_c33_coat);
  c44_coat ~ normal(0.0, prior_sd_c44_coat);

  sigma_z ~ normal(0.0, 1.0);

  /**
   * Measured resonance modes (measfq) are normally distributed around what
   * you'd expect from an RUS calculation (calculatedfreqs) with some noise
   * (sigma).
   *
   * for loop says: if a measured frequency is reported to be negative,
   * then exclude this mode from the fit.
   */
  for(i in 1:N) {
    if(measfq[i] > 0) {
      measfq[i] ~ lognormal(log(calcfq[i]), sigma);
    }
  }
}

generated quantities {
  vector[N] yhat;
  vector[N] residual;
  vector[N] calcfq = bilayer_rus(N, lookup,
          mech_rotate(build_C_cubic(c11_sub, c12_sub, c44_sub), cu2qu(cu)),
          build_C_hex(c11_coat, c12_coat, c13_coat, c33_coat, c44_coat));

  for(n in 1:N) {
    residual[n] = measfq[n] - calcfq[n];
    yhat[n] = lognormal_rng(log(calcfq[n]), sigma);
  }
}
