// marginalized.stan  (no Z in sampler; compatible with older stanc)
data {
  int<lower=1> N;
  int<lower=1> I;
  int<lower=1, upper=3> Y[N, I];  // 1=No, 2=Perhaps, 3=Yes
}
parameters {
  // Gate
  real alpha0;
  real<lower=0> sigma_eps;
  matrix[N, I] eps_raw;               // std normal -> eps = sigma_eps * eps_raw

  // Traits
  vector[N] theta_dir;
  vector[N] kappa_ext;

  // Items (we'll center bA1,bB1 manually)
  vector[I] bA1_raw;
  vector[I] bB1_raw;
  vector[I] bA2;
  vector[I] bB2;
}
transformed parameters {
  vector[I] bA1;
  vector[I] bB1;

  matrix[N, I] eps;
  matrix[N, I] logit_p;
  matrix[N, I] p_gate;

  // center bA1_raw, bB1_raw without broadcasting
  real mA1;
  real mB1;
  mA1 = mean(bA1_raw);
  mB1 = mean(bB1_raw);
  for (i in 1:I) {
    bA1[i] = bA1_raw[i] - mA1;
    bB1[i] = bB1_raw[i] - mB1;
  }
  // element-wise transforms (avoid matrix broadcasting on older Stan)
  for (n in 1:N) {
    for (i in 1:I) {
      eps[n,i]     = sigma_eps * eps_raw[n,i];
      logit_p[n,i] = alpha0 + eps[n,i];
      p_gate[n,i]  = inv_logit(logit_p[n,i]);
    }
  }
}
model {
  // Priors (avoid to_vector on matrix for older Stan)
  for (n in 1:N) for (i in 1:I) eps_raw[n,i] ~ normal(0, 1);
  alpha0    ~ normal(0, 1);
  sigma_eps ~ normal(0, 1);      // half-normal via <lower=0>
  theta_dir ~ normal(0, 1);
  kappa_ext ~ normal(0, 1);
  bA1_raw   ~ normal(0, 1.5);
  bB1_raw   ~ normal(0, 1.5);
  bA2       ~ normal(0, 1.5);
  bB2       ~ normal(0, 1.5);

  // Likelihood (Z marginalized)
  for (n in 1:N) {
    for (i in 1:I) {
      real pA1; real pA2_;
      real pB1; real pB2_;
      vector[3] pA;
      vector[3] pB;
      vector[3] pmix;

      // Tree A
      pA1 = inv_logit(theta_dir[n] - bA1[i]);  // not-No
      pA2_ = inv_logit(theta_dir[n] - bA2[i]); // Yes | not-No
      pA[3] = pA1 * pA2_;                      // Yes
      pA[1] = 1 - pA1;                         // No
      pA[2] = 1 - pA[1] - pA[3];               // Perhaps

      // Tree B
      pB1 = inv_logit(kappa_ext[n] - bB1[i]);  // Extreme
      pB2_ = inv_logit(theta_dir[n] - bB2[i]); // Yes | Extreme
      pB[2] = 1 - pB1;                         // Perhaps
      pB[3] = pB1 * pB2_;                      // Yes
      pB[1] = 1 - pB[2] - pB[3];               // No

      pmix = p_gate[n,i] * pA + (1 - p_gate[n,i]) * pB;
      target += categorical_lpmf(Y[n,i] | pmix);
    }
  }
}
generated quantities {
  matrix[N, I] w_tree1;
  matrix[N, I] log_lik;

  for (n in 1:N) {
    for (i in 1:I) {
      real pA1; real pA2_;
      real pB1; real pB2_;
      vector[3] pA;
      vector[3] pB;
      vector[3] pmix;
      int y;

      pA1 = inv_logit(theta_dir[n] - bA1[i]);
      pA2_ = inv_logit(theta_dir[n] - bA2[i]);
      pA[3] = pA1 * pA2_;
      pA[1] = 1 - pA1;
      pA[2] = 1 - pA[1] - pA[3];

      pB1 = inv_logit(kappa_ext[n] - bB1[i]);
      pB2_ = inv_logit(theta_dir[n] - bB2[i]);
      pB[2] = 1 - pB1;
      pB[3] = pB1 * pB2_;
      pB[1] = 1 - pB[2] - pB[3];

      pmix = p_gate[n,i] * pA + (1 - p_gate[n,i]) * pB;

      y = Y[n,i];
      w_tree1[n,i] = inv_logit(logit_p[n,i] + log(pA[y]) - log(pB[y]));
      log_lik[n,i] = categorical_lpmf(y | pmix);
    }
  }
}
