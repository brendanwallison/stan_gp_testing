functions {
  matrix gp_pred_rng(array[] real x2,
                     matrix y1,
                     array[] real x1,
                     real rho,
                     real sigma,
                     real delta,
                     vector alpha,
                     matrix L_Omega) {
    int N1 = rows(y1);
    int D = cols(y1);
    int N2 = size(x2);
    matrix[N2, D] f2;
    {
      matrix[N1, N1] L_K;
      matrix[N2, N2] L_K_2;
      matrix[N1, D] K_div_y1;
      matrix[N1, N2] k_x1_x2;
      matrix[N1, N2] v_pred;
      matrix[N2, D] f2_mu;
      matrix[N2, N2] cov_f2;
      matrix[N2, N2] diag_delta;
      matrix[N1, N1] K;
      matrix[N2, D] eta;
      K = gp_exp_quad_cov(x1, 1.0, rho);
      real sq_sigma = square(sigma);
      for (n1 in 1:N1) {
        K[n1, n1] = K[n1, n1] + sq_sigma;
      }
      L_K = cholesky_decompose(K);
      K_div_y1 = mdivide_left_tri_low(L_K, y1);
      K_div_y1 = mdivide_right_tri_low(K_div_y1', L_K)';
      k_x1_x2 = gp_exp_quad_cov(x1, x2, 1.0, rho);
      f2_mu = (k_x1_x2' * K_div_y1);
      v_pred = mdivide_left_tri_low(L_K, k_x1_x2);
      cov_f2 = gp_exp_quad_cov(x2, 1.0, rho) - v_pred' * v_pred;
      diag_delta = diag_matrix(rep_vector(delta, N2));
      L_K_2 = cholesky_decompose(cov_f2 + diag_delta);
      // for (y in eta){
      //   y = std_normal_rng();
      // }
      for (c in 1:D){
        for (r in 1:N2){
          eta[r, c] = std_normal_rng();
        }
      }
      f2 = L_K_2 * eta * diag_pre_multiply(alpha, L_Omega)'  + f2_mu;
    return f2;
    }
  }
}
data {
  int<lower=1> N1;
  int<lower=1> D;
  array[N1] real x1;
  matrix[N1, D] y1;
  int<lower=1> N2;
  array[N2] real x2;
}
transformed data {
  vector[N1] mu = rep_vector(0, N1);
  real delta = 1e-9;
}
parameters {
  real<lower=0> rho;
  vector<lower=0>[D] alpha;
  real<lower=0> sigma;
  cholesky_factor_corr[D] L_Omega;
  matrix[N1, D] eta;
}
model {
  matrix[N1, D] f;
  {
    matrix[N1, N1] K = gp_exp_quad_cov(x1, 1.0, rho);
    matrix[N1, N1] L_K;
    real sq_sigma = square(sigma);

    // diagonal elements
    for (n1 in 1:N1) {
      K[n1, n1] = K[n1, n1] + sq_sigma;
    }

    L_K = cholesky_decompose(K);
    f = L_K * eta
        * diag_pre_multiply(alpha, L_Omega)';
  }

  L_Omega ~ lkj_corr_cholesky(3);
  rho ~ inv_gamma(5, 5);
  alpha ~ std_normal();
  sigma ~ std_normal();
  to_vector(eta) ~ std_normal();
  to_vector(y1) ~ normal(to_vector(f), sigma);
}
generated quantities {
  matrix[N1, D] f;
  {
    matrix[N1, N1] K = gp_exp_quad_cov(x1, 1.0, rho);
    matrix[N1, N1] L_K;
    real sq_sigma = square(sigma);

    // diagonal elements
    for (n1 in 1:N1) {
      K[n1, n1] = K[n1, n1] + sq_sigma;
    }

    L_K = cholesky_decompose(K);
    f = L_K * eta
        * diag_pre_multiply(alpha, L_Omega)';
  }

  
  matrix[N2, D] f2;
  matrix[N2, D] f2_with_mo_corr;
  matrix[N2, D] y2;
  f2 = gp_pred_rng(x2, y1, x1, rho, sigma, delta, alpha, L_Omega);
  // for (n2 in 1:N2) {
  //   y2[n2] = normal_rng(f2_with_mo_corr[n2], sigma);
  // }
}
