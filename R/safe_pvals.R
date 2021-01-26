## utility functions 
log_sum_exp <- function(logx, logy) {
  if (logx > logy) {
    a = logx;
  } else {
    a = logy;
  }
  if (abs(a) == Inf) {
    a = 0;
  }
  out = exp(logx - a);
  out = out + exp(logy - a);
  return(log(out) + a)
}

log1mexp <- function(a) {
  if (a >= 0 && a <= log(2)) {
    return(log(-expm1(-a)));
  } else if (a > log(2)) {
    return(log1p(-exp(-a)));
  } else {
    stop(paste0("trying to log1mexp with a = ", a))
  }
}


log_subtract <- function(x, y) {
  if (x < y) {
    stop("log_subtract:: cannot take log of (-)ve number");
  }
  return(x + log1mexp(abs(y - x)));
}

calc_p_value_safer <- function(df, vTy, nu_norm, sig, mu = 0) {
  # TAIL_THRESHOLD <- 30
  n_intervals <- dim(df)[[1]]
  # if (abs(vTy) / sqrt(nu_norm * sig) > TAIL_THRESHOLD) {
  #   cur_interval <- df[1, ]
  #   i <- 1
  #   while (cur_interval$contained == 1 && i <= n_intervals) {
  #     b = cur_interval$max_mean;
  #     i <- i + 1;
  #     cur_interval <- df[i, ]
  #   }
  #   
  #   cur_interval <- df[n_intervals, ]
  #   i <- n_intervals
  #   while (cur_interval$contained == 1 && i > 0) {
  #     a = cur_interval$min_mean;
  #     i <- i - 1
  #     cur_interval <- df[i, ]
  #   }
  #   
  #   s_bound = max(a, abs(b));
  #   k1 = abs(vTy) / sqrt(nu_norm * sig);
  #   k2 = s_bound / sqrt(nu_norm * sig);
  #   return(min(1.0, exp(0.5 * (k2 * k2 - k1 * k1)) * (k2 * k2 + 1) / (k1 * k2)))
  # }
  # 
  n1 = -Inf;
  d1 = -Inf;
  for (i in 1:n_intervals) {
    cur_interval <- df[i, ]
    if (cur_interval$contained == 1) {
      a = pnorm((cur_interval$max_mean - mu) / sqrt(nu_norm * sig), log.p = TRUE);
      b = pnorm((cur_interval$min_mean - mu) / sqrt(nu_norm * sig), log.p = TRUE);
      arg2 = log_subtract(a, b);
      d1 = log_sum_exp(d1, arg2);
      if (cur_interval$max_mean >= vTy) {
        arg2 = log_subtract(pnorm((cur_interval$max_mean - mu) / sqrt(nu_norm * sig), log.p = TRUE ),
                            pnorm((max(cur_interval$min_mean, vTy) - mu) / sqrt(nu_norm * sig), log.p = TRUE));
        n1 = log_sum_exp(n1, arg2);
      }
    }
  }
  
  return (exp(n1 - d1));
}


