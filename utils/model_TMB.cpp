#include <TMB.hpp>
#include <limits> // for numeric_limits

// Custom lchoose function using lgamma
template<class Type>
Type lchoose(Type n, Type k) {
  return lgamma(n+Type(1.0)) - lgamma(k+Type(1.0)) - lgamma(n-k+Type(1.0));
}

template<class Type>
Type objective_function<Type>::operator() () {
  // Data
  DATA_MATRIX(xj);    // observed counts: P x T
  DATA_VECTOR(times); // timepoints (length T)
  DATA_SCALAR(dt);    // time scaling factor
  
  int P = xj.rows();
  int T = xj.cols();
  
  // Parameters
  PARAMETER_VECTOR(log_x0i); // length P
  PARAMETER_VECTOR(fi);      // length P
  
  // Compute total counts per timepoint
  vector<Type> N(T);
  for (int t = 0; t < T; t++) {
    Type sum_col = 0.0;
    for (int i = 0; i < P; i++) {
      sum_col += xj(i, t);
    }
    N(t) = sum_col;
  }
  
  // Initialize max_col with -∞
  vector<Type> max_col(T);
  for (int t = 0; t < T; t++) {
    max_col(t) = -std::numeric_limits<Type>::infinity();
  }
  
  // Compute log_abundances and track max per column
  matrix<Type> log_abundance(P, T);
  for (int t = 0; t < T; t++) {
    for (int i = 0; i < P; i++) {
      Type val = log_x0i(i) + fi(i)*times(t)*dt;
      log_abundance(i, t) = val;
      if (val > max_col(t)) max_col(t) = val;
    }
  }
  
  Type nll = 0.0;
  for (int t = 0; t < T; t++) {
    if (N(t) == 0) continue;
    // Compute normalization factor
    Type sum_exp = 0.0;
    for (int i = 0; i < P; i++) {
      sum_exp += exp(log_abundance(i, t) - max_col(t));
    }
    
    for (int i = 0; i < P; i++) {
      Type k = xj(i,t);
      if (N(t) > 0) {
        Type p = exp(log_abundance(i, t) - max_col(t)) / sum_exp;
        if (p < Type(1e-15)) p = Type(1e-15);
        if (p > Type(1 - 1e-15)) p = Type(1 - 1e-15);
        
        // log binomial pmf
        Type val = lchoose(N(t), k) + k*log(p) + (N(t)-k)*log(1-p);
        nll -= val;
      }
    }
  }
  
  // Add a small penalty on fi to stabilize
  if (P > 1) {
    Type mean_f = 0.0;
    for (int i = 0; i < P; i++) mean_f += fi(i);
    mean_f /= P;
    
    Type var_f = 0.0;
    for (int i = 0; i < P; i++) {
      Type diff = fi(i)-mean_f;
      var_f += diff*diff;
    }
    var_f /= (P-1);
    nll += sqrt(var_f);
  }
  
  return nll;
}
