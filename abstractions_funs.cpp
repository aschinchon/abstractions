#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]


// This function allows to do convolutions on the edge elements of the 
// matrix since it extracts elements from the opposite rows or columns 
// to avoid losing dimensionality after convolution
int get_index(int M, int i)
{
  if (i < 0)
    return (M + i % M) % M;
  if(i >= M)
    return i % M;
  return i;
}

// Sensor stage
Rcpp::DataFrame sensor (arma::mat envM,
                        Rcpp::DataFrame parF,
                        double FL,
                        double FR,
                        double RA,
                        double SO) {
  
  int m = envM.n_rows;
  int n = envM.n_cols;
  int npart = parF.nrows();
  
  Rcpp::IntegerVector x = parF["x"];
  Rcpp::IntegerVector y = parF["y"];
  Rcpp::NumericVector h = parF["h"];
  
  // create a vector to store heading
  Rcpp::NumericVector hnew(npart);

  for(int i = 0; i < npart; i++) {
    
    // Sensory stage
    int Fx = get_index(n, (int)round(x[i] + SO * cos(h[i])));
    int Fy = get_index(m, (int)round(y[i] + SO * sin(h[i])));
    
    int FLx = get_index(n, (int)round(x[i] + SO * (cos(h[i] + FL))));
    int FLy = get_index(m, (int)round(y[i] + SO * (sin(h[i] + FL))));
    
    int FRx = get_index(n, (int)round(x[i] + SO * (cos(h[i] + FR))));
    int FRy = get_index(m, (int)round(y[i] + SO * (sin(h[i] + FR))));
    
    double F  = envM(Fx, Fy);
    double FL = envM(FLx, FLy);
    double FR = envM(FRx, FRy);
    
    if ((F > FL) && (F > FR))
    {
      hnew[i] = h[i];
    }
    else if ((F < FL) && (F < FR))
    {
      if ( rand() % 2 == 0 )
      {
        hnew[i] = h[i] + RA;
      }
      else
      {
        hnew[i] = h[i] - RA;
      }
      
    }
    else if (FL < FR)
    {
      hnew[i] = h[i] - RA;
    }
    else if (FR < FL)
    {
      hnew[i] = h[i] + RA;
    }
    else
    {
      hnew[i] = h[i];
    }
  }
  
  Rcpp::DataFrame dest = Rcpp::DataFrame::create(Rcpp::Named("x") = x,
                                                 Rcpp::Named("y") = y,
                                                 Rcpp::Named("h") = hnew);
  
  return dest;
  
}

// Motor stage
Rcpp::DataFrame motor (Rcpp::DataFrame parF,
                       int m,
                       int n,
                       double SS){
  
  int npart = parF.nrows();

  Rcpp::IntegerVector x = parF["x"];
  Rcpp::IntegerVector y = parF["y"];
  Rcpp::NumericVector h = parF["h"];
  
  
  // create the new columns
  Rcpp::IntegerVector xnew(npart);
  Rcpp::IntegerVector ynew(npart);
  
  for(int i = 0; i < npart; i++) {
    
    xnew[i] = get_index(n, (int)round(x[i] + SS * cos(h[i])));
    ynew[i] = get_index(m, (int)round(y[i] + SS * sin(h[i])));
    
    
  }
  
  Rcpp::DataFrame dest = Rcpp::DataFrame::create(Rcpp::Named("x") = xnew,
                                                 Rcpp::Named("y") = ynew,
                                                 Rcpp::Named("h") = h);
  
  return dest;
}

// Deposition
arma::mat deposition (Rcpp::DataFrame parF,
                      double depT,
                      arma::mat envM) {
  
  int npart = parF.nrows();
  
  Rcpp::IntegerVector x = parF["x"];
  Rcpp::IntegerVector y = parF["y"];

  for(int i = 0; i < npart; i++) {
    envM(x[i], y[i]) = envM(x[i], y[i]) + depT;    
    
    
  }
    
return envM;
}

// Evaporation
arma::mat evaporate (arma::mat source,
                     double factor) {
  
  int m = source.n_rows;
  int n = source.n_cols;
  
  arma::mat dest(m, n);
  
  for (int y = 0; y < n; ++y) {
    for (int x = 0; x < m; ++x) {
      dest(x, y) = source(x, y) * (1 - factor);
    }
  }
  return dest;
};

// Gather all into a funtion called physarum
// [[Rcpp::export]]
arma::mat physarum(arma::mat envM,
                 Rcpp::DataFrame parF,
                 double decayT,
                 double FL,
                 double FR,
                 double RA,
                 double SO,
                 int SS,
                 double depT,
                 int iters)
{
  int m = envM.n_rows;
  int n = envM.n_cols;
  
  // envM = scale(envM);
  
  for (int k = 0; k < iters; k++){

    // Sensor stage
    parF = sensor (envM, parF, FL, FR, RA, SO);
    
    // Motor stage
    parF = motor (parF, m, n, SS);
    
    // Deposition stage
    envM = deposition (parF, depT, envM);

    // Evaporation
    envM = evaporate(envM, decayT);
    
  }  

  return envM;
};


