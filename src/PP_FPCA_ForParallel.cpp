#include<RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
List FPC_Kern_S(NumericVector x, NumericVector t, IntegerVector N, double h1, double h2) {
  int nx = x.size();
  int nn = N.size();
  IntegerVector idx = cumsum(N);
  idx.insert(idx.begin(), 0);
  rowvec f(nx, fill::zeros); 
  mat G(nx, nx, fill::zeros); 
  for (int k = 0; k < nn; k++) {
    NumericVector t0 = t[seq(idx[k], idx[k + 1] - 1)];
    int nt0 = t0.size();
    NumericVector z = rep(t0, nx) - rep_each(x, nt0);
    mat D1 = dnorm(z, 0, h1);
    D1.reshape(nt0, nx);
    f += sum(D1, 0);
    if (N(k) > 1) {
      mat D2 = dnorm(z, 0, h2);      
      D2.reshape(nt0, nx);  
      rowvec v = sum(D2, 0);
      G += v.t() * v - D2.t() * D2;
    }
  }
  return List::create(Named("f_mu") = f, Named("Cov_G") = G);
}

// [[Rcpp::export]]
List FPC_Kern_S_old(NumericVector x, NumericVector t, IntegerVector N, double h1, double h2){
  double tmp;
  colvec tmp1,tmp2;
  
  int nx = x.size();
  int nn = N.size();
  IntegerVector N_sum = cumsum(N); 
  N_sum.insert(N_sum.begin(),0); // insert 0 at the beginning
  vec f_mu(nx);
  mat g(nx,nx);
  
  // mean density function
  for(int i=0; i<nx; i++)
    f_mu[i] = sum(dnorm((x[i]-t)/h1))/h1;
  
  // covariance function
  for(int i=0;i<nx;i++){
    for(int j=0;j<nx;j++){
      g(i,j) = 0.0;
      if(i<=j){
        for(int k=0;k<nn;k++){
          if(N(k)>1){
            tmp1    = dnorm((x[i]-t[seq(N_sum[k],N_sum[k+1]-1)])/h2);
            tmp2    = dnorm((x[j]-t[seq(N_sum[k],N_sum[k+1]-1)])/h2);
            tmp     = sum(sum(tmp1*tmp2.t(),0))-sum(tmp1%tmp2);
            g(i,j) += tmp;
          }
        }
        g(i,j) = g(i,j)/h2/h2;
      }
    }
  }
  g = symmatu(g);
  
  return List::create(Named("f_mu")=f_mu,Named("Cov_G")=g);
}
