#include<RcppArmadillo.h>
#include<Rcpp.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
List FPC_Kern_S(NumericVector x, NumericVector t, IntegerVector N, double h1, double h2){
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




