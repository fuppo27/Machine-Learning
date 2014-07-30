#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <numeric>
#include "vector_ope.hpp"

#define PI 3.14159265358979

double standard_multivariate_normal_pdf(std::vector<double> x){
  double coef, ker;
  ker = exp(-0.5*inner_product(x.begin(), x.end(), x.begin(), 0.0));
  coef = pow(1.0/(2.0*PI), x.size()/2.0);
  return coef*ker;
}

class Gaussian_mixture{
private:
  double Alpha;
  std::vector<double> Mu1, Mu2, X;

public:
  Gaussian_mixture(){}
  Gaussian_mixture(double alpha, std::vector<double> mu1, std::vector<double> mu2, std::vector<double> x){ Alpha = alpha; Mu1 = mu1; Mu2 = mu2; X = x; }
  Gaussian_mixture(double alpha, std::vector<double> mu1, std::vector<double> mu2){ Alpha = alpha; Mu1 = mu1; Mu2 = mu2; }
  Gaussian_mixture(std::vector<double> mu1, std::vector<double> mu2, std::vector<double> x){ Alpha = 0.5; Mu1 = mu1; Mu2 = mu2; X = x; }
  Gaussian_mixture(std::vector<double> mu1, std::vector<double> mu2){ Alpha = 0.5; Mu1 = mu1; Mu2 = mu2; }

  void config_params(double alpha, std::vector<double> mu1, std::vector<double> mu2){ Alpha = alpha; Mu1 = mu1; Mu2 = mu2; }
  double calcPDF(std::vector<double> x){
    double p1 = standard_multivariate_normal_pdf(x-Mu1);
    double p2 = standard_multivariate_normal_pdf(x-Mu2);
    return Alpha*p1 + (1.0-Alpha)*p2;
  }
  double calcPDF(){
    double p1 = standard_multivariate_normal_pdf(X-Mu1);
    double p2 = standard_multivariate_normal_pdf(X-Mu2);
    return Alpha*p1 + (1.0-Alpha)*p2;
  }

  double operator () (std::vector<double> x){
    double p1 = standard_multivariate_normal_pdf(x-Mu1);
    double p2 = standard_multivariate_normal_pdf(x-Mu2);
    return Alpha*p1 + (1.0-Alpha)*p2;
  }
};


double random_normal(double mu = 0.0, double sigma = 1.0){
    double u1, u2, z1, z2;
    u1=(double)rand()/RAND_MAX;
    u2=(double)rand()/RAND_MAX;
    z1=sigma*sqrt(-2.0*log(u1))*cos(2.0*PI*u2)+mu;
    z2=sigma*sqrt(-2.0*log(u1))*sin(2.0*PI*u2)+mu;
    if(rand()%2 == 0) return z1;
    else return z2;
}

double random_uniform(double min = 0.0, double max = 1.0){
  return ((double)rand()/((double)RAND_MAX+1.0))*(max-min)+min;
}

std::vector<double> random_multivariate_normal(int dim, double mu = 0.0, double sigma = 1.0){
  std::vector<double> x(dim, 0.0);
  for(int i = 0; i < dim; i++) x[i] = random_normal(mu, sigma);
  return x;
}

std::vector<double> random_multivariate_uniform(int dim, double min = 0.0, double max = 1.0){
  std::vector<double> x(dim, 0.0);
  for(int i = 0; i < dim; i++) x[i] = random_uniform(min, max);
  return x;
}

template<class T>
std::vector<double> propose_via_Metropolis_algorithm(std::vector<double> x, T distribution){
  //param[2]分布のファンクタ
  std::vector<double> proposed = x + 0.5*random_multivariate_normal(x.size());
  /*
  distribution(proposed);
  distribution(x);
  distribution.calcPDF(proposed);
  distribution.calcPDF(x);
  */
  double ratio = distribution.calcPDF(proposed)/distribution.calcPDF(x);
  double ru = random_uniform();
  if(ru < ratio) return proposed;
  else return x;
}
