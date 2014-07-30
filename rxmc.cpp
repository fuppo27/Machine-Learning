#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <numeric>
#include "probability.hpp"

#define propose propose_via_Metropolis_algorithm

using namespace std;

int Ds, Dd, Dn, Df;
double Alpha, Mu;
vector<double> Mu1, Mu2;
vector< vector<double> > MatX;

void initValue(){
  MatX.assign(Ds, vector<double>(Dd, 0.0));
  cout << MatX.size() << " " << MatX[0].size() << endl;
}

void calRXMC(int Ds, int Dn, int Df, double alpha, double mu){
  int i, s, c1, c2;
  double e, v;
  vector<double> mu1(Dd, 0.0), mu2(Dd, 0.0);

  vector< vector<double> > x(Dn, vector<double>(Dd, 0.0));
  for(i = 0; i < Dn; i++) x[i] = random_multivariate_uniform(Dd, -mu, mu);
  vector<Gaussian_mixture> GM(Dn);
  for(i = 0, e = fabs(mu)/(Dn-1), v = 0.0; i < Dn; i++){
    mu1.assign(Dd, v);
    mu2.assign(Dd, -v);
    GM[i].config_params(alpha, mu1, mu2);
    v += e;
  }

  for(s = 0; s < Ds; s++){
    for(i = 0; i < Dn; i++) x[i] = propose(x[i], GM[i]);
    if(s%Df == 0){
      c1 = rand()%Dn;
      c2 = (c1+1)%Dn;
      if(random_uniform() < ((GM[c1].calcPDF(x[c2])*GM[c2].calcPDF(x[c1])) / (GM[c1].calcPDF(x[c1])*GM[c2].calcPDF(x[c2])))) x[c1].swap(x[c2]);
    }
    MatX[s] = x[Dn-1];
  }
}

void printValue(char *fn1){
  ofstream fout;
  fout.open(fn1);
  for(int s = 0; s < Ds; s++){
    for(int d = 0; d < Dd; d++) fout << scientific << MatX[s][d] << " ";
    fout << endl;
  }
  fout.close();
}

int main(int argc, char* argv[]){
  srand(time(0));
  Ds = atoi(argv[1]);//#simulations[10000]
  Dd = atoi(argv[2]);//次元[2, 3]
  Dn = atoi(argv[3]);//レプリカ数[5]
  Df = atoi(argv[4]);//交換頻度[20]
  Alpha = atof(argv[5]);//mixture weight[0.5]
  Mu = atof(argv[6]);//paramater basis[3]
  initValue();
  calRXMC(Ds, Dn, Df, Alpha, Mu);
  printValue(argv[7]);
}
