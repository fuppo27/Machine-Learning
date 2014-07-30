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

int Ds, Dd;
double Alpha, Mu;
vector<double> Mu1, Mu2;
vector< vector<double> > MatX;

void initValue(){
  Mu1.assign(Dd, Mu);
  Mu2.assign(Dd, -Mu);
  MatX.assign(Ds, vector<double>(Dd, 0.0));
  cout << MatX.size() << endl;
}

void calMCMC(int Ds, double alpha, vector<double> mu1, vector<double> mu2){
  Gaussian_mixture gm(alpha, mu1, mu2);
  vector<double> x(Dd, 0.0);
  //  vec_disp(x);
  for(int s = 0; s < Ds; s++){
    x = propose(x, gm);
    //    vec_disp(x);
    MatX[s] = x;
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
  Ds = atoi(argv[1]);//#simulations
  Dd = atoi(argv[2]);//次元 
  Alpha = atof(argv[3]);//mixture weight
  Mu = atof(argv[4]);//paramater basis
  initValue();
  calMCMC(Ds, Alpha, Mu1, Mu2);
  printValue(argv[5]);
}
