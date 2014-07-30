#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <numeric>

using namespace std;

template<class T>
vector<T> operator + (vector<T> vec1, vector<T> vec2){
  for(int p = 0; p < vec1.size(); p++) vec1[p] += vec2[p];
  return vec1;
}
template<class T>
vector<T> operator + (double scalar, vector<T> vec){
  for(int p = 0; p < vec.size(); p++) vec[p] += scalar;
  return vec;
}
template<class T>
vector<T> operator + (vector<T> vec, double scalar){
  for(int p = 0; p < vec.size(); p++) vec[p] += scalar;
  return vec;
}
template<class T>
vector<T> operator - (vector<T> vec1, vector<T> vec2){
  for(int p = 0; p < vec1.size(); p++) vec1[p] -= vec2[p];
  return vec1;
}
template<class T>
vector<T> operator - (double scalar, vector<T> vec){
  for(int p = 0; p < vec.size(); p++) vec[p] -= scalar;
  return vec;
}
template<class T>
vector<T> operator - (vector<T> vec, double scalar){
  for(int p = 0; p < vec.size(); p++) vec[p] -= scalar;
  return vec;
}
template<class T>
vector<T> operator * (double scalar, vector<T> vec){
  for(int p = 0; p < vec.size(); p++) vec[p] *= scalar;
  return vec;
}
template<class T>
vector<T> operator * (vector<T> vec, double scalar){
  for(int p = 0; p < vec.size(); p++) vec[p] *= scalar;
  return vec;
}

template<class T>
void vec_disp(vector<T> vec){
  //  vector<double>::iterator vp;
  //  for(vp = vec.begin(); vp < vec.end(); vp++) cout << vec[vp] << " ";
  for(int p = 0; p < vec.size(); p++) cout << scientific << vec[p] << " ";
  cout << "\n";
}

