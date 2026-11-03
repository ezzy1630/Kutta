#include <array>
#include <cmath>
#include <algorithm>
#include <vector>
#include <chrono>
#include <iostream>
#include <cstring>
#include <random>
constexpr int terms=10;
using Row=std::array<double,terms>;using Matrix=std::array<Row,terms>;
Matrix before(Matrix a) {
  Matrix q{};for(int i=0;i<terms;++i) q[i][i]=1;
  for(int iteration=0;iteration<terms*terms*16;++iteration) {
    int p=0,r=1;double scale=0;
    for(int i=0;i<terms;++i) {
      scale=std::max(scale,std::abs(a[i][i]));
      for(int j=i+1;j<terms;++j) if(std::abs(a[i][j])>std::abs(a[p][r])) {p=i;r=j;}
    }
    if(std::abs(a[p][r])<1e-13*scale) break;
    const double angle=.5*std::atan2(2*a[p][r],a[r][r]-a[p][p]);
    const double c=std::cos(angle),s=std::sin(angle),app=a[p][p],arr=a[r][r],apr=a[p][r];
    a[p][p]=c*c*app-2*c*s*apr+s*s*arr;a[r][r]=s*s*app+2*c*s*apr+c*c*arr;a[p][r]=a[r][p]=0;
    for(int i=0;i<terms;++i) {
      if(i!=p && i!=r) {const double x=a[i][p],y=a[i][r];a[i][p]=a[p][i]=c*x-s*y;a[i][r]=a[r][i]=s*x+c*y;}
      const double x=q[i][p],y=q[i][r];q[i][p]=c*x-s*y;q[i][r]=s*x+c*y;
    }
  }
  double largest=0;for(int i=0;i<terms;++i) largest=std::max(largest,a[i][i]);
  Matrix out{};
  for(int k=0;k<terms;++k) if(a[k][k]>largest*1e-11)
    for(int i=0;i<terms;++i) for(int j=0;j<terms;++j) out[i][j]+=q[i][k]*q[j][k]/a[k][k];
  return out;
}
Matrix after(Matrix a) {
  Matrix q{};for(int i=0;i<terms;++i) q[i][i]=1;
  for(int iteration=0;iteration<terms*terms*16;++iteration) {
    int p=0,r=1;double scale=0,maximum=std::abs(a[0][1]);
    for(int i=0;i<terms;++i) {
      scale=std::max(scale,std::abs(a[i][i]));
      for(int j=i+1;j<terms;++j) {const double candidate=std::abs(a[i][j]);if(candidate>maximum) {p=i;r=j;maximum=candidate;}}
    }
    if(maximum<1e-13*scale) break;
    const double angle=.5*std::atan2(2*a[p][r],a[r][r]-a[p][p]);
    const double c=std::cos(angle),s=std::sin(angle),app=a[p][p],arr=a[r][r],apr=a[p][r];
    a[p][p]=c*c*app-2*c*s*apr+s*s*arr;a[r][r]=s*s*app+2*c*s*apr+c*c*arr;a[p][r]=a[r][p]=0;
    for(int i=0;i<terms;++i) {
      if(i!=p && i!=r) {const double x=a[i][p],y=a[i][r];a[i][p]=a[p][i]=c*x-s*y;a[i][r]=a[r][i]=s*x+c*y;}
      const double x=q[i][p],y=q[i][r];q[i][p]=c*x-s*y;q[i][r]=s*x+c*y;
    }
  }
  double largest=0;for(int i=0;i<terms;++i) largest=std::max(largest,a[i][i]);
  Matrix out{};
  for(int k=0;k<terms;++k) if(a[k][k]>largest*1e-11)
    for(int i=0;i<terms;++i) for(int j=0;j<terms;++j) out[i][j]+=q[i][k]*q[j][k]/a[k][k];
  return out;
}

int main() {
 std::mt19937 rng(75);std::vector<Matrix> matrices(512);
 for(auto& a:matrices) {
  Matrix r{};for(auto& row:r) for(auto& v:row) v=(int(rng()%2001)-1000)*.001;
  for(int i=0;i<terms;++i) for(int j=0;j<terms;++j) for(int k=0;k<terms;++k) a[i][j]+=r[k][i]*r[k][j];
  const auto b=before(a),c=after(a);if(std::memcmp(&b,&c,sizeof b)) return 2;
 }
 std::cout<<"parity_matrices="<<matrices.size()<<std::endl;
 volatile double checksum=0;
 for(int rep=0;rep<5;++rep) for(int order=0;order<2;++order) {
  const bool candidate=(order+rep)%2;auto start=std::chrono::steady_clock::now();
  Matrix (*volatile operation)(Matrix)=candidate?after:before;
  for(int repeat=0;repeat<3;++repeat) for(const auto& a:matrices) {auto b=operation(a);for(const auto& row:b) for(double value:row) checksum=checksum+value;}
  std::cout<<"INVERSE_TRIAL rep="<<rep<<" candidate="<<candidate<<" ms="<<std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-start).count()<<std::endl;
 }
}
