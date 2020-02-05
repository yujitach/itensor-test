#include "itensor/all.h"
#include <iostream>
#include <sstream>

using namespace itensor;

int main(int argc,char**argv){
    if(argc!=2){
        exit(-1);
    }
    char*p=argv[1];
    std::string s=p;
    std::istringstream is(s);
    int N;
    double gam,a2;
    is>>N>>gam>>a2;
    CustomSpin sites;
    readFromFile(s+".sites",sites);
    MPS psi;
    readFromFile(s+".psi",psi);
    for(auto b=1;b<N;b++){
        psi.position(b);

        //SVD this wavefunction to get the spectrum
        //of density-matrix eigenvalues
        auto l = leftLinkIndex(psi,b);
        auto s = siteIndex(psi,b);
        auto [U,S,V] = svd(psi(b),{l,s});
        auto u = commonIndex(U,S);

        //Apply von Neumann formula
        //to the squares of the singular values
        Real SvN = 0.;
        for(auto n : range1(dim(u)))
            {
            auto Sn = elt(S,n,n);
            auto p = sqr(Sn);
            if(p > 1E-12) SvN += -p*log(p);
            }
        printfln("{%d,  %.10f},",b,SvN);
    }
    return 0;
}

