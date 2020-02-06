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
    for(auto j=1;j<N;j++){
        auto op_j = op(sites,"Sz",j);
        psi.position(j);
        auto ket=psi(j);
        auto bra=dag(prime(ket,"Site"));
        auto result=elt(bra*op_j*ket);
        printfln("{%d,  %.10f},",j,result);
    }
    return 0;
}

