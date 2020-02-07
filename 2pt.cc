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
        int i=1;
        auto op_i = op(sites,"Sz",i);
        auto op_j = op(sites,"Sz",j);
        psi.position(j);
        auto psidag = dag(psi);

        //Prime the link indices to make them distinct from
        //the original ket links
        psidag.prime("Link");

        //index linking i-1 to i:
        auto li_1 = leftLinkIndex(psi,i);

        auto C = prime(psi(i),li_1)*op_i;
        C *= prime(psidag(i),"Site");
        for(int k = i+1; k < j; ++k)
            {
            C *= psi(k);
            C *= psidag(k);
            }
        //index linking j to j+1:
        auto lj = rightLinkIndex(psi,j);

        C *= prime(psi(j),lj)*op_j;
        C *= prime(psidag(j),"Site");

        auto result = elt(C);

        printfln("{%d,  %.10f},",j,result);
    }
    return 0;
}

