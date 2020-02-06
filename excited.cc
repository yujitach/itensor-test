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
    
    auto ampo = AutoMPO(sites);
    for(int j = 1; j <= N-1; ++j){
        ampo += -1,"Sz",j,"Sz",j+1;
    }
    ampo+= -.5,"Sz",1;
    ampo+= -.5,"Sz",N;
    for(int j = 1; j <= N; ++j){
        ampo += (1+a2),"Sz",j,"Sz",j;
        ampo += -gam/2,"S+",j;
        ampo += -gam/2,"S-",j;
    }
    auto H = toMPO(ampo);

    auto en=inner(psi,H,psi);
    
    std::cout<<en<<std::endl;
    
    auto wfs = std::vector<MPS>(1);
    wfs.at(0) = psi;
    
    auto sweeps = Sweeps(60); // was 30
    sweeps.maxdim() = 10,20,100,100,200,200,300,300,400;
    sweeps.cutoff() = 1E-10;
    sweeps.niter() = 2;
    sweeps.noise() = 1E-7,1E-8,0.0;

    
    auto [en1,psi1] = dmrg(H,wfs,randomMPS(sites),sweeps,{"Silent=",true,"Weight=",20.0});

    auto enX=inner(psi1,H,psi1);
    
    std::cerr<<en1<<":"<<enX<<std::endl;
    
    std::string foo=format("%d %.10f %.10f ",N,gam,a2);
    std::ofstream ofs(foo+".excitedresult");
    ofs<<"{"<<N<<","<<gam<<","<<a2<<","<<en<<","<<en1<<"},"<<std::flush;

    writeToFile(foo+".excitedpsi",psi1);
    writeToFile(foo+".excitedsites",sites);

    return 0;
}

