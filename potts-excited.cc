#include "itensor/all.h"
#include <cmath>
#include <cstdlib>

using namespace itensor;

inline int mod(int x,int N){
    if(x>N)
        return x-N;
    return x;
}
int main(int argc,char**argv){
    int N = std::atoi(argv[1]);
    auto sites = SpinHalf(2*N,{"ConserveQNs=",false});
    //
    // Factors of 4 and 2 are to rescale
    // spin operators into Pauli matrices
    //
    Real h=10000;
    Real f=1;
    auto ampo = AutoMPO(sites);
    for(int j = 1; j <= N; ++j){
        int i=2*j-1;
        ampo+=f*0.5*2,"Sz",mod(i+1,2*N);
        ampo+=-f*0.7071067811865475*4,"Sx",mod(i+0,2*N),"Sx",mod(i+1,2*N);
        ampo+=-f*0.7071067811865475*4,"ISy",mod(i+0,2*N),"ISy",mod(i+1,2*N);
        ampo+=-f*0.5*4,"Sz",mod(i+0,2*N),"Sz",mod(i+1,2*N);
        
        ampo+=-0.375*4,"Sx",mod(i+1,2*N),"Sx",mod(i+3,2*N);
        ampo+=0.375*8,"Sx",mod(i+1,2*N),"Sz",mod(i+2,2*N),"Sx",mod(i+3,2*N);
        ampo+=-0.125*4,"Sz",mod(i+1,2*N),"Sz",mod(i+3,2*N);
        ampo+=-0.25*4,"Sz",mod(i+1,2*N),"Sz",mod(i+2,2*N);
        ampo+=-0.125*8,"Sz",mod(i+1,2*N),"Sz",mod(i+2,2*N),"Sz",mod(i+3,2*N);
        ampo+=-0.25*4,"Sz",mod(i+0,2*N),"Sz",mod(i+3,2*N);
        ampo+=-0.5*4,"Sz",mod(i+0,2*N),"Sz",mod(i+2,2*N);
        ampo+=-0.25*8,"Sz",mod(i+0,2*N),"Sz",mod(i+2,2*N),"Sz",mod(i+3,2*N);
        ampo+=0.375*8,"Sz",mod(i+0,2*N),"Sx",mod(i+1,2*N),"Sx",mod(i+3,2*N);
        ampo+=-0.375*16,"Sz",mod(i+0,2*N),"Sx",mod(i+1,2*N),"Sz",mod(i+2,2*N),"Sx",mod(i+3,2*N);
        ampo+=-0.125*8,"Sz",mod(i+0,2*N),"Sz",mod(i+1,2*N),"Sz",mod(i+3,2*N);
        ampo+=-0.25*8,"Sz",mod(i+0,2*N),"Sz",mod(i+1,2*N),"Sz",mod(i+2,2*N);
        ampo+=-0.125*16,"Sz",mod(i+0,2*N),"Sz",mod(i+1,2*N),"Sz",mod(i+2,2*N),"Sz",mod(i+3,2*N);
        
        ampo+=-h*0.25*2,"Sz",mod(i+1,2*N);
        ampo+=h*0.25*2,"Sz",mod(i+0,2*N);
        ampo+=-h*0.25*4,"Sz",mod(i+0,2*N),"Sz",mod(i+1,2*N);
        ampo+=h*0.25*4,"Sz",mod(i+0,2*N),"Sz",mod(i+0,2*N);
    }
    auto H = toMPO(ampo);
    //
    // Set the parameters controlling the accuracy of the DMRG
    // calculation for each DMRG sweep.
    //
    auto sweeps = Sweeps(60); // was 30
    sweeps.maxdim() = 10,20,100,100,200,200,300,300,400;
    sweeps.cutoff() = 1E-10;
    sweeps.niter() = 2;
    sweeps.noise() = 1E-7,1E-8,0.0;
    //
    // Begin the DMRG calculation
    // for the ground state
    //
    auto [en,psi] = dmrg(H,MPS(InitState(sites,"Up")),sweeps,{"Quiet=",true});
    auto wfs = std::vector<MPS>(1);
    wfs.at(0) = psi;

    auto [en1,psi1] = dmrg(H,wfs,MPS(InitState(sites,"Up")),sweeps,{"Quiet=",true,"Weight=",20.0});
    printfln("{%d , %.10f , %.10f  }",N,en,en1);
    return 0;
}

