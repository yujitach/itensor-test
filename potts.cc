#include "itensor/all.h"
#include <cmath>

using namespace itensor;

inline int mod(int x,int N){
    if(x>N)
        return x-N;
    return x;
}
int main(){       
    int N = 100;
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
        ampo+=f*0.7071067811865475*4,"Sy",mod(i+0,2*N),"Sy",mod(i+1,2*N);
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

    for(auto b=1;b<=2*N;b++){
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

