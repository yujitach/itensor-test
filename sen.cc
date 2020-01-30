#include "itensor/all.h"
#include <filesystem>
#include <iostream>
#include <sstream>

using namespace itensor;

inline int mod(int x,int N){
    if(x>N)
        return x-N;
    return x;
}

void go(double gam,double a2){
    int N = 40;
    auto sites = CustomSpin(N,{"ConserveQNs=",false,"2S=",3});
    //
    // Factors of 4 and 2 are to rescale
    // spin operators into Pauli matrices
    //
/*	if(argc!=3){
		exit(-1);
	}
*/	
    auto ampo = AutoMPO(sites);
    for(int j = 1; j <= N; ++j){
        ampo += -1,"Sz",j,"Sz",mod(j+1,N);
        ampo += (1+a2),"Sz",j,"Sz",j;
    } 
    for(int j = 1; j <= N; ++j){
        ampo += -gam/2,"S+",j;
        ampo += -gam/2,"S-",j;
    }
    auto H = toMPO(ampo);
    //
    // Set the parameters controlling the accuracy of the DMRG
    // calculation for each DMRG sweep. 
    //
    auto sweeps = Sweeps(20); // was 30
    sweeps.maxdim() = 10,20,100,100,200;
    sweeps.cutoff() = 1E-10;
    sweeps.niter() = 2;
    sweeps.noise() = 1E-7,1E-8,0.0;
    //
    // Begin the DMRG calculation
    // for the ground state
    //
    auto [en,psi] = dmrg(H,randomMPS(sites),sweeps,{"Silent=",true});

	auto b=N/2;
//    for(auto b=1;b<N;b++){
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
//        printfln("{%d,  %.10f},",b,SvN);
		
		
		int i=1;
		int j=N/2;
		auto op_i = op(sites,"Sz",i);
		auto op_j = op(sites,"Sz",j);

		//'gauge' the MPS to site i
		//any 'position' between i and j, inclusive, would work here
		psi.position(i); 

		//Create the bra/dual version of the MPS psi
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

		auto result = elt(C); //or eltC(C) if expecting complex
		
		std::string foo=format("%.10f %.10f ",gam,a2);
		std::ofstream ofs(foo+".result");
		ofs<<"{"<<gam<<","<<a2<<","<<SvN<<","<<result<<"},"<<std::flush;
		writeToFile(foo+".sites",sites);
		writeToFile(foo+".psi",psi);
//        printfln("{%f, %f,  %.10f},",a2,gam,SvN);
//    }
}

//auto gam=atof(argv[1]); //0.3702;
//auto a2=atof(argv[2]); // -0.0661;

namespace fs= std::filesystem;
int main(int argc,char**argv){
	for(;;){
		bool found=0;
		for(const auto& entry : fs::directory_iterator(".")){
			if(entry.path().extension()==std::string(".compute")){
				std::cerr<<entry.path()<<std::endl;
				std::istringstream is(entry.path().stem());
				double gam,a2;
				is>>gam>>a2;
				std::cerr << "processing "<<gam << " " << a2 << std::endl;
				fs::remove(entry);
				go(gam,a2);
				found=1;
				break;
			}
		}
		if(!found){
			break;
		}
	}
    return 0;
}

