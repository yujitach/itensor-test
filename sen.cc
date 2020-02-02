#include "itensor/all.h"
#include <filesystem>
#include <iostream>
#include <sstream>
#include <vector>
#include <random>
#include <iterator>

// taken from https://gist.github.com/cbsmith/5538174
// from here
template <typename RandomGenerator = std::default_random_engine>
struct random_selector
{
	//On most platforms, you probably want to use std::random_device("/dev/urandom")()
	random_selector(RandomGenerator g = RandomGenerator(std::random_device("/dev/urandom")()))
		: gen(g) {}

	template <typename Iter>
	Iter select(Iter start, Iter end) {
		std::uniform_int_distribution<> dis(0, std::distance(start, end) - 1);
		std::advance(start, dis(gen));
		return start;
	}

	//convenience function
	template <typename Iter>
	Iter operator()(Iter start, Iter end) {
		return select(start, end);
	}

	//convenience function that works on anything with a sensible begin() and end(), and returns with a ref to the value type
	template <typename Container>
	auto operator()(const Container& c) -> decltype(*begin(c))& {
		return *select(begin(c), end(c));
	}

private:
	RandomGenerator gen;
};
// to here

using namespace itensor;

inline int mod(int x,int N){
    if(x>N)
        return x-N;
    return x;
}

void go(int N,double gam,double a2){
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
    auto sweeps = Sweeps(40); // was 30
    sweeps.maxdim() = 10,20,100,100,200,200,300;
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
		
		std::string foo=format("%d %.10f %.10f ",N,gam,a2);
		std::ofstream ofs(foo+".result");
		ofs<<"{"<<N<<","<<gam<<","<<a2<<","<<SvN<<","<<result<<"},"<<std::flush;
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
		std::vector<std::filesystem::path> f;
		random_selector<> sel;
		for(const auto& entry : fs::directory_iterator(".")){
			if(entry.path().extension()==std::string(".compute")){
				f.push_back(entry.path());
			}
			found=1;
		}
		if(found){
			auto path=sel(f);
			std::istringstream is(path.stem());
			double gam,a2;
			is>>gam>>a2;
			std::cerr << "processing "<<gam << " " << a2 << std::endl;
			fs::remove(path);
			go(40,gam,a2);			
			go(60,gam,a2);			
		}else{
			break;
		}
	}
    return 0;
}

