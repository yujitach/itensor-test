#include "itensor/all.h"
#include <iostream>
#include <sstream>
#include <filesystem>
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

void go(int N,double gam,double a2){
    auto sites = CustomSpin(N,{"ConserveQNs=",false,"2S=",3});
    
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
        
    
    auto sweeps = Sweeps(60); // was 30
    sweeps.maxdim() = 10,20,100,100,200,200,300,300,400;
    sweeps.cutoff() = 1E-10;
    sweeps.niter() = 2;
    sweeps.noise() = 1E-7,1E-8,0.0;
    
    auto [en,psi] = dmrg(H,randomMPS(sites),sweeps,{"Silent=",true});

    auto wfs = std::vector<MPS>(1);
    wfs.at(0) = psi;

    auto [en1,psi1] = dmrg(H,wfs,randomMPS(sites),sweeps,{"Silent=",true,"Weight=",20.0});
    
//    auto enX=inner(psi1,H,psi1);
    
    //  std::cerr<<en1<<":"<<enX<<std::endl;
    
    std::string foo=format("%d %.10f %.10f ",N,gam,a2);
    std::ofstream ofs(foo+".excitedresult");
    ofs<<"{"<<N<<","<<gam<<","<<a2<<","<<en<<","<<en1<<"},"<<std::flush;
    
    writeToFile(foo+".psi",psi);
    writeToFile(foo+".excitedpsi",psi1);
    writeToFile(foo+".sites",sites);
}
namespace fs= std::filesystem;
int main(int argc,char**argv){
    for(;;){
        bool found=0;
        std::vector<std::filesystem::path> f;
        random_selector<> sel;
        for(const auto& entry : fs::directory_iterator(".")){
            if(entry.path().extension()==std::string(".excitedcompute")){
                f.push_back(entry.path());
            }
            found=1;
        }
        if(found){
            auto path=sel(f);
            std::istringstream is(path.stem());
            double gam,a2;
            int N;
            is>>N>>gam>>a2;
            std::cerr << "processing excited "<<N<<" "<<gam << " " << a2 << std::endl;
            fs::remove(path);
            go(N,gam,a2);
        }else{
            break;
        }
    }
    return 0;
}

