#include <iostream>
#include<vector>
#include <cctype>
#include<map>
#include<random>
#include<list>
#include<time.h>
#include<fstream>
#include <sstream>
#include<algorithm>
#include <float.h>
#include <chrono>

using namespace std;
double PI=3.14159;

#include "setup.h"
#include "fitness_landscapes.h"
#include "clone_logic.h"
#include "write.h"
#include "polyh.h"


map<vector<int>,karyotype> m;
list<karyotype> mutants;

int main(int argc, char *argv[])
{

    float sd_mutation,mean_mutation;


    string config_file_path = "config.txt";
    //string config_file_path = "/Users/taolee/Documents/GBM/runtime_ver1/ABM/primray/10/P19/fitted_landscape_config.txt";
    if(argc<2) {
        cout << "using default config file (hopefully one is there...)" << endl;
    }else{
            config_file_path = argv[1];
    }
    //load parameters from config
    parameters par(config_file_path);
    // setup for output writing
    string write_dir = make_subdir(par.output_dir);
    write_log(par.p,0,0,par.dt,write_dir);
    std::ofstream cellfile;
    std::string fname=write_dir+"/summary.txt";
    cellfile.open(fname);

    // set random seeds
    std::random_device rd;
    std::mt19937 gen(rd());

    // instantiate cell population
    vector<int> k1=par.init_kary;

    // instantiate fitness landscape
    cout << "initializing landscape...";
    par.read_landscape_file();
    fitness_landscape f;
    if(par.fitness_landscape_type=="random"){
        f.init(par.peaks,par.wavelength,par.scale,par.centre,par.deltafu);

    }
    if(par.fitness_landscape_type=="krig"){
        f.init(par.knots,par.c,par.d,par.deltafu);
    }
    if(par.fitness_landscape_type=="Predicted"){
        f.init(par.knots,par.c,par.deltafu);
    }
    if(par.fitness_landscape_type=="flat"){
        f.init(1.0);
    }

    f.setNoise(par.fitness_noise);

    cout << " type = " << f.type <<". initial clone fitness: " << f.get_fitness(k1,gen) << endl;

    //instantiate cell population
    cout << "population file: " << par.population_file << endl;
    if(par.population_file=="not_supplied"){
        instantiate_population(k1,par.init_size,m,f,gen);
    }else{
        instantiate_population(par.population_file,m,f,gen);
    }

    cout << "starting sim..." << endl;

    // iterate over generations
    for(int i=0; i<par.Nsteps; i++){
        float mean_fitness=0;

        //count how many cells
        int ncells = 0;
        auto it = m.begin();
        for (auto it = m.begin(); it != m.end();){
                ncells+=(it->second).n;
                it++;
        }


        // perform cell divisions
        float min_fitness = it->second.fitness;
        float max_fitness = min_fitness;
        ncells=0;
        for (auto& [key, value] : m){
            value.divide(par.p,par.pgd,gen,mutants,par.dt);
            mean_fitness=(value.n*value.fitness+ncells*mean_fitness)/(ncells+value.n);
            min_fitness = min(min_fitness,value.fitness);
            max_fitness = max(max_fitness,value.fitness);
            ncells+=value.n;
        }
        // assign any mis-segregated cells to their clone
        for(auto cl:mutants){
            // if the clone already exists, just add one cell
            if (auto search = m.find(cl.cn); search != m.end()){
                (search->second).n+=1;
            }else{ //if this is a new clone then generate a new entry in the map
                // the following call to get fitness needs rethought so we can
                // call the same way regardless of the fitness landscape
                //float fitness = f.get_fitness(cl.cn,cl.fitness,gen);
                float fitness  = f.get_fitness(cl.cn,gen);
                karyotype c2(cl.cn,1,fitness);
                m[cl.cn]=c2;
            }
        }
        mutants.clear();


        cout << "Ncells: " << ncells << "; mean fitness: " << mean_fitness<< "; min fitness: " << min_fitness<< "; max fitness: " << max_fitness << "; Nclones: " << m.size() << endl;
        // perform passaging if cell threshold is exceeded
        if(ncells>par.max_size){
            // write out a sample of cells to file
            float p_output = (float)par.target_output_size/(float)ncells;
            write_pop(m, i, p_output, write_dir, gen);
            // delete cells at random
            for (auto it = m.begin(); it != m.end();){
                    poisson_distribution<> dpois((it->second).n*par.bf);
                    (it->second).n = dpois(gen);
                    if((it->second).n==0){
                        it = m.erase(it);
                    }else{++it;}
            }
        } else if(i%par.pop_write_freq==0){
            float p_output = (float)par.target_output_size/(float)ncells;
            write_pop(m, i, p_output, write_dir, gen);
        }

        if(i==0){
            float p_output = (float)par.target_output_size/(float)ncells;
            write_pop(m, i, p_output, write_dir, gen);
        }


        cellfile << ncells << endl;

    }
    cellfile.close();

    return 0;
}
