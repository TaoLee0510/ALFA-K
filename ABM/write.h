#ifndef WRITE_H_INCLUDED
#define WRITE_H_INCLUDED
#include <string>

#include <sys/types.h> // required for stat.h
#include <sys/stat.h> // no clue why required -- man pages say so

using namespace std;

string make_subdir(string parent){

    bool created = false;

    int i = 0;
    string dir_name;
    while(!created){
        string rep = to_string(i);
        rep.insert(rep.begin(), 5 - rep.size(), '0');
        dir_name = parent + "/" + rep;
        const char* c = dir_name.c_str();
        mode_t nMode = 0733; // UNIX style permissions
        int nError = 0;
        #if defined(_WIN32)
            nError = mkdir(c); // can be used on Windows
        #else
            nError = mkdir(c,nMode); // can be used on non-Wind$
        #endif
        if(nError == 0){
                created=true;
        }
        ++i;
    }
    return(dir_name);
}

void write_pop(map<vector<int>,karyotype>& m, int g, float p_write, string dir, mt19937& gen){
//void write_pop(list<vector<int>>& cells, int g, int wmax, string fname){
    string gens = to_string(g);
    gens.insert(gens.begin(), 5 - gens.size(), '0');
    std::ofstream cellfile;
    std::string fnamecell=".csv";
    fnamecell = dir+"/" + gens +fnamecell;
    cellfile.open(fnamecell);
    int counter = 0;
    if(m.size()<1){
        cellfile.close();
        return;
    }

    for (auto& [key, value] : m){
        poisson_distribution<> dpois(value.n*p_write);
        int nwrite = dpois(gen);
        if(nwrite>0){
            for(const auto& j : value.cn) cellfile << j << ",";
            cellfile << nwrite << "," << value.fitness << endl;
        }
    }
    cellfile.close();
}

void write_log(float p, float sigma, float mean, float dt, string dir){
    std::ofstream cellfile;
    std::string fname=dir+"/log.txt";
    cellfile.open(fname);
    cellfile << "p,sigma,mean,dt" << endl;
    cellfile << p <<"," << sigma <<"," << mean << ","<< dt <<endl;
    cellfile.close();
}

std::ofstream make_summary_file(std::string dir){
    std::ofstream cellfile;
    std::string fname=dir+"/summary.txt";
    cellfile.open(fname);
    return cellfile;
}

void update_summary_file(int ncells, std::ofstream cellfile){
    cellfile << ncells <<endl;
}

//void write_landscape(fitness_landscape& f, string dir){

  //  std::ofstream cellfile;
    //std::string fname=dir+"/landscape.txt";
    //cellfile.open(fname);

    //int i = 0;
    //for(const auto& p:f.peaks){
      //  for(const auto& cn:p) cellfile << cn << ",";
        //cellfile << f.heights[i] << ",";
        //cellfile << f.sigmas[i] << endl;
        //i++;
    //}

    //cellfile.close();

//}

#endif // WRITE_H_INCLUDED
