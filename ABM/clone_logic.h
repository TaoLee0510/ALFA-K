#include<stdint.h>
#include <vector>
#include <random>
struct karyotype{
    int n=0;
    int nchrom;
    vector<int> cn;
    float fitness=1;
    karyotype()=default;
    karyotype(vector<int>&,int,float);
    void divide(float,float,mt19937&,list<karyotype>&);
    void divide(float,float,mt19937&,list<karyotype>&,float);
    void divide(float,float,mt19937&,list<karyotype>&,float,float,float);
};

karyotype::karyotype(vector<int> &founder, int N, float f){

    cn=founder;
    fitness = f;
    n=N;
    nchrom = accumulate(cn.begin(),cn.end(),0);

}

// all cells divide
// all cells divide
void karyotype::divide(float p,float pgd, mt19937& gen, list<karyotype>& mutants){

    binomial_distribution<> d1(static_cast<int>(cn.size()), p); // Explicit cast // isn't this an error? first arg should be n..!?

    for(int i=0;i<n;i++){
        int n_mis = d1(gen);
        if(n_mis>0){
            //if misseg happens, decrement cell count by 1
            n--;
            // vectors to store new daughters:
            vector<int> d1 = cn;
            vector<int> d2 = cn;
            // if any chromosome copy numbers go to zero, daughters not valid.
            bool d1_valid=true;
            bool d2_valid=true;
            // following is not strictly accurate since missegregations can be overwritten.
            // may be necessary to fix this.
            for(int j=0;j<n_mis;j++){
                int m = rand()%d1.size(); // choose a chromosome at random
                // following is not accurate since if cn_m1==cn[m], no mis-segregation has occurred.
                int cn_m1 = rand()%(2*cn[m]); // choose a mis-segregated copy number
                int cn_m2 = 2*cn[m]-cn_m1; // second daughter gets the inverse

                // set invalid flags if invalid cn generated:
                if(cn_m1==0) d1_valid = false;
                if(cn_m2==0) d2_valid = false;

                d1[m] = cn_m1;
                d2[m] = cn_m2;
            }
            if(d1_valid){
                    karyotype tmp(d1,1,fitness);
                    mutants.push_back(tmp);
            }
            if(d2_valid){
                    karyotype tmp(d2,1,fitness);
                    mutants.push_back(tmp);
            }
        };
    }
    n*=2; // all the cells that didn't missegregate, divide.
}

//in overloaded version, karyotype fitness is interpreted as net growth rate and
// only fitness*dt*n get to divide. This does have problems though because using the
// net growth rate means we will underestimate the number of divisions.
void karyotype::divide(float p,float pgd, mt19937& gen, list<karyotype>& mutants, float dt){


    float p_div = dt*abs(fitness);
    binomial_distribution<> d0(n, p_div);
    int n_divs = d0(gen);
    if(n_divs==0) return;
    if(fitness<0){
        n-=n_divs;
        return;
    }

    // first we will sort out the number of dividing cells that undergo WGD.
    binomial_distribution<> dgd(n_divs, pgd);
    int ngd = dgd(gen);
    if(ngd>0){
        vector<int> dgd = cn;
        for(int i = 0; i<dgd.size();i++) dgd[i]=2*dgd[i];
        karyotype tmp(dgd,ngd,fitness); // note fitness is recomputed later so the value here doesnt matter
        mutants.push_back(tmp);
    }
    n-=ngd;
    n_divs-=ngd;



    std::uniform_real_distribution<> unidis(0.0, 1.0);

    int nchecked=0;
    int n2;
    int nmis = 0;
    //cout << "ncells: " << n << "; n_divs: " << n_divs << endl;
    while(nchecked<n_divs){
        binomial_distribution<> d1((n_divs-nchecked), pow((1.0-p),nchrom)); // how many cells do not mis-segregate (again)
        int nx = d1(gen);
        if(nmis==0){
            n+=nx; // number of cells that do not mis-segregate at all
        }else{
            for(int k = 0; k<nx; k++){
                n--;
                vector<int> indices; // tells us which chromosome copies have mis-segregated
                // it is assumed the same chromosome copy cannot mis-segregate twice
                for(int i = 0; i< nmis; i++){
                    int r = rand()%nchrom;
                    while(count(indices.begin(), indices.end(), r)){
                        r = rand()%nchrom;
                    }
                    indices.push_back(r);
                }
                sort(indices.begin(),indices.end());
                vector<int> d1 = cn;
                vector<int> d2 = cn;
                bool d1_valid=true;
                bool d2_valid=true;
                int counter = 0;
                int allocated = 0;
                for(int i = 0; i < cn.size(); i++)
                {
                    counter += cn[i];
                    if(allocated < indices.size() && indices[allocated] < counter)
                    {
                        allocated++;
                        if(unidis(gen) < 0.5)
                        {
                            d1[i] += 1;
                            d2[i] -= 1;
                        } else {
                            d2[i] += 1;
                            d1[i] -= 1;
                        }
                        if(d1[i] < 1) d1_valid=false;
                        if(d2[i] < 1) d2_valid=false;
                    }
                }
                if(d1_valid)
                {
                    karyotype tmp(d1,1,fitness);
                    mutants.push_back(tmp);
                }
                if(d2_valid)
                {
                        karyotype tmp(d2,1,fitness);
                        mutants.push_back(tmp);
                }
            }
        }
        nchecked+=nx;
        //cout << "nmis: " << nmis << "; nchecked: " << nchecked << endl;
        nmis++;
    }
}

//in this version, karyotype fitness is interpreted as net growth rate but cell death is still possible for positive growth rates.
// this is required for a model with a carrying capacity to stop the model "freezing".
void karyotype::divide(float p,float pgd, mt19937& gen, list<karyotype>& mutants, float dt, float turnover,float grr){
    //float r = 1/(1+exp(-fitness*k));
    //if(r==0.5) r=0.5001; // following line undefined at r - 0.5
    //float g = fitness*r/(2*r-1);
    //float d = g-fitness;
    //float p_div = dt*g*grr;
    float p_div = (float)n*dt*grr*fitness;
    //binomial_distribution<> d0(n, p_div);
    poisson_distribution<> d0(p_div);
    int n_divs = min(n,d0(gen));
    float p_die = float(n)*dt*turnover;
    //binomial_distribution<> d1(n, p_die);
    poisson_distribution<> d1(p_die);
    int n_deaths = min(n,d1(gen));
    //if(n_divs==0) return;
    //if(fitness<0){
      //  n-=n_divs;
        //return;
    //}

    std::uniform_real_distribution<> unidis(0.0, 1.0);

    int nchecked=0;
    int n2;
    int nmis = 0;
    //cout << "ncells: " << n << "; n_divs: " << n_divs << endl;
    // this while loop is a bit odd because
    while(nchecked<n_divs){
        // approach is to find out how many of the dividing cells don't missegregate.
        binomial_distribution<> d1((n_divs-nchecked), pow((1.0-p),nchrom)); // how many cells do not mis-segregate (again)
        int nx = d1(gen);
        if(nmis==0){
            n+=nx; // number of cells that do not mis-segregate at all
        }else{
            for(int k = 0; k<nx; k++){
                n--;
                vector<int> indices; // tells us which chromosome copies have mis-segregated
                // it is assumed the same chromosome copy cannot mis-segregate twice
                for(int i = 0; i< nmis; i++){
                    int r = rand()%nchrom;
                    while(count(indices.begin(), indices.end(), r)){
                        r = rand()%nchrom;
                    }
                    indices.push_back(r);
                }
                sort(indices.begin(),indices.end());
                vector<int> d1 = cn;
                vector<int> d2 = cn;
                bool d1_valid=true;
                bool d2_valid=true;
                int counter = 0;
                int allocated = 0;
                for(int i = 0; i<cn.size(); i++){
                    counter+=cn[i];
                    if(indices[allocated]<counter){
                        allocated++;
                        if(unidis(gen)<0.5){
                            d1[i]+=1;
                            d2[i]-=1;
                        }else{
                            d2[i]+=1;
                            d1[i]-=1;
                        }
                        if(d1[i]<1) d1_valid=false;
                        if(d2[i]<1) d2_valid=false;
                    }
                }
                if(d1_valid){
                    karyotype tmp(d1,1,fitness);
                    mutants.push_back(tmp);
                }
                if(d2_valid){
                        karyotype tmp(d2,1,fitness);
                        mutants.push_back(tmp);
                }
            }
        }
        nchecked+=nx;
        //cout << "nmis: " << nmis << "; nchecked: " << nchecked << endl;
        nmis++;
    }
    n-=n_deaths;
    n=max(0,n);
}

void rank_selection(std::map<vector<int>,karyotype> &m, mt19937& gen){

        // example:
        // say there are 5 clones in m and their fitness in order of appearance in m is:
        // 4,1,2,9,3
        // we make a vector of pairs with the first item in the pair fitness, second item their order in m:
        // {{4,1}, {1,2}, {2,3}, {9,4}, {3,5}}
        // now sort them:
        // x_sort = {{1,2}, {2,3},{3,5},{4,1},{9,4}}
        // if we iterate along the i'th element of x_sort using:
        // x_rank[x_sort[i]] = i;
        // we get x_rank = {4,1,2,5,3};

        // pair vector stores fitness and the order of each clone in m
        vector<pair<float,int>> fitness_rank;
        // ranks holds the indexes of m ranked in fitness orders
        vector<int> ranks(m.size(),0);

        // populate fitness_rank vector
        int counter = 0;
        for (auto it = m.begin(); it != m.end();){
            fitness_rank.push_back(make_pair((it->second).fitness, counter));
            it++;
            counter++;
        }

        // sort by fitness
        sort(fitness_rank.begin(), fitness_rank.end());

        // populate rank vector
        for(int a = 0; a<fitness_rank.size(); a++){
            ranks[fitness_rank[a].second] = a;
        }

        // kill unfit cells based on rank, since we know
        // m.begin()+counter has rank ranks[counter]
        counter = 0;
        for (auto it = m.begin(); it != m.end();){
            float p_survival = (float)ranks[counter]/(float)m.size();
            binomial_distribution<> d1((it->second).n, p_survival);
            (it->second).n = d1(gen);
            if((it->second).n<1){
                it = m.erase(it);
            }else{
                ++it;
            }
            counter++;
        }

}

void roulette_selection(std::map<vector<int>,karyotype> &m, mt19937& gen, float min_fitness, float max_fitness){
        for (auto it = m.begin(); it != m.end();){
            float p_survival = ((it->second).fitness-min_fitness)/(max_fitness-min_fitness);
            p_survival = min(p_survival,(float)0.9);
            p_survival = max((float)0.1,p_survival);
            binomial_distribution<> d1((it->second).n, p_survival);
            (it->second).n = d1(gen);
            if((it->second).n<1){
                it = m.erase(it);
            }else{
                ++it;
            }
        }

}

void instantiate_population(string filename, map<vector<int>,karyotype>& m, fitness_landscape& fl, mt19937& gen){

   fstream fin;
   fin.open(filename, ios::in);
   std::string tmp, row;
   vector<string> words;
   char delim = ',';
   //int j = 0;
   while(std::getline(fin, row)){
     //   cout<<j<<endl;
       // j++;
        // each row is a karyotype and the final entry is the number. Everything is an int.
        words.clear();
        stringstream tmp2(row);
        while (std::getline(tmp2, tmp, delim)) {
            words.push_back(tmp);
        }
       int rowsize = static_cast<int>(words.size());
       if (words.size() > INT_MAX) {
           throw std::overflow_error("words.size() exceeds int limits");
       }
        vector<int> k;
        int N;
        float f;
        //for(int i = 0; i< (rowsize-1); i++) cout << stoi(words[i]) << ",";
        //cout << endl;
        for(int i = 0; i< (rowsize); i++) k.push_back(stoi(words[i]));
        N = stoi(words[rowsize-1]);
        f = fl.get_fitness(k,gen);

        karyotype kary(k,N,f);
        m[k] = kary;

    }

}
void instantiate_population(vector<int>& k1,int N, map<vector<int>,karyotype>& m, fitness_landscape& fl, mt19937& gen){
    float f0 = fl.get_fitness(k1,gen);
    karyotype c1(k1,N,f0);
    m[k1]=c1;
}

