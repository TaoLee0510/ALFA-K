struct fitness_landscape{
    int type; // 0 = krig, 1= random, ...
    float deltaf=0; // use to modify fitness for therapy
    std::normal_distribution<float> fitnessNoise;
    float noiseVal=0.0; // used to add random noise on a fitness landscape.
    //krig landscape variables:
    vector<vector<float>> knots;
    vector<float> c;
    vector<float> d;
    
    //random field landscape variables
    list<vector<int>> peaks;
    float wavelength;
    float scale;
    float centre;
    fitness_landscape(){};
    // krig constructor
    void init(vector<vector<float>>& k, vector<float>& cc,vector<float>& dd, float delta_f){
        type=0;
        deltaf=delta_f;
        knots=k;
        c=cc;
        d=dd;
    };
    void init(list<vector<int>>& p, float w, float s, float c ,float delta_f){
        type=1;
        deltaf=delta_f;
        peaks = p;
        wavelength = w;
        scale = s;
        centre = c;
    };// random field constructor
    void init(float fitness){
        type=2;
        scale = fitness;
    };// flat landscape constructor
    void init(vector<vector<float>>& k, vector<float>& cc, float delta_f){
        type=3;
        deltaf=delta_f;
        knots=k;
        c=cc;
    };
    
    void setNoise(float noise){
        std::normal_distribution<float> d{0, noise};
        fitnessNoise = d;
    }
    
    float get_fitness(vector<int>&, mt19937&);
    float get_grf_fitness(vector<int>&);
    float get_krig_fitness(vector<int>&);
    float get_Predicted_fitness(vector<int>&);
};

float fitness_landscape::get_Predicted_fitness(vector<int>& cn){
    float f=0;
    vector<float> xx0;
    for(const auto& i:cn) xx0.push_back((float)i);
    int chrom_number = int(xx0.size());
    for(int i = 0; i<knots.size(); i++)
    {
        int Di = 0;
        for(int j = 0; j<knots[i].size(); j++)
        {
            if(knots[i][j] == xx0[j] )
            {
                Di = Di+1;
            }
        }
        if(Di == chrom_number)
        {
            f = c[i];
            break;
        }
    }
    return(f);
}



float fitness_landscape::get_krig_fitness(vector<int>& cn){
    float f=0, xx1=0;
    vector<float> xx0;// = {1.0};
    for(const auto& i:cn) xx0.push_back((float)i);
    
    for(int i = 0; i<knots.size(); i++){
        float Di = 0;
        for(int j = 0; j<knots[i].size(); j++){
            Di+=pow((knots[i][j]-xx0[j]),2);
        }
        f+=c[i]*exp(-sqrt(Di));
    }
    f+=d[0];//xx1;
    return(f);
}

float fitness_landscape::get_grf_fitness(vector<int>&cn){
    float fitness = 0;
    for (const auto& p : peaks) {  // Use a reference to avoid copying
        float d=0;
        for(int i =0; i<cn.size(); i++){
            d+=pow((float)cn[i]-p[i],2);
        }
        d=sqrt(d);
        
        fitness+=sin(d/wavelength);
    }
    fitness*=scale;
    fitness+=centre;
    return fitness;
}

float fitness_landscape::get_fitness(vector<int>& cn, mt19937& gen){
    float noise = fitnessNoise(gen);
    if(type==0) return get_krig_fitness(cn)+deltaf+noise;
    if(type==1) return get_grf_fitness(cn)+deltaf+noise;
    if(type==2) return scale+noise;
    if(type==3) return get_Predicted_fitness(cn)+deltaf;
    cout << "unknown fitness landscape type" << endl;
    return 0.;
}
