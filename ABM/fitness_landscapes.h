struct fitness_landscape{
        int type; // 0 = krig, 1= random, ...
        float deltaf=0; // use to modify fitness for therapy
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

        float get_fitness(vector<int>&);
        float get_grf_fitness(vector<int>&);
        float get_krig_fitness(vector<int>&);
};





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
    for(const auto p:peaks){
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

float fitness_landscape::get_fitness(vector<int>& cn){
    if(type==0) return get_krig_fitness(cn)+deltaf;
    if(type==1) return get_grf_fitness(cn)+deltaf;
    if(type==2) return scale;
    cout << "unknown fitness landscape type" << endl;
    return 0.;
}
