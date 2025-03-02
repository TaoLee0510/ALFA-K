#include<stdint.h>
struct parameters{

    float p=0.001;
    float pgd=0.;
    float deltaf=0.0;
    int pop_write_freq = 1000000000;
    float deltafu=0.0;
    float wavelength=1.;
    float centre = 0.;
    float scale = 1.;
    float turnover_rate = 0.0;
    float sd_mutation,mean_mutation;
    float dt = 0.1;
    float fitness_noise=0.0;
    int Nsteps = 1000;
    int target_output_size=1000;
    float t_on = 9.;
    float t_off=3.;
    int n0 = 100000;
    int ncycles = 0;
    int treat_on = 500000;
    int treat_off=250000;
    string turnover_enabled = "no";
    string fitness_landscape_type;
    string fitness_landscape_file="not_supplied";
    string fitness_landscape_file2="not_supplied";
    string output_dir = "output";
    string population_file = "not_supplied";
    int output_gens=1000;
    int init_size=10000;
    int max_size=2000000;
    float bf=0.01;//,0.05,0.1,0.5};
    int maxchrom=8; // copy number max limit
    int G = 500;
    float p_death=0.0;
    vector<int> init_kary;

    list<vector<int>> peaks;
    vector<float> heights;
    vector<float> sigma;

    //for the polyh landscape
    vector<float> f;
    vector<float> w;
    vector<float> v;

    //for the krig landscape
    vector<vector<float>> knots;
    vector<float> c;
    vector<float> d;

    parameters(string path);
    void read_landscape_file();
    void read_landscape_file(string);
    void read_krig_file();
    void read_gaussian_file();
    void read_polyh_file();
    void read_random_file();
    void read_random_file(string);
    void read_krig_file(string);
    void read_Predicted_file();
    void read_Predicted_file(string);

};

parameters::parameters(string path){

   fstream fin;
   fin.open(path, ios::in);
   std::string tmp, row;
   vector<string> words;
   char delim = ',';

   while(std::getline(fin, row)){
        words.clear();
        stringstream tmp2(row);
        while (std::getline(tmp2, tmp, delim)) {
        words.push_back(tmp);

        }

        if(words[0]=="init_size") init_size=stoi(words[1]);
        if(words[0]=="ncycles") ncycles=stoi(words[1]);
        if(words[0]=="fitness_landscape_type") fitness_landscape_type=words[1];
        if(words[0]=="turnover_enabled") turnover_enabled=words[1];
        if(words[0]=="population_file") population_file=words[1];
        if(words[0]=="fitness_landscape_file") {
            // remove leading whitespace from string filepath:
            string s = words[1];
            s.erase(std::remove_if(s.begin(), s.end(), ::isspace),s.end());
            fitness_landscape_file=s;
        }
        if(words[0]=="fitness_landscape_file2") {
            // remove leading whitespace from string filepath:
            string s = words[1];
            s.erase(std::remove_if(s.begin(), s.end(), ::isspace),s.end());
            fitness_landscape_file2=s;
        }
        if(words[0]=="deltaf") deltaf=stof(words[1]);
        if(words[0]=="deltafu") deltafu=stof(words[1]);
        if(words[0]=="output_dir") {
            // remove leading whitespace from string filepath:
            string s = words[1];
            s.erase(std::remove_if(s.begin(), s.end(), ::isspace),s.end());
            output_dir=s;
        }
        if(words[0]=="turnover_rate") turnover_rate=stof(words[1]);
        if(words[0]=="t_on") t_on=stof(words[1]);
        if(words[0]=="t_off") t_off=stof(words[1]);
        if(words[0]=="n0") n0=stoi(words[1]);
        if(words[0]=="pop_write_freq") pop_write_freq=stoi(words[1]);
        if(words[0]=="dt") dt=stof(words[1]);
        if(words[0]=="p") p=stof(words[1]);
        if(words[0]=="pgd") pgd=stof(words[1]);
        if(words[0]=="treat_on") treat_on=stoi(words[1]);
        if(words[0]=="treat_off") treat_off=stoi(words[1]);
        if(words[0]=="wavelength") wavelength=stof(words[1]);
        if(words[0]=="scale") scale=stof(words[1]);
        if(words[0]=="centre") centre=stof(words[1]);
        if(words[0]=="fitness_noise") fitness_noise=stof(words[1]);
        if(words[0]=="Nsteps") Nsteps=stoi(words[1]);
        if(words[0]=="init_size") init_size=stof(words[1]);
        if(words[0]=="bf") bf=stof(words[1]);
        if(words[0]=="max_size") max_size=stoi(words[1]);
        if(words[0]=="init_kary"){
            for(int i=1; i<words.size(); i++){
                init_kary.push_back(stoi(words[i]));
            }
        }
    }

}

void parameters::read_landscape_file(){
    if(fitness_landscape_type=="gaussian") read_gaussian_file();
    if(fitness_landscape_type=="polyh") read_polyh_file();
    if(fitness_landscape_type=="random") read_random_file();
    if(fitness_landscape_type=="krig") read_krig_file();
    if(fitness_landscape_type=="Predicted") read_Predicted_file();
}

void parameters::read_landscape_file(string flf){
    if(fitness_landscape_type=="gaussian") read_gaussian_file();
    if(fitness_landscape_type=="polyh") read_polyh_file();
    if(fitness_landscape_type=="random") read_random_file(flf);
    if(fitness_landscape_type=="krig") read_krig_file(flf);
    if(fitness_landscape_type=="Predicted") read_Predicted_file(flf);
}

// all these read functions are better associated with their respective fitness landscapes
void parameters::read_gaussian_file(){

   fstream fin;
   fin.open(fitness_landscape_file, ios::in);
   std::string tmp, row;
   vector<string> words;
   char delim = ',';
   while(std::getline(fin, row)){

        words.clear();
        stringstream tmp2(row);
        while (std::getline(tmp2, tmp, delim)) {
            words.push_back(tmp);
        }
       int rowsize = static_cast<int>(words.size());
       if (words.size() > INT_MAX) {
           throw std::overflow_error("words.size() exceeds int limits");
       }
        vector<int> peak;
        for(int i = 0; i< (rowsize-2); i++){
            peak.push_back(stoi(words[i]));
        }
        peaks.push_back(peak);
        heights.push_back(stof(words[rowsize-2]));
        sigma.push_back(stof(words[rowsize-1]));
    }

}

void parameters::read_random_file(){
   peaks.clear();
   fstream fin;
   fin.open(fitness_landscape_file, ios::in);
   std::string tmp, row;
   vector<string> words;
   char delim = ',';
   while(std::getline(fin, row)){

        words.clear();
        stringstream tmp2(row);
        while (std::getline(tmp2, tmp, delim)) {
            words.push_back(tmp);
        }
       int rowsize = static_cast<int>(words.size());
       if (words.size() > INT_MAX) {
           throw std::overflow_error("words.size() exceeds int limits");
       }
        vector<int> peak;
        for(int i = 0; i< rowsize; i++){
            peak.push_back(stoi(words[i]));
        }
        peaks.push_back(peak);
    }

}

void parameters::read_krig_file(){
    knots.clear();
    c.clear();
    d.clear();
   fstream fin;
   fin.open(fitness_landscape_file, ios::in);
   std::string tmp, row;
   vector<string> words;
   char delim = ',';
   while(std::getline(fin, row)){
        // the last row of the input is v. Approach is to clear v and reset it every time through the loop,
        // so at the end of the loop v is correct
        d.clear();
        words.clear();
        stringstream tmp2(row);
        while (std::getline(tmp2, tmp, delim)) {
            words.push_back(tmp);
        }
       int rowsize = static_cast<int>(words.size());
       if (words.size() > INT_MAX) {
           throw std::overflow_error("words.size() exceeds int limits");
       }
        vector<float> knot;
        for(int i = 0; i< (rowsize-1); i++){
            knot.push_back(stof(words[i]));
            d.push_back(stof(words[i]));
        }
        knots.push_back(knot);
        c.push_back(stof(words[rowsize-1]));
        d.push_back(stof(words[rowsize-1]));
    }

   // remove the last elements (which are v)
   c.pop_back();
   knots.pop_back();
   //for(const auto ci:c) cout << ci << ",";
   //cout << endl;

   //for(const auto di:d) cout << di << ",";
   //cout << endl;

}

void parameters::read_Predicted_file(){
    knots.clear();
    c.clear();
   fstream fin;
   fin.open(fitness_landscape_file, ios::in);
   std::string tmp, row;
   vector<string> words;
   char delim = ',';
   while(std::getline(fin, row)){
        // the last row of the input is v. Approach is to clear v and reset it every time through the loop,
        // so at the end of the loop v is correct
        d.clear();
        words.clear();
        stringstream tmp2(row);
        while (std::getline(tmp2, tmp, delim)) {
            words.push_back(tmp);
        }
       int rowsize = static_cast<int>(words.size());
       if (words.size() > INT_MAX) {
           throw std::overflow_error("words.size() exceeds int limits");
       }
        vector<float> knot;
        for(int i = 0; i< (rowsize-1); i++){
            knot.push_back(stof(words[i]));
        }
        knots.push_back(knot);
        c.push_back(stof(words[rowsize-1]));
    }

   // remove the last elements (which are v)
   c.pop_back();
   knots.pop_back();

}

void parameters::read_random_file(string flf){
   peaks.clear();
   fstream fin;
   fin.open(flf, ios::in);
   std::string tmp, row;
   vector<string> words;
   char delim = ',';
   while(std::getline(fin, row)){

        words.clear();
        stringstream tmp2(row);
        while (std::getline(tmp2, tmp, delim)) {
            words.push_back(tmp);
        }
       int rowsize = static_cast<int>(words.size());
       if (words.size() > INT_MAX) {
           throw std::overflow_error("words.size() exceeds int limits");
       }
        vector<int> peak;
        for(int i = 0; i< rowsize; i++){
            peak.push_back(stoi(words[i]));
        }
        peaks.push_back(peak);
    }

}

void parameters::read_krig_file(string flf){
    knots.clear();
    c.clear();
    d.clear();
   fstream fin;
   fin.open(flf, ios::in);
   std::string tmp, row;
   vector<string> words;
   char delim = ',';
   while(std::getline(fin, row)){
        // the last row of the input is v. Approach is to clear v and reset it every time through the loop,
        // so at the end of the loop v is correct
        d.clear();
        words.clear();
        stringstream tmp2(row);
        while (std::getline(tmp2, tmp, delim)) {
            words.push_back(tmp);
        }
       int rowsize = static_cast<int>(words.size());
       if (words.size() > INT_MAX) {
           throw std::overflow_error("words.size() exceeds int limits");
       }
        vector<float> knot;
        for(int i = 0; i< (rowsize-1); i++){
            knot.push_back(stof(words[i]));
            d.push_back(stof(words[i]));
        }
        knots.push_back(knot);
        c.push_back(stof(words[rowsize-1]));
        d.push_back(stof(words[rowsize-1]));
    }

   // remove the last elements (which are v)
   c.pop_back();
   knots.pop_back();
   //for(const auto ci:c) cout << ci << ",";
   //cout << endl;

//   for(const auto di:d) cout << di << ",";
  // cout << endl;

}

void parameters::read_Predicted_file(string flf){
    knots.clear();
    c.clear();
   fstream fin;
   fin.open(flf, ios::in);
   std::string tmp, row;
   vector<string> words;
   char delim = ',';
   while(std::getline(fin, row)){
        // the last row of the input is v. Approach is to clear v and reset it every time through the loop,
        // so at the end of the loop v is correct
        d.clear();
        words.clear();
        stringstream tmp2(row);
        while (std::getline(tmp2, tmp, delim)) {
            words.push_back(tmp);
        }
       int rowsize = static_cast<int>(words.size());
       if (words.size() > INT_MAX) {
           throw std::overflow_error("words.size() exceeds int limits");
       }
        vector<float> knot;
        for(int i = 0; i< (rowsize-1); i++){
            knot.push_back(stof(words[i]));
        }
        knots.push_back(knot);
        c.push_back(stof(words[rowsize-1]));
    }

   // remove the last elements (which are v)
   c.pop_back();
   knots.pop_back();

}


void parameters::read_polyh_file(){

   fstream fin;
   fin.open(fitness_landscape_file, ios::in);
   std::string tmp, row;
   vector<string> words;
   char delim = ',';
   while(std::getline(fin, row)){
        // the last row of the input is v. Approach is to clear v and reset it every time through the loop,
        // so at the end of the loop v is correct
        v.clear();
        words.clear();
        stringstream tmp2(row);
        while (std::getline(tmp2, tmp, delim)) {
            words.push_back(tmp);
        }
       int rowsize = static_cast<int>(words.size());
       if (words.size() > INT_MAX) {
           throw std::overflow_error("words.size() exceeds int limits");
       }
        vector<int> peak;
        for(int i = 0; i< (rowsize-1); i++){
            peak.push_back(stoi(words[i]));
            v.push_back(stof(words[i]));
        }
        peaks.push_back(peak);
        w.push_back(stof(words[rowsize-1]));
        v.push_back(stof(words[rowsize-1]));
    }

   // remove the last elements (which are v)
   w.pop_back();
   peaks.pop_back();

}



