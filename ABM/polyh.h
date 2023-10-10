#include<math.h>
struct polyharmonic_interpolator{
    list<vector<int>> kn;
    int k;
    vector<float> f;
    vector<float> w;
    vector<float> v;

    polyharmonic_interpolator(list<vector<int>>&, vector<float>&, vector<float>&, vector<float>&, int);
    float rbf(float);
    float get_fitness(vector<int>& );
};


// as of right now I don't fancy solving linear equations in cpp so force user to supply w,v.
polyharmonic_interpolator::polyharmonic_interpolator(list<vector<int>>& knots, vector<float>& fvals,
                                                     vector<float>& wvals,vector<float>& vvals, int kk){
    kn=knots;
    k=kk;
    f=fvals;
    w=wvals;
    v=vvals;
}

float polyharmonic_interpolator::rbf(float r){
    float val;
    if(r>=1){
        if(k%2==0){
            val=pow(r,k)*log(r);
        }else{
            val = pow(r,k);
        }
    }else{
        if(k%2==0){
            val=pow(r,k-1)*log(pow(r,r));
        }else{
            val = pow(r,k);
        }
    }
    return val;
}


float polyharmonic_interpolator::get_fitness(vector<int>& pt){
    vector<int> p2;
    p2.push_back(1);
    for(int i=0; i<pt.size(); i++) p2.push_back(pt[i]);

    float fitness = 0;
    vector<float> r; // consider pre-allocating
    for(const auto knot_i:kn){
        float d_i = 0;
        for(int j = 0; j<pt.size(); j++){
            d_i+=pow(knot_i[j]-pt[j],2);
        }
        r.push_back(rbf(sqrt((float)d_i)));
    }
    for(int i = 0; i<r.size(); i++){
        fitness+=w[i]*r[i];
    }
    for(int i = 0; i<p2.size(); i++){
        fitness+=v[i]*p2[i];
    }
    return fitness;
}

