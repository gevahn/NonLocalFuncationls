#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <cmath>

#include "armadillo"

using namespace arma;
using namespace std;

float p = 2.0;

double lda_xc_hole (double k_f, double distance);
double xcp_hole(double k_fdis);
double sym_xcp_hole(vec k_f,double dis,int i, int j);
double gj_xc_hole(double rho, double distance);
double gaussian_xc_hole(double rho, double distance);

vec k_fermi (vec rho);

double wdaintegral(vector<double> rho_vec, vector<double> k_vec, vector<double> x_vec, vector<double> y_vec, vector<double> z_vec, vector<double> w_vec, double(*xc_hole)(double,double));
double f_wdaintegral(vec f,vector<double> rho_vec, vector<double> k_vec, vector<double> x_vec, vector<double> y_vec, vector<double> z_vec, vector<double> w_vec);
//double normintegral(vector<double> rho_vec, double k_vec, vector<double> x_vec, vector<double> y_vec, vector<double> z_vec, vector<double> w_vec, int j);
double symnormintegral(vector<double> rho_vec, vec k_vec, vector<double> x_vec, vector<double> y_vec, vector<double> z_vec, vector<double> w_vec, int j);
//double f_integral(vec f, vector<double> rho_vec, vector<double> k_vec, vector<double> x_vec, vector<double> y_vec, vector<double> z_vec, vector<double> w_vec, int j);
//double symnormintegral_d(vector<double> rho_vec, vec k_vec, vector<double> x_vec, vector<double> y_vec, vector<double> z_vec, vector<double> w_vec, int j);
double g_symnormintegral(vector<double> rho_vec, vec k_vec, vector<double> x_vec, vector<double> y_vec, vector<double> z_vec, vector<double> w_vec, int j, double g, double (*xc_hole)(double ,double));

double gj_normintegral(vector<double> rho_vec, double k_vec, vector<double> x_vec, vector<double> y_vec, vector<double> z_vec, vector<double> w_vec, int j);
double gj_wdaintegral(vector<double> rho_vec, vector<double> k_vec, vector<double> x_vec, vector<double> y_vec, vector<double> z_vec, vector<double> w_vec);
double gj_g_symnormintegral(vector<double> rho_vec, vec k_vec, vector<double> x_vec, vector<double> y_vec, vector<double> z_vec, vector<double> w_vec, int j, double g);

double brents(double x1, double x2, double tol,vector<double> rho_vec, vector<double> x_vec, vector<double> y_vec, vector<double> z_vec, vector<double> w_vec, int j, double(*xc_hole)(double,double));
double nwt(double x1, double x2, double tol,vector<double> rho_vec, vector<double> x_vec, vector<double> y_vec, vector<double> z_vec, vector<double> w_vec, int j);
//vec f_gcr(vec f, vec b, vector<double> rho_vec, vector<double> k_vec, vector<double> x_vec, vector<double> y_vec, vector<double> z_vec, vector<double> w_vec);
vec gcr(vec b, vector<double> rho_vec, vec k_vec, vector<double> x_vec, vector<double> y_vec, vector<double> z_vec, vector<double> w_vec, double g, double tol, double maxitr,double(*xc_hole)(double,double));

double npow(double b, double e);
int main(int argc, char** argv)
  {
    wall_clock timer;
    double n_secs;
    vector<double> rho_vec;
    vector<double> rho_avec;
    vector<double> rho_bvec;
    vector<double> rho_avec_two;
    vector<double> rho_bvec_two;
    vector<double> w_vec;
    vector<double> x_vec;
    vector<double> y_vec;
    vector<double> z_vec;
    vector<double> k_vec;
    vector<double> ka_vec;
    vector<double> kb_vec;

    cout.precision(16);

    int hole_idx = 0;



    string filename;
    if (argc == 1){
        cout << "Insert density file name: ";

//        cin >> filename;
        filename = "n2194.dat";
    }
    else {
        char *k;
        filename = argv[1];
        p = strtod(argv[2],&k);
        hole_idx = strtol(argv[3],&k,10);

    }

     // XC HOLE FUNCTION TO USE

    vector<double (*)(double,double)> xc_hole_vec;
    vector<string> hole_names;

    xc_hole_vec.push_back(&lda_xc_hole);
    hole_names.push_back("LDA");
    xc_hole_vec.push_back(&gj_xc_hole);
    hole_names.push_back("GJ");
    xc_hole_vec.push_back(&gaussian_xc_hole);
    hole_names.push_back("Gaussian");
    



    mat raw_data;
    raw_data.load(filename);
    cout << filename << endl;
    cout << hole_names[hole_idx] << endl;
    cout << "p value: " << p << endl;
    cout << raw_data.n_cols << " " << raw_data.n_rows << endl;


    vec rho_a;
    vec rho_b;
    vec weight;
    mat xyz_coord;
    xyz_coord = raw_data.cols(0, 2);
    rho_a = raw_data.col(4);
    rho_b = raw_data.col(5);
    weight = raw_data.col(3);


    vec k_f;
    vec k_f_a;
    vec k_f_b;

    // chaged K_f to be a parameter NOT FERMI WAVE VECTOR
    k_f = k_fermi(rho_a + rho_b);
    k_f_a = k_fermi(2*rho_a);
    k_f_b = k_fermi(2*rho_b);

    k_f_a.save("p0.out",raw_ascii);

    double nElec;
    double ldaX;

    nElec = dot(rho_a + rho_b, weight);
    ldaX = 0.5*(dot(pow(2 * rho_a, (4 / 3.0)), weight) + dot(pow(2 * rho_b, (4 / 3.0)), weight));
    ldaX = (pow(3 / datum::pi, (1.0 / 3))) * ldaX * (-0.75);

    cout << "number of elctrons: " << nElec << endl;
    cout << "lda: " << ldaX << endl;


    for (int i = 0; i < raw_data.n_rows; i++){
        //if (weight(i) > datum::eps * pow(10.0,0.0) && rho_a(i) > datum::eps * pow(10.0,0.0)) {
            rho_vec.push_back(rho_a(i) + rho_b(i));
            rho_avec.push_back(rho_a(i));
            rho_bvec.push_back(rho_b(i));
            rho_avec_two.push_back(2*rho_a(i));
            rho_bvec_two.push_back(2*rho_b(i));
            w_vec.push_back(weight(i));
            x_vec.push_back(xyz_coord(i, 0));
            y_vec.push_back(xyz_coord(i, 1));
            z_vec.push_back(xyz_coord(i, 2));
            k_vec.push_back(k_f(i));
            ka_vec.push_back(k_f_a(i));
            kb_vec.push_back(k_f_b(i));
        //}
    }
    cout << ka_vec.size() <<endl;
    rho_a.clear();
    rho_b.clear();

    rho_a = rho_avec;
    rho_b = rho_bvec;

    k_f_a.clear();
    k_f_b.clear();
    k_f.clear();

    k_f = k_fermi(rho_a + rho_b);
    k_f_a = k_fermi(2*rho_a);
    k_f_b = k_fermi(2*rho_b);



    int n = k_f_a.n_elem;

    timer.tic();


    float gj_wpa0 = wdaintegral(rho_avec, ka_vec, x_vec, y_vec, z_vec, w_vec,xc_hole_vec[hole_idx]);
    float wpa0 = wdaintegral(rho_avec, ka_vec, x_vec, y_vec, z_vec, w_vec,lda_xc_hole);
    cout << "WPA0 exchange: " << wpa0 << endl;

    double error_sum = 0;
    int error_sum_num = 0;
    double lower_bracket = 0.00001;
    double upper_bracket = 20;

//    for (int i = 0; i < 10 ; i++){
//        cout << g_symnormintegral(rho_avec, i*ones(ka_vec.size()) , x_vec, y_vec, z_vec, w_vec, i, 0, xc_hole) << endl;
//    }

    for (int i = 0; i < n ; i++){
        if (g_symnormintegral(rho_avec, lower_bracket*ones(ka_vec.size()) , x_vec, y_vec, z_vec, w_vec, i, 0, xc_hole_vec[hole_idx]) * g_symnormintegral(rho_avec, upper_bracket*ones(ka_vec.size()) , x_vec, y_vec, z_vec, w_vec,i,0,xc_hole_vec[hole_idx]) > 0){
               error_sum_num++;
        }
    }


    cout << "non bracketed roots: " << error_sum_num << endl;
    error_sum_num = 0;

//ka_vec[0] = brents(0.0001,1,0.00001,rho_avec, x_vec, y_vec, z_vec, w_vec, 2700 );
    for (int i = 0; i < n ; i++){
        //ka_vec[i] = nwt(ka_vec[i],100,0.0001,rho_avec, x_vec, y_vec, z_vec, w_vec, i );
        ka_vec[i] = brents(0,15,0.00001,rho_avec, x_vec, y_vec, z_vec, w_vec, i, xc_hole_vec[hole_idx] );
        if (ka_vec[i] == 1){

            ka_vec[i] = 0;
        }
        else
            if(abs(g_symnormintegral(rho_avec, ka_vec, x_vec, y_vec, z_vec, w_vec, i,0,xc_hole_vec[hole_idx])) > 0.0001)
                error_sum_num++;
        error_sum += abs(g_symnormintegral(rho_avec, ka_vec, x_vec, y_vec, z_vec, w_vec, i,0,xc_hole_vec[hole_idx]));
        kb_vec[i] = ka_vec[i];
    //cout << ka_vec[i] <<endl;



    }
    cout << "number of sum rule error points: " << error_sum_num << endl;
    cout << "sum rule error: " << error_sum << endl;

    k_f_a = ka_vec;
    k_f_a.save("p0.out",raw_ascii);


//    ka_vec[0] = nwt(ka_vec[0],15,0.0001,rho_avec, x_vec, y_vec, z_vec, w_vec, 0 );
//    for (int i = 1; i < n ; i++){
//        if (ka_vec[i] < 10 * ka_vec[i-1])
//            ka_vec[i] = nwt(ka_vec[i-1],15,0.0001,rho_avec, x_vec, y_vec, z_vec, w_vec, i );
//        else
//            ka_vec[i] = nwt(ka_vec[i],15,0.0001,rho_avec, x_vec, y_vec, z_vec, w_vec, i );
//        kb_vec[i] = ka_vec[i];
//    }


    float gj_wpa1 = wdaintegral(rho_avec, ka_vec, x_vec, y_vec, z_vec, w_vec,xc_hole_vec[hole_idx]);
    float wpa1 = wdaintegral(rho_avec, ka_vec, x_vec, y_vec, z_vec, w_vec,lda_xc_hole);
    cout << "GJ wda 0p: " << gj_wpa0 << endl;
    cout << "LDA wda 0p: " << wpa0 << endl;
    cout << "GJ wda 1p: " << gj_wpa1 << endl;
    cout << "LDA wda 1p: " << wpa1 << endl;




//    n = 500;
    int itr = 0;
    int p = 20;
    int m = 0;
    vector <double> ne;
    double ns;
    double tmp;
    double delx = 0.000001;
    double tau;
    double trust = 1;
    double max_trust = 5;
    double eta = 0.1;
    double lam;
    double alpha;
    mat C = zeros<mat>(n,p);
    mat D = zeros<mat>(n,p);
    mat R = zeros<mat>(p,p);
//    mat S = zeros<mat>(p,p);
    vec S;
    mat W = zeros<mat>(p,p);
    mat tmp_mat = zeros<mat>(p,p);
    mat tmp_vec;

    vec g;
    vec s;
    vec s_c;
    vec y;
    vec x;
    vec f;
    vec grad;
    vec b;

    mat B2 = zeros<mat>(p);

    double rho;
//    int q = 2;
    g.zeros(n);
    grad.zeros(n);
    S.ones(n);
    y.zeros(n);
    x.zeros(n);
    f.ones(n);
    b.zeros(n);
//    tmp_vec.zeros(n);

    tmp_vec.zeros(n);
    for (int i=0;i < n;i++){
        x(i) = ka_vec[i];
    }
    x.save("p1.out",raw_ascii);
//    for (int i=0;i < n;i++){
//        b(i)=-1 * f_integral(f, rho_avec, ka_vec, x_vec, y_vec, z_vec, w_vec, i);
//    }



//    while ((abs(ne[itr]) > 0.01) && (itr < 100)){
//            f = gcr(f, b, rho_avec, ka_vec, x_vec, y_vec, z_vec, w_vec);

//            for (int i=0;i < n;i++){
//                b(i)=-1 * f_integral(f, rho_avec, ka_vec, x_vec, y_vec, z_vec, w_vec, i);
//            }
//            ne.push_back(norm(b));
//            itr++;
//            cout << "itr: " << itr << " error: " << ne[itr] << endl;
//    }

    lam = 0;
    double wpa2 = wdaintegral(rho_avec, ka_vec, x_vec, y_vec, z_vec, w_vec,xc_hole_vec[hole_idx]);
    double dellam = 0.1;
    double tol = 0.0001;
    double maxitr = 20;
    while (lam < 1){

        itr = 0;

        for (int i=0;i < n;i++){
            b(i)= g_symnormintegral(rho_avec, x, x_vec, y_vec, z_vec, w_vec, i,1,xc_hole_vec[hole_idx]);
        }
        cout << "new lambda: " << lam << " lam = 1 error: " << norm(b) << endl;

        cout << "new lambda: " << lam << " wpa2: " << wpa2 << endl;

        for (int i=0;i < n;i++){
            b(i)=-1 * g_symnormintegral(rho_avec, x, x_vec, y_vec, z_vec, w_vec, i,lam,xc_hole_vec[hole_idx]);
        }
        ne.clear();
        ne.push_back(norm(b));
        cout << "new lambda: " << lam << " error: " << ne[itr] << endl;
        S = x;


        while ((abs(ne[itr]) > tol) && (itr < maxitr)){
                S = gcr(b, rho_avec, S, x_vec, y_vec, z_vec, w_vec, (lam), 0.0001*ne[itr] ,400,xc_hole_vec[hole_idx]);

                for (int i=0;i < n;i++){
                    b(i)=-1 * g_symnormintegral(rho_avec, S, x_vec, y_vec, z_vec, w_vec, i, lam, xc_hole_vec[hole_idx]);
                }
                ne.push_back(norm(b));
                itr++;
                cout <<"lambda: " << lam << " itr: " << itr << " error: " << ne[itr] << endl;
        }
        if (abs(ne[itr]) > tol){
            lam = lam - dellam;
            dellam = dellam / 2;
            maxitr = maxitr + 5;
//            tol  = tol / 2;
        }
        else
        {
//            tol = tol * 1.2;
            x = S;
        }
        lam = lam + dellam;
        for (int i = 0;i < x.n_elem;i++){
            ka_vec[i] = x(i);
        }
        wpa2 = wdaintegral(rho_avec, ka_vec, x_vec, y_vec, z_vec, w_vec,xc_hole_vec[hole_idx]);
    }
    lam = 1;
    for (int i=0;i < n;i++){
        b(i)=-1 * g_symnormintegral(rho_avec, x, x_vec, y_vec, z_vec, w_vec, i,lam,xc_hole_vec[hole_idx]);
    }

    tol = 0.01;
    maxitr = 20;
    itr = 0;
    cout <<"lambda: " << lam << " itr: " << itr << " error: " << ne[itr] << endl;
    while ((abs(ne[itr]) > tol) && (itr < maxitr)){
            S = gcr(b, rho_avec, S, x_vec, y_vec, z_vec, w_vec, (lam), 0.01, 20,xc_hole_vec[hole_idx]);

            for (int i=0;i < n;i++){
                b(i)=-1 * g_symnormintegral(rho_avec, S, x_vec, y_vec, z_vec, w_vec, i, lam,xc_hole_vec[hole_idx]);
            }
            ne.push_back(norm(b*rho_a.diag()));
            itr++;
//            cout <<"lambda: " << lam << " itr: " << itr << " error: " << ne[itr] << " vec error:" << norm(b*rho_a.diag()) << endl;
            cout <<"lambda: " << lam << " itr: " << itr << " error: " << ne[itr] << endl;
            x = S;
            for (int i = 0;i < x.n_elem;i++){
                ka_vec[i] = x(i);
            }
            wpa2 = wdaintegral(rho_avec, ka_vec, x_vec, y_vec, z_vec, w_vec,xc_hole_vec[hole_idx]);
            cout <<"lambda: " << lam << " wpa2: " << wpa2 << endl;
    }

    b.save("error.out",raw_ascii);    
////    x.ones(n);
////    mat A;
////    A.randn(n,n);


//    for (int i=0;i < n;i++){
//        g(i)=symnormintegral(rho_avec, x, x_vec, y_vec, z_vec, w_vec, i);
//    }
//    y = g;
////    g = A*x-S;
//    ne.push_back(norm(g));
//    cout << "itr: " << itr << " error: " << ne[itr] <<endl;
//    mat B = eye(n,n);
//    for (int i = 0;i < n; i++){
//        tmp_vec = x;
//        tmp_vec(i) += delx;
//        tmp = symnormintegral(rho_avec, tmp_vec, x_vec, y_vec, z_vec, w_vec, i);
//        B(i,i) = (tmp - g(i)) / delx;
//    }

//    while ((abs(ne[itr]) > 0.01) && (itr < 100)){
//        s = -B*g;
//        ns = norm(s);

//        alpha = 1;
//        while (alpha >0 && norm(y) >= norm(g)){
//            alpha = alpha - 0.1;
//            for (int i=0;i < n;i++){
//                y(i) = symnormintegral(rho_avec, x + alpha * s, x_vec, y_vec, z_vec, w_vec, i);
//            }
//            cout << "alpha: " <<alpha << " norm y:" <<norm(y) <<endl;


//        }

//        x = x + alpha * s;
//        for (int i=0;i < n;i++){
//            y(i) = symnormintegral(rho_avec, x, x_vec, y_vec, z_vec, w_vec, i);
//        }

//        y = y - g;
//        g = y + g;
//        itr++;
//        ne[itr] = norm(g);
//        B = B + (y - B * s)*s.t()/(ns*ns);

//        y = g - y;
//        for (int i = 0;i < n; i++){
//            tmp_vec = f;
//            tmp_vec(i) += delx;
//            tmp = symnormintegral(rho_avec, tmp_vec, x_vec, y_vec, z_vec, w_vec, i);
//            B(i,i) = (tmp - g(i)) / delx;
//        }

//        if (cond(B2) > 10000000000000){
//            cout << "B2 singular at " << itr << endl;
//            break;
//        }

//        cout << "itr: " << itr << " error: " << ne[itr] <<" svec: " << norm(s) <<endl;

//    }
//    double f_wpa2 = 0.5*(f_wdaintegral(f,rho_avec, ka_vec, x_vec, y_vec, z_vec, w_vec) + wdaintegral(rho_bvec, kb_vec, x_vec, y_vec, z_vec, w_vec));

//    itr = 0;
//    ne.clear();
//    f.ones(n);
//    for (int i=0;i < n;i++){
//        g(i)=f_integral(f, rho_avec, ka_vec, x_vec, y_vec, z_vec, w_vec, i);
//    }
//    for (int i = 0;i < n; i++){
//        tmp_vec = f;
//        tmp_vec(i) += delx;
//        tmp = f_integral(tmp_vec, rho_avec, ka_vec, x_vec, y_vec, z_vec, w_vec, i);
//        grad(i) = (tmp - g(i)) / delx;
//    }

//    ne.push_back(sqrt(dot(g,g)));
//    while ((ne[itr] > 0.01) && (itr < 100)){
//        cout << "itr: " << itr << " error: " << ne[itr] <<" xvec: " << norm(f) <<endl;
//        B2 = eye(p,p)-D.t()*C;
////        cout << B2 << " ";
//        if (cond(B2) > 10000000000000){
//            cout << "B2 singular at " << itr << endl;
//            break;
//        }
//        s = C * ( solve(B2 , D.t()*g )) + g;

//        ns = sqrt(dot(s,s));
//        if (ns <= 0){
//            cout <<"s is a zero vector" << endl;
//            break;
//        }
//        f = f + s;

//        for (int i=0;i < n;i++){
//            y(i) = f_integral(f, rho_avec, ka_vec, x_vec, y_vec, z_vec, w_vec, i);
//        }
//        for (int i = 0;i < n; i++){
//            tmp_vec = f;
//            tmp_vec(i) += delx;
//            tmp = f_integral(tmp_vec, rho_avec, ka_vec, x_vec, y_vec, z_vec, w_vec, i);
//            grad(i) = (tmp - y(i)) / delx;
//        }
//        y = y - g;
//        g = y + g;
//        itr++;

//        ne.push_back(sqrt(dot(g,g)));

//        if (m==p){
//            qr_econ(D,R,D);
//            C = C * R.t();
//            svd_econ(C,S,W,C);
//            C = C * diagmat(S);
//            D = D * W;
//            C.cols(0,p-2) = C.cols(1,p-1);
//            C.col(p-1) = zeros(n,1);
//            D.cols(0,p-2) = D.cols(1,p-1);
//            D.col(p-1) = zeros(n,1);
//            m = m - 1;
//        }

//        if (m == 0) C.col(m) = (y + s) / ns;
//        else C.col(m) = (y + s - C.cols(0,m-1) * D.cols(0,m-1).t() * s) / ns;

//        D.col(m) = s / ns;
//        m = m + 1;


//    }
    for (int i = 0;i < x.n_elem;i++){
        ka_vec[i] = x(i);
        kb_vec[i] = x(i);
    }
    wpa2 = wdaintegral(rho_avec, ka_vec, x_vec, y_vec, z_vec, w_vec,xc_hole_vec[hole_idx]);
    x.save("p2.out",raw_ascii);

    cout << "wda 0p: " << gj_wpa0 << endl;
    cout << "wda 1p: " << gj_wpa1 << endl;
//    cout << "wda 2p_f: " << f_wpa2 <<endl;
    cout << "wda 2p: " << wpa2 <<endl;

    n_secs = timer.toc();

    cout << "time: " << n_secs << endl;



//    ofstream logfile;
//    fstream.precision(16);
//    logfile.open("log.txt", ios::out | ios::app);
//    logfile << filename << " " << nElec << " " << ldaX << " " << wpa0 << " " << wpa1 << endl;
//    logfile.close();

  return 0;
  }

double lda_xc_hole (double k_f, double distance){

    double value;
    double k_fdis = k_f*distance;
    if (abs(k_fdis) <= datum::eps*1.5) return -1;
    value = 3 * (sin(k_fdis) / (k_fdis*k_fdis*k_fdis) - cos(k_fdis) / (k_fdis*k_fdis));
    value = -value * value;
    return value;
}
double xcp_hole(double k_f, double dis){

    double value;
    double k_fdis = k_f * dis;

    if (abs(k_f * dis) <= datum::eps*1.5) return 0;

    value = -(18*sin(k_fdis)*(sin(k_fdis)-k_fdis*cos(k_fdis))/(pow(k_f,5.0)*pow(dis,4.0)))+(54*(sin(k_fdis)-k_fdis*cos(k_fdis))*(sin(k_fdis)-k_fdis*cos(k_fdis)))/(pow(k_f,7.0)*pow(dis,6.0));

    return value;
}
double gj_xc_hole(double rho, double distance){
    double I1 = 10.2098; // magic number from constraints I1 = 4*pi*-Gamma(-2/5)/5
    double I2 = 7.42444; // magic number from constraints I2 = 4*pi*-Gamma(-3/5)/5
    double x_energy = -3.0 / 4.0 * npow(3.0 / datum::pi * rho, 1.0 / 3.0);
    //        double C = -1 / ( pow(lambda , 3.0) * rho * I2);
    double C = -1;
    double f;
    if (abs(x_energy) > 0.000000001) {
        double lambda = -I1 / (2 * I2 * x_energy);
        double u = distance / lambda;
        if (abs(u) > 0.0000001){
            f = C * ( 1 - exp(- 1 / pow(u , 6.0)));
            return f;
        }
        else
            return C;
    }
    else {

        return C;
    }

//    return -(1-exp(-1 / pow(rho*distance,5.0)));


}

double gaussian_xc_hole(double rho, double distance){
    return -(exp(-pow(rho*distance,4.0)));
}

vec k_fermi (vec rho){

    vec value = rho;
    value = pow((3 * datum::pi * datum::pi * rho), 1.0 / 3.0);
    return value;
}

double wdaintegral(vector<double> rho_vec, vector<double> k_vec, vector<double> x_vec, vector<double> y_vec, vector<double> z_vec, vector<double> w_vec, double(*xc_hole)(double,double)){

    double value = 0;
    double temp;
    double distance;
    double k_f;

    for (int i = 0; i < w_vec.size(); i++){
        for (int j = i + 1; j < w_vec.size(); j++){
            if (i == j) continue;

            distance = sqrt((x_vec[i] - x_vec[j])*(x_vec[i] - x_vec[j]) + (y_vec[i] - y_vec[j])*(y_vec[i] - y_vec[j]) + (z_vec[i] - z_vec[j])*(z_vec[i] - z_vec[j]));

            k_f = npow(0.5*(npow(k_vec[i], p) + npow(k_vec[j], p)), 1.0 / p);

//            temp = (std::sin(k_fdis) / (k_fdis*k_fdis*k_fdis) - std::cos(k_fdis) / (k_fdis*k_fdis));
//            temp = temp*temp;

            temp = xc_hole(k_f, distance);

            value += rho_vec[i] * rho_vec[j] * temp / distance * w_vec[i] * w_vec[j];

            if ((i % 100 == 0) && (j == w_vec.size()))
                cout << "* ";
        }
    }
    return value * 2;
}
double f_wdaintegral(vec f,vector<double> rho_vec, vector<double> k_vec, vector<double> x_vec, vector<double> y_vec, vector<double> z_vec, vector<double> w_vec){

    double value = 0;
    double temp;
    double distance;
    double k_fdis;

    for (int i = 0; i < w_vec.size(); i++){
        for (int j = i + 1; j < w_vec.size(); j++){
            if (i == j) continue;

            distance = sqrt((x_vec[i] - x_vec[j])*(x_vec[i] - x_vec[j]) + (y_vec[i] - y_vec[j])*(y_vec[i] - y_vec[j]) + (z_vec[i] - z_vec[j])*(z_vec[i] - z_vec[j]));

            k_fdis = npow(0.5*(npow(k_vec[i], p) + npow(k_vec[j], p)), 1.0 / p) * distance;

            temp = (std::sin(k_fdis) / (k_fdis*k_fdis*k_fdis) - std::cos(k_fdis) / (k_fdis*k_fdis));
            temp = temp*temp;

            value -= rho_vec[i] * rho_vec[j] * f(i) * temp * f(j) / distance * w_vec[i] * w_vec[j];

            if ((i % 100 == 0) && (j == w_vec.size()))
                cout << "* ";
        }
    }
    return value * 18;
}



double brents(double x1, double x2, double tol,vector<double> rho_vec, vector<double> x_vec, vector<double> y_vec, vector<double> z_vec, vector<double> w_vec, int j,double (*xc_hole)(double,double)){

    double a = x1;
    double b = x2;
    double c = x2;
    double d, e, p, q, r, s, xm, tol1;
    double fa = g_symnormintegral(rho_vec, a*ones(rho_vec.size()), x_vec, y_vec, z_vec, w_vec, j, 0 ,xc_hole);
    double fb = g_symnormintegral(rho_vec, b*ones(rho_vec.size()), x_vec, y_vec, z_vec, w_vec, j, 0 ,xc_hole);
    double fc;

//    if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))

//cout << fa << " " << fb << " " << endl;
    fc = fb;

    for (int i=0; i < 1000000 ; i++){
//        cout << "b: " << b << endl;
//        cout << "error: " << fb << endl;
        if (( fb > 0.0 && fc > 0.0) || (fb <0.0 && fc < 0.0)){
            c = a;
            fc = fa;
            e = b-a;
            d = b-a;
        }
        if (abs(fc) < abs(fb)){
            a = b;
            b = c;
            c = a;
            fa = fb;
            fb = fc;
            fc = fa;
        }
        tol1 = 2.0 * datum::eps * abs(b) + 0.5 * tol;
        xm = 0.5 * (c - b);
        if (abs(xm) <= tol1 || fb == 0.0) {
//        if (abs(fb) < tol ){
            return b;
        }
        if (abs(e) > tol1 && abs(fa) > abs(fb)){
            s = fb /fa;
            if (a == c){
                p = 2.0 * xm * s;
                q = 1.0 - s;
            }
            else{
                q = fa / fc;
                r = fb / fc;
                p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0));
                q = (q - 1.0) * (r - 1.0) * (s - 1.0);
            }
            if (p > 0.0) q = -q;
            p = abs(p);
            double min1 = 3.0 * xm * q - abs(tol1 * q);
            double min2 = abs(e * q);
            if (2.0 * p < (min1 < min2 ? min1 : min2)){
                e = d;
                d = p / q;
            }
            else {
                d = xm;
                e = d;
            }
        }
        else {
            d =xm;
            e =d;
        }
        a = b;
        fa = fb;
        if (abs(d) > tol1)
            b +=d;
        else
            if ( xm > 0 ) b += tol1;
            else if (xm < 0) b -= tol1;
        fb = g_symnormintegral(rho_vec, b*ones(rho_vec.size()), x_vec, y_vec, z_vec, w_vec, j,0,xc_hole);
    }
    return -1;

}
double nwt(double x1, double x2, double tol,vector<double> rho_vec, vector<double> x_vec, vector<double> y_vec, vector<double> z_vec, vector<double> w_vec, int j){

//    int n = rho_vec.size();
//    double y = gj_normintegral(rho_vec, x1, x_vec, y_vec, z_vec, w_vec, j);
//    double yp;
//    double s;
//    double distance;
//    int itr = 0;
//    for (int  i = 0; i < n; i++){
//        distance = sqrt((x_vec[i] - x_vec[j])*(x_vec[i] - x_vec[j]) + (y_vec[i] - y_vec[j])*(y_vec[i] - y_vec[j]) + (z_vec[i] - z_vec[j])*(z_vec[i] - z_vec[j]));
//        yp += w_vec[i]*rho_vec[i]*xcp_hole(x1,distance);
//    }
//    while (itr < 20 && abs(y) > 0.0001){
//        s = -y/yp;
//        if (abs(s) < 5)
//            x1 = x1 - y/yp;
//        else
//            if (s > 0)
//                x1 = x1 + 5;
//        else
//                x1 = x1 - 5;
//        y = g_symnormintegral(rho_vec, x1, x_vec, y_vec, z_vec, w_vec, j,1,gj_xc_hole);g_symnormintegral()
//        for (int  i = 0; i < n; i++){
//            distance = sqrt((x_vec[i] - x_vec[j])*(x_vec[i] - x_vec[j]) + (y_vec[i] - y_vec[j])*(y_vec[i] - y_vec[j]) + (z_vec[i] - z_vec[j])*(z_vec[i] - z_vec[j]));
//            yp += w_vec[i]*rho_vec[i]*xcp_hole(x1,distance);
//        }
//        itr++;

//    }
//    cout << "error: " << y << " itr: " << itr << endl;
    return x1;

}
double sym_xcp_hole(vec k_vec,double dis,int i, int j){
    double p = 1.0;
    double value;
    double k_fdis = npow(0.5*(npow(k_vec(i), p) + npow(k_vec(j), p)), 1.0 / p) * dis;
    value = (1/(pow(dis,6.0)))*3456.0*pow(k_vec(j),(p - 1))* pow((pow(k_vec(j),p) + pow(k_vec(i),p)),(-1 - 6/p)) *pow(sin(k_fdis)-dis*k_fdis * cos(k_fdis),2) - 1/pow(dis,6)* 1152.0* ((pow(k_vec(j),p) + pow(k_vec(i),p)),(-6/p)) *(-dis*k_fdis* cos(k_fdis) + sin(k_fdis))* (0.25* pow(dis,2.0)* pow(k_vec(j),p-1)*pow((pow(k_vec(j),p) + pow(k_vec(i),p)),(-1 + 2/p))*sin(k_fdis));

    return value;
}


double f_integral(vec f, vector<double> rho_vec, vector<double> k_vec, vector<double> x_vec, vector<double> y_vec, vector<double> z_vec, vector<double> w_vec, int j,double(*xc_hole)(double,double)){

    double f_x = 0;
    double tmp_k;
    double distance;
    double k_f;

    for (int i = 0; i < k_vec.size(); i++){
        if ((j == i) || (k_vec[i] < 0.0000000001) || (k_vec[j] < 0.0000000001))
            f_x += f(i)*rho_vec[i] * (-1) * w_vec[i];
        else {
            distance = sqrt((x_vec[i] - x_vec[j])*(x_vec[i] - x_vec[j]) + (y_vec[i] - y_vec[j])*(y_vec[i] - y_vec[j]) + (z_vec[i] - z_vec[j])*(z_vec[i] - z_vec[j]));

            k_f = npow(0.5*(npow(k_vec[i], p) + npow(k_vec[j], p)), 1.0 / p) * distance;

            f_x += f(i)*rho_vec[i] * xc_hole(k_f, distance) * w_vec[i];
        }
    }
    f_x = f(j)*f_x;
    f_x = f_x + 1.0;

    return f_x;
}
//vec f_gcr(vec f, vec b, vector<double> rho_vec, vector<double> k_vec, vector<double> x_vec, vector<double> y_vec, vector<double> z_vec, vector<double> w_vec){
//    vec r = b;

//    vec old_p = b;
//    vec new_p = b;
//    vec old_Ap;
//    vec new_Ap;
//    vec del_f;
//    double eps = 0.00001;
//    int n = b.size();
//    del_f.zeros(b.size());
//    old_Ap.zeros(b.size());
//    new_Ap.zeros(b.size());
//    int itr = 0;
//    double tol = 0.001;
//    double alpha;
//    double beta;
//    new_p = r;
//    for (int i = 0;i < n;i++){
//        new_Ap(i) = 1 / eps * (f_integral(f + eps * new_p, rho_vec, k_vec, x_vec, y_vec, z_vec, w_vec, i) + b(i));
//    }
//    new_p = new_p / norm(new_Ap);
//    new_Ap = new_Ap / norm(new_Ap);
//    alpha = dot(r,new_Ap);
//    del_f = del_f + alpha * new_p;
//    cout << "itr: " << itr << " gcr error: " << norm(r)/norm(b) << endl;
//    r = r - alpha * new_Ap;
//    cout << "itr: " << itr << " gcr error: " << norm(r)/norm(b) << endl;
//    while (norm(r)/norm(b) > tol && itr < 10){
//        old_p = new_p;
//        old_Ap = new_Ap;
//        new_p = r;
//        for (int i=0;i < n;i++){
//            new_Ap(i) = f_integral(f + new_p, rho_vec, k_vec, x_vec, y_vec, z_vec, w_vec, i) + b(i);
//        }
//        beta = dot(old_Ap,new_Ap);
//        new_p = new_p - beta * old_p;
//        new_Ap = new_Ap - beta * old_Ap;
//        new_p = new_p / norm(new_Ap);
//        new_Ap = new_Ap / norm(new_Ap);
//        alpha = dot(r,new_Ap);
//        del_f = del_f + alpha * new_p;
//        r = r - alpha * new_Ap;
//        itr++;
//        cout << "itr: " << itr << " gcr error: " << norm(r)/norm(b) << endl;
//    }
//    return f + del_f;


//}
vec gcr(vec b, vector<double> rho_vec, vec k_vec, vector<double> x_vec, vector<double> y_vec, vector<double> z_vec, vector<double> w_vec, double g, double tol, double maxitr,double(*xc_hole)(double,double)){
    vec r = b;

    vec old_p = b;
    vec new_p = b;
    vec old_Ap;
    vec new_Ap;
    vec del_k;
    double eps = 0.000001;
    int n = b.size();
    del_k.zeros(b.size());
    old_Ap.zeros(b.size());
    new_Ap.zeros(b.size());
    int itr = 0;

    double alpha;
    double beta;
    new_p = r;
    for (int i = 0;i < n;i++){
        new_Ap(i) = 1 / eps * (g_symnormintegral(rho_vec, k_vec + new_p * eps, x_vec, y_vec, z_vec, w_vec, i ,g,xc_hole) + b(i));
    }
    new_p = new_p / norm(new_Ap);
    new_Ap = new_Ap / norm(new_Ap);
    alpha = dot(r,new_Ap);
    del_k = del_k + alpha * new_p;
//    cout << "itr: " << itr << " gcr error: " << norm(r)/norm(b) << endl;
    r = r - alpha * new_Ap;
//    cout << "itr: " << itr << " gcr error: " << norm(r)/norm(b) << endl;
    while (norm(r)/norm(b) > tol && itr < maxitr){
        old_p = new_p;
        old_Ap = new_Ap;
        new_p = r;
        for (int i=0;i < n;i++){
            new_Ap(i) = 1 / eps * (g_symnormintegral(rho_vec, k_vec + new_p * eps, x_vec, y_vec, z_vec, w_vec, i ,g,xc_hole) + b(i));
        }
        beta = dot(old_Ap,new_Ap);
        new_p = new_p - beta * old_p;
        new_Ap = new_Ap - beta * old_Ap;
        new_p = new_p / norm(new_Ap);
        new_Ap = new_Ap / norm(new_Ap);
        alpha = dot(r,new_Ap);
        del_k = del_k + alpha * new_p;
        r = r - alpha * new_Ap;
        itr++;
//        cout << "itr: " << itr << " gcr error: " << norm(r)/norm(b) << endl;
    }

    alpha = 1;
//    vec init_b = b;
//    for (int i=0;i < n;i++){
//        b(i)=symnormintegral(rho_vec, k_vec + alpha * del_k, x_vec, y_vec, z_vec, w_vec, i);
//    }
//            cout <<"alpha: " <<alpha << " error: " << norm(b)/norm(init_b) << endl;
//    while ((norm(b) >= norm(init_b)) && (norm(alpha * del_k) > 0.001) && (alpha >= 0)){
//        for (int i=0;i < n;i++){
//            b(i)=symnormintegral(rho_vec, k_vec + alpha * del_k, x_vec, y_vec, z_vec, w_vec, i);
//        }

//        cout <<"alpha: " <<alpha << " error: " << norm(b)/norm(init_b) << endl;
//        alpha = alpha + 0.005;
//    }

    //    if (max(del_k) > 0.001)
//        del_k = del_k / max(del_k) * 0.001;

    cout << "itr: " << itr << " gcr error: " << norm(r)/norm(b) << endl;
    return k_vec + alpha * del_k;


}
double g_symnormintegral(vector<double> rho_vec, vec k_vec, vector<double> x_vec, vector<double> y_vec, vector<double> z_vec, vector<double> w_vec, int j, double g, double (*xc_hole)(double ,double)){

    double f_x = 0;
    double tmp_k;
    double distance;
    double k_f;

    for (int i = 0; i < k_vec.size(); i++){
//        if ((j == i) || (k_vec(i) < 0.0000000001) || (k_vec(j) < 0.0000000001))
//            f_x += rho_vec[i] * (-1) * w_vec[i];
//        else {
            distance = sqrt((x_vec[i] - x_vec[j])*(x_vec[i] - x_vec[j]) + (y_vec[i] - y_vec[j])*(y_vec[i] - y_vec[j]) + (z_vec[i] - z_vec[j])*(z_vec[i] - z_vec[j]));

            k_f =(1 - g) * k_vec(j) + g * (npow(0.5*(npow(k_vec(i), p) + npow(k_vec(j), p)), 1.0 / p));


            f_x += rho_vec[i] * xc_hole(k_f, distance) * w_vec[i];
//        }
    }
    f_x = f_x + 1.0;

    return f_x;

}
//double normintegral(vector<double> rho_vec, double k_vec, vector<double> x_vec, vector<double> y_vec, vector<double> z_vec, vector<double> w_vec, int j){

//    double f_x = 0;
//    double tmp_k;
//    double distance;
//    double k_fdis;

//    for (int i = 0; i < w_vec.size(); i++){
////        if ((j == i) || (k_vec < 0.0000000001))
////            f_x += rho_vec[i] * (-1) * w_vec[i];
////        else
////        {
//            distance = sqrt((x_vec[i] - x_vec[j])*(x_vec[i] - x_vec[j]) + (y_vec[i] - y_vec[j])*(y_vec[i] - y_vec[j]) + (z_vec[i] - z_vec[j])*(z_vec[i] - z_vec[j]));


//            f_x += rho_vec[i] * xc_hole(k_vec, distance) * w_vec[i];
////        }
//    }
//    f_x = f_x + 1.0;

//    return f_x;

//}
//double symnormintegral(vector<double> rho_vec, vec k_vec, vector<double> x_vec, vector<double> y_vec, vector<double> z_vec, vector<double> w_vec, int j){

//    double f_x = 0;
//    double tmp_k;
//    double distance;
//    double k_f;

//    for (int i = 0; i < k_vec.size(); i++){
////        if ((j == i) || (k_vec(i) < 0.0000000001) || (k_vec(j) < 0.0000000001))
////            f_x += rho_vec[i] * (-1) * w_vec[i];
////        else {
//            distance = sqrt((x_vec[i] - x_vec[j])*(x_vec[i] - x_vec[j]) + (y_vec[i] - y_vec[j])*(y_vec[i] - y_vec[j]) + (z_vec[i] - z_vec[j])*(z_vec[i] - z_vec[j]));

//            k_f = npow(0.5*(npow(k_vec(i), p) + npow(k_vec(j), p)), 1.0 / p);

//            f_x += rho_vec[i] * xc_hole(k_f, distance) * w_vec[i];

////        }
//    }
//    f_x = f_x + 1.0;

//    return f_x;

//}
double npow(double b,double e){
    if ( b < 0.0) return -pow(-b,e);
    else return pow(b,e);
}
double gj_normintegral(vector<double> rho_vec, double k_vec, vector<double> x_vec, vector<double> y_vec, vector<double> z_vec, vector<double> w_vec, int j){

    double f_x = 0;
    double distance;

    for (int i = 0; i < w_vec.size(); i++){
//        if ((j == i) || (k_vec < 0.0000000001))
//            f_x += rho_vec[i] * (-1) * w_vec[i];
//        else
//        {
            distance = sqrt((x_vec[i] - x_vec[j])*(x_vec[i] - x_vec[j]) + (y_vec[i] - y_vec[j])*(y_vec[i] - y_vec[j]) + (z_vec[i] - z_vec[j])*(z_vec[i] - z_vec[j]));

            f_x += rho_vec[i] * gj_xc_hole(k_vec , distance) * w_vec[i];


//        }
    }
    f_x = f_x + 1.0;

    return f_x;

}
double gj_g_symnormintegral(vector<double> rho_vec, vec k_vec, vector<double> x_vec, vector<double> y_vec, vector<double> z_vec, vector<double> w_vec, int j, double g){

    double f_x = 0;
    double tmp_k;
    double distance;
    double k_f;

    for (int i = 0; i < k_vec.size(); i++){
//        if ((j == i) || (k_vec(i) < 0.0000000001) || (k_vec(j) < 0.0000000001))
//            f_x += rho_vec[i] * (-1) * w_vec[i];
//        else {
            distance = sqrt((x_vec[i] - x_vec[j])*(x_vec[i] - x_vec[j]) + (y_vec[i] - y_vec[j])*(y_vec[i] - y_vec[j]) + (z_vec[i] - z_vec[j])*(z_vec[i] - z_vec[j]));

//            k_f =(1 - g) * k_vec(j) + g * (npow(0.5*(npow(k_vec(i), p) + npow(k_vec(j), p)), 1.0 / p));
            k_f = g * 0.5 * (gj_xc_hole(k_vec[i], distance) + gj_xc_hole(k_vec[j], distance)) + (1 - g) * gj_xc_hole(k_vec[j] , distance) ;

            f_x += rho_vec[i] * k_f * w_vec[i];

//        }
    }
    f_x = f_x + 1.0;

    return f_x;

}
double gj_wdaintegral(vector<double> rho_vec, vector<double> k_vec, vector<double> x_vec, vector<double> y_vec, vector<double> z_vec, vector<double> w_vec){

    double value = 0;
    double temp;
    double distance;
    double k_f;

    for (int i = 0; i < w_vec.size(); i++){
        for (int j = i + 1; j < w_vec.size(); j++){
            if (i == j) continue;

            distance = sqrt((x_vec[i] - x_vec[j])*(x_vec[i] - x_vec[j]) + (y_vec[i] - y_vec[j])*(y_vec[i] - y_vec[j]) + (z_vec[i] - z_vec[j])*(z_vec[i] - z_vec[j]));

            k_f = npow(0.5*(npow(k_vec[i], p) + npow(k_vec[j], p)), 1.0 / p);


            temp = gj_xc_hole(k_f,distance);

            value += rho_vec[i] * rho_vec[j] * temp / distance * w_vec[i] * w_vec[j];

            if ((i % 100 == 0) && (j == w_vec.size()))
                cout << "* ";
        }
    }
    return value * 2;
}
