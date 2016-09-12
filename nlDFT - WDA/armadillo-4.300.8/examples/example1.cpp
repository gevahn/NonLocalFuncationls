#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <cmath>

#include "armadillo"

using namespace arma;
using namespace std;

double xc_hole (double k_fdis);
double xcp_hole(double k_fdis);
double r12(double X1);
vec k_fermi (vec rho);

double integral(vector<double> rho_vec, vector<double> k_vec, vector<double> x_vec, vector<double> y_vec, vector<double> z_vec, vector<double> w_vec);

int main(int argc, char** argv)
  {
	wall_clock timer;
	double n_secs;
	vector<double> rho_vec;
	vector<double> rho_avec;
	vector<double> rho_bvec;
	vector<double> w_vec;
	vector<double> x_vec;
	vector<double> y_vec;
	vector<double> z_vec;
	vector<double> k_vec;
	vector<double> ka_vec;
	vector<double> kb_vec;
	
	cout.precision(10);
	
	
	
	string filename;
	if (argc == 1){
		cout << "Insert density file name: ";
		
		cin >> filename;
	}
	else {
		filename = argv[1];
	}
	
	
	mat raw_data;
	raw_data.load( filename + ".dat");
	cout << filename << endl;
	cout << raw_data.n_cols << " " << raw_data.n_rows << endl;


	vec rho_a;
	vec rho_b;
	vec weight;
	mat xyz_coord;
	xyz_coord = raw_data.cols(0, 2);
	rho_a = raw_data.col(4);
	rho_b = raw_data.col(5);
	weight = raw_data.col(3);

	double distance;
	vec k_f;
	vec k_f_a;
	vec k_f_b;

	double nElec;
	double ldaX;

	nElec = dot(rho_a + rho_b, weight);
	ldaX = 0.5*(dot(pow(2 * rho_a, (4 / 3.0)), weight) + dot(pow(2 * rho_b, (4 / 3.0)), weight));
	ldaX = (pow(3 / datum::pi, (1.0 / 3))) * ldaX * (-0.75);

	cout << "number of elctrons: " << nElec << endl;
	cout << "lda: " << ldaX << endl;

	k_f = k_fermi(rho_a + rho_b);
	k_f_a = k_fermi(2*rho_a);
	k_f_b = k_fermi(2*rho_b);

	for (int i = 0; i < raw_data.n_rows; i++){
		if (weight(i) > datum::eps * 10) {
			rho_vec.push_back(rho_a(i) + rho_b(i));
			rho_avec.push_back(2*rho_a(i));
			rho_bvec.push_back(2*rho_b(i));
			w_vec.push_back(weight(i));
			x_vec.push_back(xyz_coord(i, 0));
			y_vec.push_back(xyz_coord(i, 1));
			z_vec.push_back(xyz_coord(i, 2));
			k_vec.push_back(k_f(i));
			ka_vec.push_back(k_f_a(i));
			kb_vec.push_back(k_f_b(i));
		}
	}
	

	double value; double k_fdis;
	double k_f_adis;
	double k_f_bdis;

	timer.tic();
	
	int block_size = 100;
	int n_block = w_vec.size() / block_size;
	/*n_block = 10;*/

	float p = 5.0;

	
	float wpa0 = 0.5*(integral(rho_avec, ka_vec, x_vec, y_vec, z_vec, w_vec) + integral(rho_bvec, kb_vec, x_vec, y_vec, z_vec, w_vec));
	
	cout << wpa0 << endl;

	double f_x = 1;
	double fp_x;
	double x_0 = ka_vec[0];
	double max_i;
	
	double max = 0;
	float trust = 1.0;
	for (int j = 0; j< w_vec.size(); j++){
	/*int j = 6863;*/
		int itr = 0;
		x_0 = ka_vec[j];
		/*cout << "k_f initial: " << x_0 << endl;*/
		f_x = 0.1;
		if (rho_avec[j] < 0.00001) continue;
		while ((abs(f_x) > 0.0001) && (itr < 100)){
			f_x = 0;
			fp_x = 0;
			for (int i = 0; i < raw_data.n_rows; i++){
				if (j == i) continue;
				distance = sqrt((x_vec[i] - x_vec[j])*(x_vec[i] - x_vec[j]) + (y_vec[i] - y_vec[j])*(y_vec[i] - y_vec[j]) + (z_vec[i] - z_vec[j])*(z_vec[i] - z_vec[j]));
				
				f_x += rho_vec[i] * xc_hole(x_0*distance)*w_vec[i];
				fp_x += rho_vec[i] * xcp_hole(x_0*distance)*w_vec[i];
			}
			f_x = f_x + 1.0;
			
			if (f_x / fp_x > trust)
				x_0 = x_0 - trust;
			else if (f_x / fp_x < -trust)
				x_0 = x_0 + trust;
			else
				x_0 = x_0 - f_x / fp_x;

			/*ofstream nfile;
			nfile.open("nlog.txt", ios::out | ios::app );
			nfile << itr << " " << f_x << " " << x_0 << endl;*/
			itr++;

		}
		/*cout << "k_f final: " << x_0 << endl;
		cout << "itr:  " << itr << endl;*/
		if (is_finite(x_0)){
			if (abs(x_0) > 50) x_0 = k_f(j);
			if (x_0 < 0) x_0 = -x_0;
			k_vec[j] = x_0;
			k_f(j) = f_x;
			if (abs(f_x) > max){
				max = abs(f_x);
				max_i = j;
			}
		}

		
	}
		
		cout << "error: " << max << endl;
	cout << "location: " << max_i << endl;

	//for (int i = 0; i < w_vec.size(); i++){
	//	kb_vec[i] = ka_vec[i];
	//	if (!is_finite(ka_vec[i]))
	//		cout << i << endl;
	//	/*k_f(i) = ka_vec[i];*/
	//}

	k_f.save("kf.out", raw_ascii);

	
	float wpa1 = 0.5*(integral(rho_avec, ka_vec, x_vec, y_vec, z_vec, w_vec) + integral(rho_bvec, kb_vec, x_vec, y_vec, z_vec, w_vec));

	
	
	
	n_secs = timer.toc();

	cout << "time: " << n_secs << endl;

	

	
	cout << "wda 0p: " << wpa0 << endl;
	cout << "wda 1p: " << wpa1 << endl;
	
	
	
	
	mat out;
	out << nElec << endr << ldaX << endr << wpa0 << endr << wpa1 << endr;
	out.save(filename + ".out", raw_ascii);
	

	ofstream logfile;
	logfile.open("log.txt", ios::out | ios::app);
	logfile << filename << " " << nElec << " " << ldaX << " " << wpa0 << " " << wpa1 << endl;
	int test;
	cin >> test;


  return 0;
  }

double xc_hole (double k_fdis){

	double value;
	if (k_fdis <= datum::eps*1.5) return 0;
	value = 3 * (sin(k_fdis) / (k_fdis*k_fdis*k_fdis) - cos(k_fdis) / (k_fdis*k_fdis));
	value = -value * value;
	return value;
}
double xcp_hole(double k_fdis){

	double value;
	if (k_fdis <= datum::eps*1.5) return 0;

	value = -27 * (sin(k_fdis)*pow(sin(k_fdis) - k_fdis*cos(k_fdis), 2.0)*(3 * k_fdis*cos(k_fdis) + (k_fdis*k_fdis - 3)*sin(k_fdis))) / pow(k_fdis, 10.0);
	
	return value;
}
vec k_fermi (vec rho){

	vec value = rho;
	value = pow((3 * datum::pi * datum::pi * rho), 1.0 / 3.0);
	return value;
}

double integral(vector<double> rho_vec, vector<double> k_vec, vector<double> x_vec, vector<double> y_vec, vector<double> z_vec, vector<double> w_vec){
	float p = 5.0;
	double value = 0;
	double temp;
	double distance;
	double k_fdis;

	for (int i = 0; i < w_vec.size(); i++){
		for (int j = i + 1; j < w_vec.size(); j++){
			if (i == j) continue;

			distance = sqrt((x_vec[i] - x_vec[j])*(x_vec[i] - x_vec[j]) + (y_vec[i] - y_vec[j])*(y_vec[i] - y_vec[j]) + (z_vec[i] - z_vec[j])*(z_vec[i] - z_vec[j]));
				
			k_fdis = std::pow(0.5*(std::pow(k_vec[i], p) + std::pow(k_vec[j], p)), 1.0 / p) * distance;

			temp = (std::sin(k_fdis) / (k_fdis*k_fdis*k_fdis) - std::cos(k_fdis) / (k_fdis*k_fdis));
			temp = temp*temp;

			value -= rho_vec[i] * rho_vec[j] * temp / distance * w_vec[i] * w_vec[j];
									
					
		}
	}
	return value * 9 / 2;
}