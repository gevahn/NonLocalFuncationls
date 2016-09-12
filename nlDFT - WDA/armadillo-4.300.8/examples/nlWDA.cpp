#include <iostream>
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

	mat raw_data;
	raw_data.load("194.dat");
	cout << raw_data.n_cols << " " << raw_data.n_rows << endl;

	vec rho_a;
	vec rho_b;
	vec weight;
	mat xyz_coord;
	xyz_coord = raw_data.cols(0, 2);
	rho_a = raw_data.col(4);
	rho_b = raw_data.col(4);
	weight = raw_data.col(3);

	double distance;
	vec k_f;
	vec k_f_a;
	vec k_f_b;

	double integral = 0;

	k_f = k_fermi(rho_a + rho_b);
	k_f_a = k_fermi(rho_a);
	k_f_b = k_fermi(rho_b);

	for (int i = 0; i < raw_data.n_rows; i++){
		rho_vec.push_back(rho_a(i) + rho_b(i));
		rho_avec.push_back(rho_a(i));
		rho_bvec.push_back(rho_b(i));
		w_vec.push_back(weight(i));
		x_vec.push_back(xyz_coord(i, 0));
		y_vec.push_back(xyz_coord(i, 1));
		z_vec.push_back(xyz_coord(i, 2));
		k_vec.push_back(k_f(i));
		ka_vec.push_back(k_f_a(i));
		kb_vec.push_back(k_f_b(i));
	}
	

	double value; double k_fdis;
	double k_f_adis;
	double k_f_bdis;

	timer.tic();
	/*int block_size = 10;
	int n_block = raw_data.n_rows / block_size;
	for (int iblock = 0; iblock < n_block; iblock++){
		for (int jblock = 0; jblock < n_block; jblock++){			
			for (int i = iblock * block_size; i < iblock * block_size + block_size; i++){
				if (k_vec.at(i) < 0.000000001) continue;		
				for (int j = jblock * block_size; j < jblock * block_size + block_size; j++){
					if (i == j) continue;
							
					distance = sqrt((x_vec.at(i) - x_vec.at(j))*(x_vec.at(i) - x_vec.at(j)) + (y_vec.at(i) - y_vec.at(j))*(y_vec.at(i) - y_vec.at(j)) + (z_vec.at(i) - z_vec.at(j))*(z_vec.at(i) - z_vec.at(j)));
					if (distance <= datum::eps * 10) continue;
					k_fdis = k_vec.at(i) * distance;
					if (k_fdis <= datum::eps*1.5) return 0;
					else{

						value = 3 * (sin(k_fdis) / (k_fdis*k_fdis*k_fdis) - cos(k_fdis) / (k_fdis*k_fdis));
						value = -value * value;
					}
					integral += rho_vec.at(i) * rho_vec.at(j) * value / distance * w_vec.at(i) * w_vec.at(j);
					if (!is_finite(integral))
						break;
				}
			}					
		}
	}*/
	int block_size = 100;
	int n_block = raw_data.n_rows / block_size;
	/*n_block = 10;*/

	float p = 6.0;

	/*int j = 1;*/
	double f_x = 1;
	double fp_x;
	double x_0;
	double max_i;
	
	double max = 0;
	for (int j = 0; j<raw_data.n_rows; j++){
		int itr = 0;
		x_0 = ka_vec[j];
		f_x = 1;
		while ((abs(f_x) > 0.0001) && (itr < 30)){
			f_x = 0;
			fp_x = 0;
			for (int i = 0; i < raw_data.n_rows; i++){
				distance = sqrt((x_vec[i] - x_vec[j])*(x_vec[i] - x_vec[j]) + (y_vec[i] - y_vec[j])*(y_vec[i] - y_vec[j]) + (z_vec[i] - z_vec[j])*(z_vec[i] - z_vec[j]));
				if (distance == 0) continue;
				f_x += rho_avec[i] * xc_hole(x_0*distance)*w_vec[i];
				fp_x += rho_avec[i] * xcp_hole(x_0*distance)*w_vec[i];
			}
			f_x = f_x + 1;
			if (f_x / fp_x > 10)
				x_0 = x_0 - 10;
			else if (f_x / fp_x < -10)
				x_0 = x_0 + 10;
			else
				x_0 = x_0 - f_x / fp_x;

			itr++;
		}
		if (is_finite(x_0)){
			ka_vec[j] = x_0;
			if (abs(f_x) > max){
				max = abs(f_x);
				max_i = j;
			}
		}

		
	}
	cout << max << endl;

	for (int i = 0; i < raw_data.n_rows; i++){
		kb_vec[i] = ka_vec[i];
		if (!is_finite(ka_vec[i]))
			cout << i << endl;
	}
	
	for (int iblock = 0; iblock < n_block; iblock++){
		for (int jblock = 0; jblock < n_block; jblock++){
			for (int i = iblock * block_size; i < iblock * block_size + block_size; i++){
				if (ka_vec[i] < 0.000000001) continue;
				if (kb_vec[i] < 0.000000001) continue;
				for (int j = jblock * block_size; j < jblock * block_size + block_size; j++){
					if (i == j) continue;

					distance = sqrt((x_vec[i] - x_vec[j])*(x_vec[i] - x_vec[j]) + (y_vec[i] - y_vec[j])*(y_vec[i] - y_vec[j]) + (z_vec[i] - z_vec[j])*(z_vec[i] - z_vec[j]));
					k_f_adis = std::pow(0.5*(std::pow(ka_vec[i], p) + std::pow(ka_vec[j], p)), 1.0 / p) * distance;
					k_f_bdis = std::pow(0.5*(std::pow(kb_vec[i], p) + std::pow(kb_vec[j], p)), 1.0 / p) * distance;
					/*k_f_adis = ka_vec[i] * distance;
					k_f_bdis = kb_vec[i] * distance;*/
					if (k_f_adis <= datum::eps*10) continue;
					if (k_f_bdis <= datum::eps * 10) continue;
					value = 3 * (std::sin(k_f_adis) / (k_f_adis*k_f_adis*k_f_adis) - std::cos(k_f_adis) / (k_f_adis*k_f_adis));
					value = -value * value;
					
					integral += rho_avec[i] * rho_avec[j] * value / distance * w_vec[i] * w_vec[j];

					value = 3 * (std::sin(k_f_bdis) / (k_f_bdis*k_f_bdis*k_f_bdis) - std::cos(k_f_bdis) / (k_f_bdis*k_f_bdis));
					value = -value * value;

					integral += rho_bvec[i] * rho_bvec[j] * value / distance * w_vec[i] * w_vec[j];

					if (!is_finite(integral)){
						cout << distance << endl;
						cout << value << endl;
						cout << k_f_adis << endl;
						cout << k_f_bdis << endl;
						cout << ka_vec[i] << endl;
						cout << ka_vec[j] << endl;
						cout << kb_vec[i] << endl;
						cout << kb_vec[j] << endl;
						break;
						
					}
						
				}
			}
		}
	}
	double nElec;
	double ldaX;
	
	n_secs = timer.toc();

	cout << "time: " << n_secs << endl;

	nElec = dot(rho_a + rho_b, weight);
	ldaX = dot(pow(rho_a + rho_b, (4 / 3.0)), weight);
	ldaX = (pow(3 / datum::pi, (1.0 / 3))) * ldaX * (-0.75);
	
	cout << "number of elctrons: " << nElec << endl;
	cout << "lda: " << ldaX << endl;

	
	cout << "wda 0p: " << 0.5*integral << endl;
	
	
	/*vec rho;
	rho = (rho_a + rho_b);
	mat R;
	mar W;
	R = rho*rho.t();
	W = weight*weight.t();*/
	/*timer.tic();
	
	for (int i = 0; i < raw_data.n_rows; i++){
		if (k_f(i) < 0.000000001) continue;
		for (int j = 0; j < raw_data.n_rows; j++){
			if (i == j) continue;
			distance = norm(xyz_coord.row(i) - xyz_coord.row(j));
			if (distance <= datum::eps * 1.5) continue;
			k_fdis = k_f(i) * distance;
			if (k_fdis <= datum::eps*1.5) return 0;
			else{
				
				value = 3 * (sin(k_fdis) / (k_fdis*k_fdis*k_fdis) - cos(k_fdis) / (k_fdis*k_fdis));
				value = -value * value;
			}
			integral += (rho_a(i) + rho_b(i)) * (rho_a(j) + rho_b(j)) * value / distance * weight(i)*weight(j);
			if (!is_finite(integral))
				break;
		}
	}
	n_secs = timer.toc();
	cout << n_secs << endl;
	cout << 0.5*integral << endl;*/

	mat out;
	out << nElec << endr << ldaX << endr << 0.5*integral << endr;
	out.save("302.out", raw_ascii);
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
	value = pow((3 * pow(datum::pi, 2.0) * rho),1.0/3.0);
	return value;
}