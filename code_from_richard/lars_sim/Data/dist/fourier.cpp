/*
Program to calculate the fourier transform (and maybe the structure factor as well?)

*/


#include <vector>
#include <list>
#include <map>
#include <valarray>
#include <numeric>
#include <functional>
#include <random>
#include <iostream>
#include <cassert>
#include <fstream>
#include <iostream>
#include <cmath>
#include <string> 
#include <iomanip>
#include <sstream>
#include <complex>
#include <fstream>

const unsigned L = 100;
const double pi = std::acos(-1);

using namespace std;

std::vector<std::string> split(std::string text, char delim) {
    std::string line;
    std::vector<std::string> vec;
    std::stringstream ss(text);
    while(std::getline(ss, line, delim)) {
        vec.push_back(line);
    }
    return vec;
}

vector<double> coor(int n){
    vector<double> c(2);

    c[0] = n%L - 50; // x
    c[1] = n/L - 50; // y

    return c;
}

vector<double> rezcoor(int n){
    vector <double> c(2);

    c[0] = -pi/double(L) + 2.0*pi/double(L*L) * (n%L);
    c[1] = -pi/double(L) + 2.0*pi/double(L*L) * n/L;

    return c;
}

vector<double> rezcoor2(int n, int i){
    vector <double> c(2);

    c[0] = -pi/double(L) + 2.0*pi/double(L*L) * n;
    c[1] = -pi/double(L) + 2.0*pi/double(L*L) * i;

    return c;
}



vector<double> coor2(int n, int i){
    vector<double> c(2);

    c[0] = n - 50; // x
    c[1] = i - 50; // y

    return c;
}


double scalar(vector<double> a, vector<double> b){
    int l = a.size();
    double sca = 0;
    for (unsigned i = 0; i < l; i++){
        sca += a[i] * b[i];
    }
    return sca;
}



int main(){
    unsigned oo = 0;

    vector<double> r(L*L);
    double res[L][L];

    for (unsigned n = 0; n < L; n++){
        for (unsigned k = 0; k < L; k++){
            res[n][k] = 0;
        }
    }

    double p[L][L];

    for (unsigned n = 0; n < L; n++){
        for (unsigned k = 0; k < L; k++){
            p[n][k] = 0;
        }
    }
    std::ifstream file("square_dens_3.txt");
    if (file.is_open()) {
        std::string line;
        while (std::getline(file, line) && oo == 0) {
            std::vector<std::string> density = split(line, ' ');
            vector<double> dens(density.size());
            transform(density.begin(), density.end(), dens.begin(),
                [](string const& val) {return stod(val);});


            double rho = 0;

            for (const auto& m : dens) rho += m;

            rho = rho / (double(dens.size()));


            //vector<complex<double>> fourier(dens.size());
            double re[L][L];
            double im[L][L];
            for (unsigned n = 0; n < L; n++){
                for (unsigned k = 0; k < L; k++){
                    re[n][k] = 0;
                    im[n][k] = 0;
                }
            }

            /*
            complex<double> j(0,1);
            
            for (unsigned n = 0; n < L*L; n++){
                for (unsigned i = 0; i < L*L; i++){
                        fourier[n] += dens[i]*exp(-j * scalar(rezcoor(n), coor(i)));
                }
            }
            

           complex<double> fourier[L][L];

           for (unsigned n = 0; n < L; n++){
                for (unsigned k = 0; k < L; k++){
                    fourier[n][k] = 0;
                }
            }
            */

            double x[L];
            double q[L];
            for (int i = 0; i < L; i++){
                x[i] = -50 + i;
                q[i] = L*(- pi / double(L) + double(2*pi)/double(L*L) * i);
            }

            
            int count = 0;

           double d[L][L];
           for (int n = 0; n < L; n++){
                for (int k = 0; k < L; k++){
                    /*
                    if (sqrt((x[n] + 1)*(x[n] + 1) + (x[k] + 1)*(x[k] + 1)) < 5){
                        d[n][k] = 1;
                        count++;
                    }
                    else d[n][k] = 0;
                    */
                   d[n][k] = (dens[n*L + k]-rho);
                }
            }




            for (int n = 0; n < L; n++){
                for (int k = 0; k < L; k++){
                    for (int l = 0; l < L; l++){
                        for (int i = 0; i < L; i++){
                            re[n][k] += d[l][i] * cos(q[n] * x[l] + q[k] * x[i]);
                            im[n][k] -= d[l][i] * sin(q[n] * x[l] + q[k] * x[i]);
                            //fourier[n][k] += d[l][i] * exp(-j * (q[n] * x[l] + q[k] * x[i]));
                            
                        }
                    }
                }
            }
            
            /*
            for (unsigned n = 0; n < L*L; n++){
                r[n] += pow(abs(fourier[n]), 2);
                //p[n] += abs(fourier[n]);
            }
            */

            for (unsigned n = 0; n < L; n++){
                for (unsigned k = 0; k < L; k++){
                    //p[n][k] = 1.0 / double(count*count)*pow(abs(fourier[n][k]), 2);
                    res[n][k] += pow(re[n][k], 2) + pow(im[n][k], 2);//re[n][k]*re[n][k] + im[n][k]*im[n][k];
                    
                    //res[n][k] *= 1.0 / double(count*count);
                    //cout << p[n][k] << " " << res[n][k] << " ";
                    //res[n][k] = d[n][k];
                }
                //cout << endl;
            }

            oo = 1;

        }
        file.close();
    }

    ofstream outfile;
    outfile.open("fourier_sq_3.txt");
    //for (const auto& m : r) outfile << m << " ";
    for (unsigned n = 0; n < L; n++){
                for (unsigned k = 0; k < L; k++){
                    outfile << res[n][k] << " ";
        }
    }
    outfile << endl;

    //for (const auto& m : r) outfile1 << m << " ";
    //outfile1 << endl;

}

