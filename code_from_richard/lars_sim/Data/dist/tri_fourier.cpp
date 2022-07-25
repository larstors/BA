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



int main(){
    unsigned oo = 0;

    vector<double> r(L*L);
    double res[L][L];

    for (unsigned n = 0; n < L; n++){
        for (unsigned k = 0; k < L; k++){
            res[n][k] = 0;
        }
    }

    double x[L*L];
    double y[L];
    double qy[L*L];
    double qx[L];
    for (int i = 0; i < L; i++){
        for (int k = 0; k < L; k++) x[L*i + k] = k - 0.5*i - 25;
        y[i] = -sqrt(3.0)/2.0 * 50 + sqrt(3.0)/2.0 * i;

        qx[i] = 4*(-pi + 2.0*pi*double(i)/double(L));
        for (int k = 0; k < L; k++) qy[L*i + k] = 4*(-pi + 2.0*pi*double(i)/double(L));
    }

    double p[L][L];

    for (unsigned n = 0; n < L; n++){
        for (unsigned k = 0; k < L; k++){
            p[n][k] = 0;
        }
    }
    ofstream wow;
    wow.open("tri_wtf.txt");
    std::ifstream file("tri_dens_3.txt");
    if (file.is_open()) {
        std::string line;
        while (std::getline(file, line) && oo == 0) {
            std::vector<std::string> density = split(line, ' ');
            vector<double> dens(density.size());
            transform(density.begin(), density.end(), dens.begin(),
                [](string const& val) {return stod(val);});


            double rho = 0;

            for (const auto& m : dens) rho += m;


            rho = rho / double(dens.size());


            double re[L][L];
            double im[L][L];
            for (unsigned n = 0; n < L; n++){
                for (unsigned k = 0; k < L; k++){
                    re[n][k] = 0;
                    im[n][k] = 0;
                }
            }


            complex<double> j(0,1);
            
            complex<double> fourier[L][L];

           for (unsigned n = 0; n < L; n++){
                for (unsigned k = 0; k < L; k++){
                    fourier[n][k] = 0;
                }
            }

            
            int count = 0;

           double d[L][L];
           for (int n = 0; n < L; n++){
                for (int k = 0; k < L; k++){
                   d[n][k] = (dens[n*L + k] - rho);
                   wow << dens[n*L + k] << " ";
                }
            }

            wow << endl;

            for (int n = 0; n < L; n++){
                for (int k = 0; k < L; k++){
                    for (int l = 0; l < L; l++){
                        for (int i = 0; i < L; i++){
                            
                            re[n][k] += d[l][i] * cos(qx[k] * x[l*L + i] + qy[n*L + k] * y[l]);
                            im[n][k] -= d[l][i] * sin(qx[k] * x[l*L + i] + qy[n*L + k] * y[l]);
                            fourier[n][k] += d[l][i] * exp(-j * (qx[k] * x[l*L + i] + qy[n*L + k] * y[l]));
                        }
                    }
                }
            }
            

            for (unsigned n = 0; n < L; n++){
                for (unsigned k = 0; k < L; k++){
                    res[n][k] += pow(re[n][k], 2) + pow(im[n][k], 2);
                    //cout << res[n][k] << " " << pow(abs(fourier[n][k]), 2) << endl;
                    //if (pow(re[n][k], 2) + pow(im[n][k], 2) > 4700) cout << pow(re[n][k], 2) << " " << pow(im[n][k], 2) << endl;
                }
                //cout << endl;
            }

            oo = 1;

        }
        file.close();
    }

    ofstream outfile, coor;
    outfile.open("fourier_tr_3.txt");
    coor.open("tri_coor.txt");
    for (unsigned n = 0; n < L; n++){
                for (unsigned k = 0; k < L; k++){
                    outfile << res[n][k] << " ";
        }
    }
    outfile << endl;

    unsigned count = 0;
    for (const auto& m : y){
        for (unsigned o = 0; o < L; o++){
            coor << m << " ";
        }
    }
    coor << endl;

    count = 0;

    for (const auto& m : x){
        coor << m << " ";
    }
    coor << endl;

    count = 0;
    
    for (const auto& m : qy){
        coor << m << " ";
    }
    coor << endl;
    
    
    for (unsigned o = 0; o < L; o++){
        for (const auto& m : qx){
            coor << m << " ";
        }
    }
    coor << endl;
    

}

