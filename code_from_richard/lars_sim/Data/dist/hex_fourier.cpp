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
const double c30 = cos(pi/6.0);
const double s30 = sin(pi/6.0);

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


vector<double> trafo(int n, int i){
    vector <double> coor;
    int y = n/L;
    int x = 2 * (n%L) + i - y;
    x += 1;
    coor.push_back(x*c30);
    if (y%2 == 1){
        if (x%2 == 1) coor.push_back(y * (1 + s30));
        else coor.push_back((y * 3 + 1)/2);
    }
    else {
        if (x%2 == 1) coor.push_back(1.5 * y + 0.5);
        else coor.push_back(1.5 * y);
    }

    return coor;
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

    double x[L*L][2];
    double y[L*L][2];
    double qy[L];
    double qx[L];

    for (int n = 0; n < L*L; n++){
        for (int i = 0; i < 2; i++){
            x[n][i] = trafo(n, i)[0];
            y[n][i] = trafo(n, i)[1];
        }
    }
    for (int n = 0; n < L; n++){
        qx[n] = 3*(-pi + 2.0*pi*double(n)/double(L));
        qy[n] = 3*(-pi + 2.0*pi*double(n)/double(L));
    }

    double p[L][L];

    for (unsigned n = 0; n < L; n++){
        for (unsigned k = 0; k < L; k++){
            p[n][k] = 0;
        }
    }
    std::ifstream file("hex_dens_3.txt");
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

            /*
            complex<double> j(0,1);
            
            complex<double> fourier[L][L];

           for (unsigned n = 0; n < L; n++){
                for (unsigned k = 0; k < L; k++){
                    fourier[n][k] = 0;
                }
            }
            */
            
            int count = 0;

           double d[L*L][2];
           for (int n = 0; n < L*L; n++){
                for (int k = 0; k < 2; k++){
                   d[n][k] = (dens[2*n + k] - rho);
                }
            }


            for (int n = 0; n < L; n++){
                for (int k = 0; k < L; k++){
                    for (int l = 0; l < L*L; l++){
                        for (int i = 0; i < 2; i++){
                            re[n][k] += d[l][i] * cos(qx[k] * x[l][i] + qy[n] * y[l][i]);
                            im[n][k] -= d[l][i] * sin(qx[k] * x[l][i] + qy[n] * y[l][i]);
                            //fourier[n][k] += d[l][i] * exp(-j * (qx[k] * x[l*L + i] + qy[n*L + k] * y[l]));
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

            //oo = 1;

        }
        file.close();
    }

    ofstream outfile, coor;
    outfile.open("fourier_hx_3.txt");
    coor.open("hex_coor.txt");
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

