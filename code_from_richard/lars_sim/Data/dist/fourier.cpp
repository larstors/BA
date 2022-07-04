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

vector<double> coor(unsigned n){
    vector<double> c(2);

    c[0] = n%L - 50; // x
    c[1] = n/L - 50; // y

    return c;
}

vector<double> rezcoor(unsigned n){
    vector <double> c(2);

    c[0] = -pi/double(L) + 2.0*pi/double(L*L) * (n%L);
    c[1] = -pi/double(L) + 2.0*pi/double(L*L) * n/L;

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

    vector<double> r(L*L);
    //vector<double> p(L*L);
    std::ifstream file("square_dens_1.txt");
    if (file.is_open()) {
        std::string line;
        while (std::getline(file, line)) {
            std::vector<std::string> density = split(line, ' ');
            vector<double> dens(density.size());
            transform(density.begin(), density.end(), dens.begin(),
                [](string const& val) {return stod(val);});

            vector<complex<double>> fourier(dens.size());
            
            complex<double> j(0,1);
            for (unsigned n = 0; n < L*L; n++){
                for (unsigned i = 0; i < L*L; i++){
                        fourier[n] += dens[i]*exp(-j * scalar(rezcoor(n), coor(i)));
                }
            }

            
            for (unsigned n = 0; n < L*L; n++){
                r[n] += pow(abs(fourier[n]), 2);
            }



        }
        file.close();
    }

    ofstream outfile;
    outfile.open("fourier_test.txt");
    for (const auto& m : r) outfile << m << " ";
    outfile << endl;
}

