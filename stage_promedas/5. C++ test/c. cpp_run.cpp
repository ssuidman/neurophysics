//  Precompiling:
//  g++ -std=c++11 -o d.\ cpp_run.o c.\ cpp_run.cpp
//  Running the file for 3 different cases:
//  ./d.\ cpp_run.o 1
//  ./d.\ cpp_run.o 2
//  ./d.\ cpp_run.o 3
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <chrono>

using namespace std;

int main(int argc, char* argv[]) {
    int cases = stoi(argv[1]);

    vector<double> previn;
    vector<vector<double> > pfmin;
    vector<double> pfminneg;
    vector<double> posterior_julia;
    vector<double> posterior;
    double dt_julia;

    ifstream previnFile;
    ifstream pfminFile;
    ifstream pfminnegFile;
    ifstream posterior_File;
    ifstream dtFile;

    if (cases == 1){
        cout << "Reading case 1" << endl;
        previnFile.open("../variables/cpp/patient404_case_1_previn.csv");
        pfminFile.open("../variables/cpp/patient404_case_1_pfmin.csv");
        pfminnegFile.open("../variables/cpp/patient404_case_1_pfminneg.csv");
        posterior_File.open("../variables/cpp/patient404_case_1_posterior_BF.csv");
        dtFile.open("../variables/cpp/patient404_case_1_dt_exp_sum_log.csv");
    } else if (cases == 2){
        cout << "Reading case 2" << endl;
        previnFile.open("../variables/cpp/patient404_case_2_previn.csv");
        pfminFile.open("../variables/cpp/patient404_case_2_pfmin.csv");
        pfminnegFile.open("../variables/cpp/patient404_case_2_pfminneg.csv");
        posterior_File.open("../variables/cpp/patient404_case_2_posterior_BF.csv");
        dtFile.open("../variables/cpp/patient404_case_2_dt_exp_sum_log.csv");
    } else if (cases == 3){
        cout << "Reading case 3" << endl;
        previnFile.open("../variables/cpp/patient404_case_3_previn.csv");
        pfminFile.open("../variables/cpp/patient404_case_3_pfmin.csv");
        pfminnegFile.open("../variables/cpp/patient404_case_3_pfminneg.csv");
        posterior_File.open("../variables/cpp/patient404_case_3_posterior_BF.csv");
        dtFile.open("../variables/cpp/patient404_case_3_dt_exp_sum_log.csv");
    }

    string line;
    while (getline(previnFile, line, '\n')) {previn.push_back(stold(line));}
    while (getline(pfminFile, line)) {
        vector<double> row;
        istringstream ss(line); // Use a stringstream to split the line
        string cell;
        while (getline(ss, cell, ',')) {row.push_back(stod(cell));}
        pfmin.push_back(row);
    }
    while (getline(pfminnegFile, line, '\n')) {pfminneg.push_back(stod(line));}
    while (getline(posterior_File, line, '\n')) {
        posterior_julia.push_back(stold(line));
        posterior.push_back(stold(line)); // to have a vector of the same type and size as posterior_julia
    }
    while (getline(dtFile, line, '\n')) {dt_julia = stod(line);}

    int m = (!pfmin.empty()) ? pfmin.size() : 0;
    int n = previn.size();
    vector<vector<double> > prev(n + 1, previn);
    for (int i = 0; i < n; i++) {prev[i + 1][i] = 1.0;}
    vector<double> prevend(n + 1, 0.0);
    prevend[0] = 1.0;
    for (int i = 1; i <= n; i++) {prevend[i] = previn[i - 1];}
    vector<vector<double> > prevminneg(n + 1, vector<double>(n, 0.0));
    for (int i = 0; i <= n; i++) {
        for (int j = 0; j <= n; j++) {prevminneg[i][j] = (!pfminneg.empty()) ? (prev[i][j] * pfminneg[j]) : (prev[i][j]);}
    }
    vector<double> pfplus(n + 1, 0.0);
    vector<int> myset;
    double dt_cpp = 0.0;

    cout << "Implementing quickscore algorithm..." << endl;
    auto start = chrono::steady_clock::now();
    double sumExp;
    double product;
    for (int i = 0; i < (1 << m); i++) {
        vector<int> v(m, 0);
        for (int j = 0; j < m; j++) {
            v[j] = (i >> j) & 1;
        }
        myset.clear();
        for (int j = 0; j < m; j++) {
            if (v[j] == 1) {
                myset.push_back(j);
            }
        }
        // cout << "i=" << i << endl;
        for (int j = 0; j < n+1; j++) {
            sumExp = 0.0;
            for (int k = 0; k < n; k++) {
                product = 1.0;
                for (int l : myset) {product *= pfmin[l][k];}
                sumExp += log(1e-50 + product * prevminneg[j][k] + 1-prev[j][k]);
            }
            pfplus[j] += pow(-1, myset.size()) * exp(sumExp);
        }
    }
    auto end = chrono::steady_clock::now();
    dt_cpp = chrono::duration_cast<chrono::milliseconds>(end - start).count() / 1000.0;

    vector<double> P_joint(n, 0.0);
    for (int i = 1; i < n+1; i++) {
        P_joint[i - 1] = pfplus[i] * previn[i - 1];
        posterior[i-1] = P_joint[i-1] / pfplus[0];
    }

    double maxDifference = 0.0;
    for (int i = 0; i < posterior.size(); i++) {
        double diff = abs(static_cast<double>(posterior_julia[i]) - posterior[i]);
        maxDifference = max(maxDifference, diff);
    }


    // for (const double& value : prevend) {std::cout << value << " ";} cout << endl;
    // cout << "previn:\t\t" << previn.size() << endl;
    // cout << "pfmin:\t\t(" << pfmin.size() << "," << pfmin[1].size() << ")" << endl;
    // cout << "pfminneg:\t" << pfminneg.size() << endl;
    // cout << "posterior_:\t" << posterior.size() << endl;
    // cout << "prev:\t\t(" << prev.size() << "," << prev[773].size() << ")" << endl;
    // cout << "prevend:\t" << prevend.size() << endl;
    // cout << "prevminneg:\t(" << prev.size() << "," << prev[773].size() << ")" << endl;
    cout << "Size of \n\tfloat:\t\t " << sizeof(float) << " bytes" << endl;
    cout << "\tdouble:\t\t " << sizeof(double) << " bytes" << endl;
    cout << "\tlong double:\t " << sizeof(long double) << " bytes" << endl << endl;
    cout << "###############################################" << endl;
    cout << "######## Results for patient 404 (m=7) ########" << endl;
    cout << "###############################################" << endl;
    int q = 4;
    cout << "posterior: \t\t P_joint:\n";
    for (int i = 0; i < q; i++) {std::cout << "\t" << posterior[i] << "\t\t" << P_joint[i] << endl;} cout << "\t..." << endl;
    for (int i = 0; i < q; i++) {std::cout << "\t" << posterior[posterior.size()-q+i] << "\t\t" << P_joint[P_joint.size()-q+i] << endl;};
    cout << "Running time:\n\t" << "C++:\t" << dt_cpp << " sec\n\tJulia:\t" << dt_julia << " sec" << endl;
    cout << "max-abs difference of posterior C++ and 'prod BF':\n\t" << maxDifference << endl;

    ofstream output_file;
    if (cases == 1) {
        cout << "Writing case 1" << endl;
        output_file.open("../variables/cpp/cpp_output_case_1.csv");
    } else if (cases == 2) {
        cout << "Writing case 2" << endl;
        output_file.open("../variables/cpp/cpp_output_case_2.csv");
    } else if (cases == 3) {
        cout << "Writing case 3" << endl;
        output_file.open("../variables/cpp/cpp_output_case_3.csv");
    }

    output_file << "Posterior,P_joint,dt" << endl;
    for (size_t i = 0; i < posterior.size(); ++i) {
        if (i==0) {output_file << posterior[i] << "," << P_joint[i] << "," << dt_cpp << endl;}
        else {output_file << posterior[i] << "," << P_joint[i] << endl;}
    }
    output_file.close();
    cout << "Saved 'posterior', 'P_joint' and 'dt' to CSV files" << endl;

    return 0;
}

