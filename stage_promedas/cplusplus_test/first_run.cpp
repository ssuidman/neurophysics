#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <chrono>

using namespace std;

int main() {
    // Retrieve variables from CSV files
    vector<double> previn;
    vector<vector<double> > pfmin;
    vector<double> pfminneg;
    vector<long double> posterior_;

    ifstream previnFile("patient404_case_1_previn.csv");
    ifstream pfminFile("patient404_case_1_pfmin.csv");
    ifstream pfminnegFile("patient404_case_1_pfminneg.csv");
    ifstream posterior_File("patient404_case_1_posterior.csv");

    string line;

    while (getline(previnFile, line, '\n')) {
        previn.push_back(stold(line));
    }

    while (getline(pfminFile, line)) {
        vector<double> row;
        istringstream ss(line); // Use a stringstream to split the line
        string cell;
        while (getline(ss, cell, ',')) {
            row.push_back(stod(cell));
        }
        pfmin.push_back(row);
    }


    while (getline(pfminnegFile, line, '\n')) {
        pfminneg.push_back(stod(line));
    }

    while (getline(posterior_File, line, '\n')) {
        posterior_.push_back(stold(line));
    }

    int m = (!pfmin.empty()) ? pfmin.size() : 0;
    int n = previn.size();

    vector<vector<double> > prev(n + 1, previn);
    for (int i = 0; i < n; i++) {
        prev[i + 1][i] = 1.0;
    }

    vector<double> prevend(n + 1, 0.0);
    prevend[0] = 1.0;
    for (int i = 1; i <= n; i++) {
        prevend[i] = previn[i - 1];
    }

    vector<vector<double> > prevminneg(n + 1, vector<double>(n, 0.0));
    for (int i = 0; i <= n; i++) {
        for (int j = 0; j <= n; j++) {
            prevminneg[i][j] = (!pfminneg.empty()) ? (prev[i][j] * pfminneg[j]) : (prev[i][j]);
        }
    }

    vector<double> pfplus(n + 1, 0.0);
    vector<int> myset;
    double t = 0.0;

    // for (const vector<double>& row : prevminneg) {for (double value : row) {cout << value << " ";} cout << "\n\n" << endl;}    
    // for (const double& value : prevend) {std::cout << value << " ";} cout << endl;
    cout << "previn:\t\t" << previn.size() << endl;
    cout << "pfmin:\t\t(" << pfmin.size() << "," << pfmin[1].size() << ")" << endl;
    cout << "pfminneg:\t" << pfminneg.size() << endl;
    cout << "posterior_:\t" << posterior_.size() << endl;
    cout << "prev:\t\t(" << prev.size() << "," << prev[1].size() << ")" << endl;
    cout << "prevend:\t" << prevend.size() << endl;
    cout << "prevminneg:\t(" << prev.size() << "," << prev[1].size() << ")" << endl;
    cout << "################### I think variables are okay up till here ###################" << endl;

    auto start = chrono::steady_clock::now();

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
        // cout << "v = \t\t[" << v[0] << "," << v[1] << "," << v[2] << "," << v[3] << "," << v[4] << "," << v[5] << "," << v[6] << "]" << endl;
        // cout << "myset = \t["; for (int i = 0; i < myset.size(); i++) {cout << myset[i]; if (i < myset.size() - 1) {cout << ",";}} cout << "]." << endl;
        
        // HIER VERDER!!!!!!!!!
        double sumExp = 0.0;
        for (int j = 0; j <= n; j++) {
            double product = 1.0;
            for (int k : myset) {
                product *= pfmin[k][j];
            }
            sumExp += log(1e-50 + (product * prevminneg[j][0] + (1 - prevend[j])));
        }
        pfplus[i] += pow(-1, myset.size()) * exp(sumExp);
    }

    auto end = chrono::steady_clock::now();
    t = chrono::duration_cast<chrono::milliseconds>(end - start).count() / 1000.0;

    vector<double> P_joint(n, 0.0);
    for (int i = 1; i <= n; i++) {
        P_joint[i - 1] = pfplus[i] * previn[i - 1];
    }

    // for (const double& value : pfplus) {std::cout << value << " ";} cout << endl;
    cout << "Time: " << t << endl;
    // double posterior = P_joint[0] / pfplus[0];
    
    // double maxDifference = 0.0;
    // for (size_t i = 0; i < posterior_.size(); i++) {
    //     double diff = std::abs(static_cast<double>(posterior_[i]) - posterior);
    //     maxDifference = std::max(maxDifference, diff);
    // }
    // cout << "Maximum difference: " << maxDifference << endl;


    return 0;
}
