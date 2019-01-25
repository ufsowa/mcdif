#include <stdio.h>
#include <vector>
#include <fstream>
#include <iostream>

using namespace std;

ifstream clusters_file;
ifstream ECI_file;

int main() {

clusters_file.open("clusters.out");
ECI_file.open("eci.out");

if (!clusters_file) {
    cerr << "Unable to open file clusters.out";
    exit(1);   // call system to stop
}

if (!ECI_file) {
    cerr << "Unable to open file eci.out";
    exit(1);   // call system to stop
}

// while (clusters_file >> x) {
//    sum = sum + x;
//}

int x;
while (ECI_file >> x) {
    cout << x;
}
return 0;
}