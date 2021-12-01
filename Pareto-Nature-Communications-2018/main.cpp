#include <iostream>
#include <vector>
#include <algorithm>
#include "HyperNet.h"
#include "DP_Algorithm.h"
#include "time.h"
#include "math.h"

using namespace std;

// The arguments that are relevant in main
double epsilon;
std::string file_path;
std::string river_network;
bool maximize_c1;
bool maximize_c2;
unsigned int seed;
int batch_size;
size_t num_threads;

bool parse_arguments(int argc, char **argv) {
    for (int i = 1; i < argc; i++) {
        std::cout<<argv[i]<<std::endl;
        if (std::string(argv[i]) == "-criteria") {
            // Expect two max/min flags, one for each criteria
            if (std::string(argv[i+1]).find("max")==0) {
                maximize_c1 = true;
            } else if (std::string(argv[i+1]).find("min")==0){
                maximize_c1 = false;
            } else {
                std::cout << "-criteria max/min max/min = reads in whether to maximize or minimize each criteria" << std::endl;
                exit(-1);
            }
            i++;
            if (std::string(argv[i+1]).find("max")==0) {
                maximize_c2 = true;
            } else if (std::string(argv[i+1]).find("min")==0){
                maximize_c2 = false;
            } else {
                std::cout << "-criteria max/min max/min = reads in whether to maximize or minimize each criteria" << std::endl;
                exit(-1);
            }
            i++;

        } else if (std::string(argv[i]) == "-basin") {
            river_network = std::string(argv[++i]);
        } else if (std::string(argv[i]) == "-path") {
            file_path = std::string(argv[++i]);
        } else if (std::string(argv[i]) == "-seed") {
            seed = static_cast<unsigned int>(stoul(argv[++i], nullptr, 0));
        } else if (std::string(argv[i]) == "-epsilon") {
            epsilon = stod(argv[++i], nullptr);
        } else if (std::string(argv[i]) == "-batch"){
            batch_size = stoi(argv[++i], nullptr, 0);
        } else if (std::string(argv[i]) == "-thread") {
            num_threads = stoi(argv[++i], nullptr);
        } else {
            std::cout << "General Parameters" << std::endl;
            std::cout << "-criteria max/min max/min = reads in whether to maximize or minimize each criteria" << std::endl;
            std::cout << "-basin STR = the name of the river basin that we are working on" << std::endl;
            std::cout << "-path FILE = the file with the data for the river basin" << std::endl;
            std::cout << "-seed N = the random seed to use in the experiment (NOTE: if this is not given a random seed is chosen)" << std::endl;
            std::cout << "-epsilon N = the epsilon that will be used for rounding in the experiment" << std::endl;
            std::cout << "-batch N = The batch size used for dynamically processing partial policies in nlogn algorithm" << std::endl;
            exit(-1);
        }
    }
    return argc > 1;
}

int main (int argc, char **argv) {
    // Set the standard values for parameters -- for when run in local environment
    epsilon = 0.05;

    seed = (unsigned int) time(nullptr);
    river_network = "Amazon_May";
    batch_size = 100000;
    file_path="AB_energy_ghg20.txt";
    maximize_c1 = true;
    maximize_c2 = false;
    num_threads= 6;
    parse_arguments(argc, argv);

    // Set the seed!
    srand(seed);

    HyperNet net(river_network, file_path, epsilon, maximize_c1, maximize_c2);

    DP_Algorithm dp(net, net.root, maximize_c1, maximize_c2, epsilon, seed, batch_size, num_threads);

    dp.run_experiment();

    return 0;
}