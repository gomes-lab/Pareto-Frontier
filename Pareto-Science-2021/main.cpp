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
bool use_n_squared;
std::string file_path;
std::string river_newtwork;
std::vector<std::string> relevant_criteria;
unsigned int seed;
int batch_size;
size_t num_threads;
bool use_binary_tree;
bool split_values;
bool balance_tree;

bool parse_arguments(int argc, char **argv) {
    for (int i = 1; i < argc; i++) {
        if (std::string(argv[i]) == "-criteria") {
            int num_criteria = stoi(std::string(argv[++i]), 0); // Number of criteria to read in
            relevant_criteria.clear();
            for (int num = 0; num < num_criteria; num++) {
                relevant_criteria.emplace_back(argv[++i]);
            }
        } else if (std::string(argv[i]) == "-basin") {
            river_newtwork = std::string(argv[++i]);
        } else if (std::string(argv[i]) == "-path") {
            file_path = std::string(argv[++i]);
        } else if (std::string(argv[i]) == "-seed") {
            seed = static_cast<unsigned int>(stoul(argv[++i], nullptr, 0));
        } else if (std::string(argv[i]) == "-epsilon") {
            epsilon = stod(argv[++i], nullptr);
        } else if (std::string(argv[i]) == "-nlogn") {
            use_n_squared = false;
        } else if (std::string(argv[i]) == "-nsquare") {
            use_n_squared = true;
        } else if (std::string(argv[i]) == "-batch"){
            batch_size = stoi(argv[++i], nullptr, 0);
        } else if (std::string(argv[i]) == "-binary") {
            use_binary_tree = true;
        } else if (std::string(argv[i]) == "-directed") {
            use_binary_tree = false;
        } else if (std::string(argv[i]) == "-split") {
            split_values = true;
        } else if (std::string(argv[i]) == "-nsplit") {
            split_values = false;
        } else if (std::string(argv[i]) == "-balance") {
            balance_tree = true;
        } else if (std::string(argv[i]) == "-chain_r") {
            balance_tree = false;
        } else if (std::string(argv[i]) == "-thread") {
            num_threads = stoi(argv[++i], nullptr);
        } else {
            std::cout << "General Parameters" << std::endl;
            std::cout << "-criteria N STR STR = reads in N relevant criteria" << std::endl;
            std::cout << "-basin STR = the name of the river basin that we are working on" << std::endl;
            std::cout << "-path FILE = the file with the data for the river basin" << std::endl;
            std::cout << "-seed N = the random seed to use in the experiment (NOTE: if this is not given a random seed is chosen)" << std::endl;
            std::cout << "-epsilon N = the epsilon that will be used for rounding in the experiment" << std::endl;
            std::cout << "-base N = the base used for calculation of seismic energy risk (i.e. energy * base ^ risk)" <<std::endl;
            std::cout << "-nlogn = Use the divide and conquer algorithm for determining the non-dominated solutions for each node" << std::endl;
            std::cout << "-nsquare = Use the n^2 algorithm for determining the non-dominated solutions for each node" << std::endl;
            std::cout << "-batch N = The batch size used for dynamically processing partial policies in nlogn algorithm" << std::endl;
            std::cout << "-binary = Use binary tree structure" << std::endl;
            std::cout << "-directed = Use original directed hypernode graph" << std::endl;
            std::cout << "-split = split binary tree values between intermediates" << std::endl;
            std::cout << "-nsplit = assign local values of zero for all intermediates" << std::endl;
            std::cout << "-balance = construct a balanced binary tree equivalent" << std::endl;
            std::cout << "-chain_r = assign all intermediates to the right children in binary EQ" << std::endl;
            exit(-1);
        }
    }
    return argc > 1;
}


int main (int argc, char **argv) {
    // Set the standard values for parameters -- for when run in local environment
    epsilon =0.05;

    seed = (unsigned int) time(nullptr);
    use_n_squared = false;
    use_binary_tree = true;
    river_newtwork = "AB";
    batch_size = 10000000;
    split_values = true;
    balance_tree = false;
    file_path="AB_cppinput_2021_Feb.txt";
    relevant_criteria = {"energy","connectivity"};

    num_threads=6;

    parse_arguments(argc, argv);
    // Set the seed!
    srand(seed);

    HyperNet net(river_newtwork, file_path, split_values, balance_tree, epsilon);

    DP_Algorithm dp(net, net.root, relevant_criteria, epsilon, seed, use_n_squared, batch_size, use_binary_tree, num_threads);



    bool to_file = false;
    bool to_graph = false;
    bool experiment = true;

    vector<double> epsilons = {0, 0.001, 0.01, 0.025, 0.05, 0.01};



    if (experiment) {
        dp.run_experiment();
    } else {
        // Get time
        clock_t t;
        t = clock();
        //dp.build_DP_table_recursive_tree(net.root_node, true);
        dp.build_DP(true);
        //dp.build_DP_table_recursive(net.root);
        t = clock() - t;


        if (to_file) {
            // Output to file
            ofstream output("results.txt");
            if (output.is_open()) {
                output << "It took me " << t << "clicks (" << ((float) t) / CLOCKS_PER_SEC << "seconds).\n";
                dp.print_DP_File(output);
                //dp.print_DP_Output(output);
                output.close();
            }
        } else {
            printf("It took me %lu clicks (%f seconds).\n", t, ((float) t) / CLOCKS_PER_SEC);
            dp.print_DP_Output();
        }
    }

    return 0;
}