//
// Created by Jonathan Gomes Selman on 7/6/17.
//

#ifndef AMAZON_PROJECT_DP_ALGORITHM_H
#define AMAZON_PROJECT_DP_ALGORITHM_H

#include <stdlib.h>
#include "HyperNet.h"
#include "Pareto_Solution.h"
#include <stack>
#include <cmath>
#include <thread>
#include <algorithm>
#include <mutex>
#include <iomanip>
#include <unordered_set>
#include "ThreadPool.h"

#define MAX_CRITERIA 2

/*
 * Compares solutions based on a given dimension.
 * Ties are broken by comparing solutions lexicographically.
 * Used to sort partial policies
 */
struct Compare_dimensions {
    explicit Compare_dimensions(int dimension, int total_dimensions) {
        this->dimension = dimension;
        this->total_dimensions = total_dimensions;
    }
    bool operator () (Pareto_Solution* sol1, Pareto_Solution* sol2) {
        // Because arrays are 0 based index we subtract from the dimension
        for (int i = dimension - 1; i >= 0; i--) {
            // If the solutions are equal in the given dimension, try next dimension
            if (sol1->criteria_2[i].second == sol2->criteria_2[i].second) {
                continue;
            }
            return sol1->criteria_2[i].second > sol2->criteria_2[i].second; // Sort in descending order
        }
        // Resolve ties in all lower dimensions by compare iteratively across the higher dimensions.
        // Note: we will never have a tie in all of the dimensions (pre-condition to the algorithm)
        int tie = dimension;
        while (tie <= total_dimensions - 1) {
            if (sol1->criteria_2[tie].second == sol2->criteria_2[tie].second) {
                tie++;
            } else {
                return sol1->criteria_2[tie].second > sol2->criteria_2[tie].second;
            }
        }
        return false;
    }

    int dimension;
    int total_dimensions;
};

struct EqualSolutions{
    bool operator()(Pareto_Solution const* a, Pareto_Solution const* b) const {
        return *a == *b;
    }
};

struct SolutionHash {
    std::size_t operator()(const Pareto_Solution* solution) const {
        // Compute individual hash values each criterion's rounded value
        // Hash approach courtesy of ---
        // http://stackoverflow.com/a/1646913/126995
        size_t res = 17;
        for (int i = 0; i < solution->getNum_criteria(); i++) {
            res = res * 31 + std::hash<double>{}(solution->criteria_2[i].second);
        }
        return res;
    }
};


class DP_Algorithm {
public:

    /*
     * Construct a DP_Algorithm object from which you can run experiments on specified
     * river networks.
     */
    DP_Algorithm(const HyperNet &net, int root, bool maximize_c1, bool maximize_c2, double epsilon,
                 unsigned int seed, int batch_size, size_t num_threads);


    /*
     * Calculate the k values for each node.
     * -------------------------------------
     * k values are calculated by taking the local value for each criteria and
     * multiplying by the constant K_EPSILON.
     *      k_c1 = K_EPSILON * associated_c1
     */
    void calc_theoretical_ks();

    void calc_theoretical_ks_tree(HyperTreeNode* node);

    /*
     * Given a value, return the rounded value.
     * ----------------------------------------
     * Return: floor(value / k) * k
     */
    double round_with_k(double value, double k);

    /* -------------------------------------------------------------------------------------------------------------
     * -------------------------------------------------------------------------------------------------------------
     *
     * Bellow are different methods for performing the dynamic programming algorithm
     * to compute the Pareto-Optimal solutions for a given node in the Hyper-tree.
     * ------------------------------------------------------------------------
     * In each algorithm, we consider all the possible decisions for building dams
     * between the node's children. Then For each possible build decision,
     * we analyze the resulting partial policy created by the combination of Pareto_Optimal
     * solutions from the children nodes. At each node we only keep the resulting
     * Pareto_Optimal partial policies for the given parent node.
     *
     * -------------------------------------------------------------------------------------------------------------
     * -------------------------------------------------------------------------------------------------------------
     */


    /*
     * Implements the dynamic programming algorithm with an nlogn method
     * for determining non-dominated solutions. Additionally, batched pruning
     * is used to allow for parallelization and prevent memory overflows.
     * Specifically, rather than pruning all candidate solutions at once for a given
     * node, we prune the candidate solutions in batches of size specified by the user.
     */
    void computeDP_nlogn_batching(HyperTreeNode *node);

    /*
     * Helper method used by computerDP_nlogn_batching to add a new partial policy
     */
    void add_nlogn_solution(std::vector<Pareto_Solution *> &grouping, std::vector<int> &decisionVec,
                                HyperTreeNode *node,
                                std::unordered_set<Pareto_Solution *, SolutionHash, EqualSolutions> &solution_set,
                                int &num_generated, int &batch_num);

    /*
     * Nlogn Algorithm for determining the non-dominated set (i.e. set of Pareto
     * Optimal solutions) S', given a set of potential solutions S with 2 criteria.
     * ----------------------------------------------------------------------------
     * Algorithm:
     *     1) Sort the candidate solutions by the first dimension in descending order
     *     2) Iterate through the solutions, keeping track of the maximum value in
     *        the second dimension that you have seen.
     *        2a) Keep solutions that have a second dimension greater than the maximum
     *            seen so far.
     *        2b) Prune solutions otherwise
     * ----------------------------------------------------------------------------
     * Solutions returned are gauranteed to by Pareto Optimal solutions w/r to the set S
     */
    std::vector<Pareto_Solution* > L2D(std::vector<Pareto_Solution* >& all_possible);


    /*
     * Initialize the criteria array for each new partial policy
     * with the local values for node 'node'.
     */
    void initialize_criteria(HyperTreeNode* node, std::pair<double, double> *criteria);


    /*
     * Updates the criteria array by adding the appropriate values
     * from a specific child of node 'node'. The variable 'pareto_decision'
     * represents the partial policy being considered from the given child
     * at 'child_index'. 'decision' indicates whether the dam connecting
     * the child to 'node' is being built or not in the current partial policy.
     *
     * Note: This version is used by 'computeDP_n_squared_generic and add_nlogn_solution'
     */
    void update_criteria_new(HyperTreeNode *node, Pareto_Solution *pareto_decision, std::pair<double, double> *criteria,
                                 int decision, int child_index, double &c1_built, double &c2_built);

    /*
     * Assign the rounding constant k for a given node based on the dams that are built at that node.
     */
    void assign_node_k(HyperTreeNode *node, double c1_built, double c2_built);

    /*
     * After a new partial policy is created at a given 'node',
     * each criteria value is rounded using its
     * corresponding rounding constant k.
     *
     * Note: This version is used by 'computeDP_n_squared_generic and add_nlogn_solution'
     */
    void round_criteria(HyperTreeNode* node, std::pair<double, double> *criteria);

    /*
     * Recursively traverse the tree (non-binary) with a post-order traversal calculating
     * the Pareto_Opt_List for each parent node.
     * -----------------------------------------
     * 1) Recurse on children
     * 2) Compute DP for parent
     */
    void build_DP(bool round_values);

    /*
     * Recursively traverse the binary-tree with a post-order traversal calculating
     * the Pareto_Opt_List for each parent node.
     * -----------------------------------------
     * 1) Recurse on left/right children
     * 2) Compute DP for parent
     */
    void build_DP_table_recursive_tree(HyperTreeNode* root, bool round_values);

    // --------------------------------------------------------------------------------


    // --------------------------------------------------------------------------------
    // Functions to print output of experiments
    // --------------------------------------------------------------------------------

    void print_DP_Output(std::ostream &stream = std::cout); // If no output stream is provided then prints  to console

    void print_dams_generic_tree_vec(HyperTreeNode* node, std::ostream &out = std::cout);

    void print_dams_helper(Pareto_Solution *result, std::ostream &out, HyperTreeNode *node);

    // --------------------------------------------------------------------------------

    /*
     * Runs an experiment for a given epsilon and river basin.
     * Outputs all results to a file with a standardized format.
     * Keeps track of several experiment metrics such as the execution
     * time and time spent doing particular tasks.
     */
    void run_experiment();

private:
    HyperNet net; // Hyper node network with the all the data
    int root; // The id of the root node in the hyper net

    // Used to track how much progress
    double node_processed_counter;

    bool maximize_c1;
    bool maximize_c2;

    // Used to calculate the theoretical k values
    // See the function 'calc_theoretical_ks' to see how K_EPSILON is used
    double k_epsilon;

    unsigned int random_seed; // Keeps track of the random seed used for the last DP run

    int batch_size; // Used to dynamically process the possible Pareto-Solutions in effort to limit memory

    // Used for experimentation on algorithm complexity
    // Do we need these???
    //double time_sorting;
    //double time_generating;
    //double time_nlogn;
    //double time_copying;

    // Used for comparing binary vs. non-binary
    unsigned long max_policies_considered;
    unsigned long total_policies_considered;
    unsigned long total_pruned_policies;

    // Fields used for parallel batch
    int K;
    std::mutex mtx;
    std::vector<std::unordered_set<Pareto_Solution *, SolutionHash, EqualSolutions> >  solution_sets;
    ThreadPool pool;
    std::vector< std::future<int> > results;

    int batch_number;

};


#endif //AMAZON_PROJECT_DP_ALGORITHM_H
