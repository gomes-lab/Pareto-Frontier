//
// Created by Jonathan Gomes Selman on 7/6/17.
//

#ifndef AMAZON_PROJECT_DP_ALGORITHM_H
#define AMAZON_PROJECT_DP_ALGORITHM_H

#include <stdlib.h>
#include "HyperNet.h"
#include "Pareto_Opt_List.h"
#include <stack>
#include <cmath>
#include <thread>
#include <algorithm>
#include <mutex>
#include <iomanip>
#include <unordered_set>
#include "ThreadPool.h"

#define MAX_CRITERIA 8

/*
 * Used to compare solutions based on a given dimension.
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
        // We will never have solutions equal in all dimensions so now we compare across all dimensions
        // Resolve tie in all lower dimensions by comparing upper dimensions
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
    static constexpr double EQUALITY_EPSILON = 0.00001;
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
    DP_Algorithm(const HyperNet &net, int root, const std::vector<std::string> &relevant_criteria, double epsilon,
                 unsigned int seed, bool use_n_squared, int batch_size, bool use_binary_tree, size_t num_threads); // Constructor
    // Constructor


    /*
     * Calculate the k values for each node.
     * -------------------------------------
     * The k values are calculated by taking the local value for each criteria and
     * multiplying by the constant K_EPSILON.
     *      k_energy = K_EPSILON * associated_energy
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
     *
     * This is the most basic version of the dynamic programming algorithm. This implementation
     * uses the orignal, non-binary tree and an n_squared algorithm for determining the
     * Pareto-optimal partial policies at each node.
     *
     * Note: this version has not been updated to include the new rounding scheme for
     * energy and seismic (i.e. where we determine rounding constant k at each node
     * by considering the dams that are actually built)!
     */
    void computeDP_n_squared(int node);

    /*
     * Initialize the criteria array for each new partial policy
     * with the local values for node 'node'. For example, the
     * array location representing connectivity is initialized
     * to the local connectivity of the node and energy is initialized
     * to 0.
     *
     * Note: This version is used by 'computeDP_n_squared'
     */
    void initialize_criteria(int node, std::pair<double, double> *criteria);

    /*
     * Updates the criteria array by adding the appropriate values
     * from a specific child of node 'node'. The variable 'pareto_decision'
     * represents the partial policy being considered from the given child
     * at 'child_index'. 'decision' indicates whether the dam connecting
     * the child to 'node' is being built or not in the current partial policy.
     *
     * Note: This version is used by 'computeDP_n_squared'
     */
    void update_criteria(int node, Pareto_Solution *pareto_decision, std::pair<double, double> *criteria,
                         int decision, int child_index);

    /*
     * After a new partial policy is created at a given 'node',
     * each criteria value is rounded using its
     * corresponding rounding constant k.
     *
     * Note: This version is used by 'computeDP_n_squared'
     */
    void round_criteria(int node, std::pair<double, double> *criteria);

    /*
     * Implements the dynamic programming algorithm with a divide-and-conquer method
     * for determining non-dominated solutions as well as batched pruning.
     */
    void computeDP_nlogn_batching(HyperTreeNode *node);

    /*
     * Helper method used by computerDP_nlogn_batching to add a new partial policy
     */
    void add_nlogn_solution(std::vector<Pareto_Solution *> &grouping, std::vector<int> &decisionVec,
                                HyperTreeNode *node,
                                std::unordered_set<Pareto_Solution *, SolutionHash, EqualSolutions> &solution_set,
                                int &num_generated, clock_t &t, int &batch_num);



    /*
     * Implementation of dynamic programming algorithm that uses a binary tree
     * representation of the network and the naive n_squared approach for determining
     * non-dominated solutions
     */
    void computeDP_n_squared_generic(HyperTreeNode* node, bool round_values);

    /*
     * Initialize the criteria array for each new partial policy
     * with the local values for node 'node'. For example, the
     * array location representing connectivity is initialized
     * to the local connectivity of the node and energy is initialized
     * to 0.
     *
     * Note: This version is used by 'computeDP_n_squared_generic' and add_nlogn_solution
     */
    void initialize_criteria(HyperTreeNode* node, std::pair<double, double> *criteria);

    /*
     * Updates the criteria array by adding the appropriate values
     * from a specific child of node 'node'. The variable 'pareto_decision'
     * represents the partial policy being considered from the given child
     * at 'child_index'. 'decision' indicates whether the dam connecting
     * the child to 'node' is being built or not in the current partial policy.
     *
     * Note: This is an old version used before new rounding scheme
     */
    void update_criteria(HyperTreeNode* node, Pareto_Solution *pareto_decision, std::pair<double, double> *criteria,
                            int decision, int child_index);

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
                                 int decision, int child_index, double &energy_built, double &seismic_built,
                                 double &biodiversity_built, double &ghg_built, double &dor_built, double &population_built);

    /*
     * Used in the new rounding scheme to assign the rounding constant k for a given
     * node based on the dams that are built at that node. See function for the
     * full formula for calculating a nodes seismic and energy rounding constants.
     */
    void assign_node_k(HyperTreeNode *node, double energy_built, double seismic_built, double biodiversity_built,
                                double ghg_built, double dor_built, double population_built);

    /*
     * After a new partial policy is created at a given 'node',
     * each criteria value is rounded using its
     * corresponding rounding constant k.
     *
     * Note: This version is used by 'computeDP_n_squared_generic and add_nlogn_solution'
     */
    void round_criteria(HyperTreeNode* node, std::pair<double, double> *criteria);

    /*
     * Recursively traverse the tree with a post-order traversal calculating
     * the Pareto_Opt_List for each parent node.
     * -----------------------------------------
     * 1) Recurse on children
     * 2) Compute DP for parent
     */
    void build_DP_table_recursive(int root); // Probably want to save the root!!

    void build_DP(bool round_values);

    /*
     * Recursively traverse the binary-tree with a post-order traversl calculating
     * the Pareto_Opt_List for each parent node.
     * -----------------------------------------
     * 1) Recurse on left/right children
     * 2) Compute DP for parent
     */
    void build_DP_table_recursive_tree(HyperTreeNode* root, bool round_values);

    // Below are verious methods to print out information based on each experimental run
    // ---------------------------------------------------------------------------------
    void print_DP_File(std::ofstream &file);

    void print_DP_Output(std::ostream &stream = std::cout); // If no output stream is provided then prints  to console

    void print_dams_non_binary_list(int node, std::ostream &out = std::cout);

    void print_dams_generic_tree_list(HyperTreeNode *node, std::ostream &out = std::cout);

    void print_dams_generic_tree_vec(HyperTreeNode* node, std::ostream &out = std::cout);

    void print_dams_helper(Pareto_Solution *result, std::ostream &out, HyperTreeNode *node);

    void print_dams_helper(Pareto_Solution *result, std::ostream &out, int node);
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



    //static const bool use_n_squared = false;
    bool use_n_squared;
    bool use_binary_tree;

    // Used to calculate the theoretical k values
    // See the function 'calc_theoretical_ks' to see how K_EPSILON is used
    double k_epsilon;

    std::vector<std::string> relevant_criteria;

    unsigned int random_seed; // Keeps track of the random seed used for the last DP run

    int batch_size; // Used to dynamically process the possible Pareto-Solutions in effort to limit memory

    // Used for experimentation on algorithm complexity
    double time_sorting;
    double time_generating;
    double time_nlogn;
    double time_copying;

    // Used for comparing binary vs. non-binary
    unsigned long max_policies_considered;
    unsigned long total_policies_considered;
    unsigned long total_pruned_policies;

    // Maps: node id --> Pareto_Opt_List
    // Represents the list of Pareto optimal pairs for a given hyper node
    std::unordered_map<int, Pareto_Opt_List> pareto_opt_table;

    // Used for non-binary
    std::unordered_map<int, Frontier_List> pareto_partial_table;

    // Likely we should just add this to the parameters of each node!!!!!!!!!!!!
    std::unordered_map<int, double> energies_below; // Gives max energy for each node --- if all dams before were built
    std::unordered_map<int, double> min_energy; // Min energy -- don't build any new dams
    std::unordered_map<int, double> max_conn; // Max conn -- no dams * min con is always just the con of the node

    /*
     * Method to take the cartesian product of multiple sets.
     * Returns a vector of vectors, where each internal vector represents
     * a selection of objects, each from a different initial set.
     * Uses the algorithm for the python function itertools.product.
     * Does not store intermediate values when calculating the cartesian
     * product for more than 2 sets.
     */
    template <typename T>
    std::vector<std::vector<T> > cartesian_product_generic(const std::vector<std::vector<T> >& table);


    struct marked_elements{
        Pareto_Opt_Node* node;
        bool inferior;
    };


    // See the CPAIOR paper for in-depth details and pseudo-code describing the
    // implementation of these functions, used for determining the non-dominated
    // solutions in a given set of solutions.
    // -------------------------------------------------------------------------
    // Calls divide_and_conquer helper or L2D to determine the Pareto-Optimal solutions
    std::vector<Pareto_Solution* > divide_and_conquer(std::vector<Pareto_Solution* >& curr_set, int dimension);

    std::vector<Pareto_Solution* > divide_and_conquer_helper(std::vector<Pareto_Solution *> &curr_set,
                                                             unsigned long low, unsigned long high, int dimension);

    /*
     * @precondition - The vector curr_set is sorted by dimension
     */
    std::vector<Pareto_Solution* > marry(std::vector<Pareto_Solution *> &curr_set, unsigned long low,
                                         unsigned long high, int dimension);

    /*
     * @precondition - The vector curr_set is sorted by the second to last dimension
     */
    std::vector<Pareto_Solution* > marry_2d(std::vector<Pareto_Solution* >& curr_set);

    std::vector<Pareto_Solution* > L2D(std::vector<Pareto_Solution* >& all_possible);
    // -------------------------------------------------------------------------

    // Fields used for parallel batch
    int K;
    std::mutex mtx;
    std::vector<std::unordered_set<Pareto_Solution *, SolutionHash, EqualSolutions> >  solution_sets;
    ThreadPool pool;
    std::vector< std::future<int> > results;

    int batch_number;

};


#endif //AMAZON_PROJECT_DP_ALGORITHM_H
