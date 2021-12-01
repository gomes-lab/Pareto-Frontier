//
// Created by Jonathan Gomes Selman on 7/6/17.
//

#ifndef AMAZON_PROJECT_HYPERNET_H
#define AMAZON_PROJECT_HYPERNET_H

#include <vector>
#include <unordered_map>
#include "string"
#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>
#include "Pareto_Solution.h"

// Represents an edge in the tree
// Used by the adjacency list implementation
struct Edge {
    int endNode;
    int damId;
};

// Represents a Hyper Node in the tree
struct HyperNode {
    int id;

    // Individual k values used for rounding
    // The value of k represents the "bucket size" used for rounding
    //      If k = 5 --> the number 7 gets rounded down to 5
    double k_c1;
    double k_c2;
    double c1;
    double c2;
};

// Represents a Dam in the Tree
struct Dam{
    double c1;
    double c2;
    int status;
    int id;

    // Holds the possible decisions for the dam
    //      {0,1} = planned dam that we are analyzing
    //      {1} = dam is already built
    //      {0} = this is a dummy dam that exists in the binary-tree representation of the graph and is never built!
    std::vector<int> decision;
};


// Represents a Node in the binary-tree representation of the network
struct HyperTreeNode{
    HyperNode node_data; // The actual data for the tree node

    // Children
    HyperTreeNode* left;
    HyperTreeNode* right;

    // Each child is connected by an specific dam
    Dam leftD;
    Dam rightD;

    bool intermediate; // Indicates whether the node is an intermediate node --- created to make the tree binary

    std::vector<Pareto_Solution* > vec_frontier; // Used for nlogn algorithm and for generic experimentation



    HyperTreeNode(const HyperNode& node, bool is_intermediate) {
        node_data = node;
        intermediate = is_intermediate;

        left = nullptr;
        right = nullptr;
    }

    ~HyperTreeNode() {
        for (int i = 0; i < vec_frontier.size(); i++) {
            delete(vec_frontier[i]);
        }
    }

    HyperTreeNode(const HyperTreeNode& src) {
        node_data = src.node_data;
        left = src.left;
        right = src.right;
        leftD = src.leftD;
        rightD = src.rightD;
        intermediate = src.intermediate;
        vec_frontier = src.vec_frontier;
    }

};



class HyperNet {
public:

    /*
     * Constructor used by the most current version
     */
    HyperNet(const std::string &basin, const std::string &file, double epsilon, bool maximize_c1, bool maximize_c2);

    virtual ~HyperNet(); // Destructor

    /*
    * Follow the rule of 3:
    *      Destructor
    *      Copy constructor
    *      Assignment operator
    */
    void destructor_helper(HyperTreeNode* node);
    HyperNet(const HyperNet& src); // Copy constructor
    HyperNet& operator=(const HyperNet& src); // Assignment operator

    /*
     * Read the adjacent list from the input file
     */
    void readAdjList(std::ifstream& input, int num_edges);

    /*
     * Read info for each dam
     */
    void readDamInfo(std::ifstream& input, int num_dams);

    /*
     * Given a Hyper-Tree, transform the tree into a binary tree.
     * Uses helper methods to either create a balanced or unbalanced tree.
     */
    void generateBinaryTree();

    /*
     * Print the Hyper-Tree in-order
     */
    void print_adjacent_tree(int node);

    /*
     * Print the binary Hyper-Tree in-order
     */
    void print_binary(HyperTreeNode* root);

    /*
     * Calculates the rounding constants for the
     * preliminary rounding of criteria
     * ------------------------------------------------
     * These rounding constants are determined by the minimal
     * criteria values for criteria.
     * -------------------------------------------
     * The formula is:
     *      epsilon * minimum_criteria_val / 2
     */
    void calculate_k();


    /*
     * Round all of the dam criteria
     * based on the respective rounding constants. This
     * is the preliminary rounding phase before the
     * dynamic programming algorithm is run.
     */
    void round_criteria(int node);

    /*
     * Helper method to round specific values:
     *      floor ( c1 / k_factor) * k_factor
     * Note: If the k_factor is zero, the value of
     * c1 is returned.
     */
    double round_kp(double cr, double k_factor);



    // Data structures for different attributes of Hyper Tree
    std::unordered_map<int, std::vector<Edge> > adj_lists; // Adjacency list of edges for each node
    std::vector<HyperNode> nodes; // Info about each node
    std::unordered_map<int, Dam> dams; // Info for dams

    // Basic data about the tree
    int num_nodes;
    int num_edges;
    int num_dams;
    bool maximize_c1;
    bool maximize_c2;
    int root; // Defines the root of the tree

    std::vector<std::string> relevant_criteria; // Stores the criteria that we will use for the analysis

    HyperTreeNode* root_node;
    // Count for number of nodes in the binary tree
    int num_nodes_binary;

    // Data file path
    std::string filename;

    // River network name. Example: Maranon
    std::string river_network;

    // minimum value of criteria 1
    double min_c1;
    // minimum value of criteria 1
    double min_c2;

    // Epsilon accuracy
    double epsilon;
    // The rounding constants
    double k_c1;
    double k_c2;
    // The names of criteria
    std::string c1_name;
    std::string c2_name;


private:

    /*
     * Calculates the number of intermediate nodes that will be needed to represent
     * a parent node with a given number of children in a binary tree. The number
     * of intermediate nodes is num_children - 2.
     */
    int calculate_number_intermediate_nodes(int num_children);

    // Helper method to perform the operations of the defualt constructor
    void construct();

    /*
     * Recursively generates a binary tree representation of the Hyper-Tree.
     * Algorithm:
     *      1) Make the current node's left child a new node with the node_data
     *      equal to the child at children_index within the original tree.
     *      2) Recurse on the left child --- children_index = 0
     *      3) Check how many children are left to be added to the tree for the given node
     *          - If there are more than one to add we must an intermediate:
     *          Make the node's right child a new intermediate node with the same data as itself
     *          - If there is only one child to add:
     *          Make the node's right child a new node with the node_data equal to
     *          the child at children_index + 1 within the original tree
     *      4) Recurse on the right child
     *          - If right child is intermediate, children_index += 2
     *          - Else children_index = 0
     *
     *      *** NOTE ***
     *      Whenever a non-intermediate node is added, calculate the number of intermediates
     *      that are needed to represent the node. The connectivity of that node will then
     *      be split between it and its intermediates (i.e. conn = node.conn / (num_inter + 1) )
     */
    void generate_binary_tree_helper(HyperTreeNode* node, int children_index);

    /*
     * Used to copy the binary-tree into a new HyperNet object
     */
    void copy_constructor_helper(HyperTreeNode*& newTree, HyperTreeNode* oldTree);

    /*
     * Updates the criteria of the HyperNode to reflect the split
     * that is being made when introducing intermediary nodes.
     * For example, if the number of intermediates is 5, we
     * will update a connectivity value of 12 as follows:
     *      connectivity = 12 / (5 + 1)
     * NOTE: we add one to the number of intermediates to include the
     * original parent.
     */
    void split_node_data(HyperNode& data, int num_intermediates);

    /*
     * Used to read from the data file. For each line of the file
     * that includes data for each node, this method reads a single
     * criterion value into the current HyperNode. The string variable
     * criterion represents the current criterion that is being read in.
     */

    void read_dam_criterion(int i, std::stringstream& reader, Dam& dam);

    void initialize_nodes(int num_nodes);
    /*
     * This method tests to see if a line is a comment line. If the line is a comment line then it is
     * skipped until a non-comment line is reached. This function then consumes the flag characters
     * of that line and returns with the stringstream holding the contents of that line.
     */
    void read_comments(std::ifstream& input, std::stringstream& reader);

};


#endif //AMAZON_PROJECT_HYPERNET_H
