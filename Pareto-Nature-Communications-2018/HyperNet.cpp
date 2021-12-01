//
// Created by Jonathan Gomes Selman on 7/6/17.
//
// Refer to the file describing the input file for hyper_net.
// File: format-dams-v1_xiaojian.docx, found in cpp_inputs file

#include "HyperNet.h"

HyperNet::HyperNet(const std::string &basin, const std::string &file, double epsilon, bool maximize_c1, bool maximize_c2)
        : river_network(basin), filename(file), maximize_c1(maximize_c1), maximize_c2(maximize_c2), epsilon(epsilon){
    construct();
}

void HyperNet::construct() {
    std::ifstream input;
    input.open(this->filename);

    if (input.fail()) { // Check to make sure the file was properly opened
        std::cerr << "error opening file" << std::endl;
        std::cerr << "Check to make sure the file path is correct" << std::endl;
        exit(1);
    }

    std::stringstream reader; // Used for processing lines

    // First line
    read_comments(input, reader); // Check for comments and process the flag character
    reader >> this->num_nodes >> this->num_edges ;
    this->num_dams = this->num_edges;

    // Dam-criterion line
    read_comments(input, reader);
    reader >> this->c1_name;
    reader >> this->c2_name;

    std::cout << this->c1_name<<std::endl;
    std::cout << this->c2_name<<std::endl;

    // Dam Data
    readDamInfo(input, num_dams);
    initialize_nodes(num_nodes);
    // Edge List
    readAdjList(input, num_edges);

    // After reading in the data we want to update the energies for the dams
    // Using the rounding scheme proposed for AAAI paper. This is
    // the pre-rounding step.
    calculate_k();
    round_criteria(root);

    root_node = nullptr;

    generateBinaryTree();
}


HyperNet::~HyperNet() {
    destructor_helper(root_node);
    root_node = nullptr;
}

HyperNet::HyperNet(const HyperNet &src): num_nodes(src.num_nodes), num_edges(src.num_edges), num_dams(src.num_dams),
                                                root(src.root), adj_lists(src.adj_lists),
                                                nodes(src.nodes), dams(src.dams), maximize_c1(src.maximize_c1),
                                                maximize_c2(src.maximize_c2), relevant_criteria(src.relevant_criteria),
                                                filename(src.filename), river_network(src.river_network),
                                                num_nodes_binary(src.num_nodes_binary), epsilon(src.epsilon),
                                                min_c1(src.min_c1), min_c2(src.min_c2), c1_name(src.c1_name), c2_name(src.c2_name) {
    copy_constructor_helper(root_node, src.root_node);
}

void HyperNet::copy_constructor_helper(HyperTreeNode *&newTree, HyperTreeNode *oldTree) {
    // Copy the binary oldTree into newTree
    if (oldTree != nullptr) {
        newTree = new HyperTreeNode(*oldTree);
        // Recurse
        copy_constructor_helper(newTree->left, oldTree->left);
        copy_constructor_helper(newTree->right, oldTree->right);
    }
}

HyperNet &HyperNet::operator=(const HyperNet &src) {
    // Generate new tree
    HyperTreeNode* tmpRoot;
    copy_constructor_helper(tmpRoot, src.root_node);

    // Copy instance variables
    num_nodes = src.num_nodes;
    num_edges = src.num_edges;
    num_dams = src.num_dams;
    root = src.root;
    adj_lists = src.adj_lists;
    nodes = src.nodes;
    dams = src.dams;
    maximize_c1=src.maximize_c1;
    maximize_c2=src.maximize_c2;
    relevant_criteria = src.relevant_criteria;
    filename = src.filename;
    river_network = src.river_network;
    num_nodes_binary = src.num_nodes_binary;
    epsilon = src.epsilon;
    min_c1 = src.min_c1;
    min_c2 = src.min_c2;
    c1_name = src.c1_name;
    c2_name = src.c2_name;
    // Free old tree
    destructor_helper(root_node);

    this->root_node = tmpRoot;
    return *this;
}

void HyperNet::destructor_helper(HyperTreeNode* node) {
    // Free tree with post-order traversal
    if (node != nullptr) {
        destructor_helper(node->left);
        destructor_helper(node->right);
        delete node;
    }
}


void HyperNet::readAdjList(std::ifstream& input, int num_edges) {
    // Read in the root line
    std::stringstream reader;
    read_comments(input, reader); // Check for comments and process the root flag
    reader >> root; // Get root id

    // Read through each edge
    for (int i = 0; i < num_edges; i++) {
        // Get each line to be processed
        read_comments(input, reader); // Check for comments and process the root flag

        int parent;
        Edge edge_to_child;
        reader >> parent >> edge_to_child.endNode >> edge_to_child.damId;

        // Put the new edge into the adjacency list map
        adj_lists[parent].push_back(edge_to_child);


        // Split the criteria of the dam between parent and child
        nodes[edge_to_child.endNode].c1 += dams[edge_to_child.damId].c1 / 2;
        nodes[parent].c1 += dams[edge_to_child.damId].c1 / 2;

        nodes[edge_to_child.endNode].c2 += dams[edge_to_child.damId].c2 / 2;
        nodes[parent].c2 += dams[edge_to_child.damId].c2 / 2;

    }
}



void HyperNet::readDamInfo(std::ifstream& input, int num_dams) {
    // Initialize the min_c1 and min_c2 to -1 to signal that we have not read in any dam data
    min_c1 = -1;
    min_c2 = -1;

    dams.reserve(num_dams);
    for (int i = 0; i < num_dams; i++){
        Dam newDam;

        // Get each line to be processed
        std::stringstream reader;
        read_comments(input, reader); // Check for comments and process the root flag

        reader >> newDam.id;
        std::cout<< "dam: "<<newDam.id<<std::endl;
        for (int i=0;i<=2;i++) {
            read_dam_criterion(i, reader, newDam);
        }

        // Do we need this??? Are the comments correct??

        // Keeps track of the minimum dam criteria values
        if (newDam.c1!=0) {
            if (min_c1 == -1) {
                min_c1 = newDam.c1;
            } else {
                min_c1 = std::fmin(min_c1, newDam.c1);
            }
        }

        // Keeps track of minimum seismic risk
        if (newDam.c2!= 0) {
            if (min_c2 == -1) {
                min_c2 = newDam.c2;
            } else {
                min_c2 = std::fmin(min_c2, newDam.c2);
            }
        }


        // Record the dam status (i.e. planned or already built)
        if (newDam.status == 0) {
            newDam.decision = {0, 1};
        } else if (newDam.status == 1) {
            newDam.decision = {1};
        } else {
            newDam.decision = {0};
        }
        this->dams[newDam.id] = newDam;
    }
}

void HyperNet::initialize_nodes(int num_nodes) {
    // Initialize the size of the node vector
    nodes.resize(num_nodes);
    for (int i = 0; i < num_nodes; i++) {
        HyperNode newNode;

        newNode.id=i; // Read in the node id

        this->nodes[newNode.id] = newNode;
    }
}

void HyperNet::calculate_k() {
    k_c1 = epsilon * min_c1 / 2.0;
    k_c2 = epsilon * min_c2 / 2.0;
}

void HyperNet::round_criteria(int node) {
    for (Edge& e: adj_lists[node]){
        // Round the edge value
        dams[e.damId].c1 = round_kp(dams[e.damId].c1, k_c1);
        dams[e.damId].c2 = round_kp(dams[e.damId].c2, k_c2);
        round_criteria(e.endNode);
    }
}

double HyperNet::round_kp(double cr, double k_factor) {
    // Avoid the case of divide by zero
    if (k_factor == 0) {
        return cr;
    }
    return std::floor(cr / k_factor) * k_factor;
}


void HyperNet::read_dam_criterion(int i, std::stringstream &reader, Dam &dam) {
    if (i == 0) {
        reader >> dam.c1;
    } else if (i == 1) {
        reader >> dam.c2;
    } else if (i == 2) {
        reader >> dam.status;
    } else {
        std::cerr << "incorrect dam_criterion entered" << std::endl;
        exit(1);
    }
}

void HyperNet::read_comments(std::ifstream &input, std::stringstream &reader) {
    // Get the line to be checked
    std::string line;
    getline(input, line);

    // Clear and re-assign the contents of the stringstream
    reader.str(line);
    reader.clear();
    // Test the first characters to see if the line starts with '%' and is thus a comment.
    // Keep getting new lines until all comments are processed.
    // In the end the flag character(s) that appear at the beginning of the line
    // containing the data to be parsed is/are consumed.
    std::string dummy;
    reader >> dummy;
    while (dummy.compare(0, 1, "%") == 0) {
        getline(input, line);
        reader.str(line);
        reader.clear();
        reader >> dummy;
    }
}

void HyperNet::generateBinaryTree() {
    // Initialize before generating
    num_nodes_binary = 0;

   // Allow for either splitting data values or not between intermediary nodes
    HyperNode root_data = nodes[root];

    int num_intermediates_root = calculate_number_intermediate_nodes(adj_lists[root].size());
    split_node_data(root_data, num_intermediates_root);


    root_node = new HyperTreeNode(root_data, false);
    generate_binary_tree_helper(root_node, 0);


}

int HyperNet::calculate_number_intermediate_nodes(int num_children) {
    return num_children - 2;
}

void HyperNet::generate_binary_tree_helper(HyperTreeNode *node, int children_index) {
    // We are on a new node so increment count
    num_nodes_binary++;

    // If we have reached a leaf stop
    if (adj_lists[node->node_data.id].size() == 0) {
        return;
    }

    // Make the left child be the child at current children_index
    Edge currEdge = adj_lists[node->node_data.id][children_index];

    // Allow for splitting data or not
    HyperNode left_data = nodes[currEdge.endNode];

    int num_intermediates_left = calculate_number_intermediate_nodes(adj_lists[currEdge.endNode].size());
    split_node_data(left_data, num_intermediates_left);

    node->left = new HyperTreeNode(left_data, false);

    node->leftD = dams[currEdge.damId];

    //Recurse on left
    generate_binary_tree_helper(node->left, 0);

    // Check if we need to add an intermediate for the right child
    if (adj_lists[node->node_data.id].size() - (children_index) > 2) {
        // Add a new intermediate node --- NOTE: the intermediate has the same node_data as its parent
        // We are splitting so it is the same data as parent
        node->right = new HyperTreeNode(node->node_data, true);

        // Create dummy dam connecting the intermediate to its parent
        node->rightD.c1= 0;
        node->rightD.c2 = 0;
        node->rightD.status = 2;
        node->rightD.decision = {0}; // The dam connecting a parent to an intermediate is never built

        // Recurse
        generate_binary_tree_helper(node->right, children_index + 1);
    } else if (adj_lists[node->node_data.id].size() > 1){ // Make sure there is a second child
        // Add a new 'parent' node to the right (i.e. not an intermediate)
        Edge rightEdge = adj_lists[node->node_data.id][children_index + 1];

        HyperNode right_data = nodes[rightEdge.endNode];

        int num_intermediates_right = calculate_number_intermediate_nodes(adj_lists[rightEdge.endNode].size());
        split_node_data(right_data, num_intermediates_right);


        node->right = new HyperTreeNode(right_data, false);
        node->rightD = dams[rightEdge.damId];

        // Recurse on the new right child
        generate_binary_tree_helper(node->right, 0);
    }
}

void HyperNet::print_adjacent_tree(int node) {
    for (Edge& e: adj_lists[node]) {
        print_adjacent_tree(e.endNode);
    }
    std::cout << nodes[node].id << std::endl;
}

void HyperNet::print_binary(HyperTreeNode *node) {
    if (node == nullptr) {
        return;
    }

    print_binary(node->left);
    print_binary(node->right);

    std::cout << node->node_data.id << " is intermediate: " << node->intermediate << std::endl;
}


void HyperNet::split_node_data(HyperNode &data, int num_intermediates) {
    if (num_intermediates > 0) {

        data.c1 = data.c1 / (num_intermediates + 1);
        data.c2 = data.c2 / (num_intermediates + 1);
    }
}


