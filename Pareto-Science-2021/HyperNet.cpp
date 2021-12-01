//
// Created by Jonathan Gomes Selman on 7/6/17.
//
// Refer to the file describing the input file for hyper_net.
// File: format-dams-v1_xiaojian.docx, found in cpp_inputs file

#include "HyperNet.h"

HyperNet::HyperNet(const std::string &basin, const std::string &file, bool split_values, bool balance_tree,
                   double epsilon)
        : river_network(basin), filename(file),
          split_values(split_values), balance_tree(balance_tree), epsilon(epsilon){
    construct();
}

HyperNet::HyperNet(const std::vector<std::string> &relevant_criteria) : relevant_criteria(relevant_criteria) {
    construct();
}

void HyperNet::construct() {
    std::ifstream input;
    input.open(this->filename);

    if (input.fail()) { // Check to make sure the file was properly opened
        std::cerr << "error opening file" << std::endl;
        exit(1);
    }

    std::stringstream reader; // Used for processing lines

    // First line
    read_comments(input, reader); // Check for comments and process the flag character
    reader >> this->num_nodes >> this->num_edges >> this->num_criteria; // Assign values from first line
    reader >> this->num_criteria_dam;
    this->num_dams = this->num_edges;

    // Dam-criterion line
    read_comments(input, reader); // Check for comments and process the flag characters
    dam_criteria.resize(this->num_criteria_dam);
    for (int i = 0; i < this->num_criteria_dam; i++) {
        reader >> this->dam_criteria[i];
    }

    // Node-criterion line
    read_comments(input, reader); // Check for comments and process the flag characters
    node_criteria.resize(this->num_criteria);
    for (int i = 0; i < this->num_criteria; i++) {
        reader >> this->node_criteria[i];
    }

    // Dam Data
    readDamInfo(input, num_dams);
    // Node Data
    readNodeData(input, num_nodes);
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
                                            num_criteria(src.num_criteria), root(src.root), adj_lists(src.adj_lists),
                                                nodes(src.nodes), dams(src.dams),
                                                num_criteria_dam(src.num_criteria_dam), node_criteria(src.node_criteria),
                                                dam_criteria(src.dam_criteria), relevant_criteria(src.relevant_criteria),
                                                filename(src.filename), river_network(src.river_network),
                                                num_nodes_binary(src.num_nodes_binary),
                                                split_values(src.split_values), balance_tree(src.balance_tree),
                                                min_dam_energy(src.min_dam_energy), epsilon(src.epsilon),
                                                min_seismic(src.min_seismic), min_biodiversity(src.min_biodiversity), min_ghg(src.min_ghg)
                                                , min_dor(src.min_dor), min_population(src.min_population){
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
    num_criteria = src.num_criteria;
    root = src.root;
    adj_lists = src.adj_lists;
    nodes = src.nodes;
    dams = src.dams;
    num_criteria_dam = src.num_criteria_dam;
    node_criteria = src.node_criteria;
    dam_criteria = src.dam_criteria;
    relevant_criteria = src.relevant_criteria;
    filename = src.filename;
    river_network = src.river_network;
    num_nodes_binary = src.num_nodes_binary;
    split_values = src.split_values;
    balance_tree = src.balance_tree;
    min_dam_energy = src.min_dam_energy;
    epsilon = src.epsilon;
    min_seismic = src.min_seismic;
    min_biodiversity = src.min_biodiversity;
    min_ghg = src.min_ghg;
    min_dor = src.min_dor;
    min_population = src.min_population;

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


        // Split the energy of the dam between parent and child
        nodes[edge_to_child.endNode].energy_associated += dams[edge_to_child.damId].energy / 2;
        nodes[parent].energy_associated += dams[edge_to_child.damId].energy / 2;

        // Do the same for the seismic risk energy of the dam

        nodes[edge_to_child.endNode].seismic_associated += dams[edge_to_child.damId].seismic / 2;
        nodes[parent].seismic_associated += dams[edge_to_child.damId].seismic / 2;

        nodes[edge_to_child.endNode].biodiversity_associated += dams[edge_to_child.damId].biodiversity / 2;
        nodes[parent].biodiversity_associated += dams[edge_to_child.damId].biodiversity / 2;

        nodes[edge_to_child.endNode].ghg_associated += dams[edge_to_child.damId].ghg / 2;
        nodes[parent].ghg_associated += dams[edge_to_child.damId].ghg / 2;

        nodes[edge_to_child.endNode].dor_associated += dams[edge_to_child.damId].dor / 2;
        nodes[parent].dor_associated += dams[edge_to_child.damId].dor / 2;

        nodes[edge_to_child.endNode].population_associated += dams[edge_to_child.damId].population / 2;
        nodes[parent].population_associated += dams[edge_to_child.damId].population / 2;

    }
}


void HyperNet::readNodeData(std::ifstream& input, int num_nodes) {
    // Initialize the size of the node vector
    nodes.resize(num_nodes);
    for (int i = 0; i < num_nodes; i++) {
        HyperNode newNode;

        // Get each line to be processed
        std::stringstream reader;
        read_comments(input, reader); // Check for comments and process the root flag

        reader >> newNode.id; // Read in the node id
        for (const std::string& criterion: this->node_criteria) { // Read the criterion in order dictated by file
            read_criterion(criterion, reader, newNode);
        }
        // Zero these out because we will be adding to them later
        newNode.energy_associated = 0;
        newNode.seismic_associated = 0;
        newNode.biodiversity_associated = 0;
        newNode.ghg_associated = 0;
        newNode.dor_associated = 0;
        newNode.population_associated = 0;
        this->nodes[newNode.id] = newNode;
    }
}

void HyperNet::read_criterion(const std::string &criterion, std::stringstream &reader, HyperNode &newNode) {
    if (criterion == "connectivity") {
        reader >> newNode.connectivity;
    } else if (criterion == "sediment") {
        reader >> newNode.sediment;
    } else {
        std::cerr << "incorrect node criterion entered" << std::endl;
        exit(1);
    }
}


void HyperNet::readDamInfo(std::ifstream& input, int num_dams) {
    // Initialize the min_dam_energy to -1 to signal that we have not read in any dam data
    min_dam_energy = -1;
    min_seismic = -1;
    min_biodiversity = -1;
    min_ghg = -1;
    min_dor = -1;
    min_population = -1;

    dams.reserve(num_dams);
    for (int i = 0; i < num_dams; i++){
        Dam newDam;

        // Get each line to be processed
        std::stringstream reader;
        read_comments(input, reader); // Check for comments and process the root flag

        reader >> newDam.id;
        for (const std::string& criterion: this->dam_criteria) {
            read_dam_criterion(criterion, reader, newDam);
        }

        // Keeps track of the minimum dam energy
        if (min_dam_energy == -1) {
            min_dam_energy = newDam.energy;
        } else {
            min_dam_energy = std::fmin(min_dam_energy, newDam.energy);
        }

        // Keeps track of minimum seismic risk
        if (newDam.seismic!= 0) {
            if (min_seismic == -1) {
                min_seismic = newDam.seismic;
            } else {
                min_seismic = std::fmin(min_seismic, newDam.seismic);
            }
        }

        if (newDam.biodiversity!= 0) {
            if (min_biodiversity == -1) {
                min_biodiversity = newDam.biodiversity;
            } else {
                min_biodiversity = std::fmin(min_biodiversity, newDam.biodiversity);
            }
        }

        if (newDam.ghg!= 0) {
            if (min_ghg == -1) {
                min_ghg = newDam.ghg;
            } else {
                min_ghg = std::fmin(min_ghg, newDam.ghg);
            }
        }

        if (newDam.dor!= 0) {
            if (min_dor == -1) {
                min_dor = newDam.dor;
            } else {
                min_dor = std::fmin(min_dor, newDam.dor);
            }
        }

        if (newDam.population!= 0) {
            if (min_population == -1) {
                min_population = newDam.population;
            } else {
                min_population = std::fmin(min_population, newDam.population);
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

void HyperNet::calculate_k() {
    dam_rounding_k = epsilon * min_dam_energy / 2.0;
    seismic_k = epsilon * min_seismic / 2.0;
    biodiversity_k = epsilon * min_biodiversity / 2.0;
    dor_k = epsilon * min_dor / 2.0;
    ghg_k = epsilon * min_ghg / 2.0;
    population_k = epsilon * min_population / 2.0;
}

void HyperNet::round_criteria(int node) {
    for (Edge& e: adj_lists[node]){
        // Round the edge value
        dams[e.damId].energy = round_energy(dams[e.damId].energy, dam_rounding_k);
        dams[e.damId].seismic = round_energy(dams[e.damId].seismic, seismic_k);
        dams[e.damId].biodiversity = round_energy(dams[e.damId].biodiversity, biodiversity_k);
        dams[e.damId].ghg = round_energy(dams[e.damId].ghg, ghg_k);
        dams[e.damId].dor = round_energy(dams[e.damId].dor, dor_k);
        dams[e.damId].population = round_energy(dams[e.damId].population, population_k);
        round_criteria(e.endNode);
    }
}

double HyperNet::round_energy(double energy, double k_factor) {
    // Avoid the case of divide by zero
    if (k_factor == 0) {
        return energy;
    }
    return std::floor(energy / k_factor) * k_factor;
}



void HyperNet::read_dam_criterion(const std::string &criterion, std::stringstream &reader, Dam &dam) {
    if (criterion == "energy") {
        reader >> dam.energy;
    } else if (criterion == "sed_trap") {
        reader >> dam.sediment;
    } else if (criterion == "seismic") {
        reader >> dam.seismic;
    } else if (criterion == "biodiversity") {
        reader >> dam.biodiversity;
    } else if (criterion == "ghg") {
        reader >> dam.ghg;
    } else if (criterion == "dor") {
        reader >> dam.dor;
    } else if (criterion == "population") {
        reader >> dam.population;
    } else if (criterion == "status") {
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

    if (balance_tree) { // Build a balanced tree
        generate_binary_tree_balanced(root_node, adj_lists[root], 0, (int) adj_lists[root].size() - 1, false, nodes[root]);
    } else {
        // Allow for either splitting data values or not between intermediary nodes
        HyperNode root_data = nodes[root];
        if (split_values) {
            int num_intermediates_root = calculate_number_intermediate_nodes(adj_lists[root].size());
            split_node_data(root_data, num_intermediates_root);
        }

        root_node = new HyperTreeNode(root_data, false);
        generate_binary_tree_helper(root_node, 0);
    }

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
    if (split_values) {
        int num_intermediates_left = calculate_number_intermediate_nodes(adj_lists[currEdge.endNode].size());
        split_node_data(left_data, num_intermediates_left);
    }
    node->left = new HyperTreeNode(left_data, false);

    node->leftD = dams[currEdge.damId];

    //Recurse on left
    generate_binary_tree_helper(node->left, 0);

    // Check if we need to add an intermediate for the right child
    if (adj_lists[node->node_data.id].size() - (children_index) > 2) {
        // Add a new intermediate node --- NOTE: the intermediate has the same node_data as its parent

        // Don't split data values (i.e. give zero for all intermediate values)
        if (!split_values) {
            HyperNode intermediate;
            intermediate.seismic_associated = 0;
            intermediate.biodiversity_associated = 0;
            intermediate.ghg_associated = 0;
            intermediate.dor_associated = 0;
            intermediate.population_associated = 0;
            intermediate.energy_associated = 0;
            intermediate.sediment = 0;
            intermediate.connectivity = 0;
            intermediate.id = node->node_data.id;
            node->right = new HyperTreeNode(intermediate, true);
        } else {
            // We are splitting so it is the same data as parent
            node->right = new HyperTreeNode(node->node_data, true);
        }

        // Create dummy dam connecting the intermediate to its parent
        node->rightD.energy = 0;
        node->rightD.seismic = 0;
        node->rightD.biodiversity = 0;
        node->rightD.ghg = 0;
        node->rightD.dor = 0;
        node->rightD.population = 0;
        node->rightD.sediment = 0;
        node->rightD.status = 2;
        node->rightD.decision = {0}; // The dam connecting a parent to an intermediate is never built

        // Recurse
        generate_binary_tree_helper(node->right, children_index + 1);
    } else if (adj_lists[node->node_data.id].size() > 1){ // Make sure there is a second child
        // Add a new 'parent' node to the right (i.e. not an intermediate)
        Edge rightEdge = adj_lists[node->node_data.id][children_index + 1];

        HyperNode right_data = nodes[rightEdge.endNode];
        if (split_values) {
            int num_intermediates_right = calculate_number_intermediate_nodes(adj_lists[rightEdge.endNode].size());
            split_node_data(right_data, num_intermediates_right);
        }

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

    std::cout << node->node_data.id << " is intermediate: " << node->intermediate << " connectivity: " << node->node_data.connectivity << std::endl;
}


void HyperNet::split_node_data(HyperNode &data, int num_intermediates) {
    if (num_intermediates > 0) {
        data.connectivity = data.connectivity / (num_intermediates + 1);
        data.sediment = data.sediment / (num_intermediates + 1);
        data.energy_associated = data.energy_associated / (num_intermediates + 1);
        data.seismic_associated = data.seismic_associated / (num_intermediates + 1);
        data.biodiversity_associated = data.biodiversity_associated / (num_intermediates + 1);
        data.ghg_associated = data.ghg_associated / (num_intermediates + 1);
        data.dor_associated = data.dor_associated / (num_intermediates + 1);
        data.population_associated = data.population_associated / (num_intermediates + 1);
    }
}

void
HyperNet::generate_binary_tree_balanced(HyperTreeNode *&node, std::vector<Edge>& children,
                                        int low_ch_bound, int high_ch_bound, bool intermediate,
                                        HyperNode& node_data) {
    num_nodes_binary++;
    if (children.empty()) { // Leaf Node
        node = new HyperTreeNode(node_data, false);
        return;
    }

    // Assign the value for the node
    if (!intermediate) {
        node = new HyperTreeNode(node_data, false); // Make new node from the actual tree;
        if (split_values) {
            split_node_data(node->node_data, calculate_number_intermediate_nodes((int)children.size()));
        }
    } else { // Intermediate node
        if (split_values) {
            node = new HyperTreeNode(node_data, true);
        } else { // Give the intermediate values of 0
            HyperNode intermediate;
            intermediate.seismic_associated = 0;
            intermediate.biodiversity_associated = 0;
            intermediate.ghg_associated = 0;
            intermediate.dor_associated = 0;
            intermediate.population_associated = 0;
            intermediate.energy_associated = 0;
            intermediate.sediment = 0;
            intermediate.connectivity = 0;
            intermediate.id = node->node_data.id;
            node = new HyperTreeNode(intermediate, true);
        }
    }

    // Distribute the children
    if (high_ch_bound - low_ch_bound > 1) { // More than two children left so split and divide
        int median = (low_ch_bound + high_ch_bound) / 2;
        // Make left intermediate
        // Create dummy dam connecting the intermediate to its parent
        node->leftD.status = 2;
        node->leftD.decision = {0}; // The dam connecting a parent to an intermediate is never built
        generate_binary_tree_balanced(node->left, children, low_ch_bound, median, true, node->node_data);
        // Check if right is intermediate or not
        if (high_ch_bound - (median + 1) == 0) { // Not intermediate
            Edge right_edge = children[high_ch_bound];
            node->rightD = dams[right_edge.damId]; // Assign the right dam
            generate_binary_tree_balanced(node->right, adj_lists[right_edge.endNode],
                                            0, (int)adj_lists[right_edge.endNode].size() - 1, false, nodes[right_edge.endNode]);
        } else { // Need intermediate
            node->rightD.status = 2;
            node->rightD.decision = {0}; // The dam connecting a parent to an intermediate is never built
            generate_binary_tree_balanced(node->right, children, median + 1, high_ch_bound, true, node->node_data);
        }
    } else {
        // Add real left child
        Edge left = children[low_ch_bound];
        node->leftD = dams[left.damId];
        generate_binary_tree_balanced(node->left, adj_lists[left.endNode], 0, (int)adj_lists[left.endNode].size() - 1,
                                        false, nodes[left.endNode]);
        // Check for right
        if (high_ch_bound - low_ch_bound > 0) {
            Edge right_edge = children[high_ch_bound];
            node->rightD = dams[right_edge.damId]; // Assign the right dam
            generate_binary_tree_balanced(node->right, adj_lists[right_edge.endNode], 0,
                                          (int) adj_lists[right_edge.endNode].size() - 1, false, nodes[right_edge.endNode]);
        }
    }

}

