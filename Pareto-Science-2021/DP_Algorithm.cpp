//
// Created by Jonathan Gomes Selman on 7/6/17.
//

#include "DP_Algorithm.h"

DP_Algorithm::DP_Algorithm(const HyperNet &net, int root, const std::vector<std::string> &relevant_criteria, double epsilon,
                           unsigned int seed, bool use_n_squared, int batch_size, bool use_binary_tree, size_t num_threads)
        : net(net), root(root), relevant_criteria(relevant_criteria), k_epsilon(epsilon),
          random_seed(seed), use_n_squared(use_n_squared), batch_size(batch_size), use_binary_tree(use_binary_tree),
          K(num_threads),pool(K){

    // Initialize counters
    this->total_pruned_policies = 0;
    this->total_policies_considered = 0;
    this->max_policies_considered = 0;

    // Initialize timers
    this->time_copying = 0;
    this->time_generating = 0;
    this->time_sorting = 0;
    this->time_nlogn = 0;

    // Calculate rounding constants
    calc_theoretical_ks();
    calc_theoretical_ks_tree(this->net.root_node);

    for (int i = 0; i < K; ++i){
        solution_sets.push_back(std::unordered_set<Pareto_Solution *, SolutionHash, EqualSolutions>());
    }

}

void DP_Algorithm::calc_theoretical_ks() {
    for (auto &node : net.nodes) {
        node.k_conn = k_epsilon * node.connectivity;
        node.k_sediment = k_epsilon * node.sediment;

        // NOTE: For the current rounding scheme these values cannot be
        // pre-calculated as they depend on the number of dams built etc.
        // However, the old scheme depends on this. Likely this can be removed!
        node.k_seismic = k_epsilon * node.seismic_associated;
        node.k_energy = k_epsilon * node.energy_associated;
        node.k_biodiversity = k_epsilon * node.biodiversity_associated;
        node.k_ghg = k_epsilon * node.ghg_associated;
        node.k_dor = k_epsilon * node.dor_associated;
        node.k_population = k_epsilon * node.population_associated;
    }
}

void DP_Algorithm::calc_theoretical_ks_tree(HyperTreeNode* node) {
    if (node != nullptr) {
        node->node_data.k_conn = k_epsilon * node->node_data.connectivity;
        node->node_data.k_sediment = k_epsilon * node->node_data.sediment;

        // NOTE: For the current rounding scheme these values cannot be
        // pre-calculated as they depend on the number of dams built etc.
        // However, the old scheme depends on this. Likely this can be removed!
        node->node_data.k_energy = k_epsilon * node->node_data.energy_associated;
        node->node_data.k_seismic = k_epsilon * node->node_data.seismic_associated;
        node->node_data.k_biodiversity = k_epsilon * node->node_data.biodiversity_associated;
        node->node_data.k_ghg = k_epsilon * node->node_data.ghg_associated;
        node->node_data.k_dor = k_epsilon * node->node_data.dor_associated;
        node->node_data.k_population = k_epsilon * node->node_data.population_associated;

        // Recurse
        calc_theoretical_ks_tree(node->left);
        calc_theoretical_ks_tree(node->right);
    }
}


void DP_Algorithm::computeDP_n_squared(int node) {
    std::vector<Edge> childrenEdges = net.adj_lists[node]; // Get the edges connecting to the children

    // Matrix where each row represents the Pareto_Solutions for a child.
    std::vector<std::vector<Pareto_Solution* > > children_pairs;
    children_pairs.reserve(net.adj_lists[node].size()); // We know what the size will be
    for (Edge& e: net.adj_lists[node]) {
        children_pairs.push_back({});
        for (Pareto_Solution* curr = pareto_partial_table[e.endNode].get_head(); curr != nullptr; curr = curr->next) {
            children_pairs.back().push_back(curr);
        }
    }

    // For each dam assign the possible build decisions
    std::vector<std::vector<int> > decisions;
    decisions.reserve(childrenEdges.size());
    for (Edge e: childrenEdges) {
        decisions.push_back(net.dams[e.damId].decision);
    }

    std::cout << "Start - Node: " << node << std::endl;
    unsigned long num_considered = 1;

    // Check for leaf node
    if (childrenEdges.empty()) {
        // List of criteria values
        std::pair<double, double> criteria[MAX_CRITERIA];
        initialize_criteria(node, criteria);

        // Book keeping
        std::cout << "Considering: 1" << std::endl;
        total_policies_considered++;

        this->pareto_partial_table[node].add_n_squared(criteria, node, {}, {}, (int) relevant_criteria.size(), false);
    } else {
        // Get all possible decision vectors and children pairs
        std::vector<std::vector<Pareto_Solution* > > cartesian_pairs = cartesian_product_generic(children_pairs);
        std::vector<std::vector<int> > decisionVectors = cartesian_product_generic(decisions);

        // Book keeping
        num_considered = cartesian_pairs.size() * decisionVectors.size();
        std::cout << "Considering: " << num_considered << std::endl;
        total_policies_considered += num_considered;
        max_policies_considered = std::max(max_policies_considered, num_considered);

        // Loop through every possible partial policy
        for (auto &decisionVec: decisionVectors) {
            for (auto &grouping: cartesian_pairs) {
                std::pair<double, double> criteria[MAX_CRITERIA];
                initialize_criteria(node, criteria);

                // Add the corresponding criteria for each of the children
                for (int i = 0; i < decisionVec.size(); i++) {
                    update_criteria(node, grouping[i], criteria, decisionVec[i], i);
                }
                round_criteria(node, criteria); // Round the criteria values based on local ks

                this->pareto_partial_table[node].add_n_squared(criteria, node, decisionVec, grouping,
                                                   (int) relevant_criteria.size(), false);
            }
        }
    }

    // Book keeping
    std::cout << "End - kept: " << this->pareto_partial_table[node].getNum_solutions() << "\tdropped: " <<
              (num_considered - pareto_partial_table[node].getNum_solutions()) << std::endl;
    total_pruned_policies += num_considered - pareto_partial_table[node].getNum_solutions();
}

void DP_Algorithm::initialize_criteria(int node, std::pair<double, double> *criteria) {
    for (int i = 0; i < relevant_criteria.size(); i++) {
        if (relevant_criteria[i] == "connectivity") {
            criteria[i] = {net.nodes[node].connectivity, net.nodes[node].connectivity};
        } else if (relevant_criteria[i] == "sediment") {
            criteria[i] = {net.nodes[node].sediment, net.nodes[node].sediment};
        } else if (relevant_criteria[i] == "energy" || relevant_criteria[i] == "seismic" || relevant_criteria[i] == "biodiversity" ||
                relevant_criteria[i] == "ghg" || relevant_criteria[i] == "dor" || relevant_criteria[i] == "population"  ) {
            criteria[i] = {0, 0};
        } else {
            std::cerr << "Something went wrong when adding criteria values to leaf node" << std::endl;
            exit(1);
        }
    }
}

void DP_Algorithm::update_criteria(int node, Pareto_Solution *pareto_decision, std::pair<double, double> *criteria,
                                   int decision, int child_index) {
    if (decision == 0) { // We are not building the dam
        for (int i = 0; i < relevant_criteria.size(); i++) {
            // Note we may not need this if statement, as when no dam is built we always collect value from child
            if(relevant_criteria[i] == "connectivity" ||
               relevant_criteria[i] == "sediment" ||
               relevant_criteria[i] == "energy" ||
               relevant_criteria[i] == "seismic" ||
               relevant_criteria[i] == "biodiversity" ||
               relevant_criteria[i] == "ghg" ||
               relevant_criteria[i] == "dor" ||
               relevant_criteria[i] == "population" ) {
                // Add the true connectivity / sediment / energy of child
                criteria[i].first += pareto_decision->criteria_2[i].first;
                // Add the rounded connectivity / sediment / energy
                criteria[i].second += pareto_decision->criteria_2[i].second;
            }
        }
    } else { // Dam is being built
        // Calculations needed to determine sediment retention and energy gained by dam
        double sed_trap = net.dams[net.adj_lists[node][child_index].damId].sediment;
        double added_energy = net.dams[net.adj_lists[node][child_index].damId].energy;
        double seismic = net.dams[net.adj_lists[node][child_index].damId].seismic;
        double biodiversity = net.dams[net.adj_lists[node][child_index].damId].biodiversity;
        double ghg = net.dams[net.adj_lists[node][child_index].damId].ghg;
        double dor = net.dams[net.adj_lists[node][child_index].damId].dor;
        double population = net.dams[net.adj_lists[node][child_index].damId].population;

        for (int i = 0; i < relevant_criteria.size(); i++) {
            if(relevant_criteria[i] == "sediment") {
                // Add the sediment percentage that passes through dam
                criteria[i].first += pareto_decision->criteria_2[i].first * (1 - sed_trap);
                // Add the rounded sediment percentage that passes through dam
                criteria[i].second += pareto_decision->criteria_2[i].second * (1 - sed_trap);
            } else if (relevant_criteria[i] == "energy") {
                // Add energy gained from built dam and energy from child
                criteria[i].first += added_energy + pareto_decision->criteria_2[i].first;
                // Add energy gained from built dam and rounded energy from child
                criteria[i].second += added_energy + pareto_decision->criteria_2[i].second;
            } else if (relevant_criteria[i] == "seismic") { // We should check this
                // Add the risk associated with the dam in its true value and that accumulated by child
                criteria[i].first += seismic + pareto_decision->criteria_2[i].first;;
                // Add the negative of the risk associated so we can maximize negative and that from child
                criteria[i].second += -1 * seismic + pareto_decision->criteria_2[i].second;
            } else if (relevant_criteria[i] == "biodiversity") { // We should check this
                // Add the risk associated with the dam in its true value and that accumulated by child
                criteria[i].first += biodiversity + pareto_decision->criteria_2[i].first;;
                // Add the negative of the risk associated so we can maximize negative and that from child
                criteria[i].second += -1 * biodiversity + pareto_decision->criteria_2[i].second;
            } else if (relevant_criteria[i] == "ghg") { // We should check this
                // Add the risk associated with the dam in its true value and that accumulated by child
                criteria[i].first += ghg + pareto_decision->criteria_2[i].first;;
                // Add the negative of the risk associated so we can maximize negative and that from child
                criteria[i].second += -1 * ghg + pareto_decision->criteria_2[i].second;
            } else if (relevant_criteria[i] == "dor") { // We should check this
                // Add the risk associated with the dam in its true value and that accumulated by child
                criteria[i].first += dor + pareto_decision->criteria_2[i].first;;
                // Add the negative of the risk associated so we can maximize negative and that from child
                criteria[i].second += -1 * dor + pareto_decision->criteria_2[i].second;
            } else if (relevant_criteria[i] == "population") { // We should check this
                // Add the risk associated with the dam in its true value and that accumulated by child
                criteria[i].first += population + pareto_decision->criteria_2[i].first;;
                // Add the negative of the risk associated so we can maximize negative and that from child
                criteria[i].second += -1 * population + pareto_decision->criteria_2[i].second;
            }
        }
    }
}

void DP_Algorithm::round_criteria(int node, std::pair<double, double> *criteria) {
    for (int i = 0; i < relevant_criteria.size(); i++) {
        if (relevant_criteria[i] == "connectivity") {
            criteria[i].second = round_with_k(criteria[i].second, net.nodes[node].k_conn);
        } else if (relevant_criteria[i] == "sediment") {
            criteria[i].second = round_with_k(criteria[i].second, net.nodes[node].k_sediment);
        } else if (relevant_criteria[i] == "energy") {
            criteria[i].second = round_with_k(criteria[i].second, net.nodes[node].k_energy);
        } else if (relevant_criteria[i] == "seismic") {
            criteria[i].second = round_with_k(criteria[i].second, net.nodes[node].k_seismic);
        } else if (relevant_criteria[i] == "biodiversity") {
            criteria[i].second = round_with_k(criteria[i].second, net.nodes[node].k_biodiversity);
        } else if (relevant_criteria[i] == "ghg") {
            criteria[i].second = round_with_k(criteria[i].second, net.nodes[node].k_ghg);
        } else if (relevant_criteria[i] == "dor") {
            criteria[i].second = round_with_k(criteria[i].second, net.nodes[node].k_dor);
        } else if (relevant_criteria[i] == "population") {
            criteria[i].second = round_with_k(criteria[i].second, net.nodes[node].k_population);
        }
    }
}

void DP_Algorithm::computeDP_n_squared_generic(HyperTreeNode *node, bool round_values) {
    std::vector<std::vector<Pareto_Solution* > > children_pairs;

    // LOOK AT THIS -- look into the pareto-pair for each solution and the decision vector. Should probably make
    // Always of size 2 even if only 1 child exists!

    // Add Pareto_pairs from left child
    if (node->left != nullptr) {
        children_pairs.emplace_back(); //Create new empty list that we will add to
        for (Pareto_Solution *curr = node->left->local_frontier.get_head(); curr != nullptr; curr = curr->next) {
            children_pairs.back().push_back(curr);
        }
    } else {
        //children_pairs.push_back({nullptr});
    }

    // Add Pareto_pairs from right child
    if (node->right != nullptr) {
        children_pairs.emplace_back(); //Create new empty list that we will add to
        for (Pareto_Solution *curr = node->right->local_frontier.get_head(); curr != nullptr; curr = curr->next) {
            children_pairs.back().push_back(curr);
        }
    } else {
        //children_pairs.push_back({nullptr});
    }

    // For each dam assign the possible build decisions
    // Note if there is only one child then it will be the left child!!!!
    std::vector<std::vector<int> > decisions;
    if (node->left != nullptr) {
        decisions.push_back(node->leftD.decision);
    } else {
        //decisions.push_back({-1});
    }
    if (node->right != nullptr) {
        decisions.push_back(node->rightD.decision);
    } else {
        //decisions.push_back({-1});
    }

    // We can pre-compute this
    std::vector<std::vector<Pareto_Solution* > > cartesian_pairs = cartesian_product_generic(children_pairs);
    std::cout << "Start - Node: " << node->node_data.id << std::endl;
    unsigned long num_considered = 1;

    // Leaf node in the tree
    if (children_pairs.empty()) {
        // List of criteria values
        std::pair<double, double> criteria[MAX_CRITERIA];
        initialize_criteria(node, criteria);

        std::cout << "Considering: 1" << std::endl;
        total_policies_considered++;

        node->local_frontier.add_n_squared(criteria, node->node_data.id, {}, {}, (int) relevant_criteria.size(), false);
    } else {
        // Get all possible decision vectors
        std::vector<std::vector<int> > decisionVectors = cartesian_product_generic(decisions);
        // Used to track the number considered solutions
        num_considered = cartesian_pairs.size() * decisionVectors.size();
        std::cout << "Considering: " << num_considered << std::endl;
        total_policies_considered += num_considered;
        max_policies_considered = std::max(max_policies_considered, num_considered);

        // Loop through all possible dam decisions and
        // combinations of Pareto-Optimal policies from
        // the children.
        for (auto &decisionVec: decisionVectors) {
            bool calc_built = false;
            bool printed = false;
            for (auto &grouping: cartesian_pairs) {
                std::pair<double, double> criteria[MAX_CRITERIA];
                initialize_criteria(node, criteria);

                // Keep track of the energy built
                double energy_built = 0;
                // Keep track of the seismic risk taken
                double seismic_built = 0;
                double biodiversity_built = 0;
                double ghg_built = 0;
                double dor_built = 0;
                double population_built = 0;


                // Add the corresponding values for each of the Pareto pairs from the children
                for (int i = 0; i < decisionVec.size(); i++) {
                    // So that there is a decision for the left and right
                    //if (decisionVec[i] != -1) {
                    update_criteria_new(node, grouping[i], criteria, decisionVec[i], i, energy_built, seismic_built, biodiversity_built, ghg_built, dor_built, population_built);
                    //}
                }
                // Re-calculate the k used for energy and seismic for the given node
                if (!calc_built) {
                    assign_node_k(node, energy_built, seismic_built, biodiversity_built, ghg_built, dor_built, population_built);
                    calc_built = true;
                }

                // Round values
                if (round_values) {
                    round_criteria(node, criteria); // Round the criteria values based on local ks
                }
                // Use strict compare, n_squared algorithm
                node->local_frontier.add_n_squared(criteria, node->node_data.id, decisionVec, grouping,
                                                   (int) relevant_criteria.size(), false);
            }
        }
    }


    std::cout << "End - kept: " << node->local_frontier.getNum_solutions() << "\tdropped: " <<
              (num_considered - node->local_frontier.getNum_solutions()) << std::endl;
    total_pruned_policies += num_considered - node->local_frontier.getNum_solutions();
}

void DP_Algorithm::computeDP_nlogn_batching(HyperTreeNode *node) {

    std::cout << "Start - Node: " << node->node_data.id;
    unsigned long num_considered = 1;
    this->batch_number = 1;
    std::unordered_set<Pareto_Solution*, SolutionHash, EqualSolutions> all_solution_set;
    // Check if we are a leaf node
    if (node->left == nullptr && node->right == nullptr) {
        // List of criteria values
        std::pair<double, double> criteria[MAX_CRITERIA];

        initialize_criteria(node, criteria);

        std::cout << "\tConsidering: 1" << std::endl;
        total_policies_considered++;

        node->vec_frontier.push_back(new Pareto_Solution(node->node_data.id, criteria, {}, {}, (int) relevant_criteria.size()));
    } else { // Not leaf node


        clock_t t = clock();

        if (node->left != nullptr && node->right != nullptr) { // Has both children
            // Calculate the number that will be considered
            num_considered = node->left->vec_frontier.size() * node->right->vec_frontier.size() *
                                    node->leftD.decision.size() * node->rightD.decision.size();
            std::cout << "Considering: " << num_considered << std::endl;
            total_policies_considered += num_considered;
            max_policies_considered = std::max(max_policies_considered, num_considered);

            // Take all the combinations of the dam placements
            // Note we go *backwards* because when doing the cartesian
            // product it does it backwards
            if (node->right->vec_frontier.size()>2*K) {


                int n=node->right->vec_frontier.size()/K;
                std::cout<<"right vec frontier size: "<< node->right->vec_frontier.size()<<std::endl;
                std::cout<<"left vec frontier size: "<< node->left->vec_frontier.size()<<std::endl;
                for (int i = 0; i < K; ++i){
                    results.emplace_back(pool.enqueue([i,this,node,all_solution_set,n] {
                        // mtx.lock();
                        auto t = clock();
                        // printf("[%d] Start\n",i );
                        auto it1 = node->right->vec_frontier.begin();

                        std::advance(it1, i * n);

                        auto it2 = node->right->vec_frontier.begin();
                        if (i==K-1){
                            it2= node->right->vec_frontier.end();
                        }
                        else {

                            std::advance(it2, (i + 1) * n );
                        }
                        std::vector<Pareto_Solution* > partial_right(it1,it2);

                        int num_generated= 0;

                        for (int right_dam_plan: node->rightD.decision) {
                            for (int left_dam_plan: node->leftD.decision) {
                                // Create decision vector
                                std::vector<int> decision = {left_dam_plan, right_dam_plan};
                                // Take all combinations of the left and right tables
                                for (Pareto_Solution *currR: partial_right) {
                                    for (Pareto_Solution *currL: node->left->vec_frontier) {
                                        // Create pair

                                        std::vector<Pareto_Solution *> pair = {currL, currR};

                                        // Call method to generate new solution
                                        add_nlogn_solution(pair, decision, node, solution_sets[i],num_generated, t,
                                                           this->batch_number);

                                    }
                                }

                            }
                        }
                        return i;
                    }
                    ) );
                }
                for(auto && result: results) {
                    int i = result.get();

                }
                results.clear();

                for(int i = 0; i < K; ++i) {
                    all_solution_set.insert(solution_sets[i].begin(),solution_sets[i].end());
                    std::cout<<"solution size: "<< all_solution_set.size()<<std::endl;
                    solution_sets[i].clear();
                }
            }
            else {
                int num_generated = 0;
                for (int right_dam_plan: node->rightD.decision) {
                    for (int left_dam_plan: node->leftD.decision) {
                        // Create decision vector
                        std::vector<int> decision = {left_dam_plan, right_dam_plan};
                        // Take all combinations of the left and right tables
                        for (Pareto_Solution *currR: node->right->vec_frontier) {
                            for (Pareto_Solution *currL: node->left->vec_frontier) {
                                // Create pair
                                std::vector<Pareto_Solution *> pair = {currL, currR};

                                // Call method to generate new solution
                                add_nlogn_solution(pair, decision, node, all_solution_set, num_generated, t,
                                                   batch_number);
                            }
                        }

                    }
                }
            }
        } else if (node->left != nullptr) { // Just left
            // Track policies
            int num_generated = 0;
            num_considered = node->left->vec_frontier.size() * node->leftD.decision.size();
            std::cout << "Considering: " << num_considered << std::endl;
            total_policies_considered += num_considered;
            max_policies_considered = std::max(max_policies_considered, num_considered);
            // Loop through left dam
            for (int left_dam_plan: node->leftD.decision) {
                std::vector<int> decision = {left_dam_plan, -1};

                for (Pareto_Solution* currL: node->left->vec_frontier) {
                    std::vector<Pareto_Solution *> pair = {currL, nullptr};

                    // Call method to generate new solution
                    add_nlogn_solution(pair, decision, node, all_solution_set, num_generated, t,
                                       batch_number);
                }
            }
        } else { // Note this can never happen!!!!! Because we always put child as left first
            // Track policies
            int num_generated = 0;
            num_considered = node->right->vec_frontier.size() * node->rightD.decision.size();
            std::cout << "Considering: " << num_considered << std::endl;
            total_policies_considered += num_considered;
            max_policies_considered = std::max(max_policies_considered, num_considered);
            // Loop through left dam
            for (int right_dam_plan: node->rightD.decision) {
                std::vector<int> decision = {-1, right_dam_plan};

                for (Pareto_Solution *currR: node->right->vec_frontier) {
                    std::vector<Pareto_Solution *> pair = {nullptr, currR};

                    // Call method to generate new solution
                    add_nlogn_solution(pair, decision, node, all_solution_set, num_generated, t,
                                       batch_number);

                }
            }
        }
        t = clock() - t;
        time_generating += ((double)t)/CLOCKS_PER_SEC;

        // Copying
        t = clock();
        // Add the solutions from the hash set to a vector.
        // No ties exist
        std::vector<Pareto_Solution* > solution_set(all_solution_set.begin(), all_solution_set.end());
        t = clock() - t;
        time_copying += ((double)t)/CLOCKS_PER_SEC;

        std::cout << "Num-actually considered: " << solution_set.size() << std::endl << std::endl;

        // Time nlogn
        t = clock();
        // Perform the divide and conquer
        node->vec_frontier = divide_and_conquer(solution_set, (int) relevant_criteria.size());
        t = clock() - t;
        time_nlogn += ((double)t)/CLOCKS_PER_SEC;
    }
    std::cout << "End - kept: " << node->vec_frontier.size() << "\tdropped: " <<
              (num_considered - node->vec_frontier.size()) << std::endl;
    total_pruned_policies += num_considered - node->vec_frontier.size();
}

void
DP_Algorithm::add_nlogn_solution(std::vector<Pareto_Solution *> &grouping, std::vector<int> &decisionVec,
                                 HyperTreeNode *node,
                                 std::unordered_set<Pareto_Solution *, SolutionHash, EqualSolutions> &solution_set,
                                 int &num_generated, clock_t &t, int &batch_num) {
    std::pair<double, double> criteria[MAX_CRITERIA];
    initialize_criteria(node, criteria);

    // Keep track of the energy built
    double energy_built = 0;
    // Keep track of the seismic risk taken
    double seismic_built = 0;
    double biodiversity_built = 0;
    double ghg_built = 0;
    double dor_built = 0;
    double population_built = 0;

    // Add the corresponding values for each of the Pareto pairs from the children
    for (int i = 0; i < decisionVec.size(); i++) {
        if (decisionVec[i] != -1) {
            update_criteria_new(node, grouping[i], criteria, decisionVec[i], i, energy_built, seismic_built, biodiversity_built, ghg_built, dor_built, population_built);
        }
    }
    // Later we should not update this everytime!
    assign_node_k(node, energy_built, seismic_built, biodiversity_built, ghg_built, dor_built, population_built);

    round_criteria(node, criteria); // Round the criteria values based on local ks

    auto adding = new Pareto_Solution(node->node_data.id, criteria, decisionVec, grouping,
                                      (int) relevant_criteria.size());
    auto find = solution_set.find(adding);
    if (find != solution_set.end()) { // Solution tie exists!
        int rand_num = rand() % 2; // Only want to do if epsilon != 0 i.e. we are rounding!!!!

        if (rand_num == 0) { // Remove the current and replace it with the new node
            Pareto_Solution* old = *find;

            solution_set.erase(find);
            delete (old);
            solution_set.insert(adding);
        } else {
            delete(adding);
        }
    } else {
        solution_set.insert(adding);
    }

    num_generated++;
    // Lets do a trick to dynamically prune
    // Process a batch
    if (batch_size != 0 && num_generated >= batch_size) {
        std::cout << "Batch: " << batch_num << std::endl;

        t = clock() - t;
        mtx.lock();
        time_generating += ((double)t)/CLOCKS_PER_SEC;
        mtx.unlock();
        t = clock();
        std::vector<Pareto_Solution* > batched_solutions(solution_set.begin(), solution_set.end());
        t = clock() - t;
        mtx.lock();
        time_copying += ((double)t)/CLOCKS_PER_SEC;
        mtx.unlock();
        t = clock();
        std::vector<Pareto_Solution* > temp_answers = divide_and_conquer(batched_solutions, (int) relevant_criteria.size());
        t = clock() - t;
        mtx.lock();
        time_nlogn += ((double)t)/CLOCKS_PER_SEC;
        mtx.unlock();

        t = clock();
        solution_set = std::unordered_set<Pareto_Solution*, SolutionHash, EqualSolutions> (temp_answers.begin(), temp_answers.end());
        t = clock() - t;
        mtx.lock();
        time_copying += ((double)t)/CLOCKS_PER_SEC;
        mtx.unlock();
        t = clock();
        num_generated = 0;
        mtx.lock();
        batch_num++;
        mtx.unlock();
    }

}


void DP_Algorithm::initialize_criteria(HyperTreeNode *node, std::pair<double, double> *criteria) {
    for (int i = 0; i < relevant_criteria.size(); i++) {
        if (relevant_criteria[i] == "connectivity") {
            criteria[i] = {node->node_data.connectivity, node->node_data.connectivity};
        } else if (relevant_criteria[i] == "sediment") {
            criteria[i] = {node->node_data.sediment, node->node_data.sediment};
        } else if (relevant_criteria[i] == "energy" || relevant_criteria[i] == "seismic" || relevant_criteria[i] == "biodiversity" ||
        relevant_criteria[i] == "ghg" || relevant_criteria[i] == "dor" || relevant_criteria[i] == "population") {
            criteria[i] = {0, 0};
        } else {
            std::cerr << "Something went wrong when adding criteria values to leaf node" << std::endl;
            exit(1);
        }
    }
}



// This is the updated version that keeps track of the energy we build so that we can round correctly
void DP_Algorithm::update_criteria_new(HyperTreeNode *node, Pareto_Solution *pareto_decision, std::pair<double, double> *criteria,
                                       int decision, int child_index, double &energy_built, double &seismic_built, double &biodiversity_built,
                                       double &ghg_built, double &dor_built, double &population_built) {
    if (decision == 0) { // We are not building the dam
        for (int i = 0; i < relevant_criteria.size(); i++) {
            // Note we may not need this if statement, as when no dam is built we always collect value from child
            if(relevant_criteria[i] == "connectivity" ||
               relevant_criteria[i] == "sediment" ||
               relevant_criteria[i] == "energy" ||
               relevant_criteria[i] == "seismic" ||
               relevant_criteria[i] == "biodiversity" ||
               relevant_criteria[i] == "ghg" ||
               relevant_criteria[i] == "dor" ||
               relevant_criteria[i] == "population" ) {
                // Add the true connectivity / sediment / energy of child
                criteria[i].first += pareto_decision->criteria_2[i].first;
                // Make sure this is positive
                criteria[i].second += std::abs(pareto_decision->criteria_2[i].second);
            }
        }
    } else { // Dam is being built
        // Calculations needed to determine sediment retention and energy gained by dam
        double sed_trap;
        double added_energy;
        double added_seismic;
        double added_biodiversity;
        double added_ghg;
        double added_dor;
        double added_population;
        if (child_index == 0) { // Left child
            sed_trap = node->leftD.sediment;
            added_energy = node->leftD.energy;
            added_seismic = node->leftD.seismic;
            added_biodiversity = node->leftD.biodiversity;
            added_ghg = node->leftD.ghg;
            added_dor = node->leftD.dor;
            added_population = node->leftD.population;
        } else { // Right child
            sed_trap = node->rightD.sediment;
            added_energy = node->rightD.energy;
            added_seismic = node->rightD.seismic;
            added_biodiversity = node->rightD.biodiversity;
            added_ghg = node->rightD.ghg;
            added_dor = node->rightD.dor;
            added_population = node->rightD.population;
        }
        // Keep track of the energy built
        energy_built += added_energy;
        seismic_built += added_seismic;
        biodiversity_built += added_biodiversity;
        ghg_built += added_ghg;
        dor_built += added_dor;
        population_built += added_population;


        for (int i = 0; i < relevant_criteria.size(); i++) {
            if (relevant_criteria[i] == "sediment") {
                // Add the sediment percentage that passes through dam
                criteria[i].first += pareto_decision->criteria_2[i].first * (1 - sed_trap);
                // Add the rounded sediment percentage that passes through dam
                criteria[i].second += pareto_decision->criteria_2[i].second * (1 - sed_trap);
            } else if (relevant_criteria[i] == "energy") {
                // Add energy gained from built dam and energy from child
                criteria[i].first += added_energy + pareto_decision->criteria_2[i].first;
                // Add energy gained from built dam and rounded energy from child
                criteria[i].second += added_energy + pareto_decision->criteria_2[i].second;
            } else if (relevant_criteria[i] == "seismic") { // We should check this
                // Add the risk associated with the dam in its true value and that accumulated by child
                criteria[i].first += added_seismic + pareto_decision->criteria_2[i].first;

                // We first work with positive values then make negative after rounding
                criteria[i].second += added_seismic + std::abs(pareto_decision->criteria_2[i].second);
            } else if (relevant_criteria[i] == "biodiversity") { // We should check this
                // Add the risk associated with the dam in its true value and that accumulated by child
                criteria[i].first += added_biodiversity + pareto_decision->criteria_2[i].first;

                // We first work with positive values then make negative after rounding
                criteria[i].second += added_biodiversity + std::abs(pareto_decision->criteria_2[i].second);
            } else if (relevant_criteria[i] == "ghg") { // We should check this
                // Add the risk associated with the dam in its true value and that accumulated by child
                criteria[i].first += added_ghg + pareto_decision->criteria_2[i].first;

                // We first work with positive values then make negative after rounding
                criteria[i].second += added_ghg + std::abs(pareto_decision->criteria_2[i].second);
            } else if (relevant_criteria[i] == "dor") { // We should check this
                // Add the risk associated with the dam in its true value and that accumulated by child
                criteria[i].first += added_dor + pareto_decision->criteria_2[i].first;

                // We first work with positive values then make negative after rounding
                criteria[i].second += added_dor + std::abs(pareto_decision->criteria_2[i].second);
            } else if (relevant_criteria[i] == "population") { // We should check this
                // Add the risk associated with the dam in its true value and that accumulated by child
                criteria[i].first += added_population + pareto_decision->criteria_2[i].first;

                // We first work with positive values then make negative after rounding
                criteria[i].second += added_population + std::abs(pareto_decision->criteria_2[i].second);
            }
        }
    }
}

void DP_Algorithm::assign_node_k(HyperTreeNode *node, double energy_built, double seismic_built, double biodiversity_built, double ghg_built, double dor_built, double population_built) {
    node->node_data.k_energy = std::floor(energy_built / net.min_dam_energy) * net.min_dam_energy * k_epsilon / 2.0;
    node->node_data.k_seismic = std::floor(seismic_built / net.min_seismic) * net.min_seismic * k_epsilon / 2.0;
    node->node_data.k_biodiversity = std::floor(biodiversity_built / net.min_biodiversity) * net.min_biodiversity * k_epsilon / 2.0;
    node->node_data.k_ghg = std::floor(ghg_built / net.min_ghg) * net.min_ghg * k_epsilon / 2.0;
    node->node_data.k_dor = std::floor(dor_built / net.min_dor) * net.min_dor * k_epsilon / 2.0;
    node->node_data.k_population = std::floor(population_built / net.min_population) * net.min_population * k_epsilon / 2.0;
}

void DP_Algorithm::round_criteria(HyperTreeNode *node, std::pair<double, double> *criteria) {
    for (int i = 0; i < relevant_criteria.size(); i++) {
        if (relevant_criteria[i] == "connectivity") {
            criteria[i].second = round_with_k(criteria[i].second, node->node_data.k_conn);
        } else if (relevant_criteria[i] == "sediment") {
            criteria[i].second = round_with_k(criteria[i].second, node->node_data.k_sediment);
        } else if (relevant_criteria[i] == "energy") {
            criteria[i].second = round_with_k(criteria[i].second, node->node_data.k_energy);
        } else if (relevant_criteria[i] == "seismic") {
            // In this new scheme we first round the positive values then go back to negative
            criteria[i].second = -1 * round_with_k(criteria[i].second, node->node_data.k_seismic);
        } else if (relevant_criteria[i] == "biodiversity") {
            // In this new scheme we first round the positive values then go back to negative
            criteria[i].second = -1 * round_with_k(criteria[i].second, node->node_data.k_biodiversity);
        } else if (relevant_criteria[i] == "ghg") {
            // In this new scheme we first round the positive values then go back to negative
            criteria[i].second = -1 * round_with_k(criteria[i].second, node->node_data.k_ghg);
        } else if (relevant_criteria[i] == "dor") {
            // In this new scheme we first round the positive values then go back to negative
            criteria[i].second = -1 * round_with_k(criteria[i].second, node->node_data.k_dor);
        } else if (relevant_criteria[i] == "population") {
            // In this new scheme we first round the positive values then go back to negative
            criteria[i].second = -1 * round_with_k(criteria[i].second, node->node_data.k_population);
        }
    }
}


template <typename T>
std::vector<std::vector<T>> DP_Algorithm::cartesian_product_generic(const std::vector<std::vector<T>> &table) {
    unsigned long num_sets = table.size();
    std::vector<unsigned long> set_lengths; // Get the lengths for each row in the table
    set_lengths.reserve(num_sets);
    std::vector<unsigned long> num_lists_remaining;
    num_lists_remaining.reserve(num_sets + 1);
    num_lists_remaining.push_back(1);
    for(unsigned long i = 0; i < num_sets; i++) {
        unsigned long list_size = table[i].size();
        set_lengths.push_back(list_size);
        num_lists_remaining.push_back(list_size * num_lists_remaining[i]);
    }
    unsigned long num_products = num_lists_remaining[num_lists_remaining.size() - 1]; // Number of products that will be produced
    std::cout << "Cartesian size: " << num_products << std::endl;
    num_lists_remaining.pop_back();

    std::vector<std::vector<T> > result;
    result.reserve(num_products);

    for (unsigned long p = 0; p < num_products; p++) {
        std::vector<T> prod_set; // Build up each result
        prod_set.reserve(num_sets);

        for (unsigned long i = 0; i < num_sets; i++) {
            auto j = p / num_lists_remaining[i] % set_lengths[i];
            prod_set.push_back(table[i][j]);
        }
        result.push_back(prod_set);
    }
    return result;
}


double DP_Algorithm::round_with_k(double value, double k) {
    if (k == 0) {
        return value;
    }

    return std::floor(value / k) * k;
}


void DP_Algorithm::build_DP_table_recursive(int root) {
    for (Edge& e: net.adj_lists[root]) {
        build_DP_table_recursive(e.endNode);
    }

    computeDP_n_squared(root);
}


void DP_Algorithm::build_DP(bool round_values) {
    node_processed_counter = 0;
    build_DP_table_recursive_tree(net.root_node, round_values);
}


void DP_Algorithm::build_DP_table_recursive_tree(HyperTreeNode *root, bool round_values) {
    if (root != nullptr) {

        build_DP_table_recursive_tree(root->left, round_values);
        build_DP_table_recursive_tree(root->right, round_values);

        // Helpful tracker
        node_processed_counter++;
        std::cout << "Nodes processed so far: " << node_processed_counter << std::endl;

        if (use_n_squared) {
            computeDP_n_squared_generic(root, round_values);
        } else {
            computeDP_nlogn_batching(root);
        }

    }
}


void DP_Algorithm::print_DP_File(std::ofstream &file) {
    int num_sol = 0;
    file << "Node ID: " << net.root_node->node_data.id << "\n";
    for (Pareto_Opt_Node* n_ptr: net.root_node->vec_table) {
        num_sol++;
        //file << "[conn: " << n_ptr->connectivity << ", energy " << n_ptr->energy <<  "]\n";
        file << "[conn: " << n_ptr->connectivity << ", energy: " << n_ptr->energy << ", sediment: " << n_ptr->sediment <<  "]\n";
    }
    file <<  "\nNumber Pareto Solutions: " << num_sol << "\n";
}


void DP_Algorithm::print_DP_Output(std::ostream &stream) {
    // Output num_solutions
    stream << "num_solutions: ";
    if (use_binary_tree) {
        if (use_n_squared) {
            stream << net.root_node->local_frontier.getNum_solutions() << "\n";
        } else {
            stream << net.root_node->vec_frontier.size() << "\n";
        }
    } else {
        if (use_n_squared) {
            stream << this->pareto_partial_table[net.root].getNum_solutions() << "\n";
        }
    }
    // Give data about pruning
    if (use_binary_tree) {
        stream << "# pruning steps (# nodes): " << net.num_nodes_binary << "\n";
    } else {
        stream << "# pruning steps (# nodes): " << net.num_nodes << "\n";
    }
    stream << "Max policies considered: " << max_policies_considered << "\n";
    stream << "Policies considered: " << total_policies_considered << "\n";
    stream << "Pruned policies: " << total_pruned_policies << "\n";

    // Output the epsilon
    stream << "epsilon: " << this->k_epsilon << "\n";
    // Batch
    stream << "batch size: " << this->batch_size << "\n";
    // Output the relevant criteria
    stream << "criteria: ";
    for (int i = 0; i < relevant_criteria.size(); i++) {
        if (i != 0) {
            stream << ", ";
        }
        stream << relevant_criteria[i];
    }
    stream << "\n";

    // Output the info for each solution
    if (use_binary_tree) {
        if (use_n_squared) {
            print_dams_generic_tree_list(net.root_node, stream);
        } else {
            print_dams_generic_tree_vec(net.root_node, stream);
        }
    } else {
        if (use_n_squared) {
            print_dams_non_binary_list(net.root, stream);
        }
    }
}

void DP_Algorithm::print_dams_non_binary_list(int node, std::ostream &out) {
    for (Pareto_Solution *result = pareto_partial_table[node].get_head(); result != nullptr; result = result->next) {
        print_dams_helper(result, out, node);
    }
}

void DP_Algorithm::print_dams_helper(Pareto_Solution *result, std::ostream &out, int node) {
    std::stack<std::pair<int , Pareto_Solution* > > paths;
    paths.emplace(node, result);

    // Explore the "policy" to see what dams are built for the given solution
    std::vector<int> dams_built;
    while (!paths.empty()) {
        int decision_node = paths.top().first;
        Pareto_Solution* solution = paths.top().second;
        paths.pop();

        for (int i = 0; i < solution->dam_decisions.size(); i++) {
            int decision = solution->dam_decisions[i];
            if (decision == 1) {
                dams_built.push_back(net.adj_lists[decision_node][i].damId);
            }
            paths.emplace(net.adj_lists[decision_node][i].endNode, solution->pareto_decisions[i]);
        }
    }

    // Print the solutions as described by the ouput file format -- ***.word

    // Print the criteria values
    for (int i = 0; i < relevant_criteria.size(); i++) {
        out << result->criteria_2[i].first << ", ";
    }
    // Print number of dam
    out << dams_built.size() << ",";
    // Sort dam ids
    std::sort(dams_built.begin(), dams_built.end());
    // Print the dam ids
    for (int dam_id : dams_built) {
        out << " " << dam_id;
    }

    out << "\n";
}


void DP_Algorithm::print_dams_generic_tree_vec(HyperTreeNode *node, std::ostream &out) {
    for (Pareto_Solution *result: node->vec_frontier) {
        print_dams_helper(result, out, node);
    }
}


void DP_Algorithm::print_dams_generic_tree_list(HyperTreeNode *node, std::ostream &out) {
    for (Pareto_Solution *result = node->local_frontier.get_head(); result != nullptr; result = result->next) {
        print_dams_helper(result, out, node);
    }
}

void DP_Algorithm::print_dams_helper(Pareto_Solution *result, std::ostream &out, HyperTreeNode *node) {
    std::stack<std::pair<HyperTreeNode* , Pareto_Solution* > > paths;
    paths.emplace(node, result);

    // Explore the "policy" to see what dams are built for the given solution
    std::vector<int> dams_built;
    while (!paths.empty()) {
        HyperTreeNode* decision_node = paths.top().first;
        Pareto_Solution* decision = paths.top().second;
        paths.pop();

        if (decision_node->left != nullptr) { // Check if a left child exists
            if (decision->dam_decisions[0] == 1) { // Left dam is built
                dams_built.push_back(decision_node->leftD.id);
            }
            // Continue search
            paths.emplace(decision_node->left, decision->pareto_decisions[0]); // Left child with the Pareto-Solution used
        }
        if (decision_node->right != nullptr) {
            if (decision->dam_decisions[1] == 1) { // Right dam is built
                dams_built.push_back(decision_node->rightD.id);
            }
            // Continue the search
            paths.emplace(decision_node->right, decision->pareto_decisions[1]); // Right child with the Pareto-Solution used
        }
    }

    // Print the solutions as described by the ouput file format -- ***.word

    // Print the criteria values
    for (int i = 0; i < relevant_criteria.size(); i++) {
        out << result->criteria_2[i].first << ", ";
    }
    // Print number of dam
    out << dams_built.size() << ",";
    // Sort dam ids
    std::sort(dams_built.begin(), dams_built.end());
    // Print the dam ids
    for (int dam_id : dams_built) {
        out << " " << dam_id;
    }

    out << "\n";
}


std::vector<Pareto_Solution *> DP_Algorithm::L2D(std::vector<Pareto_Solution *> &all_possible) {
    // Sort by the 1st criteria (i.e. the criteria at index 1)
    clock_t t = clock();
    std::sort(all_possible.begin(), all_possible.end(), Compare_dimensions(2, (int)this->relevant_criteria.size()));
    t = clock() - t;
    time_sorting += ((double)t)/CLOCKS_PER_SEC;

    std::vector<Pareto_Solution* > solutions;
    double second_criteria_max;

    // Put the first element into the solution set and update _max_sec;
    solutions.push_back(all_possible[0]);
    second_criteria_max = all_possible[0]->criteria_2[0].second; // Get the last criteria
    int deleted = 0;
    for (unsigned long curr_index = 1; curr_index < all_possible.size(); curr_index++) {
        // If the last criteria is larger than the max so far update and include solution
        if (all_possible[curr_index]->criteria_2[0].second > second_criteria_max) {
            solutions.push_back(all_possible[curr_index]);
            second_criteria_max = all_possible[curr_index]->criteria_2[0].second;
        } else { // Free the memory for the dominated solution
            delete(all_possible[curr_index]);
            deleted++;
        }
    }
    std::cout << "Num deleted: " << deleted << std::endl;
    return solutions;
}

std::vector<Pareto_Solution *>
DP_Algorithm::divide_and_conquer(std::vector<Pareto_Solution *> &curr_set, int dimension) {
    if (dimension == 2) {
        return L2D(curr_set);
    }

    // Sort the vector by the first criteria and then call helper function
    std::sort(curr_set.begin(), curr_set.end(), Compare_dimensions(dimension, (int)this->relevant_criteria.size()));
    return divide_and_conquer_helper(curr_set, 0, curr_set.size() - 1, dimension);
}

std::vector<Pareto_Solution *>
DP_Algorithm::divide_and_conquer_helper(std::vector<Pareto_Solution *> &curr_set, unsigned long low, unsigned long high,
                                        int dimension) {
    if (high - low == 0) {
        // We want to return the single value
        return std::vector<Pareto_Solution* > (curr_set.begin() + low, curr_set.begin() + high + 1);
    }
    unsigned long midpoint = (high + low) / 2;
    std::vector<Pareto_Solution* > superior = divide_and_conquer_helper(curr_set, low, midpoint, dimension);
    std::vector<Pareto_Solution* > inferior = divide_and_conquer_helper(curr_set, midpoint + 1, high, dimension);

    // Take the union of superior and inferior, marking the inferior points, then sort by dimension - 1
    std::vector<Pareto_Solution* > sup_inf;
    sup_inf.reserve(superior.size() + inferior.size());
    // Copy those from superior
    for (Pareto_Solution* sol: superior) {
        sol->inferior = false;
        sup_inf.push_back(sol);
    }
    // Copy from inferior
    for (Pareto_Solution* sol: inferior) {
        sol->inferior = true;
        sup_inf.push_back(sol);
    }
    // Sort by d - 1 before calling marry
    std::sort(sup_inf.begin(), sup_inf.end(), Compare_dimensions(dimension - 1, (int)this->relevant_criteria.size()));

    std::vector<Pareto_Solution* > keep_inferior = marry(sup_inf, 0, sup_inf.size() - 1, dimension - 1);
    // Take the union of superior and the inferior points kept
    superior.insert(superior.end(), keep_inferior.begin(), keep_inferior.end());
    return superior;

}

std::vector<Pareto_Solution *>
DP_Algorithm::marry(std::vector<Pareto_Solution *> &curr_set, unsigned long low, unsigned long high, int dimension) {
    if (dimension == 2) { // We can use marry_2D to quickly find the optimal pairs
        return marry_2d(curr_set);
    }
    if (curr_set.size()==0){
        return curr_set;
    }
    // Create divide plane and count the inferior and superior in each section
    unsigned long median = (low + high) / 2;
    // We need to keep the superior points from the first half for future use
    std::vector<Pareto_Solution* > sup_first;
    // Lets count the number in each as a base case!
    unsigned long sup1 = 0, inf1 = 0, sup2 = 0, inf2 = 0;
    // First half
    for (unsigned long i = low; i <= median; i++) {

        if (curr_set[i]->inferior) {
            inf1++;
        } else {
            sup1++;
            sup_first.push_back(curr_set[i]);
        }

    }
    // Second half
    for (unsigned long i = median + 1; i <= high; i++) {
        if (curr_set[i]->inferior) {
            inf2++;
        }else {
            sup2++;
        }
    }
    // Check base cases!!
    if (inf1 + inf2 == high - low + 1) { // All inferior so we return all!
        return std::vector<Pareto_Solution* > (curr_set.begin() + low, curr_set.begin() + high + 1);
    } else if (sup1 + sup2 == high - low + 1) { // All superior so return none
        return {};
    }
    else if (inf1 == 0 && sup2 == 0) { // All superior points remain superior and inf remain inferior -- drop dimension
        std::vector<Pareto_Solution* > immediate_drop(curr_set.begin() + low, curr_set.begin() + high + 1);
        // Sort and drop dimension on curr set
        std::sort(immediate_drop.begin(), immediate_drop.end(), Compare_dimensions(dimension - 1, (int)this->relevant_criteria.size()));

        return marry(immediate_drop, 0, immediate_drop.size() - 1, dimension - 1);
    } else if (sup1 == 0 && inf2 == 0){ // Inferior points are in first half and superior in second
        return std::vector<Pareto_Solution* > (curr_set.begin() + low, curr_set.begin() + median + 1);
    }

    std::vector<Pareto_Solution* > inf_first = marry(curr_set, low, median, dimension);
    std::vector<Pareto_Solution* > inf_second = marry(curr_set, median + 1, high, dimension);

    // Now we must merge the sup_first with inf_second by looking at d - 1
    // Take the union of the two and sort
    sup_first.insert(sup_first.end(), inf_second.begin(), inf_second.end());
    // Sort by d - 1 before calling marry
    std::sort(sup_first.begin(), sup_first.end(), Compare_dimensions(dimension - 1, (int)this->relevant_criteria.size()));

    std::vector<Pareto_Solution* > keep_inf_second = marry(sup_first, 0, sup_first.size() - 1, dimension - 1);
    // Take union of inf_first and keep_int_second
    inf_first.insert(inf_first.end(), keep_inf_second.begin(), keep_inf_second.end());
    return inf_first;
}

std::vector<Pareto_Solution *> DP_Algorithm::marry_2d(std::vector<Pareto_Solution *> &curr_set) {
    // Solutions to be saved
    std::vector<Pareto_Solution* > saved_solutions;
    // Used to track the max first criteria seen so far
    double max_first = -1;
    for (Pareto_Solution* sol: curr_set) {
        // compare max_first or the first criteria
        if (sol->criteria_2[0].second > max_first) {
            // Keep solution if is inferior
            if (sol->inferior) {
                saved_solutions.push_back(sol);
            } else {
                max_first = sol->criteria_2[0].second;
            }
        } else if (sol->inferior) { // Free the memory pointed to by sol because it is dominated
            delete(sol);
        }
    }

    return saved_solutions;
}


void DP_Algorithm::run_experiment() {
    // Run the experiment
    // Monitor wall and CPU time
    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    clock_t t = clock();
    if (use_binary_tree) {
        build_DP(true);
    } else {
        build_DP_table_recursive(net.root);
    }
    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    t = clock() - t;

    // Output times!
    std::cout << "Time generating: " << time_generating << std::endl;
    std::cout << "Time copying: " << time_copying << std::endl;
    std::cout << "Time sorting: " << time_sorting << std::endl;
    std::cout << "Time nlogn: " << time_nlogn << std::endl;

    // Create file name
    std::ostringstream create_file;
    // Append tree structure
    if (use_binary_tree) {
        if (net.balance_tree) {
            create_file << "binary-balanced_";
        } else {
            create_file << "binary-un-balanced_";
        }
        if (net.split_values) {
            create_file << "split_";
        } else {
            create_file << "zero_";
        }
    } else {
        create_file << "directed_";
    }
    // Append algorithm type
    if (use_n_squared) {
        create_file << "nsquared_";
    } else {
        create_file << "nlogn_";
    }
    // Add the river network
    create_file << net.river_network;
    for (const std::string& str: relevant_criteria) {
        create_file << "_" << str;
    }
    // Add the epsilon
    create_file << "_e_" << k_epsilon;

    create_file << "_K_" << K;

    // Add the batching value
    create_file << "_b_" << batch_size;

    // Get the time
    time_t rawtime;
    struct tm * timeinfo;
    char buffer [80];

    time (&rawtime);
    timeinfo = localtime (&rawtime);

    strftime (buffer,80,"%c",timeinfo);
    std::string time_stamp = buffer;
    // Get rid of the spaces
    std::replace(time_stamp.begin(), time_stamp.end(), ' ', '_');
    std::replace(time_stamp.begin(), time_stamp.end(), ':', '_');
    create_file << "_" << time_stamp << ".sol";

    // Print the result of the trail
    std::ofstream output(create_file.str());
    if (output.is_open()) {
        // Print date and time for trial
        output << "Date/time: " << time_stamp << "\n";

        // File used
        output << "Data file: " << net.filename << "\n";

        // Set the precision
        std::ios::fmtflags old_settings = output.flags(); //save previous format flags
        long old_precision = output.precision();

        // NEW
        std::chrono::duration<float> duration = t2 - t1;
        output << std::fixed << std::setprecision(4) << "Wall time: " << duration.count() << " seconds.\n";

        output << std::fixed << std::setprecision(4) << "CPU time: " << ((double)t)/CLOCKS_PER_SEC << " seconds.\n";


        // Output to exp file time
        std::cout << "Wall time experiment: " << duration.count() << std::endl;
        std::cout << "CPU time experiment: " << ((double)t)/CLOCKS_PER_SEC << std::endl;


        //Restore precision
        output.flags(old_settings);
        output.precision(old_precision);
        // Print the seed
        output << "seed: " << this->random_seed << "\n";

        print_DP_Output(output);
        output.close();
    }
}