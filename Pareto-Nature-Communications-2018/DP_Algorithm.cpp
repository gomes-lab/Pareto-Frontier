//
// Created by Jonathan Gomes Selman on 7/6/17.
//

#include "DP_Algorithm.h"

DP_Algorithm::DP_Algorithm(const HyperNet &net, int root, bool maximize_c1, bool maximize_c2, double epsilon, unsigned int seed, int batch_size, size_t num_threads)
        : net(net), root(root), maximize_c1(maximize_c1), maximize_c2(maximize_c2), k_epsilon(epsilon),
          random_seed(seed), batch_size(batch_size), K(num_threads),pool(K){

    // Initialize counters
    this->total_pruned_policies = 0;
    this->total_policies_considered = 0;
    this->max_policies_considered = 0;

    // Initialize timers
    //this->time_copying = 0;
    //this->time_generating = 0;
    //this->time_sorting = 0;
    //this->time_nlogn = 0;

    // Calculate rounding constants
    calc_theoretical_ks();
    calc_theoretical_ks_tree(this->net.root_node);

    for (int i = 0; i < K; ++i){
        solution_sets.push_back(std::unordered_set<Pareto_Solution *, SolutionHash, EqualSolutions>());
    }

}

void DP_Algorithm::calc_theoretical_ks() {
    for (auto &node : net.nodes) {

        node.k_c1 = k_epsilon * node.c1;
        node.k_c2 = k_epsilon * node.c2;

    }
}

void DP_Algorithm::calc_theoretical_ks_tree(HyperTreeNode* node) {
    if (node != nullptr) {

        node->node_data.k_c1 = k_epsilon * node->node_data.c1;
        node->node_data.k_c2 = k_epsilon * node->node_data.c2;

        // Recurse
        calc_theoretical_ks_tree(node->left);
        calc_theoretical_ks_tree(node->right);
    }
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

        node->vec_frontier.push_back(new Pareto_Solution(node->node_data.id, criteria, {}, {}, 2));
    } else { // Not leaf node

        if (node->left != nullptr && node->right != nullptr) { // Has both children
            // Calculate the number that will be considered
            num_considered = node->left->vec_frontier.size() * node->right->vec_frontier.size() *
                                    node->leftD.decision.size() * node->rightD.decision.size();
            std::cout << "\tConsidering: " << num_considered << std::endl;
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
                                        add_nlogn_solution(pair, decision, node, solution_sets[i],num_generated,
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
                                add_nlogn_solution(pair, decision, node, all_solution_set, num_generated,
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
            std::cout << "\tConsidering: " << num_considered << std::endl;
            total_policies_considered += num_considered;
            max_policies_considered = std::max(max_policies_considered, num_considered);
            // Loop through left dam
            for (int left_dam_plan: node->leftD.decision) {
                std::vector<int> decision = {left_dam_plan, -1};

                for (Pareto_Solution* currL: node->left->vec_frontier) {
                    std::vector<Pareto_Solution *> pair = {currL, nullptr};

                    // Call method to generate new solution
                    add_nlogn_solution(pair, decision, node, all_solution_set, num_generated,
                                       batch_number);
                }
            }
        } else { // Note this can never happen!!!!! Because we always put child as left first
            // Track policies
            int num_generated = 0;
            num_considered = node->right->vec_frontier.size() * node->rightD.decision.size();
            std::cout << "\tConsidering: " << num_considered << std::endl;
            total_policies_considered += num_considered;
            max_policies_considered = std::max(max_policies_considered, num_considered);
            // Loop through left dam
            for (int right_dam_plan: node->rightD.decision) {
                std::vector<int> decision = {-1, right_dam_plan};

                for (Pareto_Solution *currR: node->right->vec_frontier) {
                    std::vector<Pareto_Solution *> pair = {nullptr, currR};

                    // Call method to generate new solution
                    add_nlogn_solution(pair, decision, node, all_solution_set, num_generated,
                                       batch_number);

                }
            }
        }

        // Add the solutions from the hash set to a vector.
        // No ties exist
        std::vector<Pareto_Solution* > solution_set(all_solution_set.begin(), all_solution_set.end());

        std::cout << "Num-actually considered: " << solution_set.size() << std::endl;

        // Perform the divide and conquer
        node->vec_frontier = L2D(solution_set);
    }
    std::cout << "End - kept: " << node->vec_frontier.size() << "\tdropped: " <<
              (num_considered - node->vec_frontier.size()) << std::endl << std::endl;
    total_pruned_policies += num_considered - node->vec_frontier.size();
}

void
DP_Algorithm::add_nlogn_solution(std::vector<Pareto_Solution *> &grouping, std::vector<int> &decisionVec,
                                 HyperTreeNode *node,
                                 std::unordered_set<Pareto_Solution *, SolutionHash, EqualSolutions> &solution_set,
                                 int &num_generated, int &batch_num) {
    std::pair<double, double> criteria[MAX_CRITERIA];
    initialize_criteria(node, criteria);

    // Keep track of the existing values of c1 and c2
    double c1_built = 0;
    double c2_built = 0;

    // Add the corresponding values for each of the Pareto pairs from the children
    for (int i = 0; i < decisionVec.size(); i++) {
        if (decisionVec[i] != -1) {
            update_criteria_new(node, grouping[i], criteria, decisionVec[i], i, c1_built, c2_built);
        }
    }
    // Later we should not update this everytime!
    assign_node_k(node, c1_built, c2_built);

    round_criteria(node, criteria); // Round the criteria values based on local ks

    auto adding = new Pareto_Solution(node->node_data.id, criteria, decisionVec, grouping, 2);
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

        std::vector<Pareto_Solution* > batched_solutions(solution_set.begin(), solution_set.end());

        std::vector<Pareto_Solution* > temp_answers = L2D(batched_solutions);

        solution_set = std::unordered_set<Pareto_Solution*, SolutionHash, EqualSolutions> (temp_answers.begin(), temp_answers.end());
        num_generated = 0;
        mtx.lock();
        batch_num++;
        mtx.unlock();
    }

}


void DP_Algorithm::initialize_criteria(HyperTreeNode *node, std::pair<double, double> *criteria) {
    criteria[0]={0, 0};
    criteria[1]={0, 0};

}


// This is the updated version that keeps track of the criteria we build so that we can round correctly
void DP_Algorithm::update_criteria_new(HyperTreeNode *node, Pareto_Solution *pareto_decision, std::pair<double, double> *criteria,
                                       int decision, int child_index, double &c1_built, double &c2_built) {
    if (decision == 0) { // We are not building the dam

        criteria[0].first += pareto_decision->criteria_2[0].first;
        criteria[0].second += std::abs(pareto_decision->criteria_2[0].second);

        criteria[1].first += pareto_decision->criteria_2[1].first;
        criteria[1].second += std::abs(pareto_decision->criteria_2[1].second);


    } else { // Dam is being built

        double added_c1;
        double added_c2;
        if (child_index == 0) { // Left child

            added_c1 = node->leftD.c1;
            added_c2 = node->leftD.c2;
        } else { // Right child
            added_c1 = node->rightD.c1;
            added_c2 = node->rightD.c2;
        }
        // Keep track of the existing criteria
        c1_built += added_c1;
        c2_built += added_c2;

        if (maximize_c1){
            criteria[0].first += added_c1 + pareto_decision->criteria_2[0].first;
            criteria[0].second += added_c1 + pareto_decision->criteria_2[0].second;
        } else {
            criteria[0].first += added_c1 + pareto_decision->criteria_2[0].first;
            criteria[0].second += added_c1 + std::abs(pareto_decision->criteria_2[0].second);
        }

        if (maximize_c2){
            criteria[1].first += added_c2 + pareto_decision->criteria_2[1].first;
            criteria[1].second += added_c2 + pareto_decision->criteria_2[1].second;
        } else {
            criteria[1].first += added_c2 + pareto_decision->criteria_2[1].first;
            criteria[1].second += added_c2 + std::abs(pareto_decision->criteria_2[1].second);
        }
    }
}

void DP_Algorithm::assign_node_k(HyperTreeNode *node, double c1_built, double c2_built) {
    node->node_data.k_c1 = std::floor(c1_built / net.min_c1) * net.min_c1 * k_epsilon / 2.0;
    node->node_data.k_c2 = std::floor(c2_built / net.min_c2) * net.min_c2 * k_epsilon / 2.0;
}

void DP_Algorithm::round_criteria(HyperTreeNode *node, std::pair<double, double> *criteria) {
    if (maximize_c1){
        criteria[0].second = round_with_k(criteria[0].second, node->node_data.k_c1);
    } else {
        criteria[0].second = -1 * round_with_k(criteria[0].second, node->node_data.k_c1);
    }

    if (maximize_c2){
        criteria[1].second = round_with_k(criteria[1].second, node->node_data.k_c2);
    } else {
        criteria[1].second = -1 * round_with_k(criteria[1].second, node->node_data.k_c2);
    }
}


double DP_Algorithm::round_with_k(double value, double k) {
    if (k == 0) {
        return value;
    }

    return std::floor(value / k) * k;
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

        computeDP_nlogn_batching(root);
    }
}


void DP_Algorithm::print_DP_Output(std::ostream &stream) {
    // Output num_solutions
    stream << "num_solutions: ";
    stream << net.root_node->vec_frontier.size() << "\n";

    // Give data about pruning
    stream << "# pruning steps (# nodes): " << net.num_nodes_binary << "\n";

    stream << "Max policies considered: " << max_policies_considered << "\n";
    stream << "Policies considered: " << total_policies_considered << "\n";
    stream << "Pruned policies: " << total_pruned_policies << "\n";

    // Output the epsilon
    stream << "epsilon: " << this->k_epsilon << "\n";
    // Batch
    stream << "batch size: " << this->batch_size << "\n";
    // Output the relevant criteria
    stream << "criteria: "<< this->net.c1_name << ", " << this->net.c2_name << "\n";


    // Output the info for each solution
    print_dams_generic_tree_vec(net.root_node, stream);


}


void DP_Algorithm::print_dams_generic_tree_vec(HyperTreeNode *node, std::ostream &out) {
    for (Pareto_Solution *result: node->vec_frontier) {
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
    for (int i = 0; i < 2; i++) {
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
    std::sort(all_possible.begin(), all_possible.end(), Compare_dimensions(2, 2));

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




void DP_Algorithm::run_experiment() {
    // Run the experiment
    // Monitor wall and CPU time
    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    clock_t t = clock();
    build_DP(true);

    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    t = clock() - t;

    // Create file name
    std::ostringstream create_file;

    create_file << "Two_criteria_";

    // Add the river network
    create_file << net.river_network<< "_" <<this->net.c1_name<< "_" <<this->net.c2_name;

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