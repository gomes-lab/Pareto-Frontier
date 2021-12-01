//
// Created by Jonathan Gomes Selman on 7/6/17.
//

#include "Pareto_Opt_List.h"

Pareto_Opt_List::Pareto_Opt_List() {
    // Create Dummy header
    this->head = new Pareto_Opt_Node(-1, -1, -1, -1, -1, {}, {});

    // Set random seed
    // Random number generator used to choose between pairs with same values
    srand((unsigned int)time(nullptr));
}

void Pareto_Opt_List::add(double new_connectivity, double new_energy, int nodeID, std::vector<int> dam_decisions,
                          std::vector<Pareto_Opt_Node* > pareto_decisions) {
    // Call to add helper with con_compare and energy_compare set to the true values
    add_helper(new_connectivity, new_energy, new_connectivity, new_energy, nodeID, dam_decisions, pareto_decisions);
}

/*void
Pareto_Opt_List::add(double new_connectivity, double new_energy, HyperTreeNode *node, std::vector<int> dam_decisions,
                     std::vector<Pareto_Opt_Node *> pareto_decisions) {
    add_helper(new_connectivity, new_energy, new_connectivity, new_energy, node, dam_decisions, pareto_decisions);
}*/

/*void Pareto_Opt_List::add_round(double new_connectivity, double new_energy, double con_round, double energy_round,
                                HyperTreeNode *node, std::vector<int> dam_decisions,
                                std::vector<Pareto_Opt_Node *> pareto_decisions) {
    add_helper(con_round, energy_round, new_connectivity, new_energy, node, dam_decisions, pareto_decisions);
}*/

void Pareto_Opt_List::add_round(double new_connectivity, double new_energy, double con_round, double energy_round,
                          int nodeID, std::vector<int> dam_decisions,
                          std::vector<Pareto_Opt_Node* > pareto_decisions) {
    // Call to add helper rounded values
    add_helper(con_round, energy_round, new_connectivity, new_energy, nodeID, dam_decisions, pareto_decisions);
}

Pareto_Opt_List::~Pareto_Opt_List() {
    destructor_helper();
}

void Pareto_Opt_List::destructor_helper() {
    Pareto_Opt_Node* curr = head;
    while (curr != nullptr) {
        Pareto_Opt_Node* trash = curr;
        curr = curr->next;
        delete  trash;
    }
}

/*
 * This function is primarily used for traversing the
 * list by continually getting next nodes until a nullptr
 * is reached.
 */
Pareto_Opt_Node *Pareto_Opt_List::getHead() const {
    // Skip over the Dummy header!
    return this->head->next;
}

Pareto_Opt_List::Pareto_Opt_List(const Pareto_Opt_List &src) {
    // We need to create a copy of each node
    this->head = new Pareto_Opt_Node(-1, -1, -1, -1, -1, {}, {});
    Pareto_Opt_Node *copy = this->head;
    for (Pareto_Opt_Node *curr = src.getHead(); curr != nullptr; curr = curr->next) {
        // Create a new copy of the node - utilize default copy for node as all the member variable have copy functions
        copy->next = new Pareto_Opt_Node(*curr);
        copy = copy->next;
    }
}

Pareto_Opt_List &Pareto_Opt_List::operator=(const Pareto_Opt_List &src) {
    // Generate new list
    Pareto_Opt_Node *tmpHead = new Pareto_Opt_Node(-1, -1, -1, -1, -1, {}, {});
    Pareto_Opt_Node *copy = tmpHead;
    for (Pareto_Opt_Node *curr = src.getHead(); curr != nullptr; curr = curr->next) {
        copy->next = new Pareto_Opt_Node(*curr);
        // new Pareto_Opt_Node(curr->connectivity, curr->energy);
        copy = copy->next;
    }

    // Delete old list
    destructor_helper();

    this->head = tmpHead;
    return *this;
}

void Pareto_Opt_List::add_helper(double con_compare, double energy_compare, double con_orig, double energy_orig, int nodeID,
                            std::vector<int> &dam_decisions, std::vector<Pareto_Opt_Node *> &pareto_decisions) {
    Pareto_Opt_Node *curr = this->head->next;
    Pareto_Opt_Node *prev = this->head;

    // Add the new element based on a list sorted by connectivity
    // Note this will have to change later with multiple criteria!!!!!!!!!!!!!
    while (curr != nullptr) {
        if (curr->connectivity_compare <= con_compare) {
            if (curr->energy_compare <= energy_compare) {
                if (curr->connectivity_compare == con_compare && curr->energy_compare == energy_compare) { // Choose random if they are the same
                    int rand_num = rand() % 2;
                    //std::mt19937 mt(rd());
                    //std::uniform_int_distribution<int> dist(0, 1);
                    //int rand_num = dist(mt);
                    if (rand_num == 0) { // Remove the current
                        prev->next = curr->next;
                        delete(curr);
                    } else {
                        return; // Don't added
                    }
                } else { // Remove dominated pair
                    prev->next = curr->next;
                    delete(curr);
                }
            } else if (curr->connectivity_compare == con_compare
                       && energy_compare < curr->energy_compare) {
                return; // Special case to avoid domination
            } else { // Advance the search in the list
                prev = curr;
            }
            curr = curr->next;
        } else {
            if (curr->energy_compare < energy_compare) { // Make sure the added node is not dominated
                Pareto_Opt_Node *newNode = new Pareto_Opt_Node(con_orig, energy_orig, nodeID, con_compare, energy_compare, dam_decisions, pareto_decisions);
                prev->next = newNode;
                newNode->next = curr;
            }
            return; // If we have reached this point we either added or did not and want to return
        }
    }
    // If we have reached this point then we must add at the end of the list
    prev->next = new Pareto_Opt_Node(con_orig, energy_orig, nodeID, con_compare, energy_compare, dam_decisions, pareto_decisions);
}


void
Pareto_Opt_List::add_multi(double true_con, double true_energy, double true_sed, double round_con, double round_energy,
                           double round_sed, int nodeID, std::vector<int> dam_decisions,
                           std::vector<Pareto_Opt_Node *> pareto_decisions) {
    Pareto_Opt_Node *add_node = new Pareto_Opt_Node(true_con, true_energy, true_sed, round_con, round_energy,
                                                    round_sed, nodeID, dam_decisions, pareto_decisions);
    Pareto_Opt_Node *curr = this->head->next;
    Pareto_Opt_Node *prev = this->head;

    while (curr != nullptr) {
        int compare = add_node->is_dominant(curr);
        if (compare == 2) { // They are the same so random decide which to keep
            //int rand_num = rand() % 2;

            // Generate randomly distributed number
            //std::mt19937 mt(rd());
            //std::uniform_int_distribution<int> dist(0, 1);
            //int rand_num = dist(mt);
            int rand_num = rand() % 2;

            if (rand_num == 0) { // Remove the current and replace it with the new node
                prev->next = add_node;
                add_node->next = curr->next;
                delete(curr);
            } else {
                delete(add_node);
            }
            return;
        } else if (compare == 1) { // Dominates so remove
            prev->next = curr->next;
            delete(curr);
        } else if (compare == -1) { // New is dominated so don't add
            delete(add_node);
            return;
        } else { // Non-dominant so we move forward in list
            prev = curr;
        }
        curr = curr->next;
    }
    // If we are here then its time to add
    prev->next = add_node;

}

/*void Pareto_Opt_List::add_helper(double con_compare, double energy_compare, double con_orig, double energy_orig,
                                 HyperTreeNode *node, std::vector<int> &dam_decisions,
                                 std::vector<Pareto_Opt_Node *> &pareto_decisions) {
    Pareto_Opt_Node *curr = this->head->next;
    Pareto_Opt_Node *prev = this->head;

    // Add the new element based on a list sorted by connectivity
    // Note this will have to change later with multiple criteria!!!!!!!!!!!!!
    while (curr != nullptr) {
        if (curr->connectivity_compare <= con_compare) {
            if (curr->energy_compare <= energy_compare) {
                if (curr->connectivity_compare == con_compare && curr->energy_compare == energy_compare) { // Choose random if they are the same
                    //int rand_num = rand() % 2;
                    std::mt19937 mt(rd());
                    std::uniform_int_distribution<int> dist(0, 1);
                    int rand_num = dist(mt);
                    if (rand_num == 0) { // Remove the current
                        prev->next = curr->next;
                    } else {
                        return; // Don't added
                    }
                } else { // Remove dominated pair
                    prev->next = curr->next;
                }
            } else if (curr->connectivity_compare == con_compare
                       && energy_compare < curr->energy_compare) {
                return; // Special case to avoid domination
            } else { // Advance the search in the list
                prev = curr;
            }
            curr = curr->next;
        } else {
            if (curr->energy_compare < energy_compare) { // Make sure the added node is not dominated
                Pareto_Opt_Node *newNode = new Pareto_Opt_Node(con_orig, energy_orig, node, con_compare, energy_compare, dam_decisions, pareto_decisions);
                prev->next = newNode;
                newNode->next = curr;
            }
            return; // If we have reached this point we either added or did not and want to return
        }
    }
    // If we have reached this point then we must add at the end of the list
    prev->next = new Pareto_Opt_Node(con_orig, energy_orig, node, con_compare, energy_compare, dam_decisions, pareto_decisions);
}*/







