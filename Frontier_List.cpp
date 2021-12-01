//
// Created by Jonathan Gomes Selman on 8/6/17.
//

#include "Frontier_List.h"

Frontier_List::Frontier_List() {
    // Create a dummy header for simplicity of traversal
    // and for simplicity of adding and deleting values.
    head = new Pareto_Solution(-1, {}, {}, {});
    num_solutions = 0;
}

Frontier_List::~Frontier_List() {
    destructor_helper();
}

void Frontier_List::destructor_helper() {
    Pareto_Solution* curr = head;
    // Free memory associated with each solution
    while (curr != nullptr) {
        Pareto_Solution* trash = curr;
        curr = curr->next;
        delete  trash;
    }
}

Frontier_List::Frontier_List(const Frontier_List &src) {
    // We need to create a copy of each node
    this->head = new Pareto_Solution(-1, {}, {}, {});
    Pareto_Solution *copy = this->head;
    for (Pareto_Solution *curr = src.get_head(); curr != nullptr; curr = curr->next) {
        // Create a new copy of the node - utilize default copy for node as all the member variable have copy functions
        copy->next = new Pareto_Solution(*curr);
        copy = copy->next;
    }
    this->num_solutions = src.num_solutions;
}

Frontier_List &Frontier_List::operator=(const Frontier_List &src) {
    // Generate new list
    Pareto_Solution *tmpHead = new Pareto_Solution(-1, {}, {}, {});
    Pareto_Solution *copy = tmpHead;
    for (Pareto_Solution *curr = src.get_head(); curr != nullptr; curr = curr->next) {
        copy->next = new Pareto_Solution(*curr);
        copy = copy->next;
    }

    // Delete old list
    destructor_helper();

    // Set the head to the new list
    this->head = tmpHead;
    this->num_solutions = src.num_solutions;
    return *this;
}

Pareto_Solution *Frontier_List::get_head() const{
    // We want to skip over the dummy head
    return this->head->next;
}

void Frontier_List::add_n_squared(const std::vector<std::pair<double, double> > &criteria, int nodeID,
                                  const std::vector<int> &dam_decisions,
                                  const std::vector<Pareto_Solution *> &pareto_decisions) {
    Pareto_Solution *add_node = new Pareto_Solution(nodeID, criteria, dam_decisions, pareto_decisions);
    Pareto_Solution *curr = this->head->next;
    Pareto_Solution *prev = this->head;

    while (curr != nullptr) {
        int compare = add_node->is_dominant(curr);
        if (compare == 2) { // They are the same so random decide which to keep --- May only want to do this when rounding!!!!!!
            // Generate randomly distributed number
            int rand_num = rand() % 2; // Only want to do if epsilon != 0 i.e. we are rounding!!!!

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
            num_solutions--;
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
    num_solutions++;
    prev->next = add_node;

}

int Frontier_List::getNum_solutions() const {
    return num_solutions;
}

void Frontier_List::add_n_squared(const std::pair<double, double> *criteria, int nodeID,
                                  const std::vector<int> &dam_decisions,
                                  const std::vector<Pareto_Solution *> &pareto_decisions, int num_solutions,
                                  bool use_strict_compare) {
    Pareto_Solution *add_node = new Pareto_Solution(nodeID, criteria, dam_decisions, pareto_decisions, num_solutions);
    Pareto_Solution *curr = this->head->next;
    Pareto_Solution *prev = this->head;


    while (curr != nullptr) {
        int compare;
        if (use_strict_compare) {
            compare = add_node->is_dominant_strong(curr);
        } else {
            compare = add_node->is_dominant_2(curr);
        }

        if (compare == 2) { // They are the same so random decide which to keep --- May only want to do this when rounding!!!!!!
            int rand_num = rand() % 2; // Only want to do if epsilon != 0 i.e. we are rounding!!!!

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
            this->num_solutions--;
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
    this->num_solutions++;
    prev->next = add_node;
}

void Frontier_List::add_n_squared_created(Pareto_Solution *add_node, bool strict_compare) {
    Pareto_Solution *curr = this->head->next;
    Pareto_Solution *prev = this->head;


    while (curr != nullptr) {
        int compare;
        // Will fix later
        if (strict_compare) {
            compare = add_node->is_dominant_true(curr);
        } else {
            compare = add_node->is_dominant_true(curr);
        }

        if (compare == 2) { // They are the same so random decide which to keep --- May only want to do this when rounding!!!!!!
            int rand_num = rand() % 2; // Only want to do if epsilon != 0 i.e. we are rounding!!!!

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
            this->num_solutions--;
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
    this->num_solutions++;
    prev->next = add_node;
}


