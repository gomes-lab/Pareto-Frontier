//
// Created by Jonathan Gomes Selman on 7/6/17.
//

#include "Pareto_Opt_Node.h"

Pareto_Opt_Node::Pareto_Opt_Node(double connectivity, double energy, int nodeID, double con_compare, double energy_compare,
                                 std::vector<int> dam_decisions,
                                 std::vector<Pareto_Opt_Node* > pareto_decisions) {
    this->connectivity = connectivity;
    this->energy = energy;
    //this->dam_decisions = dam_decisions;
    this->dam_decisions = std::move(dam_decisions);
    //this->pareto_decisions = pareto_decisions;
    this->pareto_decisions = std::move(pareto_decisions);
    this->next = nullptr;
    this->nodeID = nodeID;
    this->connectivity_compare = con_compare;
    this->energy_compare = energy_compare;
}

/*Pareto_Opt_Node::Pareto_Opt_Node(double connectivity, double energy, HyperTreeNode *node, double con_compare,
                                 double energy_compare, std::vector<int> dam_decisions,
                                 std::vector<Pareto_Opt_Node *> pareto_decisions) {
    this->next = nullptr;
    this->dam_decisions = std::move(dam_decisions);
    this->pareto_decisions = std::move(pareto_decisions); // Trying out the move operator, because the copied values are only used once
    this->connectivity = connectivity;
    this->energy = energy;
    this->node = node;
    this->connectivity_compare = con_compare;
    this->energy_compare = energy_compare;
}*/

Pareto_Opt_Node::Pareto_Opt_Node(double connectivity, double energy, double sediment, double con_compare,
                                 double energy_compare, double sediment_compare, int id,
                                 std::vector<int> dam_decisions, std::vector<Pareto_Opt_Node *> pareto_decisions) {
    this->next = nullptr;
    //this->dam_decisions = dam_decisions;
    this->dam_decisions = std::move(dam_decisions);
    //this->pareto_decisions = pareto_decisions;
    this->pareto_decisions = std::move(pareto_decisions); // Trying out the move operator, because the copied values are only used once
    this->connectivity = connectivity;
    this->energy = energy;
    this->nodeID = id;
    this->connectivity_compare = con_compare;
    this->energy_compare = energy_compare;
    this->sediment = sediment;
    this->sediment_compare = sediment_compare;

}

int Pareto_Opt_Node::is_dominant(Pareto_Opt_Node *compare) {
    bool dominate = true;
    bool is_dominated = true;
    // Compare connectivity
    if (connectivity_compare < compare->connectivity_compare) {
        dominate = false;
    } else if (connectivity_compare > compare->connectivity_compare) {
        is_dominated = false;
    }

    // Compare energy
    if (energy_compare < compare->energy_compare) {
        dominate = false;
    } else if (energy_compare > compare->energy_compare) {
        is_dominated = false;
    }

    // Compare sediment
    if (sediment_compare < compare->sediment_compare) {
        dominate = false;
    } else if (sediment_compare > compare->sediment_compare) {
        is_dominated = false;
    }

    if (dominate && is_dominated) { // If is both dominating and dominated then they are the same
        return 2;
    } else if(!dominate && !is_dominated) { // Neither dominates
        return 0;
    } else if (dominate) {
        return 1; // Data dominates the compare data
    } else {
        return -1; // Is dominated
    }
}

/*bool Pareto_Opt_Node::operator==(const Pareto_Opt_Node &rhs) const {
    return std::tie(connectivity_compare, energy_compare, sediment_compare) ==
           std::tie(rhs.connectivity_compare, rhs.energy_compare, rhs.sediment_compare);
}

bool Pareto_Opt_Node::operator!=(const Pareto_Opt_Node &rhs) const {
    return !(rhs == *this);
}

bool Pareto_Opt_Node::operator<(const Pareto_Opt_Node &rhs) const {
    return std::tie(connectivity_compare, energy_compare, sediment_compare) <
           std::tie(rhs.connectivity_compare, rhs.energy_compare, rhs.sediment_compare);
}

bool Pareto_Opt_Node::operator>(const Pareto_Opt_Node &rhs) const {
    return rhs < *this;
}

bool Pareto_Opt_Node::operator<=(const Pareto_Opt_Node &rhs) const {
    return !(rhs < *this);
}

bool Pareto_Opt_Node::operator>=(const Pareto_Opt_Node &rhs) const {
    return !(*this < rhs);
}*/

