//
// Created by Jonathan Gomes Selman on 7/6/17.
//

#ifndef AMAZON_PROJECT_PARETO_OPT_NODE_H
#define AMAZON_PROJECT_PARETO_OPT_NODE_H


#include <vector>


class Pareto_Opt_Node {
public:
    // Values will be likely represented differently later on --- PROBABLY IN A DICTIONARY!
    // These values are the true connectivity for a given Pareto-pair
    double connectivity;
    double energy;
    double sediment;
    int nodeID;
    //HyperTreeNode* node;

    // These values are the rounded values that are used to compare Pareto-pairs
    double connectivity_compare;
    double energy_compare;
    double sediment_compare;

    std::vector<int> dam_decisions; // Array of 0 and 1 showing if each dam was built or not

    bool inferior;

    // Represents an array of pointers to the Pareto_Opt_Node that was used for each
    // child to generate the new Pareto_Opt_Node
    std::vector<Pareto_Opt_Node*> pareto_decisions;

    Pareto_Opt_Node* next; // Pointer to next node in the list

    Pareto_Opt_Node(double connectivity, double energy, int nodeID, double con_compare, double energy_compare,
                    std::vector<int> dam_decisions,
                    std::vector<Pareto_Opt_Node* > pareto_decisions);

    Pareto_Opt_Node(double connectivity, double energy, double sediment, double con_compare, double energy_compare,
                    double sediment_compare, int id,
                    std::vector<int> dam_decisions,
                    std::vector<Pareto_Opt_Node* > pareto_decisions);

    /*Pareto_Opt_Node(double connectivity, double energy, HyperTreeNode* node, double con_compare, double energy_compare,
                    std::vector<int> dam_decisions,
                    std::vector<Pareto_Opt_Node* > pareto_decisions);*/

    /*bool operator<(const Pareto_Opt_Node &rhs) const;

    bool operator>(const Pareto_Opt_Node &rhs) const;

    bool operator<=(const Pareto_Opt_Node &rhs) const;

    bool operator>=(const Pareto_Opt_Node &rhs) const;

    bool operator==(const Pareto_Opt_Node &rhs) const;

    bool operator!=(const Pareto_Opt_Node &rhs) const;*/


    int is_dominant(Pareto_Opt_Node* compare);
};


#endif //AMAZON_PROJECT_PARETO_OPT_NODE_H
