//
// Created by Jonathan Gomes Selman on 7/6/17.
//

#ifndef AMAZON_PROJECT_PARETO_OPT_LIST_H
#define AMAZON_PROJECT_PARETO_OPT_LIST_H

#include "Pareto_Opt_Node.h"
#include <iostream>
#include <algorithm>
#include <cmath>
#include "time.h"
#include <random>

class Pareto_Opt_List {
public:
    Pareto_Opt_List(); // Constructor
    virtual ~Pareto_Opt_List(); // Destructor

    /*
     * Follow the rule of 3:
     *      Destructor
     *      Copy constructor
     *      Assignment operator
     */
    void destructor_helper();
    Pareto_Opt_List(const Pareto_Opt_List& src); // Copy constructor
    Pareto_Opt_List& operator=(const Pareto_Opt_List& src); // assignment operator

    /*
     * Add a new Pareto-pair to the Pareto List
     */
    void add(double new_connectivity, double new_energy, int nodeID,
             std::vector<int> dam_decisions,
             std::vector<Pareto_Opt_Node* > pareto_decisions); // Note the representation of the variables will change


    /*void add(double new_connectivity, double new_energy, HyperTreeNode *node,
             std::vector<int> dam_decisions,
             std::vector<Pareto_Opt_Node *> pareto_decisions);*/

    /**
     * Add a new Pareto-pair to the list with rounded values
     * @param new_connectivity the true connectivity for the given node
     * @param new_energy the true energy for the given node
     * @param con_round rounded connectivity used for comparison
     * @param energy_round rounded energy
     */
    void add_round(double new_connectivity, double new_energy, double con_round, double energy_round, int nodeID,
                   std::vector<int> dam_decisions,
                   std::vector<Pareto_Opt_Node* > pareto_decisions);

    /*void add_round(double new_connectivity, double new_energy, double con_round, double energy_round, HyperTreeNode* node,
                   std::vector<int> dam_decisions,
                   std::vector<Pareto_Opt_Node* > pareto_decisions)*/

    void add_multi(double true_con, double true_energy, double true_sed, double round_con, double round_energy,
                    double round_sed, int nodeID, std::vector<int> dam_decisions,
                   std::vector<Pareto_Opt_Node* > pareto_decisions);


    Pareto_Opt_Node * getHead()const; // Get the start of the list so we can loop


private:
    Pareto_Opt_Node* head;

    // Used for random generation
    //std::random_device rd;


    /**
     * Helper method to add a new pareto node
     * @param con_compare connectivity that will be used to compare nodes
     * @param energy_compare energy that will be used to compare nodes
     * @param con_orig true value
     * @param energy_orig true value
     */
    void add_helper(double con_compare, double energy_compare, double con_orig, double energy_orig, int nodeID,
                    std::vector<int> &dam_decisions,
                    std::vector<Pareto_Opt_Node* > &pareto_decisions);

    /*void add_helper(double con_compare, double energy_compare, double con_orig, double energy_orig, HyperTreeNode* node,
                    std::vector<int> &dam_decisions,
                    std::vector<Pareto_Opt_Node* > &pareto_decisions);*/
};


#endif //AMAZON_PROJECT_PARETO_OPT_LIST_H
