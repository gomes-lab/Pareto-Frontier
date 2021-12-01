//
// Created by Jonathan Gomes Selman on 8/6/17.
//

#ifndef AMAZON_PROJECT_FRONTIER_LIST_H
#define AMAZON_PROJECT_FRONTIER_LIST_H


#include "Pareto_Solution.h"
#include <random>
#include <iostream>

class Frontier_List {
public:
    Frontier_List(); // Constructor
    virtual ~Frontier_List(); // Destructor

    /*
    * Follow the rule of 3:
    *      Destructor
    *      Copy constructor
    *      Assignment operator
    */
    void destructor_helper();
    Frontier_List(const Frontier_List& src); // Copy constructor
    Frontier_List& operator=(const Frontier_List& src); // assignment operator

    /*
     * Get the first element in the list. Note, we skip over the dummy header
     */
    Pareto_Solution* get_head() const;

    /*
     * Add a new solution using std::vector to hold the new criteria values
     */
    void add_n_squared(const std::vector<std::pair<double, double> > &criteria, int nodeID,
                       const std::vector<int> &dam_decisions,
                       const std::vector<Pareto_Solution *> &pareto_decisions);

    /*
     * Add a new solution using a static array to hold the new criteria values
     */
    void add_n_squared(const std::pair<double, double> *criteria, int nodeID,
                       const std::vector<int> &dam_decisions,
                       const std::vector<Pareto_Solution *> &pareto_decisions, int num_solutions,
                       bool use_strict_compare);

    /*
     * Add a new solution that has already been constructed
     */
    void add_n_squared_created(Pareto_Solution *add_node, bool strict_compare);

    int getNum_solutions() const;

private:
    Pareto_Solution* head;
    int num_solutions;

};


#endif //AMAZON_PROJECT_FRONTIER_LIST_H
