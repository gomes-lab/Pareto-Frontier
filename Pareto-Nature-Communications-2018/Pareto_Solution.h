//
// Created by Jonathan Gomes Selman on 8/6/17.
//

#ifndef AMAZON_PROJECT_PARETO_SOLUTION_H
#define AMAZON_PROJECT_PARETO_SOLUTION_H


#include <vector>
#include "string"
#include <functional>
#include <cmath>

// Gives the maximum number of criteria that could be used in an experimental run.
// used to determine the size of the static array of criteria values
#define  MAX_CRITERIA 2

class Pareto_Solution {
public:
    // Class variables
    // --------------------------------------------------------------------------------
    int node_id;

    // Static array representation of the criteria values for a given solution.
    // A static array is initialized to the size of MAX_CRITERIA, even though
    // the number of filled indexes is equal to the number of criteria being
    // considered for a given experimental run. A static array is used to avoid
    // the overhead of dynamically allocated memory on the heap that is needed
    // when working with std::vector. The static array is allocated on the stack
    // and allows for much quicker allocation and data retrieval / modification.
    std::pair<double, double> criteria_2[MAX_CRITERIA];

    std::vector<int> dam_decisions; // Array of 0 and 1 showing if each dam was built or not

    // Represents an array of pointers to the Pareto_Solutions that was used for each
    // child to generate the new Pareto_Solution
    std::vector<Pareto_Solution*> pareto_decisions;

    // End class variables
    // --------------------------------------------------------------------------------




    /*
     * Constructor used for creating a solution utilizing the static array for storing criterion values
     */
    Pareto_Solution(int node_id, const std::pair<double, double> criteria[], const std::vector<int> &dam_decisions,
                    const std::vector<Pareto_Solution *> &pareto_decisions, int num_criteria);

    /*
     * Compares solutions based on the criteria values. Specifically, compares
     * the 'rounded_value' for each criteria. Returns true is the compared solution
     * has equal values for each rounded criterion value.
     */
    bool operator==(const Pareto_Solution &rhs) const;

    /*
     * Uses equal operator to evaluate not-equals
     */
    bool operator!=(const Pareto_Solution &rhs) const;

    /*
     * Gets the number of criteria that is actually being considered
     */
    int getNum_criteria() const;

private:
    int num_criteria;

};

/*
 * defines the hash protocol for a Pareto_Solution pointer
 */
namespace std
{
    template <>
    struct hash<Pareto_Solution* >
    {
        size_t operator()( const Pareto_Solution* sol ) const
        {
            // Compute individual hash values each criterion's rounded value
            // Hash approach courtesy of ---
            // http://stackoverflow.com/a/1646913/126995
            size_t res = 17;
            for (int i = 0; i < sol->getNum_criteria(); i++) {
                res = res * 31 + std::hash<double>{}(sol->criteria_2[i].second);
            }
            return res;
        }
    };
}


#endif //AMAZON_PROJECT_PARETO_SOLUTION_H
