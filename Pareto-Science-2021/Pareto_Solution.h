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
#define  MAX_CRITERIA 8

class Pareto_Solution {
public:
    // Class variables
    int node_id;

    // Holds the values of the criteria for a given solution.
    // The order of the criteria is determined by the order given
    // by the HyperNet.
    // Ex. if the algorithm is running on the criteria 'conn' and 'energy'
    //      criteria[0] --> connectivity
    //      criteria[1] --> energy
    // Note: std::pair is used to hold the true value and rounded value of a criteria.
    // The data is stored as <true_value, rounded_value>. If we choose not to round then
    // rounded_value is just equal to the true_value.
    std::vector<std::pair<double, double> > criteria; // NOTE: to allow generic types we should use std::variant

    // Static array representation of the criteria values for a given solution.
    // A static array is initialized to the size of MAX_CRITERIA, even though
    // the number of filled indexes is equal to the number of criteria being
    // considered for a given experimental run. A static array is used to avoid
    // the overhead of dynamically allocated memory on the heap that is needed
    // when working with std::vector. The static array is allocated on the stack
    // and allows for much quicker allocation and data retrieval / modification.
    std::pair<double, double> criteria_2[MAX_CRITERIA];

    std::vector<int> dam_decisions; // Array of 0 and 1 showing if each dam was built or not

    bool inferior; // Used in the divide-and-conquer algorithm. Signals if the solution is in the inferior or superior category

    // Represents an array of pointers to the Pareto_Solutions that was used for each
    // child to generate the new Pareto_Solution
    std::vector<Pareto_Solution*> pareto_decisions;

    Pareto_Solution* next; // Pointer to the next node in the linked list --- Only used with linked list data structure.

    //std::vector<std::string>* relevant_criteria; // Likely want to use this later to check whether certain criteria are min or maximized
    // End class variables



    /*
     * Constructor used for creating a solution utalizing the std::vector data structure
     * for storing criterion values
     */
    Pareto_Solution(int node_id, const std::vector<std::pair<double, double>> &criteria,
                    const std::vector<int> &dam_decisions, const std::vector<Pareto_Solution *> &pareto_decisions);

    /*
     * Constructor used for creating a solution utalizing the static array for storing criterion values
     */
    Pareto_Solution(int node_id, const std::pair<double, double> criteria[], const std::vector<int> &dam_decisions,
                    const std::vector<Pareto_Solution *> &pareto_decisions, int num_criteria);

    /*
     * Copy constructor --- Copy all data except pointer to next
     */
    Pareto_Solution(const Pareto_Solution& src);

    /*
     * Used with the vector implementation for storing criteria
     * Compares to Pareto_Solutions to see if one dominates the other.
     * Specifically returns:
     *      1 --> if the current solutions dominates 'compare'
     *      2 --> if the solutions are equal
     *      0 --> if the solutions do not dominate eachother
     *      -1 --> if the current solution is dominated by 'compare'
     */
    int is_dominant(Pareto_Solution* compare);

    /*
     * Used with the static array implementation for storing criteria
     * See documentation for 'is_dominant' to see implementation details
     */
    int is_dominant_2(Pareto_Solution* compare);

    /*
     * See if the current solution is strictly dominated by compare.
     * Additional flag (3) returned if solutions are equal in one category.
     *
     */
    int is_dominant_strict(Pareto_Solution* compare);

    /*
     * If ties exist between rounded criteria values, the true
     * solution values are compared.
     */
    int is_dominant_strong(Pareto_Solution* compare);

    /*
     * Compares non-rounded values of solutions
     */
    int is_dominant_true(Pareto_Solution* compare);

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
    static constexpr double EQUALITY_EPSILON = 0.0001;

};

/*
 * defines the hash protocal for a Pareto_Solution object
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
