//
// Created by Jonathan Gomes Selman on 8/6/17.
//

#include "Pareto_Solution.h"



Pareto_Solution::Pareto_Solution(int node_id, const std::pair<double, double> *criteria, const std::vector<int> &dam_decisions,
                                 const std::vector<Pareto_Solution *> &pareto_decisions,
                                 int num_criteria) :  node_id(node_id), dam_decisions(dam_decisions),
                                                      pareto_decisions(pareto_decisions), num_criteria(num_criteria){
    // Copy in criteria values
    for (int i = 0; i < num_criteria; i++) {
        this->criteria_2[i] = criteria[i];
    }

}



bool Pareto_Solution::operator==(const Pareto_Solution &rhs) const {
    // Compare the 'rounded-value' for each criterion
    for (int i = 0; i < num_criteria; i++) {
        if (this->criteria_2[i].second != rhs.criteria_2[i].second) { // Watch out for doubles equals!!
            return false;
        }
    }
    return true;
}

bool Pareto_Solution::operator!=(const Pareto_Solution &rhs) const {
    return !(rhs == *this);
}

int Pareto_Solution::getNum_criteria() const {
    return num_criteria;
}









