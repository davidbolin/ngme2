#ifndef NGME_VAR_H
#define NGME_VAR_H

#include <eigen/Dense>
#include <string>

using Eigen::VectorXd;
using std::string;

class Var {
public:
    VectorXd V;
    string var_type;


};

// 
class indepent_IG : Var {

};

#endif