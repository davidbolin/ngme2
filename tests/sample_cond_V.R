
# Inverse Gaussian (results in NIG distribution) - AR + FE + RE (it does not use h)

Eigen::VectorXd temporal = mu_gibbs.cwiseProduct(mu_gibbs)/(sigma*sigma);
Eigen::VectorXd arg_1;
arg_1.setOnes(m);
arg_1 = -arg_1;

Eigen::VectorXd arg_2 = temporal + eta * Eigen::VectorXd::Ones(temporal.rows());

Eigen::VectorXd arg_3 = eta + (K*new_W+h.cwiseProduct(mu_gibbs)).cwiseProduct((K*new_W+h.cwiseProduct(mu_gibbs)))/(sigma*sigma);

Eigen::MatrixXd gig_Gen = rGIG_cpp(arg_1, arg_2, arg_3);

V = gig_Gen;

# Gamma (results in GAL distribution) - AR + FE + RE (it does not use h)

Eigen::VectorXd temporal = mu_gibbs.cwiseProduct(mu_gibbs)/(sigma*sigma);
Eigen::VectorXd arg_1;
arg_1.setOnes(m);
arg_1 = eta-0.5*arg_1;

Eigen::VectorXd arg_2 = temporal + 2*eta * Eigen::VectorXd::Ones(temporal.rows());

Eigen::VectorXd arg_3 = (K*new_W+h.cwiseProduct(mu_gibbs)).cwiseProduct((K*new_W+h.cwiseProduct(mu_gibbs)))/(sigma*sigma);

Eigen::MatrixXd gig_Gen = rGIG_cpp(arg_1, arg_2, arg_3);

V = gig_Gen;


### Matern case

# Inverse Gaussian (results in NIG distribution) - Matern (it uses h)

Eigen::VectorXd temporal = mu_gibbs.cwiseProduct(mu_gibbs)/(sigma*sigma);
Eigen::VectorXd arg_1;
arg_1.setOnes(m);
arg_1 = -arg_1;

Eigen::VectorXd arg_2 = temporal + eta * Eigen::VectorXd::Ones(temporal.rows());

Eigen::VectorXd arg_3 = eta*(h.cwiseProduct(h)) + (K*new_W+h.cwiseProduct(mu_gibbs)).cwiseProduct((K*new_W+h.cwiseProduct(mu_gibbs)))/(sigma*sigma);

Eigen::MatrixXd gig_Gen = rGIG_cpp(arg_1, arg_2, arg_3);

V = gig_Gen;

# Gamma (results in GAL distribution) - Matern (it uses h)

Eigen::VectorXd temporal = mu_gibbs.cwiseProduct(mu_gibbs)/(sigma*sigma);
Eigen::VectorXd arg_1;
arg_1.setOnes(m);
arg_1 = eta*h-0.5*arg_1;

Eigen::VectorXd arg_2 = temporal + 2*eta * Eigen::VectorXd::Ones(temporal.rows());

Eigen::VectorXd arg_3 = (K*new_W+h.cwiseProduct(mu_gibbs)).cwiseProduct((K*new_W+h.cwiseProduct(mu_gibbs)))/(sigma*sigma);

Eigen::MatrixXd gig_Gen = rGIG_cpp(arg_1, arg_2, arg_3);

V = gig_Gen;
