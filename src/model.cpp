// #include "model.h"

// // Property of K
// VectorXd Model::grad_mu_1(const VectorXd& w,
//                           const VectorXd& V,
//                           const VectorXd& Y) 
// {
//     // VectorXd m = K.selfadjointView<Eigen::Upper>().llt().solve(-mu * h + mu * V);
    
//     VectorXd m = K.partialPivLu().solve(-mu * h + mu * V);
//     VectorXd r = pow(sigma_eps, -2) * A * K.partialPivLu().solve(-h + V) * 
//                    V.cwiseInverse().asDiagonal() * 
//                    (Y - A * K.partialPivLu().solve(-mu * h + mu * V) - A * (w - m));
//     return r;
// }

// double Model::grad_sgima_eps_1(const VectorXd& w,
//                              const VectorXd& V,
//                              const VectorXd& Y) 
// {
//     VectorXd m = K.partialPivLu().solve(-mu * h + mu * V);
//     VectorXd x = Y - A * K.partialPivLu().solve(-mu * h + mu * V) - A * (w - m);
//     double r = -n / sigma_eps + 1 / pow(sigma_eps, 3) * 
//                  x.transpose() * V.cwiseInverse().asDiagonal() * x;

//     return r;
// }

// // double Model::grad_kappa_sq_1(const VectorXd& w,
// //                             const VectorXd& V,
// //                             const VectorXd& Y) 
// // {
// //     double r = (C * K.inverse()).trace() - w.transpose() * C.transpose() * 
// //                 V.cwiseInverse().asDiagonal() * 
// //                 (K*w + (h-V) * mu);
// //     return r;
// // }

// // grad of mu for all samples V and Y
// // d * num_samples
// VectorXd Model::grad_mu(const MatrixXd& ws,
//                         const MatrixXd& Vs,
//                         const VectorXd& Y) 
// {
//     int samples = ws.cols();
//     MatrixXd grads (A.rows(), samples);

//     // compute grad for every sample
//     for (int i = 1; i < samples; i++) {
//         grads.col(i) = grad_mu_1(ws.col(i), Vs.col(i), Y);
//     }

//     return grads.rowwise().mean();
// }

// double Model::grad_sgima_eps(const MatrixXd& ws,
//                      const MatrixXd& Vs,
//                      const VectorXd& Y)
// {
//     int samples = ws.cols();
//     VectorXd grads (samples);
//     for (int i = 1; i < samples; i++) {
//         grads(i) = grad_sgima_eps_1(ws.col(i), Vs.col(i), Y);
//     }
//     return grads.mean();
// }

// double Model::grad_kappa_sq(const MatrixXd& ws,
//                      const MatrixXd& Vs,
//                      const VectorXd& Y)
// {
//     int samples = ws.cols();
//     VectorXd grads (samples);
//     for (int i = 1; i < samples; i++) {
//         grads(i) = grad_kappa_sq_1(ws.col(i), Vs.col(i), Y);
//     }
//     return grads.mean();
// }