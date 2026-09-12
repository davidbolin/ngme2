#ifndef NGME_TRANSFORM_UTILS_H
#define NGME_TRANSFORM_UTILS_H

#include <string>
#include <cmath>

namespace ngme {
namespace transforms {

/**
 * Apply a named transformation to a value
 * @param value The value to transform
 * @param trans_type The type of transformation to apply
 * @return The transformed value
 */
inline double apply_transform(double value, const std::string& trans_type) {
    if (trans_type == "exp4") {
        return std::exp(4 * value);
    } else if (trans_type == "exp2") {
        return std::exp(2 * value);
    } else if (trans_type == "tanh") {
        return (-1 + (2 * std::exp(value)) / (1 + std::exp(value)));
    } else if (trans_type == "sech") {
        // sqrt(1 - tanh(value)^2) for the "tanh" above, i.e. sech(value/2).
        // This is the AR(1) stationary standard deviation sqrt(1 - rho^2):
        // it has to be its own coefficient because the generic operator builds
        // K as a linear combination of fixed matrices. cosh overflows to inf
        // for large |value|, giving 0, which is the correct limit as rho -> +-1.
        return 1.0 / std::cosh(0.5 * value);
    } else if (trans_type == "identity") {
        return value;
    } else if (trans_type == "exp") {
        return std::exp(value);
    } else if (trans_type == "sqrt") {
        return std::sqrt(value);
    } else if (trans_type == "square") {
        return value * value;
    } else if (trans_type == "log") {
        return std::log(value);
    } else if (trans_type == "null") {
        return 1;
    } else {
        // Default to identity for unknown transformations
        return value;
    }
}

/**
 * First and second derivatives of apply_transform with respect to `value`.
 * Kept beside the transform itself so the two cannot drift apart.
 * @param order 1 or 2
 */
inline double transform_derivative(double value, const std::string& trans_type,
                                   int order) {
    if (trans_type == "exp4") {
        return (order == 1 ? 4.0 : 16.0) * std::exp(4 * value);
    } else if (trans_type == "exp2") {
        return (order == 1 ? 2.0 : 4.0) * std::exp(2 * value);
    } else if (trans_type == "tanh") {
        // t(v) = tanh(v/2):  t' = (1 - t^2)/2,  t'' = -t (1 - t^2)/2
        const double t = -1 + (2 * std::exp(value)) / (1 + std::exp(value));
        return order == 1 ? 0.5 * (1 - t * t) : -0.5 * t * (1 - t * t);
    } else if (trans_type == "sech") {
        // s(v) = sech(v/2):  s' = -s t / 2,  s'' = s (2 t^2 - 1) / 4,
        // with t = tanh(v/2). cosh overflows to inf for large |v|, giving 0,
        // which is the correct limit -- and the derivatives vanish there too.
        const double s = 1.0 / std::cosh(0.5 * value);
        const double t = std::tanh(0.5 * value);
        return order == 1 ? -0.5 * s * t : 0.25 * s * (2 * t * t - 1);
    } else if (trans_type == "identity") {
        return order == 1 ? 1.0 : 0.0;
    } else if (trans_type == "exp") {
        return std::exp(value);
    } else if (trans_type == "sqrt") {
        return order == 1 ? 0.5 / std::sqrt(value)
                          : -0.25 / (value * std::sqrt(value));
    } else if (trans_type == "square") {
        return order == 1 ? 2.0 * value : 2.0;
    } else if (trans_type == "log") {
        return order == 1 ? 1.0 / value : -1.0 / (value * value);
    } else if (trans_type == "null") {
        return 0.0; // contributes a constant 1
    } else {
        return order == 1 ? 1.0 : 0.0; // default identity, as above
    }
}

} // namespace transforms
} // namespace ngme

#endif // NGME_TRANSFORM_UTILS_H 