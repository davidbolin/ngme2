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

} // namespace transforms
} // namespace ngme

#endif // NGME_TRANSFORM_UTILS_H 