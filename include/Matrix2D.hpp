#ifndef MATRIX2D_HPP
#define MATRIX2D_HPP

#include <vector>
#include <cstddef>
#include <cassert>
#include <algorithm>

namespace schwartz {

/**
 * @brief Contiguous row-major 2D matrix container modeled after Duffy (2004, 2022).
 * Maps 2D finite-difference mesh indices (i, j) onto a single flat std::vector<double>
 * to maximize CPU cache locality during directional ADE sweeps.
 */
class Matrix2D {
public:
    Matrix2D() : rows_(0), cols_(0) {}

    Matrix2D(std::size_t rows, std::size_t cols, double init_val = 0.0)
        : rows_(rows), cols_(cols), data_(rows * cols, init_val) {}

    inline double& operator()(std::size_t i, std::size_t j) noexcept {
        assert(i < rows_ && j < cols_);
        return data_[i * cols_ + j];
    }

    inline double operator()(std::size_t i, std::size_t j) const noexcept {
        assert(i < rows_ && j < cols_);
        return data_[i * cols_ + j];
    }

    inline std::size_t rows() const noexcept { return rows_; }
    inline std::size_t cols() const noexcept { return cols_; }

    inline void fill(double val) noexcept {
        std::fill(data_.begin(), data_.end(), val);
    }

private:
    std::size_t rows_;
    std::size_t cols_;
    std::vector<double> data_;
};

} // namespace schwartz

#endif // MATRIX2D_HPP
