#pragma once 

#include "core.hpp"
#include "operators.hpp"
#include "views.hpp"

namespace kr {

template <typename T>
LU<T> Matrix<T>::lu() const {
    if (rows_ != cols_) throw std::runtime_error("Undefined");

    LU<T> result;
    result.U = *this;
    result.P = identity<T>(rows_);
    result.L = identity<T>(rows_);

    const T epsilon = static_cast<T>(1e-9);

    for(size_t i = 0; i < cols_; ++i){
        //Find pivot row index
        size_t maxId = i;
        T max = std::fabs(result.U(i, i));

        for(size_t j = i + 1; j < rows_; ++j){
            T comp = std::fabs(result.U(j, i));
            if(comp > max){ 
                max = comp; 
                maxId = j; 
            };
        }

        //Swap rows/cols if pivot isn't on the current row
        if(i != maxId){
            result.U.swapRow(i, maxId);
            result.P.swapRow(i, maxId);

            for (size_t col = 0; col < i; ++col) {
                std::swap(result.L(i, col), result.L(maxId, col));
            }

            result.rowOps++;
        }

        //Singular check
        if (std::fabs(result.U(i, i)) <= epsilon) {
            result.U(i, i) = 0;
            return result;
        }

        //Elimination
        for(size_t k = i + 1; k < rows_; ++k){
            T factor = result.U(k, i) / result.U(i, i);
            result.L(k, i) = factor;
            
            for(size_t n = i; n < cols_; ++n){
                result.U(k, n) -= result.U(i, n) * factor;
            }
        }
    }

    return result;
}

template <typename T>
T Matrix<T>::determinant() const {
    if (rows_ != cols_) {
        throw std::runtime_error("Undefined");
    }

    auto diag = lu();

    //See if det is negative or positive based off row swaps
    T det = (diag.rowOps % 2 == 0) ? static_cast<T>(1) : static_cast<T>(-1);

    for(size_t i = 0; i < cols_; ++i){
        T pivot = diag.U(i, i);

        if (pivot == 0) {
            return static_cast<T>(0);
        }

        det *= pivot;
    }

    return det;
}

}