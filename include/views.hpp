#pragma once

#include "core.hpp"
#include "operators.hpp"

namespace kr {

template <typename T>
void Matrix<T>::swapRow(size_t r1, size_t r2){
    for(size_t i = 0; i < cols_; ++i){
        std::swap(data[(r1 * cols_) + i], data[(r2 * cols_) + i]);
    }
}

template <typename T>
void Matrix<T>::swapCol(size_t c1, size_t c2){
    for(size_t i = 0; i < rows_; ++i){
        std::swap(data[(i * cols_) + c1], data[(i * cols_) + c2]);
    }
}

template <typename T>
Matrix<T>& Matrix<T>::transpose(){
    std::vector<T> result(rows_ * cols_);
    
    for(size_t i = 0; i < rows_; ++i){
        for(size_t j = 0; j < cols_; ++j){
            result[(j * rows_) + i] = data[(i * cols_) + j];
        }
    }
    
    data = std::move(result); 
    std::swap(rows_, cols_);
    
    return *this;
}

}