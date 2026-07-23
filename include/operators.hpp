#pragma once
#include "core.hpp"

namespace kr {

template <typename T>
Matrix<T> Matrix<T>::operator*(const Matrix& other){
   if(cols_ != other.rows_){
        throw std::runtime_error("Size does not match"); 
    }
    
    Matrix<T> product(rows_, other.cols_);

    for(size_t i = 0; i < rows_; ++i){
        for(size_t j = 0; j < other.cols_; ++j){
            for(size_t k = 0; k < cols_; ++k){
                product.data[(i * other.cols_) + j] +=
                    data[(i * cols_) + k] * 
                    other.data[(k * other.cols_) + j];
            }
        }
    }

    return product;
}

template <typename T>
Matrix<T> Matrix<T>::operator+(const Matrix& other){
    if(rows_ != other.rows_ || cols_ != other.cols_){
        throw std::invalid_argument("Dimension mismatch");
    }
    
    Matrix<T> sum(rows_, cols_);

    for(size_t i = 0; i < data.size(); i++){
        sum.data[i] = data[i] + other.data[i];
    }

    return sum;
}

template <typename T>
Matrix<T>& Matrix<T>::operator+=(const Matrix& other) {
    if (rows_ != other.rows_ || cols_ != other.cols_) {
        throw std::invalid_argument("Dimension mismatch");
    }
    for (size_t i = 0; i < data.size(); ++i) {
        data[i] += other.data[i];
    }
    return *this;
}

template <typename T>
Matrix<T> Matrix<T>::operator-(const Matrix& other){
    if(rows_ != other.rows_ || cols_ != other.cols_){
        throw std::invalid_argument("Dimension mismatch");    
    }
    
    Matrix<T> diff(rows_, cols_);

    for(size_t i = 0; i < data.size(); i++){
        diff.data[i] = data[i] - other.data[i];
    }

    return diff;
}

template <typename T>
Matrix<T> Matrix<T>::operator-(const T scalar){
    Matrix<T> diff(rows_, cols_);

    for(size_t i = 0; i < data.size(); i++){
        diff.data[i] = data[i] - scalar;
    }

    return diff;
}


template <typename T>
Matrix<T>& Matrix<T>::operator-=(const Matrix& other) {
    if (rows_ != other.rows_ || cols_ != other.cols_) {
        throw std::invalid_argument("Dimension mismatch");
    }
    for (size_t i = 0; i < data.size(); ++i) {
        data[i] -= other.data[i];
    }
    return *this;
}

//% is hadamard product
template <typename T>
Matrix<T> Matrix<T>::operator%(const Matrix& other){
    if(rows_ != other.rows_ || cols_ != other.cols_){
        throw std::invalid_argument("Dimension mismatch");
    }
    
    Matrix<T> product(rows_, cols_);

    for(size_t i = 0; i < data.size(); i++){
        product.data[i] = data[i] * other.data[i];
    }

    return product;
}

template <typename T>
Matrix<T>& Matrix<T>::operator%=(const Matrix& other) {
    if (rows_ != other.rows_ || cols_ != other.cols_) {
        throw std::invalid_argument("Dimension mismatch");
    }
    for (size_t i = 0; i < data.size(); ++i) {
        data[i] *= other.data[i];
    }
    return *this;
}

}