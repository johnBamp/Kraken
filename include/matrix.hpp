#pragma once

#include <initializer_list>
#include <ostream>
#include <stdexcept>
#include <vector>
#include <cmath>
#include <iostream>

namespace kr {

template <typename T>
class Matrix;

template <typename T>
struct LU {
    Matrix<T> P;
    Matrix<T> L; 
    Matrix<T> U;

    size_t rowOps = 0;
};

template <typename T>
Matrix<T> identity(size_t size){
    Matrix<T> id(size, size);

    for(size_t i = 0; i < size; ++i){
        id(i, i) = 1;
    }

    return id;
}

template <typename T>
class Matrix {
private:
    std::vector<T> data; 
    size_t rows_, cols_;

public:
    Matrix(size_t rows, size_t cols);
    Matrix();

    Matrix& operator=(std::vector<T> other);

    const T& operator()(size_t i, size_t j) const;
    T& operator()(size_t i, size_t j);

    Matrix operator*(const Matrix& other);
    Matrix operator+(const Matrix& other);
    Matrix& operator+=(const Matrix& other);
    Matrix operator-(const Matrix& other);
    Matrix& operator-=(const Matrix& other);
    Matrix operator%(const Matrix& other);
    Matrix& operator%=(const Matrix& other);

    friend std::ostream& operator<<(std::ostream& os, const Matrix& a) {
        for(size_t i = 0; i < a.data.size(); ++i){
            if(i != 0 && i % a.cols_ == 0){
                os << std::endl;
            }
            os << a.data[i] << ' ';
        }
        return os;
    }

    void swapRow(size_t r1, size_t r2);
    void swapCol(size_t c1, size_t c2);
    Matrix<T>& transpose();
    T determinant() const;
    LU<T> lu() const;
};

template <typename T>
Matrix<T>::Matrix(size_t rows, size_t cols) : rows_(rows), cols_(cols) {
    data.resize(rows * cols);
}

template <typename T>
Matrix<T>::Matrix() : rows_(0), cols_(0) {}

template <typename T>
Matrix<T>& Matrix<T>::operator=(std::vector<T> other){
    if(other.size() != data.size()) { throw std::runtime_error("Size does not match"); }
    data = other;
    return *this;
}

template <typename T>
const T& Matrix<T>::operator()(size_t i, size_t j) const{
    return data[(i * cols_) + j];
}

template <typename T>
T& Matrix<T>::operator()(size_t i, size_t j) {
    return data[(i * cols_) + j];
}

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

} // namespace kr
