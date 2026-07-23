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
    Matrix operator-(const T scalar);
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

} // namespace kr
