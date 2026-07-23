#pragma once

#include <initializer_list>
#include <ostream>
#include <stdexcept>
#include <vector>

namespace kr {

template <typename T>
class Matrix {
private:
    std::vector<T> data; 
    
    size_t rows_, cols_;

public:

    Matrix(size_t rows, size_t cols) : rows_(rows), cols_(cols) {
        data.resize(rows * cols);
    }

    Matrix& operator=(std::vector<T> other){
        if(other.size() != data.size()) {throw std::runtime_error("Size does not match"); }

        data = other;

        return *this;
    }

    const T& operator()(size_t i, size_t j) const{
        return data[(i * cols_) + j];
    }

    T& operator()(size_t i, size_t j) {
        return data[(i * cols_) + j];
    }

    Matrix operator*(const Matrix& other){
       if(rows_ != other.cols_){
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


    Matrix operator+(const Matrix& other){
        if(rows_ != other.rows_ || cols_ != other.cols_){
            throw std::runtime_error("Size does not match"); 
        }
        
        Matrix<T> sum(rows_, cols_);

        for(size_t i = 0; i < data.size(); i++){
            sum.data[i] = data[i] + other.data[i];
        }

        return sum;
    }

    Matrix operator+=(const Matrix& other){
        *this = *this + other;
        return *this;
    }

    Matrix operator-(const Matrix& other){
        if(rows_ != other.rows_ || cols_ != other.cols_){
            throw std::runtime_error("Size does not match"); 
        }
        
        Matrix<T> diff(rows_, cols_);

        for(size_t i = 0; i < data.size(); i++){
            diff.data[i] = data[i] - other.data[i];
        }

        return diff;
    }

    Matrix operator-=(const Matrix& other){
        *this = *this - other;
        return *this;
    }

    Matrix operator%(const Matrix& other){
        if(rows_ != other.rows_ || cols_ != other.cols_){
            throw std::runtime_error("Size does not match"); 
        }
        
        Matrix<T> product(rows_, cols_);

        for(size_t i = 0; i < data.size(); i++){
            product.data[i] = data[i] * other.data[i];
        }

        return product;
    }

    Matrix operator%=(const Matrix& other){
        *this = *this % other;
        return *this;
    }


    friend std::ostream& operator<<(std::ostream& os, const Matrix& a) {
        for(size_t i = 0; i < a.data.size(); i++){
            if(i != 0 && i % a.cols_ == 0){
                std::cout << std::endl;
            }
            std::cout << a.data[i] << ' ';
        }

        return os;
    }

    

};

} // namespace kr