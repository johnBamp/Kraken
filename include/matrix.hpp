#pragma once

#include <initializer_list>
#include <ostream>
#include <stdexcept>
#include <vector>

namespace kr {

template <typename T>
class Matrix;

template <typename T>
struct LU {
    Matrix<T> P;
    Matrix<T> L; 
    Matrix<T> U;

};

template <typename T>
Matrix<T> identity(size_t rows_, size_t cols_){
    Matrix<T> id(rows_, cols_);

    size_t row_;

    for(size_t i = 0; i < cols_; ++i){
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

    Matrix(size_t rows, size_t cols) : rows_(rows), cols_(cols) {
        data.resize(rows * cols);
    }

    Matrix(){}

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
        for(int i = 0; i < a.data.size(); ++i){
            if(i != 0 && i % a.cols_ == 0){
                std::cout << std::endl;
            }
            std::cout << a.data[i] << ' ';
        }

        return os;
    }

    void swapRow(size_t r1, size_t r2){
        for(int i = 0; i < cols_; ++i){
            std::swap(data[(r1 * cols_) + i], data[(r2 * cols_) + i]);
        }
    }

    Matrix<T>& Transpose(){
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
    
    LU<T> lu(){
        Matrix<T> U(rows_, cols_);
        Matrix<T> L(rows_, cols_);
        Matrix<T> P(rows_, cols_);
        
        P = identity<T>(rows_, cols_);
        L = P;

        U = *this;

        std::cout << L << std::endl;

        float epsilon = 1e-5f;

        for(size_t i = 0; i < cols_; ++i){
            T max = std::fabs(U.data[(i * cols_) + i]);

            size_t maxId = i;

            for(size_t j = i; j < rows_; ++j){
                T comp = std::fabs(U.data[(j * cols_) + i]);
                if(comp > max){ max = comp; maxId = j; };
            }

            U.swapRow(i, maxId);
            
            for(size_t row = 0; row < rows_; ++row) {
                std::swap(P.data[(row * cols_) + i], P.data[(row * cols_) + maxId]);
            }

            for(size_t col = 0; col < i; ++col) {
                std::swap(L.data[(i * cols_) + col], L.data[(maxId * cols_) + col]);
            }

            for(size_t k = i + 1; k < rows_; ++k){
                T div = U.data[(k * cols_) + i] / U.data[(i * cols_) + i];
                L.data[(k * cols_) + i] = div;

                for(size_t n = i; n < cols_; ++n){
                    U.data[(k * cols_) + n] -=  U.data[(i * cols_) + n] * div;
                }
            }

        }

        LU<T> result;

        result.U = U;
        result.L = L;
        result.P = P;

        return result;
        
    }

};

} // namespace kr