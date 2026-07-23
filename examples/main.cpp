#include <iostream>
#include <Kraken>

int main() {

    kr::Matrix<float> a(3, 2);
    kr::Matrix<float> b(2, 3);

    a = {1, 2, 
         3, 4,
         7, 8};

    
    b = {4, 3, 5,
         2, 1, 4};

    std::cout << a << std::endl;
    std::cout << b << std::endl;

    auto c = a * b;

    std::cout << c << std::endl; 

    return 0;
}