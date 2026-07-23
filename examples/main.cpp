#include <iostream>
#include <Kraken>

int main() {

    kr::Matrix<float> a(3, 3);

    a = {1, 2, 3, 
         4, 5, 6,
         7, 8, 9};

    std::cout << a << std::endl << std::endl;

    auto lu = a.lu();

    std::cout << std::endl << std::endl;

    std::cout << lu.U << std::endl << std::endl;
    std::cout << lu.L << std::endl << std::endl;

    std::cout << lu.P << std::endl << std::endl;
    
    std::cout << lu.P.Transpose() << std::endl << std::endl;
    
    std::cout << (lu.P.Transpose() * lu.L * lu.U) << std::endl << std::endl;

    return 0;
}