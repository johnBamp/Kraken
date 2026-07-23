#include <iostream>
#include <Kraken>

int main() {
    kr::Matrix<float> A(3, 3);
    kr::Matrix<float> b(3, 1);

    A = {10, -7, 0, 
         -3, 2, 6,
          5, -1, 5};

    b = {1, 2, 3};

    auto lu = A.lu();
    
    return 0;
}