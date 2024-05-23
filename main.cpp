#include "algorithms.h"
#include <iostream>

int main() {
    int n, k;
    std::cin >> k;
    n = 1 << k;
    std::vector<Complex> x(n);
    for (int i = 0; i < n; i++) {
        std::cin >> x[i];
    }
    auto x_coef = forward_dct(x);
    auto z = inverse_dct(x_coef);
    for (auto t : z) {
        std::cout << t << " ";
    }
    std::cout << std::endl;
    return 0;
}
