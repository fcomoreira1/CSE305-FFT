#include "integersmodp.h"
#include "utils.h"
#include <algorithm>
#include <iostream>
#include <optional>

IntegersModP IntegersModP::pow(IntegersModP n, int exp) {
    if (n == 0) {
        return IntegersModP(0);
    }
    exp = ((exp % (p - 1)) + p - 1) % (p - 1);
    if (exp < 0) {
        std::cerr << "Negative exponent mod p" << std::endl;
        exit(-1);
    }
    IntegersModP res(1);
    auto x = n;
    while (exp > 0) {
        if (exp % 2) {
            res = res * x;
        }
        exp = exp >> 1;
        x = x * x;
    }
    return res;
}

IntegersModP IntegersModP::primitive_root() {
    static std::optional<std::pair<IntegersModP, int>> omega;
    if (omega && omega->second == p) {
        return omega->first;
    }

    std::vector<int> prime_div = get_prime_divisors(p - 1);
    bool valid_root;
    for (int i = 2; i < p; i++) {
        IntegersModP i_modp(i);
        valid_root = true;
        for (auto q : prime_div) {
            if (pow(i_modp, (p - 1) / q).val == 1) {
                valid_root = false;
                break;
            }
        }
        if (valid_root) {
            omega = std::make_pair(i_modp, p);
            break;
        }
    }
    return omega->first;
}

static int gcd(int a, int b, int &x, int &y) {
    if (a == 0) {
        x = 0;
        y = 1;
        return b;
    }
    int x1, y1;
    int _gcd = gcd(b % a, a, x1, y1);
    x = y1 - (b / a) * x1;
    y = x1;
    return _gcd;
}

IntegersModP IntegersModP::inverse(IntegersModP n){
    int x, y;
    gcd(n.val, p, x, y);
    return x;
}
