#include "integersmodp.h"
#include "types.h"
#include <algorithm>
#include <iostream>
#include <optional>

template<int p>
IntegersModP<p> IntegersModP<p>::pow(IntegersModP<p> n, int exp) {
    if (n == 0) {
        return IntegersModP<p>(0);
    }
    exp = ((exp % (p - 1)) + p - 1) % (p - 1);
    if (exp < 0) {
        std::cerr << "Negative exponent mod p" << std::endl;
        exit(-1);
    }
    // return IntegersModP<p>(std::pow(n.val, exp));
    IntegersModP<p> res(1);
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
template IntegersModP<p> IntegersModP<p>::pow(IntegersModP<p> n, int exp);

static std::vector<int> get_prime_divisors(int n) {
    std::vector<int> divisors;
    for (int i = 2; i * i <= n; i++) {
        if (n % i == 0) {
            divisors.push_back(i);
            while (n % i == 0) {
                n /= i;
            }
        }
    }
    if (n != 1) {
        divisors.push_back(n);
    }
    return divisors;
}

template<int p>
IntegersModP<p> IntegersModP<p>::primitive_root() {
    static std::optional<IntegersModP<p>> omega;
    if (omega) {
        std::cout << "aaa";
        return *omega;
    }
    std::vector<int> prime_div = get_prime_divisors(p - 1);
    bool valid_root;
    for (int i = 2; i < p; i++) {
        IntegersModP<p> i_modp(i);
        valid_root = true;
        for (auto q : prime_div) {
            if (pow(i_modp, (p - 1) / q).val == 1) {
                valid_root = false;
                break;
            }
        }
        if (valid_root) {
            *omega = i_modp;
            break;
        }
    }
    return *omega;
}
template IntegersModP<p> IntegersModP<p>::primitive_root();
