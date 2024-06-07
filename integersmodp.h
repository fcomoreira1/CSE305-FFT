#pragma once
#include <ostream>
template <int p> class IntegersModP {
    int val;
    IntegersModP(int value) : val(value % p) {
        // Ensure non-negative value
        if (val < 0)
            val += p;
    }
    IntegersModP<p> pow(int n) {
        IntegersModP<p> res(1);
        int x = val;
        while (n > 0){
            if (n % 2) {
                res *= x;
            }
            n = n >> 1;
            x = x * x;
        }
        return res;
    }
    IntegersModP<p> &operator=(const IntegersModP<p> &other) {
        if (this != &other) {
            val = other.val;
        }
        return *this;
    }

    IntegersModP<p> operator+(const IntegersModP<p> &other) const {
        return IntegersModP<p>((val + other.val) % p);
    }

    IntegersModP<p> operator-(const IntegersModP<p> &other) const {
        return IntegersModP<p>((val - other.val + p) % p);
    }

    IntegersModP<p> operator*(const IntegersModP<p> &other) const {
        return IntegersModP<p>((val * other.val) % p);
    }

    IntegersModP<p> &operator+=(const IntegersModP<p> &other) {
        val = (val + other.val) % p;
        return *this;
    }

    IntegersModP<p> &operator-=(const IntegersModP<p> &other) {
        val = (val - other.val + p) % p;
        return *this;
    }

    IntegersModP<p> &operator*=(const IntegersModP<p> &other) {
        val = (val * other.val) % p;
        return *this;
    }

    IntegersModP<p> &operator/=(const IntegersModP<p> &other) {
        *this = *this / other;
        return *this;
    }

    // Comparison operators
    bool operator==(const IntegersModP<p> &other) const {
        return val == other.val;
    }

    bool operator!=(const IntegersModP<p> &other) const {
        return val != other.val;
    }
    friend std::ostream &operator<<(std::ostream &os,
                                    const IntegersModP<p> &obj) {
        os << obj.val;
        return os;
    }
    static IntegersModP<p> primitive_root() {

    }
};
