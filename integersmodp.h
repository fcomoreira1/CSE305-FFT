#pragma once
#include <ostream>

const int p = 17;
template <int p> class IntegersModP {
  public:
    int val;
    IntegersModP(int value) : val(value % p) {
        // Ensure non-negative value
        if (val < 0)
            val += p;
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

    // IntegersModP<p> &operator+=(const IntegersModP<p> &other) {
    //     val = (val + other.val) % p;
    //     return *this;
    // }

    // IntegersModP<p> &operator-=(const IntegersModP<p> &other) {
    //     val = (val - other.val + p) % p;
    //     return *this;
    // }

    // IntegersModP<p> &operator*=(const IntegersModP<p> &other) {
    //     val = (val * other.val) % p;
    //     return *this;
    // }

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
    static IntegersModP<p> primitive_root();
    static IntegersModP<p> pow(IntegersModP<p> n, int exp);
    static IntegersModP<p> inverse(IntegersModP<p> n) { return pow(n, p - 2); }
};
void test_primitive_root();
