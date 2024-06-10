#pragma once
#include <iostream>

class IntegersModP {
    long long val;

  public:
    static int p;
    IntegersModP() : val(0) {}
    IntegersModP(int value) : val(value % p) {
        if (val < 0) {
            val += p;
        }
    }
    IntegersModP &operator=(const IntegersModP &other) {
        if (this != &other) {
            val = other.val;
        }
        return *this;
    }
    int get_val() const { return val; }

    IntegersModP operator+(const IntegersModP &other) const {
        return IntegersModP((val + other.val) % p);
    }

    IntegersModP operator-(const IntegersModP &other) const {
        return IntegersModP((val - other.val) % p);
    }

    IntegersModP operator*(const IntegersModP &other) const {
        return IntegersModP((val * other.val) % p);
    }

    IntegersModP operator/(const IntegersModP &other) const {
        return IntegersModP((val * IntegersModP::inverse(other).val) % p);
    }

    bool operator==(const IntegersModP &other) const {
        return val == other.val;
    }

    bool operator!=(const IntegersModP &other) const {
        return val != other.val;
    }
    friend std::ostream &operator<<(std::ostream &os, const IntegersModP &obj) {
        os << obj.val;
        return os;
    }
    static IntegersModP primitive_root();
    static IntegersModP pow(IntegersModP n, int exp);
    static IntegersModP inverse(IntegersModP n);
};
void test_primitive_root();
void test_pow();
