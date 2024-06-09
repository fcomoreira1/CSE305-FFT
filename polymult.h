#include "integersmodp.h"
#include <functional>
#include <vector>
std::vector<int>
polmult_ntt(std::vector<int> P1, std::vector<int> P2,
            std::function<void(const IntegersModP *, IntegersModP *, int)> ntt,
            std::function<void(const IntegersModP *, IntegersModP *, int)> intt);
std::vector<int> polmult_baseline(std::vector<int> P1, std::vector<int> P2);
