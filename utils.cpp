#include "utils.h"

int intlog2(int n) {
    int log = 0;
    while (n > 1) {
        log++;
        n /= 2;
    }
    return log;
}
int pow2greater(int n) {
    int log = intlog2(n);
    return 1 << log == n ? n : 1 << (log + 1);
}
