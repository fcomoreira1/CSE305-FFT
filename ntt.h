#pragma once
#include "types.h"
#include "integersmodp.h"

void ntt_baseline(const IntegersModP<p> *x, IntegersModP<p> *y, int n);
void intt_baseline(const IntegersModP<p> *y, IntegersModP<p> *x, int n);
void ntt_radix2_seq_(const IntegersModP<p> *x, IntegersModP<p> *y, int n,
                     int d = 1);
void intt_radix2_seq_(const IntegersModP<p> *y, IntegersModP<p> *x, int n,
                      int d = 1);
void ntt_radix2_seq(const IntegersModP<p> *x, IntegersModP<p> *y, int n);
void intt_radix2_seq(const IntegersModP<p> *y, IntegersModP<p> *x, int n);
