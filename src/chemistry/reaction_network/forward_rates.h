#ifndef FORWARD_RATES_H
#define FORWARD_RATES_H


typedef double (*ForwardRateFunction)(double temperature);

extern const ForwardRateFunction *const forward_rate_table;

/* Evaluate before exponentiation so detailed balance cannot form 0 * inf. */
double chemistry_log_forward_rate_12(double temperature);


#endif
