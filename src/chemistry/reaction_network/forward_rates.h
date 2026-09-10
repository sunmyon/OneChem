#ifndef FORWARD_RATES_H
#define FORWARD_RATES_H

typedef double (*ForwardRateFunction)(double temperature);

extern const ForwardRateFunction forward_rate_table[];

#endif
