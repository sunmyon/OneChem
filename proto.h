
int compute_abundances_d(int mode, int ithis, double a_start, double a_end,int *flag_implicit, int flag_print);
int H2line_tabulate(void);
int InitChem(void);

void make_rate_tables_LW_radiation(void);
void get_Hminus_PD_rate(double Teff, double *Hminus_pd_rate, double *G0);

int compute_gamma(int ithis);
double compute_temp(int ithis);

double timestep_flag(double dt_pre, int flag_chem);
double timestep(double dt_pre);
int density_update(int, int, double, double atime, double a_end, double *);
int density_update_phys(int, double dt_phys, double * dt_density);
int density_isobaric(int, double, double, double *);
double JeansMass(double T,double nH,double mu);
int compute_Jeans_length(int ithis, double T_gas);
double compute_freefall_time(double rho_gas);
double compute_collapse_time(double SHR, double rho_gas);
double update_phys(int ithis, double dt, double SHR0);

void endrun(int ierr);
