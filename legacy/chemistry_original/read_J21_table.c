#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#include "../../src/standalone/allvars.h"
#include "../tags.h"
#include "chemistry.h"

#if defined(CHEMISTRY_SC) && defined(TABLE_J21_EVOLUTION)

static int compare_time(const void *a, const void *b);
static int n_table;

struct Time_J21{
  double time;
  double J21;
} *table_J21;

void init_J21_table(void){
  if(ThisTask == 0){
    FILE *fp;
    if(!(fp = fopen(All.Table_J21_evolution, "r"))){
      printf("Can't open the file, %s.\n", All.Table_J21_evolution);
      exit(1);
    }

    char buf[100];
    int count = 0;
    while(fgets(buf, sizeof(buf), fp) != NULL){
      int snapid;
      double redshift, J21, pos[3];
      int n = sscanf(buf, "%d %le %le %le %le %le", &snapid, &redshift, &J21, pos, pos+1, pos+2);
      if(n == 6)
	count++;
    }

    n_table = count;
    table_J21 = (struct Time_J21 *) mymalloc("table_J21", n_table * sizeof(struct Time_J21));

    rewind(fp);
    count = 0;
    while(fgets(buf, sizeof(buf), fp) != NULL){
      int snapid;
      double redshift, J21, pos[3];
      int n = sscanf(buf, "%d %le %le %le %le %le", &snapid, &redshift, &J21, pos, pos+1, pos+2);
      if(n == 6){
	table_J21[count].time = 1./(redshift + 1.);
	table_J21[count].J21  = J21;
	
	count++;
      }
    }
    
    fclose(fp);

    qsort(table_J21, n_table, sizeof(struct Time_J21), compare_time);
  }

  if(NTask > 1){
    MPI_Bcast(&n_table, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if(ThisTask != 0)
      table_J21 = (struct Time_J21 *) mymalloc("table_J21", n_table * sizeof(struct Time_J21));
    
    MPI_Bcast(table_J21, n_table*sizeof(struct Time_J21), MPI_BYTE, 0, MPI_COMM_WORLD);
  }
  
}


void set_J21(void){
  double time = All.Time;
  
  int i_t=0;
  double t_new = table_J21[i_t].time;
  while(t_new < time){
    i_t++;
    t_new = table_J21[i_t].time;

    if(i_t == n_table-1) break;
  }
  
  double J21=0.;
  if(i_t > 0){
    int i_t_left   = i_t - 1;

    if(i_t_left < 0) i_t_left = 0;
    if(i_t == n_table - 1) i_t_left = i_t - 1;

    double c_t = 0.;
    if(i_t > i_t_left)
      c_t = (table_J21[i_t].time - time) / (table_J21[i_t].time - table_J21[i_t_left].time);

    if(c_t > 1.) c_t = 1.;
    if(c_t < 0.) c_t = 0.;
    
    J21 = c_t * table_J21[i_t_left].J21 +  (1 - c_t) * table_J21[i_t].J21;  

    mpi_printf("it=%d ntable=%d ct=%g J21=%g\n", i_t, n_table, c_t, J21);
  }

  mpi_printf("J21 is %g\n", J21);
  
  All.J21 = J21;
}

static int compare_time(const void *a, const void *b)
{
  if(((struct Time_J21 *) a)->time < (((struct Time_J21 *) b)->time))
    return -1;

  if(((struct Time_J21 *) a)->time > (((struct Time_J21 *) b)->time))
    return 1;

  return 0;
}
#endif
