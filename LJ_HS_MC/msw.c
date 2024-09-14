#include "mpi.h"
#include "unistd.h"
#include "stdio.h"
 
int main (int argc, char *argv[])
{
#define MAXCMD 1024
   int i, lencmd;
   char cmd[MAXCMD];
 
   MPI_Init(&argc, &argv);

   lencmd = 0;
   for (i = 1; i < argc; i++) {
     lencmd += (strlen(argv[i]) + 1);
   }
   if (lencmd > (MAXCMD - 1)) {
     printf("command too long\n");
     return(0);
   }

   strcpy(cmd, "");
   for (i = 1; i < argc; i++) {
     strcat(cmd, argv[i]);
     strcat(cmd, " ");
   }

/*
   printf("cmd=%s\n", cmd);
*/

   system(cmd);
 
   MPI_Finalize();
 
   return(0);
}
