**==store.spg  processed by SPAG 4.52O  at 18:49 on  6 Jun 1996
      SUBROUTINE STORE(Dr)
c
c     writes configuration to disk
c
      IMPLICIT NONE
      INCLUDE 'parameter.inc'
      INCLUDE 'conf.inc'
      INCLUDE 'system.inc'
      INTEGER i
      DOUBLE PRECISION Dr
 
	open(11,file='iniconf.dat')
      WRITE (11, *) BOX
      WRITE (11, *) NPART
      WRITE (11, *) Dr
      DO i = 1, NPART
         WRITE (11, *) X(i), Y(i), Z(i)
      END DO
      close(11)

      RETURN
      END
