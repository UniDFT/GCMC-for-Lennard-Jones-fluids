**==toterg.spg  processed by SPAG 4.52O  at 18:49 on  6 Jun 1996
      SUBROUTINE TOTERG(Ener, Vir)
c
c     ---calculates total energy and virial
c
      IMPLICIT NONE
      INCLUDE 'parameter.inc'
      INCLUDE 'conf.inc'
      INCLUDE 'potential.inc'
      INCLUDE 'system.inc'
 
      DOUBLE PRECISION xi, yi, zi, Ener, eni, CORU, viri, Vir, rho,
     &                 dx, dy, dz, r2
      INTEGER i, jb
 
      Ener = 0
      Vir = 0
      DO i = 1, NPART - 1
         xi = X(i)
         yi = Y(i)
         zi = Z(i)
         jb = i + 1
c       ---calculate energy particle i with particel j=jb,naprt
         CALL ENERI(xi, yi, zi, i, jb, eni, viri)
         Ener = Ener + eni
         Vir = Vir + viri
      END DO
c     ---add tail corrections
      IF (TAILCO) THEN
         rho = NPART/(BOX**3)
         Ener = Ener + NPART*CORU(RC, rho)
      END IF

c      ---external potential
               DO i = 1, NPART
	            dx = X(i) - HBOX
	            dy = Y(i) - HBOX
	            dz = Z(i) - HBOX

	            r2 = dx*dx + dy*dy + dz*dz

                  IF (r2.LT.Rcav2) THEN
				Ener = 1.e6
                  END IF
               END DO

      RETURN
      END
 
