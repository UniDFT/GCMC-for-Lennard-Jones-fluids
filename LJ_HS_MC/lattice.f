**==lattice.spg  processed by SPAG 4.52O  at 18:49 on  6 Jun 1996
 
      SUBROUTINE LATTICE(rho)
c     ---place `npart' particles on a lattice with density 'rho'
      IMPLICIT NONE
      INCLUDE 'parameter.inc'
      INCLUDE 'conf.inc'
      INCLUDE 'system.inc'
	INCLUDE 'potential.inc'

      INTEGER i, j, k, itel, n
      DOUBLE PRECISION dx, dy, dz, del, r2, rho
 
      del = (BOX**3)**(1.D0/3.D0)
      n = INT(NPART**(1.D0/3.D0)) + 1
      IF (n.EQ.0) n = 1
      del = del/DBLE(n)
      itel = 0
      dx = -del
      DO i = 1, n
         dx = dx + del
         IF (dx.GT.BOX) dx = dx - BOX
         dy = -del
         DO j = 1, n
            dy = dy + del
            IF (dy.GT.BOX) dy = dy - BOX
            dz = -del
            DO k = 1, n
               dz = dz + del
               IF (dz.GT.BOX) dz = dz - BOX
               IF (itel.LT.NPART) THEN
                  itel = itel + 1
                  X(itel) = dx
                  Y(itel) = dy
                  Z(itel) = dz
               END IF
            END DO
         END DO
      END DO

	i=0
           DO j = 1, itel
            dx = X(j) - HBOX
            dy = Y(j) - HBOX
            dz = Z(j) - HBOX

            r2 = dx*dx + dy*dy + dz*dz

                  IF (r2.GT.Rcav2) THEN
				i=i+1
                  X(i) = X(j)
                  Y(i) = Y(j)
                  Z(i) = Z(j)
                  END IF

            END DO

	Npart=i

	rho=Npart/(BOX**3-4.*pi*Rcav**3/3.)

      WRITE (*, *)  
      WRITE (*, 99001) Npart

      RETURN
99001 FORMAT (' Initialisation on lattice: ', /, i10, 
     &        ' particles placed on a lattice')
      END
