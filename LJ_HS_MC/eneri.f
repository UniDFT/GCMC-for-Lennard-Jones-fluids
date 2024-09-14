**==eneri.spg  processed by SPAG 4.52O  at 18:49 on  6 Jun 1996
      SUBROUTINE ENERI(Xi, Yi, Zi, I, Jb, En, Vir)
      IMPLICIT NONE
      INCLUDE 'parameter.inc'
      INCLUDE 'conf.inc'
      INCLUDE 'system.inc'
 
      DOUBLE PRECISION Xi,Yi,Zi,En,dx,dy,dz, r2, Vir, virij, enij
      INTEGER I, j, Jb
c
      En = 0.D0
      Vir = 0.D0
      DO j = Jb, NPART
         IF (j.NE.I) THEN
            dx = Xi - X(j)
            dy = Yi - Y(j)
            dz = Zi - Z(j)
c           ---minimum immage
            IF (dx.GT.HBOX) THEN
               dx = dx - BOX
            ELSE
               IF (dx.LT.-HBOX) dx = dx + BOX
            END IF
            IF (dy.GT.HBOX) THEN
               dy = dy - BOX
            ELSE
               IF (dy.LT.-HBOX) dy = dy + BOX
            END IF
            IF (dz.GT.HBOX) THEN
               dz = dz - BOX
            ELSE
               IF (dz.LT.-HBOX) dz = dz + BOX
            END IF
            r2 = dx*dx + dy*dy + dz*dz
c           ---calculate energy and virial pair i,j
            CALL ENER(enij, virij, r2)
            En = En + enij
            Vir = Vir + virij
         END IF
      END DO
      Vir = Vir

c      ---external potential
            dx = Xi - HBOX
            dy = Yi - HBOX
            dz = Zi - HBOX

            r2 = dx*dx + dy*dy + dz*dz

                  IF (r2.LT.Rcav2) THEN
				En = 1.e6
                  END IF
      RETURN
      END
