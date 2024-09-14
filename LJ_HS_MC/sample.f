      SUBROUTINE SAMPLE(ibeg,switch,En, Vir, rhoav)
c
c      write quantities (pressure and energy) to file
c
c
c  Ener (input) : total energy
c  Vir  (input) : total virial
c  grL=0.5sigma, grH=3sigma for sampling g(r)
c
      IMPLICIT NONE
      INCLUDE 'parameter.inc'
      INCLUDE 'conf.inc'
      INCLUDE 'system.inc'
      INCLUDE 'potential.inc'
	INTEGER switch, nhgr, nsample, i, j, ig,ibeg,Np
      DOUBLE PRECISION En, enp, Vir, press, delg, delgi, CORU, CORP
      DOUBLE PRECISION grL,grh,grH2,dx, dy, dz, r2, rx, dV
      PARAMETER (nhgr=250, grL=0.0, grH=3.0)
	DOUBLE PRECISION r(nhgr),g(nhgr),enpx, pressx , rhob, rhoav,avNp
      SAVE nsample, delg, delgi, g, grH2, enp, press, np

      IF (Switch.EQ.0) THEN

c        ---Initialize
	   grH2=(grH+Rcav)**2
         delg = (grH-grL)/nhgr
         delgi = 1./delg

	 if(ibeg.lt.2) then
	   nsample = 0
	   Np = 0
	   enp = 0
	   press = 0

         DO i = 1, nhgr
         g(i) = 0
         END DO
	  else

          OPEN(67,FILE='intermid')
          read(67,*) nsample, np, enp, press
		DO i = 1, nhgr
          read(67,*) g(i)
          END DO
          CLOSE(67)

	  endif

      ELSEIF (Switch.EQ.1) THEN

	   nsample = nsample + 1
c        ---sample average energy and pressure
         np = np + NPART
         enp = enp + En
         press = press + Vir

c        ---sample radial distribution function
               DO j = 1, NPART
				dx = X(j) - HBOX
				dy = Y(j) - HBOX
				dz = Z(j) - HBOX

                  r2 = dx*dx + dy*dy + dz*dz

                  IF (r2.LE.grH2) THEN

                     rx = dSQRT(r2) - Rcav
                     ig = INT(rx*delgi) + 1
                     g(ig) = g(ig) + 1
                  END IF
               END DO

               IF (MOD(nsample,10).EQ.0) THEN
				avNp = float(Np/nsample)
				enpx = enp/float(nsample)
				pressx = 1. + press/(3.*temp*float(nsample))

				open(66,file='out1')
				write(66,*) nsample/10,avNp,enpx,pressx
				close(66)

                 OPEN(67,FILE='intermid')
                 WRITE(67,*) nsample, Np, enp, press
		       DO i = 1, nhgr
                 WRITE(67,*) g(i)
                 END DO
                 CLOSE(67)
			 ENDIF

      ELSE IF (Switch.EQ.2) THEN

c        ---radial distribution function
            DO i = 1, nhgr
               r(i) = Rcav + grL +(i-0.5)*delg
               dV = 4.*PI*((r(i)+delg/2)**3-(r(i)-delg/2)**3)/3.
               g(i) = g(i)/dV/float(nsample)
            END DO

	   avNp = float(Np/nsample)
	   rhob = g(nhgr)
         enp = enp/avNp/float(nsample)
         press = 1. + press/(3.*avNp*temp*float(nsample))
	   rhoav=avNp/(Box**3-4.*pi*Rcav**3/3.)


c        ---output results
		 OPEN (7,FILE='OUTPUT')
		 WRITE(7,*) 'number of samples', nsample
		 WRITE(7,*) 'Average number of particles', avNp
		 WRITE(7,*) 'average density', rhoav 
		 WRITE(7,103) box, npart,Temp,rhob
		 WRITE (7, 101) enp, press
		 WRITE (7, 104) CORU(RC, rhoav),CORP(RC, rhoav)/(rhoav*temp) 
		 DO i=1,nhgr
		 WRITE (7, 102)  r(i), g(i)
		 ENDDO
	ENDIF
101	format(1x,'Potential Energy, kT=',e10.4,2x,/1x,
     &         'Compressibility factor=',f6.4)
102	format(1x,f6.3,2x,f6.4)
103	format(1x,'box_L=',f6.2,/1x,'Npart=',I5,/1x,'TEMP=',f6.3,
     &        /1x,'rhob=',f6.4)
104	format(1x,'Potential Energy(tail), kT=',e10.4,2x,/1x,
     &         'Compressibility factor (tail)=',f8.4)
      RETURN
      END
