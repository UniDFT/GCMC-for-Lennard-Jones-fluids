      SUBROUTINE 
     &  READDAT(ibeg, Equil, Prod, Nsamp, Ndispl, Dr, Nexch, Iseed)
c
C     reads input data and model parameters
c
      IMPLICIT NONE
      INCLUDE 'parameter.inc'
      INCLUDE 'system.inc'
      INCLUDE 'potential.inc'
      INCLUDE 'conf.inc'
      INCLUDE 'grand.inc'
      INTEGER ibeg, Equil, Prod, i, Ndispl, Nsamp, Nexch, Iseed
      DOUBLE PRECISION eps, sig, CORU, CORP, vir, rho, Dr, pid
 
 
c     ---read simulation data
c
c input parameters
c      ibeg    = 0 start simulation from lattice
c              else continue sumlation from disk
c      equil   = number of equilibration cycles
c      prod    = number of produktion cycles
c      nsample = sample frequency; sample after nsample cycles
c      iseed   = seed random number
c      dr      = maximum displacement (only for ibeg=0)
c      ndispl  = number of attempts to displace a particle per cycle
c      nexch   = number of attempts to exchange particles per cycle
c      npart   = total number of particles
c      temp    = temperture
c      rho     = initial density  (only for ibeg=0)
c      pid     = ideal gas pressure reservoir (transfered to chemical potential)
c
c     ---read model parameters
	open (5, file='input.dat')
      READ (5, *)
      READ (5, *) ibeg, Equil, Prod, Nsamp, Iseed
      READ (5, *)
      READ (5, *) Dr, Rcav
      READ (5, *)
      READ (5, *) Ndispl, Nexch
      READ (5, *)
      READ (5, *) NPART, TEMP, rho, pid
      READ (5, *)
      READ (5, *) TAILCO, SHIFT
      READ (5, *)
      READ (5, *) eps, sig, MASS, RC
	close(5)

c     ---calculate parameters:
      PI = 4.D0*ATAN(1.D0)
      BETA = 1.D0/TEMP
	Rcav2=Rcav**2

c     ---grand canonical coupled to an ideal gas bath with pressure: pid
c        chemical potential bath: mu^b = mu^0 + ln(beta*pid)/beta
c                                      = ln (beta*pid*Lamda^3)/beta
c        zz is defined as         zz   = exp(beta*mu^b)/Lamda^3
c                                      = beta*pid
c        excess chemical pot.     muex = mu^b -mu^id
c                                      = mu^0 + ln(beta*pid)/beta - mu^0 - ln(rho)
c                                      = ln(zz)/beta - ln <rho>
      ZZ = BETA*pid

      IF (NPART.GT.NPMax) THEN
         WRITE (6, *) ' ERROR: number of particles too large'
         STOP
      END IF

c     ---read/generate configuration
      IF (ibeg.EQ.0) THEN
c        ---generate configuration form lattice
         BOX = (NPART/rho+4.*pi*Rcav**3/3.)**(1.D00/3.D00)
	   HBOX = 0.5D0*BOX
         CALL LATTICE(rho)
      ELSE
         WRITE (6, *) 
         WRITE (6, *) ' read conf from disk '

	   OPEN (11,file='iniconf.dat')
         READ (11, *) BOX
	   HBOX = 0.5D0*BOX
         READ (11, *) NPART
         READ (11, *) Dr
         DO i = 1, NPART
            READ (11, *) X(i), Y(i), Z(i)
         END DO
	   CLOSE(11)
         rho = NPART/(BOX**3-4.*pi*Rcav**3/3.)
      END IF

c     ---write input data
      WRITE (6, 99001) Equil, Prod, Nsamp
      WRITE (6, 99002) Ndispl, Dr, Nexch
      WRITE (6, 99003) NPART, TEMP, pid, rho, BOX
      WRITE (6, 99004) eps, sig, MASS
c     ---calculate cut-off radius potential
      RC = MIN(RC, HBOX)
      RC2 = RC*RC
      EPS4 = 4.D0*eps
      EPS48 = 48.D0*eps
      SIG2 = sig*sig
      IF (SHIFT) THEN
c     ---calculate energy of the shift (dirty trick)
         ECUT = 0.D0
         CALL ENER(ECUT, vir, RC2)
         WRITE (6, 99005) RC, ECUT
      END IF
      IF (TAILCO) THEN
         WRITE (6, 99006) RC, CORU(RC, rho), CORP(RC, rho)
      END IF

      RETURN
99001 FORMAT ('  Number of equilibration cycles             :', i10, /, 
     &        '  Number of production cycles                :', i10, /, 
     &        '  Sample frequency                           :', i10, /)
99002 FORMAT ('  Number of att. to displ. a part. per cycle :', i10, /, 
     &        '  Maximum displacement                       :', f10.3, 
     &        /, '  Number of att. to exchange part. per cycle :', i10, 
     &        //)
99003 FORMAT ('  Number of particles                        :', i10, /, 
     &        '  Temperature                                :', f10.3, 
     &        /, '  Ideal gas pressure reservoir               :', 
     &        f10.3, /, '  Density                                    :'
     &        , f10.3, /, 
     &        '  Box length                                 :', f10.3, 
     &        /)
99004 FORMAT ('  Model parameters: ', /, '     epsilon: ', f5.3, /, 
     &        '     sigma  : ', f5.3, /, '     mass   : ', f5.3)
99005 FORMAT (' Simulations with TRUNCATED AND SHIFTED potential: ', /, 
     &        ' Potential truncated at :', f10.3, /, 
     &        ' Energy shif            :', f10.6, //)
99006 FORMAT (' Simulations with tail correction: ', /, 
     &        ' Potential truncated at  Rc =:', f10.3, /, 
     &        ' Tail corrections:   energy = ', f10.3, ' pressure ', 
     &        f10.3, /, /)
99007 FORMAT (' Requested density: ', f5.2, 
     &        ' different from density on disk: ', f5.2, /, 
     &        ' Rescaling of coordinates!')
      END
