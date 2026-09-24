C---------------------------------------------------------------------------------------------
C This code defines and calculates parameters needed commonly by the forward and los routines
C author -  Pratik P Joshi [JUN 28, 2023]
C contact ppjoshi2@illinois.edu for questions
C --------------------------------------------------------------------------------------------

      SUBROUTINE global_parameters(msis_p1x,msis_p1,msis_p2x,msis_p2,
     &                             ap_histx,ap_hist,
!      &                             rtbkgx,rtbkgy,rtbkgz,rtbkg,
     &                             thbkgx,thbkgy,thbkgz,thbkg,
     &                             zgridx,zgrid)
      IMPLICIT REAL*8 (A-H,O-Z)
C     MSIS parameters (STL not used here)
C     thermospheric atomic hydrogen parameters
      INTEGER msis_p1x
      INTEGER msis_p1(msis_p1x)
      INTEGER msis_p2x
      REAL*8 msis_p2(msis_p2x)
      INTEGER ap_histx
      REAL*4 ap_hist(ap_histx)
!       INTEGER rtbkgx,rtbkgy,rtbkgz
!       REAL*8 rtbkg(rtbkgx,rtbkgy,rtbkgz)
      INTEGER thbkgx,thbkgy,thbkgz
      REAL*8 thbkg(thbkgx,thbkgy,thbkgz)
      INTEGER zgridx
      REAL*8 zgrid(zgridx)

C       SPEED OF LIGHT
        PARAMETER (SPEEDC = 2.9979D10)
C       BOLTZMANN CONSTANT
        PARAMETER (BOLTZ  = 1.3806D-16)
C       INTEGRATED ABSORPTION CONSTANT
        PARAMETER (ABSIC  = 2.654D-2)
C       ATOMIC MASS
        PARAMETER (AMASS  = 1.6738D-24)
C       FIXED
        PARAMETER (IKNT  = 76)
C       Used to be JKNT
        PARAMETER (NPT  =  1)

        DIMENSION HUT(61), WUT(61), OxUT(61),O2UT(61),TUT(61),Z(61)
        DIMENSION I_STEP(IKNT),IFIBO(IKNT)

        ! REAL, DIMENSION(1:IKNT,1:NPT) :: nHC,nHB,nO2,Tn

        REAL*8 THERMO_MSIS(7,61),TSAT,DSAT,ABSCSX
        REAL*8 DEXO,TEXO, T_EXO, FLUX_C
        REAL*4 WAVELN, FNUMBER, BRATIO, ABS_N2, ABS_O2
        REAL*4 GLAT, GLONG, UT, STL, F107, F107A, AP(7)
        REAL*4 SATT, SATD, SCALN, ADJT, BASE, TOP
        CHARACTER*13 CHOOZ
        INTEGER IDAY,IYEAR, cnt

        REAL*8 Z_GRID(IKNT)

C        COMMON/PARA_MSIS/ GLAT,GLONG,IDAY,IYEAR,UT,STL,IAPH,AP,F107,
C     &                    F107A
        COMMON/H_PROFILE/ D_EXO,FLUX,D_MAX,ALT_JNT,ALT_MAX,T_EXO,
     &                    MOD_MSIS

        COMMON/PARA_RT/ LINE_LABEL,WAVELN,FNUMBER,BRATIO,ABS_N2,ABS_O2
        COMMON /EXOS/  GM,CONL,PLANETR,RBASE,RC,RUPR,RP,TEXO,DEXO,VELT,
     &               RADPF,CENTER,ABSCSX,FTSAT,FDSAT,IGEO
        COMMON/PARA_EXOS/ IGEOM,SATT,SATD,SCALN,ADJT
        COMMON/FINE_GRID/ THERMO_MSIS,BASE,TOP,ITHERM
        COMMON /NMBR/  PI,RTPI,PID2,OFFSET
        COMMON IFIBO,I_STEP,CHOOZ
        COMMON/SOLAR_FLUX/ FLUX_C

C       Define constants
C       G-CONSTANT * PLANET MASS
        DATA GM       / 3.9898D20/
C       MEAN PLANET RADIUS
        DATA PLANETR  / 6371.0D5 /
C       GENERIC OFFSET VALUE
        DATA OFFSET   / 1.0D-6   /

        LINE_LABEL = 1
        WAVELN = 1215.67
        FNUMBER = 0.4162
        BRATIO = 1
        ABS_N2 = 0.00E+00
        ABS_O2 = 8.70E-21
        ITHERM = 60
        CHOOZ = 'LYMAN ALPHA'

        PID2  = ASIN(1.0D0)
        PI    = 2.0D0 * PID2
        RTPI  = SQRT(PI)

        CONL  = GM * AMASS / BOLTZ

        ! nHC = rtbkg(1,:,:)
        ! nO2 = rtbkg(2,:,:)
        ! Tn = rtbkg(3,:,:)

C      Read MSIS parameters
       IDAY = msis_p1(1)
       IYEAR = msis_p1(2)
       UT = msis_p2(1)
       IAPH = msis_p1(3)
       AP = ap_hist
       GLAT = msis_p2(2)
       GLONG = msis_p2(3)
       F107 = msis_p2(4)
       F107A = msis_p2(5)

C       Read thermospheric atomic hydrogen parameters
        MOD_MSIS = msis_p1(4)
        D_EXO = msis_p2(6)
        zexo = msis_p2(11)
        FLUX = msis_p2(12)
        D_MAX = msis_p2(13)
        IGEO = msis_p1(5)
        SATT = msis_p2(14)
        SATD = msis_p2(15)
        ALT_MAX = msis_p2(16)

        IGEOM=IGEO

C       Read H density low resolution grid (LR)
        Z_GRID = zgrid

        TSAT  = DBLE(SATT)
        DSAT  = DBLE(SATD)

C       Read MSIS temperatures and densities
        THERMO_MSIS = thbkg(:,1,:)

        ADJT    = 0.0
        ALT_JNT = 112.5d0
        SCALN = SNGL(D_EXO / THERMO_MSIS(4,ITHERM+1))

C Define thermospheric grid to define background atmosphere
C  for Lyman alpha, O2 photoabsorption is the main absorber channel but
C  the "black level" is some distance beneath the region where thermospheric
C  temperatures rise rapidly.  the QPIVOT approach gives a reasonably fine
C  reference grid for capturing both the rapid fall-off in O2 densities
C  above 80 km and the temperature ramp-up at altitudes above 100 km while
C  providing a monotonically increasing step size at higher altitudes where
C there is little absorption, temperature is effectively constant, and
C  atomic hydrogen densities decrease slowly with altitude.
        IF (LINE_LABEL.EQ.1) THEN
C       base at 74 km
          BASE   =  74.0
          QPIVOT = 102.0
          QFCTR  = 0.118418
C          TOP    = BASE + 376.0     392.0  exobase at 466 km
          TOP    = zexo
C          PRINT*,'TOP',TOP
        END IF

C Get thermospheric H density and temporature on the thermosphere grid
C        CALL H_DENSITY
C        THERMO_MSIS(4,ITHERM+1) = D_EXO

        TEXO  = THERMO_MSIS(3,ITHERM+1)
        DEXO  = THERMO_MSIS(4,ITHERM+1)
        T_EXO = TEXO

C Define RT grid resolution using user-defined input ZGRID
        I_STEP=Z_GRID - Z_GRID(1)+1

        IF (IGEO.EQ.1) THEN
          IF (DSAT.LT.0.0D0) THEN
            TSAT = TEXO
            DSAT = DEXO
          END IF
          FTSAT = TSAT/TEXO
          FDSAT = DSAT/DEXO
        ELSE IF (IGEO.EQ.0) THEN
          FTSAT = TSAT * PLANETR
          FDSAT = 0.0D0
          DSAT  = ABS(DSAT)
        END IF
        RBASE = PLANETR + DBLE(BASE)*1.0D5
        RC    = PLANETR + DBLE(TOP) *1.0D5
        RUPR  = RBASE + DBLE(I_STEP(IKNT))*1.0D5
        RADPF = DBLE(2.91 * (1.0 + 0.002*(F107-65.0)))
        RP    = SQRT(GM/0.1774D0/RADPF)

C        PRINT*,'F107',F107

C       Get doppler profile at exobase temperature
C       VELT   . . . MOST PROBABLE SPEED
C       CENTER . . . LINE CENTER SCATTERING CROSS SECTION
C       ABSCSX . . . ABSORBER PHOTOABSORPTION CROSS SECTION

        VELT   = SQRT(2.0d0 * BOLTZ / AMASS) * SQRT(TEXO)
        CENTER = ABSIC * DBLE(FNUMBER*WAVELN) * 1.0d-8 / RTPI / VELT
        ABSCSX = DBLE(ABS_O2)

        END SUBROUTINE global_parameters
C END FILE global_parameters.F