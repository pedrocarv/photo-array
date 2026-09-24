C---------------------------------------------------------------------------------------------
C This code defines and calculates parameters needed commonly by the forward and los routines
C author -  Pratik P Joshi [JUN 28, 2023]
C contact ppjoshi2@illinois.edu for questions
C --------------------------------------------------------------------------------------------

        SUBROUTINE global_parameters
        IMPLICIT REAL*8 (A-H,O-Z)
        PARAMETER (SPEEDC = 2.9979D10)         ! SPEED OF LIGHT
        PARAMETER (BOLTZ  = 1.3806D-16)        ! BOLTZMANN CONSTANT
        PARAMETER (ABSIC  = 2.654D-2)          ! INTEGRATED ABSORPTION CONSTANT
        PARAMETER (AMASS  = 1.6738D-24)        ! ATOMIC MASS
        PARAMETER (IKNT  = 32)                 ! FIXED
        PARAMETER (JKNT  =  32)                ! MUST BE EVEN
        DIMENSION HUT(61), WUT(61), OxUT(61),O2UT(61),TUT(61),Z(61)
        DIMENSION I_STEP(IKNT),IFIBO(IKNT)
        DIMENSION rtbkg(IKNT,JKNT,3),thermobkg(7,JKNT,61)
C        integer, parameter :: M = 40, N =19
        REAL, DIMENSION(1:IKNT,1:JKNT) :: nHC,nHB,nO2,Tn

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
C        COMMON/RT_GRID/ nHC,nHB,nO2,Tn
        COMMON/SOLAR_FLUX/ FLUX_C

C Define constants
        DATA GM       / 3.9898D20/         ! G-CONSTANT * PLANET MASS
        DATA PLANETR  / 6371.0D5 /         ! MEAN PLANET RADIUS
        DATA OFFSET   / 1.0D-6   /         ! GENERIC OFFSET VALUE

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

        open(unit=10, file = 'rt_bkg.txt', status = 'old')
        do k=1,3
            cnt=0
            do j =1+((k-1)*JKNT),k*JKNT
                cnt=cnt+1
                read(10,*) (rtbkg(i,cnt,k), i = 1, IKNT)
            end do
        end do
        close(10)
        nHC = rtbkg(:,:,1)
        nO2 = rtbkg(:,:,2)
        Tn = rtbkg(:,:,3)

        open(unit=10, file = 'thermo_bkg.txt', status = 'old')
        do k=1,7
            cnt=0
            do j =1+((k-1)*JKNT),k*JKNT
                cnt=cnt+1
                read(10,*) (thermobkg(k,cnt,i), i = 1, ITHERM+1)
            end do
        end do
        close(10)
        THERMO_MSIS = thermobkg(:,1,:)

C Read thermospheric atomic hydrogen parameters
        OPEN(98, FILE='INPUT_HDEN.txt', STATUS='OLD')
        READ(98,*) MOD_MSIS, D_EXO, FLUX, D_MAX, FLUX_C
!C981   FORMAT (I1, E7.1E2, E8.2E2, E8.2E2)
        READ(98,*) IGEO, SATT, SATD, ALT_MAX, F107

        IGEOM=IGEO

        OPEN(9, FILE='alt_file_LR.DAT', STATUS='OLD')
        DO I = 1,IKNT
          READ(9,*) Z_GRID(I)
        END DO
        CLOSE(9)
C        PRINT *, 'ZGRID=',Z_GRID
        TSAT  = DBLE(SATT)
        DSAT  = DBLE(SATD)

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
        BASE   =  74.0         ! base at 74 km
        QPIVOT = 102.0
        QFCTR  = 0.118418
C        TOP    = BASE + 376.0!392.0  ! exobase at 466 km
        TOP    = THERMO_MSIS(1,ITHERM+1)
        PRINT*,"TOP", TOP
        END IF

        TEXO  = THERMO_MSIS(3,ITHERM+1)
        DEXO  = THERMO_MSIS(4,ITHERM+1)
C        PRINT*,'EXO',TEXO,DEXO
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

C Get doppler profile at exobase temperature
C VELT   . . . MOST PROBABLE SPEED
C CENTER . . . LINE CENTER SCATTERING CROSS SECTION
C ABSCSX . . . ABSORBER PHOTOABSORPTION CROSS SECTION

        VELT   = SQRT(2.0d0 * BOLTZ / AMASS) * SQRT(TEXO)
        CENTER = ABSIC * DBLE(FNUMBER*WAVELN) * 1.0d-8 / RTPI / VELT
        ABSCSX = DBLE(ABS_O2)

        END

