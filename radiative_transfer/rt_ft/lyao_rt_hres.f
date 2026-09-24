C Geocoronal H-Lya forward model, originally developed by James Bishop,
C modified by Jianqi Qin at the University of Illinois: V16-18
C modified by Pratik Joshi at ECE Illinois [Mar2021]

C =============================================================================

C simple driver for batch execution:  lyao_rt
      SUBROUTINE FT_MAIN(exopt,exoptx,cdens,cdensx,
     &                   th,thx,thy,
     &                   br,brx,bch,bchx,chn,chnx,
     &                   iknt_all,iknt_allx,iknt_ally,
     &                   ijknt_all,ijknt_allx,ijknt_ally,ijknt_allz,
     &                   msis_p1x,msis_p1,msis_p2x,msis_p2,
     &                   ap_histx,ap_hist,
    !  &                   rtbkgx,rtbkgy,rtbkgz,rtbkg,
     &                   thbkg,thbkgx,thbkgy,thbkgz,
     &                   zgrid,zgridx)
        IMPLICIT REAL*8 (A-H,O-Z)
        INTEGER exoptx
        REAL*8 exopt(exoptx)
        INTEGER cdensx
        REAL*8 cdens(cdensx)
        INTEGER thx,thy
        REAL*8 th(thx,thy)
        INTEGER brx
        REAL*8 br(brx)
        INTEGER bchx
        REAL*8 bch(bchx)
        INTEGER chnx
        REAL*8 chn(chnx)
        INTEGER iknt_allx,iknt_ally
        REAL*8 iknt_all(iknt_allx,iknt_ally)
        INTEGER ijknt_allx,ijknt_ally,ijknt_allz
        REAL*8 ijknt_all(ijknt_allx,ijknt_ally,ijknt_allz)
        INTEGER msis_p1x
        INTEGER msis_p1(msis_p1x)
        INTEGER msis_p2x
        REAL*8 msis_p2(msis_p2x)
        INTEGER ap_histx
        REAL*4 ap_hist(ap_histx)
        ! INTEGER rtbkgx,rtbkgy,rtbkgz
        ! REAL*8 rtbkg(rtbkgx,rtbkgy,rtbkgz)
        INTEGER thbkgx,thbkgy,thbkgz
        REAL*8 thbkg(thbkgx,thbkgy,thbkgz)
        INTEGER zgridx
        REAL*8 zgrid(zgridx)
Cf2py intent(in) exoptx
Cf2py intent(out) exopt
Cf2py depend(exoptx) exopt
Cf2py intent(in) cdensx
Cf2py intent(out) cdens
Cf2py depend(cdensx) cdens
Cf2py intent(in) thx,thy
Cf2py intent(out) th
Cf2py depend(thx,thy) th
Cf2py intent(in) brx
Cf2py intent(out) br
Cf2py depend(brx) br
Cf2py intent(in) bchx
Cf2py intent(out) bch
Cf2py depend(bchx) bch
Cf2py intent(in) chnx
Cf2py intent(out) chn
Cf2py depend(chnx) chn
Cf2py intent(in) iknt_allx,iknt_ally
Cf2py intent(out) iknt_all
Cf2py depend(iknt_allx,iknt_ally) iknt_all
Cf2py intent(in) ijknt_allx,ijknt_ally,ijknt_allz
Cf2py intent(out) ijknt_all
Cf2py depend(ijknt_allx,ijknt_ally,ijknt_allz) ijknt_all
Cf2py intent(in) msis_p1x
Cf2py intent(in) msis_p1
Cf2py intent(in) msis_p2x
Cf2py intent(in) msis_p2
Cf2py intent(in) ap_histx
Cf2py intent(in) ap_hist
Cf2py intent(in) rtbkgx,rtbkgy,rtbkgz
Cf2py intent(in) rtbkg
Cf2py intent(in) thbkgx,thbkgy,thbkgz
Cf2py intent(in) thbkg
Cf2py intent(in) zgridx
Cf2py intent(in) zgrid

      CALL global_parameters(msis_p1x,msis_p1,msis_p2x,msis_p2,
     &                       ap_histx,ap_hist,
    !  &                       rtbkgx,rtbkgy,rtbkgz,rtbkg,
     &                       thbkgx,thbkgy,thbkgz,thbkg,
     &                       zgridx,zgrid)
      CALL LYAO_RT(exoptx,exopt,cdensx,cdens,
     &             thx,thy,th,
     &             brx,br,bchx,bch,chnx,chn,
     &             iknt_allx,iknt_ally,iknt_all,
     &             ijknt_allx,ijknt_ally,ijknt_allz,ijknt_all)
      END SUBROUTINE FT_MAIN

C =============================================================================
C =============================================================================
C                   |                                |                        |
C  L Y A O _ R T    |  ATMOSPHERE & SOURCE ARRAYS    |     VERSION: 5.2       |
C  (Alpha-Omega)    |  DISTRIBUTION VERSION          |        JULY 2001       |
C                   |                                |                        |
C =============================================================================
C                                                                             |
C  SUBROUTINES NEEDED:                                                        |
C     CORONA  . . . ANALYTIC EXOSPHERE MODELS                                 |
C     ZSUN    . . . POINT-TO-SUN LOS PROPAGATOR                               |
C     ZONE    . . . ZONE-TO-ZONE LOS PROPAGATOR                               |
C     TRANS   . . . ISOTHERMAL TRANSMISSION FUNCTION INTERPOLATION ROUTINE    |
C     MATRIX_INV. . MATRIX INVERSION PROCEDURES FROM NUMERICAL RECIPES        |
C     GAUSS   . . . GAUSS-LEGENDRE POINTS & WEIGHTS                           |
C     SPEED   . . . SPEED POINTS & WEIGHTS                                    |
C                                                                             |
C =============================================================================
C
C CODE TO CALCULATE SOURCE FUNCTIONS FOR RESONANTLY SCATTERED SOLAR LYMAN
C   SERIES EMISSION LINES IN A SPHERICALLY SYMMETRIC NONISOTHERMAL
C   THERMOSPHERE & EXOSPHERE.
C BASED ON ALGORITHM OF ANDERSON AND HORD [1977].
C
C VERSION FOR GENERAL USE IN LYMAN SERIES LINE AIRGLOW MODELING:
C   -- NONISOTHERMAL TRANSPORT INCLUDED
C   -- "FINE SCALE" THERMOSPHERE GRID EVALUATED IN SUBROUTINE GET_MSIS
C   -- DEFAULT OPTION FOR MSIS [H](z):  CORRECTED FOR DIFFUSIVE FLOW
C      WITH MESOSPHERIC SOURCE REGION
C   -- CHOICE OF ORIGINAL CHAMBERLAIN OR MODIFIED ANALYTIC FORMULATION
C      FOR EXOSPHERIC EXTENSION
C   -- N2 & O2 OPACITY COMBINED AS A SINGLE ABSORBER SPECIES ("O2")
C   -- A NUMBER OF MINOR FIX-UPS & MODIFICATIONS MADE (E.G., OUTPUT FILE
C      FORMATS) TO REMOVE EXTRANEOUS diff RESULTS WRT OTHER VERSIONS

      SUBROUTINE LYAO_RT(exoptx,exopt,cdensx,cdens,
     &                   thx,thy,th,
     &                   brx,br,bchx,bch,chnx,chn,
     &                   iknt_allx,iknt_ally,iknt_all,
     &                   ijknt_allx,ijknt_ally,ijknt_allz,ijknt_all)

C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        IMPLICIT REAL*8 (A-H,O-Z)
        INTEGER exoptx
        REAL*8 exopt(exoptx)
        INTEGER cdensx
        REAL*8 cdens(cdensx)
        INTEGER thx,thy
        REAL*8 th(thx,thy)
        INTEGER brx
        REAL*8 br(brx)
        INTEGER bchx
        REAL*8 bch(bchx)
        INTEGER chnx
        REAL*8 chn(chnx)
        INTEGER iknt_allx,iknt_ally
        REAL*8 iknt_all(iknt_allx,iknt_ally)
        INTEGER ijknt_allx,ijknt_ally,ijknt_allz
        REAL*8 ijknt_all(ijknt_allx,ijknt_ally,ijknt_allz)
Cf2py intent(in) exoptx
Cf2py intent(out) exopt
Cf2py depend(exoptx) exopt
Cf2py intent(in) cdensx
Cf2py intent(out) cdens
Cf2py depend(cdensx) cdens
Cf2py intent(in) thx,thy
Cf2py intent(out) th
Cf2py depend(thx,thy) th
Cf2py intent(in) brx
Cf2py intent(out) br
Cf2py depend(brx) br
Cf2py intent(in) bchx
Cf2py intent(out) bch
Cf2py depend(bchx) bch
Cf2py intent(in) chnx
Cf2py intent(out) chn
Cf2py depend(chnx) chn
Cf2py intent(in) iknt_allx,iknt_ally
Cf2py intent(out) iknt_all
Cf2py depend(iknt_allx,iknt_ally) iknt_all
Cf2py intent(in) ijknt_allx,ijknt_ally,ijknt_allz
Cf2py intent(out) ijknt_all
Cf2py depend(ijknt_allx,ijknt_ally,ijknt_allz) ijknt_all

      PARAMETER (IKNT  =  76)
C     MUST BE EVEN
      PARAMETER (JKNT  =  64)
C     IKNT * JKNT
      PARAMETER (MDIM  = 4864)
C     OPTIONS: 4 8 16 24 32 64 96
      PARAMETER (LMU   =  24)
C     MUST BE EVEN
      PARAMETER (LPHI  =  16)
C     OPTIONS: 4 8 16
      PARAMETER (LPROF =   8)
C     EXOBASE INDEX IN RADIAL BINNING
      PARAMETER (IEXO  =  14)

C SPEED OF LIGHT
      PARAMETER (SPEEDC = 2.9979D10)
C BOLTZMANN CONSTANT
      PARAMETER (BOLTZ  = 1.3806D-16)
C INTEGRATED ABSORPTION CONSTANT
      PARAMETER (ABSIC  = 2.654D-2)
C ATOMIC MASS
      PARAMETER (AMASS  = 1.6738D-24)

      LOGICAL TRIP(IKNT,JKNT),NONIT,EXOLOS
      CHARACTER*13 CHOOZ
      CHARACTER*15 HISTORY
      CHARACTER*13 ITORNIT

      REAL*4 WAVELN, FNUMBER, BRATIO, ABS_N2, ABS_O2
      REAL*4 GLAT, GLONG, UT, STL, F107, F107A, AP(7)
      REAL*4 SATT, SATD, SCALN, ADJT, BASE, TOP
      INTEGER IDAY,IYEAR

      DIMENSION IFIBO(IKNT)
      DIMENSION BRAD(IKNT+1),RADN(IKNT),I_STEP(IKNT)
      DIMENSION BCHI(JKNT+1),CHIN(JKNT)
      DIMENSION THERMO(5,61),HDNS(IKNT),O2DNS(IKNT),TEMP(IKNT)
      DIMENSION TSOL(IKNT,JKNT),HOMO(MDIM,MDIM),SOURCE(MDIM)
      DIMENSION LUI(MDIM),ALUD(MDIM,MDIM),ASRC(MDIM),SRC(IKNT,JKNT)
      DIMENSION ODLC(2,0 : LPROF,IKNT,JKNT,2),ODO2(2,IKNT,JKNT,2),
     &          STORE(LMU,LPHI,IKNT,JKNT),OSUN(3)
      DIMENSION SPDPT(LPROF),WSPD(LPROF)
      DIMENSION THERMO_MSIS(7,61)
      DIMENSION HUT1(27), TUT1(27), O2UT1(27), Z1(27)

      COMMON/PARA_RT/   LINE_LABEL,WAVELN,FNUMBER,BRATIO,ABS_N2,ABS_O2
C      COMMON/PARA_MSIS/ GLAT,GLONG,IDAY,IYEAR,UT,STL,IAPH,AP,F107,F107A
      COMMON/PARA_EXOS/ IGEOM,SATT,SATD,SCALN,ADJT
      COMMON/H_PROFILE/ D_EXO,FLUX,D_MAX,ALT_JNT,ALT_MAX,T_EXO,MOD_MSIS
      COMMON/FINE_GRID/ THERMO_MSIS,BASE,TOP,ITHERM

      COMMON /EXOS/  GM,CONL,PLANETR,RBASE,RC,RUPR,RP,TEXO,DEXO,VELT,
     &               RADPF,CENTER,ABSCSX,FTSAT,FDSAT,IGEO

      COMMON IFIBO,I_STEP,CHOOZ

      COMMON /NMBR/  PI,RTPI,PID2,OFFSET
      COMMON /GAUS04/WW04(04),XX04(04)
      COMMON /GAUS08/WW08(08),XX08(08)
      COMMON /GAUS16/WW16(16),XX16(16)
      COMMON /GAUS24/WW24(24),XX24(24)
      COMMON /GAUS32/WW32(32),XX32(32)
      COMMON /GAUS64/WW64(64),XX64(64)
      COMMON /GAUS96/WW96(96),XX96(96)
      COMMON /SPW04/ SWGT04(04),SPNT04(04)
      COMMON /SPW08/ SWGT08(08),SPNT08(08)
      COMMON /SPW16/ SWGT16(16),SPNT16(16)

C      G-CONSTANT * PLANET MASS
C      DATA GM       / 3.9898D20/
C      MEAN PLANET RADIUS
C      DATA PLANETR  / 6371.0D5 /
C      GENERIC OFFSET VALUE
C      DATA OFFSET   / 1.0D-6   /

C      OPEN (10,FILE='LYAO_source.DAT',STATUS='UNKNOWN')

C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C *************************************************************************
C * BASIC INFORMATION  (ALTITUDES IN KM)                                  *
C *                                                                       *
C *  NONIT . . . LOGICAL FLAG FOR NONISOTHERMAL EVALUATION OF SOURCE      *
C *              FUNCTIONS (STRICTLY .TRUE. IN THIS VERSION)              *
C *              (LPROF = NUMBER OF SPEED POINTS USED IN THE EVALUATION)  *
C *  LINE_LABEL  INDEX SPECIFYING LYMAN SERIES LINE OF INTEREST:          *
C *               1 - LYMAN ALPHA                                         *
C *  BASE  . . . ALTITUDE OF BASE OF RADIATING REGION (FIXED)             *
C *  TOP   . . . TOP OF THERMOSPHERE OR EXOBASE (FIXED: BASE + 392 KM)    *
C *  ITHERM. . . NUMBER OF LAYERS IN REFERENCE THERMOSPHERE MODEL         *
C *              (QFCTR,QPIVOT CONTROL THE NONUNIFORM SPACING OF THE      *
C *              MSIS REFERENCE LEVELS, MUST BE SPECIFIED FOR EACH        *
C *              EMISSION OF INTEREST)                                    *
C *  THERMO. . . ARRAY CONTAINING "FINE GRID" THERMOSPHERIC QUANTITIES    *
C *              (ALL PURE ABSORBERS - O2 & N2 - MERGED INTO A COMPOSITE  *
C *              "O2" DENSITY)                                            *
C *                                                                       *
C *  GLAT  . . . GEODETIC LATITUDE OF EVALUATION (DEGREES)                *
C *  GLONG . . . GEOGRAPHIC LONGITUDE OF EVALUATION (DEGREES)             *
C *  IDAY  . . . OBSERVATION UT DAY NUMBER                                *
C *  IYEAR . . . OBSERVATION YEAR                                         *
C *  UT    . . . UNIVERSAL TIME OF EVALUATION                             *
C *  STL   . . . LOCAL APPARENT SOLAR TIME OF EVALUATION                  *
C *  IAPH  . . . OPTION SWITCH FOR AP-HISTORY MODE                        *
C *              IF 1, NEED TO SPECIFY AP ARRAY ELEMENTS 2-7              *
C *              (SET AP ARRAY ELEMENTS 2-7 TO 0.0 OTHERWISE)             *
C *  AP(1) . . . DAILY (AVERAGE) AP                                       *
C *  F107  . . . SOLAR F(10.7) FOR PREVIOUS 24 HOURS                      *
C *  F107A . . . SOLAR F(10.7) 81-DAY AVERAGE                             *
C *                                                                       *
C *  IGEO  . . . FLAG FOR GEOCORONAL MODEL:                               *
C *              0 - ORIGINAL CHAMBERLAIN [1963] MODEL                    *
C *              1 - MODIFIED GEOCORONAL [BISHOP 1991] MODEL              *
C *  TSAT  . . . SATELLITE COMPONENT "EXOBASE" T-PARAMETER FOR IGEO=1     *
C *              SATELLITE CRITICAL RADIUS (RE) FOR IGEO=0                *
C *  DSAT  . . . SATELLITE COMPONENT "EXOBASE" n-PARAMETER FOR IGEO=1     *
C *              IF VALUE < 0 IS ENTERED, EVAPORATIVE CASE IS EVALUATED   *
C *              IF VALUE = 0, BALLISTIC CASE IS EVALUATED                *
C *              NOT USED FOR IGEO=0                                      *
C *                                                                       *
C *************************************************************************

C GENERAL STUFF
      NONIT = .TRUE.
      CONL  = GM * AMASS / BOLTZ

      TSAT  = DBLE(SATT)
      DSAT  = DBLE(SATD)

C EFFECTIVE ZERO FOR TRANSMISSION FUNCTION
      EFZERO = EXP(-100.0D0)

C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C "FINE GRID" THERMOSPHERE MODEL & RELATED QUANTITIES

C WRITE OUT FOR FUTURE REFERENCE
      IF (IAPH.EQ.0) HISTORY = 'NO AP HISTORY  '
      IF (IAPH.EQ.1) HISTORY = 'AP HISTORY MODE'
      ITORNIT = 'NONISOTHERMAL'

      IF (MOD_MSIS .LE. 0 .AND. THERMO_MSIS(4,1) .LT. 1.0D3) THEN
        THERMO_MSIS(4,1) = 1.0D3
      END IF
      DO 102 IBIN = 1,ITHERM+1
C       z    (km)
        THERMO(1,IBIN) = THERMO_MSIS(1,IBIN)
C       r    (cm)
        THERMO(2,IBIN) = THERMO_MSIS(2,IBIN)
C       T    (K)
        THERMO(3,IBIN) = THERMO_MSIS(3,IBIN)
C       [H]  (cm^-3)
        THERMO(4,IBIN) = THERMO_MSIS(4,IBIN)
C       [O2] (cm^-3)
        THERMO(5,IBIN) = THERMO_MSIS(5,IBIN)
C       [N2] (cm^-3)
     &                   + (ABS_N2/ABS_O2)*THERMO_MSIS(6,IBIN)
C        WRITE(10,12)(THERMO(JJ,IBIN),JJ=1,5)
 102  CONTINUE

        th(1,:) = THERMO(1,:)
        th(2,:) = THERMO(2,:)
        th(3,:) = THERMO(3,:)
        th(4,:) = THERMO(4,:)
        th(5,:) = THERMO(5,:)

      BRANCH = DBLE(BRATIO)

C WRITE OUT EXOSPHERE & OPTICAL QUANTITIES
      exopt(1) = IGEO
      exopt(2) = TEXO
      exopt(3) = DEXO
      exopt(4) = TSAT
      exopt(5) = DSAT
      exopt(6) = RBASE
      exopt(7) = RC
      exopt(8) = RUPR
      exopt(9) = RP
      exopt(10) = RADPF
      exopt(11) = VELT
      exopt(12) = CENTER
      exopt(13) = ABSCSX
      exopt(14) = BRANCH

C EVALUATING ZENITH COLUMN DENSITIES
C . . . MSIS THERMOSPHERE:  PLANE PARALLEL, QUASI-ISOTHERMAL TREATMENT
      COLHT = 0.0D0
      COLOT = 0.0D0
      DO 201 ISUM = 1,ITHERM
        DZ  = THERMO(2,ISUM+1) - THERMO(2,ISUM)
        H1N = THERMO(4,ISUM)
        O2N = THERMO(5,ISUM)
        H1H = DZ / LOG(H1N/THERMO(4,ISUM+1))
        O2H = DZ / LOG(O2N/THERMO(5,ISUM+1))
        H1C = H1N*H1H * (1.0D0 - EXP(-DZ/H1H))
        O2C = O2N*O2H * (1.0D0 - EXP(-DZ/O2H))

        COLHT = COLHT + H1C
        COLOT = COLOT + O2C
 201  CONTINUE
      COLTOT = COLHT + (ABSCSX/CENTER)*COLOT
      AUPR  = CONL / TEXO / RC
      ALWR  = CONL / TEXO / RUPR
      COLHE = 0.0D0
      DO 202 ISUM = 1,16
        ALOC = ((AUPR-ALWR) * XX16(ISUM) + AUPR+ALWR) / 2.0D0
        WAL  = (AUPR-ALWR)  * WW16(ISUM) / 2.0D0
        RLOC = CONL / TEXO / ALOC
        CALL CORONA(RLOC,DNST)
        COLHE = COLHE + WAL * DNST / ALOC**2
 202  CONTINUE
      COLHE = COLHE * CONL / TEXO

C WRITE ZENITH COLUMN DENSITIES
      cdens(1) = COLHT
      cdens(2) = COLOT
      cdens(3) = COLTOT
      cdens(4) = COLHE

C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C PARTITIONING ZONES & ESTABLISHING GRIDS

C RADIAL GRID POINTS

      BRAD(1) = RBASE
      DO 301 II = 1,IKNT
        BRAD(II+1) = BRAD(1) + DBLE(I_STEP(II))*1.0D5
 301  CONTINUE
C . . . DEFINE "NOMINAL" CENTROID RADII
      DO 302 II = 1,IKNT
        RADN(II) = (BRAD(II) + BRAD(II+1)) / 2.0D0
 302  CONTINUE


C SOLAR ANGLE GRID POINTS

      BCHI(1) = 0.0D0
      BCHI(2) = ASIN(RBASE/RUPR)
      XLOC    = COS(BCHI(2))
      STEP    = XLOC / DBLE(JKNT/2 - 1)
      DO 303 JJ = 3,JKNT/2 + 1
        XLOC = XLOC - STEP
        BCHI(JJ) = ACOS(XLOC)
 303  CONTINUE

      DO 304 JJ = JKNT/2 + 1,JKNT
        JR = JKNT + 1 - JJ
        BCHI(JJ+1) = PI - BCHI(JR)
 304  CONTINUE
C . . . DEFINE "NOMINAL" CENTROID SOLAR ANGLES
      DO 305 JJ = 1,JKNT
        CHIN(JJ) = (BCHI(JJ) + BCHI(JJ+1)) / 2.0D0
 305  CONTINUE

C WRITE OUT RADIAL AND SOLAR ANGLE GRIDS
        br = BRAD
        bch = BCHI

C WRITE OUT "NOMINAL" CENTROID RADII AND SOLAR ANGLES
      iknt_all(:,1) = RADN
      chn = CHIN

C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C ZONE DENSITIES & TEMPERATURES:  SPHERICALLY SYMMETRIC

C MSIS THERMOSPHERE:  PLANE PARALLEL, QUASI-ISOTHERMAL TREATMENT
      ZSOFAR = BRAD(1)
      HCOLM  = 0.0D0
      O2COLM = 0.0D0
      TWHYD  = 0.0D0
      II = 1
      DO 401 IBIN = 1,ITHERM
        INXT  = IBIN + 1
        DZ    = THERMO(2,INXT) - THERMO(2,IBIN)
        H1N   = THERMO(4,IBIN)
        O2N   = THERMO(5,IBIN)
        H1H   = DZ / LOG(H1N/THERMO(4,INXT))
        O2H   = DZ / LOG(O2N/THERMO(5,INXT))
        H1C   = H1N*H1H * (1.0D0 - EXP(-DZ/H1H))
        O2C   = O2N*O2H * (1.0D0 - EXP(-DZ/O2H))
        TPL   = THERMO(3,IBIN)
        TGRAD = (THERMO(3,INXT) - TPL) / DZ
        ZTEST = ZSOFAR + DZ
        IF ((ZTEST.GT.BRAD(II+1)).OR.(IBIN.EQ.ITHERM)) THEN
 411      CONTINUE
          DR  = BRAD(II+1) - ZSOFAR
          H1R = H1N*H1H * (1.0D0 - EXP(-DR/H1H))
          O2R = O2N*O2H * (1.0D0 - EXP(-DR/O2H))
          TPR = TPL + DR * TGRAD
          HCOLM  = HCOLM  + H1R
          O2COLM = O2COLM + O2R
          TWHYD  = TWHYD  + H1R * (TPL + TPR) / 2.0D0
          HDNS(II)  = HCOLM  / (BRAD(II+1)-BRAD(II))
          O2DNS(II) = O2COLM / (BRAD(II+1)-BRAD(II))
          TEMP(II)  = TWHYD  / HCOLM
          IF (II.EQ.IEXO) GO TO 402
          ZSOFAR = BRAD(II+1)

          II = II + 1
          IF (ZTEST.GT.BRAD(II+1)) THEN
            HCOLM  = 0.0D0
            O2COLM = 0.0D0
            TWHYD  = 0.0D0
            H1N    = H1N * EXP(-DR/H1H)
            O2N    = O2N * EXP(-DR/O2H)
            TPL    = TPL + DR * TGRAD
            H1C    = H1C - H1R
            O2C    = O2C - O2R
            GO TO 411
          END IF
          HCOLM  = H1C - H1R
          O2COLM = O2C - O2R
          TWHYD  = HCOLM * (TPR + THERMO(3,INXT)) / 2.0D0
          ZSOFAR = ZTEST
        ELSE
          HCOLM  = HCOLM  + H1C
          O2COLM = O2COLM + O2C
          TWHYD  = TWHYD  + H1C * (TPL + THERMO(3,INXT)) / 2.0D0
          ZSOFAR = ZTEST
        END IF
 401  CONTINUE
 402  CONTINUE

C LOCAL EXOSPHERE
      DO 403 II = IEXO+1,IKNT
        AUPR = CONL / TEXO / BRAD(II)
        ALWR = CONL / TEXO / BRAD(II+1)
        SUM  = 0.0D0
        DO 413 ISUM = 1,8
          ALOC = ((AUPR-ALWR) * XX08(ISUM) + AUPR+ALWR) / 2.0D0
          WAL  = (AUPR-ALWR)  * WW08(ISUM) / 2.0D0
          RLOC = CONL / TEXO / ALOC
          CALL CORONA(RLOC,DNST)
          SUM  = SUM + WAL * DNST * (RLOC**3) / ALOC
 413    CONTINUE
        DNOM = BRAD(II+1)**3 - BRAD(II)**3
        HDNS(II)  = 3.0D0*SUM/DNOM
        O2DNS(II) = 0.0D0
        TEMP(II)  = TEXO

403   CONTINUE

C WRITE OUT DENSITY & TEMPERATURE ARRAYS
      iknt_all(:,2) = HDNS
      iknt_all(:,3) = O2DNS
      iknt_all(:,4) = TEMP

C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C CONSTRUCTING SOURCE ARRAYS

C CONVERT THERMOSPHERIC DENSITIES TO DENSITY LOGARITHMS
      DO 112 IBIN = 1,ITHERM+1
        THERMO(4,IBIN) = LOG(THERMO(4,IBIN))
        THERMO(5,IBIN) = LOG(THERMO(5,IBIN))
 112  CONTINUE

C ESTABLISH SPEED POINT ARRAYS
      DO 500 LP = 1,LPROF
        IF (LPROF.EQ.4) THEN
          SPDPT(LP) = SPNT04(LP)
          WSPD(LP)  = SWGT04(LP)
        ELSE IF (LPROF.EQ.8) THEN
          SPDPT(LP) = SPNT08(LP)
          WSPD(LP)  = SWGT08(LP)
        ELSE IF (LPROF.EQ.16) THEN
          SPDPT(LP) = SPNT16(LP)
          WSPD(LP)  = SWGT16(LP)
        END IF
 500  CONTINUE

C CALCULATING SOLAR TRANSMISSION FUNCTIONS
      DO 501 II = 1,IKNT
      DO 501 JJ = 1,JKNT
        MROW = (JJ-1)*IKNT + II
        BR1  = BRAD(II)
        BR2  = BRAD(II+1)
        BX1  = SIN(BCHI(JJ))
        BX2  = SIN(BCHI(JJ+1))
 
        IF ((BCHI(JJ+1).LE.PID2).OR.(RBASE.LE.BR1*BX2)) THEN
C . . . CENTROID POINTS
          RADS = (BRAD(II) + BRAD(II+1)) / 2.0D0
          CHIS = (BCHI(JJ) + BCHI(JJ+1)) / 2.0D0
          CALL ZSUN(ITHERM,THERMO,RADS,CHIS,OSUN,
     &              LPROF,SPDPT,WSPD,NONIT,EXOLOS)
          IF (EXOLOS) THEN
            TSOL(II,JJ) = EXP(-OSUN(2)) * TRANS(OSUN(1))
          ELSE
            TSOL(II,JJ) = EXP(-OSUN(2)) * OSUN(3)
          END IF
          IF (TSOL(II,JJ).LT.EFZERO) TSOL(II,JJ) = EFZERO

        ELSE IF ((BR1*BX2.LT.RBASE).AND.(RBASE.LT.BR1*BX1)) THEN
C . . . CENTROID POINTS
          RADS = (BRAD(II) + BRAD(II+1)) / 2.0D0
          XSM  = PI - ASIN(RBASE/RADS)
          IF (BCHI(JJ+1).LT.XSM) XSM = BCHI(JJ+1)

          CHIS = (BCHI(JJ) + XSM) / 2.0D0
          XS1  = PI - ASIN(RBASE/BR1)
          XS2  = PI - ASIN(RBASE/BR2)
          IF (BCHI(JJ+1).LT.XS2) XS2 = BCHI(JJ+1)
          VOLILL = (BR2**3 - BR1**3) * (COS(BCHI(JJ)) - COS(XS1)) +
     &             (BR2**3) * (COS(XS1) - COS(XS2)) +
     &             (RBASE**3) * (1.0D0/TAN(XS2) - 1.0D0/TAN(XS1))
          VOLTOT = (BR2**3 - BR1**3) * (COS(BCHI(JJ))-COS(BCHI(JJ+1)))
          FRAC   = VOLILL / VOLTOT
          CALL ZSUN(ITHERM,THERMO,RADS,CHIS,OSUN,
     &              LPROF,SPDPT,WSPD,NONIT,EXOLOS)
          IF (EXOLOS) THEN
            TSOL(II,JJ) = FRAC * EXP(-OSUN(2)) * TRANS(OSUN(1))
          ELSE
            TSOL(II,JJ) = FRAC * EXP(-OSUN(2)) * OSUN(3)
          END IF
          IF (TSOL(II,JJ).LT.EFZERO) TSOL(II,JJ) = EFZERO

        ELSE IF ((BR1*BX1.LE.RBASE).AND.(RBASE.LT.BR2*BX1)) THEN
C . . . CENTROID POINTS
          RSM  = RBASE / BX1
          RADS = (RSM + BRAD(II+1)) / 2.0D0
          XSM  = PI - ASIN(RBASE/RADS)
          IF (BCHI(JJ+1).LT.XSM) XSM = BCHI(JJ+1)
          CHIS = (BCHI(JJ) + XSM) / 2.0D0
          XS2  = PI - ASIN(RBASE/BR2)
          IF (BCHI(JJ+1).LT.XS2) XS2 = BCHI(JJ+1)
          VOLILL = (BR2**3) * (COS(BCHI(JJ)) - COS(XS2)) +
     &             (RBASE**3) * (1.0D0/TAN(XS2) - 1.0D0/TAN(BCHI(JJ)))
          VOLTOT = (BR2**3 - BR1**3) * (COS(BCHI(JJ))-COS(BCHI(JJ+1)))
          FRAC   = VOLILL / VOLTOT
          CALL ZSUN(ITHERM,THERMO,RADS,CHIS,OSUN,
     &              LPROF,SPDPT,WSPD,NONIT,EXOLOS)
          IF (EXOLOS) THEN
            TSOL(II,JJ) = FRAC * EXP(-OSUN(2)) * TRANS(OSUN(1))
          ELSE
            TSOL(II,JJ) = FRAC * EXP(-OSUN(2)) * OSUN(3)
          END IF
          IF (TSOL(II,JJ).LT.EFZERO) TSOL(II,JJ) = EFZERO
 
C COMPLETELY IN SHADOW
        ELSE IF (BR2*BX1.LE.RBASE) THEN
          RADS = (BRAD(II) + BRAD(II+1)) / 2.0D0
          CHIS = (BCHI(JJ) + BCHI(JJ+1)) / 2.0D0
          TSOL(II,JJ) = 0.0D0
        END IF

C STUFF INTO 1-D ARRAY
        SOURCE(MROW) = TSOL(II,JJ)
 501  CONTINUE

C WRITE OUT SOLAR TRANSMISSION FUNCTIONS
      ijknt_all(:,:,1) = TSOL

C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C CONSTRUCTING HOMOGENEOUS MATRIX (USING NOMINAL CENTROID POINTS)

C DEFINING THE HOMOGENEOUS MATRIX
C . . . MANY ENTRIES ARE ZERO, SO BEGIN WITH A ZERO MATRIX
      DO 1000 MROW = 1,MDIM
      DO 1000 MCOL = 1,MDIM
        HOMO(MROW,MCOL) = 0.0D0
 1000 CONTINUE

C AZIMUTHAL WEIGHT (CONSTANT)
      WPHI = 2.0D0 * PI / DBLE(LPHI)

C LAYING OUT INTEGRALS OVER ORIENTATION, FOR EACH BIN
C (INVOKING SYMMETRY IN CHI=PI/2 PLANE: JKNT MUST BE EVEN)
      DO 1001 II = 1,IKNT
      DO 1001 JJ = 1,JKNT/2
        MROW = (JJ-1)*IKNT + II
        R1   = RADN(II)
        X1   = CHIN(JJ)
        FCTREF = TEXO/TEMP(II)

C . . . SUMMATION SWEEPS: THETA IS LOCAL ZENITH ANGLE AND
C       PHI IS AZIMUTH ABOUT LOCAL ZENITH, ZERO IN SOLAR DIRECTION
C       (INVOKING SYMMETRY BETWEEN EAST & WEST IN PHI INTEGRAL)
        DO 1101 KMU = 1,LMU
          IF (LMU.EQ.8) THEN
            THETA = ACOS(XX08(KMU))
          ELSE IF (LMU.EQ.16) THEN
            THETA = ACOS(XX16(KMU))
          ELSE IF (LMU.EQ.24) THEN
            THETA = ACOS(XX24(KMU))
          ELSE IF (LMU.EQ.32) THEN
            THETA = ACOS(XX32(KMU))
          ELSE IF (LMU.EQ.64) THEN
            THETA = ACOS(XX64(KMU))
          ELSE IF (LMU.EQ.96) THEN
            THETA = ACOS(XX96(KMU))
          END IF
          DO 1102 KPHI = 1,LPHI
            PHI = PI * DBLE(2*KPHI-1) / DBLE(2*LPHI)

            CALL ZONE(IKNT,JKNT,BRAD,BCHI,HDNS,O2DNS,TEMP,
     &                II,JJ,R1,X1,THETA,PHI,LPROF,SPDPT,NONIT,
     &                TRIP,ODLC,ODO2,IEXIT,JEXIT,ITHERM,THERMO)

C . . . NOW CONVERTING THE LINE-CENTER OPTICAL DEPTHS TO
C       LINE-INTEGRATED TRANSMISSION FUNCTIONS, USING PREVIOUSLY
C       COMPUTED INTERPOLATION TABLE.
C       INCLUDES NONISOTHERMAL EVALUATION
            DO 1103 III = 1,IKNT
            DO 1103 JJJ = 1,JKNT
              IF (TRIP(III,JJJ)) THEN
                RHO = (ODO2(2,III,JJJ,1)-ODO2(1,III,JJJ,1))/
     &                 (ODLC(2,0,III,JJJ,1)-ODLC(1,0,III,JJJ,1))
                ABS11 = EXP(-ODO2(1,III,JJJ,1))
                ABS21 = EXP(-ODO2(2,III,JJJ,1))
                TAU1  = ODLC(1,0,III,JJJ,1)
                TAU2  = ODLC(2,0,III,JJJ,1)
                TLC11 = 0.0D0
                TLC21 = 0.0D0
                DO LP = 1,LPROF
                  XPNT  = (FCTREF - 1.0D0) * (SPDPT(LP)**2)
                  TLC11 = TLC11 + WSPD(LP) * EXP(-XPNT) *
     &                    EXP(-ODLC(1,LP,III,JJJ,1))
                  TLC21 = TLC21 + WSPD(LP) * EXP(-XPNT) *
     &                    EXP(-ODLC(2,LP,III,JJJ,1))
                END DO
                TLC11 = TLC11 * SQRT(FCTREF) * 2.0D0 / RTPI
                TLC21 = TLC21 * SQRT(FCTREF) * 2.0D0 / RTPI
                HOLDIT = ABS11*TLC11 - ABS21*TLC21
                IF (RHO.GT.OFFSET) THEN
                  SUM = 0.0D0
                  DO 1104 KK = 1,16
                    TAU  = ((TAU2-TAU1)*XX16(KK)+TAU2+TAU1)/2.0D0
                    WTAU = (TAU2-TAU1)*WW16(KK)/2.0D0
                    TABS = EXP(-(ODO2(1,III,JJJ,1) + RHO*(TAU-TAU1)))
                    FCTR  = SQRT(TEXO/TEMP(III))
                    TSCAT = 0.0D0
                    DO LP = 1,LPROF
                      XPNT  = (FCTR * SPDPT(LP))**2
                      TAULP = (TAU-TAU1) * FCTR * EXP(-XPNT)
                      XPNT  = (FCTREF - 1.0D0) * (SPDPT(LP)**2)
                      TSCAT = TSCAT + WSPD(LP) * EXP(-XPNT) *
     &                        EXP(-(TAULP + ODLC(1,LP,III,JJJ,1)))
                    END DO
                    TSCAT = TSCAT * SQRT(FCTREF) * 2.0D0 / RTPI
                    SUM = SUM + RHO * WTAU * TSCAT * TABS
 1104             CONTINUE
                  HOLDIT = HOLDIT - SUM
                END IF

                IF (ODLC(2,0,III,JJJ,2).GT.0.0D0) THEN
                  RHO = (ODO2(2,III,JJJ,2)-ODO2(1,III,JJJ,2))/
     &                   (ODLC(2,0,III,JJJ,2)-ODLC(1,0,III,JJJ,2))

                  ABS12 = EXP(-ODO2(1,III,JJJ,2))
                  ABS22 = EXP(-ODO2(2,III,JJJ,2))
                  TAU1  = ODLC(1,0,III,JJJ,2)
                  TAU2  = ODLC(2,0,III,JJJ,2)
                  TLC12 = 0.0D0
                  TLC22 = 0.0D0
                  DO LP = 1,LPROF
                    XPNT  = (FCTREF - 1.0D0) * (SPDPT(LP)**2)
                    TLC12 = TLC12 + WSPD(LP) * EXP(-XPNT) *
     &                      EXP(-ODLC(1,LP,III,JJJ,2))
                    TLC22 = TLC22 + WSPD(LP) * EXP(-XPNT) *
     &                      EXP(-ODLC(2,LP,III,JJJ,2))
                  END DO

                  TLC12 = TLC12 * SQRT(FCTREF) * 2.0D0 / RTPI
                  TLC22 = TLC22 * SQRT(FCTREF) * 2.0D0 / RTPI
                  HOLDIT = HOLDIT + ABS12*TLC12 - ABS22*TLC22
                  IF (RHO.GT.OFFSET) THEN
                    SUM = 0.0D0
                    DO 1105 KK = 1,16
                      TAU  = ((TAU2-TAU1)*XX16(KK)+TAU2+TAU1)/2.0D0
                      WTAU = (TAU2-TAU1)*WW16(KK)/2.0D0
                      TABS = EXP(-(ODO2(1,III,JJJ,2) + RHO*(TAU-TAU1)))
                      FCTR  = SQRT(TEXO/TEMP(III))
                      TSCAT = 0.0D0
                      DO LP = 1,LPROF
                        XPNT  = (FCTR * SPDPT(LP))**2
                        TAULP = (TAU-TAU1) * FCTR * EXP(-XPNT)
                        XPNT  = (FCTREF - 1.0D0) * (SPDPT(LP)**2)
                        TSCAT = TSCAT + WSPD(LP) * EXP(-XPNT) *
     &                          EXP(-(TAULP + ODLC(1,LP,III,JJJ,2)))
                      END DO
                      TSCAT = TSCAT * SQRT(FCTREF) * 2.0D0 / RTPI
                      SUM = SUM + RHO * WTAU * TSCAT * TABS
 1105               CONTINUE
                    HOLDIT = HOLDIT - SUM
                  END IF
                END IF
                IF (HOLDIT.LE.0.0D0) HOLDIT = 0.0D0
                STORE(KMU,KPHI,III,JJJ) = HOLDIT
              ELSE 
                STORE(KMU,KPHI,III,JJJ) = 0.0D0
              END IF
 1103       CONTINUE

 1102     CONTINUE
 1101   CONTINUE

C . . . NOW TO PERFORM SUMMATIONS OVER UMU AND PHI, STORING RESULTS IN
C       HOMOGENEOUS MATRIX
      DO 1201 III = 1,IKNT
      DO 1201 JJJ = 1,JKNT
        MCOL = (JJJ-1)*IKNT + III
        SUM1 = 0.0D0
        DO 1202 KMU = 1,LMU
          IF (LMU.EQ.8) THEN
            WUMU = WW08(KMU)
          ELSE IF (LMU.EQ.16) THEN
            WUMU = WW16(KMU)
          ELSE IF (LMU.EQ.24) THEN
            WUMU = WW24(KMU)
          ELSE IF (LMU.EQ.32) THEN
            WUMU = WW32(KMU)
          ELSE IF (LMU.EQ.64) THEN
            WUMU = WW64(KMU)
          ELSE IF (LMU.EQ.96) THEN
            WUMU = WW96(KMU)
          END IF
          SUM2 = 0.0D0
          DO 1203 KPHI = 1,LPHI
            SUM2 = SUM2 + WPHI * STORE(KMU,KPHI,III,JJJ)
 1203     CONTINUE
          SUM1 = SUM1 + WUMU * SUM2
 1202   CONTINUE
        HOMO(MROW,MCOL) = SUM1
 1201 CONTINUE
 
 1001 CONTINUE

C . . . NOW TO FOLD MATRIX, TO GET CHI > PI/2 ENTRIES (SPHERICAL SYMMETRY)
      DO 2001 II = 1,IKNT
      DO 2001 JJ = 1,JKNT/2
        JQ = (JKNT+1) - JJ
        MM = (JJ-1)*IKNT + II
        MQ = (JQ-1)*IKNT + II
        DO 2002 III = 1,IKNT
        DO 2002 JJJ = 1,JKNT
          JJQ = (JKNT+1) - JJJ
          NN  = (JJJ-1)*IKNT + III
          NQ  = (JJQ-1)*IKNT + III
          HOMO(MQ,NQ) = HOMO(MM,NN)
 2002   CONTINUE
 2001 CONTINUE

C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C MATRIX INVERSION

C RECASTING HOMOGENEOUS MATRIX INTO COEFFICIENT MATRIX
      DO 3001 MROW = 1,MDIM
      DO 3001 MCOL = 1,MDIM
        HOMO(MROW,MCOL) = -BRANCH*HOMO(MROW,MCOL)/(4.0D0*PI)
        IF (MROW.EQ.MCOL) HOMO(MROW,MCOL) = 1.0D0 + HOMO(MROW,MCOL)
 3001 CONTINUE

C PREPARING TO CALL INVERSION ROUTINES (COPY ARRAYS TO BE DESTROYED)
      DO 3002 MROW = 1,MDIM
        ASRC(MROW) = SOURCE(MROW)
      DO 3002 MCOL = 1,MDIM
        ALUD(MROW,MCOL) = HOMO(MROW,MCOL)
 3002 CONTINUE

C NOW TO INVERT COEFFICIENT MATRIX
      CALL LUDCMP(ALUD,MDIM,LUI,D)

      CALL LUBKSB(ALUD,MDIM,LUI,ASRC)
      CALL MPROVE(HOMO,ALUD,MDIM,LUI,SOURCE,ASRC)

C STUFF INTO 2-D ARRAY
      DO 3003 II = 1,IKNT
      DO 3003 JJ = 1,JKNT
        MROW = (JJ-1)*IKNT + II
        SRC(II,JJ) = ASRC(MROW)
        IF (SRC(II,JJ).LT.EFZERO) SRC(II,JJ) = EFZERO
 3003 CONTINUE

C WRITE OUT TOTAL SOURCE FUNCTIONS
      ijknt_all(:,:,2) = SRC

      RETURN
      END

C =============================================================================

