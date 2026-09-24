C 1022 lines
C =============================================================================
C                   |                                  |                      |
C  L Y A O _ L O S  |  APPARENT COLUMN EMISSION RATES  |     VERSION: 5.2     |
C   (Alpha-Omega)   |  DISIRIBUTION VERSION            |        JULY 2001     |
C                   |                                  |                      |
C =============================================================================
C                                                                             |
C  SUBROUTINES NEEDED:                                                        |
C     CORONA  . . . ANALYTIC EXOSPHERE MODELS                                 |
C     ZSUN    . . . POINT-TO-SUN LOS PROPAGATOR                               |
C     TRANS   . . . ISOTHERMAL TRANSMISSION FUNCTION INTERPOLATION ROUTINE    |
C     SHADOW  . . . SHADOW LOCATION AND TRANSIT:                              |
C                   GRAZE  : ITERATION ROUTINE TO DETERMINE IF SHADOW         |
C                            CROSSINGS OCCUR                                  |
C                   DARKST : ITERATION ROUTINE TO FIND POINT OF DEEPEST       |
C                            PENETRATION INTO SHADOW                          |
C                   CROSS  : ITERATION ROUTINE TO LOCATE SHADOW CROSSINGS     |
C     DLOCAL  . . . ROUTINE TO ESTIMATE LOCAL DENSITIES                       |
C     STEP_TAU. . . EVALUATES INCREMENTS TO LINE-OF-SIGHT OPTICAL DEPTHS OVER |
C                   ZONES AND DARK SEGMENTS                                   |
C     INLINE_TAU. . INCREMENTS LINE-OF-SIGHT OPTICAL DEPTHS IN TANDEM WITH    |
C                   SOURCE FUNCTION EVALUATIONS                               |
C     GAUSS   . . . GAUSS-LEGENDRE POINTS & WEIGHTS                           |
C     SPEED   . . . SPEED POINTS & WEIGHTS                                    |
C                                                                             |
C =============================================================================
C
C CODE TO CALCULATE RESONANCE RADIATION INTENSITIES ALONG SPECIFIED LINES
C   OF SIGHT (LOSs) USING lyao_rt.f SOURCE FUNCTIONS.
C BASED ON ALGORITHM OF ANDERSON AND HORD [1977, EQN. (36) IN PARTICULAR].
C
C VERSION FOR USE IN LYMAN LINE SERIES AIRGLOW MODELING:
C   -- ASSUMES OBSERVATIONS TAKEN FROM WITHIN THE MEDIUM
C   -- RESONANCE LINES FROM GROUND STATE ONLY (NO FLUORESCENCE FEATURES)

      SUBROUTINE LYAO_LOS(exoptx,exopt,cdensx,cdens,
     &                    thx,thy,th,
     &                    brx,br,bchx,bch,chnx,chn,
     &                    rhotx,rhoty,rhot,
     &                    tsx,tsy,tsz,ts,
     &                    LINE_LABEL,MAXLOS,LPROF,WAVELN,
     &                    SPDPT,WSPD,ILOS,ROBS,SZA,ZNTH,
     &                    AZI,TPALT,ACER)

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
        INTEGER rhotx,rhoty
        REAL*8 rhot(rhotx,rhoty)
        INTEGER tsx,tsy,tsz
        REAL*8 ts(tsx,tsy,tsz)

C     MUST BE SAME AS IN SOURCE CODE
      PARAMETER (IKNT   = 32)
C     MUST BE SAME AS IN SOURCE CODE
      PARAMETER (JKNT   = 32)
C     OPTIONS: 4 8 16 NOW PASSED
CCCCCCPARAMETER (LPROF  =  8)
C     CUTOFF LOS ABSORBER OPTICAL DEPTH
      PARAMETER (TLIMIT = 10.0D0)

C     SPEED OF LIGHT
      PARAMETER (SPEEDC = 2.9979D10)
C     BOLTZMANN CONSTANT
      PARAMETER (BOLTZ  = 1.3806D-16)
C     ATOMIC MASS
      PARAMETER (AMASS  = 1.6738D-24)

      LOGICAL PASS,ILLUM,NONIT,EXOLOS
      LOGICAL :: DEBUG = .FALSE.
      CHARACTER*13 CHOOZ
      CHARACTER*15 HISTORY
      CHARACTER*13 ITORNIT

      DIMENSION ROBS(MAXLOS), SZA(MAXLOS), ZNTH(MAXLOS), AZI(MAXLOS)
      DIMENSION TPALT(MAXLOS), ACER(3,MAXLOS)
      DIMENSION BRAD(IKNT+1),RADN(IKNT)
      DIMENSION BCHI(JKNT+1),CHIN(JKNT)
      DIMENSION THERMO(5,61),HDNS(IKNT),O2DNS(IKNT),TEMP(IKNT)
      DIMENSION SOL(IKNT,JKNT),SRC(IKNT,JKNT),OSUN(3)
      DIMENSION SPDPT(LPROF),WSPD(LPROF)
      DIMENSION TAULPS(0 : 16),TAULPM(0 : 16),TLPLOC(0 : 16)
      DIMENSION XSMS(IKNT+2),D0SMS(IKNT+2,JKNT),D2SMS(IKNT+2,JKNT),
     &          YSMS(IKNT+2),D2YSMS(IKNT+2)

       REAL*8 THERMO_MSIS(7,61)
       REAL*4 BASE, TOP

      COMMON/FINE_GRID/ THERMO_MSIS,BASE,TOP,ITHERM
      COMMON /GAUS08/WW08(08),XX08(08)
      COMMON /SPW04/ SWGT04(04),SPNT04(04)
      COMMON /SPW08/ SWGT08(08),SPNT08(08)
      COMMON /SPW16/ SWGT16(16),SPNT16(16)

C     G-CONSTANT * PLANET MASS
      DATA GM       / 3.9898D20/
C     MEAN PLANET RADIUS
      DATA PLANETR  / 6371.0D5 /
C     GENERIC OFFSET VALUE
      DATA OFFSET   / 1.0D-6   /

C GENERAL STUFF
      NONIT = .TRUE.
      CONL  = GM * AMASS / BOLTZ

      PID2  = ASIN(1.0D0)
      PI    = 2.0D0 * PID2
      RTPI  = SQRT(PI)

C     single scattering coefficients
      ANSCAT_A = 0.9219D0
      ANSCAT_B = 0.2344D0

C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        THERMO = th

      DO 112 IBIN = 1,ITHERM+1
        THERMO(4,IBIN) = LOG(THERMO(4,IBIN))
        THERMO(5,IBIN) = LOG(THERMO(5,IBIN))
 112  CONTINUE

C Read EXOSPHERE & OPTICAL QUANTITIES
       IGEO = exopt(1)
       TEXO = exopt(2)
       DEXO = exopt(3)
       TSAT = exopt(4)
       DSAT = exopt(5)
       RBASE = exopt(6)
       RC = exopt(7)
       RUPR = exopt(8)
       RP = exopt(9)
       RADPF = exopt(10)
       VELT = exopt(11)
       CENTER = exopt(12)
       ABSCSX = exopt(13)
       BRANCH = exopt(14)

       IF (IGEO.EQ.1) THEN
        FTSAT = TSAT/TEXO
        FDSAT = DSAT/DEXO
      ELSE IF (IGEO.EQ.0) THEN
        FTSAT = TSAT * PLANETR
        FDSAT = 0.0D0
      END IF

C Read EFFECTIVE ZENITH COLUMN DENSITIES
      COLHT = cdens(1)
      COLOT = cdens(2)
      COLTOT = cdens(3)
      COLHE = cdens(4)

C Read BOUNDARY GRID POINTS
       BRAD = br
       BCHI = bch

C Read NOMINAL CENTROID POINTS
        RADN = rhot(:,1)
        CHIN = chn

C Read DENSITY & TEMPERATURE ARRAYS
        HDNS = rhot(:,2)
        O2DNS = rhot(:,3)
        TEMP = rhot(:,4)

C Read SOLAR SOURCE ARRAY
      SOL = ts(:,:,1)

C Read TOTAL SOURCE ARRAY
      SRC = ts(:,:,2)

C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C REPLACE BRANCHING RATIO USED IN SOURCE FUNCTION EVALUATIONS WITH
C BRANCHING RATIO FOR FLUORESCENT SCATTERING

C H-alpha photons escape freely; use Ly-beta->H-alpha branching ratio
      BRANCH = 0.118D0 ! Balmer-Alpha B-Ratio (n=3->2 decay, 11.8% of Ly-beta absorptions)

C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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

C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C SET UP INTERPOLATION ARRAYS FOR MULTIPLE SCATTERING SOURCE FUNCTIONS

      CALL SMS_SPLINE(IKNT,JKNT,BRAD,RADN,SOL,SRC,XSMS,D0SMS,D2SMS)

C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C LOOPING THROUGH LINES OF SIGHT

C ZSUN SHADOW OFFSET (REDUNDANT ALTITUDE GRID VIOLATION PROTECTION)
C   SOLILLUMCHK . . . Solar Illumination Check
C                     effective placement of terminator shadow height,
C                     slightly offset from nominal placement to protect
C                     against ZSUN grid violations.
      SOLILLUMCHK = BRAD(1) * (1.0D0 + OFFSET)

C LINES OF SIGHT ARE INDEPENDENT: SPREAD THEM OVER OPENMP THREADS
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,4) DEFAULT(NONE)
!$OMP& SHARED(ILOS,ROBS,SZA,ZNTH,AZI,TPALT,ACER,RBASE,RUPR,PI,PID2,
!$OMP&        PLANETR,BRAD,BCHI,CHIN,THERMO,ITHERM,LPROF,SPDPT,WSPD,
!$OMP&        NONIT,XSMS,D0SMS,D2SMS,BRANCH,CENTER,ANSCAT_A,ANSCAT_B,
!$OMP&        SOLILLUMCHK,/GAUS08/)
!$OMP& PRIVATE(AA,ANSCAT,ANUM,BB,CP1,CP2,CSCAT,CUMU,CX1,CX2,DELT2,
!$OMP&         DELTA,DL,DLWR,DNOM,DSHAD,DSHAD2,DSMU,DUPR,DWAS,ICROSS,
!$OMP&         IDEL,IEXIT,II,III,IIN,IL,ILLUM,ILOOP,IPHASE,IPHASX,
!$OMP&         ISHAD,ISPLN,ITRIP,JBOX1,JBOX2,JJ,JJJ,JJN,JLOOP,LP,OEXIT,
!$OMP&         OMEGA,OMEGL,OMEGR,OMEGX,OSKIM,OTOT,PASS,PHI,R1,R2,R2R,
!$OMP&         RAT,RATIOX,RL,RLWR,RSINXK,RSMU,SMS,SMS1,SMS2,SP1,SP2,
!$OMP&         SSOL,SUMU,SX1,SX2,TAUO2M,TAUO2S,THET2,THETA,THETR,
!$OMP&         TO2LOC,WDL,X1,X2,X2X,XEXIT,XL,XSKIM,XTOT,ZETA,HWAS,
!$OMP&         O2WAS,TWAS,HLOC,O2LOC,TLOC,OGRAZ,PHIG,ODEEP,EXOLOS,
!$OMP&         OSUN,TLPLOC,TAULPS,TAULPM,YSMS,D2YSMS)
      DO 3333 LOOK = 1,ILOS

C LOCATION AND LOOK DIRECTION IN SOLAR COORDINATES
C protections added, required by use of spherical trig relations
C keep in mind that AZI is an ``interior'' angle
C   R1    . . . RADIUS AT OBSERVATION LOCATION
      R1 = ROBS(LOOK)

      IF (R1.LT.(RBASE*1.001D0)) THEN
        R1 = RBASE * 1.001D0
      ELSE IF (R1.GT.(RUPR*0.999D0)) THEN
        R1 = RUPR * 0.999D0
      END IF
C   X1    . . . SOLAR ZENITH ANGLE AT OBSERVATION LOCATION
      X1 = SZA(LOOK) * PI / 180.0D0
      IF (X1.LT.(0.001D0*PI)) THEN
        X1 = 0.001D0 * PI
      ELSE IF (X1.GT.(0.999D0*PI)) THEN
        X1 = 0.999D0 * PI
      END IF
C   THETA . . . LOCAL ZENITH ANGLE OF LOS
      THETA = ZNTH(LOOK) * PI / 180.0D0
C      IF (THETA.LT.(0.001D0*PI)) THEN
C        THETA = 0.001D0 * PI
C      ELSE IF (THETA.GT.(0.999D0*PI)) THEN
C        THETA = 0.999D0 * PI
C      END IF
C   PHI   . . . AZIMUTH OF LOS RELATIVE TO SOLAR DIRECTION
      PHI = AZI(LOOK) * PI / 180.0D0
      IF (PHI.LT.0.0D0) PHI = -PHI
      IF (PHI.LT.(0.001D0*PI)) THEN
        PHI = 0.001D0 * PI
      ELSE IF ((PHI.LE.PID2).AND.(PHI.GE.(PID2-0.001D0*PI))) THEN
        PHI = PID2 - 0.001D0 * PI
      ELSE IF ((PHI.GT.PID2).AND.(PHI.LE.(PID2+0.001D0*PI))) THEN
        PHI = PID2 + 0.001D0 * PI
      ELSE IF (PHI.GT.(0.999D0*PI)) THEN
        PHI = 0.999D0 * PI
      END IF

C DETERMINING INITIAL BINS
      DO 3011 ILOOP = 1,IKNT
        IF ((R1.GE.BRAD(ILOOP)).AND.(R1.LT.BRAD(ILOOP+1))) THEN
          II = ILOOP
          GO TO 3012
        END IF
 3011 CONTINUE
 3012 CONTINUE
      DO 3021 JLOOP = 1,JKNT
        IF ((X1.GE.BCHI(JLOOP)).AND.(X1.LT.BCHI(JLOOP+1))) THEN
          JJ = JLOOP
          GO TO 3022
        END IF
 3021 CONTINUE
 3022 CONTINUE

C LINE OF SIGHT REFERENCE QUANTITIES
C   RSMU  . . . TANGENT POINT RADIUS
C   DSMU  . . . DISTANCE ALONG LOS TO TANGENT POINT
C   XSKIM . . . "SKIMMING" LATITUDE OF LOS GREAT CIRCLE
C   XTOT  . . . ASYMPTOTIC VALUE OF CHI
C   OSKIM . . . DISPLACEMENT ANGLE TO SKIMMING LATITUDE
C   OTOT  . . . ASYMPTOTIC VALUE OF OMEGA
C   IPHASE. . . PHASE OF PHI2, FOR COSINE EVALUATIONS
C   PASS  . . . FLAG FOR XSKIM PASSAGE

      CX1  = COS(X1)
      SX1  = SIN(X1)
      CUMU = COS(THETA)
      SUMU = SIN(THETA)
      CP1  = COS(PHI)
      SP1  = SIN(PHI)
      RSMU = R1*SUMU
      DSMU = R1*ABS(CUMU)
      IF (PHI.LE.PID2) THEN
        XSKIM  = ASIN(SX1*SP1)
        IPHASE = -1
      ELSE
        XSKIM  = PI - ASIN(SX1*SP1)
        IPHASE = +1
      END IF
      XTOT  = ACOS(CX1*CUMU + SX1*SUMU*CP1)
      ANUM  = CX1*COS(XSKIM)
      DNOM  = 1.0D0 - SX1*SIN(XSKIM)*SP1

      RAT=ANUM/DNOM
      IF (RAT>1.0) THEN
        RAT=1.0
      ELSE IF (RAT<-1.0) THEN
        RAT=-1.0
      END IF
      OSKIM = ACOS(RAT)

      OTOT  = THETA
      IF (OSKIM.LT.OTOT) THEN
        PASS = .TRUE.
      ELSE
        PASS = .FALSE.
      END IF
C . . . TANGENT POINT ALTITUDES FOR OUTPUT
      IF (ZNTH(LOOK).GT.90.0D0) THEN
        TPALT(LOOK) = (RSMU - PLANETR)*1.0D-5
      ELSE
        TPALT(LOOK) = (ROBS(LOOK) - PLANETR)*1.0D-5
      END IF

C SCATTERING FACTOR FOR LOS (NOTE THAT SCATTERING ANGLE IS XTOT)
      CSCAT  = COS(XTOT)
      ANSCAT = ANSCAT_A + ANSCAT_B * (CSCAT**2)

C SHADOW ZONE
      DSHAD  = 0.0D0
      DSHAD2 = 0.0D0
      IF ((X1.GT.PID2).AND.((R1*SX1).LT.RBASE)) THEN
        ILLUM = .FALSE.
      ELSE
        ILLUM = .TRUE.
      END IF
C . . . NOW DETERMINE LOS RADIAL EXIT BOUNDARY
C   OEXIT . . . DISPLACEMENT ANGLE TO EXIT
C   XEXIT . . . SOLAR ANGLE AT EXIT
      IF ((CUMU.LT.0.0D0).AND.(RSMU.LT.RBASE)) THEN
        IEXIT = 1

        RAT=RSMU/BRAD(IEXIT)
C       RAT outside ASIN limits gives NaN in Radiance
        IF (RAT>1.0) THEN
          RAT=1.0
        ELSE IF (RAT<-1.0) THEN
          RAT=-1.0
        END IF

        OEXIT = THETA - (PI - ASIN(RAT))
      ELSE
        IEXIT = IKNT+1

        RAT=RSMU/BRAD(IEXIT)
C       RAT outside ASIN limits gives NaN in Radiance
        IF (RAT>1.0) THEN
          RAT=1.0
        ELSE IF (RAT<-1.0) THEN
          RAT=-1.0
        END IF

        OEXIT = THETA - ASIN(RAT)
      END IF
      XEXIT = ACOS(CX1*COS(OEXIT)+SX1*SIN(OEXIT)*CP1)
C . . . NOW DETERMINE DISTANCES TO CROSSINGS.
      IF (.NOT.ILLUM) THEN
        IF ((IEXIT.EQ.1).OR.(XEXIT.GE.BCHI(JKNT))) THEN
          ISHAD = 0
        ELSE
          ISHAD = 1
          CALL CROSS(RBASE,RSMU,R1,THETA,CX1,SX1,CP1,ILLUM,OEXIT,
     &               ISHAD,DSHAD,DSHAD2)
        END IF
      ELSE

        ISHAD = 0
        IF (IEXIT.EQ.1) THEN
          IF (XEXIT.GT.PID2) THEN
            ISHAD = 1
            CALL CROSS(RBASE,RSMU,R1,THETA,CX1,SX1,CP1,ILLUM,OEXIT,
     &                 ISHAD,DSHAD,DSHAD2)
          END IF
        ELSE IF (XEXIT.GE.BCHI(JKNT)) THEN
          ISHAD = 1
          CALL CROSS(RBASE,RSMU,R1,THETA,CX1,SX1,CP1,ILLUM,OEXIT,
     &               ISHAD,DSHAD,DSHAD2)
        ELSE IF (((X1+THETA).GT.PI).AND.
     &           (.NOT.((X1.LE.PID2).AND.(PHI.LE.PID2)))) THEN
          RSINXK = RSMU * SIN(XSKIM)
          ZETA   = X1 + THETA - PI
          BB     = R1 * SX1 / SIN(ZETA)
          AA     = BB * COS(ZETA) - R1 * CX1
          IF ((AA.LT.RUPR).AND.(RSINXK.LT.RBASE)) THEN
            CALL GRAZE(AA,BB,RBASE,R1,X1,THETA,CX1,SX1,CUMU,SUMU,
     &                 OGRAZ,PHIG)
            IF (PHI.GT.PHIG) THEN
              ISHAD = 2
              CALL DARKST(OGRAZ,R1,X1,THETA,CX1,SX1,CUMU,SUMU,CP1,ODEEP)
              CALL CROSS(RBASE,RSMU,R1,THETA,CX1,SX1,CP1,ILLUM,ODEEP,
     &                   ISHAD,DSHAD,DSHAD2)
            END IF
          END IF
        END IF
      END IF

C + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
C PROPAGATION ALONG LINE OF SIGHT

C INITIALIZE LOS ARRAYS
      DO 3031 IL = 1,3
        ACER(IL,LOOK) = 0.0D0
 3031 CONTINUE
      TAUO2S = 0.0D0
      TAUO2M = 0.0D0
      DO LP = 0,LPROF
        TAULPS(LP) = 0.0D0
        TAULPM(LP) = 0.0D0
      END DO

C SEARCH FOR EDGE OF INITIAL BIN ALONG SPECIFIED LINE OF SIGHT
      ITRIP = 1
C . . . RADIAL AND SOLAR ANGLE BIN
      IF ((CUMU.GT.0.0D0).OR.(RSMU.GT.BRAD(II))) THEN
        R2R   = BRAD(II+1)
        RAT=RSMU/R2R
C       RAT outside ASIN limits gives NaN in Radiance
        IF (RAT>1.0) THEN
          RAT=1.0
        ELSE IF (RAT<-1.0) THEN
          RAT=-1.0
        END IF

        THETR = ASIN(RAT)
        III   = II + 1
      ELSE
        R2R   = BRAD(II)
        RAT=RSMU/R2R
C       RAT outside ASIN limits gives NaN in Radiance
        IF (RAT>1.0) THEN
          RAT=1.0
        ELSE IF (RAT<-1.0) THEN
          RAT=-1.0
        END IF

        THETR = PI - ASIN(RAT)
        III   = II - 1
      END IF
      OMEGR = THETA - THETR
      IF (PASS) THEN
        IF (IPHASE.EQ.+1) THEN
          IF (XSKIM.GT.BCHI(JJ+1)) THEN
            X2X    = BCHI(JJ+1)
            IPHASX = IPHASE
            JJJ    = JJ + 1
          ELSE IF (XTOT.LT.BCHI(JJ)) THEN
            X2X    = BCHI(JJ)
            IPHASX = -IPHASE
            JJJ    = JJ - 1
          ELSE
            X2X    = XTOT
            IPHASX = -IPHASE
            JJJ    = JJ
          END IF
        ELSE IF (IPHASE.EQ.-1) THEN
          IF (XSKIM.LT.BCHI(JJ)) THEN
            X2X    = BCHI(JJ)
            IPHASX = IPHASE
            JJJ    = JJ - 1
          ELSE IF (XTOT.GT.BCHI(JJ+1)) THEN
            X2X    = BCHI(JJ+1)
            IPHASX = -IPHASE
            JJJ    = JJ + 1
          ELSE
            X2X    = XTOT
            IPHASX = -IPHASE
            JJJ    = JJ
          END IF
        END IF
      ELSE
        IPHASX = IPHASE
        IF (XTOT.GT.BCHI(JJ+1)) THEN
          X2X = BCHI(JJ+1)
          JJJ = JJ + 1
        ELSE IF (XTOT.LT.BCHI(JJ)) THEN
          X2X = BCHI(JJ)
          JJJ = JJ - 1
        ELSE
          X2X = XTOT
          JJJ = JJ
        END IF
      END IF
      SX2   = SIN(X2X)
      CX2   = COS(X2X)
      SP2   = SP1*SX1/SX2
      IF (SP2.GT.1.0D0) SP2 = 1.0D0
      CP2   = IPHASX * SQRT(1.0D0 - SP2**2)
      ANUM  = CX1*CX2 - SX1*SX2*CP1*CP2
      DNOM  = 1.0D0 - SX1*SX2*SP1*SP2

      RAT=ANUM/DNOM
      IF (RAT>1.0) THEN
        RAT=1.0
      ELSE IF (RAT<-1.0) THEN
        RAT=-1.0
      END IF
      OMEGX = ACOS(RAT)

C . . . "NEAREST" CROSSING
      IF (OMEGX.LT.OMEGR) THEN
        OMEGA = OMEGX
        THET2 = THETA - OMEGA
        R2    = RSMU/SIN(THET2)
        X2    = X2X
        III   = II
      ELSE
        OMEGA = OMEGR
        THET2 = THETR
        R2    = R2R
        X2    = ACOS(CX1*COS(OMEGA)+SX1*SIN(OMEGA)*CP1)
        JJJ   = JJ
      END IF
      IF ((OMEGA.GE.OSKIM).AND.(PASS)) THEN
        IPHASE = -IPHASE
        PASS   = .FALSE.
      END IF
      DELTA = SQRT(R2**2 + R1**2 - 2.0D0*R2*R1*COS(OMEGA))

C -----------------------------------------------------------------------------
C NOW TO INTEGRATE ALONG LOS INSIDE BIN

      IF ((II.GE.1).AND.(II.LE.IKNT)) THEN

C SINGLY-SCATTERED COMPONENT
        SSOL = 0.0D0

        ICROSS = 0
        IF (ISHAD.EQ.0) THEN
          IF (.NOT.ILLUM) GO TO 1003
          DUPR = DELTA
          DLWR = 0.0D0
        ELSE IF (ISHAD.EQ.1) THEN
          IF (.NOT.ILLUM) THEN
            IF (DELTA.LE.DSHAD) THEN
              CALL STEP_TAU(ITHERM,THERMO,R1,X1,THETA,DELTA,0.0D0,
     &                      LPROF,SPDPT,NONIT,TO2LOC,TLPLOC)
              TAUO2S = TO2LOC
              DO LP = 0,LPROF
                TAULPS(LP) = TLPLOC(LP)
              END DO
              GO TO 1003
            END IF
            DUPR = DELTA
            DLWR = DSHAD
          ELSE
            IF (DELTA.LE.DSHAD) THEN
              DUPR = DELTA
            ELSE
              DUPR = DSHAD
              ICROSS = 1
            END IF
            DLWR = 0.0D0
          END IF
        ELSE IF (ISHAD.EQ.2) THEN
          IF (DELTA.LE.DSHAD) THEN
            DUPR = DELTA
          ELSE
            DUPR = DSHAD
            ICROSS = 1
            IF (DELTA.GT.DSHAD2) ICROSS = 2
          END IF
          DLWR = 0.0D0
        END IF
        IF (DLWR.GT.0.0D0) THEN
          CALL STEP_TAU(ITHERM,THERMO,R1,X1,THETA,DLWR,0.0D0,
     &                  LPROF,SPDPT,NONIT,TO2LOC,TLPLOC)
          TAUO2S = TO2LOC
          DO LP = 0,LPROF
            TAULPS(LP) = TLPLOC(LP)
          END DO
        END IF
 1002   CONTINUE
        RLWR = SQRT(R1**2 + DLWR**2 - 2.0D0*R1*DLWR*COS(PI-THETA))
        DWAS = DLWR
        CALL DLOCAL(ITHERM,THERMO,RLWR,X1,HWAS,O2WAS,TWAS)
        TO2LOC = 0.0D0
        DO LP = 0,LPROF
          TLPLOC(LP) = 0.0D0
        END DO
        DO 1001 IDEL = 1,8
          DL  = ((DUPR-DLWR)*XX08(IDEL) + DUPR+DLWR) / 2.0D0
          WDL =  (DUPR-DLWR)*WW08(IDEL) / 2.0D0
          RL  = SQRT(R1**2 + DL**2 - 2.0D0*R1*DL*COS(PI-THETA))

          RAT=RSMU/RL
C         RAT outside ASIN limits gives NaN in Radiance
          IF (RAT>1.0) THEN
            RAT=1.0
          ELSE IF (RAT<-1.0) THEN
            RAT=-1.0
          END IF

          IF ((CUMU.LT.0.0D0).AND.(DL.LT.DSMU)) THEN
            OMEGL = THETA - (PI - ASIN(RAT))
          ELSE
            OMEGL = THETA - ASIN(RAT)
          END IF
          XL  = ACOS(CX1*COS(OMEGL)+SX1*SIN(OMEGL)*CP1)
          CALL DLOCAL(ITHERM,THERMO,RL,X1,HLOC,O2LOC,TLOC)

          IF (((RL*SIN(XL)).LE.SOLILLUMCHK).AND.(XL.GE.PID2)) THEN
            SSOL = 0.0D0
          ELSE
            CALL ZSUN(ITHERM,THERMO,RL,XL,OSUN,
     &                LPROF,SPDPT,WSPD,NONIT,EXOLOS)
            IF (EXOLOS) THEN
              SSOL = BRANCH*CENTER*HLOC*EXP(-OSUN(2))*TRANS(OSUN(1))
            ELSE
              SSOL = BRANCH*CENTER*HLOC*EXP(-OSUN(2))*OSUN(3)
            END IF
          END IF

          CALL INLINE_TAU(DWAS,HWAS,O2WAS,TWAS,
     &                    DL  ,HLOC,O2LOC,TLOC,
     &                    LPROF,SPDPT,NONIT,TO2LOC,TLPLOC)
C H-alpha photons travel freely to observer: TABS=1, TSCAT=1
    !       TABS   = EXP(-(TAUO2S + TO2LOC))
    !       FCTREF = TEXO/TLOC
    !       TSCAT  = 0.0D0
    !       DO LP = 1,LPROF
    !         XPNT  = (FCTREF - 1.0D0) * (SPDPT(LP)**2)
    !         TSCAT = TSCAT + WSPD(LP) *
    !  &              EXP(-(TAULPS(LP)+TLPLOC(LP))) * EXP(-XPNT)
    !       END DO
    !       TSCAT  = TSCAT * SQRT(FCTREF) * 2.0D0 / RTPI
          ACER(1,LOOK) = ACER(1,LOOK) + WDL*ANSCAT*SSOL ! *TABS*TSCAT
 1001   CONTINUE
        CALL STEP_TAU(ITHERM,THERMO,R1,X1,THETA,DUPR,DLWR,
     &                LPROF,SPDPT,NONIT,TO2LOC,TLPLOC)
        TAUO2S = TAUO2S + TO2LOC
        DO LP = 0,LPROF
          TAULPS(LP) = TAULPS(LP) + TLPLOC(LP)
        END DO
        IF (ICROSS.EQ.1) THEN
          ICROSS = 0
          CALL STEP_TAU(ITHERM,THERMO,R1,X1,THETA,DELTA,DSHAD,
     &                  LPROF,SPDPT,NONIT,TO2LOC,TLPLOC)
          TAUO2S = TAUO2S + TO2LOC
          DO LP = 0,LPROF
            TAULPS(LP) = TAULPS(LP) + TLPLOC(LP)
          END DO
        ELSE IF (ICROSS.EQ.2) THEN
          ICROSS = 0
          CALL STEP_TAU(ITHERM,THERMO,R1,X1,THETA,DSHAD2,DSHAD,
     &                  LPROF,SPDPT,NONIT,TO2LOC,TLPLOC)
          TAUO2S = TAUO2S + TO2LOC
          DO LP = 0,LPROF
            TAULPS(LP) = TAULPS(LP) + TLPLOC(LP)
          END DO
          DUPR = DELTA
          DLWR = DSHAD2
          GO TO 1002
        END IF

 1003   CONTINUE
C MULTIPLY-SCATTERED COMPONENT

        DWAS = 0.0D0
        CALL DLOCAL(ITHERM,THERMO,R1,X1,HWAS,O2WAS,TWAS)
        TO2LOC = 0.0D0
        DO LP = 0,LPROF
          TLPLOC(LP) = 0.0D0
        END DO
        DO 1004 IDEL = 1,8
          DL  = (XX08(IDEL) + 1.0D0) * DELTA / 2.0D0
          WDL =  WW08(IDEL) * DELTA / 2.0D0
          RL  = SQRT(R1**2 + DL**2 - 2.0D0*R1*DL*COS(PI-THETA))

          RAT=RSMU/RL
C         RAT outside ASIN limits gives NaN in Radiance
          IF (RAT>1.0) THEN
            RAT=1.0
          ELSE IF (RAT<-1.0) THEN
            RAT=-1.0
          END IF

          IF ((CUMU.LT.0.0D0).AND.(DL.LT.DSMU)) THEN
            OMEGL = THETA - (PI - ASIN(RAT))
          ELSE
            OMEGL = THETA - ASIN(RAT)
          END IF
          XL  = ACOS(CX1*COS(OMEGL)+SX1*SIN(OMEGL)*CP1)
          CALL DLOCAL(ITHERM,THERMO,RL,X1,HLOC,O2LOC,TLOC)

          IF ((XL.GT.CHIN(1)).AND.(XL.LT.CHIN(JKNT))) THEN
            IF (XL.LE.CHIN(JJ)) THEN
              JBOX1 = JJ - 1
              JBOX2 = JJ
            ELSE
              JBOX1 = JJ
              JBOX2 = JJ + 1
            END IF
            DO ISPLN = 1,IKNT+2
              YSMS(ISPLN)   = D0SMS(ISPLN,JBOX1)
              D2YSMS(ISPLN) = D2SMS(ISPLN,JBOX1)
            END DO

            CALL SPLINT_LYAO(XSMS,YSMS,D2YSMS,IKNT+2,RL,SMS1)
            DO ISPLN = 1,IKNT+2
              YSMS(ISPLN)   = D0SMS(ISPLN,JBOX2)
              D2YSMS(ISPLN) = D2SMS(ISPLN,JBOX2)
            END DO
            CALL SPLINT_LYAO(XSMS,YSMS,D2YSMS,IKNT+2,RL,SMS2)
            RATIOX = (XL - CHIN(JBOX1))/(CHIN(JBOX2) - CHIN(JBOX1))
            SMS    = BRANCH*CENTER*HLOC*EXP(SMS1 + RATIOX*(SMS2-SMS1))
          ELSE
            DO ISPLN = 1,IKNT+2
              YSMS(ISPLN)   = D0SMS(ISPLN,JJ)
              D2YSMS(ISPLN) = D2SMS(ISPLN,JJ)
            END DO
            CALL SPLINT_LYAO(XSMS,YSMS,D2YSMS,IKNT+2,RL,SMS)
            SMS    = BRANCH*CENTER*HLOC*EXP(SMS)
          END IF

          CALL INLINE_TAU(DWAS,HWAS,O2WAS,TWAS,
     &                    DL  ,HLOC,O2LOC,TLOC,
     &                    LPROF,SPDPT,NONIT,TO2LOC,TLPLOC)
C H-alpha photons travel freely to observer: TABS=1, TSCAT=1
    !       TABS   = EXP(-TO2LOC)
    !       FCTREF = TEXO/TLOC
    !       TSCAT  = 0.0D0
    !       DO LP = 1,LPROF
    !         XPNT  = (FCTREF - 1.0D0) * (SPDPT(LP)**2)
    !         TSCAT = TSCAT + WSPD(LP) *
    !  &              EXP(-TLPLOC(LP)) * EXP(-XPNT)
    !       END DO
    !       TSCAT  = TSCAT * SQRT(FCTREF) * 2.0D0 / RTPI
          ACER(2,LOOK) = ACER(2,LOOK) + WDL*SMS ! *TABS*TSCAT
 1004   CONTINUE
        CALL STEP_TAU(ITHERM,THERMO,R1,X1,THETA,DELTA,0.0D0,
     &                LPROF,SPDPT,NONIT,TO2LOC,TLPLOC)
        TAUO2M = TO2LOC
        DO LP = 0,LPROF
          TAULPM(LP) = TLPLOC(LP)
        END DO

C TOTAL INTENSITY
        ACER(3,LOOK) = ACER(1,LOOK) + ACER(2,LOOK)

      END IF
C -----------------------------------------------------------------------------
C     exiting the volume
      IF ((III.LT.1).OR.(III.GT.IKNT)) GO TO 2222
C     effective end-of-the-line
      IF (TAUO2M.GT.TLIMIT) GO TO 2222

 1111 CONTINUE

C SUBSEQUENT STEPS THROUGH BINS ALONG LINE OF SIGHT
      ITRIP = ITRIP + 1
C . . . RADIAL AND SOLAR ANGLE BIN
      IF ((COS(THET2).GT.0.0D0).OR.(RSMU.GT.BRAD(III))) THEN
        R2R   = BRAD(III+1)

        RAT=RSMU/R2R
C       RAT outside ASIN limits gives NaN in Radiance
        IF (RAT>1.0) THEN
          RAT=1.0
        ELSE IF (RAT<-1.0) THEN
          RAT=-1.0
        END IF

        THETR = ASIN(RAT)
        IIN   = III + 1
      ELSE
        R2R   = BRAD(III)

        RAT=RSMU/R2R
C       RAT outside ASIN limits gives NaN in Radiance
        IF (RAT>1.0) THEN
          RAT=1.0
        ELSE IF (RAT<-1.0) THEN
          RAT=-1.0
        END IF

        THETR = PI - ASIN(RAT)
        IIN   = III - 1
      END IF
      OMEGR = THETA - THETR

      IF (PASS) THEN
        IF (IPHASE.EQ.+1) THEN
          IF (XSKIM.GT.BCHI(JJJ+1)) THEN
            X2X    = BCHI(JJJ+1)
            IPHASX = IPHASE
            JJN    = JJJ + 1
          ELSE IF (XTOT.LT.BCHI(JJJ)) THEN
            X2X    = BCHI(JJJ)
            IPHASX = -IPHASE
            JJN    = JJJ - 1
          ELSE
            X2X    = XTOT
            IPHASX = -IPHASE
            JJN    = JJJ
          END IF
        ELSE IF (IPHASE.EQ.-1) THEN
          IF (XSKIM.LT.BCHI(JJJ)) THEN
            X2X    = BCHI(JJJ)
            IPHASX = IPHASE
            JJN    = JJJ - 1
          ELSE IF (XTOT.GT.BCHI(JJJ+1)) THEN
            X2X    = BCHI(JJJ+1)
            IPHASX = -IPHASE
            JJN    = JJJ + 1
          ELSE
            X2X    = XTOT
            IPHASX = -IPHASE
            JJN    = JJJ
          END IF
        END IF

      ELSE
        IPHASX = IPHASE
        IF (XTOT.GT.BCHI(JJJ+1)) THEN
          X2X = BCHI(JJJ+1)
          JJN = JJJ + 1
        ELSE IF (XTOT.LT.BCHI(JJJ)) THEN
          X2X = BCHI(JJJ)
          JJN = JJJ - 1
        ELSE
          X2X = XTOT
          JJN = JJJ
        END IF
      END IF
      SX2   = SIN(X2X)
      CX2   = COS(X2X)
      SP2   = SP1*SX1/SX2
      IF (SP2.GT.1.0D0) SP2 = 1.0D0
      CP2   = IPHASX * SQRT(1.0D0 - SP2**2)
      ANUM  = CX1*CX2 - SX1*SX2*CP1*CP2
      DNOM  = 1.0D0 - SX1*SX2*SP1*SP2

      RAT=ANUM/DNOM
      IF (RAT>1.0) THEN
        RAT=1.0
      ELSE IF (RAT<-1.0) THEN
        RAT=-1.0
      END IF
      OMEGX = ACOS(RAT)

      IF (OMEGX.LT.OMEGR) THEN
        OMEGA = OMEGX
        THET2 = THETA - OMEGA
        R2    = RSMU/SIN(THET2)
        X2    = X2X
        IIN   = III
      ELSE
        OMEGA = OMEGR
        THET2 = THETR
        R2    = R2R
        X2    = ACOS(CX1*COS(OMEGA)+SX1*SIN(OMEGA)*CP1)
        JJN   = JJJ
      END IF
      IF ((OMEGA.GE.OSKIM).AND.(PASS)) THEN
        IPHASE = -IPHASE
        PASS   = .FALSE.
      END IF
      DELT2 = SQRT(R2**2 + R1**2 - 2.0D0*R2*R1*COS(OMEGA))

C -----------------------------------------------------------------------------
C NOW TO INTEGRATE ALONG LOS INSIDE BIN

      IF ((III.GE.1).AND.(III.LE.IKNT)) THEN

C SINGLY-SCATTERED COMPONENT
        SSOL = 0.0D0

        ICROSS = 0
        IF (ISHAD.EQ.0) THEN
          IF (.NOT.ILLUM) GO TO 2003
          DUPR = DELT2
          DLWR = DELTA
        ELSE IF (ISHAD.EQ.1) THEN
          IF (.NOT.ILLUM) THEN
            IF (DELT2.LE.DSHAD) THEN
              CALL STEP_TAU(ITHERM,THERMO,R1,X1,THETA,DELT2,DELTA,
     &                      LPROF,SPDPT,NONIT,TO2LOC,TLPLOC)
              TAUO2S = TAUO2S + TO2LOC
              DO LP = 0,LPROF
                TAULPS(LP) = TAULPS(LP) + TLPLOC(LP)
              END DO
              GO TO 2003
            END IF
            DUPR = DELT2
            DLWR = MAX(DSHAD,DELTA)
          ELSE
            IF (DELTA.GE.DSHAD) GO TO 2003
            IF (DELT2.LE.DSHAD) THEN
              DUPR = DELT2
            ELSE
              DUPR = DSHAD
              ICROSS = 1
            END IF
            DLWR = DELTA
          END IF
        ELSE IF (ISHAD.EQ.2) THEN
          IF (DELTA.LT.DSHAD) THEN
            IF (DELT2.LE.DSHAD) THEN
              DUPR = DELT2
            ELSE
              DUPR = DSHAD
              ICROSS = 1
              IF (DELT2.GT.DSHAD2) ICROSS = 2
            END IF
            DLWR = DELTA
          ELSE IF (DELT2.LT.DSHAD2) THEN
            CALL STEP_TAU(ITHERM,THERMO,R1,X1,THETA,DELT2,DELTA,
     &                    LPROF,SPDPT,NONIT,TO2LOC,TLPLOC)
            TAUO2S = TAUO2S + TO2LOC
            DO LP = 0,LPROF
              TAULPS(LP) = TAULPS(LP) + TLPLOC(LP)
            END DO
            GO TO 2003
          ELSE
            DUPR = DELT2
            DLWR = MAX(DSHAD2,DELTA)
          END IF
        END IF
        IF (DLWR.GT.DELTA) THEN
          CALL STEP_TAU(ITHERM,THERMO,R1,X1,THETA,DLWR,DELTA,
     &                  LPROF,SPDPT,NONIT,TO2LOC,TLPLOC)
          TAUO2S = TAUO2S + TO2LOC
          DO LP = 0,LPROF
            TAULPS(LP) = TAULPS(LP) + TLPLOC(LP)
          END DO
        END IF

 2002   CONTINUE
        RLWR = SQRT(R1**2 + DLWR**2 - 2.0D0*R1*DLWR*COS(PI-THETA))
        DWAS = DLWR
        CALL DLOCAL(ITHERM,THERMO,RLWR,X1,HWAS,O2WAS,TWAS)
        TO2LOC = 0.0D0
        DO LP = 0,LPROF
          TLPLOC(LP) = 0.0D0
        END DO
        DO 2001 IDEL = 1,8
          DL  = ((DUPR-DLWR)*XX08(IDEL) + DUPR+DLWR) / 2.0D0
          WDL =  (DUPR-DLWR)*WW08(IDEL) / 2.0D0
          RL  = SQRT(R1**2 + DL**2 - 2.0D0*R1*DL*COS(PI-THETA))

          RAT=RSMU/RL
C         RAT outside ASIN limits gives NaN in Radiance
          IF (RAT>1.0) THEN
            RAT=1.0
          ELSE IF (RAT<-1.0) THEN
            RAT=-1.0
          END IF

          IF ((CUMU.LT.0.0D0).AND.(DL.LT.DSMU)) THEN
            OMEGL = THETA - (PI - ASIN(RAT))
          ELSE
            OMEGL = THETA - ASIN(RAT)
          END IF
          XL  = ACOS(CX1*COS(OMEGL)+SX1*SIN(OMEGL)*CP1)
          CALL DLOCAL(ITHERM,THERMO,RL,X1,HLOC,O2LOC,TLOC)

          IF (((RL*SIN(XL)).LE.SOLILLUMCHK).AND.(XL.GE.PID2)) THEN
            SSOL = 0.0D0
          ELSE
            CALL ZSUN(ITHERM,THERMO,RL,XL,OSUN,
     &                LPROF,SPDPT,WSPD,NONIT,EXOLOS)
            IF (EXOLOS) THEN
              SSOL = BRANCH*CENTER*HLOC*EXP(-OSUN(2))*TRANS(OSUN(1))
            ELSE
              SSOL = BRANCH*CENTER*HLOC*EXP(-OSUN(2))*OSUN(3)
            END IF
          END IF

          CALL INLINE_TAU(DWAS,HWAS,O2WAS,TWAS,
     &                    DL  ,HLOC,O2LOC,TLOC,
     &                    LPROF,SPDPT,NONIT,TO2LOC,TLPLOC)
C H-alpha photons travel freely to observer: TABS=1, TSCAT=1
    !       TABS   = EXP(-(TAUO2S + TO2LOC))
    !       FCTREF = TEXO/TLOC
    !       TSCAT  = 0.0D0
    !       DO LP = 1,LPROF
    !         XPNT  = (FCTREF - 1.0D0) * (SPDPT(LP)**2)
    !         TSCAT = TSCAT + WSPD(LP) *
    !  &              EXP(-(TAULPS(LP)+TLPLOC(LP))) * EXP(-XPNT)
    !       END DO
    !       TSCAT  = TSCAT * SQRT(FCTREF) * 2.0D0 / RTPI
          ACER(1,LOOK) = ACER(1,LOOK) + WDL*ANSCAT*SSOL ! *TABS*TSCAT
 2001   CONTINUE
        CALL STEP_TAU(ITHERM,THERMO,R1,X1,THETA,DUPR,DLWR,
     &                LPROF,SPDPT,NONIT,TO2LOC,TLPLOC)
        TAUO2S = TAUO2S + TO2LOC
        DO LP = 0,LPROF
          TAULPS(LP) = TAULPS(LP) + TLPLOC(LP)
        END DO
        IF (ICROSS.EQ.1) THEN
          ICROSS = 0
          CALL STEP_TAU(ITHERM,THERMO,R1,X1,THETA,DELT2,DSHAD,
     &                  LPROF,SPDPT,NONIT,TO2LOC,TLPLOC)
          TAUO2S = TAUO2S + TO2LOC
          DO LP = 0,LPROF
            TAULPS(LP) = TAULPS(LP) + TLPLOC(LP)
          END DO
        ELSE IF (ICROSS.EQ.2) THEN
          ICROSS = 0
          CALL STEP_TAU(ITHERM,THERMO,R1,X1,THETA,DSHAD2,DSHAD,
     &                  LPROF,SPDPT,NONIT,TO2LOC,TLPLOC)
          TAUO2S = TAUO2S + TO2LOC
          DO LP = 0,LPROF
            TAULPS(LP) = TAULPS(LP) + TLPLOC(LP)
          END DO
          DUPR = DELT2
          DLWR = DSHAD2
          GO TO 2002
        END IF

 2003   CONTINUE
C MULTIPLY-SCATTERED COMPONENT

        RLWR = SQRT(R1**2 + DELTA**2 - 2.0D0*R1*DELTA*COS(PI-THETA))
        DWAS = DELTA
        CALL DLOCAL(ITHERM,THERMO,RLWR,X1,HWAS,O2WAS,TWAS)
        TO2LOC = 0.0D0
        DO LP = 0,LPROF
          TLPLOC(LP) = 0.0D0
        END DO
        DO 2004 IDEL = 1,8
          DL  = ((DELT2-DELTA)*XX08(IDEL) + DELT2+DELTA) / 2.0D0
          WDL =  (DELT2-DELTA)*WW08(IDEL) / 2.0D0
          RL  = SQRT(R1**2 + DL**2 - 2.0D0*R1*DL*COS(PI-THETA))

          RAT=RSMU/RL
C         RAT outside ASIN limits gives NaN in Radiance
          IF (RAT>1.0) THEN
            RAT=1.0
          ELSE IF (RAT<-1.0) THEN
            RAT=-1.0
          END IF

          IF ((CUMU.LT.0.0D0).AND.(DL.LT.DSMU)) THEN
            OMEGL = THETA - (PI - ASIN(RAT))
          ELSE
            OMEGL = THETA - ASIN(RAT)
          END IF
          XL  = ACOS(CX1*COS(OMEGL)+SX1*SIN(OMEGL)*CP1)
          CALL DLOCAL(ITHERM,THERMO,RL,X1,HLOC,O2LOC,TLOC)

          IF ((XL.GT.CHIN(1)).AND.(XL.LT.CHIN(JKNT))) THEN
            IF (XL.LE.CHIN(JJJ)) THEN
              JBOX1 = JJJ - 1
              JBOX2 = JJJ
            ELSE
              JBOX1 = JJJ
              JBOX2 = JJJ + 1
            END IF
            DO ISPLN = 1,IKNT+2
              YSMS(ISPLN)   = D0SMS(ISPLN,JBOX1)
              D2YSMS(ISPLN) = D2SMS(ISPLN,JBOX1)
            END DO
            CALL SPLINT_LYAO(XSMS,YSMS,D2YSMS,IKNT+2,RL,SMS1)
            DO ISPLN = 1,IKNT+2
              YSMS(ISPLN)   = D0SMS(ISPLN,JBOX2)
              D2YSMS(ISPLN) = D2SMS(ISPLN,JBOX2)
            END DO
            CALL SPLINT_LYAO(XSMS,YSMS,D2YSMS,IKNT+2,RL,SMS2)
            RATIOX = (XL - CHIN(JBOX1))/(CHIN(JBOX2) - CHIN(JBOX1))
            SMS    = BRANCH*CENTER*HLOC*EXP(SMS1 + RATIOX*(SMS2-SMS1))
          ELSE
            DO ISPLN = 1,IKNT+2
              YSMS(ISPLN)   = D0SMS(ISPLN,JJJ)
              D2YSMS(ISPLN) = D2SMS(ISPLN,JJJ)
            END DO
            CALL SPLINT_LYAO(XSMS,YSMS,D2YSMS,IKNT+2,RL,SMS)
            SMS    = BRANCH*CENTER*HLOC*EXP(SMS)
          END IF

          CALL INLINE_TAU(DWAS,HWAS,O2WAS,TWAS,
     &                    DL  ,HLOC,O2LOC,TLOC,
     &                    LPROF,SPDPT,NONIT,TO2LOC,TLPLOC)
C H-alpha photons travel freely to observer: TABS=1, TSCAT=1
    !       TABS   = EXP(-(TAUO2M + TO2LOC))
    !       FCTREF = TEXO/TLOC
    !       TSCAT  = 0.0D0
    !       DO LP = 1,LPROF
    !         XPNT  = (FCTREF - 1.0D0) * (SPDPT(LP)**2)
    !         TSCAT = TSCAT + WSPD(LP) *
    !  &              EXP(-(TAULPM(LP)+TLPLOC(LP))) * EXP(-XPNT)
    !       END DO
    !       TSCAT  = TSCAT * SQRT(FCTREF) * 2.0D0 / RTPI
          ACER(2,LOOK) = ACER(2,LOOK) + WDL*SMS ! *TABS*TSCAT

 2004   CONTINUE
        CALL STEP_TAU(ITHERM,THERMO,R1,X1,THETA,DELT2,DELTA,
     &                LPROF,SPDPT,NONIT,TO2LOC,TLPLOC)
        TAUO2M = TAUO2M + TO2LOC
        DO LP = 0,LPROF
          TAULPM(LP) = TAULPM(LP) + TLPLOC(LP)
        END DO

C TOTAL INTENSITY
        ACER(3,LOOK) = ACER(1,LOOK) + ACER(2,LOOK)

      END IF
C -----------------------------------------------------------------------------
C     exiting the volume
      IF ((IIN.LT.1).OR.(IIN.GT.IKNT)) GO TO 2222
C     effective end-of-the-line
      IF (TAUO2M.GT.TLIMIT) GO TO 2222
      DELTA = DELT2
      III   = IIN
      JJJ   = JJN
      GO TO 1111

 2222 CONTINUE

C + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +

 3333 CONTINUE
!$OMP END PARALLEL DO

C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

C OUTPUT OF APPARENT COLUMN EMISSION RATES (ACERs)

C CONVERTING TO ABSOLUTE APPARENT COLUMN EMISSION RATES (RAYLEIGHS)
C keep in mind that the line-center solar flux is in per-Angstrom units,
C so that the wavelength doppler width must be in Angstroms.

      DPPLR  = WAVELN * (VELT / SPEEDC)
      ABFCTR = RTPI * DPPLR / 1.0D6
C      ABFCTR = RTPI * DPPLR
      DO LOOK = 1,ILOS
        ACER(1,LOOK) = ABFCTR * ACER(1,LOOK)
        ACER(2,LOOK) = ABFCTR * ACER(2,LOOK)
        ACER(3,LOOK) = ABFCTR * ACER(3,LOOK)
      END DO

      RETURN
      END
C =============================================================================
