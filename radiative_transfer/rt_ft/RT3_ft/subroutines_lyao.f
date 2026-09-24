C 1857
C=============================================================================
C LIBRARY OF SUBROUTINES NEEDED BY LYMAN LINE RADIANCE MODELING CODES
C
C     rt     los     SUBROUTINE
C
C     x       x      ZSUN
C     x              ZONE
C     x       x      TRANS
C             x      SHADOW
C             x      DLOCAL
C             x      STEP_TAU
C             x      INLINE_TAU
C             x      SMS_SPLINE
C     x              matrix inversion routines   (LUDCMP+LUBKSB+MPROVE)
C     x       x      gaussian quadrature points
C     x       x      "speed" quadrature points
C
C =============================================================================
C =============================================================================

C POINT-TO-SUN LINE-OF-SIGHT PROPAGATION:
C FOR DIRECT SOLAR SOURCE FUNCTION EVALUATIONS
C . . . OSUN(1):  OPTICAL DEPTH AT LINE CENTER FOR SCATTERING
C . . . OSUN(2):  OPTICAL DEPTH FOR ABSORPTION
C NOTE THAT IN THIS APPLICATION, THETA IS IDENTICALLY EQUAL TO CHI AND
C PHI IS IDENTICALLY ZERO, SO (IN THE NOTATION OF ZONE) THET2=X2;
C HENCE THETA & THET2 ARE ELIMINATED IN FAVOR OF X1 & X2.

      SUBROUTINE ZSUN(ITHERM,THERMO,R1,X1,OSUN,
     &                LPROF,SPDPT,WSPD,NONIT,EXOLOS)
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL NONIT,EXOLOS
      DIMENSION THERMO(5,ITHERM+1),OSUN(3)
      DIMENSION SPDPT(LPROF),WSPD(LPROF),TAULP(0 : 16)
      COMMON /EXOS/  GM,CONL,PLANETR,RBASE,RC,RUPR,RP,TEXO,DEXO,VELT,
     &               RADPF,CENTER,ABSCSX,FTSAT,FDSAT,IGEO
      COMMON /NMBR/  PI,RTPI,PID2,OFFSET
      COMMON /GAUS16/WW16(16),XX16(16)
C      COMMON/FINE_GRID/ THERMO_MSIS,BASE,TOP,ITHERM1

C        PRINT*,'--------'
C        PRINT*,GM,CONL,PLANETR,RBASE,RC,RUPR,RP,RADPF
C        PRINT*,TEXO,DEXO,VELT,CENTER
C        PRINT*,ABSCSX,FTSAT,FDSAT,IGEO
C        PRINT*,THERMO_MSIS
C        PRINT*,BASE,TOP,ITHERM
C        PRINT*,PI,RTPI,PID2,OFFSET

C BASIC REFERENCE QUANTITIES
      RSMU   = R1 * SIN(X1)
      EXOLOS = .TRUE.

C -----------------------------------------------------------------------------
C LOS DOES NOT INTERSECT THERMOSPHERE: DAYSIDE EXOSPHERE (ISOTHERMAL)

      IF ((X1.LE.PID2).AND.(R1.GE.RC)) THEN
        XUPR = X1
        XLWR = ASIN(RSMU/RUPR)
        TAULP(0) = 0.0D0
        DO 101 ISUN = 1,16
          CHI  = ((XUPR-XLWR)*XX16(ISUN) + XUPR+XLWR) / 2.0D0
          WCHI = (XUPR-XLWR) *WW16(ISUN) / 2.0D0
          RR   = RSMU / SIN(CHI)
          CALL CORONA(RR,DNST)
          TAULP(0) = TAULP(0) + WCHI * DNST * RR / SIN(CHI)
 101    CONTINUE
        OSUN(1) = CENTER * TAULP(0)
        OSUN(2) = 0.0D0

      ELSE IF (RSMU.GE.RC) THEN
        XUPR = X1
        XLWR = PID2
        TAULP(0) = 0.0D0
        DO 201 ISUN = 1,16
          CHI  = ((XUPR-XLWR)*XX16(ISUN) + XUPR+XLWR) / 2.0D0
          WCHI = (XUPR-XLWR) *WW16(ISUN) / 2.0D0
          RR   = RSMU / SIN(CHI)
          CALL CORONA(RR,DNST)
          TAULP(0) = TAULP(0) + WCHI * DNST * RR / SIN(CHI)
 201    CONTINUE
        TAULP(0) = 2.0D0 * TAULP(0)
        XUPR = PI - X1
        XLWR = ASIN(RSMU/RUPR)
        DO 202 ISUN = 1,16
          CHI  = ((XUPR-XLWR)*XX16(ISUN) + XUPR+XLWR) / 2.0D0
          WCHI = (XUPR-XLWR) *WW16(ISUN) / 2.0D0
          RR   = RSMU / SIN(CHI)
          CALL CORONA(RR,DNST)
          TAULP(0) = TAULP(0) + WCHI * DNST * RR / SIN(CHI)
 202    CONTINUE
        OSUN(1) = CENTER * TAULP(0)
        OSUN(2) = 0.0D0

C -----------------------------------------------------------------------------

C THERMOSPHERE CONTRIBUTION
      ELSE IF (R1.LT.RC) THEN
        TAULP(0) = 0.0D0
        IF (NONIT) THEN
          DO LP = 1,LPROF
            TAULP(LP) = 0.0D0
          END DO
          EXOLOS = .FALSE.
        END IF
        IHI = ITHERM + 1
        ILO = 1
 301    IF ((IHI-ILO).GT.1) THEN
          II = (IHI + ILO) / 2
          IF (THERMO(2,II).GT.R1) THEN 
            IHI = II
          ELSE
            ILO = II
          END IF
          GO TO 301
        END IF
        II = ILO
C . . . ESTIMATE LOCAL DENSITIES
        DZ  = THERMO(2,II+1) - THERMO(2,II)
        DR  = (R1 - THERMO(2,II)) / DZ
        XPT = THERMO(4,II) + DR * (THERMO(4,II+1)-THERMO(4,II))
        H1N = EXP( XPT )
        XPT = THERMO(5,II) + DR * (THERMO(5,II+1)-THERMO(5,II))
        O2N = EXP( XPT )
        TPL = THERMO(3,II) + DR * (THERMO(3,II+1)-THERMO(3,II))
C . . . SEARCH FOR EDGE OF INITIAL BIN ALONG SPECIFIED LINE OF SIGHT
        ITRIP = 0
        IF (COS(X1).GT.0.0D0) THEN
          R2  = THERMO(2,II+1)
          X2  = ASIN(RSMU/R2)
          III = II + 1
          H1Z = EXP( THERMO(4,III) )
          O2Z = EXP( THERMO(5,III) )
          TPZ = THERMO(3,III)

          XUPR = X1
          XLWR = X2
          HC0 = (R1-R2)/LOG(H1Z/H1N)
          OC0 = (R1-R2)/LOG(O2Z/O2N)
          TAULP(0) = 0.0D0
          TAUO2 = 0.0D0
          DO 169 ISUN = 1,16
            CHI  = ((XUPR-XLWR)*XX16(ISUN) + XUPR+XLWR) / 2.0D0
            WCHI = (XUPR-XLWR) *WW16(ISUN) / 2.0D0
            RR   = RSMU / SIN(CHI)
            DNST = H1N*EXP(-(RR-R1)/HC0)
            ONST = O2N*EXP(-(RR-R1)/OC0)
            TAULP(0) = TAULP(0) + WCHI * DNST * RR / SIN(CHI)
            TAUO2 = TAUO2 + WCHI * ONST * RR / SIN(CHI)
 169      CONTINUE
          TAULP(0) = CENTER * TAULP(0)
          TAUO2 = ABSCSX * TAUO2

        ELSE IF (RSMU.GE.THERMO(2,II)) THEN

          R2  = THERMO(2,II+1)
          X2  = ASIN(RSMU/R2)
          III = II + 1
          DR  = (RSMU - THERMO(2,II)) / DZ
          XPT = THERMO(4,II) + DR * (THERMO(4,III)-THERMO(4,II))
          H1N = EXP( XPT )
          XPT = THERMO(5,II) + DR * (THERMO(5,III)-THERMO(5,II))
          O2N = EXP( XPT )
          TPL = THERMO(3,II) + DR * (THERMO(3,III)-THERMO(3,II))
          H1Z = EXP( THERMO(4,III) )
          O2Z = EXP( THERMO(5,III) )
          TPZ = THERMO(3,III)

          XUPR = X1
          XLWR = X2
          HC0 = (RSMU-R2)/LOG(H1Z/H1N)
          OC0 = (RSMU-R2)/LOG(O2Z/O2N)
          TAULP(0) = 0.0D0
          TAUO2 = 0.0D0
          DO 170 ISUN = 1,16
            CHI  = ((XUPR-XLWR)*XX16(ISUN) + XUPR+XLWR) / 2.0D0
            WCHI = (XUPR-XLWR) *WW16(ISUN) / 2.0D0
            RR   = RSMU / SIN(CHI)
            DNST = H1N*EXP(-(RR-RSMU)/HC0)
            ONST = O2N*EXP(-(RR-RSMU)/OC0)
            TAULP(0) = TAULP(0) + WCHI * DNST * RR / SIN(CHI)
            TAUO2 = TAUO2 + WCHI * ONST * RR / SIN(CHI)
 170      CONTINUE
          TAULP(0) = CENTER * TAULP(0)
          TAUO2 = ABSCSX * TAUO2

        ELSE

          R2  = THERMO(2,II)
          X2  = PI - ASIN(RSMU/R2)
          III = II - 1
          H1Z = EXP( THERMO(4,II) )
          O2Z = EXP( THERMO(5,II) )
          TPZ = THERMO(3,II)

          XUPR = X1
          XLWR = X2
          HC0 = (R1-R2)/LOG(H1Z/H1N)
          OC0 = (R1-R2)/LOG(O2Z/O2N)
          TAULP(0) = 0.0D0
          TAUO2 = 0.0D0
          DO 171 ISUN = 1,16
            CHI  = ((XUPR-XLWR)*XX16(ISUN) + XUPR+XLWR) / 2.0D0
            WCHI = (XUPR-XLWR) *WW16(ISUN) / 2.0D0
            RR   = RSMU / SIN(CHI)
            DNST = H1N*EXP(-(RR-R1)/HC0)
            ONST = O2N*EXP(-(RR-R1)/OC0)
            TAULP(0) = TAULP(0) + WCHI * DNST * RR / SIN(CHI)
            TAUO2 = TAUO2 + WCHI * ONST * RR / SIN(CHI)
 171      CONTINUE
          TAULP(0) = CENTER * TAULP(0)
          TAUO2 = ABSCSX * TAUO2

        END IF

        OSUN(1)  = TAULP(0)
        OSUN(2)  = TAUO2

        IF (NONIT) THEN
          TLOC = (TPL + TPZ) / 2.0D0
          TREF = TLOC
          FCTR = SQRT(TEXO/TLOC)
          DO LP = 1,LPROF
            XPNT = (SPDPT(LP)*FCTR)**2
            TAULP(LP) = TAULP(0) * FCTR * EXP(-XPNT)
          END DO
        END IF
        IF (III.GT.ITHERM) GO TO 304

 303    CONTINUE
        ITRIP = ITRIP + 1
        H1N = H1Z
        O2N = O2Z
        TPL = TPZ
        IF (COS(X2).GT.0.0D0) THEN

          RPRE = R2
          XPRE = X2
          HPRE = H1Z
          OPRE = O2Z

          R2  = THERMO(2,III+1)
          X2  = ASIN(RSMU/R2)
          IIN = III + 1
          H1Z = EXP( THERMO(4,IIN) )
          O2Z = EXP( THERMO(5,IIN) )
          TPZ = THERMO(3,IIN)

          XUPR = XPRE
          XLWR = X2
          HC0 = (RPRE-R2)/LOG(H1Z/HPRE)
          OC0 = (RPRE-R2)/LOG(O2Z/OPRE)
          TAULP(0) = 0.0D0
          TAUO2 = 0.0D0
          DO 172 ISUN = 1,16
            CHI  = ((XUPR-XLWR)*XX16(ISUN) + XUPR+XLWR) / 2.0D0
            WCHI = (XUPR-XLWR) *WW16(ISUN) / 2.0D0
            RR   = RSMU / SIN(CHI)
            DNST = HPRE*EXP(-(RR-RPRE)/HC0)
            ONST = OPRE*EXP(-(RR-RPRE)/OC0)
            TAULP(0) = TAULP(0) + WCHI * DNST * RR / SIN(CHI)
            TAUO2 = TAUO2 + WCHI * ONST * RR / SIN(CHI)
 172      CONTINUE
          TAULP(0) = CENTER * TAULP(0)
          TAUO2 = ABSCSX * TAUO2

        ELSE IF (RSMU.GE.THERMO(2,III)) THEN

          XPRE = X2

          R2  = THERMO(2,III+1)
          X2  = ASIN(RSMU/R2)
          IIN = III + 1
          DZ  = THERMO(2,IIN) - THERMO(2,III)
          DR  = (RSMU - THERMO(2,III)) / DZ
          XPT = THERMO(4,III) + DR * (THERMO(4,IIN)-THERMO(4,III))
          H1N = EXP( XPT )
          XPT = THERMO(5,III) + DR * (THERMO(5,IIN)-THERMO(5,III))
          O2N = EXP( XPT )
          TPL = THERMO(3,III) + DR * (THERMO(3,IIN)-THERMO(3,III))
          H1Z = EXP( THERMO(4,IIN) )
          O2Z = EXP( THERMO(5,IIN) )
          TPZ = THERMO(3,IIN)

          XUPR = XPRE
          XLWR = X2
          HC0 = (RSMU-R2)/LOG(H1Z/H1N)
          OC0 = (RSMU-R2)/LOG(O2Z/O2N)
          TAULP(0) = 0.0D0
          TAUO2 = 0.0D0
          DO 173 ISUN = 1,16
            CHI  = ((XUPR-XLWR)*XX16(ISUN) + XUPR+XLWR) / 2.0D0
            WCHI = (XUPR-XLWR) *WW16(ISUN) / 2.0D0
            RR   = RSMU / SIN(CHI)
            DNST = H1N*EXP(-(RR-RSMU)/HC0)
            ONST = O2N*EXP(-(RR-RSMU)/OC0)
            TAULP(0) = TAULP(0) + WCHI * DNST * RR / SIN(CHI)
            TAUO2 = TAUO2 + WCHI * ONST * RR / SIN(CHI)
 173      CONTINUE
          TAULP(0) = CENTER * TAULP(0)
          TAUO2 = ABSCSX * TAUO2

        ELSE

          RPRE = R2
          XPRE = X2
          HPRE = H1Z
          OPRE = O2Z

          R2  = THERMO(2,III)
          X2  = PI - ASIN(RSMU/R2)
          IIN = III - 1
          H1Z = EXP( THERMO(4,III) )
          O2Z = EXP( THERMO(5,III) )
          TPZ = THERMO(3,III)

          XUPR = XPRE
          XLWR = X2
          HC0 = (RPRE-R2)/LOG(H1Z/HPRE)
          OC0 = (RPRE-R2)/LOG(O2Z/OPRE)
          TAULP(0) = 0.0D0
          TAUO2 = 0.0D0
          DO 174 ISUN = 1,16
            CHI  = ((XUPR-XLWR)*XX16(ISUN) + XUPR+XLWR) / 2.0D0
            WCHI = (XUPR-XLWR) *WW16(ISUN) / 2.0D0
            RR   = RSMU / SIN(CHI)
            DNST = HPRE*EXP(-(RR-RPRE)/HC0)
            ONST = OPRE*EXP(-(RR-RPRE)/OC0)
            TAULP(0) = TAULP(0) + WCHI * DNST * RR / SIN(CHI)
            TAUO2 = TAUO2 + WCHI * ONST * RR / SIN(CHI)
 174      CONTINUE
          TAULP(0) = CENTER * TAULP(0)
          TAUO2 = ABSCSX * TAUO2

        END IF

        OSUN(1)  = OSUN(1) + TAULP(0)
        OSUN(2)  = OSUN(2) + TAUO2
        IF (NONIT) THEN
          TLOC = (TPL + TPZ) / 2.0D0
          FCTR = SQRT(TEXO/TLOC)
          DO LP = 1,LPROF
            XPNT = (SPDPT(LP)*FCTR)**2
            TAULP(LP) = TAULP(LP) + TAULP(0) * FCTR * EXP(-XPNT)
          END DO
        END IF
        IF (IIN.GT.ITHERM) GO TO 304
        III   = IIN
        GO TO 303
 304    CONTINUE

        XUPR = X2
        XLWR = ASIN(RSMU/RUPR)
        TAULP(0) = 0.0D0
        DO 305 ISUN = 1,16
          CHI  = ((XUPR-XLWR)*XX16(ISUN) + XUPR+XLWR) / 2.0D0
          WCHI = (XUPR-XLWR) *WW16(ISUN) / 2.0D0
          RR   = RSMU / SIN(CHI)
          CALL CORONA(RR,DNST)
          TAULP(0) = TAULP(0) + WCHI * CENTER * DNST * RR / SIN(CHI)
 305    CONTINUE
        OSUN(1) = OSUN(1) + TAULP(0)
        IF (NONIT) THEN
          DO LP = 1,LPROF
            XPNT = SPDPT(LP)**2
            TAULP(LP) = TAULP(LP) + TAULP(0) * EXP(-XPNT)
          END DO
        END IF

C -----------------------------------------------------------------------------
C LOS PASSES THROUGH THERMOSPHERE:  EVALUATE THERMOSPHERIC
C CONTRIBUTION BY PROPAGATING FROM TERMINATOR TO DAYSIDE EXOBASE,
C XX REFERS TO THE TERMINATOR (I.E., MTRH).

      ELSE
        TREF = TEXO
        TAULP(0) = 0.0D0
        IF (NONIT) THEN
          DO LP = 1,LPROF
            TAULP(LP) = 0.0D0
          END DO
          EXOLOS = .FALSE.
        END IF
C . . . FIND RADIAL BIN CONTAINING MTRH
        XX = PID2
        IHI = ITHERM + 1
        ILO = 1
 401    IF ((IHI-ILO).GT.1) THEN
          II = (IHI + ILO) / 2
          IF (THERMO(2,II).GT.RSMU) THEN 
            IHI = II
          ELSE
            ILO = II
          END IF
          GO TO 401
        END IF
        II = ILO
C . . . ESTIMATE LOCAL DENSITIES
        DZ  = THERMO(2,II+1) - THERMO(2,II)
        DR  = (RSMU - THERMO(2,II)) / DZ
        XPT = THERMO(4,II) + DR * (THERMO(4,II+1)-THERMO(4,II))
        H1N = EXP( XPT )
        XPT = THERMO(5,II) + DR * (THERMO(5,II+1)-THERMO(5,II))
        O2N = EXP( XPT )
        TPL = THERMO(3,II) + DR * (THERMO(3,II+1)-THERMO(3,II))
        ITRIP = 0
        R2  = THERMO(2,II+1)
        III = II + 1
        H1Z = EXP( THERMO(4,III) )
        O2Z = EXP( THERMO(5,III) )
        TPZ = THERMO(3,III)

        XUPR = XX
        XLWR = ASIN(RSMU/R2)
        HC0 = (RSMU-R2)/LOG(H1Z/H1N)
        OC0 = (RSMU-R2)/LOG(O2Z/O2N)
        TAULP(0) = 0.0D0
        TAUO2 = 0.0D0
        DO 175 ISUN = 1,16
          CHI  = ((XUPR-XLWR)*XX16(ISUN) + XUPR+XLWR) / 2.0D0
          WCHI = (XUPR-XLWR) *WW16(ISUN) / 2.0D0
          RR   = RSMU / SIN(CHI)
          DNST = H1N*EXP(-(RR-RSMU)/HC0)
          ONST = O2N*EXP(-(RR-RSMU)/OC0)
          TAULP(0) = TAULP(0) + WCHI * DNST * RR / SIN(CHI)
          TAUO2 = TAUO2 + WCHI * ONST * RR / SIN(CHI)
 175    CONTINUE
        TAULP(0) = CENTER * TAULP(0)
        TAUO2 = ABSCSX * TAUO2

        OSUN(1)  = TAULP(0)
        OSUN(2)  = TAUO2
        IF (NONIT) THEN
          TLOC = (TPL + TPZ) / 2.0D0
          FCTR = SQRT(TEXO/TLOC)
          DO LP = 1,LPROF
            XPNT = (SPDPT(LP)*FCTR)**2
            TAULP(LP) = TAULP(0) * FCTR * EXP(-XPNT)
          END DO
        END IF
        IF (III.GT.ITHERM) GO TO 404
 403    CONTINUE

        XPRE = ASIN(RSMU/R2)
        RPRE = R2

        ITRIP = ITRIP + 1
        H1N = H1Z
        O2N = O2Z
        TPL = TPZ
        R2  = THERMO(2,III+1)
        IIN = III + 1
        H1Z = EXP( THERMO(4,IIN) )
        O2Z = EXP( THERMO(5,IIN) )
        TPZ = THERMO(3,IIN)

        XUPR = XPRE
        XLWR = ASIN(RSMU/R2)
        HC0 = (RPER-R2)/LOG(H1Z/H1N)
        OC0 = (RPRE-R2)/LOG(O2Z/O2N)
        TAULP(0) = 0.0D0
        TAUO2 = 0.0D0
        DO 176 ISUN = 1,16
          CHI  = ((XUPR-XLWR)*XX16(ISUN) + XUPR+XLWR) / 2.0D0
          WCHI = (XUPR-XLWR) *WW16(ISUN) / 2.0D0
          RR   = RSMU / SIN(CHI)
          DNST = H1N*EXP(-(RR-RPRE)/HC0)
          ONST = O2N*EXP(-(RR-RPRE)/OC0)
          TAULP(0) = TAULP(0) + WCHI * DNST * RR / SIN(CHI)
          TAUO2 = TAUO2 + WCHI * ONST * RR / SIN(CHI)
 176    CONTINUE
        TAULP(0) = CENTER * TAULP(0)
        TAUO2 = ABSCSX * TAUO2

        OSUN(1)  = OSUN(1) + TAULP(0)
        OSUN(2)  = OSUN(2) + TAUO2
        IF (NONIT) THEN
          TLOC = (TPL + TPZ) / 2.0D0
          FCTR = SQRT(TEXO/TLOC)
          DO LP = 1,LPROF
            XPNT = (SPDPT(LP)*FCTR)**2
            TAULP(LP) = TAULP(LP) + TAULP(0) * FCTR * EXP(-XPNT)
          END DO
        END IF
        IF (IIN.GT.ITHERM) GO TO 404
        III   = IIN
        GO TO 403
 404    CONTINUE
        OSUN(1) = 2.0D0 * OSUN(1)
        OSUN(2) = 2.0D0 * OSUN(2)
        IF (NONIT) THEN
          DO LP = 1,LPROF
            TAULP(LP) = 2.0D0 * TAULP(LP)
          END DO
        END IF

        XUPR = X1
        XLWR = PI - ASIN(RSMU/RC)
        TAULP(0) = 0.0D0
        DO 405 ISUN = 1,16
          CHI  = ((XUPR-XLWR)*XX16(ISUN) + XUPR+XLWR) / 2.0D0
          WCHI = (XUPR-XLWR) *WW16(ISUN) / 2.0D0
          RR   = RSMU / SIN(CHI)
          CALL CORONA(RR,DNST)
          TAULP(0) = TAULP(0) + WCHI * CENTER * DNST * RR / SIN(CHI)
 405    CONTINUE

        XUPR = ASIN(RSMU/RC)
        XLWR = ASIN(RSMU/RUPR)
        DO 406 ISUN = 1,16
          CHI  = ((XUPR-XLWR)*XX16(ISUN) + XUPR+XLWR) / 2.0D0
          WCHI = (XUPR-XLWR) *WW16(ISUN) / 2.0D0
          RR   = RSMU / SIN(CHI)
          CALL CORONA(RR,DNST)
          TAULP(0) = TAULP(0) + WCHI * CENTER * DNST * RR / SIN(CHI)
 406    CONTINUE
        OSUN(1) = OSUN(1) + TAULP(0)
        IF (NONIT) THEN
          DO LP = 1,LPROF
            XPNT = SPDPT(LP)**2
            TAULP(LP) = TAULP(LP) + TAULP(0) * EXP(-XPNT)
          END DO
        END IF

      END IF

C -----------------------------------------------------------------------------
C FINAL OPTICAL DEPTHS

      OSUN(3) = 0.0D0
      IF (.NOT.EXOLOS) THEN
        FCTREF = TEXO/TREF
        DO 501 LP = 1,LPROF
          XPNT = (FCTREF - 1.0D0) * (SPDPT(LP)**2)
          OSUN(3) = OSUN(3) + WSPD(LP) * EXP(-TAULP(LP)) * EXP(-XPNT)
 501    CONTINUE
        OSUN(3) = OSUN(3) * SQRT(FCTREF) * 2.0D0 / RTPI
      END IF

      RETURN
      END

C =============================================================================
C =============================================================================

C ZONE-TO-ZONE PROPAGATION
C FOR HOMOGENEOUS MATRIX ELEMENT EVALUATIONS.

      SUBROUTINE ZONE(IKNT,JKNT,BRAD,BCHI,HDNS,O2DNS,TEMP,
     &                II,JJ,R1,X1,THETA,PHI,LPROF,SPDPT,NONIT,
     &                TRIP,ODLC,ODO2,IEXIT,JEXIT,ITHERM,THERMO)
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL   NONIT,TRIP(IKNT,JKNT),PASS
      DIMENSION THERMO(5,ITHERM+1)
      DIMENSION BRAD(IKNT+1),BCHI(JKNT+1),HDNS(IKNT),O2DNS(IKNT),
     &          TEMP(IKNT),SPDPT(LPROF) 
      DIMENSION ODLC(2,0 : LPROF,IKNT,JKNT,2),ODO2(2,IKNT,JKNT,2)
      COMMON /EXOS/  GM,CONL,PLANETR,RBASE,RC,RUPR,RP,TEXO,DEXO,VELT,
     &               RADPF,CENTER,ABSCSX,FTSAT,FDSAT,IGEO
      COMMON /NMBR/  PI,RTPI,PID2,OFFSET
      COMMON /GAUS16/WW16(16),XX16(16)

C INITIALIZING STORAGE ARRAYS, MOST ENTRIES OF WHICH WILL BE ZERO/FALSE
      DO 999 III   = 1,IKNT
      DO 999 JJJ   = 1,JKNT
        TRIP(III,JJJ) = .FALSE.
      DO 999 IN_OUT = 1,2
      DO 999 KROSS = 1,2
        ODO2(IN_OUT,III,JJJ,KROSS) = 0.0D0
      DO 999 LP    = 0,LPROF
        ODLC(IN_OUT,LP,III,JJJ,KROSS) = 0.0D0
 999  CONTINUE

C LINE OF SIGHT REFERENCE QUANTITIES
C   RSMU  . . . MINIMUM TANGENT RAY HEIGHT
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
      OSKIM = ACOS(ANUM/DNOM)

      OTOT  = THETA
      IF (OSKIM.LT.OTOT) THEN
        PASS = .TRUE.
      ELSE
        PASS = .FALSE.
      END IF

C -----------------------------------------------------------------------------
C SEARCH FOR EDGE OF INITIAL BIN ALONG SPECIFIED LINE OF SIGHT

      ITRIP = 0
C . . . RADIAL BIN DETERMINATION
      IF ((CUMU.GT.0.0D0).OR.(RSMU.GT.BRAD(II))) THEN
        R2R   = BRAD(II+1)
        THETR = ASIN(RSMU/R2R)
        III   = II + 1
      ELSE
        R2R   = BRAD(II)
        THETR = PI - ASIN(RSMU/R2R)
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
      CP2   = IPHASX * SQRT(1.0D0 - SP2**2)
      ANUM  = CX1*CX2 - SX1*SX2*CP1*CP2
      DNOM  = 1.0D0 - SX1*SX2*SP1*SP2
      OMEGX = ACOS(ANUM/DNOM)

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

      DELTA  = SQRT(R2**2 + R1**2 - 2.0D0*R2*R1*COS(OMEGA))

      COMGA = (R1**2+DELTA**2-R2**2)/(2.0D0*R1*DELTA)
      H_TAU = 0.0D0
      O2_TAU = 0.0D0
      DO 169 ISUN = 1,16
        HCHI  = (XX16(ISUN) + 1.0D0) * DELTA / 2.0D0
        WCHI = DELTA * WW16(ISUN) / 2.0D0
        RR   = SQRT(R1**2+HCHI**2-2.0D0*R1*HCHI*COMGA)

        IF (RR .GT. RC) THEN
          CALL CORONA(RR,DNST)
          ONST = 0.0D0
        ELSE
          IHI = ITHERM + 1
          ILO = 1
 170      IF ((IHI-ILO).GT.1) THEN
            IKK = (IHI + ILO) / 2
            IF (THERMO(2,IKK).GT.RR) THEN
              IHI = IKK
            ELSE
              ILO = IKK
            END IF
            GO TO 170
          END IF
          IKK = ILO

          DZ  = THERMO(2,IKK+1) - THERMO(2,IKK)
          DR  = (RR - THERMO(2,IKK)) / DZ
          XPT = THERMO(4,IKK) + DR * (THERMO(4,IKK+1)-THERMO(4,IKK))
          DNST = EXP( XPT )
          XPT = THERMO(5,IKK) + DR * (THERMO(5,IKK+1)-THERMO(5,IKK))
          ONST = EXP( XPT )
        END IF
        H_TAU = H_TAU + WCHI * DNST
        O2_TAU = O2_TAU + WCHI * ONST
 169  CONTINUE

      H_TAU = CENTER * H_TAU
      O2_TAU = ABSCSX * O2_TAU

      FCTR   = SQRT(TEXO/TEMP(II))
      ODO2(2,II,JJ,1)   = O2_TAU
      ODLC(2,0,II,JJ,1) = H_TAU
      IF (NONIT) THEN
        DO LP = 1,LPROF
          XPNT = (FCTR * SPDPT(LP))**2
          ODLC(2,LP,II,JJ,1) = H_TAU * FCTR * EXP(-XPNT)
        END DO
      END IF
      TRIP(II,JJ) = .TRUE.
      IF ((III.LT.1).OR.(III.GT.IKNT)) GO TO 222
      ODO2(1,III,JJJ,1)   = ODO2(2,II,JJ,1)
      ODLC(1,0,III,JJJ,1) = ODLC(2,0,II,JJ,1)
      IF (NONIT) THEN
        DO LP = 1,LPROF
          ODLC(1,LP,III,JJJ,1) = ODLC(2,LP,II,JJ,1)
        END DO
      END IF

C -----------------------------------------------------------------------------

 111  CONTINUE

      ITRIP = ITRIP + 1
C . . . RADIAL BIN DETERMINATION
      IF ((COS(THET2).GT.0.0D0).OR.(RSMU.GT.BRAD(III))) THEN
        R2R   = BRAD(III+1)
        THETR = ASIN(RSMU/R2R)
        IIN   = III + 1
      ELSE
        R2R   = BRAD(III)
        THETR = PI - ASIN(RSMU/R2R)
        IIN   = III - 1
      END IF
      OMEGR = THETA - THETR

C . . . . . . HERE PASS=.T., SO XSKIM PASSAGE IS STILL AHEAD
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
C . . . . . . HERE PASS=.F., SO CHI IS MONOTONICALLY APPROACHING XTOT
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
      CP2   = IPHASX * SQRT(1.0D0 - SP2**2)
      ANUM  = CX1*CX2 - SX1*SX2*CP1*CP2
      DNOM  = 1.0D0 - SX1*SX2*SP1*SP2
      OMEGX = ACOS(ANUM/DNOM)

      RPRE = R2

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
      DELT2  = SQRT(R2**2 + R1**2 - 2.0D0*R2*R1*COS(OMEGA))

      DELTX = (DELT2-DELTA)
      COMGA = (RPRE**2+DELTX**2-R2**2)/(2.0D0*RPRE*DELTX)
      H_TAU = 0.0D0
      O2_TAU = 0.0D0
      DO 171 ISUN = 1,16
        HCHI  = (XX16(ISUN) + 1.0D0) * DELTX / 2.0D0
        WCHI = DELTX * WW16(ISUN) / 2.0D0
        RR   = SQRT(RPRE**2+HCHI**2-2.0D0*RPRE*HCHI*COMGA)

        IF (RR .GT. RC) THEN
          CALL CORONA(RR,DNST)
          ONST = 0.0D0
        ELSE
          IHI = ITHERM + 1
          ILO = 1
 172      IF ((IHI-ILO).GT.1) THEN
            IKK = (IHI + ILO) / 2
            IF (THERMO(2,IKK).GT.RR) THEN
              IHI = IKK
            ELSE
              ILO = IKK
            END IF
            GO TO 172
          END IF
          IKK = ILO

          DZ  = THERMO(2,IKK+1) - THERMO(2,IKK)
          DR  = (RR - THERMO(2,IKK)) / DZ
          XPT = THERMO(4,IKK) + DR * (THERMO(4,IKK+1)-THERMO(4,IKK))
          DNST = EXP( XPT )
          XPT = THERMO(5,IKK) + DR * (THERMO(5,IKK+1)-THERMO(5,IKK))
          ONST = EXP( XPT )
        END IF
        H_TAU = H_TAU + WCHI * DNST
        O2_TAU = O2_TAU + WCHI * ONST
 171  CONTINUE
      H_TAU = CENTER * H_TAU
      O2_TAU = ABSCSX * O2_TAU

      FCTR   = SQRT(TEXO/TEMP(III))
      KROSS  = 1
      IF (TRIP(III,JJJ)) KROSS = 2
      ODO2(2,III,JJJ,KROSS)   = O2_TAU + ODO2(1,III,JJJ,KROSS)
      ODLC(2,0,III,JJJ,KROSS) = H_TAU  + ODLC(1,0,III,JJJ,KROSS)
      IF (NONIT) THEN
        DO LP = 1,LPROF
          XPNT = (FCTR * SPDPT(LP))**2
          ODLC(2,LP,III,JJJ,KROSS) = H_TAU * FCTR * EXP(-XPNT)
     &                               + ODLC(1,LP,III,JJJ,KROSS)
        END DO
      END IF
      TRIP(III,JJJ) = .TRUE.
      IF ((IIN.LT.1).OR.(IIN.GT.IKNT)) GO TO 222
      KROSN = 1
      IF (TRIP(IIN,JJN)) KROSN = 2
      ODO2(1,IIN,JJN,KROSN)   = ODO2(2,III,JJJ,KROSS)
      ODLC(1,0,IIN,JJN,KROSN) = ODLC(2,0,III,JJJ,KROSS)
      IF (NONIT) THEN
        DO LP = 1,LPROF
          ODLC(1,LP,IIN,JJN,KROSN) = ODLC(2,LP,III,JJJ,KROSS)
        END DO
      END IF
      DELTA = DELT2
      III   = IIN
      JJJ   = JJN
      GO TO 111

C -----------------------------------------------------------------------------

 222  CONTINUE
      IF (ITRIP.EQ.0) THEN
        IEXIT = II
        JEXIT = JJ
      ELSE
        IEXIT = III
        JEXIT = JJJ
      END IF

      RETURN
      END

C =============================================================================
C =============================================================================

C INTERPOLATION ROUTINE FOR LINE-INTEGRATED TRANSMISSION FUNCTIONS
C (LINEAR INTERPOLATION WRT LINE-CENTER OPTICAL DEPTH)

      FUNCTION TRANS(XX)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION TBLX(101),TBLY(101)

      DATA TBLX / -4.6052D0, -4.5131D0, -4.4210D0, -4.3289D0, -4.2368D0,
     & -4.1447D0, -4.0525D0, -3.9604D0, -3.8683D0, -3.7762D0,
     & -3.6841D0, -3.5920D0, -3.4999D0, -3.4078D0, -3.3157D0,
     & -3.2236D0, -3.1315D0, -3.0394D0, -2.9473D0, -2.8552D0,
     & -2.7631D0, -2.6710D0, -2.5789D0, -2.4868D0, -2.3947D0,
     & -2.3026D0, -2.2105D0, -2.1184D0, -2.0263D0, -1.9342D0,
     & -1.8421D0, -1.7500D0, -1.6579D0, -1.5658D0, -1.4737D0,
     & -1.3816D0, -1.2894D0, -1.1973D0, -1.1052D0, -1.0131D0,
     & -9.2103D-1,-8.2893D-1,-7.3683D-1,-6.4472D-1,-5.5262D-1,
     & -4.6052D-1,-3.6841D-1,-2.7631D-1,-1.8421D-1,-9.2103D-2,
     & -1.1102D-16,9.2103D-2, 1.8421D-1, 2.7631D-1, 3.6841D-1,
     &  4.6052D-1, 5.5262D-1, 6.4472D-1, 7.3683D-1, 8.2893D-1,
     &  9.2103D-1, 1.0131D0,  1.1052D0,  1.1973D0,  1.2894D0,
     &  1.3816D0,  1.4737D0,  1.5658D0,  1.6579D0,  1.7500D0,
     &  1.8421D0,  1.9342D0,  2.0263D0,  2.1184D0,  2.2105D0,
     &  2.3026D0,  2.3947D0,  2.4868D0,  2.5789D0,  2.6710D0,
     &  2.7631D0,  2.8552D0,  2.9473D0,  3.0394D0,  3.1315D0,
     &  3.2236D0,  3.3157D0,  3.4078D0,  3.4999D0,  3.5920D0,
     &  3.6841D0,  3.7762D0,  3.8683D0,  3.9604D0,  4.0525D0,
     &  4.1447D0,  4.2368D0,  4.3289D0,  4.4210D0,  4.5131D0, 4.6052D0 /
      DATA TBLY / 9.9296D-1, 9.9228D-1, 9.9154D-1, 9.9073D-1, 9.8984D-1,
     & 9.8887D-1, 9.8780D-1, 9.8663D-1, 9.8535D-1, 9.8395D-1,
     & 9.8242D-1, 9.8074D-1, 9.7891D-1, 9.7690D-1, 9.7470D-1,
     & 9.7230D-1, 9.6968D-1, 9.6681D-1, 9.6367D-1, 9.6025D-1,
     & 9.5651D-1, 9.5243D-1, 9.4799D-1, 9.4314D-1, 9.3785D-1,
     & 9.3209D-1, 9.2583D-1, 9.1902D-1, 9.1162D-1, 9.0358D-1,
     & 8.9486D-1, 8.8542D-1, 8.7519D-1, 8.6414D-1, 8.5221D-1,
     & 8.3935D-1, 8.2551D-1, 8.1064D-1, 7.9469D-1, 7.7763D-1,
     & 7.5943D-1, 7.4004D-1, 7.1945D-1, 6.9766D-1, 6.7466D-1,
     & 6.5048D-1, 6.2516D-1, 5.9875D-1, 5.7134D-1, 5.4302D-1,
     & 5.1393D-1, 4.8422D-1, 4.5407D-1, 4.2367D-1, 3.9326D-1,
     & 3.6306D-1, 3.3334D-1, 3.0433D-1, 2.7630D-1, 2.4946D-1,
     & 2.2404D-1, 2.0021D-1, 1.7811D-1, 1.5783D-1, 1.3943D-1,
     & 1.2288D-1, 1.0816D-1, 9.5146D-2, 8.3732D-2, 7.3766D-2,
     & 6.5091D-2, 5.7548D-2, 5.0986D-2, 4.5265D-2, 4.0265D-2,
     & 3.5878D-2, 3.2018D-2, 2.8610D-2, 2.5593D-2, 2.2916D-2,
     & 2.0537D-2, 1.8418D-2, 1.6529D-2, 1.4843D-2, 1.3336D-2,
     & 1.1988D-2, 1.0782D-2, 9.7011D-3, 8.7322D-3, 7.8630D-3,
     & 7.0829D-3, 6.3823D-3, 5.7528D-3, 5.1869D-3, 4.6780D-3,
     & 4.2202D-3, 3.8082D-3, 3.4372D-3, 3.1030D-3, 2.8020D-3, 2.5307D-3/

C -----------------------------------------------------------------------------

      IF (XX.LE.0.01D0) THEN
        TRANS = 1.0D0 - XX / 1.414213562D0
      ELSE IF (XX.GE.100.0D0) THEN
        TRANS = 1.0D0 / XX / SQRT(3.1415926536D0 * LOG(XX))
      ELSE
        XL    = LOG(XX)
        IBIN  = INT((XL+4.605170186D0)/0.09210340372D0) + 1
C        PRINT*,XX,XL,BIN
        DEL   = TBLY(IBIN+1) - TBLY(IBIN)
        XREF  = EXP(TBLX(IBIN))
        DENOM = EXP(TBLX(IBIN+1)) - XREF
        TRANS = TBLY(IBIN) + DEL*(XX-XREF)/DENOM
      END IF
      RETURN
      END

C =============================================================================
C =============================================================================

C ITERATIVE PROCEDURE TO DETERMINE THETA THAT GRAZES THE PLANETARY SHADOW, TREATED AS A RIGHT CYLINDER.  MUST HAVE X1+THETA > PI.

      SUBROUTINE GRAZE(AA,BB,RBASE,R1,X1,THETA,CX1,SX1,CUMU,SUMU,
     &                 OGRAZ,PHIG)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON / NMBR/  PI,RTPI,PID2,OFFSET

C NOW TO ENTER INTO ITERATIVE PROCEDURE
      OG1 = THETA - ASIN(R1 * SUMU / AA)
 123  CONTINUE
        XX   = PI - ASIN(RBASE*SIN(THETA-OG1)/R1/SUMU)
        RR   = RBASE / SIN(XX)
        CCSQ = R1**2 + RBASE**2 * (1.0D0/TAN(XX)**2 - 1.0D0)
     &         - 2.0D0*R1*RBASE*CX1/TAN(XX)
        COSO = (R1**2 + RR**2 - CCSQ) / (2.0D0*R1*RR)
        OG2  = ACOS(COSO)
        TEST = (OG2 - OG1) / OG2
        IF (ABS(TEST).LT.OFFSET) GO TO 234
        OG1  = OG2
        GO TO 123
 234  CONTINUE
      XX   = PI - ASIN(RBASE*SIN(THETA-OG2)/R1/SUMU)
      RR   = RBASE / SIN(XX)
      CCSQ = R1**2 + RBASE**2 * (1.0D0/TAN(XX)**2 - 1.0D0)
     &       - 2.0D0*R1*RBASE*CX1/TAN(XX)
      EESQ = AA**2 + RR**2 + 2.0D0*AA*RR*COS(XX)
      GG   = (BB**2 + CCSQ - EESQ) / (2.0D0*BB*SQRT(CCSQ))
      PHIG = PI - ACOS((GG - CUMU**2)/SUMU**2)
      OGRAZ= OG2

      RETURN
      END

C =============================================================================
C =============================================================================

C ITERATIVE PROCEDURE TO DETERMINE LOS DISPLACEMENT ANGLE TO POINT
C OF DEEPEST PENETRATION INTO SHADOW

      SUBROUTINE DARKST(OGRAZ,R1,X1,THETA,CX1,SX1,CUMU,SUMU,CP1,ODEEP)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON / NMBR/  PI,RTPI,PID2,OFFSET

C NOW TO ENTER INTO ITERATIVE PROCEDURE
      OG1 = OGRAZ
 123  CONTINUE
        XX   = ACOS(CX1*COS(OG1) + SX1*SIN(OG1)*CP1)
        RR   = R1 * SUMU / SIN(THETA-OG1)
        FFSQ = R1**2 + (RR*COS(XX))**2 - 2.0D0*R1*RR*COS(XX)*CX1
        SS   = RR * SIN(XX)
        COSO = (R1**2 + RR**2 - (FFSQ - SS**2)) / (2.0D0*R1*RR)
        OG2  = ACOS(COSO)
        TEST = (OG2 - OG1) / OG2
        IF (ABS(TEST).LT.OFFSET) GO TO 234
        OG1 = OG2
        GO TO 123
 234  CONTINUE
      ODEEP = OG2

      RETURN
      END

C =============================================================================
C =============================================================================

C SIMPLE-MINDED WAY TO LOCATE SHADOW INTERSECTION POINTS.

      SUBROUTINE CROSS(RBASE,RSMU,R1,THETA,CX1,SX1,CP1,ILLUM,OS1,
     &                 ISHAD,DSHAD1,DSHAD2)
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL ILLUM
      COMMON / NMBR/  PI,RTPI,PID2,OFFSET

      DSHAD1 = 0.0D0
      DSHAD2 = 0.0D0

      OMAX = OS1
      OMIN = 0.0D0
 123  CONTINUE
        OS2  = (OMAX + OMIN) / 2.0D0
        RS   = RSMU / SIN(THETA-OS2)
        COSX = CX1*COS(OS2) + SX1*SIN(OS2)*CP1
        XS   = ACOS(COSX)
        TEST = (RS * SIN(XS) - RBASE) / RBASE
        IF (ABS(TEST).LT.OFFSET) GO TO 234
        IF (.NOT.ILLUM) TEST = -TEST
        IF (TEST.LT.0.0D0) THEN
C . . . ON THE "FAR SIDE"
          OMAX = OS2
        ELSE
C . . . ON THE "NEAR SIDE"
          OMIN = OS2
        END IF
        GO TO 123
 234  CONTINUE
      DSHAD1 = SQRT(R1**2 + RS**2 - 2.0D0*R1*RS*COS(OS2))
      IF (ISHAD.EQ.2) THEN
        RDEEP = RSMU / SIN(THETA-OS1)
        DDEEP = SQRT(R1**2 + RDEEP**2 - 2.0D0*R1*RDEEP*COS(OS1))
        DSHAD2 = DDEEP + (DDEEP - DSHAD1)
      END IF

C . . . NOW TO OFFSET DSHAD FROM "TRUE" VALUE BY 10 M (ROUND-OFF PROTECTION)
C ***** offset protection increased to 500 m (JBishop 6Feb96)
      IF (ILLUM) THEN
        DSHAD1 = DSHAD1 - 500.0D2
        IF (DSHAD1.LE.OFFSET) DSHAD1 = OFFSET
        IF (ISHAD.EQ.2) DSHAD2 = DSHAD2 + 500.0D2
      ELSE
        DSHAD1 = DSHAD1 + 500.0D2
      END IF
      RETURN
      END

C =============================================================================
C =============================================================================

C DETERMINATION OF LOCAL DENSITIES & TEMPERATURES USING MSIS90 TABLE

      SUBROUTINE DLOCAL(ITHERM,THERMO,RL,HLOC,O2LOC,TLOC)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION THERMO(5,ITHERM+1)
      COMMON /EXOS/  GM,CONL,PLANETR,RBASE,RC,RUPR,RP,TEXO,DEXO,VELT,
     &               RADPF,CENTER,ABSCSX,FTSAT,FDSAT,IGEO
      COMMON /NMBR/  PI,RTPI,PID2,OFFSET

C  EXOSPHERE
      IF (RL.GT.(RC*(1.0d0 + OFFSET))) THEN
        CALL CORONA(RL,HLOC)
        O2LOC = 0.0D0
        TLOC  = TEXO

C  WANDERED BELOW DESIGNATED PERMITTED REGION (roundoff errors)
      ELSE IF (RL.LE.THERMO(2,1)) THEN
        HLOC  = EXP( THERMO(4,1) )
        O2LOC = EXP( THERMO(5,1) )
        TLOC  = THERMO(3,1)

C  THERMOSPHERE
      ELSE
C . . . FIND CURRENT RADIAL BIN
        IHI = ITHERM + 1
        ILO = 1
 101    IF ((IHI-ILO).GT.1) THEN
          II = (IHI + ILO) / 2
          IF (THERMO(2,II).GT.RL) THEN 
            IHI = II
          ELSE
            ILO = II
          END IF
          GO TO 101
        END IF
C . . . INTERPOLATE
        RATIO = (RL-THERMO(2,ILO)) / (THERMO(2,IHI)-THERMO(2,ILO))
        XPNT  = THERMO(4,ILO) + RATIO * (THERMO(4,IHI)-THERMO(4,ILO))
        HLOC  = EXP( XPNT )
        XPNT  = THERMO(5,ILO) + RATIO * (THERMO(5,IHI)-THERMO(5,ILO))
        O2LOC = EXP( XPNT )
        TLOC  = THERMO(3,ILO) + RATIO * (THERMO(3,IHI)-THERMO(3,ILO))

      END IF
      RETURN
      END

C =============================================================================
C =============================================================================

C ROUTINE TO EVALUATE INCREMENTS TO LINE-OF-SIGHT OPTICAL DEPTHS
C . . . TO2LOC:  LOCAL INCREMENT TO O2 PHOTOABSORPTION OPTICAL DEPTH
C . . . TLPLOC:  LOCAL INCREMENT TO H  SCATTERING OPTICAL DEPTH
C                TLPLOC(LP = 0): LINE CENTER TAU INCREMENT AT TEXO
C                TLPLOC(LP > 0): TAU INCREMENT AT DISPLACEMENT SPDPT(LP)
C                                AT LOCAL TEMPERATURE

      SUBROUTINE STEP_TAU(ITHERM,THERMO,R1,THETA,D2,D1,
     &                    LPROF,SPDPT,NONIT,TO2LOC,TLPLOC)
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL NONIT
      DIMENSION THERMO(5,ITHERM+1),SPDPT(LPROF),TLPLOC(0 : LPROF)
      COMMON /EXOS/  GM,CONL,PLANETR,RBASE,RC,RUPR,RP,TEXO,DEXO,VELT,
     &               RADPF,CENTER,ABSCSX,FTSAT,FDSAT,IGEO
      COMMON /NMBR/  PI,RTPI,PID2,OFFSET
      COMMON /GAUS08/WW08(08),XX08(08)

      TO2LOC = 0.0D0
      DO 101 LP = 0,LPROF
        TLPLOC(LP) = 0.0D0
 101  END DO

      DO 201 ITAU = 1,8
        DTL  = ((D2-D1)*XX08(ITAU) + D2+D1) / 2.0D0
        WDTL =  (D2-D1)*WW08(ITAU) / 2.0D0
        RTL  = SQRT(R1**2 + DTL**2 - 2.0D0*R1*DTL*COS(PI-THETA))
        CALL DLOCAL(ITHERM,THERMO,RTL,HLOC,O2LOC,TLOC)
        TO2LOC    = TO2LOC    + WDTL * ABSCSX * O2LOC
        TLPLOC(0) = TLPLOC(0) + WDTL * CENTER * HLOC
        IF (NONIT) THEN
          FCTR    = SQRT(TEXO/TLOC)
          CLUSTER = (WDTL*CENTER*HLOC)*FCTR
          DO 202 LP = 1,LPROF
            XPNT = (SPDPT(LP)*FCTR)**2
            TLPLOC(LP) = TLPLOC(LP) + CLUSTER * EXP(-XPNT)
 202      END DO
        END IF
 201  END DO

      RETURN
      END

C =============================================================================
C =============================================================================

C ROUTINE TO INCREMENTALLY EVALUATE LINE-OF-SIGHT OPTICAL DEPTHS
C . . . TO2LOC:  LOCAL INCREMENT TO O2 PHOTOABSORPTION OPTICAL DEPTH
C . . . TLPLOC:  LOCAL INCREMENT TO H  SCATTERING OPTICAL DEPTH
C                TLPLOC(LP = 0): LINE CENTER TAU INCREMENT AT TEXO
C                TLPLOC(LP > 0): TAU INCREMENT AT DISPLACEMENT SPDPT(LP)
C                                AT LOCAL TEMPERATURE

      SUBROUTINE INLINE_TAU(DEL1,HDEN1,O2DEN1,TEMP1,
     &                      DEL2,HDEN2,O2DEN2,TEMP2,
     &                      LPROF,SPDPT,NONIT,TO2LOC,TLPLOC)
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL NONIT
      DIMENSION SPDPT(LPROF),TLPLOC(0 : LPROF)
      COMMON /EXOS/  GM,CONL,PLANETR,RBASE,RC,RUPR,RP,TEXO,DEXO,VELT,
     &               RADPF,CENTER,ABSCSX,FTSAT,FDSAT,IGEO
      COMMON /NMBR/  PI,RTPI,PID2,OFFSET

      HLOC  = (HDEN1  + HDEN2)  / 2.0d0
      O2LOC = (O2DEN1 + O2DEN2) / 2.0d0
      TLOC  = (TEMP1  + TEMP2)  / 2.0d0
      TO2LOC    = TO2LOC    + (DEL2-DEL1) * ABSCSX * O2LOC
      TLPLOC(0) = TLPLOC(0) + (DEL2-DEL1) * CENTER * HLOC
      IF (NONIT) THEN
        FCTR    = SQRT(TEXO/TLOC)
        CLUSTER = ((DEL2-DEL1)*CENTER*HLOC)*FCTR
        DO 202 LP = 1,LPROF
          XPNT = (SPDPT(LP)*FCTR)**2
          TLPLOC(LP) = TLPLOC(LP) + CLUSTER * EXP(-XPNT)
 202    END DO
      END IF
      DEL1   = DEL2
      HDEN1  = HDEN2
      O2DEN1 = O2DEN2
      TEMP1  = TEMP2

      RETURN
      END

C =============================================================================
C =============================================================================

C SOURCE FUNCTION INTERPOLATION SET-UP ROUTINE.
C PASSES BACK:
C   XSMS  . . . PADDED RADIAL INTERPOLATION ARRAY.
C   D0SMS . . . 2-D ARRAY OF THE LOGARITHMS OF THE MULTIPLE SCATTERING
C               SOURCE FUNCTIONS AT THE POINTS DEFINED BY THE XSMS
C               AND CHIN ARRAYS.
C   D2SMS . . . SECOND (RADIAL) DERIVATIVES OF THE MULTIPLE SCATTERING
C               SOURCE FUNCTION LOGARITHMS AT THE POINTS DEFINED BY
C               THE XSMS AND CHIN ARRAYS.

      SUBROUTINE SMS_SPLINE(IKNT,JKNT,BRAD,RADN,SOL,SRC,
     &                      XSMS,D0SMS,D2SMS)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NMAX = 12000)      ! must be at least as large as IKNT+2
      DIMENSION BRAD(IKNT+1),RADN(IKNT),SOL(IKNT,JKNT),SRC(IKNT,JKNT)
      DIMENSION XSMS(IKNT+2),D0SMS(IKNT+2,JKNT),D2SMS(IKNT+2,JKNT)
      DIMENSION YSMS(NMAX),D2YSMS(NMAX)
      COMMON /NMBR/  PI,RTPI,PID2,OFFSET

c set up array of radial points
      DO 101 ISPLN = 1,IKNT
        XSMS(ISPLN+1) = RADN(ISPLN)
 101  END DO
c . . . lowermost point:  bottom boundary radius
      XSMS(1) = BRAD(1)
c . . . uppermost point:  top boundary radius
      XSMS(IKNT+2) = BRAD(IKNT+1)
       
C        PRINT*,'----'
C        PRINT*,BRAD
C        PRINT*,XSMS

c set up array of logarithms of multiple scattering source function
      DO 201 JJ = 1,JKNT
        DO 202 ISPLN = 1,IKNT
c currently have SRC & SOL archived to the LYAO_source.DAT file retaining
c five significant digits.  possibility arises of SRC-SOL=0 in optically
c thin situations, leading to NaN results.  protection trigger is set
c one decade down (current value of OFFSET is 1.0D-6). (JBishop 10Jul96)
          DELSMS = (SRC(ISPLN,JJ) - SOL(ISPLN,JJ)) / SRC(ISPLN,JJ)
          IF (DELSMS.GT.OFFSET) THEN
            D0SMS(ISPLN+1,JJ) = LOG(SRC(ISPLN,JJ) - SOL(ISPLN,JJ))
          ELSE
            D0SMS(ISPLN+1,JJ) = LOG(OFFSET * SRC(ISPLN,JJ))
          END IF
 202    END DO
c . . . lowermost point:  extrapolate to bottom boundary radius
        slope = (D0SMS(3,JJ)-D0SMS(2,JJ)) / (XSMS(3)-XSMS(2))
        D0SMS(1,JJ) = D0SMS(2,JJ) + (XSMS(1)-XSMS(2))*SLOPE
c . . . uppermost point:  extrapolate to top boundary radius
        slope = (D0SMS(IKNT+1,JJ)-D0SMS(IKNT,JJ)) /
     &          (XSMS(IKNT+1)-XSMS(IKNT))
        D0SMS(IKNT+2,JJ) = D0SMS(IKNT+1,JJ) +
     &                     (XSMS(IKNT+2)-XSMS(IKNT+1))*SLOPE
 201  END DO

c now to break into individual 1-D arrays, obtain second derivatives,
c then repack into 2-D array
      DO 301 JJ = 1,JKNT
        DO 302 ISPLN = 1,IKNT+2
          YSMS(ISPLN) = D0SMS(ISPLN,JJ)
 302    END DO
        CALL SPLINE_LYAO(XSMS,YSMS,IKNT+2,1.0D30,1.0D30,D2YSMS)
        DO 303 ISPLN = 1,IKNT+2
          D2SMS(ISPLN,JJ) = D2YSMS(ISPLN)
 303    END DO
 301  END DO

C        PRINT*,'----'
C        PRINT*,XSMS
C        PRINT*,D0SMS
C        PRINT*,D2SMS

      RETURN
      END

C -----------------------------------------------------------------------------

C SPLINE SUBROUTINE, MODIFIED TO DOUBLE PRECISION
      SUBROUTINE SPLINE_LYAO(X,Y,N,YP1,YPN,Y2)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NMAX=12000)
      DIMENSION X(N),Y(N),Y2(N),U(NMAX)
      IF (YP1.GT.0.99D30) THEN
        Y2(1)=0.0D0
        U(1) =0.0D0
      ELSE
        Y2(1)=-0.5D0
        U(1) =(3.0D0/(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1)
      ENDIF
      DO 11 I=2,N-1
        SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1))
        P  =SIG*Y2(I-1)+2.0D0
        Y2(I)=(SIG-1.0D0)/P
        U(I) =(6.0D0*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1))
     *       /(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P
11    CONTINUE
      IF (YPN.GT.0.99D30) THEN
        QN=0.0D0
        UN=0.0D0
      ELSE
        QN=0.50D0
        UN=(3.0D0/(X(N)-X(N-1)))*(YPN-(Y(N)-Y(N-1))/(X(N)-X(N-1)))
      ENDIF
      Y2(N)=(UN-QN*U(N-1))/(QN*Y2(N-1)+1.0D0)
      DO 12 K=N-1,1,-1
        Y2(K)=Y2(K)*Y2(K+1)+U(K)
12    CONTINUE
      RETURN
      END

C -----------------------------------------------------------------------------

      SUBROUTINE SPLINT_LYAO(XA,YA,Y2A,N,X,Y)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XA(N),YA(N),Y2A(N)
      KLO=1
      KHI=N
1     IF (KHI-KLO.GT.1) THEN
        K=(KHI+KLO)/2
        IF(XA(K).GT.X)THEN
          KHI=K
        ELSE
          KLO=K
        ENDIF
      GOTO 1
      ENDIF
      H=XA(KHI)-XA(KLO)
      IF (H.EQ.0.0D0) STOP 'Bad XA input.'
      A=(XA(KHI)-X)/H
      B=(X-XA(KLO))/H
      Y=A*YA(KLO)+B*YA(KHI)+
     *      ((A**3-A)*Y2A(KLO)+(B**3-B)*Y2A(KHI))*(H**2)/6.0D0
      RETURN
      END

C =============================================================================
C =============================================================================
C MATRIX INVERSION SUBROUTINES

      SUBROUTINE LUDCMP(A,N,INDX,D)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NMAX=12000,TINY=1.0D-20)
      DIMENSION A(N,N),INDX(N),VV(NMAX)
      D = 1.0D0
      DO 12 I = 1,N
        AAMAX = 0.0D0
        DO 11 J = 1,N
          IF (ABS(A(I,J)).GT.AAMAX) AAMAX = ABS(A(I,J))
 11     CONTINUE
        IF (AAMAX.EQ.0.0D0) STOP 'Singular matrix.'
        VV(I) = 1.0D0/AAMAX
 12   CONTINUE
      DO 19 J = 1,N
        IF (J.GT.1) THEN
          DO 14 I = 1,J-1
            SUM = A(I,J)
            IF (I.GT.1) THEN
              DO 13 K = 1,I-1
                SUM = SUM - A(I,K)*A(K,J)
 13           CONTINUE
              A(I,J) = SUM
            ENDIF
 14       CONTINUE
        ENDIF
        AAMAX = 0.0D0
        DO 16 I = J,N
          SUM = A(I,J)
          IF (J.GT.1) THEN
            DO 15 K = 1,J-1
              SUM = SUM - A(I,K)*A(K,J)
 15         CONTINUE
            A(I,J) = SUM
          ENDIF
          DUM = VV(I)*ABS(SUM)
          IF (DUM.GE.AAMAX) THEN
            IMAX  = I
            AAMAX = DUM
          ENDIF
 16     CONTINUE
        IF (J.NE.IMAX) THEN
          DO 17 K = 1,N
            DUM = A(IMAX,K)
            A(IMAX,K) = A(J,K)
            A(J,K) = DUM
 17       CONTINUE
          D = -D
          VV(IMAX) = VV(J)
        ENDIF
        INDX(J) = IMAX
        IF(J.NE.N) THEN
          IF (A(J,J).EQ.0.0D0) A(J,J) = TINY
          DUM = 1.0D0/A(J,J)
          DO 18 I = J+1,N
            A(I,J) = A(I,J)*DUM
 18       CONTINUE
        ENDIF
 19   CONTINUE
      IF (A(N,N).EQ.0.0D0) A(N,N) = TINY
      RETURN
      END

C -----------------------------------------------------------------------------

      SUBROUTINE LUBKSB(A,N,INDX,B)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(N,N),INDX(N),B(N)
      II = 0
      DO 12 I = 1,N
        LL    = INDX(I)
        SUM   = B(LL)
        B(LL) = B(I)
        IF (II.NE.0) THEN
          DO 11 J = II,I-1
            SUM = SUM - A(I,J)*B(J)
 11       CONTINUE
        ELSE IF (SUM.NE.0.0D0) THEN
          II = I
        ENDIF
        B(I) = SUM
 12   CONTINUE
      DO 14 I = N,1,-1
        SUM = B(I)
        IF (I.LT.N) THEN
          DO 13 J = I+1,N
            SUM = SUM - A(I,J)*B(J)
 13       CONTINUE
        ENDIF
        B(I) = SUM/A(I,I)
 14   CONTINUE
      RETURN
      END

C -----------------------------------------------------------------------------

      SUBROUTINE MPROVE(A,ALUD,N,INDX,B,X)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NMAX=12000)
      DIMENSION A(N,N),ALUD(N,N),INDX(N),B(N),X(N),R(NMAX)
      DO 12 I = 1,N
        SDP = -B(I)
        DO 11 J = 1,N
          SDP = SDP + A(I,J)*X(J)
 11     CONTINUE
        R(I) = SDP
 12   CONTINUE
      CALL LUBKSB(ALUD,N,INDX,R)
      DO 13 I = 1,N
        X(I) = X(I) - R(I)
 13   CONTINUE
      RETURN
      END

C =============================================================================
C =============================================================================
C GAUSSIAN POINTS & WEIGHTS

      BLOCK DATA GAUSS04
      IMPLICIT REAL*8 (A-Z)
      COMMON /GAUS04/ WW04(4),XX04(4)
      DATA WW04( 1),XX04( 1) /  3.47854853E-01 , -8.61136317E-01 /
      DATA WW04( 2),XX04( 2) /  6.52145147E-01 , -3.39981049E-01 /
      DATA WW04( 3),XX04( 3) /  6.52145147E-01 ,  3.39981049E-01 /
      DATA WW04( 4),XX04( 4) /  3.47854853E-01 ,  8.61136317E-01 /
      END

      BLOCK DATA GAUSS08
      IMPLICIT REAL*8 (A-Z)
      COMMON /GAUS08/ WW08(8),XX08(8)
      DATA WW08( 1),XX08( 1) /  1.01228535E-01 , -9.60289836E-01 /
      DATA WW08( 2),XX08( 2) /  2.22381040E-01 , -7.96666503E-01 /
      DATA WW08( 3),XX08( 3) /  3.13706636E-01 , -5.25532424E-01 /
      DATA WW08( 4),XX08( 4) /  3.62683773E-01 , -1.83434635E-01 /
      DATA WW08( 5),XX08( 5) /  3.62683773E-01 ,  1.83434635E-01 /
      DATA WW08( 6),XX08( 6) /  3.13706636E-01 ,  5.25532424E-01 /
      DATA WW08( 7),XX08( 7) /  2.22381040E-01 ,  7.96666503E-01 /
      DATA WW08( 8),XX08( 8) /  1.01228535E-01 ,  9.60289836E-01 /
      END

      BLOCK DATA GAUSS16
      IMPLICIT REAL*8 (A-Z)
      COMMON /GAUS16/ WW16(16),XX16(16)
      DATA WW16( 1),XX16( 1) /  2.71524601E-02 , -9.89400923E-01 /
      DATA WW16( 2),XX16( 2) /  6.22535236E-02 , -9.44575012E-01 /
      DATA WW16( 3),XX16( 3) /  9.51585099E-02 , -8.65631223E-01 /
      DATA WW16( 4),XX16( 4) /  1.24628969E-01 , -7.55404413E-01 /
      DATA WW16( 5),XX16( 5) /  1.49595991E-01 , -6.17876232E-01 /
      DATA WW16( 6),XX16( 6) /  1.69156522E-01 , -4.58016783E-01 /
      DATA WW16( 7),XX16( 7) /  1.82603419E-01 , -2.81603545E-01 /
      DATA WW16( 8),XX16( 8) /  1.89450607E-01 , -9.50125083E-02 /
      DATA WW16( 9),XX16( 9) /  1.89450607E-01 ,  9.50125083E-02 /
      DATA WW16(10),XX16(10) /  1.82603419E-01 ,  2.81603545E-01 /
      DATA WW16(11),XX16(11) /  1.69156522E-01 ,  4.58016783E-01 /
      DATA WW16(12),XX16(12) /  1.49595991E-01 ,  6.17876232E-01 /
      DATA WW16(13),XX16(13) /  1.24628969E-01 ,  7.55404413E-01 /
      DATA WW16(14),XX16(14) /  9.51585099E-02 ,  8.65631223E-01 /
      DATA WW16(15),XX16(15) /  6.22535236E-02 ,  9.44575012E-01 /
      DATA WW16(16),XX16(16) /  2.71524601E-02 ,  9.89400923E-01 /
      END

      BLOCK DATA GAUSS24
      IMPLICIT REAL*8 (A-Z)
      COMMON /GAUS24/ WW24(24),XX24(24)
      DATA WW24( 1),XX24( 1) /  1.23412302E-02 , -9.95187223E-01 /
      DATA WW24( 2),XX24( 2) /  2.85313893E-02 , -9.74728584E-01 /
      DATA WW24( 3),XX24( 3) /  4.42774370E-02 , -9.38274562E-01 /
      DATA WW24( 4),XX24( 4) /  5.92985861E-02 , -8.86415541E-01 /
      DATA WW24( 5),XX24( 5) /  7.33464807E-02 , -8.20001960E-01 /
      DATA WW24( 6),XX24( 6) /  8.61901641E-02 , -7.40124166E-01 /
      DATA WW24( 7),XX24( 7) /  9.76186544E-02 , -6.48093641E-01 /
      DATA WW24( 8),XX24( 8) /  1.07444271E-01 , -5.45421481E-01 /
      DATA WW24( 9),XX24( 9) /  1.15505666E-01 , -4.33793515E-01 /
      DATA WW24(10),XX24(10) /  1.21670470E-01 , -3.15042675E-01 /
      DATA WW24(11),XX24(11) /  1.25837460E-01 , -1.91118866E-01 /
      DATA WW24(12),XX24(12) /  1.27938196E-01 , -6.40568957E-02 /
      DATA WW24(13),XX24(13) /  1.27938196E-01 ,  6.40568957E-02 /
      DATA WW24(14),XX24(14) /  1.25837460E-01 ,  1.91118866E-01 /
      DATA WW24(15),XX24(15) /  1.21670470E-01 ,  3.15042675E-01 /
      DATA WW24(16),XX24(16) /  1.15505666E-01 ,  4.33793515E-01 /
      DATA WW24(17),XX24(17) /  1.07444271E-01 ,  5.45421481E-01 /
      DATA WW24(18),XX24(18) /  9.76186544E-02 ,  6.48093641E-01 /
      DATA WW24(19),XX24(19) /  8.61901641E-02 ,  7.40124166E-01 /
      DATA WW24(20),XX24(20) /  7.33464807E-02 ,  8.20001960E-01 /
      DATA WW24(21),XX24(21) /  5.92985861E-02 ,  8.86415541E-01 /
      DATA WW24(22),XX24(22) /  4.42774370E-02 ,  9.38274562E-01 /
      DATA WW24(23),XX24(23) /  2.85313893E-02 ,  9.74728584E-01 /
      DATA WW24(24),XX24(24) /  1.23412302E-02 ,  9.95187223E-01 /
      END

      BLOCK DATA GAUSS32
      IMPLICIT REAL*8 (A-Z)
      COMMON /GAUS32/ WW32(32),XX32(32)
      DATA WW32( 1),XX32( 1) /  7.01860990E-03 , -9.97263849E-01 /
      DATA WW32( 2),XX32( 2) /  1.62743945E-02 , -9.85611498E-01 /
      DATA WW32( 3),XX32( 3) /  2.53920648E-02 , -9.64762270E-01 /
      DATA WW32( 4),XX32( 4) /  3.42738628E-02 , -9.34906065E-01 /
      DATA WW32( 5),XX32( 5) /  4.28358987E-02 , -8.96321177E-01 /
      DATA WW32( 6),XX32( 6) /  5.09980582E-02 , -8.49367619E-01 /
      DATA WW32( 7),XX32( 7) /  5.86840920E-02 , -7.94483781E-01 /
      DATA WW32( 8),XX32( 8) /  6.58222213E-02 , -7.32182145E-01 /
      DATA WW32( 9),XX32( 9) /  7.23457932E-02 , -6.63044274E-01 /
      DATA WW32(10),XX32(10) /  7.81938955E-02 , -5.87715745E-01 /
      DATA WW32(11),XX32(11) /  8.33119228E-02 , -5.06899893E-01 /
      DATA WW32(12),XX32(12) /  8.76520947E-02 , -4.21351284E-01 /
      DATA WW32(13),XX32(13) /  9.11738798E-02 , -3.31868589E-01 /
      DATA WW32(14),XX32(14) /  9.38443989E-02 , -2.39287362E-01 /
      DATA WW32(15),XX32(15) /  9.56387222E-02 , -1.44471958E-01 /
      DATA WW32(16),XX32(16) /  9.65400860E-02 , -4.83076647E-02 /
      DATA WW32(17),XX32(17) /  9.65400860E-02 ,  4.83076647E-02 /
      DATA WW32(18),XX32(18) /  9.56387222E-02 ,  1.44471958E-01 /
      DATA WW32(19),XX32(19) /  9.38443989E-02 ,  2.39287362E-01 /
      DATA WW32(20),XX32(20) /  9.11738798E-02 ,  3.31868589E-01 /
      DATA WW32(21),XX32(21) /  8.76520947E-02 ,  4.21351284E-01 /
      DATA WW32(22),XX32(22) /  8.33119228E-02 ,  5.06899893E-01 /
      DATA WW32(23),XX32(23) /  7.81938955E-02 ,  5.87715745E-01 /
      DATA WW32(24),XX32(24) /  7.23457932E-02 ,  6.63044274E-01 /
      DATA WW32(25),XX32(25) /  6.58222213E-02 ,  7.32182145E-01 /
      DATA WW32(26),XX32(26) /  5.86840920E-02 ,  7.94483781E-01 /
      DATA WW32(27),XX32(27) /  5.09980582E-02 ,  8.49367619E-01 /
      DATA WW32(28),XX32(28) /  4.28358987E-02 ,  8.96321177E-01 /
      DATA WW32(29),XX32(29) /  3.42738628E-02 ,  9.34906065E-01 /
      DATA WW32(30),XX32(30) /  2.53920648E-02 ,  9.64762270E-01 /
      DATA WW32(31),XX32(31) /  1.62743945E-02 ,  9.85611498E-01 /
      DATA WW32(32),XX32(32) /  7.01860990E-03 ,  9.97263849E-01 /
      END

      BLOCK DATA GAUSS64
      IMPLICIT REAL*8 (A-Z)
      COMMON /GAUS64/ WW64(64),XX64(64)
      DATA WW64( 1),XX64( 1) /  1.78328075E-03 , -9.99305069E-01 /
      DATA WW64( 2),XX64( 2) /  4.14703321E-03 , -9.96340096E-01 /
      DATA WW64( 3),XX64( 3) /  6.50445791E-03 , -9.91013348E-01 /
      DATA WW64( 4),XX64( 4) /  8.84675980E-03 , -9.83336270E-01 /
      DATA WW64( 5),XX64( 5) /  1.11681391E-02 , -9.73326802E-01 /
      DATA WW64( 6),XX64( 6) /  1.34630483E-02 , -9.61008787E-01 /
      DATA WW64( 7),XX64( 7) /  1.57260299E-02 , -9.46411371E-01 /
      DATA WW64( 8),XX64( 8) /  1.79517157E-02 , -9.29569185E-01 /
      DATA WW64( 9),XX64( 9) /  2.01348234E-02 , -9.10522163E-01 /
      DATA WW64(10),XX64(10) /  2.22701747E-02 , -8.89315426E-01 /
      DATA WW64(11),XX64(11) /  2.43527032E-02 , -8.65999401E-01 /
      DATA WW64(12),XX64(12) /  2.63774693E-02 , -8.40629280E-01 /
      DATA WW64(13),XX64(13) /  2.83396728E-02 , -8.13265324E-01 /
      DATA WW64(14),XX64(14) /  3.02346572E-02 , -7.83972383E-01 /
      DATA WW64(15),XX64(15) /  3.20579298E-02 , -7.52819896E-01 /
      DATA WW64(16),XX64(16) /  3.38051617E-02 , -7.19881833E-01 /
      DATA WW64(17),XX64(17) /  3.54722142E-02 , -6.85236335E-01 /
      DATA WW64(18),XX64(18) /  3.70551273E-02 , -6.48965478E-01 /
      DATA WW64(19),XX64(19) /  3.85501534E-02 , -6.11155331E-01 /
      DATA WW64(20),XX64(20) /  3.99537422E-02 , -5.71895659E-01 /
      DATA WW64(21),XX64(21) /  4.12625633E-02 , -5.31279445E-01 /
      DATA WW64(22),XX64(22) /  4.24735136E-02 , -4.89403158E-01 /
      DATA WW64(23),XX64(23) /  4.35837246E-02 , -4.46366012E-01 /
      DATA WW64(24),XX64(24) /  4.45905589E-02 , -4.02270168E-01 /
      DATA WW64(25),XX64(25) /  4.54916283E-02 , -3.57220173E-01 /
      DATA WW64(26),XX64(26) /  4.62847948E-02 , -3.11322868E-01 /
      DATA WW64(27),XX64(27) /  4.69681844E-02 , -2.64687151E-01 /
      DATA WW64(28),XX64(28) /  4.75401655E-02 , -2.17423648E-01 /
      DATA WW64(29),XX64(29) /  4.79993895E-02 , -1.69644415E-01 /
      DATA WW64(30),XX64(30) /  4.83447611E-02 , -1.21462822E-01 /
      DATA WW64(31),XX64(31) /  4.85754684E-02 , -7.29931220E-02 /
      DATA WW64(32),XX64(32) /  4.86909561E-02 , -2.43502930E-02 /
      DATA WW64(33),XX64(33) /  4.86909561E-02 ,  2.43502930E-02 /
      DATA WW64(34),XX64(34) /  4.85754684E-02 ,  7.29931220E-02 /
      DATA WW64(35),XX64(35) /  4.83447611E-02 ,  1.21462822E-01 /
      DATA WW64(36),XX64(36) /  4.79993895E-02 ,  1.69644415E-01 /
      DATA WW64(37),XX64(37) /  4.75401655E-02 ,  2.17423648E-01 /
      DATA WW64(38),XX64(38) /  4.69681844E-02 ,  2.64687151E-01 /
      DATA WW64(39),XX64(39) /  4.62847948E-02 ,  3.11322868E-01 /
      DATA WW64(40),XX64(40) /  4.54916283E-02 ,  3.57220173E-01 /
      DATA WW64(41),XX64(41) /  4.45905589E-02 ,  4.02270168E-01 /
      DATA WW64(42),XX64(42) /  4.35837246E-02 ,  4.46366012E-01 /
      DATA WW64(43),XX64(43) /  4.24735136E-02 ,  4.89403158E-01 /
      DATA WW64(44),XX64(44) /  4.12625633E-02 ,  5.31279445E-01 /
      DATA WW64(45),XX64(45) /  3.99537422E-02 ,  5.71895659E-01 /
      DATA WW64(46),XX64(46) /  3.85501534E-02 ,  6.11155331E-01 /
      DATA WW64(47),XX64(47) /  3.70551273E-02 ,  6.48965478E-01 /
      DATA WW64(48),XX64(48) /  3.54722142E-02 ,  6.85236335E-01 /
      DATA WW64(49),XX64(49) /  3.38051617E-02 ,  7.19881833E-01 /
      DATA WW64(50),XX64(50) /  3.20579298E-02 ,  7.52819896E-01 /
      DATA WW64(51),XX64(51) /  3.02346572E-02 ,  7.83972383E-01 /
      DATA WW64(52),XX64(52) /  2.83396728E-02 ,  8.13265324E-01 /
      DATA WW64(53),XX64(53) /  2.63774693E-02 ,  8.40629280E-01 /
      DATA WW64(54),XX64(54) /  2.43527032E-02 ,  8.65999401E-01 /
      DATA WW64(55),XX64(55) /  2.22701747E-02 ,  8.89315426E-01 /
      DATA WW64(56),XX64(56) /  2.01348234E-02 ,  9.10522163E-01 /
      DATA WW64(57),XX64(57) /  1.79517157E-02 ,  9.29569185E-01 /
      DATA WW64(58),XX64(58) /  1.57260299E-02 ,  9.46411371E-01 /
      DATA WW64(59),XX64(59) /  1.34630483E-02 ,  9.61008787E-01 /
      DATA WW64(60),XX64(60) /  1.11681391E-02 ,  9.73326802E-01 /
      DATA WW64(61),XX64(61) /  8.84675980E-03 ,  9.83336270E-01 /
      DATA WW64(62),XX64(62) /  6.50445791E-03 ,  9.91013348E-01 /
      DATA WW64(63),XX64(63) /  4.14703321E-03 ,  9.96340096E-01 /
      DATA WW64(64),XX64(64) /  1.78328075E-03 ,  9.99305069E-01 /
      END

      BLOCK DATA GAUSS96
      IMPLICIT REAL*8 (A-Z)
      COMMON /GAUS96/ WW96(96),XX96(96)
      DATA WW96( 1),XX96( 1) /  7.96792039E-04 , -9.99689519E-01 /
      DATA WW96( 2),XX96( 2) /  1.85396080E-03 , -9.98364389E-01 /
      DATA WW96( 3),XX96( 3) /  2.91073183E-03 , -9.95981872E-01 /
      DATA WW96( 4),XX96( 4) /  3.96455452E-03 , -9.92543876E-01 /
      DATA WW96( 5),XX96( 5) /  5.01420256E-03 , -9.88054097E-01 /
      DATA WW96( 6),XX96( 6) /  6.05854532E-03 , -9.82517242E-01 /
      DATA WW96( 7),XX96( 7) /  7.09647080E-03 , -9.75939155E-01 /
      DATA WW96( 8),XX96( 8) /  8.12687725E-03 , -9.68326807E-01 /
      DATA WW96( 9),XX96( 9) /  9.14867129E-03 , -9.59688306E-01 /
      DATA WW96(10),XX96(10) /  1.01607703E-02 , -9.50032711E-01 /
      DATA WW96(11),XX96(11) /  1.11621022E-02 , -9.39370334E-01 /
      DATA WW96(12),XX96(12) /  1.21516045E-02 , -9.27712440E-01 /
      DATA WW96(13),XX96(13) /  1.31282294E-02 , -9.15071428E-01 /
      DATA WW96(14),XX96(14) /  1.40909422E-02 , -9.01460648E-01 /
      DATA WW96(15),XX96(15) /  1.50387213E-02 , -8.86894524E-01 /
      DATA WW96(16),XX96(16) /  1.59705635E-02 , -8.71388495E-01 /
      DATA WW96(17),XX96(17) /  1.68854799E-02 , -8.54959011E-01 /
      DATA WW96(18),XX96(18) /  1.77825019E-02 , -8.37623537E-01 /
      DATA WW96(19),XX96(19) /  1.86606795E-02 , -8.19400311E-01 /
      DATA WW96(20),XX96(20) /  1.95190813E-02 , -8.00308764E-01 /
      DATA WW96(21),XX96(21) /  2.03567967E-02 , -7.80369043E-01 /
      DATA WW96(22),XX96(22) /  2.11729407E-02 , -7.59602368E-01 /
      DATA WW96(23),XX96(23) /  2.19666436E-02 , -7.38030672E-01 /
      DATA WW96(24),XX96(24) /  2.27370691E-02 , -7.15676785E-01 /
      DATA WW96(25),XX96(25) /  2.34833993E-02 , -6.92564547E-01 /
      DATA WW96(26),XX96(26) /  2.42048409E-02 , -6.68718338E-01 /
      DATA WW96(27),XX96(27) /  2.49006338E-02 , -6.44163430E-01 /
      DATA WW96(28),XX96(28) /  2.55700368E-02 , -6.18925869E-01 /
      DATA WW96(29),XX96(29) /  2.62123402E-02 , -5.93032360E-01 /
      DATA WW96(30),XX96(30) /  2.68268660E-02 , -5.66510439E-01 /
      DATA WW96(31),XX96(31) /  2.74129622E-02 , -5.39388120E-01 /
      DATA WW96(32),XX96(32) /  2.79700067E-02 , -5.11694193E-01 /
      DATA WW96(33),XX96(33) /  2.84974109E-02 , -4.83457983E-01 /
      DATA WW96(34),XX96(34) /  2.89946143E-02 , -4.54709411E-01 /
      DATA WW96(35),XX96(35) /  2.94610895E-02 , -4.25478995E-01 /
      DATA WW96(36),XX96(36) /  2.98963450E-02 , -3.95797640E-01 /
      DATA WW96(37),XX96(37) /  3.02999150E-02 , -3.65696847E-01 /
      DATA WW96(38),XX96(38) /  3.06713767E-02 , -3.35208535E-01 /
      DATA WW96(39),XX96(39) /  3.10103334E-02 , -3.04364949E-01 /
      DATA WW96(40),XX96(40) /  3.13164257E-02 , -2.73198813E-01 /
      DATA WW96(41),XX96(41) /  3.15893292E-02 , -2.41743162E-01 /
      DATA WW96(42),XX96(42) /  3.18287574E-02 , -2.10031316E-01 /
      DATA WW96(43),XX96(43) /  3.20344567E-02 , -1.78096876E-01 /
      DATA WW96(44),XX96(44) /  3.22062038E-02 , -1.45973712E-01 /
      DATA WW96(45),XX96(45) /  3.23438235E-02 , -1.13695852E-01 /
      DATA WW96(46),XX96(46) /  3.24471630E-02 , -8.12974945E-02 /
      DATA WW96(47),XX96(47) /  3.25161181E-02 , -4.88129854E-02 /
      DATA WW96(48),XX96(48) /  3.25506143E-02 , -1.62767451E-02 /
      DATA WW96(49),XX96(49) /  3.25506143E-02 ,  1.62767451E-02 /
      DATA WW96(50),XX96(50) /  3.25161181E-02 ,  4.88129854E-02 /
      DATA WW96(51),XX96(51) /  3.24471630E-02 ,  8.12974945E-02 /
      DATA WW96(52),XX96(52) /  3.23438235E-02 ,  1.13695852E-01 /
      DATA WW96(53),XX96(53) /  3.22062038E-02 ,  1.45973712E-01 /
      DATA WW96(54),XX96(54) /  3.20344567E-02 ,  1.78096876E-01 /
      DATA WW96(55),XX96(55) /  3.18287574E-02 ,  2.10031316E-01 /
      DATA WW96(56),XX96(56) /  3.15893292E-02 ,  2.41743162E-01 /
      DATA WW96(57),XX96(57) /  3.13164257E-02 ,  2.73198813E-01 /
      DATA WW96(58),XX96(58) /  3.10103334E-02 ,  3.04364949E-01 /
      DATA WW96(59),XX96(59) /  3.06713767E-02 ,  3.35208535E-01 /
      DATA WW96(60),XX96(60) /  3.02999150E-02 ,  3.65696847E-01 /
      DATA WW96(61),XX96(61) /  2.98963450E-02 ,  3.95797640E-01 /
      DATA WW96(62),XX96(62) /  2.94610895E-02 ,  4.25478995E-01 /
      DATA WW96(63),XX96(63) /  2.89946143E-02 ,  4.54709411E-01 /
      DATA WW96(64),XX96(64) /  2.84974109E-02 ,  4.83457983E-01 /
      DATA WW96(65),XX96(65) /  2.79700067E-02 ,  5.11694193E-01 /
      DATA WW96(66),XX96(66) /  2.74129622E-02 ,  5.39388120E-01 /
      DATA WW96(67),XX96(67) /  2.68268660E-02 ,  5.66510439E-01 /
      DATA WW96(68),XX96(68) /  2.62123402E-02 ,  5.93032360E-01 /
      DATA WW96(69),XX96(69) /  2.55700368E-02 ,  6.18925869E-01 /
      DATA WW96(70),XX96(70) /  2.49006338E-02 ,  6.44163430E-01 /
      DATA WW96(71),XX96(71) /  2.42048409E-02 ,  6.68718338E-01 /
      DATA WW96(72),XX96(72) /  2.34833993E-02 ,  6.92564547E-01 /
      DATA WW96(73),XX96(73) /  2.27370691E-02 ,  7.15676785E-01 /
      DATA WW96(74),XX96(74) /  2.19666436E-02 ,  7.38030672E-01 /
      DATA WW96(75),XX96(75) /  2.11729407E-02 ,  7.59602368E-01 /
      DATA WW96(76),XX96(76) /  2.03567967E-02 ,  7.80369043E-01 /
      DATA WW96(77),XX96(77) /  1.95190813E-02 ,  8.00308764E-01 /
      DATA WW96(78),XX96(78) /  1.86606795E-02 ,  8.19400311E-01 /
      DATA WW96(79),XX96(79) /  1.77825019E-02 ,  8.37623537E-01 /
      DATA WW96(80),XX96(80) /  1.68854799E-02 ,  8.54959011E-01 /
      DATA WW96(81),XX96(81) /  1.59705635E-02 ,  8.71388495E-01 /
      DATA WW96(82),XX96(82) /  1.50387213E-02 ,  8.86894524E-01 /
      DATA WW96(83),XX96(83) /  1.40909422E-02 ,  9.01460648E-01 /
      DATA WW96(84),XX96(84) /  1.31282294E-02 ,  9.15071428E-01 /
      DATA WW96(85),XX96(85) /  1.21516045E-02 ,  9.27712440E-01 /
      DATA WW96(86),XX96(86) /  1.11621022E-02 ,  9.39370334E-01 /
      DATA WW96(87),XX96(87) /  1.01607703E-02 ,  9.50032711E-01 /
      DATA WW96(88),XX96(88) /  9.14867129E-03 ,  9.59688306E-01 /
      DATA WW96(89),XX96(89) /  8.12687725E-03 ,  9.68326807E-01 /
      DATA WW96(90),XX96(90) /  7.09647080E-03 ,  9.75939155E-01 /
      DATA WW96(91),XX96(91) /  6.05854532E-03 ,  9.82517242E-01 /
      DATA WW96(92),XX96(92) /  5.01420256E-03 ,  9.88054097E-01 /
      DATA WW96(93),XX96(93) /  3.96455452E-03 ,  9.92543876E-01 /
      DATA WW96(94),XX96(94) /  2.91073183E-03 ,  9.95981872E-01 /
      DATA WW96(95),XX96(95) /  1.85396080E-03 ,  9.98364389E-01 /
      DATA WW96(96),XX96(96) /  7.96792039E-04 ,  9.99689519E-01 /
      END

C =============================================================================
C =============================================================================
C SPEED POINTS & WEIGHTS

      BLOCK DATA SPEED04
      IMPLICIT REAL*8 (A-Z)
      COMMON /SPW04/ SWGT04(4),SPNT04(4)
      DATA SWGT04(1), SPNT04(1) / 0.3253029998D0,  0.1337764470D0  /
      DATA SWGT04(2), SPNT04(2) / 0.4211071019D0,  0.6243246902D0  /
      DATA SWGT04(3), SPNT04(3) / 0.1334425004D0,  1.3425378256D0  /
      DATA SWGT04(4), SPNT04(4) / 0.6374323486D-2, 2.2626644770D0  /
      END

      BLOCK DATA SPEED08
      IMPLICIT REAL*8 (A-Z)
      COMMON /SPW08/ SWGT08(8),SPNT08(8)
      DATA SWGT08(1), SPNT08(1) / 0.1341091885D0,  0.5297864393D-1 /
      DATA SWGT08(2), SPNT08(2) / 0.2683307545D0,  0.2673983722D0  /
      DATA SWGT08(3), SPNT08(3) / 0.2759533980D0,  0.6163028842D0  /
      DATA SWGT08(4), SPNT08(4) / 0.1574482826D0,  1.0642463121D0  /
      DATA SWGT08(5), SPNT08(5) / 0.4481410992D-1, 1.5888558623D0  /
      DATA SWGT08(6), SPNT08(6) / 0.5367935756D-2, 2.1839211531D0  /
      DATA SWGT08(7), SPNT08(7) / 0.2020636491D-3, 2.8631338837D0  /
      DATA SWGT08(8), SPNT08(8) / 0.1192596927D-5, 3.6860071627D0  /
      END

      BLOCK DATA SPEED16
      IMPLICIT REAL*8 (A-Z)
      COMMON /SPW16/ SWGT16(16),SPNT16(16)
      DATA SWGT16(1), SPNT16(1) / 0.5052463202D-1, 0.1975365846D-1 /
      DATA SWGT16(2), SPNT16(2) / 0.1136085569D0,  0.1028022452D0  /
      DATA SWGT16(3), SPNT16(3) / 0.1629212923D0,  0.2473976695D0  /
      DATA SWGT16(4), SPNT16(4) / 0.1835628011D0,  0.4466962260D0  /
      DATA SWGT16(5), SPNT16(5) / 0.1654386378D0,  0.6930737203D0  /
      DATA SWGT16(6), SPNT16(6) / 0.1165724906D0,  0.9794041703D0  /
      DATA SWGT16(7), SPNT16(7) / 0.6199969610D-1, 1.2997893213D0  /
      DATA SWGT16(8), SPNT16(8) / 0.2391970962D-1, 1.6498542404D0  /
      DATA SWGT16(9) ,SPNT16(9) / 0.6409914424D-2, 2.0268081522D0  /
      DATA SWGT16(10),SPNT16(10)/ 0.1135695311D-2, 2.4294504916D0  /
      DATA SWGT16(11),SPNT16(11)/ 0.1252862213D-3, 2.8582665285D0  /
      DATA SWGT16(12),SPNT16(12)/ 0.7950495720D-5, 3.3157692750D0  /
      DATA SWGT16(13),SPNT16(13)/ 0.2590007619D-6, 3.8073771168D0  /
      DATA SWGT16(14),SPNT16(14)/ 0.3611549140D-8, 4.3436063455D0  /
      DATA SWGT16(15),SPNT16(15)/ 0.1537677916D-10,4.9463772040D0  /
      DATA SWGT16(16),SPNT16(16)/ 0.8674204452D-14,5.6750179340D0  /
      END

C =============================================================================
