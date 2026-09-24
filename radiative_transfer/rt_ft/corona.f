C =============================================================================

C EXOSPHERE DENSITY SUBROUTINE *** COMBINED *** (14 Sept 93)
C THIS CODE COMPUTES TOTAL DENSITY USING EITHER THE ORIGINAL
C CHAMBERLAIN [1963] ANALYTIC EXOSPHERE FORMULATION OR THE
C MODIFIED VERSION [BISHOP 1991]:
C       SPHERICALLY UNIFORM EXOBASE TEMPERATURE & DENSITY
C       EXOPAUSE HANDLING OF SATELLITE COMPONENT
C       NO ROTATION OR WINDS
C INCLUDES BOLTZMANN (OR "BAROMETRIC") FACTOR.
C SELECTION FLAG IS IGEO:
C     0 . . . ORIGINAL CHAMBERLAIN MODEL WITH SATELLITE CRITICAL
C             RADIUS (PASSED IN COMMON BLOCK AS FTSAT)
C     1 . . . MODIFIED MODEL WITH EXOPAUSE, "EVAPORATIVE" SATELLITE
C             POPULATION PARAMETERS, AND SOLAR IONIZATION DECAY.

C SPEED-UP (SEP 2026): THE SOLAR ANGLE ONLY ENTERS THROUGH THE EXOBASE
C DENSITY DENS(ANG), SO DNST = DENS(ANG) * G(RR).  G(RR) IS TABULATED ONCE
C PER MODEL BY CORONA_INIT (CALLED FROM GLOBAL_PARAMETERS) AND INTERPOLATED
C WITH A CUBIC SPLINE; OUTSIDE THE TABLE RANGE IT IS EVALUATED DIRECTLY BY
C CORONA_G (THE ORIGINAL FORMULATION).

      SUBROUTINE CORONA(RR, ANG, DNST)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NCTAB = 2049)
      COMMON /EXOS/ GM,CONL,PLANETR,RBASE,RC,RUPR,RP,TEXO,DEXO,PHZ_SHFT,
     &              TOGGLE,VELT,RADPF,CENTER,ABSCSX,FTSAT,FDSAT,IGEO
      COMMON /NMBR/ PI,RTPI,PID2,OFFSET
      COMMON /CORTAB/ ALTAB,AHTAB,GTAB(NCTAB),G2TAB(NCTAB),IFTAB,LOGTAB

C Uncomment for original functionality
      ! DENS = DEXO

C Exobase density varying with solar angle
      d_shift = PHZ_SHFT - (PI / 2)
      DENS = DEXO + TOGGLE * DEXO * 0.7 *(SIN(ANG + d_shift))/2

      AA = CONL/(RR*TEXO)
      IF ((IFTAB.EQ.0).OR.(AA.LT.ALTAB).OR.(AA.GE.AHTAB)) THEN
        CALL CORONA_G(RR, G)
      ELSE
C . . . TABLE NODES ARE UNIFORM IN T, WHERE AA = AL + (AH-AL)*(1-COS(PI*T))/2
        X  = (AA-ALTAB)/(AHTAB-ALTAB)
        T  = ACOS(1.0D0 - 2.0D0*X) / PI * DBLE(NCTAB-1)
        I  = MIN(INT(T) + 1, NCTAB-1)
        B  = T - DBLE(I-1)
        A  = 1.0D0 - B
        G  = A*GTAB(I) + B*GTAB(I+1) +
     &       ((A**3-A)*G2TAB(I) + (B**3-B)*G2TAB(I+1)) / 6.0D0
        IF (LOGTAB.EQ.1) G = EXP(G)
      END IF
      DNST = DENS * G
      RETURN
      END

C -----------------------------------------------------------------------------

C BUILDS THE G(RR) TABLE FOR THE CURRENT EXOSPHERE PARAMETERS IN /EXOS/.
C NODES ARE CLUSTERED AT BOTH ENDS (COSINE MAPPING IN AA = CONL/(RR*TEXO))
C SO THE SQRT BEHAVIOUR OF THE PARTITION FUNCTIONS AT THE EXOBASE AND THE
C EXOPAUSE IS INTERPOLATED SMOOTHLY.  THE ENDPOINTS STOP JUST SHORT OF THE
C BRANCH SWITCHES IN CORONA_G, WHICH ARE DISCONTINUOUS.

      SUBROUTINE CORONA_INIT
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NCTAB = 2049)
      COMMON /EXOS/ GM,CONL,PLANETR,RBASE,RC,RUPR,RP,TEXO,DEXO,PHZ_SHFT,
     &              TOGGLE,VELT,RADPF,CENTER,ABSCSX,FTSAT,FDSAT,IGEO
      COMMON /NMBR/ PI,RTPI,PID2,OFFSET
      COMMON /CORTAB/ ALTAB,AHTAB,GTAB(NCTAB),G2TAB(NCTAB),IFTAB,LOGTAB
      DIMENSION XT(NCTAB)

      IFTAB = 0
      TEMP  = TEXO
      AHTAB = CONL/(RC*TEMP) * (1.0D0 - OFFSET) * (1.0D0 - 1.0D-9)
C . . . TABLE REACHES A LITTLE BEYOND THE TOP OF THE RT GRID
      ALTAB = CONL/(1.01D0*RUPR*TEMP)
      IF (IGEO.EQ.1) THEN
C . . . . . . G DROPS TO ZERO AT THE EXOPAUSE RP
        ALTAB = MAX(ALTAB, CONL/(RP*TEMP) * (1.0D0 + 1.0D-9))
      ELSE IF (IGEO.EQ.0) THEN
C . . . . . . KINK AT THE SATELLITE CRITICAL RADIUS: NO TABLE IF IN RANGE
        AS = CONL/(FTSAT*TEMP)
        IF ((AS.GT.ALTAB).AND.(AS.LT.AHTAB)) RETURN
      ELSE
        RETURN
      END IF
      IF (ALTAB.GE.AHTAB) RETURN

      LOGTAB = 1
      DO I = 1,NCTAB
        XT(I) = DBLE(I-1)
        T  = DBLE(I-1) / DBLE(NCTAB-1)
        AA = ALTAB + (AHTAB-ALTAB) * (1.0D0 - COS(PI*T)) / 2.0D0
        CALL CORONA_G(CONL/(AA*TEMP), GTAB(I))
        IF (GTAB(I).LE.0.0D0) LOGTAB = 0
      END DO
      IF (LOGTAB.EQ.1) THEN
        DO I = 1,NCTAB
          GTAB(I) = LOG(GTAB(I))
        END DO
      END IF
C . . . NATURAL CUBIC SPLINE ON THE UNIT-SPACED NODE INDEX
      CALL SPLINE_LYAO(XT, GTAB, NCTAB, 1.0D31, 1.0D31, G2TAB)
      IFTAB = 1
      RETURN
      END

      BLOCK DATA CORTAB0
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NCTAB = 2049)
      COMMON /CORTAB/ ALTAB,AHTAB,GTAB(NCTAB),G2TAB(NCTAB),IFTAB,LOGTAB
      DATA IFTAB / 0 /
      END

C -----------------------------------------------------------------------------

C ORIGINAL CORONA FORMULATION WITH UNIT EXOBASE DENSITY: DNST = G(RR)

      SUBROUTINE CORONA_G(RR, DNST)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (GG   = 0.8862269D0)
      COMMON /EXOS/ GM,CONL,PLANETR,RBASE,RC,RUPR,RP,TEXO,DEXO,PHZ_SHFT,
     &              TOGGLE,VELT,RADPF,CENTER,ABSCSX,FTSAT,FDSAT,IGEO
      COMMON /NMBR/ PI,RTPI,PID2,OFFSET
C     SOLAR IONIZATION LIFETIME
      DATA TSOL   / 1.10D6 /

      DENS = 1.0D0
      ! t_shift = PHZ_SHFT + (PI / 2)
      ! TEMP = TEXO + TOGGLE * TEXO * 0.3 *(SIN(ANG + t_shift))/2
      TEMP = TEXO
C
C +++ ORIGINAL MODEL +++
C
C Comments refer to parameters in Chamberlain 1963
      IF (IGEO.EQ.0) THEN
C RS is the critical satellite raidal distance.
C FTSAT=TSAT*Re
        RS  = FTSAT
C AC,AS,AA are potential enegies inversely proportional
C to RC, RS, RR, and represented by lambdaC, lambdaS,
C and lambda, respectively
        AC  = CONL/(RC*TEMP)
        AS  = CONL/(RS*TEMP)
        AA  = CONL/(RR*TEMP)
        IF (AA.GE.(AC*(1.0d0 - OFFSET))) THEN
          DNST = DENS * EXP(AA-AC)
        ELSE IF (AA.LE.0.0D0) THEN
          DNST = 0.0D0
        ELSE
C Calculate psi1
        YA = (AA**2)/(AC+AA)
        GA = GAMMA16(AA)
        GC = GAMMA16(AA-YA)
        CCOEF = SQRT(AC**2-AA**2) * EXP(-YA) / AC
C Calculate ballistic partition function
        BALST = (GA - CCOEF*GC)*2.0D0/RTPI
C Calculate escape partition function
        ESCAP = (1.0D0 - CCOEF)*GG/RTPI - BALST/2.0D0
C Calculate bound partition function
C Satellite partition function = BOUND-BALST
        BOUND = 2.0D0*GA/RTPI
C If radial distance exceeds critical satellite distance
        IF (AA.LE.AS) THEN
          YS = (AA**2)/(AS+AA)
          GS = GAMMA16(AA-YS)
          SCOEF = SQRT(AS**2-AA**2) * EXP(-YS) / AS
          SATEL = (GA - SCOEF*GS)*2.0D0/RTPI - BALST
          PART  = BALST + SATEL + ESCAP
C If radial distance exceeds critical satellite distance
        ELSE
          PART  = BOUND + ESCAP
        END IF
C Calculate density using total partition function
        DNST = PART * DENS * EXP(AA-AC)
      END IF
C
C +++ MODIFIED MODEL +++
C
      ELSE IF (IGEO.EQ.1) THEN
      AC  = CONL/(RC*TEMP)
      AR  = CONL/(RP*TEMP)
      AA  = CONL/(RR*TEMP)
      IF (AA.GE.(AC*(1.0d0 - OFFSET))) THEN
        DNST = DENS * EXP(AA-AC)
      ELSE IF (AA.LE.AR) THEN
        DNST = 0.0D0
      ELSE
        Y1 = (AA**2)/(AC+AA)
        YA = AA-AR
        YB = AA - AR*AC/(AC+AR)
        YE = (AA**2)/(AR+AA)
        GA = GAMMA16(YA)
        GB = GAMMA16(YB-Y1)
        GS = GAMMA16((YB-Y1)/FTSAT)
        ZA = ZAMMA16(YE-YA)
        ZB = ZAMMA16(YE-YB)
        ZS = ZAMMA16((YE-YB)/FTSAT)
        ACOEF = 2.0D0 * SQRT(AA**2-AR**2) * EXP(-YE) / AR / RTPI
        CCOEF = 2.0D0 * SQRT(AC**2-AA**2) * EXP(-Y1) / AC / RTPI
        BALST = 2.0D0*GA/RTPI + ACOEF*(ZA-ZB) - CCOEF*GB
        ESCAP = (GG-GA)/RTPI - ACOEF*(ZA-ZB)/2.0D0 - CCOEF*(GG-GB)/2.0D0
        PART1 = (BALST+ESCAP) * EXP(AA-AC)
        SATEL = ACOEF * ZS * EXP((FTSAT-1.0D0)*YE/FTSAT) +
     &          CCOEF * GS * EXP((FTSAT-1.0D0)*Y1/FTSAT)
        DK    = 2.97D6 * ASIN(1.0D0 - RC/RR) / SQRT(RR) / RADPF
        SATEL = SATEL * EXP(-DK/TSOL)
        PART2 = FDSAT * SATEL * EXP((AA-AC)/FTSAT)
        DNST  = DENS * (PART1 + PART2)
      END IF
C +++
      END IF
      RETURN
      END
 
C -----------------------------------------------------------------------------
 
C  EVALUATION OF THE INCOMPLETE GAMMA FUNCTION
      FUNCTION GAMMA16(AAA)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON / GAUS16 / WW16(16),XX16(16)
C  INTEGRATION ROUTINE
      ROOT = SQRT(AAA)
      ERF  = 0.0D0
      DO 100 JJ = 1,16
        XPNT= ((1.0D0+XX16(JJ))*ROOT/2.0D0)**2
        ERF = ERF + WW16(JJ)*EXP(-XPNT)*ROOT/2.0D0
 100  CONTINUE
      GAMMA16 = ERF - ROOT*EXP(-AAA)
      RETURN
      END

C -----------------------------------------------------------------------------

C  EVALUATION OF THE ZAMMA FUNCTION
      FUNCTION ZAMMA16(AAA)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON / GAUS16 / WW16(16),XX16(16)
C INTEGRATION ROUTINE
      ZAMMA16 = 0.0D0
      DO 100 JJ = 1,16
      XXX   = (1.0D0+XX16(JJ))*AAA/2.0D0
      ZAMMA16 = ZAMMA16 + WW16(JJ)*EXP(XXX)*SQRT(XXX)*AAA/2.0D0
 100  CONTINUE
      RETURN
      END

C =============================================================================
