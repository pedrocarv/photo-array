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

      SUBROUTINE CORONA(RR,DNST)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (GG   = 0.8862269D0)
      COMMON /EXOS/ GM,CONL,PLANETR,RBASE,RC,RUPR,RP,TEXO,DEXO,VELT,
     &              RADPF,CENTER,ABSCSX,FTSAT,FDSAT,IGEO
      COMMON /NMBR/ PI,RTPI,PID2,OFFSET
      DATA TSOL   / 1.10D6 /                ! SOLAR IONIZATION LIFETIME
C
C +++ ORIGINAL MODEL +++
C
      IF (IGEO.EQ.0) THEN
        RS  = FTSAT
        AC  = CONL/(RC*TEXO)
        AS  = CONL/(RS*TEXO)
        AA  = CONL/(RR*TEXO)
        IF (AA.GE.(AC*(1.0d0 - OFFSET))) THEN
          DNST = DEXO * EXP(AA-AC)
        ELSE IF (AA.LE.0.0D0) THEN
          DNST = 0.0D0
        ELSE
        YA = (AA**2)/(AC+AA)
        GA = GAMMA16(AA)
        GC = GAMMA16(AA-YA)
        CCOEF = SQRT(AC**2-AA**2) * EXP(-YA) / AC
        BALST = (GA - CCOEF*GC)*2.0D0/RTPI
        ESCAP = (1.0D0 - CCOEF)*GG/RTPI - BALST/2.0D0
        BOUND = 2.0D0*GA/RTPI
        IF (AA.LE.AS) THEN
          YS = (AA**2)/(AS+AA)
          GS = GAMMA16(AA-YS)
          SCOEF = SQRT(AS**2-AA**2) * EXP(-YS) / AS
          SATEL = (GA - SCOEF*GS)*2.0D0/RTPI - BALST
          PART  = BALST + SATEL + ESCAP
        ELSE
          PART  = BOUND + ESCAP
        END IF
        DNST = PART * DEXO * EXP(AA-AC)
      END IF
C
C +++ MODIFIED MODEL +++
C
      ELSE IF (IGEO.EQ.1) THEN
      AC  = CONL/(RC*TEXO)
      AR  = CONL/(RP*TEXO)
      AA  = CONL/(RR*TEXO)
      IF (AA.GE.(AC*(1.0d0 - OFFSET))) THEN
        DNST = DEXO * EXP(AA-AC)
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
        DNST  = DEXO * (PART1 + PART2)
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
