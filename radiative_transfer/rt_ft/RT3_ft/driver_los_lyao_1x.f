c example driver based on STP 78-1 subsolar limb profile
C modified by Pratik Joshi [JUL2023]
C -----------------------------------------------------------------------------

      PROGRAM DRIVER
      IMPLICIT REAL*8 (A-H,O-Z)
C     sequential    1x 1024     NFI 1024*1024   WFI 512*512
C     parallel 8    1x  128     NFI  128*1024   WFI 128*512
C     parallel 16   1x   64     NFI   64*1024   WFI  64*512
C     parallel 32   1x   32     NFI   32*1024   WFI  32*512
C     inversion     6 (INV6   MAXLOS = 6)
C     NFI 1024*1024, WFI 512*512
      PARAMETER (MAXLOS = 1*1024)
      PARAMETER (LPROF  =   8)

      DIMENSION ROBS(MAXLOS), SZA(MAXLOS), ZNTH(MAXLOS), AZI(MAXLOS)
      DIMENSION TPALT(MAXLOS), ACER(3,MAXLOS), THETLS(MAXLOS)
      DIMENSION SPDPT(LPROF), WSPD(LPROF)

      REAL*8 FLUX_C

      COMMON/SOLAR_FLUX/ FLUX_C

      PLANETR = 6371.0D5
      ILOS = MAXLOS
      LINE_LABEL = 1
      WAVELN = 1215.67

C      OPEN(98, FILE='solar_flux.txt', STATUS='OLD')
C      READ(98,*) FLUX_C
C      CLOSE(98)

      OPEN(96, FILE='orbitinfo.txt',STATUS='OLD')
C      READ(96,*)
      DO JK = 1,MAXLOS
        READ(96,*) ROBS(JK),SZA(JK),ZNTH(JK),AZI(JK)
      END DO
      CLOSE(96)

      DO JK = 1,MAXLOS
         ROBS(JK) = ROBS(JK)*1.0D5
         IF (ZNTH(JK) .GT. 180.0D0) THEN
             ZNTH(JK) = 360.0D0-ZNTH(JK)
         END IF
      END DO


C -----------------------------------------------------------------------------
C the big call

C      CALL GET_MSIS(LINE_LABEL)
      CALL global_parameters
      CALL LYAO_LOS(LINE_LABEL,MAXLOS,LPROF,WAVELN,
     &              SPDPT,WSPD,ILOS,ROBS,SZA,ZNTH,AZI,TPALT,ACER)

C -----------------------------------------------------------------------------
C      PRINT*, "FLUX_C", FLUX_C

      OPEN(10, FILE='lyao_los.DAT', STATUS='unknown')

C 10   FORMAT(1X, 2I4, 0P1F9.2, 1P1E11.2)
 11   FORMAT(1X,1P1E10.3,0P3F15.9,5X,G22.16,1X,G22.16,1X,G22.16)

      DO LK = 1,ILOS
        FOURPI_I = ACER(3,LK) * FLUX_C
        RATIO    = ACER(1,LK) / ACER(3,LK)
        WRITE(10,11) ROBS(LK), SZA(LK), ZNTH(LK), AZI(LK),
     &              TPALT(LK), FOURPI_I, RATIO
      END DO

      CLOSE(10)

      END
