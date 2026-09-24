c example driver based on STP 78-1 subsolar limb profile
C modified by Pratik Joshi [JUL2023]
C -----------------------------------------------------------------------------

      SUBROUTINE FT_MAIN(los,losx,losy,sf,
     &                   orbit,orbx,orby,
     &                   exopt,exoptx,cdens,cdensx,
     &                   th,thx,thy,
     &                   br,brx,bch,bchx,chn,chnx,
     &                   rhot,rhotx,rhoty,
     &                   ts,tsx,tsy,tsz,
     &                   msis_p1x,msis_p1,msis_p2x,msis_p2,
     &                   ap_histx,ap_hist,
!      &                   rtbkgx,rtbkgy,rtbkgz,rtbkg,
     &                   thbkg,thbkgx,thbkgy,thbkgz,
     &                   zgrid,zgridx)
        IMPLICIT REAL*8 (A-H,O-Z)
        INTEGER losx,losy
        REAL*8 los(losx,losy)
        REAL*8 sf
        INTEGER orbx,orby
        REAL*8 orbit(orbx,orby)
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
        INTEGER msis_p1x
        INTEGER msis_p1(msis_p1x)
        INTEGER msis_p2x
        REAL*8 msis_p2(msis_p2x)
        INTEGER ap_histx
        REAL*4 ap_hist(ap_histx)
      !   INTEGER rtbkgx,rtbkgy,rtbkgz
      !   REAL*8 rtbkg(rtbkgx,rtbkgy,rtbkgz)
        INTEGER thbkgx,thbkgy,thbkgz
        REAL*8 thbkg(thbkgx,thbkgy,thbkgz)
        INTEGER zgridx
        REAL*8 zgrid(zgridx)
Cf2py intent(in) losx,losy
Cf2py intent(out) los
Cf2py depend(losx,losy) los
Cf2py intent(in) sf
Cf2py intent(in) orbx,orby
Cf2py intent(in) orbit
Cf2py intent(in) exoptx
Cf2py intent(in) exopt
Cf2py intent(in) cdensx
Cf2py intent(in) cdens
Cf2py intent(in) thx,thy
Cf2py intent(in) th
Cf2py intent(in) brx
Cf2py intent(in) br
Cf2py intent(in) bchx
Cf2py intent(in) bch
Cf2py intent(in) chnx
Cf2py intent(in) chn
Cf2py intent(in) rhotx,rhoty
Cf2py intent(in) rhot
Cf2py intent(in) tsx,tsy,tsz
Cf2py intent(in) ts
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

C     sequential    1x 1024     NFI 1024*1024   WFI 512*512
C     parallel 8    1x  128     NFI  128*1024   WFI 128*512
C     parallel 16   1x   64     NFI   64*1024   WFI  64*512
C     parallel 32   1x   32     NFI   32*1024   WFI  32*512
C     inversion     6 (INV6   MAXLOS = 6)
C     NFI 1024*1024, WFI 512*512
      PARAMETER (MAXLOS = 512)
      PARAMETER (LPROF  =   8)

C       FOURPI_I(MAXLOS), RATIO(MAXLOS) are only for wrapper output
      DIMENSION FOURPI_I(MAXLOS), RATIO(MAXLOS)
      DIMENSION ROBS(MAXLOS), SZA(MAXLOS), ZNTH(MAXLOS), AZI(MAXLOS)
      DIMENSION TPALT(MAXLOS), ACER(3,MAXLOS), THETLS(MAXLOS)
      DIMENSION SPDPT(LPROF), WSPD(LPROF)

      REAL*8 FLUX_C

      COMMON/SOLAR_FLUX/ FLUX_C

      PLANETR = 6371.0D5
      ILOS = MAXLOS
      LINE_LABEL = 1
      WAVELN = 1215.67

C READ SOLAR_FLUX
      FLUX_C = sf

C READ ORBI INFO
      ROBS = orbit(:,1)
      SZA = orbit(:,2)
      ZNTH = orbit(:,3)
      AZI = orbit(:,4)

      DO JK = 1,MAXLOS
         ROBS(JK) = ROBS(JK)*1.0D5
         IF (ZNTH(JK) .GT. 180.0D0) THEN
             ZNTH(JK) = 360.0D0-ZNTH(JK)
         END IF
      END DO


C -----------------------------------------------------------------------------

      CALL global_parameters(msis_p1x,msis_p1,msis_p2x,msis_p2,
     &                       ap_histx,ap_hist,
!      &                       rtbkgx,rtbkgy,rtbkgz,rtbkg,
     &                       thbkgx,thbkgy,thbkgz,thbkg,
     &                       zgridx,zgrid)
      CALL LYAO_LOS(exoptx,exopt,cdensx,cdens,
     &              thx,thy,th,
     &              brx,br,bchx,bch,chnx,chn,
     &              rhotx,rhoty,rhot,
     &              tsx,tsy,tsz,ts,
     &              LINE_LABEL,MAXLOS,LPROF,WAVELN,
     &              SPDPT,WSPD,ILOS,ROBS,SZA,ZNTH,AZI,TPALT,ACER)

C -----------------------------------------------------------------------------
       DO LK = 1,ILOS
         FOURPI_I(LK) = ACER(3,LK) * FLUX_C
         RATIO(LK)    = ACER(1,LK) / ACER(3,LK)
       END DO

C WRITE OUTPUTS to los_1x array
        los(:,1) = ROBS
        los(:,2) = SZA
        los(:,3) = ZNTH
        los(:,4) = AZI
        los(:,5) = TPALT
        los(:,6) = FOURPI_I
        los(:,7) = RATIO

      END SUBROUTINE FT_MAIN
C END FILE driver_los_lyao.F
