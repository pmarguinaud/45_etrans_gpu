MODULE EUVTVD_MOD
CONTAINS
SUBROUTINE EUVTVD(KFIELD,PU,PV,PVOR,PDIV)

!**** *EUVTVD* - Compute vor/div from u and v in spectral space

!     Purpose.
!     --------
!        To compute vorticity and divergence from u and v in spectral
!       space. Input u and v from KM to NTMAX+1, output vorticity and
!       divergence from KM to NTMAX - calculation part.

!**   Interface.
!     ----------
!        CALL EUVTVD(KM,KFIELD,PEPSNM,PU,PV,PVOR,PDIV)

!        Explicit arguments :  KM - zonal wave-number
!        --------------------  KFIELD - number of fields (levels)
!                              KFLDPTR - fields pointers
!                              PEPSNM - REPSNM for wavenumber KM
!                              PU - u wind component for zonal
!                                   wavenumber KM
!                              PV - v wind component for zonal
!                                   wavenumber KM
!                              PVOR - vorticity for zonal
!                                     wavenumber KM
!                              PDIV - divergence for zonal
!                                     wavenumber KM

!     Method.  See ref.
!     -------

!     Externals.  None.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Mats Hamrud and Philippe Courtier  *ECMWF*

!     Modifications.
!     --------------
!        Original : 91-07-01
!        D. Giard : NTMAX instead of NSMAX
!        01-08-27 : R. El Khatib Fix for NPROMATR /= 0
!        03-03-03 : G. Radnoti: b-level conform mean-wind distribution
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        F. Vana + NEC 28-Apr-2009 MPI-OpenMP fix
!        D. Degrauwe  (Feb 2012): Alternative extension zone (E')
!        N. Lopes & R. El Khatib 15-Jun-2012 Scalability enhancement
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

USE TPM_DIM         ,ONLY : R
!USE TPM_FIELDS
USE TPM_DISTR       ,ONLY : D, NPRCIDS, NPRTRW, MYSETV, MYSETW, MYPROC, NPROC
USE TPM_DISTR       ,ONLY : D_NUMP,D_MYMS
USE TPMALD_GEO      ,ONLY : GALD
USE TPMALD_DISTR    ,ONLY : DALD
USE ABORT_TRANS_MOD ,ONLY : ABORT_TRANS
!

IMPLICIT NONE

INTEGER(KIND=JPIM), INTENT(IN) :: KFIELD
REAL(KIND=JPRB), INTENT(INOUT) :: PU  (:,:,:),PV  (:,:,:)
REAL(KIND=JPRB), INTENT(OUT)   :: PVOR(:,:,:),PDIV(:,:,:)

INTEGER(KIND=JPIM) :: II, IN, IR, J, JN
INTEGER(KIND=JPIM) :: KM, KMLOC, JNMAX

REAL(KIND=JPRB) :: ZKM
REAL(KIND=JPRB) :: ZIN
INTEGER(KIND=JPIM) :: JA,ITAG,ILEN,IFLD,ISND
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('EUVTVD_MOD:EUVTVD',0,ZHOOK_HANDLE)

!*       1.    COMPUTE U V FROM VORTICITY AND DIVERGENCE.
!              ------------------------------------------


!$acc parallel loop collapse(3) private(J,KMLOC,JN,IR,II,KM,ZKM) present (PVOR, PDIV, PU, PV)
DO J=1,KFIELD
  DO KMLOC=1,D_NUMP
    DO JN=1,R%NDGL+R%NNOEXTZG
      KM = D_MYMS(KMLOC)
      ZKM=REAL(KM,JPRB)*GALD%EXWN
      IR=2*J-1
      II=IR+1
      PDIV(JN,IR,KMLOC)=-ZKM*PU(JN,II,KMLOC)
      PDIV(JN,II,KMLOC)= ZKM*PU(JN,IR,KMLOC)
      PVOR(JN,IR,KMLOC)=-ZKM*PV(JN,II,KMLOC)
      PVOR(JN,II,KMLOC)= ZKM*PV(JN,IR,KMLOC)
    ENDDO
  ENDDO
ENDDO
!$acc end parallel loop

JNMAX = MAXVAL(DALD%NCPL2M)

!$acc parallel loop collapse(3) private(J,KMLOC,JN,KM,ZIN,IN) copyin (JNMAX) present (PVOR, PDIV, PU, PV)
DO J=1,2*KFIELD
  DO KMLOC=1,D_NUMP
    DO JN=1,JNMAX,2
      KM = D_MYMS(KMLOC)
      IN=(JN-1)/2
      ZIN=REAL(IN,JPRB)*GALD%EYWN
      PVOR(JN,J  ,KMLOC)=PVOR(JN  ,J,KMLOC)+ZIN*PU(JN+1,J,KMLOC)
      PVOR(JN+1,J,KMLOC)=PVOR(JN+1,J,KMLOC)-ZIN*PU(JN  ,J,KMLOC)
      PDIV(JN,J  ,KMLOC)=PDIV(JN  ,J,KMLOC)-ZIN*PV(JN+1,J,KMLOC)
      PDIV(JN+1,J,KMLOC)=PDIV(JN+1,J,KMLOC)+ZIN*PV(JN  ,J,KMLOC)
    ENDDO
  ENDDO
ENDDO
!$acc end parallel loop

IF (LHOOK) CALL DR_HOOK('EUVTVD_MOD:EUVTVD',1,ZHOOK_HANDLE)

END SUBROUTINE EUVTVD
END MODULE EUVTVD_MOD
