MODULE EPRFI1B_MOD
CONTAINS
SUBROUTINE EPRFI1B(PFFA,PSPEC,KFIELDS,KFLDPTR)

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

!USE TPM_DIM
USE TPM_DISTR
USE TPMALD_DISTR    ,ONLY : DALD, DALD_NESM0, DALD_NCPL2M
!
!**** *PRFI1* - Prepare spectral fields for inverse Legendre transform

!     Purpose.
!     --------
!        To extract the spectral fields for a specific zonal wavenumber
!        and put them in an order suitable for the inverse Legendre           .
!        tranforms.The ordering is from NSMAX to KM for better conditioning.
!        Elements 1,2 and NLCM(KM)+1 are zeroed in preparation for computing
!        u,v and derivatives in spectral space.

!**   Interface.
!     ----------
!        *CALL* *PRFI1B(...)*

!        Explicit arguments :  KM     - zonal wavenumber
!        ------------------    PIA    - spectral components for transform
!                              PSPEC  - spectral array
!                              KFIELDS  - number of fields

!        Implicit arguments :  None.
!        --------------------

!     Method.
!     -------

!     Externals.   None.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Mats Hamrud and Philippe Courtier  *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-02-01 From PRFI1B in IFS CY22R1
!        M.Hamrud      01-Oct-2003 CY28 Cleaning

!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)   :: KFIELDS
REAL(KIND=JPRB)   ,INTENT(IN)   :: PSPEC(:,:)
REAL(KIND=JPRB)   ,INTENT(OUT)  :: PFFA(:,:,:)
INTEGER(KIND=JPIM),INTENT(IN),OPTIONAL :: KFLDPTR(:)

INTEGER(KIND=JPIM) :: II, INM, IR, J, JFLD, ILCM, IOFF,IFLD
INTEGER(KIND=JPIM) :: KM, KMLOC
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

!*       1.    EXTRACT FIELDS FROM SPECTRAL ARRAYS.
!              --------------------------------------------------

IF (LHOOK) CALL DR_HOOK('EPRFI1B_MOD:EPRFI1B',0,ZHOOK_HANDLE)

PFFA = 0._JPRB

DO KMLOC = 1, D%NUMP

  KM   = D%MYMS(KMLOC)
  ILCM = DALD%NCPL2M(KM)
  IOFF = DALD%NESM0(KM)

  IF(PRESENT(KFLDPTR)) THEN
    DO JFLD=1,KFIELDS
      IR = 2*(JFLD-1)+1
      II = IR+1
      IFLD = KFLDPTR(JFLD)
      DO J=1,ILCM,2
        INM = IOFF+(J-1)*2
        PFFA(J  ,IR,KMLOC) = PSPEC(IFLD,INM  )
        PFFA(J+1,IR,KMLOC) = PSPEC(IFLD,INM+1)
        PFFA(J  ,II,KMLOC) = PSPEC(IFLD,INM+2)
        PFFA(J+1,II,KMLOC) = PSPEC(IFLD,INM+3)
      ENDDO
    ENDDO
  
  ELSE
    DO J=1,ILCM,2
      INM = IOFF+(J-1)*2
      DO JFLD=1,KFIELDS
        IR = 2*(JFLD-1)+1
        II = IR+1
        PFFA(J  ,IR,KMLOC) = PSPEC(JFLD,INM  )
        PFFA(J+1,IR,KMLOC) = PSPEC(JFLD,INM+1)
        PFFA(J  ,II,KMLOC) = PSPEC(JFLD,INM+2)
        PFFA(J+1,II,KMLOC) = PSPEC(JFLD,INM+3)
      ENDDO
    ENDDO
  
  ENDIF

ENDDO

#ifdef UNDEF
!$acc update device (PFFA)

WRITE (*, *) __FILE__, ':', __LINE__ 

!$acc serial present (PFFA) present (D_NUMP, DALD_NCPL2M, DALD_NESM0, D_MYMS) 
DO KMLOC = 1, D_NUMP

  KM   = D_MYMS(KMLOC)
  ILCM = DALD_NCPL2M(KM)
  IOFF = DALD_NESM0(KM)

print *, KMLOC, ILCM, IOFF, KM

    DO J=1,ILCM,2
      INM = IOFF+(J-1)*2
      DO JFLD=1,KFIELDS
        IR = 2*(JFLD-1)+1
        II = IR+1
PRINT *, J, IR, II, KMLOC, PFFA(J  ,IR,KMLOC), PFFA(J+1,IR,KMLOC), PFFA(J  ,II,KMLOC), PFFA(J+1,II,KMLOC)
      ENDDO
    ENDDO
  
ENDDO
!$acc end serial

#endif


IF (LHOOK) CALL DR_HOOK('EPRFI1B_MOD:EPRFI1B',1,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

END SUBROUTINE EPRFI1B
END MODULE EPRFI1B_MOD
