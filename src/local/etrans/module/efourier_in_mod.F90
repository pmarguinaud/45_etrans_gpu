MODULE EFOURIER_IN_MOD
CONTAINS
SUBROUTINE EFOURIER_IN(PREEL,KFIELDS)

!**** *FOURIER_IN* - Copy fourier data from buffer to local array

!     Purpose.
!     --------
!        Routine for copying fourier data from buffer to local array

!**   Interface.
!     ----------
!     CALL FOURIER_IN(...)

!     Explicit arguments :  PREEL - local fourier/GP array
!     --------------------  KFIELDS - number of fields
!
!     Externals.  None.
!     ----------

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 2000-04-01

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRBT

USE TPM_DISTR       ,ONLY : D, MYSETW, MYPROC, NPROC, D_NSTAGTF,D_MSTABF,D_NSTAGT0B,D_NPNTGTB0,D_NPROCM,D_NPTRLS
USE TPM_TRANS       ,ONLY : FOUBUF
USE TPM_GEOMETRY    ,ONLY : G, G_NMEN,G_NMEN_MAX
!

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN) :: KFIELDS

INTEGER(KIND=JPIM) :: KGL

REAL(KIND=JPRBT), INTENT(OUT) :: PREEL(:,:)

INTEGER(KIND=JPIM) :: JM,JF,IGLG,IPROC,IR,II,ISTA
INTEGER(KIND=JPIM) :: IBEG,IEND,IINC

!     ------------------------------------------------------------------

IF(MYPROC > NPROC/2)THEN
  IBEG=1
  IEND=D%NDGL_FS
  IINC=1
ELSE
  IBEG=D%NDGL_FS
  IEND=1
  IINC=-1
ENDIF

!$acc data &
!$acc& copyin(D_NPTRLS,D_NSTAGTF,D_MSTABF,D_NSTAGT0B,D_NPNTGTB0,G_NMEN,G_NMEN_MAX,D_NPROCM) &
!$acc& present(PREEL,FOUBUF)

!$acc parallel loop collapse(3) private(IGLG,IPROC,ISTA)
DO KGL=IBEG,IEND,IINC
   DO JM=0,G_NMEN_MAX      
      DO JF=1,KFIELDS     
         IGLG = D_NPTRLS(MYSETW)+KGL-1
         if ( JM .le. G_NMEN(IGLG)) then
            IPROC = D_NPROCM(JM)
            ISTA  = (D_NSTAGT0B(D_MSTABF(IPROC))+D_NPNTGTB0(JM,KGL))*2*KFIELDS
            PREEL(2*JF-1,2*JM+1+D_NSTAGTF(KGL)) = FOUBUF(ISTA+2*JF-1)
            PREEL(2*JF,  2*JM+1+D_NSTAGTF(KGL)) = FOUBUF(ISTA+2*JF  )
         end if
      ENDDO
   ENDDO
END DO

!$acc end data

END SUBROUTINE EFOURIER_IN
END MODULE EFOURIER_IN_MOD

