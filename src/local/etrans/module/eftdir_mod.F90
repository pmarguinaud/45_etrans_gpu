MODULE EFTDIR_MOD
CONTAINS
SUBROUTINE EFTDIR(PREEL, KFIELDS)


!**** *FTDIR - Direct Fourier transform

!     Purpose. Routine for Grid-point to Fourier transform
!     --------

!**   Interface.
!     ----------
!        CALL FTDIR(..)

!        Explicit arguments :  PREEL   - Fourier/grid-point array
!        --------------------  KFIELDS - number of fields

!     Method.
!     -------

!     Externals.  FFT992 - FFT routine
!     ----------
!

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-03-03
!        G. Radnoti 01-04-24 2D model (NLOEN=1)
!        D. Degrauwe  (Feb 2012): Alternative extension zone (E')
!        G. Mozdzynski (Oct 2014): support for FFTW transforms
!        G. Mozdzynski (Jun 2015): Support alternative FFTs to FFTW

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM, JPIB, JPRBT, JPRB

USE TPM_DISTR       ,ONLY : D, MYSETW, MYPROC, NPROC,D_NSTAGTF,D_NPTRLS
USE TPM_GEOMETRY    ,ONLY : G,G_NMEN,G_NMEN_MAX,G_NLOEN,G_NLOEN_MAX
USE TPM_FFT         ,ONLY : T, TB
#ifdef WITH_FFTW
USE TPM_FFTW        ,ONLY : TW, EXEC_FFTW
#endif
USE TPM_FFTC        ,ONLY : CREATE_PLAN_FFT
USE TPM_DIM         ,ONLY : R,R_NNOEXTZL
USE CUDA_DEVICE_MOD
!

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)  :: KFIELDS
REAL(KIND=JPRB), INTENT(INOUT) :: PREEL(:,:)
INTEGER(KIND=JPIM)  :: KGL

INTEGER(KIND=JPIM) :: IGLG,IST,ILEN,IJUMP,JJ,JF,IST1
INTEGER(KIND=JPIM) :: IOFF,IRLEN,ICLEN, ITYPE
INTEGER(KIND=JPIM) :: IPLAN_R2C
INTEGER(KIND=JPIM) :: JMAX
REAL(KIND=JPRBT)   :: SCAL

INTEGER(KIND=JPIM) :: IBEG,IEND,IINC,ISCAL
INTEGER(KIND=JPIM) :: OFFSET_VAR, IUNIT, ISIZE
integer :: istat

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

OFFSET_VAR=D_NPTRLS(MYSETW)

BLOCK
INTEGER :: II, JJ
!$acc serial
DO JJ = 1, 10
DO II = 1, 23
  PRINT *, II, JJ, PREEL (1, II + 23*(JJ-1))
ENDDO
ENDDO
!$acc end serial
ENDBLOCK

!istat = cuda_Synchronize()
DO KGL=IBEG,IEND,IINC

  ITYPE=-1
  IGLG = D_NPTRLS(MYSETW)+KGL-1

  IOFF=D_NSTAGTF(KGL)+1

  CALL CREATE_PLAN_FFT(IPLAN_R2C,-1,KN=G_NLOEN(IGLG),KLOT=SIZE (PREEL, 1))
  !$acc host_data use_device(PREEL)
  CALL EXECUTE_PLAN_FFTC(IPLAN_R2C,-1,PREEL(1, IOFF))
  !$acc end host_data
END DO

WRITE (*, *) __FILE__, ':', __LINE__ 

BLOCK
INTEGER :: II, JJ
!$acc serial
DO JJ = 1, 10
DO II = 1, 23
  PRINT *, II, JJ, PREEL (1, II + 23*(JJ-1))
ENDDO
ENDDO
!$acc end serial
ENDBLOCK

istat = cuda_Synchronize()

!$acc data &
!$acc& copy(D,D_NSTAGTF,D_NPTRLS,G_NMEN,G_NMEN_MAX,G_NLOEN,G_NLOEN_MAX,R_NNOEXTZL) &
!$acc& present(PREEL)

!$acc parallel loop collapse(3) private(JMAX,KGL,IOFF,SCAL)
DO IGLG=IBEG+OFFSET_VAR-1,IEND+OFFSET_VAR-1,IINC
   DO JJ=1, G_NLOEN_MAX+2
      DO JF=1,KFIELDS
         JMAX = G_NLOEN(IGLG)
         if (JJ .le. JMAX) then
           KGL=IGLG-OFFSET_VAR+1
           IOFF=D_NSTAGTF(KGL)+1
           SCAL = 1._JPRBT/REAL(G_NLOEN(IGLG),JPRBT)
           PREEL(JF,IOFF+JJ-1)= SCAL * PREEL(JF, IOFF+JJ-1)
         end if
      ENDDO
   ENDDO
END DO

!$acc parallel loop 
DO KGL=IBEG,IEND,IINC

   IGLG = D_NPTRLS(MYSETW)+KGL-1
   IST  = 2*(G_NMEN(IGLG)+1)+1
   ILEN = G_NLOEN(IGLG)+R_NNOEXTZL+3-IST
   
   IST1=1
   IF (G_NLOEN(IGLG)==1) IST1=0

   !$acc loop collapse(2)
   DO JJ=IST1, ILEN
      DO JF=1,KFIELDS
        PREEL(JF,IST+D_NSTAGTF(KGL)+JJ-1) = 0.0_JPRBT
      ENDDO
   ENDDO
END DO

!$acc end data

!     ------------------------------------------------------------------

END SUBROUTINE EFTDIR
END MODULE EFTDIR_MOD
