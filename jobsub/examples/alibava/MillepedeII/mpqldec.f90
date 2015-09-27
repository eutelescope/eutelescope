!> \file
!! QL decompostion.
!!
!! \author Claus Kleinwort, DESY, 2015 (Claus.Kleinwort@desy.de)
!!
!! \copyright
!! Copyright (c) 2015 Deutsches Elektronen-Synchroton,
!! Member of the Helmholtz Association, (DESY), HAMBURG, GERMANY \n\n
!! This library is free software; you can redistribute it and/or modify
!! it under the terms of the GNU Library General Public License as
!! published by the Free Software Foundation; either version 2 of the
!! License, or (at your option) any later version. \n\n
!! This library is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU Library General Public License for more details. \n\n
!! You should have received a copy of the GNU Library General Public
!! License along with this program (see the file COPYING.LIB for more
!! details); if not, write to the Free Software Foundation, Inc.,
!! 675 Mass Ave, Cambridge, MA 02139, USA.
!!
!! QL decomposition of constraints matrix by Householder transformations
!! for solution by elimination.
!!

!> QL data.
MODULE mpqldec
    USE mpdef
    IMPLICIT NONE

    INTEGER(mpi) :: npar   !< number of parameters
    INTEGER(mpi) :: ncon   !< number of constraints
    REAL(mpd), DIMENSION(:), ALLOCATABLE :: matV !< unit normals (v_i) of Householder reflectors
    REAL(mpd), DIMENSION(:), ALLOCATABLE :: matL !< lower diagonal matrix L
END MODULE mpqldec

!> Initialize QL decomposition.
!!
!! \param [in]     n  number of rows (parameters)
!! \param [in]     m  number of columns (constraints)
!!
SUBROUTINE qlini(n,m)
    USE mpqldec
    USE mpdalc

    IMPLICIT NONE
    INTEGER(mpl) :: length

    INTEGER(mpi), INTENT(IN)          :: n
    INTEGER(mpi), INTENT(IN)          :: m

    npar=n
    ncon=m
    ! allocate 
    length=npar*ncon
    CALL mpalloc(matV,length,'QLDEC: V')
    length=ncon*ncon
    CALL mpalloc(matL,length,'QLDEC: L')
END SUBROUTINE qlini

!                                                 141217 C. Kleinwort, DESY-FH1
!> QL decomposition.
!!
!! QL decomposition with Householder transformations.
!! Decompose N-By-M matrix A into orthogonal N-by-N matrix Q and a
!! N-by-M matrix containing zeros except for a lower triangular
!! M-by-M matrix L (at the bottom):
!!
!!              | 0 |
!!      A = Q * |   |
!!              | L |
!!
!! The decomposition is stored in a N-by-M matrix matV containing the unit
!! normal vectors v_i of the hyperplanes (Householder reflectors) defining Q.
!! The lower triangular matrix L is stored in the M-by-M matrix matL.
!!
!! \param [in]     a  Npar-by-Ncon matrix
!!
SUBROUTINE qldec(a)
    USE mpqldec
    USE mpdalc

    ! cost[dot ops] ~= Npar*Ncon*Ncon

    IMPLICIT NONE
    INTEGER(mpi) :: i
    INTEGER(mpl) :: ioff1
    INTEGER(mpl) :: ioff2
    INTEGER(mpl) :: ioff3
    INTEGER(mpi) :: k
    INTEGER(mpi) :: kn
    INTEGER(mpl) :: length
    REAL(mpd) :: nrm
    REAL(mpd) :: sp

    REAL(mpd), INTENT(IN)             :: a(*)

    REAL(mpd)                         :: v(npar)

    ! prepare 
    length=npar*ncon
    matV=a(1:length)
    matL=0.0_mpd

    ! Householder procedure
    DO k=ncon,1,-1
        kn=npar+k-ncon
        ! column offset
        ioff1=(k-1)*npar
        ! get column
        v(1:kn)=matV(ioff1+1:ioff1+kn)
        nrm = SQRT(dot_product(v(1:kn),v(1:kn)))
        IF (nrm == 0.0_mpd) CYCLE
        !
        IF (v(kn) >= 0.0_mpd) THEN
            v(kn)=v(kn)+nrm
        ELSE
            v(kn)=v(kn)-nrm
        END IF
        ! create normal vector
        nrm = SQRT(dot_product(v(1:kn),v(1:kn)))
        v(1:kn)=v(1:kn)/nrm
        ! transformation
        ioff2=0
        DO i=1,k
            sp=dot_product(v(1:kn),matV(ioff2+1:ioff2+kn))
            matV(ioff2+1:ioff2+kn)=matV(ioff2+1:ioff2+kn)-2.0_mpd*v(1:kn)*sp
            ioff2=ioff2+npar
        END DO
        ! store column of L
        ioff3=(k-1)*ncon
        matL(ioff3+k:ioff3+ncon)=matV(ioff1+kn:ioff1+npar)
        ! store normal vector
        matV(ioff1+1:ioff1+kn)=v(1:kn)
        matV(ioff1+kn+1:ioff1+npar)=0.0_mpd
    END DO

END SUBROUTINE qldec


!> Multiply left by Q(t).
!!
!! Multiply left by Q(t) from QL decomposition.
!!
!! \param [in,out] x    Npar-by-M matrix, overwritten with Q*X (t=false) or Q^t*X (t=true)
!! \param [in]     m    number of columns
!! \param [in]     t    use transposed of Q
!!
SUBROUTINE qlmlq(x,m,t)
    USE mpqldec

    ! cost[dot ops] ~= N*M*Nhr

    IMPLICIT NONE
    INTEGER(mpi) :: i
    INTEGER(mpl) :: ioff1
    INTEGER(mpl) :: ioff2
    INTEGER(mpi) :: j
    INTEGER(mpi) :: k
    INTEGER(mpi) :: kn
    REAL(mpd) :: sp

    REAL(mpd), INTENT(IN OUT)         :: x(*)
    INTEGER(mpi), INTENT(IN)          :: m
    LOGICAL, INTENT(IN)               :: t

    DO j=1,ncon
        k=j
        IF (t) k=ncon+1-j
        kn=npar+k-ncon
        ! column offset
        ioff1=(k-1)*npar
        ! transformation
        ioff2=0
        DO i=1,m
            sp=dot_product(matV(ioff1+1:ioff1+kn),x(ioff2+1:ioff2+kn))
            x(ioff2+1:ioff2+kn)=x(ioff2+1:ioff2+kn)-2.0_mpd*matV(ioff1+1:ioff1+kn)*sp
            ioff2=ioff2+npar
        END DO
    END DO

END SUBROUTINE qlmlq


!> Multiply right by Q(t).
!!
!! Multiply right by Q(t) from QL decomposition.
!!
!! \param [in,out] x    M-by-Npar matrix, overwritten with X*Q (t=false) or X*Q^t (t=true)
!! \param [in]     m    number of rows
!! \param [in]     t    use transposed of Q
!!
SUBROUTINE qlmrq(x,m,t)
    USE mpqldec

    ! cost[dot ops] ~= N*M*Nhr

    IMPLICIT NONE
    INTEGER(mpi) :: i
    INTEGER(mpl) :: ioff1
    INTEGER(mpl) :: iend
    INTEGER(mpi) :: j
    INTEGER(mpi) :: k
    INTEGER(mpi) :: kn
    REAL(mpd) :: sp

    REAL(mpd), INTENT(IN OUT)         :: x(*)
    INTEGER(mpi), INTENT(IN)          :: m
    LOGICAL, INTENT(IN)               :: t

    DO j=1,ncon
        k=j
        IF (.not.t) k=ncon+1-j
        kn=npar+k-ncon
        ! column offset
        ioff1=(k-1)*npar
        ! transformation
        iend=m*kn
        DO i=1,npar
            sp=dot_product(matV(ioff1+1:ioff1+kn),x(i:iend:m))
            x(i:iend:m)=x(i:iend:m)-2.0_mpd*matV(ioff1+1:ioff1+kn)*sp
        END DO
    END DO

END SUBROUTINE qlmrq


!> Similarity transformation by Q(t).
!!
!! Similarity transformation by Q from QL decomposition.
!!
!! \param [in,out] x    Npar-by-Npar matrix, overwritten with Q*X*Q^t (t=false) or Q^t*X*Q (t=true)
!! \param [in]     t    use transposed of Q
!!
SUBROUTINE qlsmq(x,t)
    USE mpqldec

    ! cost[dot ops] ~= N*N*Nhr

    IMPLICIT NONE
    INTEGER(mpi) :: i
    INTEGER(mpl) :: ioff1
    INTEGER(mpl) :: ioff2
    INTEGER(mpl) :: iend
    INTEGER(mpi) :: j
    INTEGER(mpi) :: k
    INTEGER(mpi) :: kn
    REAL(mpd) :: sp

    REAL(mpd), INTENT(IN OUT)         :: x(*)
    LOGICAL, INTENT(IN)               :: t

    DO j=1,ncon
        k=j
        IF (t) k=ncon+1-j
        kn=npar+k-ncon
        ! column offset
        ioff1=(k-1)*npar
        ! transformation
        iend=npar*kn
        DO i=1,npar
            sp=dot_product(matV(ioff1+1:ioff1+kn),x(i:iend:npar))
            x(i:iend:npar)=x(i:iend:npar)-2.0_mpd*matV(ioff1+1:ioff1+kn)*sp
        END DO
        ioff2=0
        DO i=1,npar
            sp=dot_product(matV(ioff1+1:ioff1+kn),x(ioff2+1:ioff2+kn))
            x(ioff2+1:ioff2+kn)=x(ioff2+1:ioff2+kn)-2.0_mpd*matV(ioff1+1:ioff1+kn)*sp
            ioff2=ioff2+npar
        END DO
    END DO

END SUBROUTINE qlsmq


!> Similarity transformation by Q(t).
!!
!! Similarity transformation for symmetric matrix by Q from QL decomposition.
!!
!! \param [in]     aprod    external procedure to calculate A*v
!! \param [in,out] A        symmetric Npar-by-Npar matrix A in symmetric storage mode
!!                          (V(1) = V11, V(2) = V12, V(3) = V22, V(4) = V13, ...),
!!                          overwritten with Q*A*Q^t (t=false) or Q^t*A*Q (t=true)
!! \param [in]     t        use transposed of Q
!!
SUBROUTINE qlssq(aprod,A,t)
    USE mpqldec
    USE mpdalc

    ! cost[dot ops] ~= N*N*Nhr

    IMPLICIT NONE
    INTEGER(mpi) :: i
    INTEGER(mpl) :: ioff1
    INTEGER(mpl) :: ioff2
    INTEGER(mpi) :: j
    INTEGER(mpi) :: k
    INTEGER(mpi) :: kn
    INTEGER(mpi) :: l
    INTEGER(mpl) :: length
    REAL(mpd) :: vtAv
    REAL(mpd), DIMENSION(:), ALLOCATABLE :: Av

    REAL(mpd), INTENT(IN OUT)         :: A(*)
    LOGICAL, INTENT(IN)               :: t

    INTERFACE
        SUBROUTINE aprod(n,x,y) ! y=A*x
            USE mpdef
            INTEGER(mpi), INTENT(in) :: n
            REAL(mpd), INTENT(IN)    :: x(n)
            REAL(mpd), INTENT(OUT)   :: y(n)
        END SUBROUTINE aprod
    END INTERFACE

    length=npar
    CALL mpalloc(Av,length,'qlssq: A*v')

    DO j=1,ncon
        k=j
        IF (t) k=ncon+1-j
        kn=npar+k-ncon
        ! column offset
        ioff1=(k-1)*npar
        ! A*v
        CALL aprod(npar,matV(ioff1+1:ioff1+npar),Av(1:npar))
        ! transformation
        ! diagonal block
        ! v^t*A*v
        vtAv=dot_product(matV(ioff1+1:ioff1+kn),Av(1:kn))
        ! update 
        ioff2=0
        DO i=1,kn
            ! correct with  2*(2v*vtAv*v^t - Av*v^t - (Av*v^t)^t)
            DO l=1,i
                ioff2=ioff2+1
                A(ioff2)=A(ioff2)+2.0_mpd*((2.0_mpd*matV(ioff1+i)*vtAv-Av(i))*matV(ioff1+l)-Av(l)*matV(ioff1+i))
            END DO
        END DO
        ! off diagonal block
        DO i=kn+1,npar
            ! correct with -2Av*v^t
            A(ioff2+1:ioff2+kn)=A(ioff2+1:ioff2+kn)-2.0_mpd*matV(ioff1+1:ioff1+kn)*Av(i)
            ioff2=ioff2+i
        END DO
    END DO

    CALL mpdealloc(Av)

END SUBROUTINE qlssq


!> Partial similarity transformation by Q(t).
!!
!! Partial similarity transformation for symmetric matrix by Q from QL decomposition.
!! Calculate corrections to band part of matrix.
!!
!! \param [in]     aprod    external procedure to calculate A*v
!! \param [in,out] B        band part of symmetric Npar-by-Npar matrix A in symmetric storage mode,
!!                          overwritten with band part of Q^t*A*Q (t=false) or Q^t*A*Q (t=true)
!! \param [in]     m        band width (including diagonal)
!! \param [in]     t        use transposed of Q
!!
SUBROUTINE qlpssq(aprod,B,m,t)
    USE mpqldec
    USE mpdalc

    ! cost[dot ops] ~= N*N*Nhr

    IMPLICIT NONE
    INTEGER(mpi) :: i
    INTEGER(mpl) :: ioff1
    INTEGER(mpl) :: ioff2
    INTEGER(mpi) :: j
    INTEGER(mpi) :: j2
    INTEGER(mpi) :: k
    INTEGER(mpi) :: k2
    INTEGER(mpi) :: kn
    INTEGER(mpi) :: l
    INTEGER(mpl) :: length
    INTEGER(mpi) :: mbnd
    REAL(mpd) :: vtAv
    REAL(mpd) :: vtAvp
    REAL(mpd) :: vtvp
    REAL(mpd), DIMENSION(:), ALLOCATABLE :: Av   ! A*v 

    REAL(mpd), INTENT(IN OUT)         :: B(*)
    INTEGER(mpi), INTENT(IN)          :: m
    LOGICAL, INTENT(IN)               :: t

    INTERFACE
        SUBROUTINE aprod(n,x,y) ! y=A*x
            USE mpdef
            INTEGER(mpi), INTENT(in) :: n
            REAL(mpd), INTENT(IN)    :: x(n)
            REAL(mpd), INTENT(OUT)   :: y(n)
        END SUBROUTINE aprod
    END INTERFACE

    length=npar
    length=npar*ncon
    CALL mpalloc(Av,length,'qlpssq: Av')

    mbnd=max(0,m-1) ! band width without diagonal
    ! A*V
    ioff1=0
    DO i=1,ncon
        CALL aprod(npar,matV(ioff1+1:ioff1+npar),Av(ioff1+1:ioff1+npar))
        ioff1=ioff1+npar
    END DO

    DO j=1,ncon
        k=j
        IF (t) k=ncon+1-j
        kn=npar+k-ncon
        ! column offset
        ioff1=(k-1)*npar
        ! transformation (diagonal block)
        ! diagonal block
        ! v^t*A*v
        vtAv=dot_product(matV(ioff1+1:ioff1+kn),Av(ioff1+1:ioff1+kn))
        ! update
        ioff2=0
        DO i=1,kn
            ! correct with  2*(2v*vtAv*v^t - Av*v^t - (Av*v^t)^t)
            DO l=max(1,i-mbnd),i
                ioff2=ioff2+1
                B(ioff2)=B(ioff2)+2.0_mpd*((2.0_mpd*matV(ioff1+i)*vtAv-Av(ioff1+i))*matV(ioff1+l)-Av(ioff1+l)*matV(ioff1+i))
            END DO
        END DO
        ! off diagonal block
        DO i=kn+1,npar
            ! correct with -2Av*v^t
            DO l=max(1,i-mbnd),i
                ioff2=ioff2+1
                B(ioff2)=B(ioff2)-2.0_mpd*Av(ioff1+i)*matV(ioff1+l)
            END DO
        END DO
        ! correct A*v for the remainung v
        DO j2=j+1,ncon
            k2=j2
            IF (t) k2=ncon+1-j2
            ioff2=(k2-1)*npar
            vtvp=dot_product(matV(ioff1+1:ioff1+npar),matV(ioff2+1:ioff2+npar)) ! v^t*v'
            vtAvp=dot_product(matV(ioff1+1:ioff1+npar),Av(ioff2+1:ioff2+npar)) ! v^t*(A*v')
            DO i=1,kn 
                Av(ioff2+i)=Av(ioff2+i)+2.0_mpd*((2.0_mpd*matV(ioff1+i)*vtAv-Av(ioff1+i))*vtvp-matV(ioff1+i)*vtAvp)
            END DO
            DO i=kn+1,npar 
                Av(ioff2+i)=Av(ioff2+i)-2.0_mpd*Av(ioff1+i)*vtvp
            END DO
        END DO

    END DO

    CALL mpdealloc(Av)

END SUBROUTINE qlpssq


!> Get eigenvalues.
!!
!! Get smallest and largest eingenvalue of L.
!!
!! \param [out]    emin  smallest eigenvalue
!! \param [out]    emax  largest eigenvalue
!!
SUBROUTINE qlgete(emin,emax)
    USE mpqldec

    IMPLICIT NONE
    INTEGER(mpi) :: i
    INTEGER(mpl) :: idiag

    REAL(mpd), INTENT(OUT)         :: emin
    REAL(mpd), INTENT(OUT)         :: emax

    idiag=1
    emax=matL(1)
    emin=emax
    DO i=2,ncon
        idiag=idiag+ncon+1
        IF (ABS(emax) < ABS(matL(idiag))) emax=matL(idiag)
        IF (ABS(emin) > ABS(matL(idiag))) emin=matL(idiag)
    END DO

END SUBROUTINE qlgete


!> Backward substitution.
!!
!! Get y from L^t*y=d.
!!
!! \param [in]     d  Ncon vector, resdiduals
!! \param [out]    y  Ncon vector, solution
!!
SUBROUTINE qlbsub(d,y)
    USE mpqldec

    IMPLICIT NONE
    INTEGER(mpl) :: idiag
    INTEGER(mpi) :: k

    REAL(mpd), INTENT(IN)         :: d(ncon)
    REAL(mpd), INTENT(OUT)        :: y(ncon)

    ! solve L*y=d by forward substitution
    idiag=ncon*ncon
    DO k=ncon,1,-1
        y(k)=(d(k)-dot_product(matL(idiag+1:idiag+ncon-k),y(k+1:ncon)))/matL(idiag)
        idiag=idiag-ncon-1
    END DO

END SUBROUTINE qlbsub
