MODULE deriv_MQMP4

CONTAINS
!============================================================================================
! 8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!============================================================================================

SUBROUTINE Rec_derivMQMP( matcoord, noderiv, matelem, Uaprox, nderiv, nelem, GradUaprox )

IMPLICIT NONE
INTEGER					:: nglob,nelem, nderiv
REAL(8),DIMENSION(:,:), INTENT(in)	:: matcoord!, matelem
INTEGER, DIMENSION(:,:), INTENT(in)	:: matelem
REAL(8),DIMENSION(:,:), INTENT(out)	:: GradUaprox
REAL(8),DIMENSION(:), INTENT(in)	:: Uaprox
INTEGER,DIMENSION(:), intent(in)	:: noderiv

REAL(8), ALLOCATABLE, DIMENSION(:)	:: funcaoP, vetx

INTEGER,PARAMETER	:: n=6, nn=15
REAL(8),DIMENSION(:,:)	:: matC_loc(nn,3), G(n,n)
REAL(8)			:: Rg, f(n), R(2)=0.0d0, m=2.2d0, dm=.2d0
INTEGER			:: i, j, Ng, nloc, cont, cont1, no_interp(nn)

cont1=0
matC_loc(nn,3)=0.d0
G(n,n)=0


DO i=1, nderiv


	Ng=noderiv(i)
!	CALL Raio_Domloc(matcoord,matelem,nelem,Rg,Ng,m)
!	CALL coord_Domloc(matcoord,nglob,Rg,matC_loc,nloc,Ng)
!
!	do while(nloc < 6)
!		m = m+dm
!		CALL Raio_Domloc(matcoord,matelem,nelem,Rg,Ng,m)
!		CALL coord_Domloc(matcoord,nglob,Rg,matC_loc,nloc,Ng)
!	end do
!
!	m=2.2d0
	call nos_envolta( matelem, nelem, no_interp, nn, Ng, nloc )
	call coord_loc( matcoord, no_interp, nloc, matC_loc )

	ALLOCATE( funcaoP(nloc) )
	CALL funcaoPeso( matC_loc, funcaoP, nloc )
	funcaoP = 1.d0
	CALL sistema_EqN( matC_loc, funcaoP, Uaprox, G, f, nloc )
	CALL ReduzMat( G, n )
	CALl RetroMat( G, f, n )

!	cont1=cont1+1
!	GradUaprox(cont1,1)=Ng
!	GradUaprox(cont1,2)=f(2)+f(4)*matC_loc(1,3)+f(5)*matC_loc(1,2)*2
!	GradUaprox(cont1,3)=f(3)+f(4)*matC_loc(1,2)+f(6)*matC_loc(1,3)*2
	GradUaprox(i,1) = Ng
	GradUaprox(i,2) = f(2) + f(4)*matC_loc(1,3) + f(5)*matC_loc(1,2)*2
	GradUaprox(i,3) = f(3) + f(4)*matC_loc(1,2) + f(6)*matC_loc(1,3)*2
	DEALLOCATE(funcaoP)

END DO


END SUBROUTINE Rec_derivMQMP

!============================================================================================
! 8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!============================================================================================

SUBROUTINE Rec_derivMQMP_3D( matcoord, noderiv, matelem, Uaprox, nderiv, nelem, dfdx, dfdy, dfdz, nconstb )

IMPLICIT NONE
INTEGER					:: nglob,nelem, nderiv, nconstb
REAL(8), DIMENSION(:,:), INTENT(in)	:: matcoord!, matelem
INTEGER, DIMENSION(:,:), INTENT(in)	:: matelem
REAL(8), DIMENSION(:), INTENT(in)	:: Uaprox
REAL(8), DIMENSION(:), INTENT(out)	:: dfdx, dfdy, dfdz
INTEGER, DIMENSION(:), intent(in)	:: noderiv

REAL(8), ALLOCATABLE, DIMENSION(:)	:: funcaoP, vetx

INTEGER,PARAMETER	:: nn=305
REAL(8),DIMENSION(:,:)	:: matC_loc(nn,3), G(nconstb,nconstb)
REAL(8)			:: Rg, f(nconstb), R(2)=0.0d0, m=2.2d0, dm=.2d0
INTEGER			:: i, j, Ng, nloc, cont, cont1, no_interp(nn)

cont1=0
matC_loc(nn,3)=0.d0
G=0.d0


DO i=1, nderiv

	Ng=noderiv(i)
	call nos_envolta_3D( matelem, nelem, no_interp, nn, Ng, nloc, nconstb )
	call coord_loc_3D( matcoord, no_interp, nloc, matC_loc )

	ALLOCATE( funcaoP(nloc) )
	CALL funcaoPeso_3D( matC_loc, funcaoP, nloc )

  funcaoP = 1.d0
  if( nconstb == 8)then
    CALL EqN_3D_base_l( matC_loc, no_interp, funcaoP, Uaprox, G, f, nloc )
  else
    CALL EqN_3D_base_q( matC_loc, no_interp, funcaoP, Uaprox, G, f, nloc )
  endif
	CALL ReduzMat( G, nconstb )
	CALl RetroMat( G, f, nconstb )

	if( nconstb == 8)then
		dfdx(i) = f(2) + f(5)*matC_loc(1,2) + f(6)*matC_loc(1,3) + f(8)*matC_loc(1,2)*matC_loc(1,3)
		dfdy(i) = f(3) + f(5)*matC_loc(1,1) + f(7)*matC_loc(1,3) + f(8)*matC_loc(1,1)*matC_loc(1,3)
		dfdz(i) = f(4) + f(6)*matC_loc(1,1) + f(7)*matC_loc(1,2) + f(8)*matC_loc(1,1)*matC_loc(1,2)
	else
		dfdx(i) = f(2) + f(5)*matC_loc(1,2) + f(6)*matC_loc(1,3) + f(8)*matC_loc(1,2)*matC_loc(1,3) &
									 + f(9)*matC_loc(1,1)*2.d0
		dfdy(i) = f(3) + f(5)*matC_loc(1,1) + f(7)*matC_loc(1,3) + f(8)*matC_loc(1,1)*matC_loc(1,3) &
									 + f(10)*matC_loc(1,2)*2.d0
		dfdz(i) = f(4) + f(6)*matC_loc(1,1) + f(7)*matC_loc(1,2) + f(8)*matC_loc(1,1)*matC_loc(1,2) &
									 + f(11)*matC_loc(1,3)*2.d0
	end if

	DEALLOCATE(funcaoP)

END DO


END SUBROUTINE Rec_derivMQMP_3D

!============================================================================================
! 8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!============================================================================================

SUBROUTINE nos_envolta( matelem, nelem, no_interp, nn, ng, nloc )
IMPLICIT NONE
INTEGER, INTENT(in)			:: nn, ng, nelem
INTEGER, INTENT(out)			:: nloc
INTEGER, DIMENSION(:), INTENT(out)	:: no_interp
!REAL(8), DIMENSION(:,:), INTENT(in)	:: matelem
INTEGER, DIMENSION(:,:), INTENT(in)	:: matelem

INTEGER:: l, k, n, m, cont, noelem(3)

no_interp = 0
no_interp(1) = ng
cont = 0

do l=1, nelem
	noelem = int(matelem(l,2:4))
	do m =1,3
		if( noelem(m) == ng )then
			cont  = cont + 1
			if( cont > nn )then
				write(*,*)' dimensão do vetor no_interp excedida para a derivada do nó ', ng
			end if
			do k=1, 3
				do n=1, nn
					if( no_interp(n) /= noelem(k) .and. no_interp(n) /= 0 )cycle
					no_interp(n) = noelem(k)
					exit
				end do
			end do
			exit
		end if
	end do
end do

if( cont < 5 )then
	print*,' O número de nós ao redor do nó', ng,' é < 6. Portanto o problema de '
	print*,' Mínimos quadrados e subdeterminado '
	stop
end if
nloc = cont + 1

END SUBROUTINE nos_envolta

!============================================================================================
! 8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!============================================================================================

SUBROUTINE nos_envolta_3D( matelem, nelem, no_interp, nn, ng, nloc, nconstb )
IMPLICIT NONE
INTEGER, INTENT(in)			:: nn, ng, nelem, nconstb
INTEGER, INTENT(out)			:: nloc
INTEGER, DIMENSION(:), INTENT(out)	:: no_interp
!REAL(8), DIMENSION(:,:), INTENT(in)	:: matelem
INTEGER, DIMENSION(:,:), INTENT(in)	:: matelem

INTEGER:: l, k, n, m, cont, noelem(4)

no_interp = 0
no_interp(1) = ng
cont = 0

do l=1, nelem
	noelem = matelem(l,1:4)
	do m =1,4
		if( noelem(m) == ng )then
			cont  = cont + 1
			if( cont > nn )then
				write(*,*)' dimensão do vetor no_interp excedida para a derivada do nó ', ng
			end if
			do k=1, 4
				do n=2, nn
					if( no_interp(n) /= noelem(k) .and. no_interp(n) /= 0 )cycle
					no_interp(n) = noelem(k)
					exit
				end do
			end do
			exit
		end if
	end do
end do

if( cont < nconstb-1 )then
	print*,' O número de nós ao redor do nó', ng,' é',cont,' < ',nconstb-1,'. Portanto o problema de '
	print*,' Mínimos quadrados e subdeterminado '
	nloc = cont + 1
	return
!	stop
end if

nloc = cont + 1

END SUBROUTINE nos_envolta_3D

!============================================================================================
! 8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!============================================================================================
SUBROUTINE coord_loc( matcoord, no_interp, nloc, matC_loc )
IMPLICIT NONE
INTEGER, INTENT(in):: nloc
INTEGER, DIMENSION(:), INTENT(in):: no_interp
REAL(8), DIMENSION(:, :), INTENT(in):: matcoord
REAL(8), DIMENSION(:, :), INTENT(out):: matC_loc

INTEGER:: l

do l=1 , nloc
	matC_loc(l, 1) = no_interp(l)
	matC_loc(l, 2) = matcoord(no_interp(l), 2)
	matC_loc(l, 3) = matcoord(no_interp(l), 3)
end do


END SUBROUTINE coord_loc

!============================================================================================
! 8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!============================================================================================

SUBROUTINE coord_loc_3D( matcoord, no_interp, nloc, matC_loc )
IMPLICIT NONE
INTEGER, INTENT(in):: nloc
INTEGER, DIMENSION(:), INTENT(in):: no_interp
REAL(8), DIMENSION(:, :), INTENT(in):: matcoord
REAL(8), DIMENSION(:, :), INTENT(out):: matC_loc

INTEGER:: l

do l=1 , nloc
	matC_loc(l, 1) = matcoord(no_interp(l), 1)
	matC_loc(l, 2) = matcoord(no_interp(l), 2)
	matC_loc(l, 3) = matcoord(no_interp(l), 3)
end do


END SUBROUTINE coord_loc_3D

!============================================================================================
! 8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!============================================================================================


SUBROUTINE Raio_Domloc(matcoord,matelem,nelem,Rg,Ng,m)

IMPLICIT NONE
INTEGER,INTENT(IN)::Ng,nelem
REAL(8),INTENT(INOUT)::Rg,m
REAL(8),DIMENSION(:,:),INTENT(IN)::matcoord!,matelem
INTEGER, DIMENSION(:,:), INTENT(in)	:: matelem
REAL(8)::x(3),z(3),a(3),det(2),alfa
REAL(8),PARAMETER:: pi = datan(1.d0)
INTEGER::i,j,cont,No(3)

cont=0
DO i=1,nelem
   No(1:3)=INT(matelem(i,2:4))
   IF(Ng == No(1) .or. Ng == No(2) .or. Ng == No(3))THEN
    cont=cont+1
    IF(cont==1)THEN
     x(1:3)=matcoord(No(1:3),2)
     z(1:3)=matcoord(No(1:3),3)
     a(1)=x(2)*z(3)-x(3)*z(2)
     a(2)=x(3)*z(1)-x(1)*z(3)
     a(3)=x(1)*z(2)-x(2)*z(1)
     det(1)=SUM(a(1:3))
    ELSE
     x(1:3)=matcoord(No(1:3),2)
     z(1:3)=matcoord(No(1:3),3)
     a(1)=x(2)*z(3)-x(3)*z(2)
     a(2)=x(3)*z(1)-x(1)*z(3)
     a(3)=x(1)*z(2)-x(2)*z(1)
     det(2)=SUM(a(1:3))
      IF(det(2) < det(1))THEN
        det(1)=det(2)
      END IF

    END IF
   END IF
END DO

det(1) = abs(det(1))/2.d0
alfa = dsqrt( m*pi/det(1) )
Rg=alfa*det(1)

END SUBROUTINE Raio_Domloc

!============================================================================================
! 8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!============================================================================================

SUBROUTINE coord_Domloc(matcoord,nglob,Rg,matC_loc,nloc,Ng)
IMPLICIT NONE
INTEGER,INTENT(IN)::nglob,Ng
INTEGER,INTENT(OUT)::nloc
REAL(8),INTENT(IN)::Rg
REAL(8),INTENT(IN)::matcoord(nglob,4)
REAL(8),DIMENSION(:,:),INTENT(OUT)::matC_loc
REAL(8)::raiog,x(2),z(2)
INTEGER::i

nloc=1
x(1)=matcoord(Ng,2)
z(1)=matcoord(Ng,3)
matC_loc(1,1)=Ng
matC_loc(1,2)=x(1)
matC_loc(1,3)=z(1)
DO i=1,nglob
   IF(i==Ng)CYCLE
   x(2)=matcoord(i,2)
   z(2)=matcoord(i,3)
   raiog=SQRT((x(1)-x(2))**2+(z(1)-z(2))**2)
   IF(Rg >= raiog)THEN
   nloc=nloc+1
   matC_loc(nloc,1)=i
   matC_loc(nloc,2)=x(2)
   matC_loc(nloc,3)=z(2)
   END IF
END DO
END SUBROUTINE coord_Domloc

!============================================================================================
! 8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!============================================================================================

SUBROUTINE funcaoPeso(matC_loc,funcaoP,nloc)
IMPLICIT NONE
INTEGER,INTENT(IN)::nloc
REAL(8),DIMENSION(:),INTENT(OUT)::funcaoP
REAL(8),DIMENSION(:,:),INTENT(IN)::matC_loc
REAL(8)::xmax,xmin,zmax,zmin
INTEGER::j

xmax=matC_loc(1,2)
zmax=matC_loc(1,3)
DO j=2,nloc

  xmin=matC_loc(j,2)
  zmin=matC_loc(j,3)

  IF(xmax < xmin )THEN
    xmax=xmin
  END IF

  IF(zmax < zmin)THEN
    zmax=zmin
  END IF

END DO

DO j=1,nloc
   funcaoP(j) = EXP(-(3.0)**2*( ((matC_loc(1,2)-matC_loc(j,2))/(matC_loc(1,2)-xmax))**2 + &
                              ((matC_loc(1,3)-matC_loc(j,3))/(matC_loc(1,3)-zmax))**2 ))
END DO

END SUBROUTINE funcaoPeso

!============================================================================================
! 8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!============================================================================================

SUBROUTINE funcaoPeso_3D( matC_loc, funcaoP, nloc )
IMPLICIT NONE
INTEGER,INTENT(IN)::nloc
REAL(8),DIMENSION(:),INTENT(OUT)::funcaoP
REAL(8),DIMENSION(:,:),INTENT(IN)::matC_loc
REAL(8):: xmax, xmin, ymax, ymin, zmax, zmin
INTEGER::j

xmax=(matC_loc(1,1)-matC_loc(2,1))
ymax=(matC_loc(1,2)-matC_loc(2,2))
zmax=(matC_loc(1,3)-matC_loc(2,3))

DO j=3,nloc

  xmin=(matC_loc(j,1)-matC_loc(1,1))
  ymin=(matC_loc(j,2)-matC_loc(1,2))
  zmin=(matC_loc(j,3)-matC_loc(1,3))

  IF(xmax < xmin )THEN
    xmax=xmin
  END IF

  IF(ymax < ymin )THEN
    ymax=ymin
  END IF

  IF(zmax < zmin)THEN
    zmax=zmin
  END IF

END DO

!DO j=1,nloc
!   funcaoP(j) = EXP(-(3.0)**2*( ((matC_loc(1,2)-matC_loc(j,2))/(matC_loc(1,2)-xmax))**2 + &
!                              ((matC_loc(1,3)-matC_loc(j,3))/(matC_loc(1,3)-zmax))**2 ))
!END DO

DO j=1,nloc
   funcaoP(j) = EXP(-(3.0)**2*( ((matC_loc(1,1)-matC_loc(j,1))/(matC_loc(1,1)-xmax))**2 + &
	        ((matC_loc(1,2)-matC_loc(j,2))/(matC_loc(1,2)-ymax))**2 + ((matC_loc(1,3)-matC_loc(j,3))/(matC_loc(1,3)-zmax))**2 ))
END DO

END SUBROUTINE funcaoPeso_3D

!============================================================================================
! 8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!============================================================================================

SUBROUTINE sistema_EqN( matC_loc, funcaoP, Uaprox, G, f, nloc )

IMPLICIT NONE
INTEGER,INTENT(IN)::nloc
REAL(8),DIMENSION(:),INTENT(OUT)::f
REAL(8),DIMENSION(:),INTENT(IN)::funcaoP,Uaprox
REAL(8),DIMENSION(:,:),INTENT(OUT)::G
REAL(8),DIMENSION(:,:),INTENT(IN)::matC_loc
REAL(8),DIMENSION(:,:)::A(nloc,6),B(6,nloc)
REAL(8),DIMENSION(:)::f1(nloc)
INTEGER::j

DO j=1,nloc
  f1(j)=Uaprox(INT(matC_loc(j,1)))*funcaoP(j)
  A(j,1)=1
  A(j,2)=matC_loc(j,2)
  A(j,3)=matC_loc(j,3)
  A(j,4)=matC_loc(j,2)*matC_loc(j,3)
  A(j,5)=matC_loc(j,2)**2
  A(j,6)=matC_loc(j,3)**2
END DO
B=TRANSPOSE(A)
DO j=1,nloc
  A(j,1:6)=A(j,1:6)*funcaoP(j)
END DO
f=MATMUL(B,f1)
G=MATMUL(B,A)

END SUBROUTINE sistema_EqN

!============================================================================================
! 8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!============================================================================================
SUBROUTINE EqN_3D_base_l(matC_loc, no_interp, funcaoP, Uaprox,G ,f , nloc )

IMPLICIT NONE
INTEGER, DIMENSION(:), INTENT(IN)	:: no_interp
INTEGER, INTENT(IN)			:: nloc
REAL(8), DIMENSION(:), INTENT(OUT)	:: f
REAL(8), DIMENSION(:), INTENT(IN)	:: funcaoP,Uaprox
REAL(8), DIMENSION(:,:), INTENT(OUT)	:: G
REAL(8), DIMENSION(:,:), INTENT(IN)	:: matC_loc

REAL(8), DIMENSION(:,:)	:: A(nloc,8), B(8,nloc)
REAL(8), DIMENSION(:)	:: f1(nloc)
INTEGER			:: j

DO j=1,nloc
  f1(j) = Uaprox(no_interp(j))*funcaoP(j)
  A(j,1) = 1.d0
  A(j,2) = matC_loc(j,1)
  A(j,3) = matC_loc(j,2)
  A(j,4) = matC_loc(j,3)
  A(j,5) = matC_loc(j,1)*matC_loc(j,2)
  A(j,6) = matC_loc(j,1)*matC_loc(j,3)
  A(j,7) = matC_loc(j,2)*matC_loc(j,3)
  A(j,8) = matC_loc(j,1)*matC_loc(j,2)*matC_loc(j,3)
END DO
B=TRANSPOSE(A)
DO j=1,nloc
	 A(j,1:8)=A(j,1:8)*funcaoP(j)
END DO
f=MATMUL(B,f1)
G=MATMUL(B,A)

END SUBROUTINE EqN_3D_base_l

!============================================================================================
! 8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!============================================================================================

SUBROUTINE EqN_3D_base_q(matC_loc, no_interp, funcaoP, Uaprox,G ,f , nloc)

IMPLICIT NONE
INTEGER, DIMENSION(:), INTENT(IN)	:: no_interp
INTEGER, INTENT(IN)			:: nloc
REAL(8), DIMENSION(:), INTENT(OUT)	:: f
REAL(8), DIMENSION(:), INTENT(IN)	:: funcaoP,Uaprox
REAL(8), DIMENSION(:,:), INTENT(OUT)	:: G
REAL(8), DIMENSION(:,:), INTENT(IN)	:: matC_loc

REAL(8), DIMENSION(:,:)	:: A(nloc,11), B(11,nloc)
!REAL(8), DIMENSION(:,:)	:: A(nloc,14), B(14,nloc)
REAL(8), DIMENSION(:)	:: f1(nloc)
INTEGER			:: j

DO j=1,nloc
  f1(j)=Uaprox(no_interp(j))*funcaoP(j)
  A(j,1) = 1.d0
  A(j,2) = matC_loc(j,1)
  A(j,3) = matC_loc(j,2)
  A(j,4) = matC_loc(j,3)
  A(j,5) = matC_loc(j,1)*matC_loc(j,2)
  A(j,6) = matC_loc(j,1)*matC_loc(j,3)
  A(j,7) = matC_loc(j,2)*matC_loc(j,3)
  A(j,8) = matC_loc(j,1)*matC_loc(j,2)*matC_loc(j,3)
  A(j,9) = matC_loc(j,1)**2
  A(j,10) = matC_loc(j,2)**2
  A(j,11) = matC_loc(j,3)**2
!  A(j,12) = matC_loc(j,1)**2*matC_loc(j,2)*matC_loc(j,3)
!  A(j,13) = matC_loc(j,2)**2*matC_loc(j,1)*matC_loc(j,3)
!  A(j,14) = matC_loc(j,3)**2*matC_loc(j,2)*matC_loc(j,1)
END DO
B=TRANSPOSE(A)
DO j=1,nloc
	 A(j,1:11)=A(j,1:11)*funcaoP(j)
END DO
f=MATMUL(B,f1)
G=MATMUL(B,A)

END SUBROUTINE EqN_3D_base_q

!============================================================================================
! 8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!============================================================================================

SUBROUTINE ReduzMat( G, n )
!***************************************************************************
!  Triangulariza a matriz global do método de elementos finitos
!  G     - matriz global NxN
!  n     - número de variáveis do sistema
!  Versão 1.0 04/05
!  Autor: Marcos Welby
!***************************************************************************
IMPLICIT NONE
INTEGER(4),PARAMETER::dp = KIND(1.d0)
INTEGER(4),INTENT(in)::n
real(dp),INTENT(inout)::G(n,n)
INTEGER(4)::k,p,m

DO p = 1, n-1
   G(p+1:n,p) = G(p+1:n,p)/G(p,p)
   DO k = p+1, n
      do m = p+1, n
         G(m,k) = G(m,k) - G(m,p)*G(p,k)
      enddo
      !G(p+1:n,k) = G(p+1:n,k) - G(p+1:n,p)*G(p,k)
   ENDDO
ENDDO
RETURN
END SUBROUTINE ReduzMat

!============================================================================================
! 8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!============================================================================================

SUBROUTINE RetroMat( G, f, n )
!***************************************************************************
!  Reduz o vetor fonte e executa a retrosubstituição
!  G     - matriz global NxN
!  f     - vetor fonte de entrada e solução de saída
!  n     - número de variáveis do sistema
!  Versão 1.0 04/05
!  Autor: Marcos Welby
!***************************************************************************
IMPLICIT NONE
INTEGER(4),PARAMETER::dp=KIND(1.d0)
INTEGER(4),INTENT(in)::n
real(dp),INTENT(in)::G(n,n)
real(dp),INTENT(inout)::f(n)
INTEGER(4)::k,p

DO p = 1, n-1
   do k = p+1, n
      f(k) = f(k) - G(k,p)*f(p)
   !f(p+1:n) = f(p+1:n) - G(p+1:n,p)*f(p)
   enddo
ENDDO
!
DO p = n, 1, -1
   f(p) = (f(p) - SUM(G(p,p+1:n)*f(p+1:n)))/G(p,p)
ENDDO
END SUBROUTINE RetroMat

!============================================================================================
! 8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!============================================================================================

SUBROUTINE Ident_ord_nos_segmento33_x(matcoord,vetlinha,vetx,nglob,cont1)

IMPLICIT NONE
INTEGER,INTENT(IN)::nglob,cont1
REAL(8),DIMENSION(:,:),INTENT(IN)::matcoord
REAL(8),DIMENSION(:)::vetx
INTEGER,DIMENSION(:)::vetlinha
INTEGER::i,j,num,nomax,nomin
REAL(8)::valmax,valmin

num=0
DO i=1,nglob
  IF(matcoord(i,4)==33)THEN
    num=num+1
    vetlinha(num)=INT(matcoord(i,1))
    vetx(num)=matcoord(i,2)
  END IF
END DO
num=1
DO i=1,cont1-1
  num=num+1
  valmax=vetx(i)
  nomax=vetlinha(i)
     DO j=num,cont1
       valmin=vetx(j)
       nomin=vetlinha(j)
       IF(valmax < valmin)CYCLE
       vetlinha(i)=nomin
       vetlinha(j)=nomax
       vetx(i)=valmin
       vetx(j)=valmax
       valmax=valmin
       nomax=nomin
     END DO
END DO
END SUBROUTINE Ident_ord_nos_segmento33_x

!============================================================================================
! 8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!============================================================================================

END MODULE deriv_MQMP4
