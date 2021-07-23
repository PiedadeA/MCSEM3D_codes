module modulo_auxiliar

use var_glob
use fem_3D

contains

!%%%===========================================%%%

subroutine cfx(n,x,y,z,f)
implicit none
integer,intent(in)   :: n
real(dpc),intent(in)  :: x,y,z
complex(dpc),intent(out) :: f(n)

f(1) = 0.d0 ; f(2) = 1.d0; f(3) = 0.d0;
f(4) = 0.d0 ; f(5) = y; f(6) = z;
f(7) = 0.d0 ; f(8) = x*2.d0; f(9) = 0.d0;
f(10) = 0.d0;

end subroutine

!%%%===========================================%%%

subroutine cfy(n,x,y,z,f)
implicit none
integer,intent(in)   :: n
real(dpc),intent(in)  :: x,y,z
complex(dpc),intent(out) :: f(n)

f(1) = 0.d0 ; f(2) = 0.d0; f(3) = 1.d0;
f(4) = 0.d0; f(5) = x; f(6) = 0.d0;
f(7) = z; f(8) = 0.d0; f(9) = 2.d0*y;
f(10) = 0.d0;

end subroutine

!%%%===========================================%%%

subroutine cfz(n,x,y,z,f)
implicit none
integer,intent(in)   :: n
real(dpc),intent(in)  :: x,y,z
complex(dpc),intent(out) :: f(n)

f(1) = 0.d0 ; f(2) = 0.d0; f(3) = 0.d0;
f(4) = 1.d0; f(5) = 0.d0; f(6) = x;
f(7) = y; f(8) = 0.d0; f(9) = 0.d0;
f(10) = 2.d0*z;

end subroutine

!%%%===========================================%%%

subroutine cf(n,x,y,z,f)
implicit none
integer,intent(in)   :: n
real(dpc),intent(in)  :: x,y,z
complex(dpc),intent(out) :: f(n)

f(1) = 1.d0 ; f(2) = x; f(3) = y;
f(4) = z; f(5) = x*y; f(6) = x*z;
f(7) = y*z; f(8) = x**2; f(9) = y**2;
f(10) = z**2;

end subroutine

!%%%===========================================%%%

subroutine compose_E(sol, nrcp, xyzr, Ext_in, Eyt_in, Ezt_in, Hx, Hy, Hz)
implicit none
integer,intent(in)       :: nrcp
real(dpc),intent(in)     :: xyzr(nrcp,3)
complex(dpc),intent(in)  :: sol(ngl*nnosg)
complex(dpc),intent(out) :: ext_in(nrcp), eyt_in(nrcp), ezt_in(nrcp)
complex(dpc),intent(out) :: Hx(nrcp), Hy(nrcp), Hz(nrcp)

real(dpc), allocatable   :: xe(:),ye(:),ze(:),dists(:)
integer   :: elems(10000)
real(dpc) :: dr, xm, ym, zm, xr, yr, zr
real(dpc) :: drce
integer   :: i,j,k,l,M,ii
integer,allocatable :: inds(:), indsa(:), indsaa(:)
real(dpc),allocatable       :: cxe(:), cye(:), cze(:)
complex(dpc),allocatable    :: mA(:,:), y(:), yp(:), cW(:,:), mAp(:,:), Ac(:,:), Bc(:)
complex(dpc),allocatable    :: yy(:) , yz(:)

complex(dpc),allocatable    :: aTw(:,:)
integer                     :: Ncs = 10
complex(dpc)                :: f(10), fp(10), fy(10), fz(10)
real(dpc)                   :: cc = 1.d0

dr = rad

Ext_in = 0.d0
Eyt_in = 0.d0
Ezt_in = 0.d0

Hx=Ext_in
Hy=Eyt_in
Hz=Ezt_in

allocate(dists(nnosg))
allocate(xe(nnoele),ye(nnoele),ze(nnoele))

do k = 1, nrcp

 xr = xyzr(k,1)
 yr = xyzr(k,2)
 zr = xyzr(k,3)

 do i = 1,nnosg
  !--
   xm = mat_coord(i,1);
   ym = mat_coord(i,2);
   zm = mat_coord(i,3);
  !--
   drce  = dsqrt( (xr-xm)**2 + (yr-ym)**2 + (zr-zm)**2 )
   dists(i) = drce
  !--
 end do

 M=0 ; elems=0
 do i = 1,nnosg
   if(dists(i)<dr)then
      M = M + 1
      elems(M) = i
   end if
 end do

 allocate(indsaa(M))

 allocate(mA(M,Ncs),y(M),cW(M,M),mAp(M,Ncs),yp(M))
 allocate(Ac(Ncs,Ncs),Bc(Ncs))
 allocate(cxe(M),cye(M),cze(M))
 allocate(yy(M),yz(M))

 do i = 1,M
    indsaa(i) = elems(i)
    cxe(i) = mat_coord(indsaa(i),1);
    cye(i) = mat_coord(indsaa(i),2);
    cze(i) = mat_coord(indsaa(i),3);
 end do
 
 cW = 0.d0

 do i = 1,M

    xm = mat_coord(indsaa(i),1);
    ym = mat_coord(indsaa(i),2);
    zm = mat_coord(indsaa(i),3);

    y(i)  = Sol(4*indsaa(i)-3)
    yy(i) = Sol(4*indsaa(i)-2)
    yz(i) = Sol(4*indsaa(i)-1)
    yp(i) = Sol(4*indsaa(i))    

    call cf(Ncs,xm,ym,zm,f)

    mA(i,:)  = f

    cW(i,i) = dexp( - cc**2 * ( ((xr-xm)/(xr-maxval(cxe)))**2 + ((yr-ym)/(yr-maxval(cye)))**2 &
                    + ((zr-zm)/(zr-maxval(cze)))**2 ) )
    
 end do

 allocate(aTw(Ncs,M))

 do i = 1,Ncs
    do j = 1,M   
       aTw(i,j) = sum( mA(:,i)*cW(:,j) )
    end do
 end do

!----- Ax
 do i = 1,Ncs
    do j = 1,Ncs   
       Ac(i,j) = sum( aTw(i,:)*mA(:,j) )
    end do
    Bc(i) = sum( aTw(i,:)*y(:) )
 end do

 call LU_Gauss_compacto( Ac, Bc, Ncs)
 call cf(Ncs,xr,yr,zr,f)
 Ext_in(k) = sum(f*Bc)

!----- Ay
 do i = 1,Ncs
    do j = 1,Ncs   
       Ac(i,j) = sum( aTw(i,:)*mA(:,j) )
    end do
    Bc(i) = sum( aTw(i,:)*yy(:) )
 end do

 call LU_Gauss_compacto( Ac, Bc, Ncs)
 call cf(Ncs,xr,yr,zr,f)
 Eyt_in(k) = sum(f*Bc)

!----- Az
 do i = 1,Ncs
    do j = 1,Ncs   
       Ac(i,j) = sum( aTw(i,:)*mA(:,j) )
    end do
    Bc(i) = sum( aTw(i,:)*yz(:) )
 end do

 call LU_Gauss_compacto( Ac, Bc, Ncs)
 call cf(Ncs,xr,yr,zr,f)
 Ezt_in(k) = sum(f*Bc) 

!----- Phi

 do i = 1,Ncs
    do j = 1,Ncs   
       Ac(i,j) = sum( aTw(i,:)*mA(:,j) )
    end do
    Bc(i) = sum( aTw(i,:)*yp(:) )
 end do
 call LU_Gauss_compacto( Ac, Bc, Ncs)

!--

 call cfx(Ncs,xr,yr,zr,fp)
 Ext_in(k) = Ext_in(k) + sum(fp(:)*Bc(:)) 
 Ext_in(k) = -zet0 * (Ext_in(k)) 

 call cfy(Ncs,xr,yr,zr,fp)
 Eyt_in(k) = Eyt_in(k) + sum(fp(:)*Bc(:)) 
 Eyt_in(k) = -zet0 * (Eyt_in(k)) 

 call cfz(Ncs,xr,yr,zr,fp)
 Ezt_in(k) = Ezt_in(k) + sum(fp(:)*Bc(:)) 
 Ezt_in(k) = -zet0 * (Ezt_in(k))
!-----

! print*,cdabs(Ext_in(k)),cdabs(Eyt_in(k)),cdabs(Ezt_in(k))

 deallocate(aTw)

 deallocate(yy,yz)
 deallocate(indsaa) 
 deallocate(mA,y,cW,mAp,yp)
 deallocate(Ac,Bc)
 deallocate(cxe,cye,cze)

end do ! Receptores


end subroutine

!%%%===========================================%%%

SUBROUTINE LU_Gauss_compacto( A, b, n)                               
!====================================================================
!     Fatora A em LU e armazena em A. resolver o sistema (U*x) = b e 
!     guardar a solução em b                                         
!     A    - Matriz do sistema          =============================
!     b    - Vetor fonte                - Versão 1.0   01 / 01 / 2012
!     n    - Número de variáveis        - Autor: C.M. Barriga Nunes  
!====================================================================
IMPLICIT NONE
INTEGER,INTENT(in):: n
!REAL(8), DIMENSION(:,:),INTENT(inout):: A(n,n)
!REAL(8), DIMENSION(:),INTENT(inout):: b(n)
COMPLEX(8), DIMENSION(:,:),INTENT(inout):: A(n,n)
COMPLEX(8), DIMENSION(:),INTENT(inout):: b(n)
INTEGER:: i, j, k, cont

A(2:n,1) = A(2:n,1) / A(1,1)
cont=1
DO i=2, n
 cont = cont + 1
   do k=1, i-1
   b(i) = b(i) - b(k) * A(i,k)
   enddo
   DO j=cont, n
     DO k=1, i-1
       A(i,j) = A(i,j) - A(k,j) * A(i,k)
     END DO
     IF( j == n )EXIT
     IF( i == n )EXIT
     DO k=1, i-1
       A(j+1,i) = A(j+1,i) - A(j+1,k) * A(k,i)
     END DO
     A(j+1,i) =A(j+1,i) / A(i,i)
  END DO
END DO
b(n) = b(n) / A(n,n)
DO i=n-1, 1, -1
  b(i) = (b(i) - sum(A(i,i+1:n) * b(i+1:n))) / A(i,i)
END DO
END SUBROUTINE LU_Gauss_compacto

!%%%===========================================%%%

subroutine mesh_qual(comand)
implicit none
character(len=100),intent(out) :: comand

integer :: n
real(dpc),allocatable :: F(:)

open(44,file='./parametros_modelo/frequencias_modelo1d.in',status='old',action='read')

read(44,*)n
allocate(F(n))
read(44,*)F

if(maxval(F)>=0.75d0)then 
   comand="tetgen -Qpq1.17aAemk modelo3d.poly"
else
   comand="tetgen -Qpq1.17aAemk modelo3d.poly"
end if

deallocate(F)
close(44)

end subroutine

!%%%===========================================%%%

subroutine flagsdt
implicit none

integer :: i,j,k,l,c1,c2,c3
real(dpc) :: dx,dy,r

integer :: Ndata

Ndata = nfreqs*sum(NTPP(:)*NRPP(:)) 

allocate(flagF(Ndata))
flagF = 1;

!&&&&&&&&&&&&&&&&&&&&
! c2 = 0 ;
!do l = 1, nfreqs
!	 c1 = 0; 
!	do i = 1,num_perf
!	  do j = 1,NTPP(i)
!	    c1 = c1 + 1;
!	      do k = 1,NRPP(i);
!	 
!		  dx = COOTXP(c1,1)-COORCP((i-1)*NRPP(i) + k,1) ; 
!		  dy = COOTXP(c1,2)-COORCP((i-1)*NRPP(i) + k,2) ;
!
!                  r = dsqrt( dx**2 + dy**2 )

!		  if( ( r ).ge.minper .and. &
!	              ( r ).le.maxper )then
!		    c2 = c2 + 1
!                    flagF(c2) = 1
!                  else
!		    c2 = c2 + 1                   
!                    flagF(c2) = 0
!		  end if

!	      end do
!	  end do
!	end do
!end do
!&&&&&&&&&&&&&&&&&&&&

end subroutine

!%%%===========================================%%%

subroutine arquivos_plot
implicit none

integer :: i,j,k,l,c1,c2,c3
real(dpc) :: dx,dy,r

 open(44,file='dados1.dat',status='replace',action='write')
 open(45,file='dados2.dat',status='replace',action='write')

 write(44,*)nfreqs 
 write(44,*)num_perf 
 do i = 1,num_perf
   write(44,*)NTPP(i) 
 end do 

!&&&&&&&&&&&&&&&&&&&&
 c3 = 0;
do l = 1, nfreqs
	 c1 = 0; 
	do i = 1,num_perf
	  do j = 1,NTPP(i)
	    c1 = c1 + 1; c2 = 0 ;
	      do k = 1,NRPP(i);
	 
		    c3 = c3 + 1;
		    if( flagF(c3).eq.1 ) c2 = c2 + 1
                    !c2=c2+1

	      end do
	            write(45,*) c2 
	  end do
	end do
end do
!&&&&&&&&&&&&&&&&&&&&

 close(44)
 close(45)

end subroutine

!%%%===========================================%%%

subroutine Estima_num_dados( num_fields, NDPF ) 
implicit none
integer,intent(out) :: num_fields,NDPF(:)

integer :: i,j,k,l,c1,c2,c3,c4
real(dpc) :: dx,dy,r

 num_fields = 0

!&&&&&&&&&&&&&&&&&&&&
 c3 = 0; c2 = 0;
do l = 1, nfreqs
	 c1 = 0; c4 = 0;
	do i = 1,num_perf
	  do j = 1,NTPP(i)
	    c1 = c1 + 1; 
	      do k = 1,NRPP(i);

		    c3 = c3 + 1;
		    if( flagF(c3).eq.1 )then
			 c2 = c2 + 1 ; c4 = c4 + 1
		    end if
	 
	      end do 
	  end do
	end do
	NDPF(l) = c4;
end do
!&&&&&&&&&&&&&&&&&&&&

num_fields = c2 ! Depois de uma recontagem aqui e atribuido o valor correto do numero de campos.

end subroutine

!%%%===========================================%%%

subroutine Numrec_portransm2( lp, lp_frq , trsm , xtr , ytr , nrc, ind_rec, ind_rec_proxma , ind_rec_proxme )
implicit none

integer,intent(in)  :: lp , lp_frq , trsm
real(dpc),intent(in) :: xtr , ytr

integer,intent(out) :: nrc
integer,allocatable,intent(out) :: ind_rec(:), ind_rec_proxma(:), ind_rec_proxme(:)

integer,allocatable :: ind_reca(:), ind_rec_proxmaa(:), ind_rec_proxmea(:)

integer :: i,j,k,l,c1,c2,c3,c4
real(dpc) :: dx,dy,r
integer :: lp_perf,aux

allocate(ind_reca(NRPP(lp)),ind_rec_proxmaa(NRPP(lp)),ind_rec_proxmea(NRPP(lp)))

allocate(NDPF(nfreqs)); call Estima_num_dados( aux, NDPF ) 

if(lp_frq.gt.1)then
 c1 = (lp_frq-1)*sum(NTPP(:)*NRPP(:)); 
else
 c1 = 0; 
end if

!print*,'c1=',c1

do i = 1,lp

     if(i.lt.lp)then

        do j = 1,NTPP(i)
              do k = 1,NRPP(i);
	           c1 = c1 + 1; 
              end do 
        end do 

     else 

        do j = 1,NTPP(i)

             if(j.lt.trsm)then
                do k = 1,NRPP(i);
	             c1 = c1 + 1; 
                end do 
             elseif(j==trsm)then
                c2 = 0
                do k = 1,NRPP(i);	 
	           c1 = c1 + 1; 
	           if( flagF(c1).eq.1 )then
		       c2 = c2 + 1 ;
	               ind_reca(c2) = INDRCP( (i-1)*NRPP(i) + k )
		       ind_rec_proxmaa(c2) = INDRCPX( (i-1)*NRPP(i) + k,1 )
		       ind_rec_proxmea(c2) = INDRCPX( (i-1)*NRPP(i) + k,2 )
	           end if
                 end do
                 !print*,'c2=',c2 
             end if

        end do


     end if

end do

nrc = c2;

allocate(ind_rec(nrc),ind_rec_proxma(nrc),ind_rec_proxme(nrc))

do i = 1,nrc
   ind_rec(i) = ind_reca(i)
   ind_rec_proxma(i) = ind_rec_proxmaa(i)
   ind_rec_proxme(i) = ind_rec_proxmea(i)
end do

!print*,'nrc',nrc

deallocate(ind_reca,ind_rec_proxmaa,ind_rec_proxmea)

deallocate(NDPF)

end subroutine

!%%%===========================================%%%

subroutine deriv_fbsEH( sol, n , coord , Ex, Ey, Ez, Hx, Hy, Hz )
implicit none
integer,intent(in)   :: n
real(dpc),intent(in) :: coord(n,3)
complex(dpc),intent(in)  :: sol(:) 
complex(dpc),intent(out) :: Ex(n), Ey(n), Ez(n), Hx(n), Hy(n), Hz(n)

complex(dpc),allocatable :: Ax(:),Ay(:),Az(:),phi(:)
complex(dpc),allocatable :: lAx(:),lAy(:),lAz(:),lphi(:)
real(dpc)                :: delx,dely,modu
integer :: i,j,ie,k,outl
real(dpc),allocatable :: xe(:),ye(:),ze(:)
real(dpc),allocatable :: xa(:),ya(:),za(:),b(:),c(:),a(:),d(:)
complex(dpc) :: delphix,delphiy
real(dpc) :: vol
complex(dpc),allocatable :: phis(:)

real(dpc) :: dr, xm, ym, zm, xr, yr, zr, drce
integer   :: M
real(dpc),allocatable :: dists(:), Ni(:)

allocate(dists(nelemg))

allocate(Ax(nnosg))
allocate(Ay(nnosg))
allocate(Az(nnosg))
allocate(phi(nnosg))


allocate(lAx(4))
allocate(lAy(4))
allocate(lAz(4))
allocate(lphi(4))
allocate(Ni(4))

allocate(xa(3),ya(3),za(3))
allocate(xe(4),ye(4),ze(4),b(4),c(4),a(4),d(4))

do i = 1,nnosg
	Ax(i)  = Sol( 4*i-3 + (irs-1)*ngb )
	Ay(i)  = Sol( 4*i-2 + (irs-1)*ngb )
	Az(i)  = Sol( 4*i-1 + (irs-1)*ngb )
	phi(i) = Sol( 4*i + (irs-1)*ngb   )
end do

do i = 1,n

        xr = coord(i,1)
        yr = coord(i,2)
        zr = coord(i,3)

        do j = 1,nelemg
            xe = mat_coord(mat_elem(j,:),1);
            ye = mat_coord(mat_elem(j,:),2);
            ze = mat_coord(mat_elem(j,:),3);
            xm = sum(xe)/nnoele
            ym = sum(ye)/nnoele
            zm = sum(ze)/nnoele

            drce  = dsqrt( (xr-xm)**2 + (yr-ym)**2 + (zr-zm)**2 )
            dists(j) = drce
        end do

        ie = minloc(dists,dim=1)

        xe = mat_coord(mat_elem(ie,:),1);
        ye = mat_coord(mat_elem(ie,:),2);
        ze = mat_coord(mat_elem(ie,:),3);

       !#--

        xa(1) = xe(2) ; xa(2) = xe(3) ; xa(3) = xe(4) ;
        ya(1) = ye(2) ; ya(2) = ye(3) ; ya(3) = ye(4) ; 
        za(1) = ze(2) ; za(2) = ze(3) ; za(3) = ze(4) ;  a(1) = ((-1.d0)**(2+3+4+1)) * det(xa,ya,za);

        xa(1) = xe(3) ; xa(2) = xe(4) ; xa(3) = xe(1) ;
        ya(1) = ye(3) ; ya(2) = ye(4) ; ya(3) = ye(1) ; 
        za(1) = ze(3) ; za(2) = ze(4) ; za(3) = ze(1) ;  a(2) = ((-1.d0)**(3+4+1+1)) * det(xa,ya,za);

        xa(1) = xe(4) ; xa(2) = xe(1) ; xa(3) = xe(2) ;
        ya(1) = ye(4) ; ya(2) = ye(1) ; ya(3) = ye(2) ; 
        za(1) = ze(4) ; za(2) = ze(1) ; za(3) = ze(2) ;  a(3) = ((-1.d0)**(4+1+2+1)) * det(xa,ya,za);

        xa(1) = xe(1) ; xa(2) = xe(2) ; xa(3) = xe(3) ;
        ya(1) = ye(1) ; ya(2) = ye(2) ; ya(3) = ye(3) ; 
        za(1) = ze(1) ; za(2) = ze(2) ; za(3) = ze(3) ;  a(4) = ((-1.d0)**(1+2+3+1)) * det(xa,ya,za);

        !#---

	xa(1) = 1.d0  ; xa(2) = 1.d0  ; xa(3) = 1.d0  ;
	ya(1) = ye(2) ; ya(2) = ye(3) ; ya(3) = ye(4) ; 
	za(1) = ze(2) ; za(2) = ze(3) ; za(3) = ze(4) ;  b(1) = ((-1.d0)**(2+3+4)) * det(xa,ya,za);

	xa(1) = 1.d0  ; xa(2) = 1.d0   ; xa(3) = 1.d0 ;
	ya(1) = ye(3) ; ya(2) = ye(4) ; ya(3) = ye(1) ; 
	za(1) = ze(3) ; za(2) = ze(4) ; za(3) = ze(1) ;  b(2) = ((-1.d0)**(3+4+1)) * det(xa,ya,za);

	xa(1) = 1.d0  ; xa(2) = 1.d0  ; xa(3) = 1.d0  ;
	ya(1) = ye(4) ; ya(2) = ye(1) ; ya(3) = ye(2) ; 
	za(1) = ze(4) ; za(2) = ze(1) ; za(3) = ze(2) ;  b(3) = ((-1.d0)**(4+1+2)) * det(xa,ya,za);

	xa(1) = 1.d0  ; xa(2) = 1.d0  ; xa(3) = 1.d0  ;
	ya(1) = ye(1) ; ya(2) = ye(2) ; ya(3) = ye(3) ; 
	za(1) = ze(1) ; za(2) = ze(2) ; za(3) = ze(3) ;  b(4) = ((-1.d0)**(1+2+3)) * det(xa,ya,za);

        !#---

	xa(1) = xe(2)  ; xa(2) = xe(3)  ; xa(3) = xe(4);
	ya(1) = 1.d0 ; ya(2) = 1.d0 ; ya(3) = 1.d0 ; 
	za(1) = ze(2) ; za(2) = ze(3) ; za(3) = ze(4) ;  c(1) = ((-1.d0)**(2+3+4)) * det(xa,ya,za);

	xa(1) = xe(3)  ; xa(2) = xe(4)  ; xa(3) = xe(1);
	ya(1) = 1.d0 ; ya(2) = 1.d0 ; ya(3) = 1.d0 ; 
	za(1) = ze(3) ; za(2) = ze(4) ; za(3) = ze(1) ;  c(2) = ((-1.d0)**(3+4+1)) * det(xa,ya,za);

	xa(1) = xe(4)  ; xa(2) = xe(1)  ; xa(3) = xe(2);
	ya(1) = 1.d0 ; ya(2) = 1.d0 ; ya(3) = 1.d0 ; 
	za(1) = ze(4) ; za(2) = ze(1) ; za(3) = ze(2) ;  c(3) = ((-1.d0)**(4+1+2)) * det(xa,ya,za);

	xa(1) = xe(1)  ; xa(2) = xe(2)  ; xa(3) = xe(3);
	ya(1) = 1.d0 ; ya(2) = 1.d0 ; ya(3) = 1.d0 ; 
	za(1) = ze(1) ; za(2) = ze(2) ; za(3) = ze(3) ;  c(4) = ((-1.d0)**(1+2+3)) * det(xa,ya,za);

        !#---

        xa(1) = xe(2)  ; xa(2) = xe(3)  ; xa(3) = xe(4);
        ya(1) = ye(2) ; ya(2) = ye(3) ; ya(3) = ye(4); 
        za(1) = 1.d0 ; za(2) = 1.d0 ; za(3) = 1.d0;      d(1) = ((-1.d0)**(2+3+4)) * det(xa,ya,za);

        xa(1) = xe(3)  ; xa(2) = xe(4)  ; xa(3) = xe(1);
        ya(1) = ye(3) ; ya(2) = ye(4) ; ya(3) = ye(1); 
        za(1) = 1.d0 ; za(2) = 1.d0 ; za(3) = 1.d0;      d(2) = ((-1.d0)**(3+4+1)) * det(xa,ya,za);

        xa(1) = xe(4)  ; xa(2) = xe(1)  ; xa(3) = xe(2);
        ya(1) = ye(4) ; ya(2) = ye(1) ; ya(3) = ye(2); 
        za(1) = 1.d0 ; za(2) = 1.d0 ; za(3) = 1.d0;      d(3) = ((-1.d0)**(4+1+2)) * det(xa,ya,za);

        xa(1) = xe(1)  ; xa(2) = xe(2)  ; xa(3) = xe(3);
        ya(1) = ye(1) ; ya(2) = ye(2) ; ya(3) = ye(3); 
        za(1) = 1.d0 ; za(2) = 1.d0 ; za(3) = 1.d0;      d(4) = ((-1.d0)**(1+2+3)) * det(xa,ya,za);

        !#--

        vol = vol_tetraedro(xe,ye,ze) 

        lphi = phi(mat_elem(ie,:))
        lAx  = Ax(mat_elem(ie,:))
        lAy  = Ay(mat_elem(ie,:))
        lAz  = Az(mat_elem(ie,:))

        Ni = (1.d0/(6.d0*vol))*(a + b*xr + c*yr + d*zr)        

        Ex(i) = -zet0 * ( sum( (lphi(:)*b(:)/(6.d0*vol)) + lAx(:)*Ni(:) ) )
        Ey(i) = -zet0 * ( sum( (lphi(:)*c(:)/(6.d0*vol)) + lAy(:)*Ni(:) ) )
        Ez(i) = -zet0 * ( sum( (lphi(:)*d(:)/(6.d0*vol)) + lAz(:)*Ni(:) ) )

        Hx(i) = 1.d0/(6.d0*vol) * ( sum( lAz(:)*c(:) - lAy(:)*d(:) ) )
        Hy(i) = 1.d0/(6.d0*vol) * ( sum( lAx(:)*d(:) - lAz(:)*b(:) ) )
        Hz(i) = 1.d0/(6.d0*vol) * ( sum( lAy(:)*b(:) - lAx(:)*c(:) ) )

end do


end subroutine

!%%%===========================================%%%

subroutine deriv_fbs( sol, n , ind , indpma , E )
implicit none
integer,intent(in) :: n
integer,intent(in) :: ind(n), indpma(n)
complex(dpc),intent(in)  :: sol(:) 
complex(dpc),intent(out) :: E(n)

complex(dpc),allocatable :: Ax(:),Ay(:),Az(:),phi(:)
real(dpc)                :: delx,dely,modu
integer :: i,j,ie,k,outl
real(dpc),allocatable :: xe(:),ye(:),ze(:)
real(dpc),allocatable :: xa(:),ya(:),za(:),b(:),c(:)
complex(dpc) :: delphix,delphiy
real(dpc) :: vol
complex(dpc),allocatable :: phis(:)

allocate(Ax(nnosg))
allocate(Ay(nnosg))
allocate(Az(nnosg))
allocate(phi(nnosg))

allocate(xa(3),ya(3),za(3))
allocate(xe(4),ye(4),ze(4),b(4),phis(4),c(4))

do i = 1,nnosg
	Ax(i)  = Sol( 4*i-3 + (irs-1)*ngb )
	Ay(i)  = Sol( 4*i-2 + (irs-1)*ngb )
	Az(i)  = Sol( 4*i-1 + (irs-1)*ngb )
	phi(i) = Sol( 4*i + (irs-1)*ngb   )
end do

do i = 1,n

	delx= mat_coord(indpma(i),1) - mat_coord(ind(i),1);
	dely= mat_coord(indpma(i),2) - mat_coord(ind(i),2);
	modu= dsqrt( delx**2 + dely**2 );
 
        outl=0
        do j = 1,nelemg
             do k = 1,nnoele
                if(ind(i)==mat_elem(j,k))then
                  ie = j ; outl=1 ;exit
                end if
             end do
             if(outl==1)exit
        end do

        xe = mat_coord(mat_elem(ie,:),1);
        ye = mat_coord(mat_elem(ie,:),2);
        ze = mat_coord(mat_elem(ie,:),3);

        !#---

	xa(1) = 1.d0  ; xa(2) = 1.d0  ; xa(3) = 1.d0  ;
	ya(1) = ye(2) ; ya(2) = ye(3) ; ya(3) = ye(4) ; 
	za(1) = ze(2) ; za(2) = ze(3) ; za(3) = ze(4) ;  b(1) = ((-1.d0)**(2+3+4)) * det(xa,ya,za);

	xa(1) = 1.d0  ; xa(2) = 1.d0   ; xa(3) = 1.d0 ;
	ya(1) = ye(3) ; ya(2) = ye(4) ; ya(3) = ye(1) ; 
	za(1) = ze(3) ; za(2) = ze(4) ; za(3) = ze(1) ;  b(2) = ((-1.d0)**(3+4+1)) * det(xa,ya,za);

	xa(1) = 1.d0  ; xa(2) = 1.d0  ; xa(3) = 1.d0  ;
	ya(1) = ye(4) ; ya(2) = ye(1) ; ya(3) = ye(2) ; 
	za(1) = ze(4) ; za(2) = ze(1) ; za(3) = ze(2) ;  b(3) = ((-1.d0)**(4+1+2)) * det(xa,ya,za);

	xa(1) = 1.d0  ; xa(2) = 1.d0  ; xa(3) = 1.d0  ;
	ya(1) = ye(1) ; ya(2) = ye(2) ; ya(3) = ye(3) ; 
	za(1) = ze(1) ; za(2) = ze(2) ; za(3) = ze(3) ;  b(4) = ((-1.d0)**(1+2+3)) * det(xa,ya,za);

        !#---

	xa(1) = xe(2)  ; xa(2) = xe(3)  ; xa(3) = xe(4);
	ya(1) = 1.d0 ; ya(2) = 1.d0 ; ya(3) = 1.d0 ; 
	za(1) = ze(2) ; za(2) = ze(3) ; za(3) = ze(4) ;  c(1) = ((-1.d0)**(2+3+4)) * det(xa,ya,za);

	xa(1) = xe(3)  ; xa(2) = xe(4)  ; xa(3) = xe(1);
	ya(1) = 1.d0 ; ya(2) = 1.d0 ; ya(3) = 1.d0 ; 
	za(1) = ze(3) ; za(2) = ze(4) ; za(3) = ze(1) ;  c(2) = ((-1.d0)**(3+4+1)) * det(xa,ya,za);

	xa(1) = xe(4)  ; xa(2) = xe(1)  ; xa(3) = xe(2);
	ya(1) = 1.d0 ; ya(2) = 1.d0 ; ya(3) = 1.d0 ; 
	za(1) = ze(4) ; za(2) = ze(1) ; za(3) = ze(2) ;  c(3) = ((-1.d0)**(4+1+2)) * det(xa,ya,za);

	xa(1) = xe(1)  ; xa(2) = xe(2)  ; xa(3) = xe(3);
	ya(1) = 1.d0 ; ya(2) = 1.d0 ; ya(3) = 1.d0 ; 
	za(1) = ze(1) ; za(2) = ze(2) ; za(3) = ze(3) ;  c(4) = ((-1.d0)**(1+2+3)) * det(xa,ya,za);

        !#---

        phis = phi(mat_elem(ie,:))
        vol = vol_tetraedro(xe,ye,ze) 

        delphix = 0.d0
        delphiy = 0.d0
        do j = 1,nnoele
             delphix=delphix+(phis(j)*b(j)/(6.d0*vol))
             delphiy=delphiy+(phis(j)*c(j)/(6.d0*vol))
        end do 

        E(i) = -zet0 * ( Ax(ind(i))*(delx/modu) + Ay(ind(i))*(dely/modu) + delphix*(delx/modu) + delphiy*(dely/modu)  )

end do

deallocate(xa,ya,za)
deallocate(xe,ye,ze,b,phis,c)

deallocate(Ax)
deallocate(Ay)
deallocate(Az)
deallocate(phi)

end subroutine

!%%%===========================================%%%

subroutine deriv_defin( sol, n , indpma , indpme , ind , E )
implicit none
integer,intent(in) :: n
integer,intent(in) :: indpma(n), indpme(n), ind(n)
complex(dpc),intent(in)  :: sol(:) 
complex(dpc),intent(out) :: E(n)

complex(dpc),allocatable :: Ax(:),Ay(:),Az(:),phi(:)
real(dpc)                :: delx,dely,modu
integer :: i

allocate(Ax(nnosg))
allocate(Ay(nnosg))
allocate(Az(nnosg))
allocate(phi(nnosg))

do i = 1,nnosg
	Ax(i)  = Sol( 4*i-3 + (irs-1)*ngb )
	Ay(i)  = Sol( 4*i-2 + (irs-1)*ngb )
	Az(i)  = Sol( 4*i-1 + (irs-1)*ngb )
	phi(i) = Sol( 4*i + (irs-1)*ngb   )
end do

do i = 1,n

	delx= mat_coord(indpma(i),1) - ( mat_coord(indpma(i),1) - 3.d0 ) ;
	dely= mat_coord(indpma(i),2) - mat_coord(ind(i),2);
	modu= dsqrt( delx**2 + dely**2 );

	!E(i) = -zet0 * ( ( phi(indpma(i)) - phi(indpme(i)) )/(2.d0*modu) + Ax(ind(i))*(delx/modu) + Ay(ind(i))*(dely/modu) )

	E(i) = -zet0 * ( ( phi(indpma(i)) - phi(indpme(i)) )/(2.d0*modu) + &
                       ( (Ax(indpma(i)) + Ax(indpme(i)))/2.d0 )*(delx/modu) + &
                       ( (Ay(indpma(i)) + Ay(indpme(i)))/2.d0 )*(dely/modu) )

        !write(1,*) dreal( phi(ind(i)) ), dreal( Ax(ind(i)) ), dreal( Ay(ind(i)) ), dimag( phi(ind(i)) ), &
	!	   dimag( Ax(ind(i)) ), dimag( Ay(ind(i)) );

end do

deallocate(Ax)
deallocate(Ay)
deallocate(Az)
deallocate(phi)

end subroutine

!%%%===========================================%%%

subroutine inf_memory (my_rank)  

! use ifport !if on intel compiler

implicit none

integer,intent(in) :: my_rank

character(len=30) :: pid_char
integer :: valueRSS , pid
real(8) :: rmem_MB , rmem_GB

call system_mem_usage(valueRSS)

!transgormando para GB ###

rmem_MB = (1.d0*valueRSS)/1024.d0

rmem_GB  = (1.d0*valueRSS)/(1024.d0*1024.d0)

!#########################

write(*,*)' '
write(*,'(A33,I5,A2)')'Memoria utilizada no processo',my_rank,':'
write(*,'(F10.3,A4)')rmem_GB,'GB'
write(*,*)' '
 
end subroutine

!%%%===========================================%%%

subroutine system_mem_usage(valueRSS)

 use ifport !if on intel compiler

implicit none

 integer, intent(out) :: valueRSS

 character(len=200):: filename=' '
 character(len=80) :: line
 character(len=8)  :: pid_char=' '
 integer :: pid
 logical :: ifxst

valueRSS=-1    ! return negative number if not found

!--- get process ID

pid=getpid()

write(pid_char,'(I8)') pid
filename='/proc/'//trim(adjustl(pid_char))//'/status'

!--- read system file

inquire (file=filename,exist=ifxst)
if (.not.ifxst) then
  write (*,*) 'system file does not exist'
  return
endif

open(unit=100, file=filename, action='read')
do
  read (100,'(a)',end=120) line
  if (line(1:6).eq.'VmRSS:') then
     read (line(7:),*) valueRSS
     exit
  endif
enddo

120 continue
close(100)

return
end subroutine

!%%%===========================================%%%

subroutine open_arqprim (  )
!subroutine open_arqprim ( idp )
implicit none

! character(len=100),intent(in) :: idp

 character(len=100) :: pathprim
 integer :: lenc

!pathprim = idp
lenc = len_trim(pathprim)
  
open(1111,file='./primarios/reEx_s.bin',status='old',action='read', access='direct', form='unformatted', recl=8  )
open(1112,file='./primarios/imEx_s.bin',status='old',action='read', access='direct', form='unformatted', recl=8  )
open(1113,file='./primarios/reEy_s.bin',status='old',action='read', access='direct', form='unformatted', recl=8  )
open(1114,file='./primarios/imEy_s.bin',status='old',action='read', access='direct', form='unformatted', recl=8  )
open(1115,file='./primarios/reEz_s.bin',status='old',action='read', access='direct', form='unformatted', recl=8  )
open(1116,file='./primarios/imEz_s.bin',status='old',action='read', access='direct', form='unformatted', recl=8  )

end subroutine

!%%%===========================================%%%

subroutine close_arqprim
implicit none

 close(1111);close(1112);close(1113)
 close(1114);close(1115);close(1116)

end subroutine

!%%%===========================================%%%

subroutine desaloc_vglob
implicit none

deallocate( ihet1,Inh )

end subroutine

!%%%===========================================%%%

subroutine primarios_modelo2(activ,np,nem1,ihet1,Inh,NInh,n_p,my_rank)
implicit none
include 'mpif.h'

! character(len=100),intent(in) :: idp
 integer,intent(in) :: n_p, my_rank, activ, np
 integer,intent(out)  :: nem1 , NInh
 integer,allocatable,intent(out) :: ihet1(:),Inh(:)

 character(len=100) :: pathprim
 integer :: lenc
 integer :: no1,fri,frf,err,aux
 integer,allocatable :: nd(:)
 real(dpc),allocatable :: cnd(:,:)
 integer :: i,j,k,l,m,ct,noi
 real(dpc) :: x,y,z 
 complex(dpc) :: Ex,Ey,Ez,Haux
 real(dpc) :: t1,t2
 integer :: ne1, nel, e1, lb1, ki
 integer,allocatable :: vlb(:)
 integer :: cntt

 real(dpc) :: delx, dely

!pathprim=idp
!lenc = len_trim(pathprim)

if( n_p.gt.nfreqs )then
   if(my_rank.eq.0) print*,'  '
   if(my_rank.eq.0) print*,'============== ATENCAO ! ========================='
   if(my_rank.eq.0) print*,' Programa interrompido, por favor escolha a quantidade de processos '
   if(my_rank.eq.0) print*,' menor ou igual o numero de frequencias ! '
   if(my_rank.eq.0) print*,'=================================================='
   if(my_rank.eq.0) print*,'  '
   if(n_p.gt.1) call MPI_BARRIER(MPI_COMM_WORLD,err)
   call MPI_ABORT(MPI_COMM_WORLD,err)
end if

open(10,file='./malha/'//trim(mname)//'.1.ele',status='old',action='read' )
!open(10,file=pathprim(1:lenc)//'malha/modelo.1.ele',status='old',action='read')

read(10,*)ne1
allocate(vlb(ne1))
nel = 0
do i = 1,ne1
           read(10,*)e1,aux,aux,aux,aux,lb1
           do j = 1,np                
                 if( lb1.eq.( 10*(ncm + 1) + j*10) )then
                   nel = nel + 1
                   vlb(nel) = e1
                   exit
                 end if                   
           end do            
end do 
nem1=nel ;
 close(10)
allocate( ihet1(nel) )
do i = 1, nel
     ihet1(i) = vlb(i) 
end do
deallocate(vlb)

open(30,file='./malha/'//trim(mname)//'.1.node',status='old',action='read' )
!open(30,file=pathprim(1:lenc)//'malha/modelo.1.node',status='old',action='read' )
read(30,*)no1

allocate(nd(no1),cnd(no1,3))

if(my_rank.eq.0) print*,'  '
if(my_rank.eq.0) print*,'Numero de nos da malha', no1

 ct = 0
do i = 1,no1
        read(30,*)noi,x,y,z
	if( ( (x.ge.xih).and.(x.le.xfh) ).and.( (y.ge.yih).and.(y.le.yfh) ).and.( (z.ge.zih).and.(z.le.zfh) ) )then
		ct = ct + 1
		nd(ct)=noi
		cnd(ct,1) =x ; cnd(ct,2) =y ; cnd(ct,3) =z ;
	end if
end do

if(my_rank.eq.0) print*,'  '
if(my_rank.eq.0) print*,'Numero de nos nas heterogeneidades', ct
if(my_rank.eq.0) print*,'  '

NInh = ct;

 close(30)

allocate(Inh(no1))

Inh = 0;
do i = 1,ct
	Inh(nd(i))=i
end do

if(n_p.gt.1)call MPI_BARRIER(MPI_COMM_WORLD,err)

if(activ.eq.1)then

   if(my_rank.eq.0)then
     t1 = MPI_WTIME()
     print*,'     '
     print*,'Determinando os campos primarios' 
     print*,'     '
   end if

   open(350,file='./primarios/reEx_s.bin',status='replace',action='write', access='direct', &
 form='unformatted', recl=8   )
   open(360,file='./primarios/imEx_s.bin',status='replace',action='write', access='direct', &
 form='unformatted', recl=8   )
   open(370,file='./primarios/reEy_s.bin',status='replace',action='write', access='direct', &
 form='unformatted', recl=8   )
   open(380,file='./primarios/imEy_s.bin',status='replace',action='write', access='direct', &
 form='unformatted', recl=8   )
   open(390,file='./primarios/reEz_s.bin',status='replace',action='write', access='direct', &
 form='unformatted', recl=8   )
   open(3100,file='./primarios/imEz_s.bin',status='replace',action='write', access='direct', &
 form='unformatted', recl=8  )

  !=-=-=-=
  call para_range( 1,nfreqs,n_p,my_rank,fri,frf )
  !=-=-=-=

  do i = fri,frf
	cntt = 0;
	do m = 1,num_perf

	       !=========
		!delx = COORCP( sum(NRPP(1:m-1))+1 ,1) - COORCP( sum(NRPP(1:m-1))+2 ,1)
		!dely = COORCP( sum(NRPP(1:m-1))+1 ,2) - COORCP( sum(NRPP(1:m-1))+2 ,2)
		ttap = 0.d0!datan(dely/delx); !print*,'tta:',ttap
	       !=========

		  do j = 1,NTPP(m)
		           cntt = cntt + 1;
        		   print'(a,1x,I5,1x,a,1x,I5,1x,a,1x,I5)','frequencia',i,'perfil',m,'transmissor',j
				do l = 1,ct

				        tr_x = COOTXP(cntt,1) ; tr_y = COOTXP(cntt,2) ; tr_z = COOTXP(cntt,3)

					call prim_MCSEM( ttap, Freqs(i) , tr_x , tr_y , tr_z , & 
                                  			 cnd( l , 1 ),cnd( l , 2 ),cnd( l , 3 ), & 
                                  			 Ex, Ey, Ez )

					if(m.eq.1)then
						k = (i-1)*sum( NTPP(:) )*ct + (j-1)*ct + l;
					else
						k = (i-1)*sum( NTPP(:) )*ct + sum( NTPP(1:m-1) )*ct + (j-1)*ct + l;
					end if

			                write(350,rec=k) dreal(Ex) 
			                write(370,rec=k) dreal(Ey) 
			                write(390,rec=k) dreal(Ez) 
	
			                write(360,rec=k) dimag(Ex) 
			                write(380,rec=k) dimag(Ey) 
			                write(3100,rec=k) dimag(Ez) 

				end do
		  end do ! Trnasmissor


	end do	! Perfil

  end do ! Frequencia

  close(350) ; close(360) 
  close(370) ; close(380) ; close(390) 
  close(3100)

  if(n_p.gt.1)call MPI_BARRIER(MPI_COMM_WORLD,err)

  if(my_rank.eq.0)then
     t2 = MPI_WTIME()
     !call cpu_time(t2) 
     print*,'     '
     print*,'Tempo para calcular os campos primarios',(t2-t1)/(60.d0),'minutos'
     print*,'     '
  end if

!stop

else

  if(my_rank.eq.0)then
    print*,'     '
    print*,'Campos primarios obtidos.' 
    print*,'     '
  end if

end if

deallocate(nd,cnd)

end subroutine

!%%%===========================================%%%

subroutine para_range(n1, n2, nprocs, irank, ista, iend)
implicit none

integer,intent(in)  :: n1,n2,nprocs,irank
integer,intent(out) :: ista,iend
integer :: iwork1,iwork2

iwork1 = (n2 - n1 + 1) / nprocs
iwork2 = mod(n2 - n1 + 1, nprocs)
ista = irank * iwork1 + n1 + min(irank, iwork2)
iend = ista + iwork1 - 1
if (iwork2 > irank) iend = iend + 1

end subroutine para_range

!%%%===========================================%%%


subroutine primarios_modelo(activ,np,nem1,ihet1,Inh,NInh,n_p,my_rank)
implicit none
include 'mpif.h'

! character(len=100),intent(in) :: idp
 integer,intent(in) :: n_p, my_rank, activ, np
 integer,intent(out)  :: nem1 , NInh
 integer,allocatable,intent(out) :: ihet1(:),Inh(:)

 character(len=100) :: pathprim
 integer :: lenc
 integer :: no1,fri,frf,err,aux
 integer,allocatable :: nd(:)
 real(dpc),allocatable :: cnd(:,:)
 integer :: i,j,k,l,m,ct,noi
 real(dpc) :: x,y,z 
 complex(dpc) :: Ex,Ey,Ez,Haux
 real(dpc) :: t1,t2
 integer :: ne1, nel, e1, lb1, ki
 integer,allocatable :: vlb(:),inde(:,:),nd2(:)
 integer :: cntt

 real(dpc) :: delx, dely

!pathprim=idp
!lenc = len_trim(pathprim)

if( n_p.gt.nfreqs )then
   if(my_rank.eq.0) print*,'  '
   if(my_rank.eq.0) print*,'============== ATENCAO ! ========================='
   if(my_rank.eq.0) print*,' Programa interrompido, por favor escolha a quantidade de processos '
   if(my_rank.eq.0) print*,' menor ou igual o numero de frequencias ! '
   if(my_rank.eq.0) print*,'=================================================='
   if(my_rank.eq.0) print*,'  '
   if(n_p.gt.1) call MPI_BARRIER(MPI_COMM_WORLD,err)
   call MPI_ABORT(MPI_COMM_WORLD,err)
end if

open(10,file='./malha/'//trim(mname)//'.1.ele',status='old',action='read' )
!open(10,file=pathprim(1:lenc)//'malha/modelo.1.ele',status='old',action='read')


!-- Aqui e verificado quais elementos estao na regiao de heterogeneidades --
read(10,*)ne1
allocate(vlb(ne1))
nel = 0
do i = 1,ne1
           read(10,*)e1,aux,aux,aux,aux,lb1
           do j = 1,np                
                 if( lb1.eq.( 10*(ncm + 1) + j*10) )then
                   nel = nel + 1
                   vlb(nel) = e1
                   exit
                 end if                   
           end do            
end do 
nem1=nel ;
 close(10)
allocate( ihet1(nel) ) !-- vetor que guarda os elementos que estao nas heterogeneidades --
do i = 1, nel
     ihet1(i) = vlb(i) 
end do
deallocate(vlb)

open(30,file='./malha/'//trim(mname)//'.1.node',status='old',action='read' )

read(30,*)no1

allocate(nd(no1),cnd(no1,3))
allocate(nd2(no1))
!--
!if(my_rank.eq.0) print*,'  '
!if(my_rank.eq.0) print*,'Numero de nos da malha', no1

!--
nd = 0
do i = 1,nem1
   e1 = ihet1(i)
   do j = 1,4
      nd(mat_elem(e1,j)) = 1
   end do
end do
 ct = sum(nd)
 nd2 = nd
 nd  = 0
 m   = 0
do i = 1,no1
   if( (nd2(i) == 1) )then
      m = m + 1
      nd(m) = i
   end if
end do
!--
 cnd = 0.d0
do i = 1,ct
   cnd(i,1) = mat_coord(nd(i),1) ; cnd(i,2) = mat_coord(nd(i),2) ; cnd(i,3) = mat_coord(nd(i),3) ;
end do
!--
!if(my_rank.eq.0) print*,'  '
!if(my_rank.eq.0) print*,'Numero de nos nas heterogeneidades', ct
!if(my_rank.eq.0) print*,'  '
!--
NInh = ct;
!--
 close(30)
!--
allocate(Inh(no1))
!--
Inh = 0;
do i = 1,ct
   Inh(nd(i))=i
end do
!--
if(n_p.gt.1)call MPI_BARRIER(MPI_COMM_WORLD,err)

if(activ.eq.1)then
!--
   if(my_rank.eq.0)then
     t1 = MPI_WTIME()
     print*,'     '
     print*,'Determinando os campos primarios' 
     print*,'     '
   end if
!--
   open(350,file='./primarios/reEx_s.bin',status='replace',action='write', access='direct', &
 form='unformatted', recl=8   )
   open(360,file='./primarios/imEx_s.bin',status='replace',action='write', access='direct', &
 form='unformatted', recl=8   )
   open(370,file='./primarios/reEy_s.bin',status='replace',action='write', access='direct', &
 form='unformatted', recl=8   )
   open(380,file='./primarios/imEy_s.bin',status='replace',action='write', access='direct', &
 form='unformatted', recl=8   )
   open(390,file='./primarios/reEz_s.bin',status='replace',action='write', access='direct', &
 form='unformatted', recl=8   )
   open(3100,file='./primarios/imEz_s.bin',status='replace',action='write', access='direct', &
 form='unformatted', recl=8  )
!--
  !=-=-=-=
  call para_range( 1,nfreqs,n_p,my_rank,fri,frf )
  !=-=-=-=
!--
  do i = fri,frf
	cntt = 0;
	do m = 1,num_perf

	       !=========
!		delx = COORCP( sum(NRPP(1:m-1))+1 ,1) - COORCP( sum(NRPP(1:m-1))+2 ,1)
!		dely = COORCP( sum(NRPP(1:m-1))+1 ,2) - COORCP( sum(NRPP(1:m-1))+2 ,2)
!		ttap = datan(dely/delx); !print*,'tta:',ttap
	       !=========

		  do j = 1,NTPP(m)
		           cntt = cntt + 1; itrm = cntt
        		   print'(1x,a,1x,I5,1x,a,1x,I5,1x,a,1x,I5)','frequencia',i,'perfil',m,'transmissor',j
				do l = 1,ct

				        tr_x = COOTXP(cntt,1) ; tr_y = COOTXP(cntt,2) ; tr_z = COOTXP(cntt,3)

					call prim_MCSEM( ttap, Freqs(i) , tr_x , tr_y , tr_z , & 
                                  			 cnd( l , 1 ),cnd( l , 2 ),cnd( l , 3 ), & 
                                  			 Ex, Ey, Ez )

					if(m.eq.1)then
						k = (i-1)*sum( NTPP(:) )*ct + (j-1)*ct + l;
					else
						k = (i-1)*sum( NTPP(:) )*ct + sum( NTPP(1:m-1) )*ct + (j-1)*ct + l;
					end if

			                write(350,rec=k) dreal(Ex)
			                write(370,rec=k) dreal(Ey) 
			                write(390,rec=k) dreal(Ez) 
	
			                write(360,rec=k) dimag(Ex) 
			                write(380,rec=k) dimag(Ey) 
			                write(3100,rec=k) dimag(Ez) 

				end do
		  end do ! Trnasmissor


	end do	! Perfil

  end do ! Frequencia
!--
  close(350) ; close(360) 
  close(370) ; close(380) ; close(390) 
  close(3100)

  if(n_p.gt.1)call MPI_BARRIER(MPI_COMM_WORLD,err)

  if(my_rank.eq.0)then
     t2 = MPI_WTIME()
     !call cpu_time(t2) 
     print*,'     '
     print*,'Tempo para calcular os campos primarios',(t2-t1)/(60.d0),'minutos'
     print*,'     '
  end if

!stop

else

  if(my_rank.eq.0)then
    print*,'     '
    print*,'Campos primarios obtidos.' 
    print*,'     '
  end if

end if

end subroutine

!--

subroutine primarios(ncp,E,H)
implicit none
integer,intent(in) :: ncp
complex(dpc),intent(out) :: E(ncp,3),H(ncp,3)

integer   :: i,j,k,p
real(dpc) :: r,x,y,z
real(dpc) :: xt,yt,zt,f

integer :: ncpf,it,ix,ii,iat
real(dpc),allocatable :: off(:)
complex(dpc) :: Ex,Ey,Ez, E1, E2, E3

!open(unit = 20, file ='campo_E_aresta1D.dat', status = 'replace', action = 'write' )

E = 0.d0

ncpf = sum(NTPP(:)*NRPP(:))

ii=0
do p = 1,nfreqs
 iat = 0
 f = Freqs(p)  
 do i = 1,Num_perf
   do j = 1,NTPP(i)

      iat = iat + 1 ; itrm = iat

      if(i==1)then
        it = j
      else
        it = sum(NTPP(1:i-1)) + j
      end if

      xt = COOTXP(it,1) ; yt = COOTXP(it,2) ; zt = COOTXP(it,3); 

      do k = 1,NRPP(i)

         if(i==1)then
           ix = k
         else
           ix = sum(NRPP(1:i-1)) + k
         end if

         x = COORCP(ix,1); y = COORCP(ix,2); z = COORCP(ix,3);          
         r = dsqrt( (x-xt)**2 + (y-yt)**2 + (z-zt)**2)

         ii = ii + 1

!         rhop(3)=1.d0
!         call prim_MCSEM( 0.d0,f,xt,yt,zt,x,y,z,E1,E2,E3 )
!         rhop(3)=100.d0

         call prim_MCSEM( 0.d0,f,xt,yt,zt,x,y,z,Ex,Ey,Ez )

!         Ex = Ex - E1
!         Ey = Ey - E2
!         Ez = Ez - E3
!         write(20,'(6E20.10E3)')dreal(Ex),dimag(Ex),dreal(Ey),dimag(Ey),dreal(Ez),dimag(Ez)

         E(ii,1) = Ex; E(ii,2) = Ey; E(ii,3) = Ez
         H(ii,1) = Hpx; H(ii,2) = Hpy; H(ii,3) = Hpz

      end do
   end do
 end do
end do

!stop

! close (20)

end subroutine

!--

end module
