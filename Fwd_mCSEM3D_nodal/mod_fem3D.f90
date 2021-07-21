module fem_3D

use var_glob
use Arrays_SLU
use deriv_MQMP4
use DEHx
use arrayspds

contains

!%%%===========================================%%%

subroutine prim_MCSEM( tta,f,tx,ty,tz,x,y,z,Ex,Ey,Ez )
implicit none

real(dp),intent(in) :: f,tx,ty,tz,x,y,z,tta
complex(dp),intent(inout) :: Ex,Ey,Ez

real(dp), parameter :: mu = 4.d0*pi*1d-7
real(dp), parameter :: eps = 8.85d0*1.d-12
real(dp), parameter :: Iw = 1.d0
complex(dp) :: Hx,Hy,Hz
real(dp)    :: omega , h0 , dsx
complex*16  :: zeta, neta0
real(8), dimension(:), allocatable :: sigmas, h
integer :: ncam , i 

real(8) :: xr,yr
real(8) :: Txr,Tyr

ncam = ncm

allocate( sigmas( ncam ), h( ncam ) )

do i = 1,ncam-1
 h(i) = hjp(i)
end do
h( ncam ) = 1.d300
do i = 1,ncam
 sigmas(i) = 1.d0/rhop(i)
end do
omega = 2.d0 * pi * f
neta0 = 1.d0/rho_ar
zeta = (0.d0, 1.d0) * omega * mu
h0 = Tz
dsx = mdip

!==============

Txr=Tx ; Tyr=Ty;
xr=x   ; yr = y;

call dehx_xyz_loops( Txr, Tyr, Iw, dsx, h0, ncam, h, sigmas, neta0, zeta, xr, yr, z, Ex, Ey, Ez, Hx, Hy, Hz )

!==============

Hpx = Hx
Hpy = Hy
Hpz = Hz

deallocate( sigmas, h )

end subroutine

!%%%===========================================%%%

subroutine Atualiza_Parametros(npr,P)
implicit none
integer,intent(in) :: npr
real(dpc),intent(in) :: P(npr)

integer :: l,m

vprop(1) = rho_ar
m = 1
do l = 1,ncm
  m = m + 1
  vprop(m) = rhop(l)
end do
do l = 1, npr
    m = m + 1
    vprop( m ) = P(l)
end do

end subroutine

!%%%===========================================%%%

subroutine Atualiza_Prop_malha
implicit none
integer :: k,l,erro

do k = 1, nelemg
    do l = 1, nprop
        if( vidprop(k)==vflprop(l) )then
          rho_elem(k) = vprop(l) ; exit
        end if
    end do
end do

end subroutine

!%%%===========================================%%%

subroutine Montagem_Matriz_Global_tetraedro( ngl , nnos , nel , n_nz , vet_nn )
implicit none
integer,intent(in) :: nel, n_nz , ngl , nnos
complex(dpc),intent(out) :: vet_nn(n_nz)

integer :: ie,i,j,k,l,m
real(dpc),allocatable :: xe(:),ye(:),ze(:)
complex(dpc),allocatable :: mloc(:,:)
real(dpc),allocatable :: a(:),b(:),c(:),d(:)
real(dpc),allocatable :: xa(:),ya(:),za(:) 
integer,allocatable   :: lnode(:)

real(dpc)                :: vol , rho_e , sig_e , omega , zm
complex(dpc),allocatable :: M_base(:,:) , Phi_xyz(:,:)
complex(dpc),allocatable :: Phi_x(:,:) , Phi_y(:,:) , Phi_z(:,:)  

allocate(lnode(16))
allocate(Phi_xyz(4,4),M_base(4,4))
allocate(Phi_x(4,4),Phi_y(4,4),Phi_z(4,4))
allocate(xa(3),ya(3),za(3))
allocate(xe(4),ye(4),ze(4))
allocate(a(4),b(4),c(4),d(4))

vet_nn = (0.d0,0.d0)

allocate(mloc(16,16)) ; 

omega = 2*pi*Freq ;

do ie = 1,nel

 xe = mat_coord(mat_elem(ie,:),1);
 ye = mat_coord(mat_elem(ie,:),2);
 ze = mat_coord(mat_elem(ie,:),3);

 lnode =[ mat_elem(ie,1)*4-3,mat_elem(ie,1)*4-2,mat_elem(ie,1)*4-1,mat_elem(ie,1)*4,&
	  mat_elem(ie,2)*4-3,mat_elem(ie,2)*4-2,mat_elem(ie,2)*4-1,mat_elem(ie,2)*4,&
	  mat_elem(ie,3)*4-3,mat_elem(ie,3)*4-2,mat_elem(ie,3)*4-1,mat_elem(ie,3)*4,&
          mat_elem(ie,4)*4-3,mat_elem(ie,4)*4-2,mat_elem(ie,4)*4-1,mat_elem(ie,4)*4 ]


!%%%%=====%%%

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

!%%%%=====%%%

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

!%%%%=====%%%

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

!%%%%=====%%%

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

!%%%%=====%%%

rho_e = rho_elem(ie);
vol   = dabs(vol_tetraedro(xe,ye,ze));
sig_e = 1.d0/rho_e;
mloc  = (0.d0,0.d0);
Phi_x = (0.d0,0.d0); Phi_y = (0.d0,0.d0); Phi_z = (0.d0,0.d0); 

a = a/(6*vol) ; b = b/(6*vol) ; c = c/(6*vol) ; d = d/(6*vol)

do k = 1,4
	   do l = 1,4

	      M_base(k,l) = vol * ( b(k) * b(l) + c(k) * c(l) + d(k) * d(l) ) + &
                            ip * omega * mi0 * sig_e * ( vol/20.d0 ) * ( 1.d0 + delta(k,l) )

	      Phi_x(k,l)   = ip * omega * mi0 * sig_e * vol * (b(l)/4.d0) 
	      Phi_y(k,l)   = ip * omega * mi0 * sig_e * vol * (c(l)/4.d0)
	      Phi_z(k,l)   = ip * omega * mi0 * sig_e * vol * (d(l)/4.d0)
	
              Phi_xyz(k,l) = ip * omega * mi0 * sig_e * vol * ( b(k) * b(l) + c(k) * c(l) + d(k) * d(l) )

	   end do
end do

do k = 1, 4
	do l = 1, 4

              mloc( 4*k-3,4*l-3 )  = M_base(k,l)
              mloc( 4*k-2,4*l-2 )  = M_base(k,l)
              mloc( 4*k-1,4*l-1 )  = M_base(k,l)

              mloc( 4*k-3,4*l )    = Phi_x(k,l)
              mloc( 4*k-2,4*l )    = Phi_y(k,l)
              mloc( 4*k-1,4*l )    = Phi_z(k,l)

              mloc( 4*k,4*l-3 )    = Phi_x(l,k)
              mloc( 4*k,4*l-2 )    = Phi_y(l,k)
              mloc( 4*k,4*l-1 )    = Phi_z(l,k)

              mloc( 4*k,4*l )      = Phi_xyz(k,l)

	end do
end do

!-- chamar aquia a montagem da matriz global %%%%
call set_arrays_pds(cdim,ngl,nnosg,nnz,row_inda,lnode,mloc,col_ptra,val_nz)
!--

!================================================

end do ! loop elemosntos

deallocate(lnode)
deallocate(mloc)
deallocate( Phi_xyz , M_base )
deallocate( Phi_x , Phi_y , Phi_z )
deallocate(xa,ya,za)
deallocate(xe,ye,ze)
deallocate(a,b,c,d)

end subroutine

!%%%===========================================%%%

subroutine Montagem_Vetor_Fonte_tetraedro( ngl , nel , nnode , vetF )
implicit none
integer,intent(in) :: ngl , nnode , nel
complex(dpc),intent(out) :: vetF(:)

integer :: ie,i,j,k,l,m
real(dpc),allocatable    :: xe(:),ye(:),ze(:)
complex(dpc),allocatable :: floc(:)
real(dpc),allocatable :: a(:),b(:),c(:),d(:)
real(dpc),allocatable :: xa(:),ya(:),za(:) 
integer,allocatable   :: lnode(:)

real(dpc)                :: vol , rho_e , sig_e , sig_p , del_sig , omega , zm
complex(dpc),allocatable :: Ep_x(:) , Ep_y(:) , Ep_z(:)  
integer                  :: cmel  

integer                  :: k1
real(dpc)                :: rEx,rEy,rEz,iEx,iEy,iEz

allocate(lnode(16))
allocate(Ep_x(4),Ep_y(4),Ep_z(4))
allocate(xa(3),ya(3),za(3))
allocate(xe(4),ye(4),ze(4))
allocate(a(4),b(4),c(4),d(4))

!vetF = (0.d0,0.d0);

allocate(floc(16)) ; 

do ie = 1,nel

 xe = mat_coord(mat_elem(ie,:),1);
 ye = mat_coord(mat_elem(ie,:),2);
 ze = mat_coord(mat_elem(ie,:),3);

 lnode =[ mat_elem(ie,1)*4-3,mat_elem(ie,1)*4-2,mat_elem(ie,1)*4-1,mat_elem(ie,1)*4,&
	  mat_elem(ie,2)*4-3,mat_elem(ie,2)*4-2,mat_elem(ie,2)*4-1,mat_elem(ie,2)*4,&
	  mat_elem(ie,3)*4-3,mat_elem(ie,3)*4-2,mat_elem(ie,3)*4-1,mat_elem(ie,3)*4,&
          mat_elem(ie,4)*4-3,mat_elem(ie,4)*4-2,mat_elem(ie,4)*4-1,mat_elem(ie,4)*4 ]

!%%%%=====%%%

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

!%%%%=====%%%

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

!%%%%=====%%%

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

!%%%%=====%%%

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

!%%%%=====%%%

rho_e = rho_elem(ie) ; 
vol   = dabs(vol_tetraedro (xe,ye,ze));
sig_e = 1.d0/rho_e;
floc  = (0.d0,0.d0);

a = a/(6*vol) ; b = b/(6*vol) ; c = c/(6*vol) ; d = d/(6*vol)

zm = sum(ze)/size(ze);

if( zm < z_int(1) )then 
 sig_p = (1.d0/rho_ar) ; goto 30
elseif( zm > z_int(ncm) )then 
 cmel = ncm 
 sig_p = 1.d0/rhop(cmel)
 goto 30
elseif( zm > z_int(1) .and. zm < z_int(ncm))then
 do i = 2,ncm
     if( zm > z_int(i-1) .and. zm < z_int(i)  )then
       cmel = i-1; exit
     end if 
 end do
 sig_p = 1.d0/rhop(cmel)
end if

30 continue

del_sig = sig_e - sig_p;

if( dabs(del_sig) .gt. 1.d-8 )then

!print*,1.d0/sig_e , 1.d0/sig_p , zm

do k = 1,4

 if(uprm.eq.0)then
    call prim_MCSEM( ttap, Freq , tr_x , tr_y , tr_z , & 
                   xe(k),ye(k),ze(k), & 
                   Ep_x(k),Ep_y(k),Ep_z(k) )
 else

    k1 = (ifrq-1)*sum(NTPP)*NInh + (itrm-1)*NInh + Inh( mat_elem(ie,k) )			 			 

    read(1111,rec=k1)rEx
    read(1112,rec=k1)iEx

    read(1113,rec=k1)rEy
    read(1114,rec=k1)iEy

    read(1115,rec=k1)rEz
    read(1116,rec=k1)iEz
        
    Ep_x(k) = dcmplx(rEx,iEx)
    Ep_y(k) = dcmplx(rEy,iEy)
    Ep_z(k) = dcmplx(rEz,iEz)

 end if

end do

do k = 1,4
!%%%%=======%%%%%%%%%%

  floc(4*k-3)  = del_sig * (vol/20.d0) * ( Ep_x(1)*(1 + delta(k,1)) + Ep_x(2)*(1 + delta(k,2)) + &
 Ep_x(3)*(1 + delta(k,3)) + Ep_x(4)*(1 + delta(k,4)) )
  floc(4*k-2)  = del_sig * (vol/20.d0) * ( Ep_y(1)*(1 + delta(k,1)) + Ep_y(2)*(1 + delta(k,2)) + &
 Ep_y(3)*(1 + delta(k,3)) + Ep_y(4)*(1 + delta(k,4)) )
  floc(4*k-1)  = del_sig * (vol/20.d0) * ( Ep_z(1)*(1 + delta(k,1)) + Ep_z(2)*(1 + delta(k,2)) + &
 Ep_z(3)*(1 + delta(k,3)) + Ep_z(4)*(1 + delta(k,4)) )

  floc(4*k)    =  del_sig * vol * (1.d0/4.d0) * ( ( Ep_x(1) + Ep_x(2) + Ep_x(3) + Ep_x(4) )*b(k)  &
               + ( Ep_y(1) + Ep_y(2) + Ep_y(3) + Ep_y(4) )*c(k) &
               + ( Ep_z(1) + Ep_z(2) + Ep_z(3) + Ep_z(4) )*d(k) )

!%%%%=======%%%%%%%%%%
end do

else

  floc = (0.d0,0.d0) ;

end if

!%%% chamar aquia a montagem do vetor fonte global %%%%

do l = 1,16
    vetF( lnode(l) + (irs-1)*ngb ) = vetF( lnode(l) + (irs-1)*ngb ) + floc( l )
end do

!================================================

end do ! loop elemosntos

deallocate(lnode)
deallocate(floc)
deallocate(Ep_x,Ep_y,Ep_z)
deallocate(xa,ya,za)
deallocate(xe,ye,ze)
deallocate(a,b,c,d)

end subroutine

!%%%===========================================%%%

real(dpc) function vol_tetraedro ( x,y,z )
implicit none
real(dpc),intent(in) :: x(4),y(4),z(4)

vol_tetraedro = (1.d0/6.d0) * ( xij(x(2),x(1))*( yij(y(2),y(3))*zij(z(3),z(4)) - yij(y(3),y(4))*zij(z(2),z(3)) ) + &
                            xij(x(3),x(2))*( yij(y(3),y(4))*zij(z(1),z(2)) - yij(y(1),y(2))*zij(z(3),z(4)) ) + & 
                            xij(x(4),x(3))*( yij(y(1),y(2))*zij(z(2),z(3)) - yij(y(2),y(3))*zij(z(1),z(2)) ) )

end function

!%%%===========================================%%%

real(dpc) function det ( x , y , z )
implicit none
real(dpc),intent(in) :: x(3),y(3),z(3)

det = x(1)*( y(2)*z(3) - y(3)*z(2) ) + y(1)*( z(2)*x(3) - z(3)*x(2) ) + z(1)*( x(2)*y(3) - x(3)*y(2) )

end function

!%%%===========================================%%%

real(dpc) function xij(xi,xj)
implicit none
real(dpc),intent(in) :: xi,xj
  xij = xi - xj
end function
real(dpc) function yij(yi,yj)
implicit none
real(dpc),intent(in) :: yi,yj
  yij = yi - yj
end function
real(dpc) function zij(zi,zj)
implicit none
real(dpc),intent(in) :: zi,zj
  zij = zi - zj
end function

!%%%===========================================%%%

integer function delta(i,j)
implicit none
integer,intent(in) :: i,j

if(i.eq.j)then
    delta = 1;
else
    delta = 0;
end if

end function

!%%%===========================================%%%

subroutine aloc_opr (ng,nr)
implicit none
integer,intent(in) :: ng,nr

allocate(Qp(ng,nr),Hp(ng,nr))
allocate(sum_jaca(nprm,nr))

end subroutine 

!%%%===========================================%%%

subroutine desal_opr
implicit none

deallocate(Qp,Hp)
deallocate(sum_jaca)

end subroutine 

!%%%===========================================%%%


!--

end module
