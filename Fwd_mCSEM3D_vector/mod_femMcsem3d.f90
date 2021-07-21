module fem_mcsem3d

use vglob
use dehx
use arrayspds
use Arrays_SLU
use lstm

implicit none
contains

!'''''''''''''''''''''''''''''''''

!-- Campo primario

!'''''''''''''''''''''''''''''''''

!--

subroutine prim_MCSEM( tta,f,tx,ty,tz,x,y,z,rEx,rEy,Ez )
implicit none

real(db),intent(in) :: f,tx,ty,tz,x,y,z,tta
complex(db),intent(inout) :: rEx,rEy,Ez

real(db), parameter :: mu = 4.d0*pi*1d-7
real(db), parameter :: eps = 8.85d0*1.d-12
real(db), parameter :: Iw = 1.d0
complex(db) :: Hx,Hy,Hz,Ex,Ey
real(db)    :: omega , h0 , dsx
complex*16  :: zeta, neta0
real(db), dimension(:), allocatable :: sigmas, h
integer :: ncam , i 

real(db) :: xr,yr
real(db) :: Txr,Tyr

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

Txr = tx  
Tyr = ty 

xr = x  
yr = y 

call dehx_xyz_loops( Txr, Tyr, Iw, dsx, h0, ncam, h, sigmas, neta0, zeta, xr, yr, z, Ex, Ey, Ez, Hx, Hy, Hz )

 rEx = Ex 
 rEy = Ey 

!==============
 cpHx = Hx 
 cpHy = Hy
 cpHz = Hz
!==============

deallocate( sigmas, h )

end subroutine

!--

subroutine open_arqprim (  )
implicit none
  
open(1111,file='./primarios/reEx_s.bin',status='old',action='read', access='direct', form='unformatted', recl=lendt  )
open(1112,file='./primarios/imEx_s.bin',status='old',action='read', access='direct', form='unformatted', recl=lendt  )
open(1113,file='./primarios/reEy_s.bin',status='old',action='read', access='direct', form='unformatted', recl=lendt  )
open(1114,file='./primarios/imEy_s.bin',status='old',action='read', access='direct', form='unformatted', recl=lendt  )
open(1115,file='./primarios/reEz_s.bin',status='old',action='read', access='direct', form='unformatted', recl=lendt  )
open(1116,file='./primarios/imEz_s.bin',status='old',action='read', access='direct', form='unformatted', recl=lendt  )

end subroutine

!--

subroutine close_arqprim
implicit none

 close(1111);close(1112);close(1113)
 close(1114);close(1115);close(1116)

end subroutine

!--

subroutine primarios_modelo(activ,np,nem1,ihet1,Inh,NInh,n_p,my_rank)
implicit none
include 'mpif.h'

 integer,intent(in) :: n_p, my_rank, activ, np
 integer,intent(out)  :: nem1 , NInh
 integer,allocatable,intent(out) :: ihet1(:),Inh(:)

 character(len=100) :: pathprim
 integer :: lenc
 integer :: no1,fri,frf,err,aux
 integer,allocatable :: nd(:)
 real(db),allocatable :: cnd(:,:)
 integer :: i,j,k,l,m,ct,noi
 real(db) :: x,y,z 
 complex(db) :: Ex,Ey,Ez,Haux
 complex(db) :: rEx,rEy
 real(db) :: t1,t2
 integer :: ne1, nel, e1, lb1, ki
 integer,allocatable :: vlb(:),inde(:,:),nd2(:)
 integer :: cntt

 real(db) :: delx, dely

 integer  :: i1, i2  

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

!-- Aqui e verificado quais elementos estao nas heterogeneidades --

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
nem1=nel;
 close(10)
allocate( ihet1(nel) )
do i = 1, nel
     ihet1(i) = vlb(i) 
end do
deallocate(vlb)

!--

open(30,file='./malha/'//trim(mname)//'.1.edge',status='old',action='read')

read(30,*)no1

allocate(nd(no1),cnd(no1,3)) ; allocate(inde(no1,2))
allocate(nd2(no1))

!if(my_rank.eq.0) print*,'  '
!if(my_rank.eq.0) print*,'Numero de arestas na malha:',no1
!-- Aqui quais arestas estao nas heterogeneidades --
!--
do i = 1,no1
   read(30,*)noi,inde(i,1),inde(i,2)
end do
!--
nd = 0
do i = 1,nem1
   e1 = ihet1(i)
   do j = 1,6
      nd(gEdge(e1,j)) = 1
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
do i = 1,ct
   i1 = inde(nd(i),1) ; i2 = inde(nd(i),2)
   x = (mat_coord(i1,1) + mat_coord(i2,1))/2.d0
   y = (mat_coord(i1,2) + mat_coord(i2,2))/2.d0
   z = (mat_coord(i1,3) + mat_coord(i2,3))/2.d0
   cnd(i,1) = x ; cnd(i,2) = y ; cnd(i,3) = z ;
end do
!--
!if(my_rank.eq.0) print*,'  '
!if(my_rank.eq.0) print*,'Numero de arestas nas heterogeneidades:',ct
!if(my_rank.eq.0) print*,'  '

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
    print*,'Determinando os campos primarios...' 
    print*,'     '
 end if

   open(350,file='./primarios/reEx_s.bin',status='replace',action='write', access='direct', &
 form='unformatted', recl=lendt)
   open(360,file='./primarios/imEx_s.bin',status='replace',action='write', access='direct', &
 form='unformatted', recl=lendt)
   open(370,file='./primarios/reEy_s.bin',status='replace',action='write', access='direct', &
 form='unformatted', recl=lendt)
   open(380,file='./primarios/imEy_s.bin',status='replace',action='write', access='direct', &
 form='unformatted', recl=lendt)
   open(390,file='./primarios/reEz_s.bin',status='replace',action='write', access='direct', &
 form='unformatted', recl=lendt)
   open(3100,file='./primarios/imEz_s.bin',status='replace',action='write', access='direct', &
 form='unformatted', recl=lendt)

  !=-=-=-=
  call para_range( 1,nfreqs,n_p,my_rank,fri,frf )
  !=-=-=-=

  do i = fri,frf
	cntt = 0;
	do m = 1,num_perf

		  do j = 1,NTPP(m)
		           cntt = cntt + 1; itrm = cntt;
        		   print'(1x,a,1x,I5,1x,a,1x,I5,1x,a,1x,I5)','frequencia',i,'perfil',m,'transmissor',j
				do l = 1,ct

                                       tr_x = COOTXP(cntt,1) 
                                       tr_y = COOTXP(cntt,2) 
                                       tr_z = COOTXP(cntt,3)

                                       x = cnd( l , 1 ) 
                                       y = cnd( l , 2 )

					call prim_MCSEM( ttap, Freqs(i) , tr_x , tr_y , tr_z , & 
                                  			 x,y,cnd( l , 3 ), & 
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
     print*,'     '
     print*,'Tempo para calcular os campos primarios',(t2-t1)/(60.d0),'minutos'
     print*,'     '
  end if

else

  if(my_rank.eq.0)then
    print*,'     '
    print*,'Campos primarios obtidos.' 
    print*,'     '
  end if

end if

deallocate(nd,cnd)

end subroutine

!--

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

!--

subroutine Atualiza_Parametros(npr,P)
implicit none
integer,intent(in) :: npr
real(db),intent(in) :: P(npr)

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

!--

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

!-- 

subroutine compose_Esec ( xr, yr, zr, Ex, Ey, Ez, Et )
implicit none
!integer,intent(in)   :: no
real(db),intent(in)  :: xr, yr, zr
complex(db),intent(in)  :: Et(:)
complex(db),intent(out) :: Ex,Ey,Ez

integer :: elems(5000)
integer :: i,j,k,ce,M,ie

real(db)             :: dx,dy,dz,leng
real(db)             :: vol,Nx,Ny,Nz,xm,ym,zm
real(db),allocatable :: xe(:),ye(:),ze(:)
real(db),allocatable :: xa(:),ya(:),za(:) 
real(db),allocatable :: a(:),b(:),c(:),d(:),dists(:)
integer              :: i1,i2,j1,j2,aux
complex(db) :: sx,sy,sz
integer                 :: redg(nedgs,2), si

real(db) :: Pr, Pce, drce, dr

!--
complex(db),allocatable    :: fedD(:,:)
complex(db)                :: f(10), ca(10)
complex(db),allocatable    :: mA(:,:), u(:), y(:), cW(:,:), Ac(:,:), Bc(:)
complex(db),allocatable    :: yy(:), yz(:)
complex(db),allocatable    :: aTw(:,:)
integer                    :: ii, Ncs = 10
real(db)                   :: cc = 1.d0
real(db),allocatable       :: cxe(:), cye(:), cze(:)
integer                    :: ndga
integer,allocatable        :: iedg(:)
real(db)                   :: auxx, auxy, auxz
complex(db)                :: auxf(3) 

integer,allocatable        :: iedga(:)
real(db),allocatable       :: xea(:),yea(:),zea(:)
complex(db),allocatable    :: Fda(:,:)

dr = rad

!--

allocate(dists(nelemg))
allocate(xa(3),ya(3),za(3))
allocate(xe(nnoele),ye(nnoele),ze(nnoele))
allocate(a(nnoele),b(nnoele),c(nnoele),d(nnoele))

!-- procura em qual elemento esta o no, pela minima distancia ao centro do mesmo...

do i = 1,nelemg

   xe = mat_coord(mat_elem(i,:),1);
   ye = mat_coord(mat_elem(i,:),2);
   ze = mat_coord(mat_elem(i,:),3);
   xm = sum(xe)/nnoele
   ym = sum(ye)/nnoele
   zm = sum(ze)/nnoele

   drce  = dsqrt( (xr-xm)**2 + (yr-ym)**2 + (zr-zm)**2 )
   dists(i) = drce

end do

!-- Para o caso de um elemento.
!M = 1 ! Apenas um elemento, o de distancia minima
!elems(M) = minloc(dists,dim=1)

!-- Para mais elementos, dentro do raio estabelecido.

M=0
do i = 1,nelemg
   if(dists(i)<=dr)then
      M = M + 1
      elems(M) = i
   end if
end do

!print*,'M=',M !; stop

!ndga = nedgs ; nedgs = 1

allocate(fedD(M*nedgs,3),iedg(M*nedgs))
allocate(cxe(nedgs*M), cye(nedgs*M), cze(nedgs*M))
!--

!Ex = 0.d0
!Ey = 0.d0
!Ez = 0.d0

!-- Relacao n\'o aresta, localmente.

redg(1,1) = 1 ; redg(1,2) = 2
redg(2,1) = 1 ; redg(2,2) = 3
redg(3,1) = 1 ; redg(3,2) = 4
redg(4,1) = 2 ; redg(4,2) = 3
redg(5,1) = 4 ; redg(5,2) = 2
redg(6,1) = 3 ; redg(6,2) = 4

!--

!--
ii = 0
!--

do i = 1,M

 ie = elems(i)

 xe = mat_coord(mat_elem(ie,:),1);
 ye = mat_coord(mat_elem(ie,:),2);
 ze = mat_coord(mat_elem(ie,:),3);

!-- Calculando as coordenadas dos centros dos elementos.
 !xm = sum(xe)/size(xe);
 !ym = sum(ye)/size(ye);
 !zm = sum(ze)/size(ze);
!--
 !xm = xr
 !ym = yr
 !zm = zr
 !print*,xm,ym,zm,'coordenadas'
!--

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

 vol   = dabs(vol_tetraedro(xe,ye,ze));

!--
!Ex = 0.d0
!Ey = 0.d0
!Ez = 0.d0

do k = 1,nedgs

 !-- Calculando as coordenadas dos centros das arestas
 i1 = redg(k,1); i2 = redg(k,2);
 xm  = (mat_coord(mat_elem(ie,i2),1) + mat_coord(mat_elem(ie,i1),1))/2.d0 
 ym  = (mat_coord(mat_elem(ie,i2),2) + mat_coord(mat_elem(ie,i1),2))/2.d0 
 zm  = (mat_coord(mat_elem(ie,i2),3) + mat_coord(mat_elem(ie,i1),3))/2.d0 
 !--

 ii = ii + 1
 cxe(ii)=xm; cye(ii)=ym; cze(ii)=zm

 iedg(ii) = gEdge(ie,k)

 sx = 0.; sy = 0.; sz = 0.;
 do j = 1,nedgs
 !do j = 1,nedga

    i1 = redg(j,1); i2 = redg(j,2);

    dx = mat_coord(mat_elem(ie,i2),1) - mat_coord(mat_elem(ie,i1),1)
    dy = mat_coord(mat_elem(ie,i2),2) - mat_coord(mat_elem(ie,i1),2)
    dz = mat_coord(mat_elem(ie,i2),3) - mat_coord(mat_elem(ie,i1),3)

    si = (mat_elem(ie,i2) - mat_elem(ie,i1))/abs( mat_elem(ie,i2) - mat_elem(ie,i1) )

    leng = si * dsqrt( dx**2 + dy**2 + dz**2 )

    Nx    = (leng / (6 * vol)) * ( (fnodal( vol,i1,xm,ym,zm,a,b,c,d ) * b(i2))-(fnodal( vol,i2,xm,ym,zm,a,b,c,d ) * b(i1)) )

    Ny    = (leng / (6 * vol)) * ( (fnodal( vol,i1,xm,ym,zm,a,b,c,d ) * c(i2))-(fnodal( vol,i2,xm,ym,zm,a,b,c,d ) * c(i1)) )

    Nz    = (leng / (6 * vol)) * ( (fnodal( vol,i1,xm,ym,zm,a,b,c,d ) * d(i2))-(fnodal( vol,i2,xm,ym,zm,a,b,c,d ) * d(i1)) )

    sx = sx + Nx * Et( gEdge(ie,j) )
    sy = sy + Ny * Et( gEdge(ie,j) ) 
    sz = sz + Nz * Et( gEdge(ie,j) ) 
    
 end do

!--
! Ex = Ex + sx
! Ey = Ey + sy
! Ez = Ez + sz
!--
 fedD(ii,1) = sx
 fedD(ii,2) = sy
 fedD(ii,3) = sz
!--

end do
!--

end do

!Ex = Ex/M
!Ey = Ey/M
!Ez = Ez/M

!-- Recontagem dos valores de campo

!-- reordenacao

do i = 1 , M*nedgs-1
       do j = 1 , M*nedgs-i 
          if(iedg(j) > iedg(j+1))then

           aux = iedg(j)
           iedg(j) = iedg(j+1)
           iedg(j+1) = aux            

           auxx = cxe(j)
           cxe(j) = cxe(j+1)
           cxe(j+1) = auxx            

           auxy = cye(j)
           cye(j) = cye(j+1)
           cye(j+1) = auxy            

           auxz = cze(j)
           cze(j) = cze(j+1)
           cze(j+1) = auxz            

           auxf = fedD(j,:)
           fedD(j,:) = fedD(j+1,:)
           fedD(j+1,:) = auxf            

          end if
       end do 
end do

!-- recontagem

allocate(iedga(nedgs*M),xea(nedgs*M),yea(nedgs*M),zea(nedgs*M))
allocate(fda(nedgs*M,3))

ii=0
do i = 1,nedgs*M
   if(i == 1)then
      ii=ii+1
      iedga(ii) = iedg(i)
      xea(ii) = cxe(i)
      yea(ii) = cye(i)
      zea(ii) = cze(i)
      fda(ii,:) = FedD(i,:)
   else if( i < nedgs*M )then
      if(iedg(i) < iedg(i+1))then
         ii=ii+1
         iedga(ii) = iedg(i+1)
         xea(ii) = cxe(i+1)
         yea(ii) = cye(i+1)
         zea(ii) = cze(i+1)
         fda(ii,:) = FedD(i+1,:)
      end if
      if(iedg(i+1) == maxval(iedg))exit
   else
      ii=ii+1
      iedga(ii) = iedg(i)
      xea(ii) = cxe(i)
      yea(ii) = cye(i)
      zea(ii) = cze(i)
      fda(ii,:) = FedD(i,:)
   end if
end do

M = ii

deallocate(cxe,cye,cze,FedD)

allocate(cxe(M),cye(M),cze(M),FedD(M,3))

do i = 1,M
   cxe(i) = xea(i)
   cye(i) = yea(i)
   cze(i) = zea(i)
   FedD(i,:) = fda(i,:)
end do

ndga = nedgs
nedgs = 1

!--

allocate( mA(nedgs*M,Ncs), u(nedgs*M), y(nedgs*M), cW(nedgs*M,nedgs*M) )
allocate( yy(nedgs*M),yz(nedgs*M) )

 cW = 0.d0
 ca = 0.d0
 ii=0
do i = 1,M
  ie = elems(i)
  do j = 1,nedgs
    !--
     !i1 = redg(j,1); i2 = redg(j,2);
     !xm  = (mat_coord(mat_elem(ie,i2),1) + mat_coord(mat_elem(ie,i1),1))/2.d0
     !ym  = (mat_coord(mat_elem(ie,i2),2) + mat_coord(mat_elem(ie,i1),2))/2.d0 
     !zm  = (mat_coord(mat_elem(ie,i2),3) + mat_coord(mat_elem(ie,i1),3))/2.d0 
    !-- 
     !xe = mat_coord(mat_elem(ie,:),1);
     !ye = mat_coord(mat_elem(ie,:),2);
     !ze = mat_coord(mat_elem(ie,:),3);
     !xm = sum(xe)/size(xe);
     !ym = sum(ye)/size(ye);
     !zm = sum(ze)/size(ze);
    !--
      xm = cxe(i)
      ym = cye(i)
      zm = cze(i)
    !--

     ii=ii+1
     call cf(Ncs,xm,ym,zm,f)

     mA(ii,:)  = f

     y(ii)     = fedD(ii,1)
     yy(ii)    = fedD(ii,2)
     yz(ii)    = fedD(ii,3)

     cW(ii,ii) = dexp( - cc**2 * ( ((xr-xm)/(xr-maxval(cxe)))**2 + ((yr-ym)/(yr-maxval(cye)))**2 &
                 + ((zr-zm)/(zr-maxval(cze)))**2 ) )
      
  end do
end do

allocate(Ac(Ncs,Ncs), Bc(Ncs))
allocate(aTw(Ncs,M))

do i = 1,Ncs
   do j = 1,M   
      aTw(i,j) = sum( mA(:,i)*cW(:,j) )
   end do
end do

!-- Ex

do i = 1,Ncs
   do j = 1,Ncs   
      Ac(i,j) = sum( aTw(i,:)*mA(:,j) )
   end do
   Bc(i) = sum( aTw(i,:)*y(:) )
end do
call LU_Gauss_compacto( Ac, Bc, Ncs)
call cf(Ncs,xr,yr,zr,f)

Ex = sum(f*Bc)

!-- Ey

do i = 1,Ncs
   do j = 1,Ncs   
      Ac(i,j) = sum( aTw(i,:)*mA(:,j) )
   end do
   Bc(i) = sum( aTw(i,:)*yy(:) )
end do
call LU_Gauss_compacto( Ac, Bc, Ncs)
call cf(Ncs,xr,yr,zr,f)

Ey = sum(f*Bc)

!-- Ez

do i = 1,Ncs
   do j = 1,Ncs   
      Ac(i,j) = sum( aTw(i,:)*mA(:,j) )
   end do
   Bc(i) = sum( aTw(i,:)*yz(:) )
end do
call LU_Gauss_compacto( Ac, Bc, Ncs)
call cf(Ncs,xr,yr,zr,f)

Ez = sum(f*Bc)

!--

!print*,cdabs(Ex),cdabs(Ey),cdabs(Ez)

nedgs = ndga

end subroutine

!--

subroutine compose_Esec3( xr, yr, zr, Ex, Ey, Ez, Hx, Hy, Hz, Et )

implicit none

real(db),intent(in)  :: xr, yr, zr
complex(db),intent(in)  :: Et(:)
complex(db),intent(out) :: Ex,Ey,Ez
complex(db),intent(out) :: Hx,Hy,Hz

integer :: elems(5000)
integer :: i,j,k,ce,M,ie

real(db)             :: dx,dy,dz,leng
real(db)             :: vol,Nx,Ny,Nz,xm,ym,zm
real(db),allocatable :: xe(:),ye(:),ze(:)
real(db),allocatable :: xa(:),ya(:),za(:) 
real(db),allocatable :: a(:),b(:),c(:),d(:),dists(:)
integer              :: i1,i2,j1,j2,aux
complex(db) :: sx,sy,sz
integer              :: redg(nedgs,2), si

real(db) :: Pr, Pce, drce, dr

!--
complex(db),allocatable    :: fedD(:,:)
complex(db)                :: f(10), ca(10)
complex(db),allocatable    :: mA(:,:), u(:), y(:), cW(:,:), Ac(:,:), Bc(:)
complex(db),allocatable    :: yy(:), yz(:)
complex(db),allocatable    :: aTw(:,:)
integer                    :: ii, Ncs = 10
real(db)                   :: cc = 1.d0
real(db),allocatable       :: cxe(:), cye(:), cze(:)
integer                    :: ndga
integer,allocatable        :: iedg(:)
real(db)                   :: auxx, auxy, auxz
complex(db)                :: auxf(3) 

integer,allocatable        :: iedga(:)
real(db),allocatable       :: xea(:),yea(:),zea(:)
complex(db),allocatable    :: Fda(:,:)

complex(db),allocatable    :: mx(:),my(:),mz(:)
complex(db),allocatable    :: dfdx(:),dfdy(:),dfdz(:)
real(db)                   :: fw

complex(db)                :: dExdz,dExdy
complex(db)                :: dEydz,dEydx
complex(db)                :: dEzdy,dEzdx
complex(db)                :: fma(10), fme(10)
real(db)                   :: perc = 1.5


!--

dr = rad

allocate(mx(Ncs),my(Ncs),mz(Ncs))
allocate(dfdx(Ncs),dfdy(Ncs),dfdz(Ncs))

allocate(dists(nelemg))
allocate(xa(3),ya(3),za(3))
allocate(xe(nnoele),ye(nnoele),ze(nnoele))
allocate(a(nnoele),b(nnoele),c(nnoele),d(nnoele))

!-- procura em qual elemento esta o no, pela minima distancia ao centro do mesmo...

do i = 1,nelemg

   xe = mat_coord(mat_elem(i,:),1);
   ye = mat_coord(mat_elem(i,:),2);
   ze = mat_coord(mat_elem(i,:),3);
   xm = sum(xe)/nnoele
   ym = sum(ye)/nnoele
   zm = sum(ze)/nnoele

   drce  = dsqrt((xr-xm)**2 + (yr-ym)**2 + (zr-zm)**2)
   dists(i) = drce

end do

!-- Para mais elementos, dentro do raio estabelecido.

M=0
do i = 1,nelemg
   if(dists(i)<=dr)then
      M = M + 1
      elems(M) = i
   end if
end do

!--
! print*,xr,yr,zr,'coordenadas'
!--

!print*,'M=',M !; stop

allocate(fedD(M,3),iedg(M))
allocate(cxe(M),cye(M),cze(M))
!--

!-- Relacao n\'o aresta, localmente.
redg(1,1) = 1 ; redg(1,2) = 2
redg(2,1) = 1 ; redg(2,2) = 3
redg(3,1) = 1 ; redg(3,2) = 4
redg(4,1) = 2 ; redg(4,2) = 3
redg(5,1) = 4 ; redg(5,2) = 2
redg(6,1) = 3 ; redg(6,2) = 4
!--

!--
ii = 0
!--

do i = 1,M

 ie = elems(i)

 xe = mat_coord(mat_elem(ie,:),1);
 ye = mat_coord(mat_elem(ie,:),2);
 ze = mat_coord(mat_elem(ie,:),3);

!-- Calculando as coordenadas dos centros dos elementos.
 xm = sum(xe)/size(xe);
 ym = sum(ye)/size(ye);
 zm = sum(ze)/size(ze);
!--

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

 vol   = dabs(vol_tetraedro(xe,ye,ze));

 !--

 ii = ii + 1
 cxe(ii)=xm; cye(ii)=ym; cze(ii)=zm

 sx = 0.; sy = 0.; sz = 0.;
 do j = 1,nedgs

    i1 = redg(j,1); i2 = redg(j,2);

    dx = mat_coord(mat_elem(ie,i2),1) - mat_coord(mat_elem(ie,i1),1)
    dy = mat_coord(mat_elem(ie,i2),2) - mat_coord(mat_elem(ie,i1),2)
    dz = mat_coord(mat_elem(ie,i2),3) - mat_coord(mat_elem(ie,i1),3)

    si = (mat_elem(ie,i2) - mat_elem(ie,i1))/abs(mat_elem(ie,i2) - mat_elem(ie,i1))

    leng = si * dsqrt( dx**2 + dy**2 + dz**2 )

    Nx   = (leng / (6 * vol)) * ( (fnodal( vol,i1,xm,ym,zm,a,b,c,d ) * b(i2))-(fnodal( vol,i2,xm,ym,zm,a,b,c,d ) * b(i1)) )

    Ny   = (leng / (6 * vol)) * ( (fnodal( vol,i1,xm,ym,zm,a,b,c,d ) * c(i2))-(fnodal( vol,i2,xm,ym,zm,a,b,c,d ) * c(i1)) )

    Nz   = (leng / (6 * vol)) * ( (fnodal( vol,i1,xm,ym,zm,a,b,c,d ) * d(i2))-(fnodal( vol,i2,xm,ym,zm,a,b,c,d ) * d(i1)) )

    sx = sx + Nx * Et( gEdge(ie,j) )
    sy = sy + Ny * Et( gEdge(ie,j) ) 
    sz = sz + Nz * Et( gEdge(ie,j) ) 
    
 end do

!--
 fedD(ii,1) = sx
 fedD(ii,2) = sy
 fedD(ii,3) = sz
!--

end do

!--

allocate( mA(M,Ncs), u(M), y(M), cW(M,M) )
allocate( yy(M),yz(M) )

 cW = 0.d0
ca = 0.d0
ii=0
do i = 1,M
     ie = elems(i)
    !-- 
     xe = mat_coord(mat_elem(ie,:),1);
     ye = mat_coord(mat_elem(ie,:),2);
     ze = mat_coord(mat_elem(ie,:),3);
     xm = sum(xe)/size(xe);
     ym = sum(ye)/size(ye);
     zm = sum(ze)/size(ze);
    !--

     ii=ii+1
     call cf(Ncs,xm,ym,zm,f)

     mA(ii,:)  = f

     y(ii)     = fedD(ii,1)
     yy(ii)    = fedD(ii,2)
     yz(ii)    = fedD(ii,3)

     cW(ii,ii) = dexp( - cc**2 * ( ((xr-xm)/(xr-maxval(cxe)))**2 + ((yr-ym)/(yr-maxval(cye)))**2 &
                 + ((zr-zm)/(zr-maxval(cze)))**2 ) )
      
end do

allocate(Ac(Ncs,Ncs), Bc(Ncs))
allocate(aTw(Ncs,M))

do i = 1,Ncs
   do j = 1,M   
      aTw(i,j) = sum( mA(:,i)*cW(:,j) )
   end do
end do

!-- Ex

do i = 1,Ncs
   do j = 1,Ncs   
      Ac(i,j) = sum( aTw(i,:)*mA(:,j) )
   end do
   Bc(i) = sum( aTw(i,:)*y(:) )
end do
call LU_Gauss_compacto( Ac, Bc, Ncs)
call cf(Ncs,xr,yr,zr,f)

mx = Bc

Ex = sum(f*Bc)

!-- Ey

do i = 1,Ncs
   do j = 1,Ncs   
      Ac(i,j) = sum( aTw(i,:)*mA(:,j) )
   end do
   Bc(i) = sum( aTw(i,:)*yy(:) )
end do
call LU_Gauss_compacto( Ac, Bc, Ncs )
call cf(Ncs,xr,yr,zr,f)

my = Bc

Ey = sum(f*Bc)

!-- Ez

do i = 1,Ncs
   do j = 1,Ncs   
      Ac(i,j) = sum( aTw(i,:)*mA(:,j) )
   end do
   Bc(i) = sum( aTw(i,:)*yz(:) )
end do
call LU_Gauss_compacto( Ac, Bc, Ncs )
call cf(Ncs,xr,yr,zr,f)

mz = Bc

Ez = sum(f*Bc)

!-- Campo Magnetico --

fw = 2.d0 * pi * freq

call dcfdx(Ncs,xr,yr,zr,dfdx)
call dcfdy(Ncs,xr,yr,zr,dfdy)
call dcfdz(Ncs,xr,yr,zr,dfdz)

Hx = -(sum(dfdy*mz) - sum(dfdz*my))/(ip*fw*mi0)
Hy = -(sum(dfdz*mx) - sum(dfdx*mz))/(ip*fw*mi0)
Hz = -(sum(dfdx*my) - sum(dfdy*mx))/(ip*fw*mi0)


!call cf(Ncs,xr,yr,zr+perc,fma)
!call cf(Ncs,xr,yr,zr-perc,fme)
!dExdz = ( sum(fma*mx) - sum(fme*mx) )/(2.d0*perc)

!call cf(Ncs,xr,yr+perc,zr,fma)
!call cf(Ncs,xr,yr-perc,zr,fme)
!dExdy = ( sum(fma*mx) - sum(fme*mx) )/(2.d0*perc)

!call cf(Ncs,xr,yr,zr+perc,fma)
!call cf(Ncs,xr,yr,zr-perc,fme)
!dEydz = ( sum(fma*my) - sum(fme*my) )/(2.d0*perc)

!call cf(Ncs,xr+perc,yr,zr,fma)
!call cf(Ncs,xr-perc,yr,zr,fme)
!dEydx = ( sum(fma*my) - sum(fme*my) )/(2.d0*perc)

!call cf(Ncs,xr,yr+perc,zr,fma)
!call cf(Ncs,xr,yr-perc,zr,fme)
!dEzdy = ( sum(fma*mz) - sum(fme*mz) )/(2.d0*perc)

!call cf(Ncs,xr+perc,yr,zr,fma)
!call cf(Ncs,xr-perc,yr,zr,fme)
!dEzdx = ( sum(fma*mz) - sum(fme*mz) )/(2.d0*perc)

!Hx = -( dEzdy - dEydz )/(ip*fw*mi0)
!Hy = -( dExdz - dEzdx )/(ip*fw*mi0)
!Hz = -( dEydx - dExdy )/(ip*fw*mi0)

!--

!print*,cdabs(Ex),cdabs(Ey),cdabs(Ez)

!-- Campos na direcao do dipolo --

Ex = sum(f*mx) * cosd(aTx(itrm)) + sum(f*my) * sind(aTx(itrm))

end subroutine

!--

subroutine compose_Esec2 ( xr, yr, zr, Ex, Ey, Ez, Hx, Hy, Hz, Et )
! 
! Calcula o campo interpolando o valor apenas no elemento ao qual o n\'o pertence.
! 
implicit none

real(db),intent(in)  :: xr, yr, zr
complex(db),intent(in)  :: Et(:)
complex(db),intent(out) :: Ex,Ey,Ez
complex(db),intent(out) :: Hx, Hy, Hz
integer :: elems(1)
integer :: i,j,k,ce,M,ie

real(db)             :: dx,dy,dz,leng
real(db)             :: vol,Nx,Ny,Nz,xm,ym,zm
real(db),allocatable :: xe(:),ye(:),ze(:),dists(:)
real(db),allocatable :: xa(:),ya(:),za(:) 
real(db),allocatable :: a(:),b(:),c(:),d(:)
integer              :: i1,i2,j1,j2,aux
complex(db) :: sx,sy,sz
integer                 :: redg(nedgs,2), si

real(db) :: Pr, Pce, drce, dr

!--

complex(db),allocatable    :: fedD(:,:)
complex(db)                :: f(10), ca(10) 
complex(db),allocatable    :: mA(:,:), u(:), y(:), cW(:,:), Ac(:,:), Bc(:)
integer                    :: ii, Ncs = 10
real(db)                   :: cc = 1.d0
real(db),allocatable       :: cxe(:), cye(:), cze(:)
integer                    :: ndga

!--

complex(db)                :: dExdz, dExdy
complex(db)                :: dEydz, dEydx
complex(db)                :: dEzdy, dEzdx
real(db)                   :: fw

complex(db)                :: dNxdz, dNxdy
complex(db)                :: dNydz, dNydx
complex(db)                :: dNzdy, dNzdx

complex(db)                :: sdNxdz, sdNxdy
complex(db)                :: sdNydz, sdNydx
complex(db)                :: sdNzdy, sdNzdx


!--


dr = rad

allocate(dists(nelemg))
allocate(xa(3),ya(3),za(3))
allocate(xe(nnoele),ye(nnoele),ze(nnoele))
allocate(a(nnoele),b(nnoele),c(nnoele),d(nnoele))

!-- procura em qual elemento esta o no, pela minima distancia ao centro do mesmo...

do i = 1,nelemg
   xe = mat_coord(mat_elem(i,:),1);
   ye = mat_coord(mat_elem(i,:),2);
   ze = mat_coord(mat_elem(i,:),3);
   xm = sum(xe)/nnoele
   ym = sum(ye)/nnoele
   zm = sum(ze)/nnoele

   drce  = dsqrt( (xr-xm)**2 + (yr-ym)**2 + (zr-zm)**2 )
   dists(i) = drce
end do
!-- Para o caso de um elemento.
M = 1 ! Apenas um elemento, o de distancia minima
elems(M) = minloc(dists,dim=1)
!-- Para mais elementos, dentro do raio estabelecido.

Ex = 0.d0
Ey = 0.d0
Ez = 0.d0

redg(1,1) = 1 ; redg(1,2) = 2
redg(2,1) = 1 ; redg(2,2) = 3
redg(3,1) = 1 ; redg(3,2) = 4
redg(4,1) = 2 ; redg(4,2) = 3
redg(5,1) = 4 ; redg(5,2) = 2 ! 4  2
redg(6,1) = 3 ; redg(6,2) = 4

!--

do i = 1,M

 ie = elems(i)

!--
 xe = mat_coord(mat_elem(ie,:),1);
 ye = mat_coord(mat_elem(ie,:),2);
 ze = mat_coord(mat_elem(ie,:),3);
!--

 xm = xr
 ym = yr
 zm = zr

!--
! print*,xm,ym,zm,'coordenadas'
!--

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

 vol   = dabs(vol_tetraedro(xe,ye,ze));

!--

 sdNxdz= 0. ; sdNxdy= 0.
 sdNydz= 0. ; sdNydx= 0.
 sdNzdy= 0. ; sdNzdx= 0.

 sx = 0.; sy = 0.; sz = 0.;

 do j = 1,nedgs

    i1 = redg(j,1); i2 = redg(j,2);

    dx = mat_coord(mat_elem(ie,i2),1) - mat_coord(mat_elem(ie,i1),1)
    dy = mat_coord(mat_elem(ie,i2),2) - mat_coord(mat_elem(ie,i1),2)
    dz = mat_coord(mat_elem(ie,i2),3) - mat_coord(mat_elem(ie,i1),3)

    si = (mat_elem(ie,i2) - mat_elem(ie,i1))/abs( mat_elem(ie,i2) - mat_elem(ie,i1) )

    leng = si * dsqrt( dx**2 + dy**2 + dz**2 )

    Nx   = (leng / (6 * vol)) * ( (fnodal( vol,i1,xm,ym,zm,a,b,c,d ) * b(i2))-(fnodal( vol,i2,xm,ym,zm,a,b,c,d ) * b(i1)) )

    Ny   = (leng / (6 * vol)) * ( (fnodal( vol,i1,xm,ym,zm,a,b,c,d ) * c(i2))-(fnodal( vol,i2,xm,ym,zm,a,b,c,d ) * c(i1)) )

    Nz   = (leng / (6 * vol)) * ( (fnodal( vol,i1,xm,ym,zm,a,b,c,d ) * d(i2))-(fnodal( vol,i2,xm,ym,zm,a,b,c,d ) * d(i1)) )

    sx = sx + Nx * Et( gEdge(ie,j) )
    sy = sy + Ny * Et( gEdge(ie,j) ) 
    sz = sz + Nz * Et( gEdge(ie,j) ) 

    !--
    dNxdz = (leng / (6 * vol)) * ( (d(i1)*b(i2)/(6*vol)) - (d(i2)*b(i1)/(6*vol)) )
    dNxdy = (leng / (6 * vol)) * ( (c(i1)*b(i2)/(6*vol)) - (c(i2)*b(i1)/(6*vol)) )
    
    dNydz = (leng / (6 * vol)) * ( (d(i1)*c(i2)/(6*vol)) - (d(i2)*c(i1)/(6*vol)) )
    dNydx = (leng / (6 * vol)) * ( (b(i1)*c(i2)/(6*vol)) - (b(i2)*c(i1)/(6*vol)) )

    dNzdy = (leng / (6 * vol)) * ( (c(i1)*d(i2)/(6*vol)) - (c(i2)*d(i1)/(6*vol)) )
    dNzdx = (leng / (6 * vol)) * ( (b(i1)*d(i2)/(6*vol)) - (b(i2)*d(i1)/(6*vol)) )
    !--

    sdNxdz = sdNxdz + dNxdz * Et( gEdge(ie,j) )
    sdNxdy = sdNxdy + dNxdy * Et( gEdge(ie,j) ) 

    sdNydz = sdNydz + dNydz * Et( gEdge(ie,j) ) 
    sdNydx = sdNydx + dNydx * Et( gEdge(ie,j) )

    sdNzdy = sdNzdy + dNzdy * Et( gEdge(ie,j) ) 
    sdNzdx = sdNzdx + dNzdx * Et( gEdge(ie,j) ) 
    
 end do

 Ex = sx 
 Ey = sy
 Ez = sz

 dExdz= sdNxdz ; dExdy= sdNxdy
 dEydz= sdNydz ; dEydx= sdNydx
 dEzdy= sdNzdy ; dEzdx= sdNzdx

!-- Campo na direcao do dipolo --

 Ex = sx * cosd(aTx(itrm)) + sy * sind(aTx(itrm))

!--
end do

!-- Composicao das componentes do campo magnetico --
fw = 2.d0 * pi * freq
Hx = -( dEzdy - dEydz )/(ip*fw*mi0)
Hy = -( dExdz - dEzdx )/(ip*fw*mi0)
Hz = -( dEydx - dExdy )/(ip*fw*mi0)
!--

!print*,cdabs(Ex),cdabs(Ey),cdabs(Ez)

end subroutine

!--

subroutine cf(n,x,y,z,f)
implicit none
integer,intent(in)   :: n
real(db),intent(in)  :: x,y,z
complex(db),intent(out) :: f(n)

f(1) = 1.d0 ; f(2) = x; f(3) = y;
f(4) = z; f(5) = x*y; f(6) = x*z;
f(7) = y*z; f(8) = x**2; f(9) = y**2;
f(10) = z**2;

end subroutine

!--
subroutine dcfdx(n,x,y,z,f)
implicit none
integer,intent(in)   :: n
real(db),intent(in)  :: x,y,z
complex(db),intent(out) :: f(n)

f(1) = 0.d0 ; f(2) = 1.d0; f(3) = 0.d0;
f(4) = 0.d0; f(5) = y; f(6) = z;
f(7) = 0.d0; f(8) = x*2.d0; f(9) = 0.d0;
f(10) = 0.d0;

end subroutine
!--
subroutine dcfdy(n,x,y,z,f)
implicit none
integer,intent(in)   :: n
real(db),intent(in)  :: x,y,z
complex(db),intent(out) :: f(n)

f(1) = 0.d0 ; f(2) = 0.d0; f(3) = 1.d0;
f(4) = 0.d0; f(5) = x; f(6) = 0.d0;
f(7) = z; f(8) = 0.d0; f(9) = y*2.d0;
f(10) = 0.d0;

end subroutine
!--
subroutine dcfdz(n,x,y,z,f)
implicit none
integer,intent(in)   :: n
real(db),intent(in)  :: x,y,z
complex(db),intent(out) :: f(n)

f(1) = 0.d0 ; f(2) = 0.d0; f(3) = 0.d0;
f(4) = 1.d0; f(5) = 0.d0; f(6) = x;
f(7) = y; f(8) = 0.d0; f(9) = 0.d0;
f(10) = z*2.d0;

end subroutine
!--

real(db) function fnodal( vol,i,x,y,z,a,b,c,d )
implicit none
integer,intent(in)  :: i
real(db),intent(in) :: x,y,z
real(db),intent(in) :: vol,a(nnoele),b(nnoele),c(nnoele),d(nnoele)

fnodal = (1.d0/(6*vol)) * ( a(i) + b(i)*x + c(i)*y + d(i)*z )

end function

!--

real(db) function arrayA( ie,i,j,b,c,d,i1,i2,j1,j2,vol )
implicit none
integer,intent(in)  :: ie,i,j                            ! Indice do elemento e Indices das atestas 1-6.
real(db),intent(in) :: vol,b(nnoele),c(nnoele),d(nnoele) ! Coeficientes das funcoes base, volume.
integer,intent(inout)  :: i1,i2,j1,j2

integer :: aux
real(db) :: leni , lenj, delx, dely, delz
integer  :: Si,Sj

delx = mat_coord(mat_elem(ie,i2),1) - mat_coord(mat_elem(ie,i1),1)
dely = mat_coord(mat_elem(ie,i2),2) - mat_coord(mat_elem(ie,i1),2)
delz = mat_coord(mat_elem(ie,i2),3) - mat_coord(mat_elem(ie,i1),3)

si   =  (mat_elem(ie,i2)-mat_elem(ie,i1))/abs(mat_elem(ie,i2)-mat_elem(ie,i1))

leni = si * dsqrt( delx**2 + dely**2 + delz**2 )

delx = mat_coord(mat_elem(ie,j2),1) - mat_coord(mat_elem(ie,j1),1)
dely = mat_coord(mat_elem(ie,j2),2) - mat_coord(mat_elem(ie,j1),2)
delz = mat_coord(mat_elem(ie,j2),3) - mat_coord(mat_elem(ie,j1),3)

sj   =  (mat_elem(ie,j2)-mat_elem(ie,j1))/abs(mat_elem(ie,j2)-mat_elem(ie,j1))

lenj = sj * dsqrt( delx**2 + dely**2 + delz**2 )

arrayA = ( 4.d0 * (leni * lenj * vol) / ( 6.d0 * vol )**4 ) * &
         ( (( c(i1)*d(i2) - c(i2)*d(i1) ) * ( c(j1)*d(j2) - c(j2)*d(j1) ))   + &
           (( d(i1)*b(i2) - d(i2)*b(i1) ) * ( d(j1)*b(j2) - d(j2)*b(j1) ))   + &
           (( b(i1)*c(i2) - b(i2)*c(i1) ) * ( b(j1)*c(j2) - b(j2)*c(j1) )) ) 

end function

!--

real(db) function arrayB( ie,i,j,b,c,d,i1,i2,j1,j2,vol )
implicit none
integer,intent(in)  :: ie,i,j ! Indice do elemento e Indices das atestas 1-6.
real(db),intent(in) :: vol,b(nnoele),c(nnoele),d(nnoele) ! Coeficientes das funcoes base, volume.
integer,intent(inout)  :: i1,i2,j1,j2

integer :: aux
real(db) :: leni , lenj, delx, dely, delz
integer  :: Si, Sj

delx = mat_coord(mat_elem(ie,i2),1) - mat_coord(mat_elem(ie,i1),1)
dely = mat_coord(mat_elem(ie,i2),2) - mat_coord(mat_elem(ie,i1),2)
delz = mat_coord(mat_elem(ie,i2),3) - mat_coord(mat_elem(ie,i1),3)

si   =  (mat_elem(ie,i2)-mat_elem(ie,i1))/abs(mat_elem(ie,i2)-mat_elem(ie,i1))

leni = si * dsqrt( delx**2 + dely**2 + delz**2 )

delx = mat_coord(mat_elem(ie,j2),1) - mat_coord(mat_elem(ie,j1),1)
dely = mat_coord(mat_elem(ie,j2),2) - mat_coord(mat_elem(ie,j1),2)
delz = mat_coord(mat_elem(ie,j2),3) - mat_coord(mat_elem(ie,j1),3)

sj   =  (mat_elem(ie,j2)-mat_elem(ie,j1))/abs(mat_elem(ie,j2)-mat_elem(ie,j1))

lenj = sj * dsqrt( delx**2 + dely**2 + delz**2 )

arrayB = ( (leni*lenj) / (6.d0*vol)**2 ) * (vol/20.d0) * ( ( (1 + delta(i1,j1)) * fipjq( i2,j2,b,c,d ) ) - &
           ( (1 + delta(i1,j2)) * fipjq( i2,j1,b,c,d ) ) - &
           ( (1 + delta(i2,j1)) * fipjq( i1,j2,b,c,d ) ) + &
           ( (1 + delta(i2,j2)) * fipjq( i1,j1,b,c,d ) ) ) 

end function

!--

real(db) function fipjq( ip,jq,b,c,d )
implicit none
integer,intent(in)  :: ip,jq 
real(db),intent(in) :: b(nnoele),c(nnoele),d(nnoele) ! Coeficientes das funcoes base, volume.

fipjq = b(ip)*b(jq) + c(ip)*c(jq) + d(ip)*d(jq)

end function

!--

subroutine Montagem_Matriz_Global_tetraedro( ngl , nnos , nel , n_nz , vet_nn )
implicit none
integer,intent(in) :: nel, n_nz , ngl , nnos
complex(db),intent(out) :: vet_nn(n_nz)

integer :: ie,i,j,k,l,m
real(db),allocatable :: xe(:),ye(:),ze(:)
complex(db),allocatable :: mloc(:,:)
real(db),allocatable :: a(:),b(:),c(:),d(:)
real(db),allocatable :: xa(:),ya(:),za(:) 
integer,allocatable  :: lnode(:)

real(db)                :: vol , rho_e , sig_e , omega , zm
complex(db),allocatable :: M_base(:,:)
complex(db)             :: kk

complex(db),allocatable :: Al(:,:),Bl(:,:)

integer                 :: i1,i2,j1,j2

integer                 :: ci, cj

integer                 :: redg(nedgs,2)

allocate(lnode(nedgs))

allocate(xa(3),ya(3),za(3))
allocate(xe(nnoele),ye(nnoele),ze(nnoele))
allocate(a(nnoele),b(nnoele),c(nnoele),d(nnoele))

allocate(Al(nedgs,nedgs), & 
         Bl(nedgs,nedgs))

allocate(mloc(nedgs,nedgs))

vet_nn = (0.d0,0.d0)

omega = 2*pi*Freq ;

redg(1,1) = 1 ; redg(1,2) = 2
redg(2,1) = 1 ; redg(2,2) = 3
redg(3,1) = 1 ; redg(3,2) = 4
redg(4,1) = 2 ; redg(4,2) = 3
redg(5,1) = 4 ; redg(5,2) = 2
redg(6,1) = 3 ; redg(6,2) = 4

do ie = 1,nel

 xe = mat_coord(mat_elem(ie,:),1);
 ye = mat_coord(mat_elem(ie,:),2);
 ze = mat_coord(mat_elem(ie,:),3);

!--

 do i = 1,6
! do i = 1,4
    do m = ngl-1,0,-1
       lnode( ngl*i - m ) = ngl * gEdge(ie,i) - m
!       lnode( ngl*i - m ) = ngl * mat_elem(ie,i) - m
    end do
 end do

!--

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

kk = (ip*omega*mi0)*(sig_e + ip*omega*ep0)

Al=0.d0;Bl=0.d0;

mloc = 0.d0

do i = 1,nedgs
   i1 = redg(i,1) ; i2 = redg(i,2)
   do j = 1,nedgs
      j1 = redg(j,1) ; j2 = redg(j,2)
      Al(i,j)   = arrayA( ie,i,j,b,c,d,i1,i2,j1,j2,vol )
      Bl(i,j)   = arrayB( ie,i,j,b,c,d,i1,i2,j1,j2,vol )
      mloc(i,j) = Al(i,j) + kk * Bl(i,j)         
   end do
end do

!--
call set_arrays_pds(cdim,ngl,nnosg,nedgst,nnz,row_ind,lnode,mloc,col_ptr,val_nz)
!--

end do ! loop elemosntos

deallocate(lnode)
deallocate(mloc)
deallocate(xa,ya,za)
deallocate(xe,ye,ze)
deallocate(a,b,c,d)

deallocate(Al, & 
           Bl)

end subroutine

!--

subroutine Montagem_Vetor_Fonte_tetraedro( ngl , nel , nnode , vetFxyz )
implicit none
integer,intent(in) :: ngl, nnode, nel
complex(db),intent(out) :: vetFxyz(:)

integer :: ie,i,j,k,l,m
real(db),allocatable    :: xe(:),ye(:),ze(:)
complex(db),allocatable :: flocx(:),flocy(:),flocz(:),floc(:)
real(db),allocatable    :: a(:),b(:),c(:),d(:)
real(db),allocatable    :: xa(:),ya(:),za(:) 
integer,allocatable     :: lnode(:)

real(db)                :: vol , rho_e , sig_e , sig_p , del_sig , omega , zm
complex(db),allocatable :: Ep_x(:) , Ep_y(:) , Ep_z(:)  
integer                 :: cmel  

integer                 :: k1
real(db)                :: rEx,rEy,rEz,iEx,iEy,iEz

integer                 :: cont

integer                 :: i1,i2,j1,j2,aux

real(db) :: xm, ym
real(db) :: dx, dy, dz, leng

complex(db),allocatable :: Ept(:), soma

integer                 :: redg(nedgs,2)

integer                 :: Si

allocate(lnode(ngl*nedgs))
allocate(Ep_x(ngl*nedgs),Ep_y(ngl*nedgs),Ep_z(ngl*nedgs))
allocate(Ept(ngl*nedgs))

allocate(xa(3),ya(3),za(3))
allocate(xe(nnoele),ye(nnoele),ze(nnoele))
allocate(a(nnoele),b(nnoele),c(nnoele),d(nnoele))

allocate(floc(ngl*nedgs)) 

omega = 2*pi*Freq

redg(1,1) = 1 ; redg(1,2) = 2
redg(2,1) = 1 ; redg(2,2) = 3
redg(3,1) = 1 ; redg(3,2) = 4
redg(4,1) = 2 ; redg(4,2) = 3
redg(5,1) = 4 ; redg(5,2) = 2
redg(6,1) = 3 ; redg(6,2) = 4

do ie = 1,nel

 xe = mat_coord(mat_elem(ie,:),1);
 ye = mat_coord(mat_elem(ie,:),2);
 ze = mat_coord(mat_elem(ie,:),3);

!--
 do i = 1,6
! do i = 1,4
    do m = ngl-1,0,-1
       lnode( ngl*i - m ) = ngl * gEdge(ie,i) - m
!       lnode( ngl*i - m ) = ngl * mat_elem(ie,i) - m
    end do
 end do

!--

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

   do k = 1,nedgs

        i1 = redg(k,1) ; i2 = redg(k,2)

        call halfedge ( xm,ym,zm,dx,dy,dz,leng,mat_elem(ie,i1),mat_elem(ie,i2) )

      if(uprm==0)then

        ttap=0.d0
        call prim_MCSEM( ttap, Freq , tr_x , tr_y , tr_z , &
                         xm,ym,zm, &
                         Ep_x(k),Ep_y(k),Ep_z(k) )

      elseif(uprm==1)then

        k1 = (ifrq-1)*sum(NTPP)*NInh + (itrm-1)*NInh + Inh( gEdge(ie,k) )
		 			 
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

        si   =  (mat_elem(ie,i2)-mat_elem(ie,i1))/abs(mat_elem(ie,i2)-mat_elem(ie,i1))

        leng = si * leng

        Ept(k) = ( (Ep_x(k)*dx) + (Ep_y(k)*dy) + (Ep_z(k)*dz) ) / (leng)

   end do

   floc = 0.d0;

   do i = 1,nedgs
       i1 = redg(i,1) ; i2 = redg(i,2)
       soma = 0.d0    
       do j = 1,nedgs 
          j1 = redg(j,1) ; j2 = redg(j,2)
          soma = soma + arrayB( ie,i,j,b,c,d,i1,i2,j1,j2,vol ) * (Ept(j))
       end do
       floc(i) = -(ip*omega*mi0) * del_sig * soma
   end do

else

  floc  = (0.d0,0.d0)

end if

!%%% chamar aquia a montagem do vetor fonte global %%%%

do l = 1,nedgs
      vetFxyz( lnode(l) ) = vetFxyz( lnode(l) ) + floc( l )
end do

!================================================

end do ! loop elemosntos

deallocate(lnode)
deallocate(floc)
deallocate(Ep_x,Ep_y,Ep_z)
deallocate(Ept)
deallocate(xa,ya,za)
deallocate(xe,ye,ze)
deallocate(a,b,c,d)

end subroutine

!--

subroutine halfedge ( xm,ym,zm,delx,dely,delz,leng,no1,no2 )
implicit none
integer,intent(in)   :: no1,no2
real(db),intent(out) :: xm,ym,zm,delx,dely,delz,leng

xm = 0.5d0*(mat_coord(no1,1) + mat_coord(no2,1))
ym = 0.5d0*(mat_coord(no1,2) + mat_coord(no2,2))
zm = 0.5d0*(mat_coord(no1,3) + mat_coord(no2,3))

delx = mat_coord(no2,1) - mat_coord(no1,1)
dely = mat_coord(no2,2) - mat_coord(no1,2)
delz = mat_coord(no2,3) - mat_coord(no1,3)

leng = dsqrt( delx**2 + dely**2 + delz**2 )

end subroutine

!--

real(db) function vol_tetraedro ( x,y,z )
implicit none
real(db),intent(in) :: x(4),y(4),z(4)

vol_tetraedro = (1.d0/6.d0) * ( xij(x(2),x(1))*( yij(y(2),y(3))*zij(z(3),z(4)) - yij(y(3),y(4))*zij(z(2),z(3)) ) + &
                            xij(x(3),x(2))*( yij(y(3),y(4))*zij(z(1),z(2)) - yij(y(1),y(2))*zij(z(3),z(4)) ) + & 
                            xij(x(4),x(3))*( yij(y(1),y(2))*zij(z(2),z(3)) - yij(y(2),y(3))*zij(z(1),z(2)) ) )

end function

!--

real(db) function det ( x , y , z )
implicit none
real(db),intent(in) :: x(3),y(3),z(3)

det = x(1)*( y(2)*z(3) - y(3)*z(2) ) + y(1)*( z(2)*x(3) - z(3)*x(2) ) + z(1)*( x(2)*y(3) - x(3)*y(2) )

end function

!--

real(db) function xij(xi,xj)
implicit none
real(db),intent(in) :: xi,xj
  xij = xi - xj
end function

!--

real(db) function yij(yi,yj)
implicit none
real(db),intent(in) :: yi,yj
  yij = yi - yj
end function

!--

real(db) function zij(zi,zj)
implicit none
real(db),intent(in) :: zi,zj
  zij = zi - zj
end function

!--

integer function delta(i,j)
implicit none
integer,intent(in) :: i,j

if(i.eq.j)then
    delta = 1;
else
    delta = 0;
end if

end function


!--

end module
