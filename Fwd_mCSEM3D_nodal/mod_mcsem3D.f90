!include '/opt/newintel/intel/compilers_and_libraries_2019.1.144/linux/mkl/include/mkl_pardiso.f90'  ! (path da sua maquina) Se usar o compilador gfortran
include 'mkl_pardiso.f90'  ! (path da sua maquina) Se usar o compilador gfortran
module campo_MCSEM3D

use var_glob
use Arrays_SLU
use deriv_MQMP4
use Fem_3D	
use DEHx
use mkl_pardiso
use modulo_auxiliar
use arrayspds

contains

!%%%===========================================%%%

subroutine Foward_MCSEM3D(path,np,p,no,Ob,n_p,my_rank )
implicit none

include 'mpif.h'

 integer,intent(in)   :: np,no,n_p,my_rank
 character(len=100),intent(in) :: path
 real(dp),intent(in)  :: p(np)
 real(dp),intent(out) :: Ob(no)

real(dp)             :: ti,tf
integer,allocatable  :: nrc_ptsm(:,:),NDPF(:)
integer :: i,j,nc
complex(dp),allocatable :: Ex(:),Ey(:),Ez(:)
complex(dp),allocatable :: Hx(:),Hy(:),Hz(:)

Ob = 0.d0;

!----------------------------------
call inicia_parametros (  )
!----------------------------------

allocate(NDPF(nfreqs)); call Estima_num_dados( nc, NDPF ) ; deallocate(NDPF)

call primarios_modelo(dprm,np,nem1,ihet1,Inh,NInh,n_p,my_rank) ; dprm = 0;

allocate(Ex(nc),Ey(nc),Ez(nc))
allocate(Hx(nc),Hy(nc),Hz(nc))

call open_arqprim (  )

  call Campos_MCSEM3D_par( np, nc, P, Ex, Ey, Ez, Hx, Hy, Hz, n_p, my_rank )
  call observacoes ( nc,Ex,no,Ob )

call close_arqprim
deallocate(Ex,Ey,Ez)
deallocate(Hx,Hy,Hz)

call desaloc_vglob
!----------------------------------
call finaliza_parametros
!----------------------------------

end subroutine

!%%%===========================================%%%

subroutine observacoes ( nc,E,no,Obs )
implicit none
integer,intent(in)     :: nc,no
complex(8),intent(in)  :: E(nc)
real(8),intent(out)    :: Obs(no)

integer  :: i,j,nr,cont
real(8),allocatable :: phi(:),phase(:)

cont = 0
do i = 1,nc
  cont = cont + 1
  Obs(cont)=dlog10(cdabs(E(i)))
end do

end subroutine

!%%%===========================================%%%

subroutine unwrap(xw,n)
implicit none
integer,intent(in)  :: n
real(dp),intent(inout) :: xw(n)

integer :: i
real(dp) :: difference
real(dp),allocatable :: xu(:)

allocate(xu(n))
xu = xw;
do i=2,n
	difference = xw(i) - xw(i-1);
	if (difference > pi)then
		xu(i:n) = xu(i:n) - 2*pi;
	elseif (difference < -pi)then
		xu(i:n) = xu(i:n) + 2*pi;     
	end if
end do
xw = xu;
deallocate(xu)

end subroutine

!%%%===========================================%%%

subroutine para_range2(n1, n2, nprocs, irank, ista, iend)
implicit none

integer,intent(in)  :: n1,n2,nprocs,irank
integer,intent(out) :: ista,iend
integer :: iwork1,iwork2

iwork1 = (n2 - n1 + 1) / nprocs
iwork2 = mod(n2 - n1 + 1, nprocs)
ista = irank * iwork1 + n1 + min(irank, iwork2)
iend = ista + iwork1 - 1
if (iwork2 > irank) iend = iend + 1

end subroutine para_range2

!%%%===========================================%%%

subroutine Campos_MCSEM3D_par( npr , ncp , P , cEx , cEy , cEz, cHx , cHy , cHz , n_p , my_rank )
implicit none

include 'mpif.h'

integer,intent(in)      :: n_p , my_rank
integer,intent(in)      :: npr , ncp
real(dpc),intent(in)    :: P(npr)
complex(dp),intent(out) :: cEx( ncp ) , cEy( ncp ), cEz( ncp )
complex(dp),intent(out) :: cHx( ncp ) , cHy( ncp ), cHz( ncp )

integer :: i,j,k,l, cont 
complex(dpc),allocatable :: sol(:)
complex(dpc) :: E_aux , zeta0 , Expm , Eypm , Ezpm

integer :: id_deriv = 8    

real(dpc)               :: delx,dely,modu,ti,tf, tii, tff

!%%% variaveis pardiso %%%%

INTEGER ::maxfct, mnum, mtype, phase, nnz, error, msglvl, ng, nrhs
INTEGER, ALLOCATABLE :: iparm( : ),perm(:)
INTEGER :: idum(1) , iopt
complex(dpc):: ddum(1)
TYPE(MKL_PARDISO_HANDLE), ALLOCATABLE  :: pt(:)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%

real(dpc),allocatable :: redAxdx(:), redAxdy(:), redAxdz(:)
real(dpc),allocatable :: imdAxdx(:), imdAxdy(:), imdAxdz(:)
real(dpc),allocatable :: redAydx(:), redAydy(:), redAydz(:)
real(dpc),allocatable :: imdAydx(:), imdAydy(:), imdAydz(:)
real(dpc),allocatable :: redAzdx(:), redAzdy(:), redAzdz(:)
real(dpc),allocatable :: imdAzdx(:), imdAzdy(:), imdAzdz(:)
real(dpc),allocatable :: redphidx(:), redphidy(:), redphidz(:)
real(dpc),allocatable :: imdphidx(:), imdphidy(:), imdphidz(:)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%

!$$$ Variaveis da Paralelizacao $$$$$$$$$$$$$$$$$$$$$$$
integer :: fri,frf,ierr
integer,allocatable     :: jlec(:),idsc(:)
complex(dpc),allocatable :: Vjc(:) , Vjcy(:) , Vjcz(:)
complex(dpc),allocatable :: hVjc(:) , hVjcy(:) , hVjcz(:)
integer :: icp,it
integer :: contt,contr
integer :: aux
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ 

complex(dpc),allocatable :: VjJac(:)
integer,allocatable     :: jlecJac(:), idscJac(:)
integer                 :: contj

complex(dpc),allocatable :: vjj(:)

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

if(my_rank.eq.0)then
  print*,'-----------------------------'
  print*,'Iniciando o programa MCSEM3D'
  print*,'-----------------------------'
end if

call Atualiza_Parametros(npr,P)
call Atualiza_Prop_malha

!%%%%%%%%

allocate(NDPF(nfreqs)); call Estima_num_dados( aux, NDPF );

Allocate( jlec(0:n_p-1),idsc(0:n_p-1) )
Do i = 0,n_p-1
    call para_range2(1,nfreqs,n_p,i,fri,frf)
    jlec(i) = sum(NDPF(fri:frf))
    if(i.eq.0)then
	    idsc(i) = i
    else 
	    idsc(i) = sum(jlec(0:i-1))
    end if
End do
call para_range2(1,nfreqs,n_p,my_rank,fri,frf) 

allocate( vjc( sum(NDPF(fri:frf)) ) ) ; vjc = 0.d0;
allocate( vjcy( sum(NDPF(fri:frf)) ) ) ; vjcy = 0.d0;
allocate( vjcz( sum(NDPF(fri:frf)) ) ) ; vjcz = 0.d0;

allocate( hvjc( sum(NDPF(fri:frf)) ) ) ; hvjc = 0.d0;
allocate( hvjcy( sum(NDPF(fri:frf)) ) ) ; hvjcy = 0.d0;
allocate( hvjcz( sum(NDPF(fri:frf)) ) ) ; hvjcz = 0.d0;

deallocate(NDPF)

!+++++ Vetores iniciais ++++++
ti = MPI_WTIME();

call initial_arrays(cdim,nelemg,nnosg,ngl,nnzi,mat_elem,col_ptra,row_inda)

tf = MPI_WTIME(); print*,'Tempo arrays inic',(tf-ti)/60.d0,'min'
!+++++++++++++++++++++++++++++

!%%%%%%%%

 irs   = 0;
 cont  = 0; 

 cEx   = (0.d0,0.d0);
 cEy   = (0.d0,0.d0);
 cEz   = (0.d0,0.d0);
 cHx   = (0.d0,0.d0);
 cHy   = (0.d0,0.d0);
 cHz   = (0.d0,0.d0);

 contj = 0; 

do i = fri,frf;
     
     ifrq = i; ! Variavel global

     Freq = Freqs(i);

     zeta0 = ip * 2.d0*pi*Freq * mi0 ; zet0 = zeta0

     nrs = sum( NTPP );

     allocate( vfonte(nrs*ngl*nnosg) ) ; vfonte = (0.d0,0.d0);

     ti = MPI_WTIME();

     nnz=nnzi
     allocate( val_nz(nnz),row_ind(nnz),col_ptr(ngl*nnosg+1) )    
     row_ind=row_inda ; col_ptr=col_ptra;
     val_nz=0.d0;
     call Montagem_Matriz_Global_tetraedro( ngl , nnosg , nelemg , nnz , val_nz )
     call recount_arrays(ngl,nnosg,row_ind,col_ptr,val_nz,nnz)

     tf = MPI_WTIME();
     if(my_rank.eq.0)print*,'Matriz global montada.',(tf-ti)/60.d0,'min'

     ngb = ngl*nnosg;
     call Cond_Fronteira_DH( vnobord , val_nz , row_ind , col_ptr , vfonte , nbord , ngl )

     !-- Pardiso --

     ng = ngl * nnosg;
     nrhs = 1
     maxfct = 1
     mnum = 1
     ALLOCATE( iparm( 64 ) , perm(nnosg*ngl) , sol(nrs*nnosg*ngl) )
     sol=(0.d0,0.d0);
     iparm = 0
     perm = 0
     msglvl = 1
     ALLOCATE ( pt ( 64 ) )
     mtype = 6

     if(my_rank.eq.0)print*,'Memoria Antes:'
     call inf_memory (my_rank)

     ti = MPI_WTIME()

     call pardisoinit(pt,mtype,iparm)

     phase = 11 ! only reordering and symbolic factorization
     CALL pardiso (pt, maxfct, mnum, mtype, phase, ng, val_nz, col_ptr, row_ind, &
	              perm, nrhs, iparm, msglvl, vfonte, sol, error)
     if(error.ne.0) then
          print*,'Erro na fat. simb.';stop
     end if
     if(my_rank.eq.0)print*,'Fatorando...'
     !.. Factorization.
     phase = 22 ! only factorization
     CALL pardiso (pt, maxfct, mnum, mtype, phase,ng, val_nz, col_ptr, row_ind, &
	              perm, nrhs, iparm, msglvl, vfonte, sol, error)
     if(error.ne.0) then
          print*,'Erro na fat.';stop
     end if
     tf = MPI_WTIME();
     print*,''
     if(my_rank==0)print*,'$> Tempo de fatoracao',(tf-ti)/60.d0,'min'
     print*,''

     if(my_rank.eq.0)print*,'Memoria Depois:'
     call inf_memory (my_rank)

     !#####################################

     irs = 0; 
     vfonte = ( 0.d0 , 0.d0 );
     contt = 0;
     do l = 1,Num_perf
        iperf = l; ! Variavel global
       !=========
	ttap = 0.d0 !;
       !=========
       do j = 1,NTPP(l); 	  

          contt = contt + 1; irs = irs + 1;

          itrm = contt

          tr_x = COOTXP(contt,1) ; tr_y = COOTXP(contt,2) ; tr_z = COOTXP(contt,3)           

          call Montagem_Vetor_Fonte_tetraedro( ngl , nelemg , nnosg , vfonte )
          call Cond_Fronteira_DH( vnobord , val_nz , row_ind , col_ptr , vfonte , nbord , ngl )

       end do ! Transmissor
     end do ! Perfis

     nrhs = nrs
     phase = 33 ! system solution
     CALL pardiso (pt, maxfct, mnum, mtype, phase,ng, val_nz, col_ptr, row_ind, &
                        perm, nrhs, iparm, msglvl, vfonte, sol, error)
     if(error.ne.0) then
          print*,'Erro na solucao so sistema' ; stop
     end if

     !#####################################

     irs = 0;
     contt = 0;

     do l = 1,Num_perf
        iperf = l; ! Variavel global
       !=========
	ttap = 0.d0
       !=========

       do j = 1,NTPP(l);          
 
          contt = contt + 1; irs = irs + 1;

          itrm = contt

          tr_x = COOTXP(contt,1) ; tr_y = COOTXP(contt,2) ; tr_z = COOTXP(contt,3)           

	 !===================
          call Numrec_portransm2( l , i , j , tr_x , tr_y , nrcp , ind_perf , ind_perfpxma , ind_perfpxme ); 

          allocate(coor_perf(nrcp,3));
	  coor_perf(:,1) = COORCP( ind_perf(:) , 1 )
          coor_perf(:,2) = COORCP( ind_perf(:) , 2 )
          coor_perf(:,3) = COORCP( ind_perf(:) , 3 )
	 !===================


	  allocate(Ext_in(nrcp) , Exp_in(nrcp))
	  allocate(Eyt_in(nrcp) , Ezt_in(nrcp))
          Ext_in = (0.d0,0.d0) ; Exp_in = (0.d0,0.d0); 
          Eyt_in = (0.d0,0.d0) ; Ezt_in = (0.d0,0.d0); 

	  allocate(Hxt_in(nrcp) , Hxp_in(nrcp))
	  allocate(Hyt_in(nrcp) , Hzt_in(nrcp))
          Hxt_in = (0.d0,0.d0) ; Hxp_in = (0.d0,0.d0); 
          Hyt_in = (0.d0,0.d0) ; Hzt_in = (0.d0,0.d0); 

          if(interp==1)then
            call deriv_fbsEH( Sol, nrcp, coor_perf, Ext_in, Eyt_in, Ezt_in, Hxt_in, Hyt_in, Hzt_in )
          elseif(interp==2)then
            call compose_E(Sol, nrcp, coor_perf, Ext_in, Eyt_in, Ezt_in, Hxt_in, Hyt_in, Hzt_in )
          else
            print*,'Escolha errada do interpolador';stop
          end if

          do k = 1,nrcp
 

                !call prim_MCSEM( ttap , Freq , tr_x , tr_y , tr_z , & 
                !                  coor_perf(k,1),coor_perf(k,2),coor_perf(k,3), & 
                !                  Expm , Eypm, Ezpm )

                Ext_in(k) = Ext_in(k) !+ Expm
                Eyt_in(k) = Eyt_in(k) !+ Eypm
                Ezt_in(k) = Ezt_in(k) !+ Ezpm
                 
                Hxt_in(k) = Hxt_in(k) !+ Hpx
                Hyt_in(k) = Hyt_in(k) !+ Hpy
                Hzt_in(k) = Hzt_in(k) !+ Hpz

          end do

          !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

          !$$$$$$$$$$$$$$$$$$$$ Para a paralelizacao $$$$$$$$
	   do k = 1,nrcp
             cont = cont + 1
	     vjc(cont)  = Ext_in(k)
	     vjcy(cont) = Eyt_in(k)
	     vjcz(cont) = Ezt_in(k)

	     hvjc(cont)  = Hxt_in(k)
	     hvjcy(cont) = Hyt_in(k)
	     hvjcz(cont) = Hzt_in(k)
	   end do
          !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

	  deallocate(Ext_in, Exp_in)
	  deallocate(Eyt_in, Ezt_in)

	  deallocate(Hxt_in, Hxp_in)
	  deallocate(Hyt_in, Hzt_in)

          deallocate( ind_perf , ind_perfpxma , ind_perfpxme , coor_perf )

       end do ! Transmissor

     end do ! Perfis

     nrhs  = 1
     iopt  = 3
     phase = -1 ! release internal memory
     CALL pardiso (pt, maxfct, mnum, mtype, phase,ng, val_nz, col_ptr, row_ind, &
 		   perm, nrhs, iparm, msglvl, vfonte, sol, error)

     deallocate(col_ptr)
     deallocate(val_nz)
     deallocate(row_ind)
     deallocate(vfonte) 
     deallocate(iparm,perm,sol,pt)

     print'(1x,a,1x,I3,1x,a)','Campos calculado para a frequencia',i,'em todos os perfis e fontes.'

end do ! Frequencia

if(n_p.gt.1)then

    call MPI_BARRIER(MPI_COMM_WORLD,ierr)

    call MPI_ALLGATHERV( vjc , jlec(my_rank) , MPI_DOUBLE_COMPLEX , &
                     cEx , jlec , idsc , MPI_DOUBLE_COMPLEX , MPI_COMM_WORLD , ierr )

    call MPI_ALLGATHERV( vjcy , jlec(my_rank) , MPI_DOUBLE_COMPLEX , &
                     cEy , jlec , idsc , MPI_DOUBLE_COMPLEX , MPI_COMM_WORLD , ierr )

    call MPI_ALLGATHERV( vjcz , jlec(my_rank) , MPI_DOUBLE_COMPLEX , &
                     cEz , jlec , idsc , MPI_DOUBLE_COMPLEX , MPI_COMM_WORLD , ierr )


    call MPI_ALLGATHERV( hvjc , jlec(my_rank) , MPI_DOUBLE_COMPLEX , &
                     cHx , jlec , idsc , MPI_DOUBLE_COMPLEX , MPI_COMM_WORLD , ierr )

    call MPI_ALLGATHERV( hvjcy , jlec(my_rank) , MPI_DOUBLE_COMPLEX , &
                     cHy , jlec , idsc , MPI_DOUBLE_COMPLEX , MPI_COMM_WORLD , ierr )

    call MPI_ALLGATHERV( hvjcz , jlec(my_rank) , MPI_DOUBLE_COMPLEX , &
                     cHz , jlec , idsc , MPI_DOUBLE_COMPLEX , MPI_COMM_WORLD , ierr )

elseif(n_p.eq.1)then

	do i = 1,ncp
		cEx(i) = vjc(i)
		cEy(i) = vjcy(i)
		cEz(i) = vjcz(i)

		cHx(i) = hvjc(i)
		cHy(i) = hvjcy(i)
		cHz(i) = hvjcz(i)
	end do

end if

if(my_rank.eq.0)then
  print*,'-------------------------'
  print*,'Fim do programa MCSEM3D'
  print*,'-------------------------'
end if

!+++++++
deallocate(col_ptra,row_inda)
!++++++

deallocate(jlec,idsc)

deallocate(vjc)
deallocate(vjcy,vjcz)
deallocate(hvjc)
deallocate(hvjcy,hvjcz)

end subroutine

!%%%===========================================%%%

subroutine inicia_parametros (  )
implicit none

call leitura_aloc_perfiz()
call leitura_prop_modelo()
call leit_aloc_propmod_var_1D
call leit_aloc_varglob_FEM
call flagsdt

end subroutine

!%%%===========================================%%%

subroutine finaliza_parametros

call desaloca_var_perfiz
call desaloca_var_1D_var_modelo
call desaloca_var_FEM
deallocate(flagF)

end subroutine

!%%%===========================================%%%

subroutine leitura_aloc_perfiz()
implicit none

integer :: i,j,k,cont
real(dpc) :: Aux 
 
open(unit = 50, file ='./parametros_modelo/coord_Rx_perfis.in', status = 'old', action = 'read')
open(unit = 60, file ='./parametros_modelo/coord_Tx_perfis.in', status = 'old', action = 'read')

read(50,*)Num_perf

allocate( NRPP(Num_perf) ) 

do i = 1,Num_perf
          read(50,*)NRPP(i)
          do j = 1,NRPP(i)
		read(50,*)Aux
          end do  
end do 

NRCT = sum( NRPP(:) )

rewind(50)

allocate( COORCP(NRCT,3) , INDRCP(NRCT) , INDRCPX(NRCT,2) ) 

read(50,*)Num_perf
 cont = 0
do i = 1,Num_perf
          read(50,*)NRPP(i)
          do j = 1,NRPP(i)
                cont = cont + 1		
                INDRCP(cont) = cont
		read(50,*)COORCP(cont,1),COORCP(cont,2),COORCP(cont,3)
          end do  
end do

 close(50)

read(60,*)Num_perf

allocate( NTPP(Num_perf) )

do i = 1,Num_perf
       read(60,*)NTPP(i)
       do j = 1,NTPP(i)
          read(60,*)Aux
       end do
end do

NFTT = sum( NTPP(:) )  

allocate( COOTXP(NFTT,3) )

rewind(60)

read(60,*)Num_perf

 cont = 0
do i = 1,Num_perf
       read(60,*)NTPP(i)
       do j = 1,NTPP(i)
          cont = cont + 1
          read(60,*)COOTXP(cont,1),COOTXP(cont,2),COOTXP(cont,3)
       end do
end do

 close(60)

open(unit = 60, file ='./parametros_modelo/frequencias_modelo1d.in', status = 'old', action = 'read')

read(60,*)nFreqs ; allocate(freqs(nfreqs))
read(60,*)freqs

 close(60)

end subroutine

!%%%===========================================%%%

subroutine leitura_prop_modelo( )
implicit none

 integer :: i
 real(dpc) :: Aux

 open(unit = 101, file ='./parametros_modelo/heterogeneidades.in', status = 'old', action = 'read' )

 read(101,*)nhet

 allocate (rhohet(nhet))

 read(101,*)rhohet(1:nhet)
 read(101,*)rad,interp

 close(101)

end subroutine

!%%%===========================================%%%

subroutine leit_aloc_propmod_var_1D
implicit none

integer :: i,j
real(dpc) :: Aux

open( unit = 101, file ='./parametros_modelo/frequencias_modelo1d.in', status = 'old', action = 'read' )

do i = 1,2
	read(101,*)aux
end do

read(101,*)ncm

allocate( rhop(ncm) , hjp(ncm) )

read(101,*)rhop(1:ncm)

if(ncm==1)then
  hjp(ncm) = 1.d300; 
else
  read(101,*)hjp(1:ncm-1)
  hjp(ncm) = 1.d300; 
end if

allocate(z_int(ncm))

if(ncm==1)then
  z_int(1) = 0.d0
else
  z_int(1) = 0.d0
  do i = 2,ncm
    z_int(i) = z_int(i-1) + hjp(i-1)
  end do
end if

nprop = nhet + ncm + 1; 
nprm  = nhet;

allocate( vprop(nprop) , vflprop(nprop) ) 

vprop = 0.d0

do i = 1,nprop 
    vflprop( i ) = i*10 
end do  

read(101,*)mdip

read(101,*)dprm ;  uprm = 1;

 close(101)

end subroutine

!%%%===========================================%%%

subroutine leit_aloc_varglob_FEM
implicit none

integer :: i,j,k,aux,labels
integer :: flbnd = -33

open(unit = 10, file = './malha/'//trim(mname)//'.1.ele', status = 'old', action = 'read' )
open(unit = 20, file = './malha/'//trim(mname)//'.1.node', status = 'old', action = 'read')
open(unit = 30, file = './malha/'//trim(mname)//'.1.face', status = 'old', action = 'read')

read(10,*)nelemg
read(20,*)nnosg
read(30,*)nfaceg

nnoele = 4 ; ! Para o caso do tetraedro

allocate( mat_coord(nnosg,3) , mat_elem(nelemg,nnoele) , mat_face(nfaceg,nnoele) )
allocate( vidprop(nelemg) , rho_elem(nelemg) )

do i = 1,nelemg
  read(10,*)aux,mat_elem(i,:),vidprop(i)
end do
do i = 1,nnosg
  read(20,*)aux,mat_coord(i,:),labels
end do
do i = 1,nfaceg
  read(30,*)aux,mat_face(i,:)
end do

k = 0;
do i = 1,nfaceg
  if(mat_face(i,4)==flbnd)then 
    do j = 1,3
      k = k + 1
    end do
  end if 
end do
nbord = k

allocate(vnobord(nbord))

k = 0;
do i = 1,nfaceg
  if(mat_face(i,4)==flbnd)then
    do j = 1,3
       k = k + 1 ; 
       vnobord(k) = mat_face(i,j)
    end do 
  end if 
end do

deallocate(mat_face)

 close(10)
 close(20)
 close(30)

end subroutine

!%%%===========================================%%%

subroutine desaloca_var_FEM
implicit none

 deallocate( mat_coord , mat_elem , vnobord )
 deallocate( vidprop , rho_elem )

end subroutine

!%%%===========================================%%%

subroutine desaloca_var_1D_var_modelo
implicit none

 deallocate( rhop , hjp )
 deallocate( z_int )

 deallocate( vprop , vflprop ) 

end subroutine


!%%%===========================================%%%

subroutine desaloca_var_perfiz
implicit none

 deallocate( NRPP )
 deallocate( COORCP, INDRCP, INDRCPX )
 deallocate( NTPP )
 deallocate( COOTXP )
 deallocate( freqs )

end subroutine

!%%%===========================================%%%

end module
