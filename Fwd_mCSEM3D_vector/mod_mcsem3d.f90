!include '/opt/newintel/intel/compilers_and_libraries_2019.1.144/linux/mkl/include/mkl_pardiso.f90'
include 'mkl_pardiso.f90'
module mcsem_3d

use mkl_pardiso
use vglob
use fem_mcsem3d
use lstm
use Arrays_SLU

implicit none

contains

!'''''''''''''''''''''''''''''''''

!-- Rotinas de calculo do campo

!'''''''''''''''''''''''''''''''''

!--

subroutine Foward_MCSEM3D(np,p,no,Ob,n_p,my_rank)
implicit none

include 'mpif.h'

integer,intent(in)   :: np,no,n_p,my_rank
real(db),intent(in)  :: p(np)
real(db),intent(out) :: Ob(no)

real(db)             :: ti,tf
integer,allocatable  :: nrc_ptsm(:,:),NDPF(:)
integer :: i,j,nc
complex(db),allocatable :: Ex(:),Ey(:),Ez(:)
complex(db),allocatable :: Hx(:),Hy(:),Hz(:)

Ob = 0.d0; idjac = 0;

!----------------------------------
!call inicia_parametros (  )
!----------------------------------

allocate(NDPF(nfreqs)); call Estima_num_dados( nc, NDPF ) ; deallocate(NDPF)

call primarios_modelo(dprm,np,nem1,ihet1,Inh,NInh,n_p,my_rank) ; dprm = 0;

allocate(Ex(nc),Ey(nc),Ez(nc))
allocate(Hx(nc),Hy(nc),Hz(nc))
call open_arqprim (  )

 call Campos_MCSEM3D_par( np , nc , P , Ex , Ey , Ez, Hx , Hy , Hz , n_p , my_rank )
 call observacoes ( nc,Ex,no,Ob )

 Exyz(:,1) = Ex;
 Exyz(:,2) = Ey;
 Exyz(:,3) = Ez;

 Hxyz(:,1) = Hx;
 Hxyz(:,2) = Hy;
 Hxyz(:,3) = Hz;

call close_arqprim
deallocate(Ex,Ey,Ez)
deallocate(Hx,Hy,Hz)

call desaloc_vglob
!----------------------------------
!call finaliza_parametros
!----------------------------------

end subroutine

!--

subroutine Campos_MCSEM3D_par( npr , ncp , P , cEx, cEy, cEz, cHx, cHy, cHz, n_p , my_rank )
implicit none

include 'mpif.h'

integer,intent(in)      :: n_p , my_rank
integer,intent(in)      :: npr , ncp
real(db),intent(in)     :: P(npr)
complex(db),intent(out) :: cEx( ncp ), cEy( ncp ), cEz( ncp )
complex(db),intent(out) :: cHx( ncp ), cHy( ncp ), cHz( ncp )

integer                 :: i,j,k,l, cont 
complex(db),allocatable :: sol(:)
complex(db)             :: E_aux , zeta0 , Expm , Eypm, Ezpm, vcc = 0.d0

real(db) :: delx,dely,modu,ti,tf,tii,tff

complex(db),allocatable :: Exyz(:), vftxyz(:) , mL(:)

!-- variaveis pardiso 

INTEGER ::maxfct, mnum, mtype, phase, nnz, error, msglvl, ng, nrhs
INTEGER, ALLOCATABLE :: iparm( : ),perm(:)
INTEGER :: idum(1) , iopt
complex(db):: ddum(1)
TYPE(MKL_PARDISO_HANDLE), ALLOCATABLE  :: pt(:)

!-- Variaveis da Paralelizacao

integer :: fri,frf,ierr
integer,allocatable     :: jlec(:),idsc(:)
complex(db),allocatable :: Vjc(:), Vjcy(:), Vjcz(:)
complex(db),allocatable :: hVjc(:), hVjcy(:), hVjcz(:)
integer :: icp,it
integer :: contt,contr
integer :: aux
!--
complex(db),allocatable :: VjJac(:)
integer,allocatable     :: jlecJac(:), idscJac(:)
integer                 :: contj
complex(db),allocatable :: vjj(:)
!--
complex(db),allocatable  :: Ext_in(:) , Exp_in(:), Eyt_in(:), Ezt_in(:)
complex(db),allocatable  :: Hxt_in(:) , Hxp_in(:), Hyt_in(:), Hzt_in(:)
!--
complex(db)              :: dExdz,dExdy
complex(db)              :: dEydz,dEydx
complex(db)              :: dEzdy,dEzdx
complex(db)              :: iEx,iEy,iEz,auxH,Ema,Eme
real(db)                 :: xr,yr,zr,perc=2.5d0,fW
!--

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

allocate( vjc( sum(NDPF(fri:frf)) ) )  ; vjc = 0.d0;
allocate( vjcy( sum(NDPF(fri:frf)) ) ) ; vjcy = 0.d0;
allocate( vjcz( sum(NDPF(fri:frf)) ) ) ; vjcz = 0.d0;

allocate( hvjc( sum(NDPF(fri:frf)) ) )  ; hvjc = 0.d0;
allocate( hvjcy( sum(NDPF(fri:frf)) ) ) ; hvjcy = 0.d0;
allocate( hvjcz( sum(NDPF(fri:frf)) ) ) ; hvjcz = 0.d0;

deallocate(NDPF)

!-- Vetores iniciais --

ti = MPI_WTIME();
call initial_arrays(cdim,nelemg,nnosg,ngl,nnzi,mat_elem,nedgst,gEdge,col_ptra,row_inda)
tf = MPI_WTIME(); print*,'Tempo arrays inic',(tf-ti)/60.d0,'min'

!--
 cont  = 0; 
 contj = 0; 
 cEx   = (0.d0,0.d0);
 cEy   = (0.d0,0.d0);
 cEz   = (0.d0,0.d0);
 cHx   = (0.d0,0.d0);
 cHy   = (0.d0,0.d0);
 cHz   = (0.d0,0.d0);

do i = fri,frf;
     
     ifrq = i
     Freq = Freqs(i)
     zeta0 = ip * 2.d0*pi*Freq * mi0

     nrs = 1

     allocate(vftxyz(ngl*nedgst)) ; vftxyz = (0.d0,0.d0);

     ti = MPI_WTIME();

     nnz=nnzi
     allocate( val_nz(nnz),row_ind(nnz),col_ptr(ngl*nedgst+1) )    
     row_ind=row_inda ; col_ptr=col_ptra;
     val_nz=0.d0;

     print*,'nnz antes:',nnz

     call Montagem_Matriz_Global_tetraedro(ngl , nnosg , nelemg , nnz , val_nz)
     call recount_arrays(ngl,nnosg,nedgst,row_ind,col_ptr,val_nz,nnz)
     call condicao_contorno(cdim , ngl , nbord , vnobord , vcc , col_ptr , nnz , val_nz , nnosg, nedgst , Vftxyz)
     !--

     print*,'nnz depois:',nnz

     tf = MPI_WTIME();

     ngb = ngl*nedgst;

     print*,'Matriz global montada.',(tf-ti)/60.d0,'min'

     !-- Pardiso

     ng = ngl*nedgst;

     nrhs = nrs
     maxfct = 1
     mnum = 1
     allocate( iparm(64), perm(ngl*nedgst), Exyz(ngl*nedgst) )

     Exyz = (0.d0,0.d0);

     iparm    = 0

     perm = 0
     msglvl = 1
     allocate( pt(64) )
     mtype = 6

     print*,'Memoria Antes:'
     call inf_memory (my_rank)  

     call pardisoinit(pt,mtype,iparm)

     ti = MPI_WTIME()

     phase = 11 ! only reordering and symbolic factorization
     CALL pardiso (pt, maxfct, mnum, mtype, phase, ng, val_nz, col_ptr, row_ind, &
	              perm, nrhs, iparm, msglvl, vftxyz, Exyz, error)
     if(error.ne.0) then
          print*,'Erro na fat. simb.';stop
     end if

     print*,'Fatorando...'

     !.. Factorization.
     phase = 22 ! only factorization
     CALL pardiso (pt, maxfct, mnum, mtype, phase,ng, val_nz, col_ptr, row_ind, &
	              perm, nrhs, iparm, msglvl, vftxyz, Exyz, error)
     if(error.ne.0) then
          print*,'Erro na fat.';stop
     end if

     print*,'Memoria Depois:'
     call inf_memory (my_rank)  

     tf = MPI_WTIME();
     print*,''
     if(my_rank==0)print*,'$> Tempo de fatoracao',(tf-ti)/60.d0,'min'

     !-- Fatoracao Matriz Global --
     !allocate(mL(nnz))
     !call cCholesk2( ng,nnz,col_ptr,row_ind,val_nz,mL )
     !--

     contt = 0;

     do l = 1,Num_perf
        iperf = l;

       do j = 1,NTPP(l);          
 
          contt = contt+1 

          tr_x = COOTXP(contt,1) 
          tr_y = COOTXP(contt,2) 
          tr_z = COOTXP(contt,3)

          itrm = contt           

          vftxyz = (0.d0,0.d0);

          call Montagem_Vetor_Fonte_tetraedro( ngl , nelemg , nnosg , vftxyz )
          call condicao_contorno( cdim , ngl , nbord , vnobord , vcc , col_ptr , nnz , val_nz , nnosg, nedgst , Vftxyz )

          !call cpchlesky_cgrad_solv( ng,nnz,col_ptr,row_ind,val_nz,Vftxyz,mL,Exyz )

          nrhs = nrs
          phase = 33 ! system solution
          call pardiso (pt, maxfct, mnum, mtype, phase,ng, val_nz, col_ptr, row_ind, &
                        perm, nrhs, iparm, msglvl, vftxyz, Exyz, error)
          if(error.ne.0) then
             print*,'Erro na solucao so sistema' ; stop
          end if

	 !--
	  call Numrec_portransm( l , i , j , tr_x , tr_y , nrcp , ind_perf, ind_edg ); 

	  allocate(coor_perf(nrcp,3));

	  coor_perf(:,1) = COORCP(ind_perf(:),1)
          coor_perf(:,2) = COORCP(ind_perf(:),2)
          coor_perf(:,3) = COORCP(ind_perf(:),3)

	 !--
 
	  allocate(Ext_in(nrcp) , Exp_in(nrcp))
	  allocate(Eyt_in(nrcp) , Ezt_in(nrcp))
          Ext_in = (0.d0,0.d0) ; Exp_in = (0.d0,0.d0); 
          Eyt_in = (0.d0,0.d0) ; Ezt_in = (0.d0,0.d0); 

	  allocate(Hxt_in(nrcp) , Hxp_in(nrcp))
	  allocate(Hyt_in(nrcp) , Hzt_in(nrcp))
          Hxt_in = (0.d0,0.d0) ; Hxp_in = (0.d0,0.d0); 
          Hyt_in = (0.d0,0.d0) ; Hzt_in = (0.d0,0.d0); 

          do k = 1,nrcp
 
                xr = coor_perf(k,1)
                yr = coor_perf(k,2)

                call prim_MCSEM( ttap , Freq , tr_x , tr_y , tr_z , & 
                                 xr,yr,coor_perf(k,3), &
                                 Expm , Eypm, Ezpm )

                if(interp==1)then
                   call compose_Esec2( coor_perf(k,1),coor_perf(k,2),coor_perf(k,3), Ext_in(k), Eyt_in(k), Ezt_in(k), &
                                    Hxt_in(k), Hyt_in(k), Hzt_in(k), Exyz )
                elseif(interp==2)then

                   call compose_Esec3( coor_perf(k,1),coor_perf(k,2),coor_perf(k,3), Ext_in(k), Eyt_in(k), Ezt_in(k), &
                                       Hxt_in(k), Hyt_in(k), Hzt_in(k), Exyz )
                else
                   print*,'Escolha errada do interpolador';stop
                end if

                Ext_in(k) = Ext_in(k) !+ Expm ;
                Eyt_in(k) = Eyt_in(k) !+ Eypm ; 
                Ezt_in(k) = Ezt_in(k) !+ Ezpm ; 

                Hxt_in(k) = Hxt_in(k) !+ cpHx ; 
                Hyt_in(k) = Hyt_in(k) !+ cpHy ; 
                Hzt_in(k) = Hzt_in(k) !+ cpHz ; 

          end do
          !--

          !--
	  do k = 1,nrcp
             cont = cont + 1
	     vjc(cont)  = Ext_in(k)
	     vjcy(cont) = Eyt_in(k)
	     vjcz(cont) = Ezt_in(k)

	     hvjc(cont)  = Hxt_in(k)
	     hvjcy(cont) = Hyt_in(k)
	     hvjcz(cont) = Hzt_in(k)
	  end do
          !--

	  deallocate(Ext_in, Exp_in)
	  deallocate(Eyt_in, Ezt_in)

	  deallocate(Hxt_in, Hxp_in)
	  deallocate(Hyt_in, Hzt_in)

          deallocate( ind_perf, ind_edg, coor_perf )

       end do ! Transmissor

     end do ! Perfis

     nrhs  = nrs
     iopt  = 3
     phase = -1 ! release internal memory
     CALL pardiso (pt, maxfct, mnum, mtype, phase,ng, val_nz, col_ptr, row_ind, &
 		   perm, nrhs, iparm, msglvl, vftxyz, Exyz, error)

     !deallocate(mL)
     deallocate(col_ptr)
     deallocate(val_nz)
     deallocate(row_ind)
     deallocate(vftxyz,Exyz) 
     deallocate(iparm,perm,pt)

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

!--
deallocate(col_ptra,row_inda)
!--

deallocate(jlec,idsc)
deallocate(vjc)
deallocate(vjcy,vjcz)

deallocate(hvjc)
deallocate(hvjcy,hvjcz)

end subroutine

!--

!'''''''''''''''''''''''''''''''''

!-- Leitura do modelo e de seus 
!   parametros

!'''''''''''''''''''''''''''''''''

!--

subroutine inicia_parametros (  )
implicit none

integer :: no

call leitura_aloc_perfiz
call leitura_prop_modelo
call leit_aloc_propmod_var_1D
call leit_aloc_varglob_FEM

!--
no = nfreqs * sum(NRPP(:)*NTPP(:))
allocate(Flagf(no))
Flagf = 1
!--

end subroutine

!--

subroutine finaliza_parametros

call desaloca_var_perfiz
call desaloca_var_1D_var_modelo
call desaloca_var_FEM

deallocate(Flagf)

end subroutine

!--

subroutine leitura_aloc_perfiz()
implicit none

integer  :: i,j,k,cont
real(db) :: aux
integer  :: Nrct, Nftt
 
open(unit = 60, file ='./parametros_modelo/coord_Rx_perfis.in', status = 'old', action = 'read')
open(unit = 70, file ='./parametros_modelo/coord_Tx_perfis.in', status = 'old', action = 'read')

read(60,*)Num_perf

allocate( NRPP(Num_perf) )

do i = 1,Num_perf
          read(60,*)NRPP(i)
          do j = 1,NRPP(i)
		read(60,*)aux
          end do  
end do 

Nrct = sum( NRPP(:) )

rewind(60)

allocate( COORCP(Nrct,3), INDRCP(Nrct) ) 

read(60,*)aux
cont = 0
do i = 1,Num_perf
          read(60,*)aux
          do j = 1,NRPP(i)
                cont = cont + 1		

 		INDRCP(cont) = cont
		read(60,*)COORCP(cont,1),COORCP(cont,2),COORCP(cont,3)

          end do  
end do

close(60) 

read(70,*)Num_perf

allocate( NTPP(Num_perf) )

do i = 1,Num_perf
       read(70,*)NTPP(i)
       do j = 1,NTPP(i)
          read(70,*)Aux
       end do
end do

NFTT = sum( NTPP(:) )  

allocate( COOTXP(NFTT,3) )
allocate( aTx(NFTT) )

rewind(70)

read(70,*)Num_perf

aTx = 0.d0

 cont = 0
do i = 1,Num_perf
       read(70,*)NTPP(i)
       do j = 1,NTPP(i)
          cont = cont + 1
          read(70,*)COOTXP(cont,1),COOTXP(cont,2),COOTXP(cont,3)
       end do
end do

close(70)

!print*,'passou1'

end subroutine

!--

subroutine leitura_prop_modelo( )
implicit none

integer  :: i
real(db) :: aux

open(unit = 101, file ='./parametros_modelo/heterogeneidades.in', status = 'old', action = 'read' )

read(101,*)nhet

allocate(rhohet(nhet))

read(101,*)rhohet(1:nhet)

read(101,*)rad,interp

close(101)

end subroutine

!--

subroutine leit_aloc_propmod_var_1D
implicit none

integer  :: i,j
real(db) :: aux

open(unit = 101, file ='./parametros_modelo/frequencias_modelo1d.in', status = 'old', action = 'read' )

read(101,*)nFreqs
allocate(freqs(nfreqs))
read(101,*)freqs

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

read(101,*) mdip

read(101,*) dprm; uprm = 1;

close(101)

end subroutine

!--

subroutine leit_aloc_varglob_FEM
implicit none

integer :: i,j,k,aux,labels
integer :: flbnd = -1
integer,allocatable :: bd(:),lnd(:)

open(unit = 10, file = './malha/'//trim(mname)//'.1.ele', status = 'old', action = 'read' )
open(unit = 20, file = './malha/'//trim(mname)//'.1.node', status = 'old', action = 'read')
open(unit = 30, file = './malha/'//trim(mname)//'.1.face', status = 'old', action = 'read')
open(unit = 40, file = './malha/'//trim(mname)//'.1.nedge', status = 'old', action = 'read')
open(unit = 50, file = './malha/'//trim(mname)//'.1.edge', status = 'old', action = 'read')

read(10,*)nelemg
read(20,*)nnosg
read(30,*)nfaceg
read(40,*)aux
read(50,*)nedgst

allocate( mat_coord(nnosg,3) , mat_elem(nelemg,nnoele) , mat_face(nfaceg,nnoele) )
allocate( vidprop(nelemg) , rho_elem(nelemg) )
allocate( gEdge(nelemg,nedgs) )
allocate( node_edge(nedgst,2),bd(nedgst) ) ; bd = 0;
allocate(lnd(nnosg)) ; lnd = 0;

allocate(fbd(nnosg))

do i = 1,nedgst
  read(50,*)aux,node_edge(i,:),bd(i)
end do
do i = 1,nelemg
  read(10,*)aux,mat_elem(i,:),vidprop(i)
  read(40,*)gEdge(i,:)
end do
do i = 1,nnosg
  read(20,*)aux,mat_coord(i,:),lnd(i)
end do

do i = 1,nfaceg
  read(30,*)aux,mat_face(i,:)
end do

!--
lnd = 0
do i = 1,nfaceg
   if(mat_face(i,4)==lbnd)then
     do j = 1,3     
        lnd(mat_face(i,j)) = lbnd
     end do
   end if   
end do
fbd = lnd;
!--

k = 0;
do i = 1,nedgst
   if( (lnd(node_edge(i,1))==lbnd).and.(lnd(node_edge(i,2))==lbnd) )then
        k = k + 1
   end if
end do

nbord = k
!print*,'arestas na borda:',nbord

allocate(vnobord(nbord))

k = 0;
do i = 1,nedgst
   if( (lnd(node_edge(i,1))==lbnd).and.(lnd(node_edge(i,2))==lbnd) )then
      k = k + 1
      vnobord(k) = i
   end if 
end do

deallocate(mat_face)
deallocate(bd)
deallocate(lnd)

 close(10)
 close(20)
 close(30)
 close(40)
 close(50)

end subroutine

!--

subroutine desaloca_var_perfiz
implicit none

 deallocate( NRPP )
 deallocate( COORCP, INDRCP )
 deallocate( NTPP )
 deallocate( COOTXP )
 deallocate( aTx )

end subroutine

!--

subroutine desaloca_var_1D_var_modelo
implicit none

 deallocate( rhop , hjp )
 deallocate( z_int )

 deallocate( vprop , vflprop ) 

 deallocate( freqs )

end subroutine

!--

subroutine desaloca_var_FEM
implicit none

 deallocate( mat_coord , mat_elem , vnobord )
 deallocate( vidprop , rho_elem )
 deallocate( gEdge )
 deallocate( node_edge )
 deallocate( fbd )

end subroutine

!--

subroutine desaloc_vglob
implicit none

deallocate( ihet1,Inh )

end subroutine

!--

!'''''''''''''''''''''''''''''''''

!-- Rotinas extras

!'''''''''''''''''''''''''''''''''

!--

subroutine Estima_num_dados( num_fields, NDPF ) 
implicit none
integer,intent(out) :: num_fields,NDPF(:)

integer :: i,j,k,l,c1,c2,c3,c4
real(db) :: dx,dy,r

 num_fields = 0

!--
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
!--

num_fields = c2 ! Depois de uma recontagem aqui e atribuido o valor correto do numero de campos.

end subroutine

!--

subroutine Numrec_portransm( lp, lp_frq , trsm , xtr , ytr , nrc, ind_rec, ind_edg )
implicit none

integer,intent(in)  :: lp , lp_frq , trsm
real(db),intent(in) :: xtr , ytr

integer,intent(out) :: nrc
integer,allocatable,intent(out) :: ind_rec(:), ind_edg(:)

integer,allocatable :: ind_reca(:)

integer :: i,j,k,l,c1,c2,c3,c4
real(db) :: dx,dy,r
integer :: lp_perf,aux

integer :: edg

allocate(ind_reca(NRPP(lp)))

allocate(NDPF(nfreqs)); call Estima_num_dados( aux, NDPF ) 

if(lp_frq.gt.1)then
 c1 = (lp_frq-1)*sum(NTPP(:)*NRPP(:)); 
else
 c1 = 0; 
end if

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

                       if(i==1)then
	                 ind_reca(c2) = INDRCP( k )
                       else
	                 ind_reca(c2) = INDRCP( sum(NRPP(1:i-1)) + k )
                       end if

	           end if
                 end do

             end if

        end do


     end if

end do

nrc = c2;

allocate(ind_rec(nrc) , ind_edg(nrc))

do i = 1,nrc
   ind_rec(i) = ind_reca(i)
   ind_edg(i) = ind_reca(i)
end do


deallocate(ind_reca)

deallocate(NDPF)

end subroutine

!--

subroutine ind_edge( edg,no )
implicit none
integer,intent(in)  :: no
integer,intent(out) :: edg

integer :: ind,i,j
integer :: n1,n2

do i = 1,sum(NRPP)
   if(no==inds_Rx_edg(i,1))then
      ind = i
      exit
   end if
end do

n1 = inds_Rx_edg(ind,1)
n2 = inds_Rx_edg(ind,2)

do i = 1,nedgst
   if( (n1==node_edge(i,1)).and.(n2==node_edge(i,2)) )then
      edg = i
      exit
   end if
   if( (n2==node_edge(i,1)).and.(n1==node_edge(i,2)) )then
      edg = i
      exit
   end if
end do

end subroutine

!--

!'''''''''''''''''''''''''''''''''

!-- Outras rotinas

!'''''''''''''''''''''''''''''''''

!--

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

!--

subroutine arquivos_plot
implicit none

integer :: i,j,k,l,c1,c2,c3
real(db) :: dx,dy,r

 open(44,file='dados1.dat',status='replace',action='write')
 open(45,file='dados2.dat',status='replace',action='write')
 open( 446,file='rrcp.dat',status='replace',action='write')

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

		    if( flagF(c3).eq.1 )then 
                        if(i==1)then
                           dx = COORCP((k),1) - COOTXP(j,1)
                           dy = COORCP((k),2) - COOTXP(j,2)
                           r  = dsqrt( dx**2 + dy**2 )
                           write(446,*)r 
                        else 
                           dx = COORCP((sum(NRPP(1:i-1)) + k),1) - COOTXP(sum(NTPP(1:i-1)) + j,1)
                           dy = COORCP((sum(NRPP(1:i-1)) + k),2) - COOTXP(sum(NTPP(1:i-1)) + j,2)
                           r  = dsqrt( dx**2 + dy**2 )
                           write(446,*)r
                        end if
                    end if

	      end do
	            write(45,*) c2 
	  end do
	end do
end do
!&&&&&&&&&&&&&&&&&&&&

 close(44)
 close(45)
 close(446)

end subroutine

!--

subroutine observacoes ( nc,E,no,Obs )
implicit none
integer,intent(in)     :: nc,no
complex(8),intent(in)  :: E(nc)
real(8),intent(out)    :: Obs(no)

integer  :: i,j,nr,cont,cnt,ii
real(8),allocatable :: phi(:),aphi(:),aux(:)
real(8),allocatable :: p1(:), p2(:),offst(:)
integer             :: imn

allocate(phi(nc),aphi(nc))

cont = 0
!-- absoluto 
do i = 1,nc
  cont = cont + 1
  Obs(cont)=dlog10(cdabs(E(i)))
end do
!-- fase
do i = 1,nc
  phi(i) = datan2(dimag(E(i)),dreal(E(i)))
end do
!--

 cnt=0

do i = 1,nfreqs*sum(NTPP)

  allocate(aux(nrpft(i)))
  allocate(offst(nrpft(i)))

  if(i==1)then
    do j = 1,nrpft(i)
       aux(j)   = phi(j)
       offst(j) = inrpft(j)
    end do
  else
    do j = 1,nrpft(i)
       aux(j)   = phi(sum(nrpft(1:i-1)) + j)
       offst(j) = inrpft(sum(nrpft(1:i-1)) + j)
    end do
  end if

  imn = minloc(offst,1)

  allocate( p1(imn-1),p2(nrpft(i)-(imn-1)) )

  ii=0
  do j = imn-1,1,-1
     ii=ii+1
     p1(ii) = aux(j)
  end do
  ii=0
  do j = imn,nrpft(i)
     ii=ii+1
     p2(ii) = aux(j)
  end do

  call unwrap(p1,imn-1)
  call unwrap(p2,nrpft(i)-(imn-1))

  ii=0
  do j = imn-1,1,-1
     ii=ii+1
     aux(ii) = p1(j)
  end do
  do j = 1,nrpft(i)-(imn-1)
     ii=ii+1
     aux(ii) = p2(j)
  end do
  
  do j = 1,nrpft(i)
     cnt=cnt+1
     aphi(cnt) = aux(j)
  end do
 
  deallocate(p1,p2)
  deallocate(aux)
  deallocate(offst)
end do

phi = aphi 

do i = 1,nc
  cont = cont + 1
  Obs(cont) = phi(i)
end do


end subroutine

!--

subroutine unwrap(phase,n) 
!implicit none

!c---------------------------------------------------------------------- 
!c  Routine unwrap:To unwrap the phase function computed by 
!c                 fuction ATAN(HR,HI); 
!c input parameters: 
!c   PHASE:  n dimensioned real array. PHASE(0) to PHASE(n-1) are the 
!c           phase value to be unwraped. Note: the phase is in radian. 
!c    N   :  the dimension of array PHASE. 
!c output parameters: 
!c   PHASE:  n dimensioned real array, PHASE(0) to PHASE(n-1) are the 
!c           phase value have been unwraped; 
!c                                      in Chapter 2 
!c---------------------------------------------------------------------- 

integer,intent(in)   :: n
real(8),intent(out) :: phase

integer  :: k 
real(8) :: angle, anglejump, dx 

        dimension phase(0:n-1) 
        angle=4.*datan(1.d0)
        anglejump=0.d0 
        do 20 k=1,n-1 
           dx=phase(k)-(phase(k-1)-anglejump) 
           if(dabs(dx).le.angle) go to 10 
           anglejump=anglejump-sign(2.*angle,dx) 
10      phase(k)=phase(k)+anglejump 
20      continue 
        return 

!        end 

end subroutine unwrap

!--

!'''''''''''''''''''''''''''''''''

!-- Gerador de malha

!'''''''''''''''''''''''''''''''''

!--

subroutine edge_element( )

implicit none
  
  character(len=100) :: comand
  integer :: ierr

  call execute_command_line ("gfortran gera_arestas.f90 -o gedge.x", exitstat=ierr)
  call execute_command_line ("./gedge.x "//trim(mname), exitstat=ierr)

end subroutine

!--

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

!--

subroutine noise(miperc,maperc,minp,maxp,rpf,Nd,Dt)
!subroutine noise(miperc,maperc,minp,maxp,xtr,ytr,ixp,Nd,Dt)
implicit none

integer,intent(in)   :: Nd
real(8),intent(in) :: rpf(Nd)
!integer,intent(in)   :: Nd,ixp(:)
real(8),intent(in) :: miperc,maperc,minp,maxp
!real(dpc),intent(in) :: miperc,maperc,xtr,ytr,minp,maxp

complex(8),intent(inout) :: Dt(Nd)

real(8),allocatable :: re(:),im(:)
integer :: i
real(8) :: a,b,Nor,Noi,rp

allocate(re(Nd),im(Nd))

a = ( (maperc-miperc)/100.d0 )/( maxp-minp )
b = ( maperc/100.d0 ) - a * maxp

re=0.d0;im=0.d0
do i = 1,Nd

 call random_number(Nor)
 call random_number(Noi)
 Nor = 2.d0*Nor - 1.d0
 Noi = 2.d0*Noi - 1.d0
 Nor = Nor/dabs(Nor)
 Noi = Noi/dabs(Noi)

 rp = rpf(i)
 !rp = dsqrt( (mat_coord(ixp(i),1)-xtr)**2 + (mat_coord(ixp(i),2)-ytr)**2 )
 
 re(i) = (a * rp + b)*Nor
 im(i) = (a * rp + b)*Noi

end do

Dt = cmplx( dreal(Dt)*( 1.d0 + re ) , dimag(Dt)*( 1.d0 + im ) )

deallocate(re,im)

end subroutine

!--

subroutine primarios(ncp,E,H)
implicit none
integer,intent(in) :: ncp
complex(db),intent(out) :: E(ncp,3),H(ncp,3)

integer   :: i,j,k,p
real(db) :: r,x,y,z
real(db) :: xt,yt,zt,f

integer :: ncpf,it,ix,ii,iat
real(db),allocatable :: off(:)
complex(db) :: Ex,Ey,Ez, E1, E2, E3

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

      xt = COOTXP(it,1) ; yt = COOTXP(it,2) ; zt = COOTXP(it,3) ; 

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
         H(ii,1) = cpHx; H(ii,2) = cpHy; H(ii,3) = cpHz

      end do
   end do
 end do
end do

! close (20)

end subroutine

!--

subroutine ofssets(ncp,ofs)
implicit none
integer,intent(in)    :: ncp
real(db),intent(out) :: ofs(ncp)


integer   :: i,j,k
real(db) :: r,x,y,z
real(db) :: xt,yt,zt

integer :: ncpf,it,ix,ii
real(db),allocatable :: off(:)

ncpf = sum(NTPP(:)*NRPP(:))

allocate(off(ncpf))

ii=0
do i = 1,Num_perf
   do j = 1,NTPP(i)

      if(i==1)then
        it = j
      else
        it = sum(NTPP(1:i-1)) + j
      end if

      xt = COOTXP(it,1) ; yt = COOTXP(it,2) ; zt = COOTXP(it,3) ; 

      do k = 1,NRPP(i)

         if(i==1)then
           ix = k
         else
           ix = sum(NRPP(1:i-1)) + k
         end if

         x = COORCP(ix,1); y = COORCP(ix,2); z = COORCP(ix,3); 
         
         r = dsqrt( (x-xt)**2 + (y-yt)**2 + (z-zt)**2)

         if ( x < xt )then      
            r = -r
         else
            r =  r
         end if
         ii = ii + 1
         off(ii) = r
      end do
   end do
end do

open(234,file='offsets.dat',status='replace',action='write')

ii = 0
do i = 1,nfreqs
   do j = 1,ncpf
      ii = ii + 1
      ofs(ii) = off(j)
      if(flagF(ii)==1)then
         write(234,'(F20.10,2x)') off(j)
      end if
   end do
end do

ii = 0
do i = 1,nfreqs
   do j = 1,ncpf
      ii = ii + 1
      if(flagF(ii)==1)then
         write(234,'(F20.10,2x)') off(j)
      end if
   end do
end do

close(234)

end subroutine

end module
