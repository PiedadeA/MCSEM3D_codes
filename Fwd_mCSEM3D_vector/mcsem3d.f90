Program FMcsem3d
use vglob
use mcsem_3d
implicit none

include 'mpif.h'

integer              :: ierr,n_p,my_rank
real(db),allocatable :: P(:),Ob(:)
integer              :: no,np,i,j,bnd,nc
real(db),allocatable :: Jc(:,:)
character(20)        :: str
real(db)             :: t1,t2,absl
complex(db),allocatable :: E(:,:), H(:,:)

call MPI_INIT(ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD, n_p , ierr)  
call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr) 

call get_command_argument(1,mname)

!if(n_p.gt.1) call MPI_BARRIER(MPI_COMM_WORLD,ierr)

t1 = MPI_WTIME();

!--
print*,'gerando a relação elemento aresta...'
call edge_element( );
print*,'feito.'

call inicia_parametros()
!--

print*,''
print*,'Numero de variaveis do problema direto:',ngl*nedgst
print*,''

nc = nfreqs * sum(NRPP(:)*NTPP(:))
no = 2 * nc
np = nprm

call arquivos_plot

!--
allocate(nrpft(nfreqs*sum(NTPP)))
open(45,file='dados2.dat',status='old',action='read')
do i = 1,nfreqs*sum(NTPP)
   read(45,*)nrpft(i)
end do
allocate(inrpft(sum(nrpft)))
open(446,file='rrcp.dat',status='old',action='read')
do i = 1,sum(nrpft)
   read(446,*)inrpft(i)
end do
!--

allocate(Exyz(nc,3), Hxyz(nc,3)) ! Variaveis globais.
allocate(P(np),Ob(no))
allocate(Eo(nc))
allocate(ofst(nc))
call ofssets(nc,ofst)
allocate( E(nc,3), H(nc,3) )

!------------------------------------
!-- Calculando o campo primário 1D -- 
!------------------------------------

if(my_rank==0)then

  call primarios(nc,E,H)

  Exyz(:,1) = E(:,1); Exyz(:,2) = E(:,2); Exyz(:,3) = E(:,3);
  Hxyz(:,1) = H(:,1); Hxyz(:,2) = H(:,2); Hxyz(:,3) = H(:,3);
  
  open(unit = 10, file ='./output/Eprimario.dat', status = 'replace', action = 'write' )
  open(unit = 20, file ='./output/Hprimario.dat', status = 'replace', action = 'write' )

  do i = 1,nc
     write(10,'(6E20.10E3)')dreal(Exyz(i,1)),dimag(Exyz(i,1)), &
             dreal(Exyz(i,2)),dimag(Exyz(i,2)), &
             dreal(Exyz(i,3)),dimag(Exyz(i,3))

     write(20,'(6E20.10E3)')dreal(Hxyz(i,1)),dimag(Hxyz(i,1)), &
             dreal(Hxyz(i,2)),dimag(Hxyz(i,2)), &
             dreal(Hxyz(i,3)),dimag(Hxyz(i,3))
  end do 

  close(10)
  close(20)
    
end if

!----------------------------------------------------
!-- Calculando o campo secundário para o modelo 3D -- 
!----------------------------------------------------

P = rhohet

call Foward_MCSEM3D(np,P,no,Ob,n_p,my_rank)

if(my_rank.eq.0)then

  open(unit = 10, file ='./output/Esecundario.dat', status = 'replace', action = 'write' )
  open(unit = 20, file ='./output/Hsecundario.dat', status = 'replace', action = 'write' )

  do i = 1,nc
     write(10,'(6E20.10E3)')dreal(Exyz(i,1)),dimag(Exyz(i,1)), &
             dreal(Exyz(i,2)),dimag(Exyz(i,2)), &
             dreal(Exyz(i,3)),dimag(Exyz(i,3))

     write(20,'(6E20.10E3)')dreal(Hxyz(i,1)),dimag(Hxyz(i,1)), &
             dreal(Hxyz(i,2)),dimag(Hxyz(i,2)), &
             dreal(Hxyz(i,3)),dimag(Hxyz(i,3))
  end do 

  close(10)
  close(20)

end if

!--

t2 = MPI_WTIME();
if(my_rank.eq.0)then

 print*,'Tempo de execucao:'
 if((t2-t1).lt.60.d0)then
   print*,'-----------------------'
   print*,(t2-t1) ,'segundos.'
   print*,'-----------------------'
 elseif((t2-t1).ge.60.d0)then
   print*,'-----------------------'
   print*,(t2-t1)/60.d0 ,'minutos.'
   print*,'-----------------------'
 elseif((t2-t1).ge.3600.d0)then
   print*,'-----------------------'
   print*,(t2-t1)/3600.d0 ,'horas.'
   print*,'-----------------------'
 end if

end if
!--

!--
call finaliza_parametros()
!--

call MPI_FINALIZE(ierr)

end program
