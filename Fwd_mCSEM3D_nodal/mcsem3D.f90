program mcsem3D
   use var_glob
   use fem_3D
   use campo_MCSEM3D
   use modulo_auxiliar
implicit none

include 'mpif.h'

integer  :: ierr,n_p,my_rank

integer :: i,j,em,nf,np,ntr,ncp,eesc, bnd
real(dpc),allocatable :: trm(:,:),Fr(:),P(:),phaseEx(:)
complex(dpc),allocatable :: cEx(:), cEy(:), cEz(:)
complex(dpc),allocatable :: cHx(:), cHy(:), cHz(:)
complex(dpc),allocatable :: E(:,:), H(:,:)

real(dpc) :: t1,t2
character(len=100) :: path,comand,filename

real(dpc),allocatable :: jac(:,:) , Obs(:)
integer               :: no

call MPI_INIT(ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD, n_p , ierr)  
call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr) 

!## Criando a malha do modelo ##

call get_command_argument(1,mname)

if(n_p.gt.1) call MPI_BARRIER(MPI_COMM_WORLD,ierr)

!## --

t1 = MPI_WTIME();

!----------------------------------
call inicia_parametros ()
!----------------------------------

print*,''
print*,'Numero de variaveis do problema direto:',ngl*nnosg

call arquivos_plot

allocate(NDPF(nfreqs)); call Estima_num_dados( ncp, NDPF ) ; deallocate(NDPF)

np = nprm;

allocate( cEx(ncp), cEy(ncp), cEz(ncp))
allocate( cHx(ncp), cHy(ncp), cHz(ncp))
allocate( P(np) )
allocate( E(ncp,3), H(ncp,3) )


!------------------------------------
!-- Calculando o campo primário 1D -- 
!------------------------------------

if(my_rank==0)then

  call primarios(ncp,E,H)

  cEx = E(:,1); cEy = E(:,2); cEz = E(:,3);
  cHx = H(:,1); cHy = H(:,2); cHz = H(:,3);
  
  open(unit = 10, file ='./output/Eprimario.dat', status = 'replace', action = 'write' )
  open(unit = 20, file ='./output/Hprimario.dat', status = 'replace', action = 'write' )

  do i = 1,ncp
     write(10,'(6E20.10E3)')dreal(cEx(i)),dimag(cEx(i)), &
                dreal(cEy(i)),dimag(cEy(i)), &
                dreal(cEz(i)),dimag(cEz(i))

     write(20,'(6E20.10E3)')dreal(cHx(i)),dimag(cHx(i)), &
                dreal(cHy(i)),dimag(cHy(i)), &
                dreal(cHz(i)),dimag(cHz(i))
  end do 

  close(10)
  close(20)
    
end if

!----------------------------------------------------
!-- Calculando o campo secundário para o modelo 3D -- 
!----------------------------------------------------

P  = rhohet

call primarios_modelo(dprm,np,nem1,ihet1,Inh,NInh,n_p,my_rank) ; dprm = 0;

call open_arqprim (  )

call Campos_MCSEM3D_par( np , ncp , P , cEx , cEy, cEz, cHx, cHy, cHz, n_p , my_rank );

call close_arqprim

call desaloc_vglob

if(my_rank.eq.0)then

      open(unit = 10, file ='./output/Esecundario.dat', status = 'replace', action = 'write' )
      open(unit = 20, file ='./output/Hsecundario.dat', status = 'replace', action = 'write' )

      do i = 1,ncp
         write(10,'(6E20.10E3)')dreal(cEx(i)),dimag(cEx(i)), &
                dreal(cEy(i)),dimag(cEy(i)), &
                dreal(cEz(i)),dimag(cEz(i))

         write(20,'(6E20.10E3)')dreal(cHx(i)),dimag(cHx(i)), &
                dreal(cHy(i)),dimag(cHy(i)), &
                dreal(cHz(i)),dimag(cHz(i))
      end do 

      close(10)
      close(20)

end if

!----------------------------------
call finaliza_parametros
!----------------------------------

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

call MPI_FINALIZE(ierr)

end program
