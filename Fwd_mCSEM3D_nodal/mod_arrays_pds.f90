!include 'mkl_pardiso.f90'
module arrayspds

! Modulo que constroi os vetores no formato de linha suprimida, para serem usados nas rotinas do Pardiso.
! Autor: Anderson A. Piedade. 26/12/2017

integer,parameter :: dpcs=kind(1.d0) 

contains

!************************************

subroutine initial_arrays(cdim,nel,nnos,ngl,nnz,matel,npr,ic)
implicit none
integer,intent(in)  :: cdim,ngl,nel,nnos
integer,intent(out) :: nnz
integer,intent(in)  :: matel(:,:)
integer,allocatable,intent(inout) :: npr(:),ic(:)

integer,allocatable :: bnds(:),cnts(:),cntsin(:),vnn(:)
integer :: i,j,k,l,bnd,nnzi,j_f,i_f,dif,qi
integer,allocatable :: ica(:),el(:),ml(:,:)

integer,allocatable :: indcol(:,:),cnt(:),ics(:),icp(:)

nnz=0
allocate(bnds(ngl*nnos)) ; bnds=0
allocate(cnts(ngl*nnos)) ; cnts=0
allocate(cnt(ngl*nnos)) ; cnt=0

if(cdim==2)then
    allocate(ml(ngl*3,ngl*3),el(ngl*3)) ; ml=1 ; el=0
    i_f = ngl*3 ; j_f=i_f
elseif(cdim==3)then
    allocate(ml(ngl*4,ngl*4),el(ngl*4)) ; ml=1 ; el=0
    i_f = ngl*4 ; j_f=i_f
else
  print*,'Escolha a dimensao correta, 2 ou 3 !'
  stop
end if

do k =1,nel

    if( ngl==1 .and. ( cdim==2 .or. cdim==3 ) )then
      el = matel(k,:)
    elseif( cdim==2 .and. ngl==2 )then
      do j=1,3
         el( 2*j - 1 )= 2*matel(k,j)-1 ; el(2*j)= 2*matel(k,j)
      end do
    elseif( cdim==3 .and. ngl==4 )then
      do j=1,4
         el( 4*j - 3 )= 4*matel(k,j)-3 ; el(4*j - 2)= 4*matel(k,j)-2
         el( 4*j - 1 )= 4*matel(k,j)-1 ; el(4*j)= 4*matel(k,j)
      end do
    end if

    do i = 1,i_f
       do j = 1,j_f   

	   cnts(el(i)) = cnts(el(i)) + 1

       end do
    end do 

end do

if( cdim==2 )then
  qi = ngl*2*maxval(cnts/3)
  allocate(indcol(ngl*nnos,qi))
elseif( cdim==3 )then
  qi = ngl*3*maxval(cnts/3)
  allocate(indcol(ngl*nnos,qi))
end if

 cnt=0;
do k =1,nel

    if( ngl==1 .and. ( cdim==2 .or. cdim==3 ) )then
      el = matel(k,:)
    elseif( cdim==2 .and. ngl==2 )then
      do j=1,3
         el( 2*j - 1 )= 2*matel(k,j)-1 ; el(2*j)= 2*matel(k,j)
      end do
    elseif( cdim==3 .and. ngl==4 )then
      do j=1,4
         el( 4*j - 3 )= 4*matel(k,j)-3 ; el(4*j - 2)= 4*matel(k,j)-2
         el( 4*j - 1 )= 4*matel(k,j)-1 ; el(4*j)= 4*matel(k,j)
      end do
    end if

    do i = 1,i_f

       do j = 1,j_f   

           if(el(j) > el(i))then
	        cnt(el(i)) = cnt(el(i)) + 1
		indcol(el(i),cnt(el(i))) = el(j)
           end if

       end do
    end do 

end do

nnzi= ngl*nnos
do i = 1,ngl*nnos-1
   l = 0
   do j = 1,qi
     if( indcol(i,j).eq.0 )then
        exit
     else
        nnzi = nnzi + 1 ; l = l + 1;
     end if
   end do
   bnds(i) = l
end do

print*,''
print*,'Chute inicial:',nnzi

allocate(npr(ngl*nnos+1))
allocate(ica(nnzi))

npr(1)=1

nnz = 0
do i = 1,ngl*nnos-1

   nnz = nnz + 1
   ica(nnz) = i

   allocate(ics(bnds(i)))

   ics = indcol( i,1:bnds(i) )

   call oerdena(bnds(i),ics)
  
   l = bnds(i)

   call recount(l,ics)

   npr(i+1) = npr(i)+l+1
 
   do j=1,l
      nnz=nnz+1    
      ica(nnz)=ics(j)
   end do

   deallocate(ics)
   
end do

nnz=nnz+1 ; ica(nnz)=ngl*nnos;
npr(ngl*nnos + 1)= npr(ngl*nnos) + 1

allocate(ic(nnz))

ic=ica(1:nnz)

print*,'Primeira estimativa de nnz:',nnz
print*,''

deallocate(ica)
deallocate(indcol)

deallocate(el,ml)
deallocate(bnds)
deallocate(cnts,cnt)

end subroutine

!************************************

subroutine recount(n,icl)
implicit none

integer,intent(inout) :: n
integer,allocatable,intent(inout) :: icl(:)

integer :: i,j,cnt
integer,allocatable :: aux(:)

if( n==1 )then
  n=1
  icl(n)=icl(n)
elseif( n==2 )then
  if( icl(1)==icl(2) )then
    allocate(aux(1))
    aux = icl(1)
    deallocate(icl)
    allocate(icl(1))
    icl=aux ; n=1;
    deallocate(aux)
  else
    icl=icl; n=2;
  end if
else
  allocate(aux(n)); aux=0;

  cnt=0
  do i=1,n-1
        if( icl(i) /= icl(i+1) )then
          if( (i+1)==n )then
            cnt = cnt+1 ; aux(cnt)=icl(i)         
            cnt = cnt+1 ; aux(cnt)=icl(i+1)
          else
	    cnt=cnt+1 ; aux(cnt)=icl(i)	  
	  end if
	else
          if( (i+1)==n )then
            cnt = cnt+1 ; aux(cnt)=icl(i)
	  end if
        end if
  end do

  deallocate(icl)
  n=cnt  
  allocate(icl(n))
  icl=aux(1:n);
  deallocate(aux);  
end if

end subroutine

!************************************

subroutine initial_arrays2(cdim,nel,nnos,ngl,nnz,matel,npr,ic)
implicit none
integer,intent(in)  :: cdim,ngl,nel,nnos
integer,intent(out) :: nnz
integer,intent(in)  :: matel(:,:)
integer,allocatable,intent(inout) :: npr(:),ic(:)

integer,allocatable :: bnds(:),bndst(:),ini(:)
integer :: i,j,k,l,bnd,nnzi,j_f,i_f,dif,ift,cnt,aux,cnt2,mtp
integer,allocatable :: ica(:),el(:),ml(:,:)

if(cdim==2)then
 mtp=6
elseif(cdim==3)then
 mtp=25
end if

ift = ngl*(nnos+mtp*nel);

allocate(bnds(ngl*nnos),bndst(ngl*nnos)) 
allocate( ica(ift) )

print*,'chute inicial:',ift

if(cdim==2)then
    allocate(ml(ngl*3,ngl*3),el(ngl*3)) ; ml=1 ; el=0
    i_f = ngl*3 ; j_f=i_f
elseif(cdim==3)then
    allocate(ml(ngl*4,ngl*4),el(ngl*4)) ; ml=1 ; el=0
    i_f = ngl*4 ; j_f=i_f
else
  print*,'Escolha a dimensao correta, 2 ou 3 !'
  stop
end if

ica=0
cnt=0
do l = 1,ngl*nnos

        bnds=0;
	do k =1,nel

	    if( ngl==1 .and. ( cdim==2 .or. cdim==3 ) )then
	      el = matel(k,:)
	    elseif( cdim==2 .and. ngl==2 )then
	      do j=1,3
	         el( 2*j - 1 )= 2*matel(k,j)-1 ; el(2*j)= 2*matel(k,j)
	      end do
	    elseif( cdim==3 .and. ngl==4 )then
	      do j=1,4
	         el( 4*j - 3 )= 4*matel(k,j)-3 ; el(4*j - 2)= 4*matel(k,j)-2
	         el( 4*j - 1 )= 4*matel(k,j)-1 ; el(4*j)= 4*matel(k,j)
	      end do
	    end if
	    
	    do i = 1,i_f
	       do j = 1,j_f   

			if( el(i)==l )then 
		           if(el(j) >= el(i))then

				if(bnds(el(j)).eq.0)then
			                  bndst(el(i)) = bndst(el(i)) + 1
			                  bnds(el(j))  = el(j)

					  cnt=cnt+1
			                  ica(cnt) = el(j)
				else
			                  bnds(el(j))  = el(j)
				end if

		           end if
			end if

	       end do
	    end do

	end do

        allocate(ini(bndst(l)))

        if(l==1)then
          ini = ica( 1:bndst(l) )
        else 
          ini = ica( sum(bndst(1:l-1))+1:sum(bndst(1:l)) )
	end if

	call oerdena( bndst(l) , ini )	

        ica( sum(bndst(1:l-1))+1:sum(bndst(1:l)) ) = ini

        deallocate(ini)

end do

nnzi = sum(bndst) ; print*,'Primeira Estimativa de nnz: ',nnzi ; print*,''
nnz  = nnzi

allocate(ic(nnz),npr(ngl*nnos+1))

npr(1)=1
do i = 2,ngl*nnos
     npr(i) = sum(bndst(1:i-1)) + 1
end do
npr(ngl*nnos+1)=npr(ngl*nnos)+1
do i = 1,nnz
     ic(i) = ica(i)
end do

deallocate(bnds,bndst)

end subroutine

!************************************

subroutine recount_arrays(ngl,nnos,ic,npr,vl,nnz)
implicit none
integer,intent(in)    :: ngl,nnos
integer,intent(inout) :: nnz
integer,allocatable,intent(inout) :: ic(:),npr(:)

!integer,allocatable,intent(inout) :: vl(:)
complex(dpcs),allocatable,intent(inout) :: vl(:)
!real(dpcs),allocatable,intent(inout) :: vl(:)

integer :: i,j,k,cnt
integer,allocatable :: ica(:),npra(:)

!integer,allocatable :: vla(:)
complex(dpcs),allocatable :: vla(:)
!real(dpcs),allocatable :: vla(:) 

allocate(vla(nnz),ica(nnz),npra(ngl*nnos+1)) ; vla=0.d0 ; ica=0;

npra(1)=1
i = 1
nnz=0
do k = 1,ngl*nnos
   cnt=0
   do j = npr(k),npr(k+1)-1

      !if( vl(j) .ne. 0 )then
      if( vl(j) .ne. (0.d0,0.d0) )then
      !if( vl(j) .ne. 0.d0 )then

         nnz=nnz+1
         cnt=cnt+1
         vla(nnz) = vl(j)
         ica(nnz) = ic(j)

      end if

   end do
   i = i + 1; npra(i) = npra(i-1) + cnt;
end do

npra(ngl*nnos+1)=nnz+1

deallocate(vl,ic)
allocate(vl(nnz),ic(nnz))

npr=npra
do i = 1,nnz
 vl(i) = vla(i) ; ic(i) = ica(i)
end do

deallocate(vla,ica,npra)

end subroutine

!************************************

subroutine arrays_pds3(cdim,ngl,nnos,nnzi,ic,el,ml,npr,vl)
implicit none

integer,intent(in) :: cdim,ngl,nnos,nnzi,el(:)
integer,intent(in) :: ic(nnzi),npr(ngl*nnos+1)

!integer,intent(inout) :: vl(nnzi),ml(:,:)
complex(dpcs),intent(inout) :: vl(nnzi),ml(:,:)
!real(dpcs),intent(inout) :: vl(nnzi),ml(:,:)

integer :: i,j,k,i_f,j_f,ii,jj,pos

if(cdim==2)then
    i_f=ngl*3 ; j_f=i_f
elseif(cdim==3)then
    i_f=ngl*4 ; j_f=i_f
else
    print*,'Escolha a dimensao correta !';stop
end if

do i = 1,i_f
    ii = el(i)
       do j = 1,j_f
          jj = el(j)
              if(ii.le.jj)then
                 do k = npr(ii),npr(ii+1)-1               
                    if( jj.eq.ic(k) )then
                      vl(k) = vl(k) + ml(i,j)   
                      exit
                    end if
                 end do 
              end if 
       end do
end do

end subroutine

!************************************

subroutine reorg_arrays(nnzi,nnz,ngl,nnos,ic,npr)
!========
! Rotina que faz a recontagem de valores nao nulos dos vetores e os reorganizam, excuindo os valores nulos restantes.
!========
! Entradas:
! nnzi - Quantidade inicial de valores nao nulos. 
! nnz - Quantidade de valores nao nulos recontada.
! ngl - Numero de graus de liberdade.
! nnos - Numero de nos da malha.
! Saidas:
! ic - Indice de colunas dos valores nao nulos do vetor Vl, com zeros remanescentes.
! npr - vetor que guarda a posicao onde inicia os valores nao nulos de cada linha da matriz global, guardados em Vl.
implicit none

integer,intent(in)    :: nnzi,nnz,ngl,nnos
integer,allocatable,intent(inout) :: ic(:),npr(:) 

!complex(dpcs),intent(inout) :: vl(:)
!real(dpcs),intent(inout) :: vl(:)

integer :: cnt,i,nv
integer,allocatable :: icf(:) 

!integer,allocatable :: vlf(:)
!complex(dpcs),allocatable :: vlf(:)
!real(dpcs),allocatable :: vlf(:)

allocate(icf(nnz))

nv=ngl*nnos

cnt=0
do i = 1,nnzi
   if(ic(i).ne.0)then
     cnt=cnt+1
     icf(cnt)=ic(i)
   end if
end do
do i = 2,nv
     npr(i) = sum(npr(i-1:i))
end do
npr(1)=1
npr(nv+1)=npr(nv)+1

deallocate(ic)
allocate(ic(nnz))

ic = icf;

deallocate(icf)

end subroutine

!************************************

subroutine arrays_pds4( cdim,ie,nelem,ngl,nnos,bnds,el,ml,npr,ic,nnz )
!========
! Rotina responsavel por montar os vetores com os valores nao nulos da matriz global
!========
! Entradas:
! cdim - variavel que informa a rotina a dimensao do problema, se cdim=2 -> 2D  e cdim=3 -> 3D.
! ie   - indice atual do do elemento da malha.
! ngl  - Numero de graus de liberdade ; nelem - Numero de elementos ; nnos - Numero de nos da malha.
! el   - Vetor com os indices de nos do elemento atual; ml - Matriz local. 
! Saidas:
! npr  - vetor que guarda a posicao onde inicia os valores nao nulos de cada linha da matriz global.
! ic   - Indice de colunas dos valores nao nulos do vetor Vl.
! nnz  - Quantidade de valores nao nulos em Vl.
implicit none

integer,intent(in) :: cdim,ie,nelem,ngl,nnos
integer,intent(inout) :: el(:),npr(:),ic(:),nnz,bnds(:)

integer,intent(inout) :: ml(:,:)
!complex(dpcs),intent(inout) :: ml(:,:)
!real(dpcs),intent(inout) :: ml(:,:)

integer :: i_f,j_f,i,j,l
integer :: sm

if(cdim==2)then
    i_f=ngl*3 ; j_f=i_f
elseif(cdim==3)then
    i_f=ngl*4 ; j_f=i_f
else
    print*,'Escolha a dimensao correta !';stop
end if

if(ie==1)then
    npr=0
    npr(1)=1
    nnz = 0
    ic=0
end if

do i = 1,i_f
	do j = 1,j_f
           if( el(j) >= el(i) )then
                 
                  sm = ( sum(bnds(1:el(i)-1)) ) + ( (el(j)-el(i)) + 1 )

                  if( ic( sm ) == 0 .and. ml(i,j) /= 0 )then
                  !if( ic( sm ) == 0 .and. ml(i,j) /= dcmplx(0.d0,0.d0) )then
                  !if( ic( sm ) == 0 .and. ml(i,j) /= 0.d0 )then
                     nnz=nnz+1 
                     npr(el(i)+1) = npr(el(i)+1) + 1
                     ic( sm ) = el(j)
                  else
                     npr(el(i)+1) = npr(el(i)+1)
                  end if

           end if
	end do
end do

end subroutine

!************************************

subroutine solver_pardiso( ngl , nglob , nnz , Vnn , ia , ja , vfonte )
!=======
! Rotina que resolve o sistema linear via metodo direto, utilizando as rotinas do Pardiso.
!======
!Entradas:
!ngl - Graus de liberdade
!nglob - Numero de variaveis do problema
!nnz - Numero de valores nao nulos da matriz global ( somente diagonal superior )
!Vnn - Vetor com os valores nao nulos da matriz global.
!ia  - Vetor que guarda a posicao onde inicia os valores nao nulos de cada linha da matriz global.
!ja  - Indice de colunas dos valores nao nulos do vetor Vl.
!Vfonte - Vetor fonte do problema
!Saidas:
!Vfonte - A solucao do sistema linear e armazenada no mesmo vetor fonte.

!!use mkl_pardiso

implicit none

integer,intent(in)     :: ngl,nglob,nnz
integer,intent(in)     :: ia(ngl*nglob+1) , ja(nnz)

complex(dpcs),intent(in)    :: Vnn(nnz)
complex(dpcs),intent(inout) :: vfonte(ngl*nglob)
!real(dpcs),intent(in)    :: Vnn(nnz)
!real(dpcs),intent(inout) :: vfonte(ngl*nglob)

!!TYPE(MKL_PARDISO_HANDLE), ALLOCATABLE  :: pt(:)

INTEGER maxfct, mnum, mtype, phase, nrhs, error, msglvl
INTEGER error1
INTEGER, ALLOCATABLE :: iparm( : )
INTEGER i, idum(1)

COMPLEX(KIND=dpcs) ddum(1)
!REAL(KIND=dpcs) ddum(1)

complex(dpcs),allocatable :: x(:)
!Real(dpcs),allocatable :: x(:)

Real(dpcs) :: t1,t2

!%%%%%%% Pardiso %%%%%%

allocate(x(ngl*nglob))

! Initialize the solver.

call cpu_time(t1)

!%%%%

nrhs = 1 
maxfct = 1 
mnum = 1

ALLOCATE( iparm ( 64 ) )

iparm  = 0 ! Parameters initialization
error  = 0 ! initialize error flag
msglvl = 0 ! print statistical information
mtype  = 6 ! symmetric, indefinite

!!ALLOCATE ( pt ( 64 ) )

!!call pardisoinit( pt, mtype, iparm )

phase = 12 ! only reordering and symbolic factorization

!!CALL pardiso (pt, maxfct, mnum, mtype, phase, ngl*nglob, vnn, ia, ja, &
!!idum, nrhs, iparm, msglvl, ddum, ddum, error)    
!WRITE(*,*) 'Reordering completed ... '

if (error1 /= 0) then
        write(*,*) 'the following error on release stage was detected in reordering: ', error
        stop
endif

phase = 33 ! only factorization
!!CALL pardiso (pt, maxfct, mnum, mtype, phase, ngl*nglob, vnn, ia, ja, &
!!idum, nrhs, iparm, msglvl, vfonte, x, error)
!WRITE(*,*) 'Solve completed ... '

if (error1 /= 0) then
        write(*,*) 'the following error on release stage was detected in solver: ', error
        stop
endif

!.. Termination and release of memory
phase = -1 ! release internal memory
!!CALL pardiso (pt, maxfct, mnum, mtype, phase, ngl*nglob, ddum, idum, idum, &
!!idum, nrhs, iparm, msglvl, ddum, ddum, error)

if (error1 /= 0) then
        write(*,*) 'the following error on release stage was detected in release memory: ', error
        stop
endif


IF ( ALLOCATED( iparm ) )   DEALLOCATE( iparm )

!!IF ( ALLOCATED( pt ) )   DEALLOCATE( pt )

!%%%%

call cpu_time(t2)

!print*,'Tempo para resolver o sistema linear :',(t2-t1)/60.,'minutos'
!print*,'  '

vfonte = x

!%%%%%%%%%%%%%%%%%%%

deallocate(x)

end subroutine

!************************************

subroutine montagem_Vf( cdim, ngl , nos , nnos , vfl , Vf )
!====
! Rotina responsavel por montar o vetor fonte global
!====
! Entradas:
! cdim - Variavel que informa a rotina a dimensao do problema, se cdim=2 -> 2D  e cdim=3 -> 3D.
! ngl  - Numero de graus de liberdade
! nos  - Vetor com os nos do elemento corrente
! nnos - Numero de nos da malha
! vfl - Vetor fonte local
implicit none
integer,intent(in) :: ngl,nnos,cdim
integer,intent(in) :: nos(:)

complex(dpcs),intent(in) :: vfl(:) 
complex(dpcs),intent(inout) :: Vf(ngl*nnos)
!real(dpcs),intent(in) :: vfl(ngl*3) 
!real(dpcs),intent(inout) :: Vf(ngl*nnos)

integer :: i,j,i_f

if(cdim==2)then
    i_f=ngl*3 ;
elseif(cdim==3)then
    i_f=ngl*4 ;
else
    print*,'Escolha a dimensao correta !';stop
end if

do i = 1,i_f
  Vf( nos(i) ) = Vf( nos(i) ) + vfl(i)
end do

end subroutine

!************************************

subroutine condicao_contorno( cdim , ngl , nnb , n_b , vcc , ia , nnz , V_nn , nnos , Vf )
!====
!Rotina responsavel por aplicar as condicoes de contorno na matriz global e no vetor fonte global
!====
! Entradas:
! cdim - Variavel que informa a rotina a dimensao do problema, se cdim=2 -> 2D  e cdim=3 -> 3D.
! ngl - Numero de graus de liberdade
! nnm - Numero de nos da borda
! n_b - Vetor com os nos da borda
! vcc - Vetor com os valores para os nos da borda
! ia  - Vetor que guarda a posicao onde inicia os valores nao nulos de cada linha da matriz global.
! nnz - Numero de valores nao nulos da matriz global
! nnos - Numero de nos da malha
!Entrada e saida: 
! V_nn - Vetor com os indices nao nulos da matriz global 
! Vf - Vetor fonte global
implicit none
integer,intent(in)  :: ngl,nnz,nnb,nnos,cdim 
integer,intent(in)  :: n_b(nnb)
integer,intent(in)  :: ia( ngl*nnos + 1 ) 

complex(dpcs),intent(inout)  :: V_nn(nnz),Vf(ngl*nnos) 
complex(dpcs),intent(in)     :: vcc(nnb)
!real(dpcs),intent(inout)  :: V_nn(nnz),Vf(ngl*nnos) 
!real(dpcs),intent(in)     :: vcc(nnb)

integer,allocatable       :: nos(:)
complex(dpcs),allocatable :: val_borda(:)
!real(dpcs),allocatable :: val_borda(:)
integer                   :: i,j   
complex(dpcs),parameter   :: big = (1.d20,0.d0) 
!real(dpcs),parameter   :: big = 1.d20 


allocate ( nos(ngl*nnb) , val_borda(ngl*nnb) )

if(ngl.eq.1)then
  nos = n_b
  val_borda = vcc  
elseif( ngl.eq.2 .and. cdim.eq.2 )then
  do i=1,nnb
     nos(2*i-1) = 2*n_b(i)-1 
     nos(2*i)   = 2*n_b(i)
     val_borda(2*i-1) = vcc(i)
     val_borda(2*i)   = vcc(i)
  end do 
elseif( ngl.eq.4 .and. cdim.eq.3 )then
  do i=1,nnb
    nos(4*i-3)=n_b(i)*4-3;nos(4*i-2)=n_b(i)*4-2;
    nos(4*i-1)=n_b(i)*4-1;nos(4*i)=n_b(i)*4
    val_borda(4*i-3)=vcc(i);val_borda(4*i-2)=vcc(i);
    val_borda(4*i-1)=vcc(i);val_borda(4*i)=vcc(i);
  end do
end if

do i = 1,ngl*nnb
     V_nn( ia( nos(i) ) ) = big
     Vf(nos(i)) = big * val_borda(i)
end do

deallocate ( nos , val_borda )

end subroutine

!************************************

subroutine banda(cdim,ngl,nel,mel,bnd)
!====
! Rotina que determina a banda da matriz global
!====
! Entradas
! cdim - Variavel que informa a rotina a dimensao do problema, se cdim=2 -> 2D  e cdim=3 -> 3D.
! ngl - Numero de graus de liberdade
! nel - Numero de elementos da malha
! mel - Matriz de elementos da malha
! Saidas:
! bnd - Banda da matriz global
implicit none
integer,intent(in)  :: cdim,nel,ngl,mel(:,:)
integer,intent(out) :: bnd

integer,allocatable :: el(:),difs(:)
integer :: i,j,miv,mav,cn

allocate(difs(nel))

if(cdim==2)then

 allocate(el(ngl*3))

 do i = 1,nel    
    if(ngl==1)then
      el = mel(i,:)
    else
      cn=0
      do j = 1,3
         el(cn+1)=2*mel(i,j)-1;el(cn+2)=2*mel(i,j)
         cn=cn+2
      end do
    end if
    miv=minval(el) ; mav=maxval(el)
    difs(i)=mav-miv
 end do

 deallocate(el)

elseif(cdim==3)then

 allocate(el(ngl*4))

 do i = 1,nel    
    if(ngl==1)then
      el = mel(i,:)
    else
      cn=0
      do j = 1,4
         el(cn+1)=4*mel(i,j)-3;el(cn+2)=4*mel(i,j)-2
         el(cn+3)=4*mel(i,j)-1;el(cn+4)=4*mel(i,j)
         cn=cn+4
      end do
    end if
    miv=minval(el) ; mav=maxval(el)
    difs(i)=mav-miv
 end do

 deallocate(el)

end if

bnd = maxval(difs)+1

deallocate(difs)

end subroutine

!************************************

subroutine oerdena( n , a )
implicit none
integer,intent(in) :: n
integer,intent(inout) :: a(n)

integer :: aux , i , j 

do i = 1 , n-1
       do j = 1 , n-i 
          if(a(j) > a(j+1))then
           aux = a(j)
           a(j) = a(j+1)
           a(j+1) = aux            
          end if
       end do 
end do

end subroutine

!************************************

end module
