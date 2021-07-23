module arrayspds

! Modulo que constroi os vetores no formato de linha suprimida, para serem usados nas rotinas do Pardiso.
! Autor: Anderson A. Piedade. 26/12/2017

integer,parameter :: dpcs=kind(1.d0) 

contains

!************************************

subroutine initial_arrays(cdim,nel,nnos,ngl,nnz,matel,nedge,gEdge,npr,ic)
! Rotina para iniciar os vetores do pardiso

! nedge = quantidade de arestas do modelo
! gEdge = matriz coma numeracao global das arestas
! cdim  = Escolha da dimensao do problema 2 -> 2D ou 3 -> 3D;
! nel   = Numero de elementos da malha de EF
! nnos  = Numero de nos da malha
! ngl   = Nuero de graus de liberdade
! nnz   = Numero de valores nao nulos da matriz global, inicialmente.
! matel = matriz de elementos ( somente os indices dos nos de cada elemento )
! npr   = Vetor com a informacao da quantidade de elementos nao nulos em cada linha da matriz global
! ic    = indice de coluna dos elementos nao nulos de cada linha da matriz global.

implicit none
integer,intent(in)  :: cdim,ngl,nel,nnos,nedge
integer,intent(out) :: nnz
integer,intent(in)  :: matel(:,:),gEdge(:,:)
integer,allocatable,intent(inout) :: npr(:),ic(:)

integer,allocatable :: bnds(:),cnts(:),cntsin(:),vnn(:)
integer :: i,j,k,l,bnd,nnzi,j_f,i_f,dif,qi
integer,allocatable :: ica(:),el(:),ml(:,:)

integer,allocatable :: indcol(:,:),cnt(:),ics(:),icp(:)

integer :: m

nnz=0
allocate(bnds(ngl*nedge)) ; bnds=0
allocate(cnts(ngl*nedge)) ; cnts=0
allocate(cnt(ngl*nedge))  ; cnt=0
!allocate(bnds(ngl*nnos)) ; bnds=0
!allocate(cnts(ngl*nnos)) ; cnts=0
!allocate(cnt(ngl*nnos)) ; cnt=0

if(cdim==2)then
    allocate(ml(ngl*3,ngl*3),el(ngl*3)) ; ml=1 ; el=0
    i_f = ngl*3 ; j_f=i_f
elseif(cdim==3)then
    allocate(ml(ngl*6,ngl*6),el(ngl*6)) ; ml=1 ; el=0
    i_f = ngl*6 ; j_f=i_f
!    allocate(ml(ngl*4,ngl*4),el(ngl*4)) ; ml=1 ; el=0
!    i_f = ngl*4 ; j_f=i_f
else
  print*,'Escolha a dimensao correta, 2 ou 3 !'
  stop
end if

do k =1,nel

    if( ngl==1 .and. ( cdim==2 .or. cdim==3 ) )then
      el = gEdge(k,:)
      !el = matel(k,:)
    elseif( cdim==2 .and. ngl > 1 )then
      do j=1,3
         do m = ngl-1,0,-1
            el( ngl*j - m ) = ngl*gEdge(k,j) - m
            !el( ngl*j - m ) = ngl*matel(k,j) - m
         end do
!         el( 2*j - 1 )= 2*matel(k,j)-1 ; el(2*j)= 2*matel(k,j)
      end do
    elseif( cdim==3 .and. ngl > 1 )then
      do j=1,6
!      do j=1,4
         do m = ngl-1,0,-1
            el( ngl*j - m )= ngl*gEdge(k,j) - m
            !el( ngl*j - m )= ngl*matel(k,j) - m
         end do
!         el( 4*j - 3 )= 4*matel(k,j)-3 ; el(4*j - 2)= 4*matel(k,j)-2
!         el( 4*j - 1 )= 4*matel(k,j)-1 ; el(4*j)= 4*matel(k,j)
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
  allocate(indcol(ngl*nedge,qi))
!  allocate(indcol(ngl*nnos,qi))
elseif( cdim==3 )then
  qi = ngl*3*maxval(cnts/3)
  allocate(indcol(ngl*nedge,qi))
!  allocate(indcol(ngl*nnos,qi))
end if

 cnt=0;
do k =1,nel

    if( ngl==1 .and. ( cdim==2 .or. cdim==3 ) )then
      el = gEdge(k,:)
      !el = matel(k,:)
    elseif( cdim==2 .and. ngl > 1 )then
      do j=1,3
         do m = ngl-1,0,-1
            el( ngl*j - m )= ngl*gEdge(k,j) - m
!            el( ngl*j - m )= ngl*matel(k,j) - m
         end do
!         el( 2*j - 1 )= 2*matel(k,j)-1 ; el(2*j)= 2*matel(k,j)
      end do
    elseif( cdim==3 .and. ngl > 1 )then
      do j=1,6
!      do j=1,4
         do m = ngl-1,0,-1
            el( ngl*j - m )= ngl*gEdge(k,j) - m
!            el( ngl*j - m )= ngl*matel(k,j) - m
         end do
!         el( 4*j - 3 )= 4*matel(k,j)-3 ; el(4*j - 2)= 4*matel(k,j)-2
!         el( 4*j - 1 )= 4*matel(k,j)-1 ; el(4*j)= 4*matel(k,j)
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

nnzi= ngl*nedge
!nnzi= ngl*nnos
do i = 1,ngl*nedge-1
!do i = 1,ngl*nnos-1
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

allocate(npr(ngl*nedge+1))
!allocate(npr(ngl*nnos+1))
allocate(ica(nnzi))

npr(1)=1

nnz = 0
do i = 1,ngl*nedge-1
!do i = 1,ngl*nnos-1

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

nnz=nnz+1 ; ica(nnz)=ngl*nedge;
!nnz=nnz+1 ; ica(nnz)=ngl*nnos;
npr(ngl*nedge + 1)= npr(ngl*nedge) + 1
!npr(ngl*nnos + 1)= npr(ngl*nnos) + 1

allocate(ic(nnz))

ic=ica(1:nnz)

!print*,'Primeira estimativa de nnz:',nnz
!print*,''

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

subroutine recount_arrays(ngl,nnos,nedge,ic,npr,vl,nnz)
! Rotina que verifica se ainda tem valores nulos apos a montagem da matriz global

! nedge = Numero de arestas
! ngl  = Numero de graus de liberdade
! nnos = Numero de nos da malha
! npr   = Vetor com a informacao da quantidade de elementos nao nulos em cada linha da matriz global
! ic    = indice de coluna dos elementos nao nulos de cada linha da matriz global.
! vl   = Valores nao nulos da matriz global
! nnz  = Numero de valores nao nulos da matriz global apos a recontagem

implicit none
integer,intent(in)    :: ngl,nnos,nedge
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

allocate(vla(nnz),ica(nnz),npra(ngl*nedge+1)) ; vla=0.d0 ; ica=0;
!allocate(vla(nnz),ica(nnz),npra(ngl*nnos+1)) ; vla=0.d0 ; ica=0;

npra(1)=1
i = 1
nnz=0
do k = 1,ngl*nedge
!do k = 1,ngl*nnos
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

npra(ngl*nedge+1)=nnz+1
!npra(ngl*nnos+1)=nnz+1

deallocate(vl,ic)
allocate(vl(nnz),ic(nnz))

npr=npra
do i = 1,nnz
 vl(i) = vla(i) ; ic(i) = ica(i)
end do

deallocate(vla,ica,npra)

end subroutine

!************************************

subroutine set_arrays_pds(cdim,ngl,nnos,nedge,nnz,ic,el,ml,npr,vl)

! Rotina que monta o vetor com os valores nao nulos da matriz global

! nedge = Numero de arestas
! cdim  = Escolha da dimensao do problema 2 -> 2D ou 3 -> 3D;
! nel   = Numero de elementos da malha de EF
! nnos  = Numero de nos da malha
! ngl   = Nuero de graus de liberdade
! nnz   = Numero de valores nao nulos da matriz global, inicialmente.
! ic    = indice de coluna dos elementos nao nulos de cada linha da matriz global.
! el    = Vetor com os indices do elemento
! ml    = Matriz local
! npr   = Vetor com a informacao da quantidade de elementos nao nulos em cada linha da matriz global
! vl    = Valores nao nulos da matriz global

implicit none

integer,intent(in) :: cdim,ngl,nnos,nnz,el(:),nedge
integer,intent(in) :: ic(nnz),npr(ngl*nedge+1)
!integer,intent(in) :: ic(nnz),npr(ngl*nnos+1)

!integer,intent(inout) :: vl(nnz),ml(:,:)
complex(dpcs),intent(inout) :: vl(nnz),ml(:,:)
!real(dpcs),intent(inout) :: vl(nnz),ml(:,:)

integer :: i,j,k,i_f,j_f,ii,jj,pos

if(cdim==2)then
    i_f=ngl*3 ; j_f=i_f
elseif(cdim==3)then
    i_f=ngl*6 ; j_f=i_f
!    i_f=ngl*4 ; j_f=i_f
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

subroutine montagem_Vf( cdim, ngl , nos , nnos , nedge , vfl , Vf )
!====
! Rotina responsavel por montar o vetor fonte global
!====
! Entradas:
! nedge - Numero de arestas
! cdim - Variavel que informa a rotina a dimensao do problema, se cdim=2 -> 2D  e cdim=3 -> 3D.
! ngl  - Numero de graus de liberdade
! nos  - Vetor com os nos do elemento corrente
! nnos - Numero de nos da malha
! vfl - Vetor fonte local
implicit none
integer,intent(in) :: ngl,nnos,cdim,nedge
integer,intent(in) :: nos(:)

complex(dpcs),intent(in) :: vfl(:) 
complex(dpcs),intent(inout) :: Vf(ngl*nedge)
!complex(dpcs),intent(inout) :: Vf(ngl*nnos)

!real(dpcs),intent(in) :: vfl(ngl*3) 
!real(dpcs),intent(inout) :: Vf(ngl*nnos)

integer :: i,j,i_f

if(cdim==2)then
    i_f=ngl*3 ;
elseif(cdim==3)then
    i_f=ngl*6 ;
!    i_f=ngl*4 ;
else
    print*,'Escolha a dimensao correta !';stop
end if

do i = 1,i_f
  Vf( nos(i) ) = Vf( nos(i) ) + vfl(i)
end do

end subroutine

!************************************

subroutine condicao_contorno( cdim , ngl , nnb , n_b , vcc , ia , nnz , V_nn , nnos , nedge , Vf )
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
integer,intent(in)  :: ngl,nnz,nnb,nnos,cdim,nedge 
integer,intent(in)  :: n_b(nnb)
integer,intent(in)  :: ia( ngl*nedge + 1 )
!integer,intent(in)  :: ia( ngl*nnos + 1 ) 

complex(dpcs),intent(inout)  :: V_nn(nnz),Vf(ngl*nedge) 
!complex(dpcs),intent(inout)  :: V_nn(nnz),Vf(ngl*nnos) 
!complex(dpcs),intent(inout)  :: V_nn(nnz),Vf(ngl*nnos,3) 

complex(dpcs),intent(in)     :: vcc
!real(dpcs),intent(inout)  :: V_nn(nnz),Vf(ngl*nnos) 
!real(dpcs),intent(in)     :: vcc(nnb)

integer,allocatable       :: nos(:)
complex(dpcs),allocatable :: val_borda(:)
!real(dpcs),allocatable :: val_borda(:)
integer                   :: i,j   
complex(dpcs),parameter   :: big = (1.d20,0.d0) 
!real(dpcs),parameter   :: big = 1.d20 

integer                   :: m 

allocate ( nos(ngl*nnb) , val_borda(ngl*nnb) )

if(ngl.eq.1)then
  nos = n_b
  val_borda = vcc  
elseif( ngl.gt.1 .and. cdim.eq.2 )then
  do i=1,nnb
     do m = ngl-1,0,-1
        nos(ngl*i - m) = ngl*n_b(i) - m        
        val_borda(ngl*i - m) = vcc
     end do
!     nos(2*i-1) = 2*n_b(i)-1 
!     nos(2*i)   = 2*n_b(i)
!     val_borda(2*i-1) = vcc
!     val_borda(2*i)   = vcc
  end do 
elseif( ngl.gt.1 .and. cdim.eq.3 )then
  do i=1,nnb
     do m = ngl-1,0,-1
        nos(ngl*i-m) = ngl*n_b(i) - m
        val_borda(ngl*i-m) = vcc
     end do
!     nos(4*i-3)=n_b(i)*4-3;nos(4*i-2)=n_b(i)*4-2;
!     nos(4*i-1)=n_b(i)*4-1;nos(4*i)=n_b(i)*4
!     val_borda(4*i-3)=vcc;val_borda(4*i-2)=vcc;
!     val_borda(4*i-1)=vcc;val_borda(4*i)=vcc;
  end do
end if

do i = 1,ngl*nnb

     V_nn( ia( nos(i) ) ) = big
     Vf(nos(i)) = big * val_borda(i)

!     Vf(nos(i),1) = big * val_borda(i)
!     Vf(nos(i),2) = big * val_borda(i)
!     Vf(nos(i),3) = big * val_borda(i)
end do

deallocate ( nos , val_borda )

end subroutine

!************************************

subroutine banda(cdim,ngl,nel,mel,gEdge,bnd)
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
integer,intent(in)  :: cdim,nel,ngl,mel(:,:),gEdge(:,:)
integer,intent(out) :: bnd

integer,allocatable :: el(:),difs(:)
integer :: i,j,miv,mav,cn

integer :: m

allocate(difs(nel))

if(cdim==2)then

 allocate(el(ngl*3))

 do i = 1,nel    
    if(ngl==1)then
      el = gEdge(i,:)
!      el = mel(i,:)
    else
      cn=0
      do j = 1,3
         do m = ngl-1,0,-1
            cn = cn + 1
            el(cn)= ngl*gEdge(i,j) - m
!            el(cn)= ngl*mel(i,j) - m
         end do
!         el(cn+1)=2*mel(i,j)-1;el(cn+2)=2*mel(i,j)
!         cn=cn+2
      end do
    end if
    miv=minval(el) ; mav=maxval(el)
    difs(i)=mav-miv
 end do

 deallocate(el)

elseif(cdim==3)then

 allocate(el(ngl*6))
! allocate(el(ngl*4))

 do i = 1,nel    
    if(ngl==1)then
      el = mel(i,:)
    else
      cn=0
      do j = 1,6
!      do j = 1,4
        do m = ngl-1,0,-1
           cn = cn + 1
           el(cn) = ngl*gEdge(i,j) - m
!           el(cn) = ngl*mel(i,j) - m
        end do
      end do
!      do j = 1,4
!         el(cn+1)=4*mel(i,j)-3;el(cn+2)=4*mel(i,j)-2
!         el(cn+3)=4*mel(i,j)-1;el(cn+4)=4*mel(i,j)
!         cn=cn+4
!      end do
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

! Rotina que ordena um vetor

! n = Tamanho do vetor.
! a = vetor a ser ordenado.

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
