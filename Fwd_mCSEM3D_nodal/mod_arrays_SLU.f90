module Arrays_SLU

  use var_glob

implicit none

integer,parameter   :: dp = kind(1.d0)

contains

subroutine coluna_ptr2( i, inc, nglob, noele, ngl, nnz )
!==============================================================================================
!
!==============================================================================================
IMPLICIT NONE
integer, intent(in)         :: i, inc, ngl, noele, nglob
integer, intent(out)        :: nnz

integer :: nelem
integer::  j, k, cont, aux(2),ii,id, numb

nelem = nelemg

allocate( col_ptr( ngl*nglob+1 )) ; col_ptr = 0 ;


col_ptr = 100 +i*inc

nnz = sum(col_ptr(:))
aux(1) = col_ptr(1)
col_ptr(1) = 1

do ii=2, ngl*nglob
  aux(2) = col_ptr(ii)
  col_ptr(ii) = aux(1) + 1
  aux(1) = aux(1) + aux(2)
end do

col_ptr(ngl*nglob+1) = aux(1) + 1

end subroutine coluna_ptr2

!=================================================================================================
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!=================================================================================================


subroutine coluna_ptr( nbord, vnbord, nelem, melem, nglob, noele, ngl, nnz )
IMPLICIT NONE
integer, intent(in)         :: ngl,noele,nglob,nelem,nbord
integer, intent(out)        :: nnz
integer,intent(in)          :: vnbord(nbord) , melem(nelem,4) 

!integer :: nelem
integer::  j, k, cont, aux(2),ii,id, numb

!nelem = nelem
!allocate( col_ptr( ngl*nglob+1 )) ; col_ptr = 0 ; ! Alocar esse vetor fora ...

do ii=1, nglob

	numb = 0
	cont = 0
	do j=1, nelem
		do k=1, noele
			if( ii == melem(j,k) )then
				cont = cont + 1
				 exit
			end if
		end do
	end do

	do id = 1, size(vnbord)
		if( ii  == vnbord(id) )then
			do k=1, ngl
				col_ptr(ngl*ii-(k-1)) = (3 + cont)*ngl
			end do
			numb = 1
		end if
	end do

	if( numb == 0 )then
		do k=1, ngl
			col_ptr(ngl*ii-(k-1)) = (2 + cont)*ngl
		end do
	end if

end do

nnz = sum(col_ptr(:))
aux(1) = col_ptr(1)
col_ptr(1) = 1

do ii=2, ngl*nglob
  aux(2) = col_ptr(ii)
  col_ptr(ii) = aux(1) + 1
  aux(1) = aux(1) + aux(2)
end do

col_ptr(ngl*nglob+1) = aux(1) + 1

end subroutine

!======================================================================================
!88888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!======================================================================================

subroutine val_nnz( mloc, noelem, ngl, noele, col_ptr, row_ind, val_nz )
implicit none

integer, intent(in)        :: ngl, noele
integer, dimension(:)  , intent(in)  :: noelem, col_ptr
complex(dp), dimension(:,:), intent(in)  :: mloc

complex(dp), dimension(:), intent(inout) :: val_nz
integer, dimension(:), intent(inout)   :: row_ind

integer:: i, j, k, dim_ij

dim_ij = noele*ngl

do j=1, dim_ij
	do i=1, dim_ij
		do k=col_ptr(noelem(j)), col_ptr(noelem(j)+1)-1
			if( row_ind(k) == noelem(i) )then
				val_nz(k) = val_nz(k) + mloc(i,j)
				exit
			else if( row_ind(k) == 0 )then
				val_nz(k) = mloc(i,j)
				row_ind(k)= noelem(i)
				exit
			else
				cycle
			end if
		end do
	end do
end do

end subroutine

subroutine val_nnz_sim( mloc, noelem, ngl, noele, col_ptr, row_ind, val_nz )
!======================================================================================
!
!======================================================================================
implicit none

    integer, intent(in)   :: ngl, noele
    integer, dimension(:)  , intent(in)  :: noelem, col_ptr
    complex(dp), dimension(:,:), intent(in)  :: mloc

    complex(dp), dimension(:), intent(inout) :: val_nz
    integer, dimension(:), intent(inout)   :: row_ind

    integer:: i, j, k, dim_ij

    dim_ij = noele*ngl

    do j=1, dim_ij
        do i=1, dim_ij
            do k=col_ptr(noelem(j)), col_ptr(noelem(j)+1)-1
                if( row_ind(k) == noelem(i) )then
                    val_nz(k) = val_nz(k) + mloc(i,j)
                    exit
                else if( row_ind(k) == 0 .and. noelem(i) >= noelem(j) )then
                    val_nz(k) = mloc(i,j)
                    row_ind(k)= noelem(i)
                    exit
                else
                    cycle
                end if
            end do
        end do
    end do

end subroutine


!======================================================================================
!88888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!======================================================================================

subroutine eliminando_zeros_cm( nnosg , nnz , ngl )
IMPLICIT NONE
integer, intent(inout)        :: ngl, nnz , nnosg

integer, allocatable, dimension(:) :: row_ind_aux, col_ptr_aux
complex(dp), allocatable, dimension(:)  :: val_nz_aux

integer:: l(2), j, k, nzeros, ng, ncol

nzeros=0 ; l(2)=0 ; l(1)=1 ; ncol = nnosg*ngl
!print*, 'ncol', ncol

ng = 0
do j=1, ncol
	l(1) = col_ptr(j)
	l(2) = col_ptr(j+1)-1
	do k=l(1), l(2)
		if( val_nz(k) /= (0.d0, 0.d0) )then
!			nzeros = nzeros + 1
			ng = ng + 1
		end if
	end do
end do

!ng = nnz-nzeros
!print*, 'ng', ng

allocate( val_nz_aux( ng ), row_ind_aux( ng ), col_ptr_aux(ncol+1) )
val_nz_aux = 0
row_ind_aux = 0
col_ptr_aux = 1

nnz=0 ; nzeros = 0

do j=1, ncol
	l(1) = col_ptr(j)
	l(2) = col_ptr(j+1)-1
	do k=l(1), l(2)
		if( val_nz(k) /= (0.d0, 0.d0) )then
			nnz = nnz + 1
			val_nz_aux(nnz) = val_nz(k)
			row_ind_aux(nnz) = row_ind(k)
			col_ptr_aux(j+1) = col_ptr(j+1) - nzeros
		else
			nzeros = nzeros + 1
			col_ptr_aux(j+1) = col_ptr(j+1) - nzeros
		end if
	end do
end do

deallocate( row_ind, val_nz )

allocate( val_nz( nnz ), row_ind( nnz ) )

val_nz = val_nz_aux
row_ind = row_ind_aux
col_ptr = col_ptr_aux

deallocate( val_nz_aux, row_ind_aux, col_ptr_aux )

end subroutine

!======================================================================================
!88888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!======================================================================================

subroutine eliminando_zeros_cm2( nnz, ngl )
!======================================================================================
!
!======================================================================================
IMPLICIT NONE
integer, intent(inout)        :: ngl, nnz

integer, allocatable, dimension(:) :: row_ind_aux, col_ptr_aux
complex(dp), allocatable, dimension(:)  :: val_nz_aux

integer:: l(2), j, k, nzeros, ng, ncol

nzeros=0 ; l(2)=0 ; l(1)=1 ; ncol = nnosg*ngl
!print*, 'ncol', ncol

ng = 0
do j=1, ncol
	l(1) = col_ptr(j)
	l(2) = col_ptr(j+1)-1
	do k=l(1), l(2)
		if( val_nz(k) /= (0.d0, 0.d0) )then
!			nzeros = nzeros + 1
			ng = ng + 1
		end if
	end do
end do

!ng = nnz-nzeros
!print*, 'ng', ng

allocate( val_nz_aux( ng ), row_ind_aux( ng ), col_ptr_aux(ncol+1) )
val_nz_aux = 0
row_ind_aux = 0
col_ptr_aux = 1

nnz=0 ; nzeros = 0

do j=1, ncol
	l(1) = col_ptr(j)
	l(2) = col_ptr(j+1)-1
	do k=l(1), l(2)
		if( val_nz(k) /= (0.d0, 0.d0) )then
			nnz = nnz + 1
			val_nz_aux(nnz) = val_nz(k)
			row_ind_aux(nnz) = row_ind(k)
			col_ptr_aux(j+1) = col_ptr(j+1) - nzeros
		else
			nzeros = nzeros + 1
			col_ptr_aux(j+1) = col_ptr(j+1) - nzeros
		end if
	end do
end do

deallocate( row_ind, val_nz )

allocate( val_nz( nnz ), row_ind( nnz ) )

val_nz = val_nz_aux
row_ind = row_ind_aux
col_ptr = col_ptr_aux

deallocate( val_nz_aux, row_ind_aux, col_ptr_aux )

end subroutine

!======================================================================================
!88888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!======================================================================================

subroutine cond_front ( n,vnob,nnz,ngl,column,rowind,nnob,values,F )
implicit none
! n = numero de nós globais.
! nnz = número de valores não nulos, ou nulos se forem da diagonal, da MG.
! column = vetor que contém as colunas dos valores não nulos da MG.
! rowind = vetor que guarda a posição de onde se inicia uma nova linha da matriz global em values.
! matcoord = matriz de coordenadas.
! values = valores não nulos da matriz global, ou nulos se forem da diagonal.
! F = Vetor Fonte VF.
integer,parameter :: dpc = kind(1.d0)
integer,intent(in)     :: n,nnz,ngl,nnob
!real(dp),intent(in)    :: matcoord(n)
!complex(dp),intent(inout) :: Values(nnz),F(ngl*n)
!integer,intent(in)     :: column(nnz),rowind(ngl*n+1)

complex(dp),dimension(:),intent(inout) :: Values,F
integer,dimension(:),intent(in)     :: column,rowind

integer :: i,j,k,l
integer,dimension(8),intent(in) :: vnob
complex(dp) :: big = (1.d30,0.d0)

!nnob = 0
!do i = 1,n
!   if(int(matcoord(i,4)).eq.1)then
!      nnob = nnob + 1
!   end if
!end do

    !allocate(vnob(ngl*nnob))

!k = 0
!do i = 1,n
!   if( int(matcoord(i,4)).eq.1 .and. ngl.eq.2 )then
!     k = k + 1
!     vnob( 2*k - 1 ) = 2*int(matcoord(i,1)) - 1
!     vnob( 2*k ) = 2*int(matcoord(i,1))
!   elseif( int(matcoord(i,4)).eq.1 .and. ngl.eq.1 )then
!     k = k + 1
!     vnob( k ) = int(matcoord(i,1))
!   end if
!end do

do i = 1,ngl*nnob
   do j = 1, (size(rowind)-1)
        do k = rowind(j),rowind(j+1)-1
           if( (vnob(i).eq.j).and.(vnob(i).eq.column(k)) )then
              values(k) = big
           end if
        end do
   end do
 F(vnob(i)) = (0.d0,0.d0)
end do

!deallocate(vnob)

end subroutine

!======================================================================================
!88888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!======================================================================================

SUBROUTINE Cond_Fronteira_DH( no_borda, val_nz, row_ind, col_ptr, vfonte, nglob, ngl )

IMPLICIT NONE
 INTEGER, INTENT(IN)        :: ngl, nglob
! REAL(wp), DIMENSION(:,:), INTENT(IN)    :: matcoord
 COMPLEX(dp), DIMENSION(:), INTENT(INOUT)  :: vfonte, val_nz
 INTEGER, DIMENSION(:), INTENT(IN)    :: row_ind, col_ptr, no_borda
! CHARACTER(LEN=20), INTENT(IN)      :: modeling
 INTEGER          :: i, j, k


 do i=1, nglob
      do k=1, ngl
        do j=col_ptr(ngl*no_borda(i)-(k-1)), col_ptr(ngl*no_borda(i)-(k-1)+1)-1
          if( row_ind(j) == ngl*no_borda(i)-(k-1) )then
            val_nz(j) = (1.d30,0.d0)
          end if
          !vfonte(ngl*no_borda(i)-(k-1))=(0.d0,0.d0)
          if(irs.ne.0) vfonte( ngl*no_borda(i)-(k-1) + (irs-1)*ngb )=(0.d0,0.d0)
        end do
      end do
  end do



END SUBROUTINE Cond_Fronteira_DH


subroutine ordenacao( nnode, ngl )

!=================================================================================================
!
! Programa que ordena o vetor coluna comprimida no sentido crescente
!
!
!
! nnode :-------------> Nº de variáveis do sistema
! ngl :---------------> Nº de graus de liberdade do problema
!
!=================================================================================================

IMPLICIT NONE
integer, intent(in)             :: nnode, ngl

integer :: l, k, ntemp, min_ind_tmp
integer :: max_ind_tmp, min_ind, max_ind
complex(dp) :: val_tmp(2)

integer, allocatable, dimension(:):: ind_temp, ind_aux
complex(dp), allocatable, dimension(:):: val_temp
integer    :: ndim

ndim = 1
do l=1, nnode*ngl

ntemp = col_ptr(l+1) - col_ptr(l)
allocate( val_temp(ntemp), ind_temp(ntemp) )
val_temp(1:ntemp) = (0.d0, 0.d0)
ind_temp(1:ntemp) = 0
val_temp(1:ntemp) = val_nz( col_ptr(l):col_ptr(l+1)-1 )
ind_temp(1:ntemp) = row_ind( col_ptr(l):col_ptr(l+1)-1 )

do k=1, int(ntemp/2)

min_ind = minloc( ind_temp(k:ntemp-(k-1)), ndim )+(K-1)

min_ind_tmp = ind_temp(k)
ind_temp(k) = ind_temp(min_ind)
ind_temp(min_ind) = min_ind_tmp

val_tmp(1) = val_temp(k)
val_temp(k) = val_temp(min_ind)
val_temp(min_ind) = val_tmp(1)

max_ind = maxloc( ind_temp(k:ntemp-(k-1)), ndim )+(K-1)

max_ind_tmp = ind_temp(ntemp-(k-1))
ind_temp(ntemp-(k-1)) = ind_temp(max_ind)
ind_temp(max_ind) = max_ind_tmp

val_tmp(2) = val_temp(ntemp-(k-1))
val_temp(ntemp-(k-1)) = val_temp(max_ind)
val_temp(max_ind) = val_tmp(2)
end do
val_nz( col_ptr(l):col_ptr(l+1)-1 ) = val_temp
row_ind( col_ptr(l):col_ptr(l+1)-1 )= ind_temp
deallocate( val_temp, ind_temp )
end do

end subroutine

!======================================================================================
!88888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!======================================================================================

subroutine liberar_memoria
  implicit none

  deallocate(col_ptr,val_nz,row_ind)

endsubroutine liberar_memoria

end module
