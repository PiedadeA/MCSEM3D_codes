module lstm

contains

!------

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

!------

subroutine cgrad_solv( n,nnz,stl,ca,G,F,x )
implicit none
integer,intent(in)  :: n,nnz
integer,intent(in)  :: stl(n+1),ca(nnz)
real(8),intent(in)  :: G(nnz), F(n)
real(8),intent(out) :: x(n)

integer :: nit=1000,i,j,k
real(8),allocatable :: r(:),p(:),Gx(:),Gp(:)
real(8) :: alp,dotr,bta,dotro
real(8),parameter :: tol = 1.d-8
real(8) :: dif

allocate(r(n),Gx(n),p(n),Gp(n))

x = 0.d0;

call prod(n,nnz,stl,ca,G,x,Gx)

r = F - Gx

p = r

do i = 1,nit

   Gp=0.d0
   call prod(n,nnz,stl,ca,G,p,Gp)

   dotr = sum(r*r) ; dotro = dotr
   alp = dotr/sum(p*Gp)

   x = x + alp * p
   r = r - alp * Gp

   dif = dsqrt( sum((r)*r) )/dsqrt(sum((F)*F))

   if( dif < tol )then
!   if( dotr < tol )then
     print*,''
     print*,'End of execution.'
     print*,'erro:',dotr
     print*,'Iterations:',i
     print*,''
     exit
   elseif( i==nit )then
     print*,''
     print*,'End of execution, not converged.'
     print*,'erro:',dotr
     print*,'Iterations:',i
     print*,''
     exit
   end if

   dotr = sum(r*r) 
   bta = dotr/dotro

   p = r + bta * p
   
end do

deallocate(r,Gx,p,Gp)

end subroutine

!------

subroutine prod(n,nnz,stl,ca,vna,b,p)
implicit none
integer,intent(in)  :: n,nnz
integer,intent(in)  :: stl(n+1),ca(nnz)
real(8),intent(in)  :: vna(nnz), b(n)
real(8),intent(out) :: p(n)

integer :: i,j

p=0.d0
do i = 1,n
	do j = stl(i),stl(i+1)-1
		if(j==stl(i))then
			p(i) = p(i) + vna(j) * b(ca(j))
		else
			p(ca(j)) = p(ca(j)) + vna(j) * b(i)
			p(i) = p(i) + vna(j) * b(ca(j))
		end if
	end do
end do

end subroutine

!------

subroutine ccgrad_solv( n,nnz,stl,ca,G,F,x )
implicit none
integer,intent(in)  :: n,nnz
integer,intent(in)  :: stl(n+1),ca(nnz)
complex(8),intent(in)  :: G(nnz), F(n)
complex(8),intent(out) :: x(n)

integer :: nit=1000,i,j,k
complex(8),allocatable :: r(:),p(:),Gx(:),Gp(:)
complex(8) :: alp,dotr,bta,dotro
real(8),parameter :: tol = 1.d-8
real(8) :: dif

allocate(r(n),Gx(n),p(n),Gp(n))

x = 0.d0;

call cprod(n,nnz,stl,ca,G,x,Gx)

r = F - Gx

p = r

do i = 1,nit

   Gp=0.d0
   call cprod(n,nnz,stl,ca,G,p,Gp)

   dotr = sum((r)*r) ; dotro = dotr
   alp = dotr/sum((p)*Gp)

   x = x + alp * p
   r = r - alp * Gp

   dif = cdsqrt( sum(conjg(r)*r) )/cdsqrt(sum(conjg(F)*F))

   if( dif < tol )then
     print*,''
     print*,'End of execution.'
     print*,'erro:',dif
     print*,'Iterations:',i
     print*,''
     exit
   elseif( i==nit )then
     print*,''
     print*,'End of execution, not converged.'
     print*,'erro:',dif
     print*,'Iterations:',i
     print*,''
     exit
   end if

   dotr = sum((r)*r) 
   bta = dotr/dotro

   p = r + bta * p
   
end do

deallocate(r,Gx,p,Gp)

end subroutine

!------

subroutine cprod(n,nnz,stl,ca,vna,b,p)
implicit none
integer,intent(in)  :: n,nnz
integer,intent(in)  :: stl(n+1),ca(nnz)
complex(8),intent(in)  :: vna(nnz), b(n)
complex(8),intent(out) :: p(n)

integer :: i,j

p=0.d0
do i = 1,n
	do j = stl(i),stl(i+1)-1
		if(j==stl(i))then
			p(i) = p(i) + vna(j) * b(ca(j))
		else
			p(ca(j)) = p(ca(j)) + vna(j) * b(i)
			p(i) = p(i) + vna(j) * b(ca(j))
		end if
	end do
end do

end subroutine

!------

subroutine Cholesk2( N,nnz,ia,ja,va,L )
implicit none
integer,intent(in)  :: N,nnz
integer,intent(in)  :: ia(N+1),ja(nnz)
real(8),intent(in)  :: va(nnz)
real(8),intent(out) :: L(nnz)
integer,allocatable :: Mrw(:,:)

integer :: i,j,k,p,q,s
real(8) :: sm,vc,vr
integer :: col,row,q1,q2,ic,ir

integer,allocatable :: ill(:),ial(:),iold(:)
real(8),allocatable :: aa(:),ll(:)
integer             :: cnt,cnnz,ind

integer,allocatable :: qnr(:),c(:)

allocate(ill(nnz),ial(n+1),iold(nnz))
!allocate(aa(nnz),ll(nnz))

allocate(qnr(n),c(n))

print*,'Reescrevendo os indices...'
qnr=0
do i = 1,n
   do j = ia(i),ia(i+1)-1
      qnr(ja(j)) = qnr(ja(j)) + 1
   end do
   if(i==1)then
      ial(i) = 1
   else
      ial(i) = sum(qnr(1:i-1)) + 1
   end if
end do
ial(n+1) = ial(n) + qnr(n)
c=0
do i = 1,n
   do j = ia(i),ia(i+1)-1      
      ind       = ial(ja(j)) + c(ja(j))
      ill(ind)  = i
      iold(ind) = j     
      c(ja(j))  = c(ja(j)) + 1
   end do
end do
print*,'Feito!'

do i = 1,N
   do j = ia(i),ia(i+1)-1
      if(j==ia(i))then       !# Termos da diagonal
         if(ia(i)==1)then
            L(j)=dsqrt((va(j)))
         else
            col = i 
            row = ja(j)
            sm=0.d0
            do k = ial(row),ial(row+1)-1-1
                 sm = sm + L(iold(k)) * L(iold(k))
            end do
            L(j) = dsqrt((va(j) - sm))
         end if
      else                   !# Termos abaixo da diagonal
         if(ia(i)==1)then
             L(j) = va(j)/L(ia(i))
         else
             col = i 
             row = ja(j)
             sm=0.d0
             do p = ial(row),ial(row+1)-1-1
                  do q = ial(col),ial(col+1)-1-1-1
                        if( ill(p)==ill(q) )then
                           sm = sm + L(iold(p)) * L(iold(q))
                        end if
                  end do
             end do
             L(j) = ( 1.d0/L(ia(i)) )*( va(j) - sm )
         end if
      end if
   end do
end do

deallocate(qnr,c)
deallocate(ill,ial,iold)
!deallocate(aa,ll)

end subroutine

!------

subroutine cCholesk2( N,nnz,ia,ja,va,L )
implicit none
integer,intent(in)  :: N,nnz
integer,intent(in)  :: ia(N+1),ja(nnz)
complex(8),intent(in)  :: va(nnz)
complex(8),intent(out) :: L(nnz)
integer,allocatable :: Mrw(:,:)

integer :: i,j,k,p,q,s
complex(8) :: sm,vc,vr
integer :: col,row,q1,q2,ic,ir

integer,allocatable :: ill(:),ial(:),iold(:)
complex(8),allocatable :: aa(:),ll(:)
integer             :: cnt,cnnz,ind

integer,allocatable :: qnr(:),c(:)

allocate(ill(nnz),ial(n+1),iold(nnz))
!allocate(aa(nnz),ll(nnz))

allocate(qnr(n),c(n))

!fazer a quantidade em um vetor e ir organizando... Essa e a solucao...

print*,'Reescrevendo os indices...'

qnr=0
do i = 1,n
   do j = ia(i),ia(i+1)-1
      qnr(ja(j)) = qnr(ja(j)) + 1
   end do
   if(i==1)then
      ial(i) = 1
   else
      ial(i) = sum(qnr(1:i-1)) + 1
   end if
end do
ial(n+1) = ial(n) + qnr(n)
c=0
do i = 1,n
   do j = ia(i),ia(i+1)-1      
      ind       = ial(ja(j)) + c(ja(j))
      ill(ind)  = i
      iold(ind) = j     
      c(ja(j))  = c(ja(j)) + 1
   end do
end do
print*,'Feito!'

do i = 1,N
   do j = ia(i),ia(i+1)-1
      if(j==ia(i))then       !# Termos da diagonal
         if(ia(i)==1)then
            L(j) = cdsqrt((va(j)))
         else
            col = i 
            row = ja(j)
            sm=0.d0
            do k = ial(row),ial(row+1)-1-1
                 sm = sm + L(iold(k)) * (L(iold(k)))
            end do
            L(j) = cdsqrt((va(j) - sm))
         end if
      else                   !# Termos abaixo da diagonal
         if(ia(i)==1)then
             L(j) = va(j)/L(ia(i))
         else
             col = i 
             row = ja(j)
             sm=0.d0
             do p = ial(row),ial(row+1)-1-1
                  do q = ial(col),ial(col+1)-1-1
                        if( ill(p)==ill(q) )then
                           sm = sm + L(iold(p)) * (L(iold(q)))
                        end if
                  end do
             end do
             L(j) = ( 1.d0/L(ia(i)) )*( va(j) - sm )
         end if
      end if
   end do
end do

deallocate(qnr,c)
deallocate(ill,ial,iold)
!deallocate(aa,ll)

end subroutine

!------

subroutine pchlesky_cgrad_solv( n,nnz,stl,ca,G,F,L,x )
implicit none
integer,intent(in)  :: n,nnz
integer,intent(in)  :: stl(n+1),ca(nnz)
real(8),intent(in)  :: G(nnz), F(n) , L(nnz)
real(8),intent(out) :: x(n)

integer :: nit=2000,i,j,k
real(8),allocatable :: r(:),p(:),Gx(:),Gp(:)
real(8) :: alp,dotr,bta,dotro
real(8),parameter :: tol = 1.d-8

real(8),allocatable :: rr(:),aux(:),w(:)
real(8) :: dif

allocate(r(n),p(n),Gp(n),aux(n),w(n),rr(n))
!--
x = 0.d0;
!--
call retrosub_dwn (n,nnz,stl,ca,L,F,aux)
r = aux ; rr=aux
!--
p = r
!--
do i = 1,nit

   call retrosub_up (n,nnz,stl,ca,L,p,w)
   call prod(n,nnz,stl,ca,G,w,aux)
   call retrosub_dwn (n,nnz,stl,ca,L,aux,Gp)

   dotr = sum((r)*r) ; dotro = dotr
   alp = dotr/sum((p)*Gp)

   x = x + alp * p
   r = r - alp * Gp   

   !-- out --
   dif = sqrt( sum(r*(r)) )/ sqrt( sum(rr*(rr)) )
   !--

   if( dif < tol )then
     print*,''
     print*,'End of execution.'
     print*,'erro:',dif
     print*,'Iterations:',i
     print*,''
     exit
   elseif( i==nit )then
     print*,''
     print*,'End of execution, not converged.'
     print*,'erro:',dif
     print*,'Iterations:',i
     print*,''
     exit
   end if

   dotr = sum((r)*r) 
   bta = dotr/dotro

   p = r + bta * p
   
end do

aux = x
call retrosub_up (n,nnz,stl,ca,L,aux,x)

deallocate(r,p,Gp,aux,w,rr)

end subroutine

!------

subroutine matmul_lower( n,nnz,ia,ja,L,x,Lx )
implicit none
integer,intent(in) :: n,nnz
integer,intent(in) :: ia(n+1),ja(nnz)
real(8),intent(in) :: L(nnz),x(n)
real(8),intent(out) :: Lx(n)

integer :: i,j
real(8) :: soma

Lx=0.d0
do i = 1,n
   do j = ia(i),ia(i+1)-1
      Lx(ja(j)) = Lx(ja(j)) + L(j) * x(i)
   end do
end do
end subroutine

!------

subroutine matmul_upper( n,nnz,ia,ja,L,x,Lx )
implicit none
integer,intent(in) :: n,nnz
integer,intent(in) :: ia(n+1),ja(nnz)
real(8),intent(in) :: L(nnz),x(n)
real(8),intent(out) :: Lx(n)

integer :: i,j
real(8) :: soma

Lx=0.d0
do i = 1,n
   soma=0.d0
   do j = ia(i),ia(i+1)-1
      soma = soma + L(j) * x(ja(j))
   end do
   Lx(i) = soma
end do

end subroutine

!------

subroutine cpchlesky_cgrad_solv( n,nnz,stl,ca,G,F,L,x )
implicit none
integer,intent(in)  :: n,nnz
integer,intent(in)  :: stl(n+1),ca(nnz)
complex(8),intent(in)  :: G(nnz), F(n) , L(nnz)
complex(8),intent(out) :: x(n)

integer :: nit=1500,i,j,k
complex(8),allocatable :: r(:),p(:),Gx(:),Gp(:)
complex(8) :: alp,dotr,bta,dotro
real(8),parameter :: tol = 1.d-5
real(8) :: dif

complex(8),allocatable :: rr(:),LTLr(:),aux(:),w(:)

allocate(r(n),p(n),Gp(n),aux(n),w(n),rr(n))
!--
x = 0.d0;
!--
call cretrosub_dwn (n,nnz,stl,ca,L,F,aux)
r = aux ; rr = aux
!--
p = r
!--
do i = 1,nit

   call cretrosub_up (n,nnz,stl,ca,L,p,w)
   call cprod(n,nnz,stl,ca,G,w,aux)
   call cretrosub_dwn (n,nnz,stl,ca,L,aux,Gp)

   dotr = sum((r)*r) ; dotro = dotr
   alp = dotr/sum((p)*Gp)

   x = x + alp * p
   r = r - alp * Gp

   !-- out --
   dif = sqrt( sum(r*conjg(r)) )/ sqrt( sum(rr*conjg(rr)) )
   !--

   if( dif < tol )then
     print*,''
     print*,'End of execution.'
     print*,'erro:',dif
     print*,'Iterations:',i
     print*,''
     exit
   elseif( i==nit )then
     print*,''
     print*,'End of execution, not converged.'
     print*,'erro:',dif
     print*,'Iterations:',i
     print*,''
     exit
   end if

   dotr = sum((r)*r) 
   bta = dotr/dotro

   p = r + bta * p
   
end do

aux = x
call cretrosub_up (n,nnz,stl,ca,L,aux,x)

deallocate(r,p,Gp,aux,w,rr)

end subroutine

!------

subroutine biconj_cholesk_solv( n,nnz,stl,ca,G,F,L,x )
implicit none
integer,intent(in)  :: n,nnz
integer,intent(in)  :: stl(n+1),ca(nnz)
complex(8),intent(in)  :: G(nnz), F(n) , L(nnz)
complex(8),intent(out) :: x(n)

integer :: nit=1000,i,j,k
complex(8),allocatable :: r(:),p(:),Gx(:),Gp(:),rb(:),pb(:)
complex(8) :: alp,dotr,bta,dotro
real(8),parameter :: tol = 1.d-5

complex(8),allocatable :: wb(:),aux(:),w(:),rr(:)
complex(8),allocatable :: Lc(:),Gc(:)
real(8) :: dif

allocate(r(n),p(n),Gp(n),aux(n),w(n),rr(n))
allocate(rb(n),pb(n),wb(n))

allocate(Lc(nnz),Gc(nnz))

!call cCholesk2( n,nnz,stl,ca,conjg(G),Lc )
Lc = conjg(L)
Gc = conjg(G)

!--
x = 0.d0;
!--
call cretrosub_dwn (n,nnz,stl,ca,L,F,aux)
r  = aux ; rr = aux
rb = conjg(r)
!--
p  = r
pb = rb
!--
do i = 1,nit

   call cretrosub_up (n,nnz,stl,ca,L,p,w)
   call cprod(n,nnz,stl,ca,G,w,aux)
   call cretrosub_dwn (n,nnz,stl,ca,L,aux,w) ! W = Gp

   call cretrosub_up (n,nnz,stl,ca,Lc,pb,wb)
   call cprod(n,nnz,stl,ca,Gc,wb,aux)
   call cretrosub_dwn (n,nnz,stl,ca,Lc,aux,wb) ! W = G*p

   dotr = sum(conjg(rb)*r) ; dotro = dotr
   alp = dotr/sum(conjg(pb)*w)

   x = x + alp * p
   r = r - alp * w
   rb = rb - conjg(alp) * wb

   !-- out --
   dif = sqrt( sum(r*conjg(r)) )/ sqrt( sum(rr*conjg(rr)) )
   !--

   if( dif < tol )then
     print*,''
     print*,'End of execution.'
     print*,'erro:',dif
     print*,'Iterations:',i
     print*,''
     exit
   elseif( i==nit )then
     print*,''
     print*,'End of execution, not converged.'
     print*,'erro:',dif
     print*,'Iterations:',i
     print*,''
     exit
   end if

   dotr = sum(conjg(rb)*r) 
   bta = dotr/dotro

   p = r + bta * p
   pb = rb + conjg(bta) * pb
   
end do

aux = x
call cretrosub_up (n,nnz,stl,ca,L,aux,x)

deallocate(Lc,Gc)
deallocate(rb,pb)
deallocate(r,p,Gp,aux,w,wb,rr)

end subroutine

!--------

subroutine cmatmul_lower( n,nnz,ia,ja,L,x,Lx )
implicit none
integer,intent(in) :: n,nnz
integer,intent(in) :: ia(n+1),ja(nnz)
complex(8),intent(in) :: L(nnz),x(n)
complex(8),intent(out) :: Lx(n)

integer :: i,j
complex(8) :: soma

Lx=0.d0
do i = 1,n
   do j = ia(i),ia(i+1)-1
      Lx(ja(j)) = Lx(ja(j)) + L(j) * x(i)
   end do
end do
end subroutine

!------

subroutine cmatmul_upper( n,nnz,ia,ja,L,x,Lx )
implicit none
integer,intent(in) :: n,nnz
integer,intent(in) :: ia(n+1),ja(nnz)
complex(8),intent(in) :: L(nnz),x(n)
complex(8),intent(out) :: Lx(n)

integer :: i,j
complex(8) :: soma

Lx=0.d0
do i = 1,n
   soma=0.d0
   do j = ia(i),ia(i+1)-1
      soma = soma + L(j) * x(ja(j))
   end do
   Lx(i) = soma
end do

end subroutine

!------

subroutine retrosub_dwn (n,nnz,ia,ja,L,y,x)
implicit none
integer,intent(in) :: n,nnz
integer,intent(in) :: ia(n+1),ja(nnz)
real(8),intent(in) :: L(nnz),y(n)
real(8),intent(out) :: x(n)

real(8),allocatable :: sm(:)
integer :: i,j

allocate(sm(n))

sm = 0.d0
x  = 0.d0
do i = 1,n
   do j = ia(i),ia(i+1)-1
      if(j==ia(i))then
         sm(ja(j)) = sm(ja(j)) + 0.d0
         x(ja(j))  = ( y(ja(j)) - sm(ja(j)) ) / L(ia(i))
      else
         sm(ja(j)) = sm(ja(j)) + L(j) * x(i)
         x(ja(j)) = ( y(ja(j)) - sm(ja(j)) ) / L(ia(i))
      end if
   end do
end do 

deallocate(sm)

end subroutine

!------

subroutine cretrosub_dwn (n,nnz,ia,ja,L,y,x)
implicit none
integer,intent(in) :: n,nnz
integer,intent(in) :: ia(n+1),ja(nnz)
complex(8),intent(in) :: L(nnz),y(n)
complex(8),intent(out) :: x(n)

complex(8),allocatable :: sm(:)
integer :: i,j

allocate(sm(n))

sm = 0.d0
x  = 0.d0
do i = 1,n
   do j = ia(i),ia(i+1)-1
      if(j==ia(i))then
         sm(ja(j)) = sm(ja(j)) + 0.d0
         x(ja(j))  = ( y(ja(j)) - sm(ja(j)) ) / L(ia(i))
      else
         sm(ja(j)) = sm(ja(j)) + L(j) * x(i)
         x(ja(j)) = ( y(ja(j)) - sm(ja(j)) ) / L(ia(i))
      end if
   end do
end do 

deallocate(sm)

end subroutine

!------

subroutine retrosub_up (n,nnz,ia,ja,L,y,x)
implicit none
integer,intent(in) :: n,nnz
integer,intent(in) :: ia(n+1),ja(nnz)
real(8),intent(in) :: L(nnz),y(n)
real(8),intent(out) :: x(n)

real(8) :: sm
integer :: i,j

x  = 0.d0
do i = n,1,-1
      sm = 0.d0
      do j = ia(i)+1,ia(i+1)-1
        sm = sm + L(j)*x(ja(j))
      end do
      x(i) = ( y(i) - sm ) / L(ia(i))
end do 

end subroutine

!------

subroutine cretrosub_up (n,nnz,ia,ja,L,y,x)
implicit none
integer,intent(in) :: n,nnz
integer,intent(in) :: ia(n+1),ja(nnz)
complex(8),intent(in) :: L(nnz),y(n)
complex(8),intent(out) :: x(n)

complex(8) :: sm
integer :: i,j

x  = 0.d0
do i = n,1,-1
      sm = 0.d0
      do j = ia(i)+1,ia(i+1)-1
        sm = sm + L(j)*x(ja(j))
      end do
      x(i) = ( y(i) - sm ) / L(ia(i))
end do 

end subroutine

!-----

end module
