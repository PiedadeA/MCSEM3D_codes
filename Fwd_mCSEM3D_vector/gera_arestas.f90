program gera_num_arestas
implicit none

integer :: i,j,k,l,m,n,nel,aux,cnt,ntedj,nglob
integer :: ni,nj,outl
integer :: nedj=6, nnoele=4
integer,allocatable :: mel(:,:), edj(:,:), cedj(:,:), tedj(:,:)
integer,allocatable :: qtd(:), gedgs(:,:), conts(:)

integer,allocatable :: ind_edg_glob(:,:),bd(:)

integer :: n1,n2,flbd=-1

real(4) :: auxr
character(len=100) :: mname



call get_command_argument(1,mname)

open(10,file='./malha/'//trim(mname)//'.1.ele',status='old',action='read')
read(10,*)nel
allocate( mel(nel,nnoele), edj(nel,nedj), cedj(nel,12))
do i = 1,nel
   read(10,*)aux,mel(i,1:nnoele)
   cedj(i,1)  = mel(i,1); cedj(i,2)  = mel(i,2)
   cedj(i,3)  = mel(i,1); cedj(i,4)  = mel(i,3)
   cedj(i,5)  = mel(i,1); cedj(i,6)  = mel(i,4)
   cedj(i,7)  = mel(i,2); cedj(i,8)  = mel(i,3)
   cedj(i,9)  = mel(i,4); cedj(i,10) = mel(i,2)
   cedj(i,11) = mel(i,3); cedj(i,12) = mel(i,4)
end do
 close(10)

open(20,file='./malha/'//trim(mname)//'.1.edge',status='old',action='read')
read(20,*)ntedj
allocate(tedj(ntedj,2)) ; allocate(bd(ntedj)) ; bd = 0;
do i = 1,ntedj
   read(20,*)aux,tedj(i,:),bd(i)
end do
 close(20)

open(30,file='./malha/'//trim(mname)//'.1.node',status='old',action='read')
read(30,*)nglob
 close(30)

allocate(qtd(nglob)); qtd = 0;

do i = 1,ntedj
   n1 = tedj(i,1) ; n2 = tedj(i,2)
   if(n1 < n2)then
      qtd(n1) = qtd(n1) + 1
   else
      qtd(n2) = qtd(n2) + 1
   end if
!   qtd(tedj(i,1)) = qtd(tedj(i,1)) + 1
end do

!print*,'Maio Valor:',maxval(qtd)

allocate(gedgs(nglob,maxval(qtd)), conts(nglob))
allocate(ind_edg_glob(nglob,maxval(qtd)))

gedgs=0; conts=0;

do i = 1,ntedj
   n1 = tedj(i,1) ; n2 = tedj(i,2)
   if(n1 < n2)then
      conts(n1) = conts(n1) + 1
      gedgs(n1,conts(n1)) = i
      ind_edg_glob(n1,conts(n1)) = n2
   else
      conts(n2) = conts(n2) + 1
      gedgs(n2,conts(n2)) = i
      ind_edg_glob(n2,conts(n2)) = n1
   end if
!   conts(tedj(i,1)) = conts(tedj(i,1)) + 1
!   gedgs(tedj(i,1),conts(tedj(i,1))) = i
!   ind_edg_glob(tedj(i,1),conts(tedj(i,1))) = tedj(i,2)
end do

!do i = 1,nglob
!   print*,conts(i),i
!   print*,ind_edg_glob(i,:)
!   print*,gedgs(i,:)
!end do

do i = 1,nel
   m = 0
   do j = 1,nedj
      ni = cedj(i,m+1); nj = cedj(i,m+2)
      l = 0
      do k = 1,conts(ni) 
           l = l + 1
           if( nj == ind_edg_glob(ni,k) )then        
              edj(i,j)=gedgs(ni,k)
              exit
           end if
      end do
      if( l == conts(ni) )then 
         aux = ni
         ni = nj
         nj = aux
         do k = 1,conts(ni) 
            if( nj == ind_edg_glob(ni,k) )then        
               edj(i,j)=gedgs(ni,k)
               exit
            end if
         end do
      end if
      m = m + 2
   end do
end do

open(30,file='./malha/'//trim(mname)//'.1.nedge',status='replace',action='write')

write(30,*)nel
do i = 1,nel
   write(30,*)edj(i,:)
end do
!print*,maxval(edj)
 close(30)

!open(40,file='imag_Ex_IE_1hz.in',status='old',action='read')
!do i = 1,50
!   read(40,*)auxr
!   write(4040,*)auxr, 0.0, 1005.0
!end do
! close(40)

end program




































!do i = 1,nglob
!   print*,qtd(i)
!end do

!do i = 1,nel
!   m = 0
!   do k = 1,nedj
!      ni = cedj(i,m+1); nj = cedj(i,m+2)
!      do j = 1,ntedj
!         if( ( (ni==tedj(j,1)).and.(nj==tedj(j,2)) ).or.( (nj==tedj(j,1)).and.(ni==tedj(j,2)) ) )then
!              edj(i,k) = j
!              exit
!           end if
!      end do 
!      m=m+2
!   end do

!end do

!cnt=0
!do i = 1,nel

!   l=0
   
!   print*,'Edjes do elemento',i

!   do j = 1,nedj

!      if(i==1)then

!        cnt=cnt+1
!        edj(i,j) = cnt

!      else

!        ni = cedj(i,l+1) ; nj = cedj(i,l+2)

!        outl=0
!        do n = 1,i-1

!           m=0
!           do k = 1,nedj

!              if( ( (ni==cedj(n,m+1)).and.(nj==cedj(n,m+2)) ).or.( (nj==cedj(n,m+1)).and.(ni==cedj(n,m+2)) ) )then
!                  edj(i,j) = edj(n,k)
!                  outl=1
!                  exit
!              end if
!              m=m+2

!           end do    
!           if(outl.eq.1)exit

!        end do

!        if(outl==0)then
!           cnt=cnt+1
!           edj(i,j) = cnt 
!        end if

!      end if

!      l=l+2;

!   end do

!end do
