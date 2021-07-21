module vglob

integer,parameter    :: db = kind(1.d0)
real(db),parameter   :: pi  = 3.141592653589793238462643383279502884197d0
real(db),parameter   :: mi0 = (4.d-7)*pi , ep0 = 8.85d-12
complex,parameter    :: ip  = dcmplx(0.d0,1.d0)

real(db),allocatable :: mat_coord(:,:) 
integer,allocatable  :: mat_elem(:,:) , mat_face(:,:)
integer              :: nbord
integer,allocatable  :: vnobord(:)  
integer              :: ngl = 1, nelemg , nnosg , nfaceg , nnoele = 4, nedgs = 6, nedgst, lbnd = -33
integer              :: num_perf

integer,allocatable  :: NRPP(:), INDRCP(:), NTPP(:)
real(db),allocatable :: COORCP(:,:), COOTXP(:,:),atx(:)

integer              :: nnpx,nnpy,nnpz
real(db)             :: xih,xfh
real(db)             :: yih,yfh
real(db)             :: zih,zfh

integer              :: ncm
real(db),allocatable :: rhop(:), hjp(:) , z_int(:)

integer              :: nprop, nprm
integer,allocatable  :: vflprop(:)
real(db),allocatable :: vprop(:)

real(db)             :: mind, maxd, minr, maxr

integer              :: nfreqs
real(db),allocatable :: freqs(:)
integer,allocatable  :: vidprop(:)
real(db),allocatable :: rho_elem(:)

integer              :: lendt = 8

real(db)             :: tr_x, tr_y, tr_z, ttap, Freq

integer,allocatable  :: flagf(:), NDPF(:)

integer              :: cdim = 3
integer              :: nnz, nnzi
integer,allocatable  :: col_ptr(:),row_ind(:),col_ptra(:),row_inda(:)
complex(db),allocatable  :: val_nz(:),vfonte(:)

integer              :: ifrq, itrm, iperf

integer              :: nem1 , Ninh
integer,allocatable  :: Inh(:) , ihet1(:)

real(db),allocatable     :: coor_perf(:,:)

integer,allocatable      :: ind_perf(:)
integer                  :: ngb, nrcp, nrs

integer                  :: dprm, uprm

real(db)                 :: rho_ar = 1.d6

integer,allocatable      :: gEdge(:,:),node_edge(:,:)

integer,allocatable      :: inds_Rx_edg(:,:)

integer,allocatable      :: ind_edg(:)  

integer,allocatable      :: fbd(:)

integer                  :: idjac

!-- Variaveis para a Jacobiana --

complex(db),allocatable  :: jacb(:,:), jacbT(:,:), Qp(:,:), Hp(:,:), Gp(:)
complex(db)              :: sum_jac
complex(db),allocatable  :: sum_jaca(:,:)
integer,allocatable      :: Ielpj(:,:) 
integer                  :: ilp,collp

!--

real(db)                 :: rad
complex(db)              :: cpHx, cpHy, cpHz
integer                  :: interp

!--

integer,allocatable      :: nrpft(:)
real(db),allocatable     :: inrpft(:), ofst(:)

!--

complex(db),allocatable  :: Eo(:)

real(db) :: pcmin, pcmax

real(db) :: mdip

complex(db),allocatable :: Exyz(:,:), Hxyz(:,:)
real(db),allocatable    :: rhohet(:)
character(len=100)      :: mname
integer                 :: nhet

!--

end module
