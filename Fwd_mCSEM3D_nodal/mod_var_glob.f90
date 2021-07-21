module var_glob

integer,parameter         :: dpc = kind(1.d0)
real(dpc),parameter       :: pi  = 3.141592653589793238462643383279502884197d0
real(dpc),parameter       :: mi0 = (4.d-7)*pi , ep0 = 8.85d-12 , rho_ar = 1.d6
complex,parameter         :: ip  = dcmplx(0.d0,1.d0)

real(dpc),allocatable     :: mat_coord(:,:) 
integer,allocatable       :: mat_elem(:,:) , mat_face(:,:)
integer                   :: nbord
integer,allocatable       :: vnobord(:)  
integer                   :: ngl = 4 , nelemg , nnosg , nfaceg , nnoele = 4   

real(dpc)                 :: Freq , tr_x , tr_y , tr_z
real(dpc),allocatable     :: freqs(:)
integer                   :: nfreqs

integer,allocatable       :: ind_perf(:) , lab_perfs(:)
real(dpc),allocatable     :: pos_tr(:,:)
real(dpc),allocatable     :: coor_perf(:,:)
integer                   :: num_perf, nrcp 
integer,allocatable       :: NTPP(:),NRPP(:),ind_perfpxma(:),ind_perfpxme(:)
integer,allocatable       :: INDRCP(:),INDRCPX(:,:)
real(dpc),allocatable     :: COORCP(:,:),COOTXP(:,:)
integer                   :: NRCT , NFTT

real(dpc)                 :: minper , maxper
real(dpc)                 :: mipc = 0.d0 , mapc = 0.d0

integer,allocatable       :: flagF(:),NDPF(:)

integer                   :: ncm
real(dpc),allocatable     :: z_int(:)
real(dpc),allocatable     :: rhop(:) , hjp(:)

real(dpc),allocatable     :: vprop(:) , Param(:) , rho_elem(:)
integer,allocatable       :: vidprop(:) , vflprop(:)
integer                   :: nprop , nprm
integer                   :: nnpx,nnpy,nnpz

integer                   :: nnz
integer,allocatable       :: col_ptr(:),row_ind(:)
complex(dpc),allocatable  :: val_nz(:),vfonte(:)

integer                   :: nnzi,cdim=3 
integer,allocatable       :: col_ptra(:),row_inda(:)

character(len=255)       :: input_ele,input_face,input_node,input_1d

complex(dpc),allocatable  :: Ext_in(:) , Exp_in(:) , Eyt_in(:) , Ezt_in(:)
complex(dpc),allocatable  :: Hxt_in(:) , Hxp_in(:) , Hyt_in(:) , Hzt_in(:)
complex(dpc),allocatable  :: potAx(:), potAy(:), potAz(:), potPhi(:)

integer                   :: emd 
complex(dpc)              :: zet0
integer                   :: irs , nrs 
integer                   :: ngb

integer                   :: uprm , dprm
real(dpc)                 :: xih,xfh,yih,yfh,zih,zfh

integer                   :: nem1 , Ninh
integer,allocatable       :: Inh(:) , ihet1(:)

integer                   :: iperf , ifrq , itrm

complex(dpc),allocatable  :: jacb(:,:) , jacbT(:,:) , Qp(:,:) , Hp(:,:) , Gp(:)
complex(dpc)              :: sum_jac
complex(dpc),allocatable  :: sum_jaca(:,:)


integer,allocatable       :: Ielpj(:,:) 
integer                   :: ilp,collp

integer                   :: irsj
integer                   :: idjac

real(dpc)                 :: ttap, rad

complex(dpc),allocatable  :: vtjac(:)

complex(dpc)              :: Hpx, Hpy, Hpz
integer                   :: interp

character(len=100)        :: mname

real(dpc),allocatable     :: rhohet(:)
integer                   :: nhet
real(dpc)                 :: mdip

end module
