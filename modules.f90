module parm

 ! holds data for a molecule
 type molecule
  integer(4), allocatable :: iat(:) ! atom type (H=1, C=6, etc)
  real(8), allocatable :: xyz(:,:) ! cartesian coordinates (originally in au)
  real(8), allocatable :: dist(:) ! distance between all atoms
 end type molecule


 type trajectory
  integer(4), allocatable :: iat(:) ! atom type (H=1, C=6, etc)
  real(8), allocatable :: mxyz(:,:,:) ! cartesian coordinates (originally in au) per step (3,nat,nmol)
 end type 

 ! global variables
  integer(4) nat   ! number of atoms
  integer(4) nmol  ! number of molecule in traj
  integer(4) npair ! number of atom paris
  integer i,j,k,l  ! reserved integers
  integer a,b,c,d  ! reserved integers
  character(200) string

! for dummy arrays
integer, parameter:: max_dummy=50000


end module

module logic
integer, parameter :: max_files=10
logical,save:: echo 
logical, save:: do_compare
logical, save:: do_traj
logical, save:: do_single
character(200) filevec (max_files) ! file name vector

! primitive analysis
real(8) thresh_bond,thresh_ang,thresh_tor


!hbonds
real(8) thresh_hbr
real(8) thresh_hba
end module logic

module internals

! internal arrays             
integer int_nb,int_na,int_nt
real(8), allocatable :: int_bval(:)
integer, allocatable :: int_bcast(:,:)
real(8), allocatable :: int_aval(:)
integer, allocatable :: int_acast(:,:)
real(8), allocatable :: int_tval(:)
integer, allocatable :: int_tcast(:,:)



real(8), allocatable :: hb_val(:)
integer, allocatable :: hb_cast(:,:)

end module internals


module constant
! in part taken from psi4
real(8), parameter:: pi = 3.141592653589793d0
real(8), parameter:: au2ang = 0.52917720859d0
real(8), parameter:: amu2au=1.66053886E-27/9.10938215E-31
real(8), parameter:: au2cm =219474.63067d0
real(8), parameter:: au2kcal = 627.5095d0
real(8), parameter:: au2cm1 = 219474.6d0
real(8), parameter:: au2ev= 27.21138d0
real(8), parameter:: au2mhz=6.579684E9

real(8), parameter:: kb_J=1.3806504E-23 ! J/K
real(8), parameter:: kB_au=3.1668114d-6 ! Eh/K
real(8), parameter:: au2fs=0.02418884326505d0


real(8), parameter :: Planck= 6.62606896E-34 !   The Planck constant (Js) 
real(8), parameter :: AvoN= 6.02214179D23    !Avagadro's number 
real(8), parameter :: pc_c= 2.99792458E8      ! Speed of light (ms$^{-1}$) 
real(8), parameter :: bohr2m= 0.52917720859E-10 !   Bohr to meters conversion factor 
real(8), parameter :: amu2kg= 1.660538782E-27  !  Atomic mass units to kg conversion factor 

end module

module atomdata
real(8) ams(118)
real(8) rvdw(94),rcov(94)
character(2) el(95)


data EL/'H ','HE',                                             &
  'LI','BE','B ','C ','N ','O ','F ','NE',                       &
  'NA','MG','AL','SI','P ','S ','CL','AR',                       &
  'K ','CA','SC','TI','V ','CR','MN','FE','CO','NI','CU',        &
  'ZN','GA','GE','AS','SE','BR','KR',                            &
  'RB','SR','Y ','ZR','NB','MO','TC','RU','RH','PD','AG',        &
  'CD','IN','SN','SB','TE','I ','XE',                            &
  'CS','BA','LA','CE','PR','ND','PM','SM','EU','GD','TB','DY',   &
  'HO','ER','TM','YB','LU','HF','TA','W ','RE','OS','IR','PT',   &
  'AU','HG','TL','PB','BI','PO','AT','RN',                       &
  'FR','RA','AC','TH','PA','U ','NP','PU','XX'/ ! ELEM(95) BEING A DUMMY


! taken from PSI4
data ams /1.00782503207,4.00260325415,7.016004548,9.012182201,11.009305406, &
12,14.00307400478,15.99491461956,18.998403224,19.99244017542, &
22.98976928087,23.985041699,26.981538627,27.97692653246,30.973761629, &
31.972070999,34.968852682,39.96238312251,38.963706679,39.962590983, &
44.955911909,47.947946281,50.943959507,51.940507472,54.938045141,  &
55.934937475,58.933195048,57.935342907,62.929597474,63.929142222,  &
68.925573587,73.921177767,74.921596478,79.916521271,78.918337087,  &
85.910610729,84.911789737,87.905612124,88.905848295,89.904704416,  &
92.906378058,97.905408169,98.906254747,101.904349312,102.905504292 , &
105.903485715,106.90509682,113.90335854,114.903878484,119.902194676, &
120.903815686,129.906224399,126.904472681,131.904153457,132.905451932, &
137.905247237,138.906353267,139.905438706,140.907652769,141.907723297, &
144.912749023,151.919732425,152.921230339,157.924103912,158.925346757, &
163.929174751,164.93032207,165.930293061,168.93421325,173.938862089, &
174.940771819,179.946549953,180.947995763,183.950931188,186.955753109, &
191.96148069,192.96292643,194.964791134,196.966568662,201.970643011, &
204.974427541,207.976652071,208.980398734,208.982430435,210.987496271, &
222.017577738,222.01755173,228.031070292,227.027752127,232.038055325, &
231.03588399,238.050788247,237.048173444,242.058742611,243.06138108, &
247.07035354,247.07030708,251.079586788,252.082978512,257.095104724, &
258.098431319,255.093241131,260.105504,263.112547,255.107398,259.114500, &
262.122892,263.128558,265.136151,281.162061,272.153615,283.171792,283.176451, &
285.183698,287.191186,292.199786,291.206564,293.214670/


! atomic radii from Mantina, Valero, Cramer, Truhlar "Atomic radii of elements"
! Copyed by hand,  may contain typos...
! ANGSTROM
!            H       He
data rvdw /1.10d0,1.40d0, &
    ! Li     Be     B     C       N      O     F      Ne
    1.82d0,1.53d0,1.92d0,1.70d0,1.54d0,1.52d0,1.47d0,1.54d0, &
    ! Na    Mg     Al     Si     P      S       Cl     Ar
    2.27d0,1.73d0,1.84d0,2.10d0,1.80d0,1.80d0,1.75d0,1.88d0, &
    ! K      Ca     Sc     Ti     V      Cr      Mn     Fe     Co    Ni     Cu     Zn     Ga     Ge     As     Se     Br    Kr
    2.75d0,2.31d0,2.15d0,2.11d0,2.07d0,2.06d0,2.05d0,2.04d0,2.00d0,1.97d0,1.96d0,2.01d0,1.87d0,2.11d0,1.85d0,1.90d0,1.85d0,2.02d0, &
    !  Rb    Sr     Y      Zr      Nb     Mo    Tc     Ru     Rh     Pd     Ag     Cd     In    Sn      Sb      Te     I     Xe
    3.03d0,2.49d0,2.26d0,2.23d0,2.18d0,2.17d0,2.16d0,2.13d0,2.10d0,2.10d0,2.11d0,2.18d0,1.93d0,2.17d0,2.06d0,2.06d0,1.98d0,2.16d0, &
    ! Cs Ba
    3.32d0,2.68d0, &
    ! La-Lu
    2.43d0,2.42d0,2.40d0,2.46d0,2.38d0,2.36d0,2.35d0,2.34d0,2.33d0,2.31d0,2.30d0,2.29d0,2.27d0,2.26d0,2.24d0, &
    ! Hf     Ta     W       Re     Os    Ir     Pt     Au     Hg     Ti     Pb     Bi     Po     At     Rn
    2.23d0,2.22d0,2.18d0,2.16d0,2.16d0,2.13d0,2.13d0,2.23d0,2.23d0,2.11d0,2.02d0,2.07d0,1.97d0,2.02d0,2.20d0, &
    ! Fr-Pu
    3.48d0,2.83d0,2.47d0,2.45d0,2.43d0,2.41d0,2.39d0,2.43d0/

data rcov /0.32d0,0.37d0, &
    ! Li     Be     B     C       N      O     F      Ne
    1.30d0,0.99d0,0.84d0,0.75d0,0.71d0,0.64d0,0.60d0,0.62d0,  &
    ! Na    Mg     Al     Si     P      S       Cl     Ar
    1.60d0,1.40d0,1.24d0,1.14d0,1.09d0,1.04d0,1.00d0,1.01d0,  &
    ! K      Ca     Sc     Ti     V      Cr      Mn     Fe     Co    Ni     Cu     Zn     Ga     Ge     As     Se     Br    Kr
    2.00d0,1.74d0,1.59d0,1.48d0,1.44d0,1.30d0,1.29d0,1.24d0,1.18d0,1.17d0,1.22d0,1.20d0,1.23d0,1.20d0,1.20d0,1.18d0,1.17d0,1.24d0, &
    !  Rb    Sr     Y      Zr      Nb     Mo    Tc     Ru     Rh     Pd     Ag     Cd     In    Sn      Sb      Te     I     Xe
    2.15d0,1.90d0,1.78d0,1.64d0,1.56d0,1.46d0,1.38d0,1.36d0,1.34d0,1.30d0,1.36d0,1.40d0,1.42d0,1.40d0,1.40d0,1.37d0,1.32d0,1.36d0, &
    ! Cs Ba
    2.38d0,2.06d0,  &
    ! La-Lu
     1.94d0,1.84d0,1.90d0,1.73d0,1.86d0,1.85d0,1.83d0,1.82d0,1.81d0,1.80d0,1.79d0,1.77d0,1.77d0,1.78d0,1.74d0,  &
    ! Hf     Ta     W       Re     Os    Ir     Pt     Au     Hg     Ti     Pb     Bi     Po     At     Rn
    1.64d0,1.58d0,1.50d0,1.41d0,1.36d0,1.32d0,1.30d0,1.64d0,1.88d0,1.48d0,1.45d0,1.50d0,1.42d0,1.47d0,1.46d0,  &
    ! Fr-Pu
    2.42d0,2.11d0,2.01d0,1.90d0,1.84d0,1.83d0,1.80d0,1.80d0/



end module

