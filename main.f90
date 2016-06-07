!*********************************
!* A Tool for Structure Analysis *
!*********************************

program  STRUCA
use parm  ! essential parameter
use logic ! essential logic
use omp_lib
!use ls_rmsd
implicit none
integer ntr

type(molecule) mol1
type(molecule) mol2
type(trajectory) traj1



echo=.false.
do_traj=.false.
do_compare=.false.
thresh_bond=0.05d0
thresh_ang=1d0
thresh_tor=5d0


print*, '**********************'
print*, '* structure analysis *'
print*, '*                    *'
print*, '**********************'


call eval_options()
!call read_options

! do_rmsd=.true.
!do_primitives=.true.



if (do_traj) then
!call getarg(1,infile)
!call read_trajxyz(filevec(1),nat,traj1%iat,traj1%mxyz,nmol,.true.)
!allocate(traj1%mxyz(3,nat,nmol),traj1%iat(nat))
!call read_trajxyz(filevec(1),nat,traj1%iat,traj1%mxyz,nmol,.false.)
!print*,traj1%mxyz
!stop 'halt'
endif


if(do_compare) then
!  call getarg(filevec(1),infile)
!  call getarg(filevec(2),infile2)
  call tmolrd(trim(filevec(1)),.true.,1,1) 

  npair=(nat*(nat-1))/2
  allocate(mol1%xyz(3,nat),mol2%xyz(3,nat))
  allocate(mol1%dist(npair),mol2%dist(npair))
  allocate(mol1%iat(nat),mol2%iat(nat))

  call tmolrd(trim(filevec(1)),.false.,mol1%xyz,mol1%iat) 

! read custom internals
!  call read_intcoord('struca.control',mol1%xyz)

  call get_dist(mol1) 
  call tmolrd(trim(filevec(2)),.false.,mol2%xyz,mol2%iat) 
  call get_dist(mol2) 
  call quaternion_fit(nat,mol1%xyz,mol2%xyz)
  call wrxyz(mol1%iat,nat,mol2%xyz,trim(filevec(2))//'.rot.xyz')
!  call run_drmsd(mol1%xyz,mol1%dist,mol2%xyz,mol1%dist) ! needs re-work
  call get_rmsd(nat,mol1%xyz,mol2%xyz)
  call analyse_primitives(mol1,mol2)
endif


end
