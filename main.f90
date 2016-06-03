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
character(120) infile,infile2
real(8) t0,w0,t1,w1

type(molecule) mol1
type(molecule) mol2
type(trajectory) traj1

call cpu_time(t0)
w0=omp_get_wtime ()


echo=.false.
do_rmsd=.false.
do_primitives=.false.
fstr=.false.

print*, '**********************'
print*, '* structure analysis *'
print*, '*                    *'
print*, '**********************'

!call read_options

! do_rmsd=.true.
!do_primitives=.true.


!call getarg(1,infile)
!call read_trajxyz(infile,nat,traj1%iat,traj1%mxyz,nmol,.true.)
!allocate(traj1%mxyz(3,nat,nmol),traj1%iat(nat))
!call read_trajxyz(infile,nat,traj1%iat,traj1%mxyz,nmol,.false.)
!print*,traj1%mxyz
!stop 'halt'

! DO TASKS
! this NEEDS to be simplified
!if(do_rmsd) then

  call getarg(1,infile)
  call getarg(2,infile2)
  call tmolrd(trim(infile),.true.,1,1) 

  npair=(nat*(nat-1))/2
  allocate(mol1%xyz(3,nat),mol2%xyz(3,nat))
  allocate(mol1%dist(npair),mol2%dist(npair))
  allocate(mol1%iat(nat),mol2%iat(nat))

  call tmolrd(trim(infile),.false.,mol1%xyz,mol1%iat) 
  call get_dist(mol1) 
  call tmolrd(trim(infile2),.false.,mol2%xyz,mol2%iat) 
  call get_dist(mol2) 
  call quaternion_fit(nat,mol1%xyz,mol2%xyz)
  call wrxyz(mol1%iat,nat,mol2%xyz,trim(infile2)//'.rot.xyz')
!  call analyse_primitives(mol1,mol2)
!  call run_drmsd(mol1%xyz,mol1%dist,mol2%xyz,mol1%dist)
  call get_rmsd(nat,mol1%xyz,mol2%xyz)
!endif
!call getarg(3,aa)

 if(do_primitives) then
  call getarg(1,infile)
  call getarg(2,infile2)
  print*, 'TASK:  internal coords:  ',trim(infile),' +  ', trim(infile2)
  call tmolrd(trim(infile),.true.,1,1) 

  npair=(nat*(nat-1))/2
  allocate(mol1%xyz(3,nat),mol2%xyz(3,nat))
  allocate(mol1%dist(npair),mol2%dist(npair))
  allocate(mol1%iat(nat),mol2%iat(nat))
  call tmolrd(trim(infile),.false.,mol1%xyz,mol1%iat) 
  call get_dist(mol1) 
  call tmolrd(trim(infile2),.false.,mol2%xyz,mol2%iat) 
  call get_dist(mol2) 
  call analyse_primitives(mol1,mol2)
 endif

!$omp parallel
call omp_set_dynamic( .true. )
ntr=omp_get_num_threads()
!$omp end parallel






end
