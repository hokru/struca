!*********************************
!* A Tool for Structure Analysis *
!*********************************

program  STRUCA
use parm  ! essential parameter
use logic ! essential logic
!use omp_lib
!use ls_rmsd
implicit none
character(80) basename
real(8) com(3),rot1(3),rot2(3)

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
call version

call eval_options()
!call read_options

! do_rmsd=.true.
!do_primitives=.true.

!call rm_substr(trim(filevec(1)),'.xyz',basename)


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

!process molecules
  call tmolrd(trim(filevec(1)),.false.,mol1%xyz,mol1%iat) 
  call get_dist(mol1) 

  call tmolrd(trim(filevec(2)),.false.,mol2%xyz,mol2%iat) 
  call get_dist(mol2) 

!align + rmsd from fit
  call quaternion_fit(nat,mol1%xyz,mol2%xyz)
  call wrxyz(mol1%iat,nat,mol2%xyz,'rot_'//trim(filevec(2)))

!  call run_drmsd(mol1%xyz,mol1%dist,mol2%xyz,mol1%dist) ! needs re-work

! normal RMSD, must be the same as above or something is wrong!
  call get_rmsd(nat,mol1%xyz,mol2%xyz)

!set nice name
  basename=trim(filevec(2))//'_'//trim(filevec(1))

! internals
  call analyse_primitives(mol1,mol2,trim(basename))

print*,''
print*,'*******************'
print*,'* INERTIA TENSOR  *'
print*,'*******************'
! calc rotational constants
! one needs to move the molecule to the COM first!
 call getCOM(com,nat,mol1%xyz,mol1%iat,.true.)
print*,'MOLECULE 1  : ', trim(filevec(1))
 call getIntertia(nat,mol1%iat,mol1%xyz,.false.,rot1)
print*,'MOLECULE 2  : ', trim(filevec(2))
 call getCOM(com,nat,mol2%xyz,mol2%iat,.true.)
 call getIntertia(nat,mol2%iat,mol2%xyz,.false.,rot2)
print*,''
print*,'Deviation (2-1) :'
rot2=rot2-rot1
! REMINDER: mol1/mol2 are now at the com, but NOT in the principle axis frame
 write(*,'(3x,a,3(F12.5,a))') 'delta_rot:  A= ', rot2(1),' B= ',rot2(2) ,' C= ',rot2(3), ' [Mhz]'
print*,''


endif


end
