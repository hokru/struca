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
real(8) mem
character(1) flag

type(molecule) mol1
type(molecule) mol2
type(trajectory) traj1

integer natPDB

echo=.false.
do_traj=.false.
do_compare=.false.
do_frag=.false.
do_single=.false.

thresh_bond=0.05d0
thresh_ang=1d0
thresh_tor=5d0
thresh_hbr=3.0d0
!thresh_hba=120d0
thresh_hba=35d0

print*, '**********************'
print*, '* structure analysis *'
print*, '*                    *'
print*, '**********************'
call version
print*,''

!nat=natPDB('1a.pdb')
!  allocate(mol1%xyz(3,nat))
!  allocate(mol1%dist(npair))
!  allocate(mol1%iat(nat))
!print*, nat
!call readPDB('1a.pdb',nat,mol1%xyz,mol1%iat)
!stop

call eval_options()

!***********************
!* TRAJECTORY ANALYSIS *
!***********************
if (do_traj) then
 call read_trajxyz(filevec(1),nat,traj1%iat,traj1%mxyz,nmol,.true.)
 mem=(8d0*3.0d0*dble(nat)*dble(nmol))/(1024d0**2)
 print'(a,F8.1,a)', ' memory for just the trajectory : ',mem ,' Mb '
 if(mem>15000) then
   print*,' ** WARNING ** '
   print*,' memory over 15Gb required. Continue [y/n] ? '
   read(*,*) flag
   if(flag/='y') stop
 endif
 allocate(traj1%mxyz(3,nat,nmol),traj1%iat(nat))
 call read_trajxyz(filevec(1),nat,traj1%iat,traj1%mxyz,nmol,.false.)

 ! all data in memory

endif

!*************************
!*  COMPARE 2 MOLECULES  *
!*************************
if(do_compare) then

! poke dimensions
  call tmolrd(trim(filevec(1)),.true.,1,1) 
  npair=(nat*(nat-1))/2
  allocate(mol1%xyz(3,nat),mol2%xyz(3,nat))
  allocate(mol1%dist(npair),mol2%dist(npair))
  allocate(mol1%iat(nat),mol2%iat(nat))

!process molecules
print*,'mol1:'
  call tmolrd(trim(filevec(1)),.false.,mol1%xyz,mol1%iat) 
  call get_dist(mol1) 
print*,'mol2:'
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

! hbonds
call delta_hbonds(mol1,mol2)

! calc rotational constants
 call header_it
! NOTE: one needs to move the molecules to the COM first!
print*,'MOLECULE 1  : ', trim(filevec(1))
 call getCOM(com,nat,mol1%xyz,mol1%iat,.true.)
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


!******************************
!*  SINGLE MOLECULE TREATMENT *
!******************************
if(do_single) then
  call tmolrd(trim(filevec(1)),.true.,1,1)
  npair=(nat*(nat-1))/2
  allocate(mol1%xyz(3,nat))
  allocate(mol1%dist(npair))
  allocate(mol1%iat(nat))
print*,'mol:'
  call tmolrd(trim(filevec(1)),.false.,mol1%xyz,mol1%iat)
  call get_dist(mol1)
  call print_primitives(mol1)
  call hbonds(mol1)

 call header_it
 call getCOM(com,nat,mol1%xyz,mol1%iat,.true.)
 call getIntertia(nat,mol1%iat,mol1%xyz,.false.,rot1)
endif


!***********************
!*  FRAGMENT ANALYSIS  *
!***********************
if(do_frag) then
  call tmolrd(trim(filevec(1)),.true.,1,1)
  npair=(nat*(nat-1))/2
  allocate(mol1%xyz(3,nat))
  allocate(mol1%dist(npair))
  allocate(mol1%iat(nat))
print*,'mol:'
  call tmolrd(trim(filevec(1)),.false.,mol1%xyz,mol1%iat)
  call get_dist(mol1)
! get fragments

! get com distance between fragments

endif

end

subroutine header_it()
print*,''
print*,'*******************'
print*,'* INERTIA TENSOR  *'
print*,'*******************'
print*,''
end subroutine
