! ToDo
! long term goal: re-factor so that we have both the custom internals as well as the automatic constructed ones
! are in the same array. It allows better analysis and printing.


! analyse bonds,angles,dihedral deviations
! 2 task are handled:
! 1) 'all' unique bonds/angles/torsions. detailed output in *.dat files
! 2) user defined bonds/angles/torsions.
subroutine analyse_primitives(mol1,mol2)
use parm
use constant, only: au2ang
use internals
use logic, only: thresh_bond,thresh_ang,thresh_tor
implicit none
type(molecule) :: mol1,mol2
integer bond(nat,nat),kk
real(8) s,db(nat),x
real(8) ang1,ang2,anggrad1,anggrad2
integer ii,jj
integer io,tot
real(8) dummy(1000),thresh
real(8), allocatable :: tvec(:)

dummy=0
! bonding from molecule 1
call bondmatrix(mol1%xyz,mol1%iat,bond)

call read_intcoord('struca.control',mol1%xyz)


print*,'******************'
print*,'* BOND LENGTHS   *'
print*,'******************'

thresh=thresh_bond

open(newunit=io,file='bonds.dat')
k=0
s=0
kk=0
do i=1,nat-1
 do j=i+1,nat
   k=k+1
   if(i==j) cycle
   if(bond(i,j)==1) then
    kk=kk+1
    x=mol1%dist(k)-mol2%dist(k)
    s=s+x
    dummy(kk)=x
   if(abs(x)>=thresh) then
     write(io,'(a,I5,I5,F8.4,a)'),'delta_bond',i,j,x,' * '
   else
     write(io,'(a,I5,I5,F8.4)'),'delta_bond',i,j,x
   endif
  endif
 enddo
enddo

!do i=1,int_nb
! print*, int_bval(i),int_bcast(i,1:2)
!enddo

!tot=int_nb+kk
tot=kk
allocate(tvec(tot))
tvec=0d0
tvec(1:tot)=dummy(1:tot)
print*,  ' # covalent bonds ', kk
print'(a,F8.4)',' mean bond deviation [A]     : ', s/dble(kk)
print'(a,F8.4)',' max       deviation [A]     : ', maxval(tvec)
print*,  ' output: bonds.dat'
close(io)

! currently, we dont have the atom number data... :(
print'(2x,a,F6.2,a)','|deviations| above ',thresh,' A :'
print*,'(see * in bonds.dat)'
do i=1,tot
  if(abs(tvec(i))>=thresh) then
!    print '(1x,a,x,2(I4,x),F8.2)',tvec(i)
    print '(1x,x,F8.3)',tvec(i)
  endif
enddo

!call print_threshold(tot,tvec,0.1d0)

print*,'******************'
print*,'* VALENCE ANGLES *'
print*,'******************'
thresh=thresh_ang
open(newunit=io,file='angles.dat')

s=0
kk=0
do i=1,nat-1
 do j=1,nat
  ii=nat*(i-1)+j
  if(i==j) cycle
  do k=i+1,nat
   jj=nat*(j-1)+l
   if(k==j.or.k==i) cycle
!   print*, i,j,k
   if(bond(i,j)==1.and.bond(j,k)==1) then
    kk=kk+1
    s=s+x
    call  angle(mol1%xyz,i,j,k,ang1,anggrad1)
    call  angle(mol2%xyz,i,j,k,ang2,anggrad2)
    x=anggrad1-anggrad2
    if(abs(x)>=thresh) then
     write(io,'(a,3(I5,x),F8.1,a)'),'delta_angle',i,j,k,x,' * '
    else
     write(io,'(a,3(I5,x),F8.1)'),'delta_angle',i,j,k,x
    endif
   endif
  enddo
 enddo
enddo

print*,  ' # valence angles ', kk
print'(a,F8.1)',' mean angle deviation [deg]  : ', s/dble(kk)
print*,  ' output: angles.dat'
print'(2x,a,F6.2,a)','|deviations| above ',thresh,' deg are marked (*) '
close(io)
print*,'******************'
print*,'* TORSIONS       *'
print*,'******************'


thresh=thresh_tor
open(newunit=io,file='torsions.dat')

s=0
kk=0
! i-l pair loop + full j,k
! ok ???
do i=1,nat-1
  do j=1,nat
    if(i==j) cycle
    do k=1,nat
      ii=nat*(i-1)*k
      if(k==j.or.k==i) cycle
         do l=i+1,nat
           jj=nat*(j-1)*l
           if(l==k.or.l==j.or.l==i) cycle
            if(bond(i,j)==1.and.bond(j,k)==1.and.bond(k,l)==1) then
             kk=kk+1
             call dihed(mol1%xyz,i,j,k,l,ang1,anggrad1)
             call dihed(mol2%xyz,i,j,k,l,ang2,anggrad2)
             call torsionfix(anggrad1,anggrad2,x)
             x=anggrad1-anggrad2
             s=s+x
             if(abs(x)>=thresh) then
              write(io,'(a,4(I5,x),F10.1,a)'),'delta_tors',i,j,k,l,x,' * '
             else
              write(io,'(a,4(I5,x),F10.1)'),'delta_tors',i,j,k,l,x
             endif
            endif
         enddo
      enddo
   enddo
enddo

print*,  ' # torsions ', kk
print'(a,F8.1)',' mean torsion deviation [deg]: ', s/dble(kk)
print*,  ' output: torsions.dat'
print'(2x,a,F6.2,a)','|deviations| above ',thresh,' deg are marked (*) '
close(io)

end subroutine

! fixed 0/360 switches in torsion differences
subroutine torsionfix(t0,t1,dt)
use constant, only: pi
implicit none
logical qI,qIV,shift
real(8) t0,t1,dt
! Check for the case that the torsion goes from the I to the IV quadrant and adjust accordingly
! we check for "greater/less equal" since we might want to reach 0 as target value

!  IV to I
qI=.false.
qIV=.false.
shift=.false.
 if(t1.ge.0.and.t1.le.90) qI=.true.
 if(t0.ge.270.and.t0.le.360) qIV=.true.
 if(qI.and.qIV) shift=.true.

 if(shift) then
  t1=t1+360
  dt=t1-t0
  goto 999
 endif

 qI=.false.
 qIV=.false.

! I to VI
qI=.false.
qIV=.false.
shift=.false.
 if(t0.ge.0.and.t0.le.90) qI=.true.
 if(t1.ge.270.and.t1.le.360) qIV=.true.
 if(qI.and.qIV) shift=.true.

 if(shift) then
  t1=t1-360
  dt=t1-t0
  goto 999
 endif

999 continue

 if(.not.shift)  dt=t1-t0

return
end





! UNFINISHED
! read user defined bonds/angles/torsions
! modified from xopt
subroutine read_intcoord(infile,xyz0)
use constant, only:au2ang
use internals
implicit none
real(8) xyz0(3,*) ! coordinate that will serve as reference
character(80) aa,bb,cc
logical da
integer ii,jj,kk,ll
integer s2i,io
real(8) s2r,dih,dgrad,ang,dang,bl,di360,anggrad,grad2rad
real(8) bdist,dbond

character(*) infile
logical fstr
!real(8) val0()
!integer n_bond(*),n_ang(*),n_dihed(*),ires
integer iint,n_ints

! internal arrays             
!integer int_nb,int_na,int_nt
!real(8), allocatable :: int_bval(:)
!integer, allocatable :: int_bcast(:,:)
!real(8), allocatable :: int_aval(:)
!integer, allocatable :: int_acast(:,:,:)
!real(8), allocatable :: int_tval(:)
!integer, allocatable :: int_tcast(:,:,:,:)

da=.false.
inquire(file=infile,exist=da)
if(.not.da)  return
open(newunit=io,file=infile,status='old')


! determine size
int_nb=0
int_na=0
int_nt=0
do
 read(io,'(a)',end=125) aa
 if(index(aa,'$prim').ne.0) then 
  do
    read(io,'(a)',end=126) aa
    if(index(aa,'#').ne.0) cycle
    if(fstr(aa,'bond')) int_nb=int_nb+1
    if(fstr(aa,'ang'))   int_na=int_na+1
    if(fstr(aa,'dihed').or.fstr(aa,'torsion'))  int_nt=int_nt+1
   enddo
 endif
enddo
125 continue
126 continue

n_ints=int_nb+int_na+int_nt
print*,''
print*, '# custom primitives:', n_ints
print*, '  # bonds    :', int_nb
print*, '  # angles   :', int_na
print*, '  # torsions :', int_nt
print*,''

allocate(int_bval(int_nb),int_bcast(int_nb,2))
allocate(int_aval(int_na),int_acast(int_na,3))
allocate(int_tval(int_nt),int_tcast(int_nt,4))

rewind(io)
iint=0
int_nb=0
int_na=0
int_nt=0
do
 read(io,'(a)',end=124) aa
 if(index(aa,'$prim').ne.0) then 
 print*,''
 print*,' user definied internal coordinates: '
 print*,''
  do
    read(io,'(a)',end=123) aa
    if(index(aa,'#').ne.0) cycle

    if(fstr(aa,'bond')) then 
     iint=iint+1
     int_nb=int_nb+1
     call charXsplit(aa,bb,2)
      ii=s2i(bb)
     call charXsplit(aa,bb,3)
      jj=s2i(bb)
      bdist=dbond(xyz0(1,ii),xyz0(1,jj))
      int_bval(int_nb)=bdist
      int_bcast(int_nb,1)=ii
      int_bcast(int_nb,2)=jj
      write(*,'(2x,a,2I4,2x,a,F7.2)') 'bond [A]      ', ii,jj,' | value: ',bdist
    endif

    if(fstr(aa,'ang')) then 
     iint=iint+1
     int_na=int_na+1
     call charXsplit(aa,bb,2)
      ii=s2i(bb)
     call charXsplit(aa,bb,3)
      jj=s2i(bb)
     call charXsplit(aa,bb,4)
      kk=s2i(bb)
     call angle(xyz0,ii,jj,kk,ang,dang)
      int_aval(int_na)=dang
      int_acast(int_na,1)=ii
      int_acast(int_na,2)=jj
      int_acast(int_na,3)=kk
     write(*,'(2x,a,3I4,2x,a,F7.2)') 'angle [deg]   ', ii,jj,kk,' | value: ',dang
    endif

    if(fstr(aa,'dihed').or.fstr(aa,'torsion')) then 
     iint=iint+1
     int_nt=int_nt+1
     call charXsplit(aa,bb,2)
      ii=s2i(bb)
     call charXsplit(aa,bb,3)
      jj=s2i(bb)
     call charXsplit(aa,bb,4)
      kk=s2i(bb)
     call charXsplit(aa,bb,5)
      ll=s2i(bb)
     call dihed(xyz0,ii,jj,kk,ll,dih,dgrad)
     int_tval(int_nt)=dgrad
     int_tcast(int_nt,1)=ii
     int_tcast(int_nt,1)=jj
     int_tcast(int_nt,1)=kk
     int_tcast(int_nt,1)=ll
     write(*,'(2x,a,4I4,2x,a,F7.2)') 'torsion [deg] ',ii,jj,kk,ll,' | value: ',dgrad
    endif

    if(index(aa,'$').ne.0) then
     backspace(io)
     exit
    endif
  enddo
 endif

enddo
123 continue
124 continue
close(io)

! should not happen
if(iint/=n_ints) stop 'I/O error in custom internals'

print*,''
end subroutine

! print all values about thresh for an input vector
subroutine print_threshold(n,vec,thresh)
implicit none
real(8) vec(*), thresh
integer i,n

do i=1,n
  if(vec(i)>=thresh) print '(1x,F8.2)',vec(i)
enddo
end subroutine
