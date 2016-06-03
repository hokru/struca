! analyse bonds,angles,dihedral deviations
! 2 task are handled:
! 1) 'all' unique bonds/angles/torsions. detailed output in *.dat files
! 2) user defined bonds/angles/torsions.
subroutine analyse_primitives(mol1,mol2)
use parm
use constant, only: au2ang
implicit none
type(molecule) :: mol1,mol2
integer bond(nat,nat),kk
real(8) s,db(nat),x
real(8) ang1,ang2,anggrad1,anggrad2
integer ii,jj
integer io



! bonding from molecule 1
call bondmatrix(mol1%xyz,mol1%iat,bond)

print*,'******************'
print*,'* BOND LENGTHS   *'
print*,'******************'


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
   write(io,'(a,I5,I5,F8.4)'),'delta_bond',i,j,x
  endif
 enddo
enddo

print*,  ' # covalent bonds ', kk
print'(a,F8.4)',' mean bond deviation [A]     : ', s/dble(kk)
print*,  ' output: bonds.dat'
close(io)

print*,'******************'
print*,'* VALENCE ANGLES *'
print*,'******************'

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
   write(io,'(a,3(I5,x),F8.1)'),'delta_angle',i,j,k,x
   endif
  enddo
 enddo
enddo

print*,  ' # valence angles ', kk
print'(a,F8.1)',' mean angle deviation [deg]  : ', s/dble(kk)
print*,  ' output: angles.dat'
close(io)
print*,'******************'
print*,'* TORSIONS       *'
print*,'******************'
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
       !    if(ii>jj) cycle
       !     print*, i,j,k,l,ii,jj
            if(bond(i,j)==1.and.bond(j,k)==1.and.bond(k,l)==1) then
            kk=kk+1
            call dihed(mol1%xyz,i,j,k,l,ang1,anggrad1)
            call dihed(mol2%xyz,i,j,k,l,ang2,anggrad2)
            call torsionfix(anggrad1,anggrad2,x)
            x=anggrad1-anggrad2
           s=s+x
           write(io,'(a,4(I5,x),F10.1)'),'delta_tors',i,j,k,l,x
            endif
         enddo
      enddo
   enddo
enddo

print*,  ' # torsions ', kk
print'(a,F8.1)',' mean torsion deviation [deg]: ', s/dble(kk)
print*,  ' output: torsions.dat'
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
integer ires

da=.false.
inquire(file=infile,exist=da)
if(.not.da)  return
open(newunit=io,file=infile,status='old')


do
 read(io,'(a)',end=124) aa

 if(index(aa,'$prim').ne.0) then 
 print*,' user definied internal coordinates: '
! n_bond=0
! n_ang=0
! n_dihed=0
  do
    read(io,'(a)',end=123) aa
    if(index(aa,'#').ne.0) cycle

    if(index(aa,'bond').ne.0) then 
     ires=ires+1
     call charXsplit(aa,bb,2)
!      n_bond(ires,1)=s2i(bb)
!      ii=n_bond(ires,1)
     call charXsplit(aa,bb,3)
!      n_bond(ires,2)=s2i(bb)
!      jj=n_bond(ires,2)

      bdist=dbond(xyz0(1,ii),xyz0(1,jj))
!      val0(ires)=bdist
      write(*,'(2x,a,2I4,2x,a,F7.2)') 'bond [A]', ii,jj,' | value: ',bdist
    endif
    if(index(aa,'ang').ne.0) then 
     ires=ires+1
     call charXsplit(aa,bb,2)
!      n_ang(ires,1)=s2i(bb)
!      ii=n_ang(ires,1)
     call charXsplit(aa,bb,3)
!      n_ang(ires,2)=s2i(bb)
!      jj=n_ang(ires,2)
     call charXsplit(aa,bb,4)
!      n_ang(ires,3)=s2i(bb)
!      kk=n_ang(ires,3)
     call angle(xyz0,ii,jj,kk,ang,dang)

     call angle(xyz0,ii,jj,kk,ang,dang)
!     val0(ires)=ang
     write(*,'(2x,a,3I4,2x,a,F7.2)') 'angle [deg]', ii,jj,kk,' | value: ',dang
    endif
    if(index(aa,'dihed').ne.0.or.fstr(aa,'torsion')) then 
     ires=ires+1
     call charXsplit(aa,bb,2)
!      n_dihed(ires,1)=s2i(bb)
!      ii=n_dihed(ires,1)
     call charXsplit(aa,bb,3)
!      n_dihed(ires,2)=s2i(bb)
!      jj=n_dihed(ires,2)
     call charXsplit(aa,bb,4)
!      n_dihed(ires,3)=s2i(bb)
!      kk=n_dihed(ires,3)
     call charXsplit(aa,bb,5)
!      n_dihed(ires,4)=s2i(bb)
!      ll=n_dihed(ires,4)
     call dihed(xyz0,ii,jj,kk,ll,dih,dgrad)
!     val0(ires)=dih
     write(*,'(2x,a,4I4,2x,a,F7.2)') 'torsion [deg]',ii,jj,kk,ll,' | value: ',dgrad
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
print*,''
end subroutine


