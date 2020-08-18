!
!    geom_util is free software: you can redistribute it and/or modify
!    it under the terms of the GNU Lesser General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    geom_util is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public License
!    along with geom_util.  If not, see <https://www.gnu.org/licenses/>.
!
! Copyright (C) 2019 - 2020 Holger Kruse, Institute of Biophysics of the CAS, Czechia

module io_unit
integer iounit
end module


subroutine geom_util(nat,xyz,iat,options,idstring,lib_io)
use io_unit
implicit none

integer ( kind = 4 ), parameter :: dim_num = 3
real ( kind = 8 ) centroid(dim_num)

real(8), intent(in) :: xyz(3,nat)
integer(4), intent(in) :: nat
integer(4), intent(in) :: iat(nat)
character(256), intent(in) :: options(*)
character(256), intent(in) :: idstring

real(8), allocatable:: v(:,:)
integer, allocatable :: al(:)
integer :: i,j,ial, lib_io
character(120) :: ns,frmt,mode

! FOR PLANE DISTANCE
integer nat1,nat2,k
real(8), allocatable :: coord1(:,:), coord2(:,:)
!integer, allocatable :: ic1(:),ic2(:)

integer n ! polygon size

iounit=lib_io

write(iounit,*) ''
write(iounit,*) 'x-------------------------x'
write(iounit,*) 'x       geom_util lib     x'
write(iounit,*) 'x-------------------------x'
write(iounit,*) 'version: 0.4b'
write(iounit,*) ''
write(iounit,*) ' computes controid of a 3d polygon'
write(iounit,*) ' e.g. for a ring in a molecule'
write(iounit,*) ' and a fitted distance between fitted planes'
! write(iounit,*) 'usage:'
! write(iounit,*) 'options: '
! write(iounit,*) '  plane <index list plane1> <index list plane2>'
! write(iounit,*) '  centroid <int ring size> <index list of ring>'
! write(iounit,*) ''
! write(iounit,*) ' <file xyz> = xmol or tmol, prints centroid in given units!'
! write(iounit,*) ''
! write(iounit,*) '!! atom numbers must be given in order of the ring !!'
! write(iounit,*) '!! (clock -or anticlock-wise does not matter)      !!'
! write(iounit,*) ''

write(iounit,*) 'ID: '//trim(idstring)


mode=trim(options(1))
write(iounit,*) 'MODE:',trim(mode)


if(trim(mode)=='centroid') then
write(iounit,*) " CENTROID MODE"
! ring size
read(options(2),*) n
write(iounit,'(a,I2)')' --> ring size: ', n

! ring info
allocate(al(n))
call atlist(options(3),al,ial)

allocate(v(3,n))

! order atom list
!call int_bsort(n,al)
write(ns,*) n
write(frmt,'(a)') "(a,2x,"//trim(adjustl(ns))//"(I1,x))"
write(iounit,frmt) ' --> ring atoms:',al

! assign polygon vertices
do j=1,n
 do i=1,nat
  if(i==al(j)) then
   v(1:3,j)=xyz(1:3,i)
  endif
 enddo
enddo


call r8mat_transpose_print ( dim_num, n, v, '  ring coordinates:' )

! compute!
call polygon_centroid_3d ( n, v, centroid )
write(iounit,*)  'CENTROID'
write(iounit,*)  centroid

elseif(trim(mode)=='plane') then
write(iounit,*) " PLANE MODE"


allocate(al(100))


! ATOMS OF PLANE 1
al=0
call atlist(options(2),al,ial)

nat1=0
do i=1,nat
 if(any(al==i)) nat1=nat1+1
enddo
write(iounit,*) 'Atoms plane 1', nat1
allocate(coord1(3,nat1))

k=0
do i=1,nat
 if(any(al==i)) then
   k=k+1
   coord1(1:3,k)=xyz(1:3,i)
 endif
enddo

! ATOMS OF PLANE 1
al=0
call atlist(options(3),al,ial)

nat2=0
do i=1,nat
 if(any(al==i)) nat2=nat2+1
enddo
write(iounit,*) 'Atoms plane 2', nat2
allocate(coord2(3,nat2))

k=0
do i=1,nat
 if(any(al==i)) then
   k=k+1
   coord2(1:3,k)=xyz(1:3,i)
 endif
enddo


write(iounit,*) coord1
write(iounit,*) ''
write(iounit,*) coord2

call dist_plane(nat1,nat2,coord1,coord2)
endif






end



subroutine bsort(n,e)
! bubble sort the vector e
implicit none
integer i,j,k,l,nn,n,ii
real(8) e(n),tt
character(80) cc
logical order

nn=n
order=.false.
do
if(order) exit
order=.true.
 do i=1,nn-1
    if (e(i).gt.e(i+1) ) then ! swap
      tt=e(i)
      e(i)=e(i+1)
      e(i+1)=tt
      order = .false.
     endif
 enddo
nn=nn-1
enddo

return
end subroutine


subroutine int_bsort(n,e)
! bubble sort the vector e
implicit none
integer i,j,k,l,nn,n,ii
integer e(n),tt
character(80) cc
logical order

nn=n
order=.false.
do
if(order) exit
order=.true.
 do i=1,nn-1
    if (e(i).gt.e(i+1) ) then ! swap
      tt=e(i)
      e(i)=e(i+1)
      e(i+1)=tt
      order = .false.
     endif
 enddo
nn=nn-1
enddo

return
end subroutine


