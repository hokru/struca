! print internals
subroutine print_primitives(mol1)
use parm
use constant, only: au2ang
use internals
use logic, only: thresh_bond,thresh_ang,thresh_tor
use omp_lib
use atomdata,only: el
implicit none
type(molecule) :: mol1
integer bond(nat,nat),kk
real(8) s,db(nat),x,rab,xs
real(8) ang1,ang2,anggrad1,anggrad2
integer ii,jj
integer io,tot,it
integer iat(nat)
!real(8) dummy(max_dummy),thresh
!real(8), allocatable :: tvec(:)
!character(*) basename
!character(200) file_out

print*,'************************'
print*,'* PRIMITIVE ANALYSIS   *'
print*,'************************'

! bonding from molecule 1
call bondmatrix(mol1%xyz,mol1%iat,bond)

call read_intcoord('struca.control',mol1%xyz)


print*,'******************'
print*,'* BOND LENGTHS   *'
print*,'******************'


k=0
s=0
kk=0
do i=1,nat-1
 do j=i+1,nat
   k=k+1
   if(i==j) cycle
   if(bond(i,j)==1) then
     kk=kk+1
     x=mol1%dist(k)
     write(*,'(2x,a,x,I5,''['',a2,'']'',I5,''['',a2,'']'',x,F8.4)'),' bond',i,el(mol1%iat(i)),j,el(mol1%iat(j)),x
  endif
 enddo
enddo

if(int_nb>0) print*,' custom bonds: '
do i=1,int_nb
 ii=int_bcast(i,1)
 jj=int_bcast(i,2)
 x=int_bval(i)
 kk=kk+1
 write(*,'(2x,a,x,I5,''['',a2,'']'',I5,''['',a2,'']'',x,F8.4,a3)'),'  bond',ii,el(mol1%iat(ii)),jj,el(mol1%iat(jj)),x
enddo
print*,''

print*,  ' # bonds ', kk

print*,'******************'
print*,'* VALENCE ANGLES *'
print*,'******************'

s=0
kk=0
do i=1,nat-1
 do j=1,nat
  ii=nat*(i-1)+j
  if(i==j) cycle
  do k=i+1,nat
   jj=nat*(j-1)+l
   if(k==j.or.k==i) cycle
   if(bond(i,j)==1.and.bond(j,k)==1) then
    kk=kk+1
    call  angle(mol1%xyz,i,j,k,ang1,anggrad1)
    x=anggrad1
    write(*,'(2x,a,x,3(I5,''['',a2,'']''),x,F8.4)'),' angle',i,el(mol1%iat(i)),j,el(mol1%iat(j)),k,el(mol1%iat(k)),x
   endif
  enddo
 enddo
enddo


if(int_na>0) print*,' custom angles: '
do it=1,int_na
 i=int_acast(it,1)
 j=int_acast(it,2)
 k=int_acast(it,3)
 x=int_aval(i)
    kk=kk+1
 write(*,'(2x,a,x,3(I5,''['',a2,'']''),x,F8.4)'),' angle',i,el(mol1%iat(i)),j,el(mol1%iat(j)),k,el(mol1%iat(k)),x
enddo
print*,''


print*,  ' # angles ', kk

print*,'******************'
print*,'* TORSIONS       *'
print*,'******************'


s=0
kk=0
do i=1,nat-1
  do j=1,nat
    if(i==j) cycle
    if(bond(i,j)==0) cycle 
    do k=1,nat
      ii=nat*(i-1)*k
      if(k==j.or.k==i) cycle
      if(bond(j,k)==0) cycle
         do l=i+1,nat
           jj=nat*(j-1)*l
           if(l==k.or.l==j.or.l==i) cycle
            if(bond(k,l)==1) then
             kk=kk+1
             call dihed(mol1%xyz,i,j,k,l,ang1,anggrad1)
             x=anggrad1
             write(*,'(2x,a,x,4(I5,''['',a2,'']''),x,F8.4)'),' torsion',i,el(mol1%iat(i)),j,el(mol1%iat(j)),k,el(mol1%iat(k)),l,el(mol1%iat(l)),x
            endif
         enddo
      enddo
   enddo
enddo


if(int_nt>0) print*,' custom torsions: '
do it=1,int_nt
 i=int_tcast(it,1)
 j=int_tcast(it,2)
 k=int_tcast(it,3)
 l=int_tcast(it,4)
 x=int_tval(it)
 kk=kk+1
 write(*,'(2x,a,x,4(I5,''['',a2,'']''),x,F8.4)'),' torsion',i,el(mol1%iat(i)),j,el(mol1%iat(j)),k,el(mol1%iat(k)),l,el(mol1%iat(l)),x
enddo
print*,''


print*,  ' # torsions ', kk

end subroutine

