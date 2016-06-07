! according to Seok,Coutsias,Dill
! cite: DOI: 10.1002/jcc.20110

! quarternion fit,
! nat # number of atoms
! xyz1(3,nat) of molecule 1 (reference)
! xyz2(3,nat) of molecule 2 (to be aligned, overwritten with rotated coordinates)

! not sure how to reverse the order
subroutine quaternion_fit(nat,xyz1,xyz2)
implicit none
integer i,j,nat,iat(nat)
real(8) xx(nat),yy(nat)
real(8) q(4), rot(3,3), p(4,4), ev(4),coord(3,nat)
real(8) xyz1(3,nat),xyz2(3,nat)
real(8) coord1(3,nat),coord2(3,nat)
real(8) ddot
real(8) cc(3),rr(3)
real(8) xrot,yrot,zrot
real(8) c1(3),c2(3),dum1(nat),dum2(nat)
real(8) R(3,3),S(4,4),xnorm,ynorm,lam
real(8) rmsd

!coord1=xyz1
!coord2=xyz2

print*,''
print*,''
print*,' Aligning mol2 -> mol1'
print*,''
print*,'    quaternion RMSD & alignment:'
print*,'      E.A. Coutsias, C. Seok, K.A.Dill J. Comput. Chem.,2004,25,1849-1857'
print*,'      DOI:10.1002/jcc.20110'


xnorm=0d0
ynorm=0d0
! move to centroids
do i=1,3
  xx(1:nat)=xyz1(i,1:nat)
  yy(1:nat)=xyz2(i,1:nat)
  c1(i)=sum(xx(1:nat))/dble(nat)
  c2(i)=sum(yy(1:nat))/dble(nat)
  coord1(i,1:nat)=xyz1(i,1:nat)-c1(i) !centroid translated coords 1
  coord2(i,1:nat)=xyz2(i,1:nat)-c2(i) !centroid translated coords 2
  xnorm=xnorm+ddot(nat,coord1(i,:),1,coord1(i,:),1)
  ynorm=ynorm+ddot(nat,coord2(i,:),1,coord2(i,:),1)
enddo

!call centroid(3,nat,xyz1,c1)
!call centroid(3,nat,xyz2,c2)
!print*,xnorm,ynorm
!print*,c1
!print*,c2


! here we chose which molecule we want to rotate onto which
 do i=1,3
  do j=1,3
!    R(i,j) = ddot(nat,coord1(i,:),1,coord2(j,:),1)  ! 1->2
    R(i,j) = ddot(nat,coord2(i,:),1,coord1(j,:),1)   ! 2->1
  end do
 end do

 ! S matrix
 S(1, 1) = R(1, 1) + R(2, 2) + R(3, 3)
 S(2, 1) = R(2, 3) - R(3, 2)
 S(3, 1) = R(3, 1) - R(1, 3)
 S(4, 1) = R(1, 2) - R(2, 1)

 S(1, 2) = S(2, 1)
 S(2, 2) = R(1, 1) - R(2, 2) - R(3, 3)
 S(3, 2) = R(1, 2) + R(2, 1)
 S(4, 2) = R(1, 3) + R(3, 1)

 S(1, 3) = S(3, 1)
 S(2, 3) = S(3, 2)
 S(3, 3) =-R(1, 1) + R(2, 2) - R(3, 3)
 S(4, 3) = R(2, 3) + R(3, 2)

 S(1, 4) = S(4, 1)
 S(2, 4) = S(4, 2)
 S(3, 4) = S(4, 3)
 S(4, 4) =-R(1, 1) - R(2, 2) + R(3, 3) 



! diagonalize, return eigenvectors in s, eigenvalues in ev
! destroys s
call DiagSM(4,s,ev)

!print*,'ev',ev
!  quaternion for rotation corresponds to last and largest eigenvalue
q(1:4)=s(1:4,4)

! this RMSD is correct!
print*,''
rmsd=sqrt(max(0.0d0,( (xnorm+ynorm)-2.0d0*ev(4) ))/dble(nat))
write(*,'(2x,''RMSD (qfit): [ang] '',F8.4)') rmsd
print*,''

!print*,' RMSD (qfit) [A]: ',sqrt(max(0.0d0,( (xnorm+ynorm)-2.0d0*ev(4) ))/dble(nat))
!print*,'q',q

! translate to rotation matrix
call q_to_rot(q,rot)

! rot(1,1)=q(1)**2 + q(2)**2 - q(3)**2 - q(4)**2
! rot(2,2)=q(1)**2 - q(2)**2 + q(3)**2 - q(4)**2
! rot(3,3)=q(1)**2 - q(2)**2 - q(3)**2 + q(4)**2
! rot(1,2)=2d0*(q(2)*q(3)-q(1)*q(4))
! rot(1,3)=2d0*(q(2)*q(4)+q(1)*q(3))
! rot(2,1)=2d0*(q(2)*q(3)+q(1)*q(4))
! rot(3,1)=2d0*(q(2)*q(4)-q(1)*q(3))
! rot(2,3)=2d0*(q(3)*q(4)-q(1)*q(2))
! rot(3,2)=2d0*(q(3)*q(4)+q(1)*q(2))

! compared well with Dill's code
!print*,'rotation matrix',rot


! translate & rotate 
! select correct molecule
do i=1,nat
    xyz2(1:3,i)=matmul(rot,coord2(:,i))+c1  ! rotate mol 2
!    xyz1(1:3,i)=matmul(rot,coord1(:,i))+c2 ! rotate mol 1
enddo

return
end subroutine

! forms rotation matrix from quarternion
subroutine q_to_rot(q,rot)
implicit none
real(8), intent(in)    :: q(4)
real(8), intent(inout) :: rot(3,3)

rot(1,1)=q(1)**2 + q(2)**2 - q(3)**2 - q(4)**2
rot(2,2)=q(1)**2 - q(2)**2 + q(3)**2 - q(4)**2
rot(3,3)=q(1)**2 - q(2)**2 - q(3)**2 + q(4)**2
rot(1,2)=2d0 * (q(2)*q(3) - q(1)*q(4))
rot(1,3)=2d0 * (q(2)*q(4) + q(1)*q(3))
rot(2,1)=2d0 * (q(2)*q(3) + q(1)*q(4))
rot(3,1)=2d0 * (q(2)*q(4) - q(1)*q(3))
rot(2,3)=2d0 * (q(3)*q(4) - q(1)*q(2))
rot(3,2)=2d0 * (q(3)*q(4) + q(1)*q(2))
return

end subroutine
