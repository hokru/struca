subroutine unitvec(x,e)
implicit none
real(8) x(3),e(3),t(3)
t=DOT_PRODUCT(x,x)
e=x/sqrt(t)
end

subroutine veclen(x,v)
implicit none
real(8) x(3),v
v=sqrt(dot_product(x,x))
end subroutine


! distance between vectors a(3),b(3)
real(8) function rab(a,b)
implicit none
real(8), intent(in) :: a(3),b(3)
real(8) dx,dy,dz
 dx=a(1)-b(1)
 dy=a(2)-b(2)
 dz=a(3)-b(3)
 rab=sqrt(dx*dx+dy*dy+dz*dz)
end function

! distance between cartesian vectors a and b
real(8) pure function dbond(a,b)
implicit none
real(8), intent(in) :: a(3),b(3)
integer i
real(8) ab(3)
ab=a-b
dbond=sqrt( ab(1)**2 + ab(2)**2 + ab(3)**2 )
end function


subroutine veclen2(a,b,v)
implicit none
real(8) a(3),b(3),v,x(3)
x=a-b
v=sqrt(dot_product(x,x))
end subroutine


subroutine cross_prod(y,x2,x3)
implicit none
real(8) y(3),x2(3),x3(3)
  y(1) =  x2(2)*x3(3) - x3(2)*x2(3)
  y(2) = -x2(1)*x3(3) + x3(1)*x2(3)
  y(3) =  x2(1)*x3(2) - x3(1)*x2(2)
end subroutine


subroutine angle(xyz,aa,bb,cc,ang,anggrad)
use constant
implicit none
integer aa,bb,cc
real(8) xyz(3,*),v1(3),v2(3),lv1,lv2,uv1(3),uv2(3)
real(8) t,tt,ang,anggrad
real(8) cv(3),lcv

v1=xyz(1:3,bb)-xyz(1:3,aa)
v2=xyz(1:3,bb)-xyz(1:3,cc)

call cross_prod(cv,v1,v2)
call veclen(cv,lcv)
ang=atan2(lcv, dot_product(v1,v2))
anggrad=ang*180.0d0/pi

end subroutine

subroutine dihed(xyz,aa,bb,cc,dd,dih,dihgrad)
use constant
implicit none
integer aa,bb,cc,dd
real(8) b1(3),b2(3),b3(3),n1(3),n2(3)
real(8) un1(3),un2(3),ub2(3),m(3),um(3),dix,diy
real(8) dih,dihgrad,xyz(3,*)

 b1=xyz(1:3,aa)-xyz(1:3,bb)
 b2=xyz(1:3,bb)-xyz(1:3,cc)
 b3=xyz(1:3,cc)-xyz(1:3,dd)

 ! normal of the planes
 call cross_prod(n1,b1,b2)
 call cross_prod(n2,b2,b3)

 call unitvec(n1,un1)
 call unitvec(n2,un2)
 call unitvec(b2,ub2)

 call cross_prod(m,un1,ub2)
 call unitvec(m,um)

 dix=DOT_PRODUCT(un1,un2)
 diy=DOT_PRODUCT(um,un2)

 dih=atan2(diy,dix)

!  Quadrant    Angle              sin    cos    tan
!----------------------------------------------------
!  I           0    < Î± < Ï/2     > 0    > 0    > 0       
!  II          Ï/2  < Î± < Ï       > 0    < 0    < 0
!  III         Ï    < Î± < 3Ï/2    < 0    < 0    > 0
!  IV          3Ï/2 < Î± < 2Ï      < 0    > 0    < 0
! atan2(0,1) =   0
! atan2(1,0) =   pi/2
! atan2(-1,0) = -pi/2
! atan2(0,-1) =  pi



! give results in 0 to 360 degree
 if(dih<0.0d0) dih=dih+pi*2
 dihgrad=dih*180.0d0/pi
end subroutine


real(8) function di360(x)
implicit none
real(8), parameter:: pi = 3.141592653589793d0
real(8) x
if(x<0.0d0) x=x+pi*2
di360=x*180.0d0/pi
end function

real(8) function grad2rad(x)
implicit none
real(8), parameter:: pi = 3.141592653589793d0
real(8) x
grad2rad=x*pi/180.0d0
end function

real(8) function anggrad(x)
implicit none
real(8) x
real(8), parameter:: pi = 3.141592653589793d0
anggrad=x*180.0d0/pi
end function


subroutine DiagSM(xvar,mat,eig)
implicit none
integer i,j,k
real(8), allocatable :: aux(:)
integer info,lwork,xvar
real(8) ,intent(inout) :: mat(xvar,xvar)
real(8) xx(1)
real(8), intent(out) :: eig(xvar)

eig=0
call dsyev ('V','U',xvar,mat,xvar,eig,xx,-1,info)
lwork=int(xx(1))
allocate(aux(lwork))
call dsyev ('V','U',xvar,mat,xvar,eig,aux,lwork,info)
if(info/=0) print*,'Diagonalization failed !!'
end subroutine

