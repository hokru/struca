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
subroutine dist_plane(nat1,nat2,coord1,coord2)
use io_unit
implicit none
character(40) a8                                                           
character(80) atmp                                                         
integer nat1,nat2,i
real(8) aa1,aa2,bb2,bb1,cc1,cc2,dd1,dd2,rx,rsum,cosph
real(8) xyz(3,4), xx(500), s(3), t(3), u(3)                              
real(8) coord1(3,nat1), coord2(3,nat2)
real(8) xxx(500), yyy(500), zzz(500)
!real(8) ic1(nat1),ic2(nat2) 
                                                

xxx(1:nat1)=coord1(1,1:nat1)
yyy(1:nat1)=coord1(2,1:nat1)
zzz(1:nat1)=coord1(3,1:nat1)
if(nat1>1) call plfit(nat1,xxx,yyy,zzz,aa1,bb1,cc1,dd1)


xxx(1:nat2)=coord2(1,1:nat2)
yyy(1:nat2)=coord2(2,1:nat2)
zzz(1:nat2)=coord2(3,1:nat2)
if(nat2>1) call plfit(nat2,xxx,yyy,zzz,aa2,bb2,cc2,dd2)

if(nat1>1.and.nat2>2) then 
  cosph=(aa1*aa2+bb1*bb2+cc1*cc2)/sqrt((aa1**2+bb1**2+cc1**2)*(aa2**2+bb2**2+cc2**2))
 write(iounit,'(''**********************************************'')')
 write(iounit,'('' angle between planes'',f10.4,5x,a)') acos(cosph)*180.0d0/3.1415927
else
 
endif


rsum=0.0d0
if(nat2>1) then
 write(iounit,'('' distance set 1 atoms to plane 2'')')
 do i=1,nat1
    rx=(aa2*coord1(1,i)+bb2*coord1(2,i)+cc2*coord1(3,i)+dd2)/dsqrt(aa2*aa2+bb2*bb2+cc2*cc2)
    write(iounit,*) i,rx
    rsum=rsum+abs(rx)
 enddo
endif
if(nat1>1) then
 write(iounit,'('' distance set 2 atoms to plane 1'')')
 do i=1,nat2
    rx=(aa1*coord2(1,i)+bb1*coord2(2,i)+cc1*coord2(3,i)+dd1)/dsqrt(aa1*aa1+bb1*bb1+cc1*cc1)
    write(iounit,*) i,rx
    rsum=rsum+abs(rx)
 enddo
endif

if(nat1>1.and.nat2>1) write(iounit,'('' average distance between planes (Ang.)'',F12.4)') rsum/(nat1+nat2)
write(iounit,'(''**********************************************'')')

end subroutine






!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc         
!c plane with 3 points x,y,z(1)                                        c         
!c                     x,y,z(2)                                        c         
!c                     x,y,z(3)                                        c         
!c plane is aa*x + bb*y + cc*z + dd = 0                                c         
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc         
                                                                                
subroutine ebene(x,y,z,aa,bb,cc,dd)                                       
implicit real*8 (a-h,o-z)                                                 
dimension x(*), y(*), z(*), a(3,3)                                        
                                                                          
a(1,1)=y(1)                                                               
a(1,2)=z(1)                                                               
a(1,3)=1.0d0                                                              
a(2,1)=y(2)                                                               
a(2,2)=z(2)                                                               
a(2,3)=1.0d0                                                              
a(3,1)=y(3)                                                               
a(3,2)=z(3)                                                               
a(3,3)=1.0d0                                                              
                                                                          
aa=det(a)                                                                 
                                                                          
a(1,1)=x(1)                                                               
a(1,2)=z(1)                                                               
a(1,3)=1.0d0                                                              
a(2,1)=x(2)                                                               
a(2,2)=z(2)                                                               
a(2,3)=1.0d0                                                              
a(3,1)=x(3)                                                               
a(3,2)=z(3)                                                               
a(3,3)=1.0d0                                                              
                                                                          
bb=-det(a)                                                                
                                                                          
a(1,1)=x(1)                                                               
a(1,2)=y(1)                                                               
a(1,3)=1.0d0                                                              
a(2,1)=x(2)                                                               
a(2,2)=y(2)                                                               
a(2,3)=1.0d0                                                              
a(3,1)=x(3)                                                               
a(3,2)=y(3)                                                               
a(3,3)=1.0d0                                                              
                                                                          
cc=det(a)                                                                 
                                                                          
a(1,1)=x(1)                                                               
a(1,2)=y(1)                                                               
a(1,3)=z(1)                                                               
a(2,1)=x(2)                                                               
a(2,2)=y(2)                                                               
a(2,3)=z(2)                                                               
a(3,1)=x(3)                                                               
a(3,2)=y(3)                                                               
a(3,3)=z(3)                                                               
                                                                          
dd=-det(a)                                                                
return                                                                    
end    

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc         
!c determinant of 3x3 matrix                                           c         
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc         
                                                                                
real(8) function det(a)                                                    
implicit real*8 (a-h,o-z)                                                 
dimension a(3,3)                                                          
                                                                          
det=a(1,1)*a(2,2)*a(3,3)+ &
    a(1,2)*a(2,3)*a(3,1)+ &
    a(1,3)*a(2,1)*a(3,2)- &
    a(3,1)*a(2,2)*a(1,3)- &
    a(3,2)*a(2,3)*a(1,1)- &
    a(3,3)*a(2,1)*a(1,2)  
det=-1.0d0*det                                                            
return                                                                    
end     
