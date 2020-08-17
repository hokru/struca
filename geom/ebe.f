!    This file is part of geom_util.
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
! Copyright (C) Stefan Grimme, date unknow.


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c determinant of 3x3 matrix                                           c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      real*8 function det3(a)
      implicit real*8 (a-h,o-z)
      dimension a(3,3)

      det=a(1,1)*a(2,2)*a(3,3)+
     .    a(1,2)*a(2,3)*a(3,1)+
     .    a(1,3)*a(2,1)*a(3,2)-
     .    a(3,1)*a(2,2)*a(1,3)-
     .    a(3,2)*a(2,3)*a(1,1)-
     .    a(3,3)*a(2,1)*a(1,2)
      det3=-1.0d0*det
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c plane with 3 points x,y,z(1)                                        c
c                     x,y,z(2)                                        c
c                     x,y,z(3)                                        c
c plane is aa*x + bb*y + cc*z + dd = 0                                c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine febene(x,y,z,aa,bb,cc,dd)
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

      aa=det3(a)

      a(1,1)=x(1)
      a(1,2)=z(1)
      a(1,3)=1.0d0
      a(2,1)=x(2)
      a(2,2)=z(2)
      a(2,3)=1.0d0
      a(3,1)=x(3)
      a(3,2)=z(3)
      a(3,3)=1.0d0

      bb=-det3(a)

      a(1,1)=x(1)
      a(1,2)=y(1)
      a(1,3)=1.0d0
      a(2,1)=x(2)
      a(2,2)=y(2)
      a(2,3)=1.0d0
      a(3,1)=x(3)
      a(3,2)=y(3)
      a(3,3)=1.0d0

      cc=det3(a)

      a(1,1)=x(1)
      a(1,2)=y(1)
      a(1,3)=z(1)
      a(2,1)=x(2)
      a(2,2)=y(2)
      a(2,3)=z(2)
      a(3,1)=x(3)
      a(3,2)=y(3)
      a(3,3)=z(3)

      dd=-det3(a)
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c plane with 3 points x,y,z(1)                                        c
c                     x,y,z(2)                                        c
c                     x,y,z(3)                                        c
c distance to point x(i),y(i),z(i)                                    c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine pldist(i,x,y,z,s,t,u,rx)
      implicit real*8 (a-h,o-z)
      dimension x(*), y(*), z(*)
      dimension s(3), t(3), u(3)
      call febene(s,t,u,aa,bb,cc,dd)
      rx=(aa*x(i)+bb*y(i)+cc*z(i)+dd)/dsqrt(aa*aa+bb*bb+cc*cc)
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine rms(n,x,y,z,s,t,u,r)
      implicit real*8 (a-h,o-z)
      dimension x(*), y(*), z(*)
      dimension s(3), t(3), u(3)

      r=0.0d0
      do i=1,n
         call pldist(i,x,y,z,s,t,u,rx)
         r=r+rx*rx
      enddo
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine grad(n,x,y,z,s,t,u,g)
      implicit real*8 (a-h,o-z)
      dimension x(*), y(*), z(*), g(*)
      dimension s(3), t(3), u(3)
      thr=5.0d-4

      do i=1,3
      s(i)=s(i)+thr
      call rms(n,x,y,z,s,t,u,rr)
      s(i)=s(i)-2.*thr
      call rms(n,x,y,z,s,t,u,rl)
      g(i)=0.5*(rr-rl)/thr
      s(i)=s(i)+thr
      enddo

      do i=1,3
      t(i)=t(i)+thr
      call rms(n,x,y,z,s,t,u,rr)
      t(i)=t(i)-2.*thr
      call rms(n,x,y,z,s,t,u,rl)
      g(i+3)=0.5*(rr-rl)/thr
      t(i)=t(i)+thr
      enddo

      do i=1,3
      u(i)=u(i)+thr
      call rms(n,x,y,z,s,t,u,rr)
      u(i)=u(i)-2.*thr
      call rms(n,x,y,z,s,t,u,rl)
      g(i+6)=0.5*(rr-rl)/thr
      u(i)=u(i)+thr
      enddo

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine plfit(n,x,y,z,aa,bb,cc,dd)
      use io_unit
      implicit real*8 (a-h,o-z)
       real*8 bohr
       integer maxdiag
       PARAMETER (BOHR=0.52917726)
       PARAMETER (MAXDIAG=3)
      dimension x(*), y(*), z(*)
      dimension s(3), t(3), u(3), g(9)

      thr=1.d-5
      alp=0.1


      s(1)=x(1)+thr
      s(2)=x(2)-thr
      s(3)=x(3)+thr
      t(1)=y(1)+thr
      t(2)=y(2)-thr
      t(3)=y(3)+thr
      u(1)=z(1)+thr
      u(2)=z(2)-thr
      u(3)=z(3)+thr

      call rms(n,x,y,z,s,t,u,r0)
      r0=dsqrt(r0/(n-1))

      it=0
      r0=0.0d0
 10   continue
      it=it+1
      call grad(n,x,y,z,s,t,u,g)
      gn=0.0d0
      do i=1,3
         gn=gn+g(i)**2+g(i+3)**2+g(i+6)**2
         s(i)=s(i)-alp*g(i)
         t(i)=t(i)-alp*g(i+3)
         u(i)=u(i)-alp*g(i+6)
      enddo
      call rms(n,x,y,z,s,t,u,r)
      r=dsqrt(r/(n-1))
      !write(iounit,'(i4,''  rms '',d14.6,''   norm g'',d14.6)') it, r, gn

      if(r-r0.gt.0) alp=alp*0.3
      if(r-r0.lt.0) alp=alp*1.3
      if(alp.gt.0.5) alp=1.0
      if(alp.lt.0.001) alp=0.01


      dr=dabs(r-r0)
      r0=r
      if(gn.gt.1.0d-7.or.dr.gt.1.0d-7) then
         goto 10
      endif


      write(iounit,*)
      sum1=0.0d0
      sum2=0.0d0
      do i=1,n
         call pldist(i,x,y,z,s,t,u,rx)
         write(iounit,'(i4,''.atom, dev (ang.)'',f14.6)') i,rx
         sum1=sum1+    rx 
         sum2=sum2+abs(rx)
      enddo

      write(iounit,'(/''mean     dev (ang.)'',f14.6)') sum1/n
      write(iounit,'(/''mean abs dev (ang.)'',f14.6)') sum2/n

      f=1./bohr
      write(iounit,'(/'' plane coordinates (au)''/)')
      write(iounit,'('' x, y, z'',3f14.6)') s(1)*f,t(1)*f,u(1)*f
      write(iounit,'('' x, y, z'',3f14.6)') s(2)*f,t(2)*f,u(2)*f
      write(iounit,'('' x, y, z'',3f14.6)') s(3)*f,t(3)*f,u(3)*f

      call febene(s,t,u,aa,bb,cc,dd)

      write(iounit,'(/'' plane parameters A B C D'',4F14.6/)')
     & aa,bb,cc,dd

      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  test
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     implicit real*8 (a-h,o-z)
c     dimension x(100), y(100), z(100)
c     read(*,*) n
c     do i=1,n
c        read(*,*)x(i),y(i),z(i)
c     enddo
c     call plfit(n,x,y,z)
c     end
