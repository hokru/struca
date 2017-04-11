! RDF g(r)=1/Np Sum(i) SUm(k/=i) <delta (r+rk-ri)>
! PCF g(r)=1/4pir**2*1/Np  Sum(i) SUm(k/=i)  <delta (r-|rk-ri|)>


subroutine calc_rdf(traj)
use parm, only: nat, npair,trajectory
type(trajectory) traj
real(8) dr,r,rab
integer i,j


end subroutine

subroutine dist_ana(base,traj)
use logic
use parm
implicit none
type(trajectory) traj
integer io
character(*) base
character(255) fname
character(2) esym

fname=trim(base)//'_traj.dist.dat'
open(newunit=io,file=trim(fname))
write(*,*) ' DISTANCE ANALYSIS'
write(*,*) ' --> ',trim(fname)

  
do i=1,traj%nat-1
    do j=i+1,traj%nat
        if(i==j) cycle
        ! print*,i,j
        write(io,'(a,I4,x,a,a2,a,x,I4,x,a,a2,a)') 'atom pair:', i,'[',esym(traj%iat(i)),']',j,'[',esym(traj%iat(j)),']'
        adist=0d0
        maxd=-HUGE(0d0)
        mind=HUGE(0d0)
        do k=1,traj%nmol
            ra(1:3)=traj%mxyz(1:3,i,k)
            rb(1:3)=traj%mxyz(1:3,j,k)
            rab=dbond(ra,rb)
            adist=adist+rab
            if(rab>maxd) maxd=rab
            if(rab<mind) mind=rab
        enddo
        write(io,'(4x,a,4(F9.4,x))') 'avg/max/min/span :', adist/nmol,maxd,mind,abs(mind)+maxd
    enddo
enddo


close(io)
end subroutine
