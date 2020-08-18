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
write(*,*) ' ATOM-PAIR DISTANCE ANALYSIS'
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
        write(io,'(4x,a,4(F9.4,x))') 'avg/max/min/span :', adist/nmol,maxd,mind,maxd-mind
    enddo
enddo


close(io)
end subroutine


subroutine ana_bond(base,traj,ii,jj)
use logic
use parm
implicit none
type(trajectory) traj
integer io,ii,jj
character(*) base
character(255) fname
character(2) esym

fname=trim(base)//'_traj.bond.dat'
open(newunit=io,file=trim(fname))
write(*,*) ' DISTANCE ANALYSIS'
write(*,*) ' --> ',trim(fname)

write(io,'(a,I4,x,a,a2,a,x,I4,x,a,a2,a)') '# atom pair distance:', ii,'[',esym(traj%iat(ii)),']',jj,'[',esym(traj%iat(jj)),']'

do k=1,traj%nmol
    ra(1:3)=traj%mxyz(1:3,ii,k)
    rb(1:3)=traj%mxyz(1:3,jj,k)
    rab=dbond(ra,rb)
    write(io,'(4x,F9.4,x)')  rab
enddo

close(io)
end subroutine



subroutine external(base,traj)
use logic
use parm

type(trajectory) traj
integer io, p
character(*) base
character(255) fname,runme
character(2) esym
real(8) nstep

fname=trim(base)//'_traj.external.dat'
open(newunit=io,file=trim(fname),status='replace')
write(*,*) ' EXTERNAL PROGRAM ANALYSIS'
write(*,*) ' --> ',trim(fname)
runme=trim(ecommand)//' >> '//trim(fname)
print*,'External command:', runme

do i=1,traj%nmol
    call counter(i,traj%nmol)
    call wrxyz(traj%iat,traj%nat,traj%mxyz(:,:,i),'tmp.xyz',.false.)
    call execute_command_line(runme)
enddo


end subroutine


subroutine geom_util_lib(base,traj)
use logic
use parm
implicit none
type(trajectory) traj
integer io, iword
character(*) base
character(256) fname, runme, idstring
character(2) esym
real(8) nstep

fname=trim(base)//'_traj.geom.dat'
open(newunit=io,file=trim(fname),status='replace')
write(*,*) ' geom_util library ANALYSIS'
write(*,*) ' --> ',trim(fname)
print*,'options: '
print*,' -- '
do iword=1,size(options)
  print*,trim(options(iword))
enddo
print*,' -- '
print*,traj%nmol

nstep=10d0
do i=1,traj%nmol
    call counter(i,traj%nmol)
    nstep=nstep+10d0
    write(idstring,'(a,I7)') "molecule",i
    call geom_util(traj%nat,traj%mxyz(:,:,i),traj%iat,options,idstring,io)
enddo
end subroutine


subroutine counter(i,n)
use helper, only: nstep
implicit none
integer i,n
real(8) p
p=100d0*dble(i)/dble(n)
if(p>=nstep ) then
    nstep=nstep+10.0d0
    write(*,'(I3,A)') nint(p),' % done'
endif
end 