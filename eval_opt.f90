subroutine eval_options()
use logic
implicit none
integer i,maxarg
character(120), allocatable :: arg(:)
character(120) ftmp
logical fstr
real(8) s2r
integer s2i

maxarg=iargc()
if(maxarg==0) stop 'get help with: struca -h'

if(maxarg.gt.0) then

 allocate(arg(maxarg))
   do i=1,maxarg
     call getarg(i,arg(i))
   enddo

 do i=1,maxarg
  ftmp=arg(i)
  if(fstr(ftmp,'-h ')) then
   print*,'struca [options]'
   print*,'   options:'
   print*,'   -h                       this help'
   print*,'   -comp <ref. structure> <structure-to-compare>    compare 2 molecules'
!   print*,'   -traj <xyz trajectory>                           analyse trajectory'
   print*,'   -struc <structure>                               analyse singular molecule'
   print*,'   '
   print*,'   -bthr/-athr/-tthr        bond/angle/torsion thresholds'
   print*,'   '
   print*,' <structure> formats: XMOL TMOL  '
   stop
  endif
  if(fstr(ftmp,'-comp ')) then
     do_compare=.true.
    filevec(1)=arg(i+1)
    filevec(2)=arg(i+2)
  endif
  if(fstr(ftmp,'-struc ')) then
     do_single=.true.
    filevec(1)=arg(i+1)
  endif
  if(fstr(ftmp,'-traj '))  then
   do_traj=.true.
   filevec(1)=arg(i+1)
  endif
  if(fstr(ftmp,'-bthr')) thresh_bond=s2r(arg(i+1))
  if(fstr(ftmp,'-athr')) thresh_ang=s2r(arg(i+1))
  if(fstr(ftmp,'-tthr')) thresh_tor=s2r(arg(i+1))
 enddo


endif


end subroutine
