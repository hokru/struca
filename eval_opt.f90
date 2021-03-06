subroutine eval_options()
use logic
implicit none
integer i,maxarg
character(120), allocatable :: arg(:)
character(120) ftmp
logical fstr
real(8) s2r
integer s2i
integer iword,n_options

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
   print*,'   '
   print*,'   -traj <xyz trajectory>                           analyse trajectory'
   print*,'   -tdist                                           bond statistics'
   print*,'   -tpair <int> <int>                               atom-pair distance '
   print*,'   -traj <xyz trajectory>                           analyse trajectory'
   print*,'   '
   print*,'   -struc <structure>                               analyse singular molecule'
   print*,'   '
   print*,'   -bthr/-athr/-tthr        bond/angle/torsion thresholds'
   print*,'        defaults: [0.05 A / 1.0 deg / 5.0 deg ] '
   print*,'    '
   print*,'   -hb_r/-hb_a       distance/angle thresholds for Hbond detection '
   print*,'       defaults: [ 3.0 A / 120 deg ] '
   print*,'    '
   print*,'    '
   print*,'   -noH                                              exclude all H atoms '
   print*,'   -bondfactor <real>                                scaling factor for bond detection (def: 1.1) '
   print*,'    '
   print*,'   -ext <string>                                     executes <string> for every structure in a trajectory'
   print*,'                                                     a local XMOL file called tmp.xyz will be written to disk'
   print*,'    '
   print*,'   -geom <string>                                    executes geom_util library with given option string'
   print*,'                                                       options: '
   print*,'                                                       * plane <index list plane1> <index list plane2>'
   print*,'                                                       * centroid <int ring size> <index list of ring>'
   print*,'                                                     !! atom numbers must be given in order of the ring !!'
   print*,'                                                     !! (clock -or anticlock-wise does not matter)      !!'
   print*,'    '
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
  if(fstr(ftmp,'-frag ')) then
     do_frag=.true.
    filevec(1)=arg(i+1)
  endif
  if(fstr(ftmp,'-noH '))  excludeH=.true.
  if(fstr(ftmp,'-onlyEL ')) then
      onlyEL=.true.
      iel=s2i(arg(i+1))
  endif
  if(fstr(ftmp,'-traj '))  then
   do_traj=.true.
   filevec(1)=arg(i+1)
  endif
  if(fstr(ftmp,'-tdist')) then
   traj_dist=.true.
  endif
  if(fstr(ftmp,'-tpair')) then
   traj_pair=.true.
   ia=s2i(arg(i+1))
   ja=s2i(arg(i+2))
  endif
  ! primitive thr
  if(fstr(ftmp,'-bthr')) thresh_bond=s2r(arg(i+1))
  if(fstr(ftmp,'-athr')) thresh_ang=s2r(arg(i+1))
  if(fstr(ftmp,'-tthr')) thresh_tor=s2r(arg(i+1))
  ! Hbond thr
  if(fstr(ftmp,'-hb_r')) thresh_hbr=s2r(arg(i+1))
  if(fstr(ftmp,'-hb_a')) thresh_hba=s2r(arg(i+1))
  ! bond detection
  if(fstr(ftmp,'-bondfactor')) bondf=s2r(arg(i+1))
  if(fstr(ftmp,'-ext')) then
    traj_ext=.true.
    ecommand=(arg(i+1))
  endif
  if(fstr(ftmp,'-geom')) then
    geom_lib=.true.
    call cstring(arg(i+1),n_options)
    print*,'Found ',n_options,' options'
    if(n_options>5) STOP "wrong option string! Re-consider input!"
    allocate(options(n_options))
    do iword=1,n_options
      call charXsplit(arg(i+1),options(iword),iword)
      print*,trim(options(iword))
    enddo

  endif
 enddo


endif


end subroutine
