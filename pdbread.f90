! reads a PDB and guesses the atom types


integer function natPDB(filename)
implicit none
integer io,i
character(*) filename
character(200) atmp
open(newunit=io,file=filename)
i=0
do
 read (io,'(a)',end=777) atmp
 if(index(atmp,'ATOM').ne.0) i=i+1
enddo
777 continue
natPDB=i
close(io)
end function


character(2) function id_name(string)
implicit none
integer i
character(*) string
character(2) a
character(2) el(95)

data el/'h ','he',                                             &
  'li','be','b ','c ','n ','o ','f ','ne',                       &
  'na','mg','al','si','p ','s ','cl','ar',                       &
  'k ','ca','sc','ti','v ','cr','mn','fe','co','ni','cu',        &
  'zn','ga','ge','as','se','br','kr',                            &
  'rb','sr','y ','zr','nb','mo','tc','ru','rh','pd','ag',        &
  'cd','in','sn','sb','te','i ','xe',                            &
  'cs','ba','la','ce','pr','nd','pm','sm','eu','gd','tb','dy',   &
  'ho','er','tm','yb','lu','hf','ta','w ','re','os','ir','pt',   &
  'au','hg','tl','pb','bi','po','at','rn',                       &
  'fr','ra','ac','th','pa','u ','np','pu','xx'/ ! elem(95) being a dummy

string=adjustl(string)
call lower_case(string)
do i=1,10
 a=string(1:1)//' '
 if(a==el(i)) print*,a,'hiut'
enddo

do i=1,10
 a=string(1:2)
 if(a==el(i)) print*,a,'hiut'
enddo


!do i=1,len_trim(string)
! if(any(string(1:i)==el)) print*, string
! if(any(el==string(1:i))) print*, string
!enddo
end function

subroutine readPDB(filename,nat,xyz,iat)
implicit none
integer i,nat,iat(nat),io
real(8) xyz(3,nat),s2r
character(*) filename
character(200) atmp

character(4) aname
character(2) id_name
!character(3) resname

!write(*,'('' reading...'',$)')

open(newunit=io,file=filename)
do
 read(io,'(a)',end=666) atmp
 if(index(atmp,'ATOM').ne.0) then
  i=i+1
  aname=id_name(atmp(13:16))
  print*, aname
  call elem(aname,iat(i))
  xyz(1,i)=s2r(atmp(31:38))
  xyz(2,i)=s2r(atmp(39:46))
  xyz(3,i)=s2r(atmp(47:54))
 endif
enddo
666 continue
close(io)
write(*,'(5x,'' PDB file [angst]: '',a)')  trim(filename)
return
end subroutine

! OLUMNS        DATA  TYPE    FIELD        DEFINITION
! -------------------------------------------------------------------------------------
!  1 -  6        Record name   "ATOM  "
!  7 - 11        Integer       serial       Atom  serial number.
! 13 - 16        Atom          name         Atom name.
! 17             Character     altLoc       Alternate location indicator.
! 18 - 20        Residue name  resName      Residue name.
! 22             Character     chainID      Chain identifier.
! 23 - 26        Integer       resSeq       Residue sequence number.
! 27             AChar         iCode        Code for insertion of residues.
! 31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
! 39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
! 47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
! 55 - 60        Real(6.2)     occupancy    Occupancy.
! 61 - 66        Real(6.2)     tempFactor   Temperature  factor.
! 77 - 78        LString(2)    element      Element symbol, right-justified.
! 79 - 80        LString(2)    charge       Charge  on the atom.


