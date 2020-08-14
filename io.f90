!************************************************************
!* reads a turbomole (bohr) or xmol (angst)rom file.        *
!* Tests if xmol starts with "number of atoms + blank" or   *
!* directly with the coordinates.                           *
!************************************************************
!! THIS ONE READS IN ANGSTROM!
subroutine tmolrd(infile,c_nat,xyz,iat)
use parm, only: nat
use constant, only: au2ang
implicit none
real(8),intent(out)    :: xyz(3,nat)
integer(4), intent(out) :: iat(nat)
character(2) cc,ff
character(80)  atmp
character(*) infile
real(8) txyz(3,nat),xx(5)
real(8) bohr
integer tiat(nat),nn,tifrez(nat),istat,iff
integer i,j,ifrez(nat)
logical da,c_nat,echo
bohr=au2ang
i=0
tifrez=0
iff=0
echo=.true.
inquire(file=infile,exist=da)
select case (da)
case (.true.)
      if(.not.c_nat) write(*,'('' reading...'',$)')

open(unit=33,file=infile)
! test for tmol or txyz file

 read(33,'(a)') atmp ! $coord
rewind(33)
if(index(atmp,'$coord').ne.0) then

 ! count number of atoms
 do while (da)
  read(33,'(a)',end=100) atmp ! $coord
   if(index(atmp,'$coord').ne.0) cycle
   if(index(atmp,'$').ne.0) exit
   i=i+1
  enddo
 nat=i
 100 continue
 if(c_nat) then
  close(33)
  return  ! just return number of atoms
 endif
 rewind(unit=33)

 ! read TMOL file
 read(33,*) atmp ! $coord
 do j=1,nat
!    read(33,'(a)') atmp ! $coord
!    backspace(33)
!    if(index(atmp,' f ').ne.0) then
!    read(33,*) txyz(1,j),txyz(2,j),txyz(3,j),cc,ff
!     tifrez(j)=1
!    iff=iff+1
!    else ! default
     read(33,*) txyz(1,j),txyz(2,j),txyz(3,j),cc
!   endif
   call elem(cc,tiat(j))
   txyz(1:3,j)=txyz(1:3,j)*au2ang
  enddo
 if(echo) write(*,*) ' Turbomole file [bohr] :  ', trim(infile)

 close(33)

else ! txyz file
       read(33,'(a)',end=101) atmp
! check for first two lines
       call readl(atmp,xx,nn)
        if(nn.gt.1) then   ! more than one argument found, assuming they are coords
           do
            nat=nat+1
            read(33,'(a)',end=123) atmp
           enddo
            if(c_nat) then
              close(33)
              return  ! just return number of atoms
            endif
          else
            nat=idint(xx(1))
            if(c_nat) then
             close(33)
             return  ! just return number of atoms
            endif
           read(33,'(a)',end=101) atmp  !titel line
        endif
 123   if(nn.gt.1) rewind(3)
       do i=1,nat
            read(33,'(a)') atmp
            call readl(atmp,xx,nn)
            call elem(atmp,tiat(i))
            txyz(1:3,i)=xx(1:3)*1d0
!             txyz(1:3,i)=xx(1:3)*1d0/bohr
       enddo
 101  close(33)
      if(echo) write(*,'(5x,'' XYZ file [angst]: '',a)')  trim(infile)
      endif

if(maxval(tifrez,nat).eq.1) then
 if(echo) then
  write(*,'(a,x,I4,x,a)') '  found ',iff,' frozen cart. coordinates'
  if(iff.lt.50) then ! dont spam the output to much ...
   write(*,'(a,$)') '  atom nr: '
   do i=1,nat
     if(tifrez(i)==1) write(*,'(x,I2,$)') i
   enddo
   print*,''
  endif
 endif
endif

case (.false.)
  write(*,*) ' no input file <',trim(infile) ,'> found !! '
end select

do i=1,nat
xyz(1:3,i)=txyz(1:3,i)
iat(i)=tiat(i)
ifrez(i)=tifrez(i)
enddo
return

end subroutine


      SUBROUTINE READL(A1,X,N)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*(*) A1
      DIMENSION X(*)
      I=0
      IS=1
  10  I=I+1
      X(I)=READAA(A1,IS,IB,IE)
      IF(IB.GT.0 .AND. IE.GT.0) THEN
                                IS=IE
                                GOTO 10
      ENDIF
      N=I-1
      RETURN
      END


      FUNCTION READAA(A,ISTART,IEND,IEND2)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 READAA
      CHARACTER*(*) A
      NINE=ICHAR('9')
      IZERO=ICHAR('0')
      MINUS=ICHAR('-')
      IDOT=ICHAR('.')
      ND=ICHAR('D')
      NE=ICHAR('E')
      IBL=ICHAR(' ')
      IEND=0
      IEND2=0
      IDIG=0
      C1=0
      C2=0
      ONE=1.D0
      X = 1.D0
      NL=LEN(A)
      DO 10 J=ISTART,NL-1
         N=ICHAR(A(J:J))
         M=ICHAR(A(J+1:J+1))
         IF(N.LE.NINE.AND.N.GE.IZERO .OR.N.EQ.IDOT)GOTO 20
         IF(N.EQ.MINUS.AND.(M.LE.NINE.AND.M.GE.IZERO.OR. M.EQ.IDOT)) GOTO 20

   10 CONTINUE
      READAA=0.D0
      RETURN
   20 CONTINUE
      IEND=J
      DO 30 I=J,NL
         N=ICHAR(A(I:I))
         IF(N.LE.NINE.AND.N.GE.IZERO) THEN
            IDIG=IDIG+1
            IF (IDIG.GT.10) GOTO 60
            C1=C1*10+N-IZERO
         ELSEIF(N.EQ.MINUS.AND.I.EQ.J) THEN
            ONE=-1.D0
         ELSEIF(N.EQ.IDOT) THEN
            GOTO 40
         ELSE
            GOTO 60
         ENDIF
   30 CONTINUE
   40 CONTINUE
      IDIG=0
      DO 50 II=I+1,NL
         N=ICHAR(A(II:II))
         IF(N.LE.NINE.AND.N.GE.IZERO) THEN
            IDIG=IDIG+1
            IF (IDIG.GT.10) GOTO 60
            C2=C2*10+N-IZERO
            X = X /10
         ELSEIF(N.EQ.MINUS.AND.II.EQ.I) THEN
            X=-X
         ELSE
            GOTO 60
         ENDIF
   50 CONTINUE
!C
!C PUT THE PIECES TOGETHER
!C
   60 CONTINUE
      READAA= ONE * ( C1 + C2 * X)
      DO 55 J=IEND,NL
         N=ICHAR(A(J:J))
         IEND2=J
         IF(N.EQ.IBL)RETURN
   55 IF(N.EQ.ND .OR. N.EQ.NE)GOTO 57
      RETURN

   57 C1=0.0D0
      ONE=1.0D0
      DO 31 I=J+1,NL
         N=ICHAR(A(I:I))
         IEND2=I
         IF(N.EQ.IBL)GOTO 70
         IF(N.LE.NINE.AND.N.GE.IZERO) C1=C1*10.0D0+N-IZERO
         IF(N.EQ.MINUS)ONE=-1.0D0
   31 CONTINUE
   61 CONTINUE
   70 READAA=READAA*10**(ONE*C1)
      RETURN
      END
 SUBROUTINE ELEM(KEY1, NAT)
 IMPLICIT DOUBLE PRECISION (A-H,O-Z)
 CHARACTER*(*) KEY1
 CHARACTER*2 ELEMNT(94),E

 DATA ELEMNT/'h ','he',                                      &
 'li','be','b ','c ','n ','o ','f ','ne',                    &
 'na','mg','al','si','p ','s ','cl','ar',                    &
 'k ','ca','sc','ti','v ','cr','mn','fe','co','ni','cu',     &
 'zn','ga','ge','as','se','br','kr',                         &
 'rb','sr','y ','zr','nb','mo','tc','ru','rh','pd','ag',     &
 'cd','in','sn','sb','te','i ','xe',                         &
 'cs','ba','la','ce','pr','nd','pm','sm','eu','gd','tb','dy',&
 'ho','er','tm','yb','lu','hf','ta','w ','re','os','ir','pt',&
 'au','hg','tl','pb','bi','po','at','rn',                    &
 'fr','ra','ac','th','pa','u ','np','pu'/

 nat=0
 e='  '
 k=1
 DO J=1,len(key1)
    if (k.gt.2)exit
    N=ICHAR(key1(J:J))
    if(n.ge.ichar('A') .and. n.le.ichar('Z') )then
       e(k:k)=char(n+ICHAR('a')-ICHAR('A'))
       k=k+1
    endif
    if(n.ge.ichar('a') .and. n.le.ichar('z') )then
       e(k:k)=key1(j:j)
       k=k+1
    endif
 enddo
 DO I=1,107
         if(e.eq.elemnt(i))then
            NAT=I
            RETURN
         ENDIF
      ENDDO

      end

FUNCTION ESYM(I)
CHARACTER*2 ESYM
CHARACTER*2 ELEMNT(94)
DATA ELEMNT/'h ','he',                                           &
  'li','be','b ','c ','n ','o ','f ','ne',                       &
  'na','mg','al','si','p ','s ','cl','ar',                       &
  'k ','ca','sc','ti','v ','cr','mn','fe','co','ni','cu',        &
  'zn','ga','ge','as','se','br','kr',                            &
  'rb','sr','y ','zr','nb','mo','tc','ru','rh','pd','ag',        &
  'cd','in','sn','sb','te','i ','xe',                            &
  'cs','ba','la','ce','pr','nd','pm','sm','eu','gd','tb','dy',   &
  'ho','er','tm','yb','lu','hf','ta','w ','re','os','ir','pt',   &
  'au','hg','tl','pb','bi','po','at','rn',                       &
  'fr','ra','ac','th','pa','u ','np','pu'/
  ESYM=ELEMNT(I)
  RETURN
END

! read xyz-trajectory 
subroutine read_trajxyz(infile,nat,iat,mxyz,nmol,do_allocate)
implicit none
real(8) mxyz(3,nat,nmol)
logical do_allocate,fstr,debug
character(*) infile
integer i,j,nat,nmol,io,iat(nat),ierr
character(200) aa,cnat
character(2) el

open(newunit=io,file=infile,iostat=ierr)
if(ierr/=0) stop 'trajectory file not found'

debug=.true.
if(do_allocate) then
  if(debug) print*,'checking: ',trim(infile)
  read(io,'(a)') cnat
  read(cnat,*) nat
  if(debug) print'(a,I6)', ' found nat  : ',nat
  rewind(io)
  nmol=0
  do
   read(io,'(a)',end=666) aa
   if(fstr(aa,cnat)) then
     nmol=nmol+1 ! count lines with number of atoms
   endif
  enddo
  666 continue
  if(debug) print'(a,I6)',' found nmol : ',nmol
  close(io)
return
endif

! reading EVERYTHING into memory
do i=1,nmol
   read(io,'(a)',end=667) aa           ! nat
   read(io,'(a)') aa                   ! title
   do j=1,nat
     read(io,*) el,mxyz(1:3,j,i) ! xyz, element
     call elem(el,iat(j))
   enddo
enddo
667 continue
close(io)

end subroutine


!******************
!* write xyz      *
!******************
subroutine wrxyz(iat,nat,xyz,infile,do_print)
!use parm
implicit none
integer i,j,k,l,nat,iat(nat)
real(8) xyz(3,*)
integer ierr
character(*) infile
character(120) message
character(2) esym
real(8) f
logical do_print
if(do_print) print*,'writing ',trim(infile)
!f=0.5291770d0
f=1d0
open(unit=55,file=infile,iostat=ierr,status='replace')
message='cannot write '//infile
if(ierr.ne.0) call error(6,trim(message))
write(55,'(I5)') nat
write(55,*)
!write(55,'(2F16.8)')energy,gnorm
do i=1,nat
 write(55,'(a2,5x,3(F18.14,3x))') esym(iat(i)), xyz(1,i)*f,xyz(2,i)*f,xyz(3,i)*f
enddo
close(55)
end subroutine wrxyz



subroutine error(io,aa)
implicit none
integer io
character(*) aa


write(io,'(a)')
write(io,'(3x,a)')'************  ERROR  *********************'
write(io,'(4x,a)') trim(aa)
write(io,'(3x,a)')'******************************************'
write(io,'(a)')
stop

end subroutine error

