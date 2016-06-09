Pure Function to_upper (str) Result (string)

!   ==============================
!   Changes a string to upper case
!   ==============================

    Implicit None
    Character(*), Intent(In) :: str
    Character(LEN(str))      :: string

    Integer :: ic, i

    Character(26), Parameter :: cap = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    Character(26), Parameter :: low = 'abcdefghijklmnopqrstuvwxyz'

!   Capitalize each letter if it is lowecase
    string = str
    do i = 1, LEN_TRIM(str)
        ic = INDEX(low, str(i:i))
        if (ic > 0) string(i:i) = cap(ic:ic)
    end do
End Function to_upper





!********************************
!* convert a word to lower case *
!********************************
      subroutine lower_case(word)
      character (len=*) , intent(in out) :: word
      integer :: i,ic,nlen
      nlen = len(word)
      do i=1,nlen
      ic = ichar(word(i:i))
      if (ic >= 65 .and. ic < 90) word(i:i) = char(ic+32)
      end do
      end subroutine lower_case


!********************************************
!* split a string s into n separate words w *
!********************************************
subroutine charsplit(s,n,w)
implicit none
integer i,n,k
character*80, intent(in) :: s
character*80, intent(out) :: w(n)
character*80   a,aa

aa=adjustl(s)
do i=1,n
  a=aa
  k=index(a,' ')
  w(i)=trim(a(:(k-1)))
  aa=adjustl(trim(a((k+1):)))
enddo
return
end subroutine



!***************************************************
!* split a string s iand the return the x'ths word *
!***************************************************
subroutine charXsplit(s,wx,x)
implicit none
integer i,n,k,x
character*80, intent(in) :: s
character*80, intent(out) ::wx
character*80   w(20)
character*80   a,aa

aa=adjustl(s)
i=0
do while (index(aa,' ').ne.1)
  i=i+1
  a=aa
  k=index(a,' ')
  w(i)=trim(a(:(k-1)))
  aa=adjustl(trim(a((k+1):)))
!  print*,'AA',trim(aa),index(aa,' ')
  if(i.gt.50) stop '!* string split error: subroutine charXsplit *!'
enddo
wx=w(x)
!print*, trim(wx)
return
end subroutine


!***********************************
!*  count words n in string s      *
!***********************************
subroutine cstring(s,n)
implicit none
integer i,n,k,x
character(*), intent(in) :: s
character*80   a,aa
aa=adjustl(s)
i=0
do while (index(aa,' ').ne.1)
  i=i+1
  a=aa
  k=index(a,' ')
!  w(i)=trim(a(:(k-1)))
  aa=adjustl(trim(a((k+1):)))
!  print*,'AA',trim(aa),index(aa,' ')
  if(i.gt.50) stop '!* string split error: subroutine charXsplit *!'
enddo
n=i
return
end subroutine

integer function nwords(s)
implicit none
integer i,k,x
character(*), intent(in) :: s
character*80   a,aa
aa=adjustl(s)
i=0
do while (index(aa,' ').ne.1)
  i=i+1
  a=aa
  k=index(a,' ')
!  w(i)=trim(a(:(k-1)))
  aa=adjustl(trim(a((k+1):)))
!  print*,'AA',trim(aa),index(aa,' ')
  if(i.gt.50) stop '!* string split error: subroutine charXsplit *!'
enddo
nwords=i
return
end function



! String to integer
integer pure function s2i(a)
implicit none
character(*), intent(in):: a
!read(a,'(i)') s2i
read(a,*) s2i
return
end function


! String to real(8)
real(8) pure function s2r(a)
implicit none
character(*), intent(in):: a
!read(a,'(f)') s2r
read(a,*) s2r
return
end function


! read a string containing either a single (atom) number or a range or number, eg.: "5,1-3, 7-19"
! ignores spaces
! string : input string to read from
! al     : integer vector containing the atom numbers
! ial    : number of records/atom numbers
subroutine atlist(string,al,ial)
implicit none
character(*) string
character(120) aa,xx,bb
character(1) a,x
character(2) b
character(10) records(100)
logical dash,r
integer sl,sl2,c,k, i,j,s2i,ns,ne,list(10000)
integer al(*),id(120),ii,ial

dash=.false.
r=.false.
records=''
list=0
! removes spaces in string and return aa clean of spaces and resulting character length sl
call rmspace(string,aa,sl)
!write(*,'(a)') trim(aa(1:sl))
id=0
c=1
! separator positions
do  i=1,sl
a=aa(i:i)
if(a==',') then
 id(c)=i
 c=c+1
endif
enddo


!*****************************
! assign records
c=1
k=count(id.ne.0)  ! = number of separators
do i=1,k
 records(i)=aa(c:id(i)-1)
!print*,records(i),c,id(i)-1
c=id(i)+1
enddo
 ! last record
records(k+1)=aa(c:sl)
!print*,records(i),c,sl
!*****************************

c=0
! now read records into atom numbers
do i=1,k+1
ii=index(records(i),'-')
sl=len_trim(records(i))
xx=trim(records(i))
if(ii.ne.0) then
! print*,trim(xx(1:ii-1))
! print*,trim(xx(ii+1:sl))
 read(xx(1:ii-1),*) ns 
 read(xx(ii+1:sl),*) ne 
   do k=ns,ne
     c=c+1
     list(c)=k
   enddo
else
 c=c+1
 read(xx(1:sl),*) list(c)
endif
enddo
ial=count(list.ne.0)
!print*,'Atom number list debug:'
!print*,list(1:count(list.ne.0))

if(c.ne.count(list.ne.0)) stop ' error in records assignments'
al(1:ial)=list(1:ial)

! sort?


end subroutine


subroutine rmspace(aa,sout,len)
implicit none
character(*) aa
character(120) xx,sout
character(1) a
integer i,ic,sl,c,len

c=0
sl=len_trim(adjustr(aa))
do i=1,sl
  a=aa(i:i)
  ic=iachar(a)
  select case(ic)
    case(0:32)  ! space, tab, or control character
     cycle
    case(33:)
     c=c+1
     sout(c:c)=a
  end select
enddo
!print*,'c',c
!write(*,'(a)') sout(1:c)
len=c
end 


! return true if characters 'c' are found in input string 'str'
logical function fstr(str,c)
character(*), intent(in):: str,c
fstr=.false.
fstr=index(str,c).ne.0
end function



subroutine rm_substr(instr,substr,outstr)
character(*) instr,substr,outstr
integer i,s(2)

s(1)=index(instr,substr,.false.)
outstr=instr(1:s(1)-1)
!print*,instr(1:s(1)-1),' ',substr
!print*,s(1)
end subroutine

