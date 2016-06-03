!*****************
!* print timings *
!*****************
subroutine prtim(io,tt,is,string)
integer io
real*8 tt,t,sec
integer day,hour,min
character*(*) string
character*1 is

t=tt
day=idint(t/86400)
t=t-day*86400
hour=idint(t/3600)
t=t-hour*3600
min=idint(t/60)
t=t-60*min
sec=t

if(day.ne.0)then
   if(is=='w')write(io,2)trim(string),day,hour,min,sec
   if(is=='t')write(io,22)trim(string),day,hour,min,sec
   return
endif
if(hour.ne.0)then
   if(is=='w')write(io,3)trim(string),hour,min,sec
   if(is=='t')write(io,33)trim(string),hour,min,sec
   return
endif
if(min .ne.0)then
   if(is=='w')write(io,4)trim(string),min,sec
   if(is=='t')write(io,44)trim(string),min,sec
   return
endif
  if(is=='w')write(io,5)trim(string),sec
  if(is=='t')write(io,55)trim(string),sec
return

 1    format('======================================')
 2    format('wall-time for ',a,2x,i3,' d',i3,' h',i3,' m',f5.1,' s')
 3    format('wall-time for ',a,2x,i3,' h',i3,' m',f5.1,' s')
 4    format('wall-time for ',a,2x,i3,' m',f5.1,' s')
 5    format('wall-time for ',a,2x,f5.1,' s')

 22    format('cpu-time for ',a,2x,i3,' d',i3,' h',i3,' m',f5.1,' s')
 33    format('cpu-time for ',a,2x,i3,' h',i3,' m',f5.1,' s')
 44    format('cpu-time for ',a,2x,i3,' m',f5.1,' s')
 55    format('cpu-time for ',a,2x,f5.1,' s')

return
end

!SUBROUTINE PRINTS MATRIX R,WHICH IS SUPPOSED
!TO HAVE DIMENSION N,M  WHEN M IS NONZERO AND
!((N+1)*N)/2 WHEN M IS ZERO
SUBROUTINE PRMATI(IUOUT,R,N,M,HEAD)
integer R
CHARACTER*(*) HEAD
DIMENSION R(*)

WRITE(IUOUT,1001) HEAD
NKPB=6
IF(M)10,10,80

 10 CONTINUE
    IBL=N/NKPB
    IR=N-IBL*NKPB
    J1=1
    K1S=1
    KD=0
    IF(IBL.EQ.0) GO TO 50
    J2=NKPB
    DO 40 I=1,IBL
    WRITE(IUOUT,1002)(J,J=J1,J2)
    K1=K1S
    K2=K1
    KK=0
    DO 20 J=J1,J2
    WRITE(IUOUT,1003)J,(R(K),K=K1,K2)
    KK=KK+1
    K1=K1+KD+KK
 20 K2=K1+KK
    J1=J1+NKPB
    IF(J1.GT.N) RETURN
    J2=J2+NKPB
    K2=K1-1
    K1=K2+1
    K2=K1+(NKPB-1)
    K1S=K2+1
    KK=KD+NKPB
    DO 30 J=J1,N
    WRITE(IUOUT,1003)J,(R(K),K=K1,K2)
    KK=KK+1
    K1=K1+KK
 30 K2=K2+KK
 40 KD=KD+NKPB
 50 IF(IR.EQ.0) GO TO 70
    K1=K1S
    J2=J1+IR-1
    KK=0
    K2=K1
    WRITE(IUOUT,1002)(J,J=J1,J2)
    WRITE(IUOUT,1003)
    DO 60 J=J1,J2
    WRITE(IUOUT,1003)J,(R(K),K=K1,K2)
    KK=KK+1
    K1=K1+KD+KK
 60 K2=K1+KK
 70 RETURN
 80 IBL=M/NKPB
    IR=M-IBL*NKPB
    I2=0
    K2=0
    IF(IBL.EQ.0) GO TO 100
    DO 90 I=1,IBL
    I1=(I-1)*N*NKPB+1
    I2=I1+(NKPB-1)*N
    K1=K2+1
    K2=K1+(NKPB-1)
    WRITE(IUOUT,1002)(K,K=K1,K2)
    DO 90 J=1,N
    WRITE(IUOUT,1003)J,(R(IJ),IJ=I1,I2,N)
    I1=I1+1
 90 I2=I1+(NKPB-1)*N
100 IF(IR.EQ.0) GO TO 120
    I1=IBL*N*NKPB+1
    I2=I1+(IR-1)*N
    K1=K2+1
    K2=M
    WRITE(IUOUT,1002)(K,K=K1,K2)
    WRITE(IUOUT,1003)
    DO 110 J=1,N
    WRITE(IUOUT,1003)J,(R(IJ),IJ=I1,I2,N)
    I1=I1+1
    I2=I1+(IR-1)*N
110 CONTINUE
120 WRITE(IUOUT,1003)
   RETURN
1001 FORMAT(/'->>',2X,A)
1002 FORMAT(/,' ',4X,6(1X,I4,1X),/)
1003 FORMAT(' ',I4,6(1x,I4,1x))
END

