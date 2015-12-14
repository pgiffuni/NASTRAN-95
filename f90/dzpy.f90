SUBROUTINE dzpy(kb,ks,ls,i,j1,j2,nyflag,          sgr,cgr,  &
        fmach,arb,nbea,dt)
!   ***   GENERATES ROWS OF THE SUBMATRICES  DZP, DZZ  AND  DZY
!         USING  SUBROUTINE  SUBB
 
 INTEGER, INTENT(IN OUT)                  :: kb
 INTEGER, INTENT(IN OUT)                  :: ks
 INTEGER, INTENT(OUT)                     :: ls
 INTEGER, INTENT(IN OUT)                  :: i
 INTEGER, INTENT(IN)                      :: j1
 INTEGER, INTENT(IN)                      :: j2
 INTEGER, INTENT(IN)                      :: nyflag
 REAL, INTENT(IN OUT)                     :: sgr
 REAL, INTENT(IN OUT)                     :: cgr
 REAL, INTENT(IN)                         :: fmach
 REAL, INTENT(IN)                         :: arb(1)
 INTEGER, INTENT(IN)                      :: nbea(1)
 COMPLEX, INTENT(OUT)                     :: dt(1)
 COMPLEX :: sum
 
 COMMON /dlbdy/ nj1,nk1,np,nb,ntp,nbz,nby,ntz,nty,nt0,ntzs,ntys,  &
     inc,ins,inb,inas,izin,iyin,inbea1,inbea2,insbea,izb,iyb,  &
     iavr,iarb,infl,ixle,ixte,int121,int122,izs,iys,ics,iee,isg,  &
     icg,ixij,ix,idelx,ixic,ixlam,ia0,ixis1,ixis2,ia0p,iria  &
     ,inasb,ifla1,ifla2,ith1a,ith2a, ecore,next,scr1,scr2,scr3,scr4,scr5
 COMMON /zzzzzz / z(1)
 
 ndy  = 0
 nyfl = nyflag
 pi   = 3.1415926
 eps  = 0.00001
 beta = SQRT(1.0-fmach**2)
 jz   = 0
 lb   = 1
 jb   = 0
 ar   = arb(lb)
!  LS IS THE INDEX OF THE  Y  AND  Z  COORDINATES OF SENDING POINT  J --
!  LS RUNS FROM NSTRIP+1  THROUGH  NSTRIP+NBZ
 DO    j=j1,j2
   jb   = jb+1
!  JB  IS THE BODY-ELEMENT NUMBER IN BODY  LB  --  JB  RUNS FROM  1
!  THROUGH  NTZ
   jz   = jz+1
!  JZ  RUNS FROM 1  THROUGH  NBE-SUB-LB
   CALL subb(kb,ks,i,j,jb,lb,ls,ndy,nyfl,          pi,eps,  &
       sgr,cgr,ar,beta,sum,z(iria),z(idelx),z(iyb),z(izb),z(iys), z(izs),z(ix))
   dt(j)= sum
   IF  (jz == nbea(lb))    GO TO  20
   CYCLE
   20 CONTINUE
   jz   = 0
   lb   = lb+1
   ls   = ls+1
   ar   = arb(lb)
 END DO
 RETURN
END SUBROUTINE dzpy
