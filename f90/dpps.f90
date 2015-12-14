SUBROUTINE dpps(ks,i,j1,j2,sgr,cgr,ys,zs,nbaray,ncaray,dt,work)
!   ***   GENERATES ROWS OF THE  DPP  SUBMATRIX USING
!         SUBROUTINE  SUBP
 
 INTEGER, INTENT(IN OUT)                  :: ks
 INTEGER, INTENT(IN OUT)                  :: i
 INTEGER, INTENT(IN)                      :: j1
 INTEGER, INTENT(IN)                      :: j2
 REAL, INTENT(IN OUT)                     :: sgr
 REAL, INTENT(IN OUT)                     :: cgr
 REAL, INTENT(IN)                         :: ys(1)
 REAL, INTENT(IN)                         :: zs(1)
 INTEGER, INTENT(IN)                      :: nbaray(1)
 INTEGER, INTENT(IN)                      :: ncaray(1)
 COMPLEX, INTENT(OUT)                     :: dt(1)
 REAL, INTENT(IN OUT)                     :: work(1)
 
 COMPLEX :: sum
 COMMON /dlcom / np,nstrip,ntp,f,njj,next,length,  &
     inc,inb,iys,izs,iee,isg,icg, ixic,idelx,ixlam,idt,  &
     icore
 
 l    = 1
!  L IS THE PANEL NUMBER ASSOCIATED WITH SENDING   POINT  J
 ls   = 1
!  LS IS THE STRIP NUMBER ASSOCIATED WITH SENDING   POINT  J
 nbxs = nbaray(l)
 nc1  = ncaray(l)
 nbcum= nc1
 yrec = ys(ks)
 zrec = zs(ks)
 DO     j=j1,j2
   CALL subp(i,l,ls,j,sgr,cgr,yrec,zrec,sum,  &
       work(ixic),work(idelx),work(iee),work(ixlam), work(isg),work(icg),ys,zs)
   dt(j)= sum
   IF (j == j2)  CYCLE
   IF (j < nbxs)   GO TO  10
   l    = l+1
   nc1  = ncaray(l)
   nbxs = nbaray(l)
   10 CONTINUE
   IF (j < nbcum)  CYCLE
   ls   = ls+1
   nbcum= nbcum+nc1
 END DO
 RETURN
END SUBROUTINE dpps
