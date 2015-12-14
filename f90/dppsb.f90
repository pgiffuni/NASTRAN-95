SUBROUTINE dppsb(  ks,i,j1,j2,sgr,cgr,           ys,zs,nbaray,  &
        ncaray,dt,z)
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
 INTEGER, INTENT(IN)                      :: z(1)
 
 
 COMPLEX :: sum
 COMMON /dlbdy/ nj1,nk1,np,nb,ntp,nbz,nby,ntz,nty,nt0,ntzs,ntys,  &
     inc,ins,inb,inas,izin,iyin,inbea1,inbea2,insbea,izb,iyb,  &
     iavr,iarb,infl,ixle,ixte,int121,int122,izs,iys,ics,iee,isg,  &
     icg,ixij,ix,idelx,ixic,ixlam,ia0,ixis1,ixis2,ia0p,iria  &
     ,inasb,ifla1,ifla2,ith1a,ith2a, ecore,next,scr1,scr2,scr3,scr4,scr5
 
 l    = 1
!  L IS THE PANEL NUMBER ASSOCIATED WITH SENDING   POINT  J
 ls   = 1
 lsp = 0
!  LS IS THE STRIP NUMBER ASSOCIATED WITH SENDING   POINT  J
 nbxs = nbaray(l)
 nc1  = ncaray(l)
 nbcum= nc1
 yrec = ys(ks)
 zrec = zs(ks)
 DO     j=j1,j2
   CALL subpb(i,l,ls,j,sgr,cgr,yrec,zrec,sum,z(ixic),z(idelx),z(iee)  &
       ,z(ixlam),z(isg),z(icg),z(iys),z(izs),z(inas),z(inasb+lsp),  &
       z(iavr),z(izb),z(iyb),z(iarb),z(ixle),z(ixte),z(ix),nb)
   dt(j)= sum
   IF (j == j2)  CYCLE
   IF (j < nbxs)   GO TO  10
   lsp = lsp + z(inas+l-1)
   l    = l+1
   nc1  = ncaray(l)
   nbxs = nbaray(l)
   10 CONTINUE
   IF (j < nbcum)  CYCLE
   ls   = ls+1
   nbcum= nbcum+nc1
 END DO
 RETURN
END SUBROUTINE dppsb
