SUBROUTINE gend(ncaray,nbaray,ys,zs,sg,cg,dt,work,matout)
!  GENERATE THE INFLUENCE COEFFICIENT MATRIX ADPP
 
 INTEGER, INTENT(IN)                      :: ncaray(1)
 INTEGER, INTENT(IN)                      :: nbaray(1)
 REAL, INTENT(IN OUT)                     :: ys(1)
 REAL, INTENT(IN OUT)                     :: zs(1)
 REAL, INTENT(IN)                         :: sg(1)
 REAL, INTENT(IN)                         :: cg(1)
 COMPLEX, INTENT(OUT)                     :: dt(1)
 REAL, INTENT(IN OUT)                     :: work(1)
 INTEGER, INTENT(IN OUT)                  :: matout
 
 
 COMMON /amgmn / mcb(7),nrow,nd,NE,refc,fmach,rfk
 COMMON /dlcom / np,nstrip,ntp,f,njj,next,length,  &
     inc,inb,iys,izs,iee,isg,icg, ixic,idelx,ixlam,idt,  &
     icore
 
 i1 = 1
 i2 = ntp
 j1 = 1
 j2 = ntp
 
!     POSITION IN DT TO START OF THIS PART OF MATRIX
 
 idtpt = i1 + nrow
 DO  i = i1,njj
   dt(i) = (0.0,0.0)
 END DO
!     DPP LOOP
 k = 1
!     K IS THE PANEL NUMBER
 ks = 1
!     KS IS THE STRIP NUMBER
 nbxr = ncaray(k)
 DO  i = i1,i2
   sgr = sg(ks)
   cgr = cg(ks)
   CALL dpps(ks,i,j1,j2,sgr,cgr,ys,zs,nbaray,ncaray,dt(idtpt),work)
   CALL pack(dt,matout,mcb)
   IF(i == i2) CYCLE
   IF(i == nbaray(k)) k=k+1
   IF(i == nbxr) GO TO 50
   CYCLE
   50 CONTINUE
   ks = ks +1
   nbxr = nbxr + ncaray(k)
 END DO
 RETURN
END SUBROUTINE gend
