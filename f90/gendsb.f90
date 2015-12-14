SUBROUTINE gendsb(ncaray,nbaray,sg,cg,nfl,nbea1,nbea2,ifla1,  &
        ifla2,dt,dpz,dpy)
     
 INTEGER, INTENT(IN)                      :: ncaray(1)
 INTEGER, INTENT(IN)                      :: nbaray(1)
 REAL, INTENT(IN)                         :: sg(1)
 REAL, INTENT(IN)                         :: cg(1)
 INTEGER, INTENT(IN)                      :: nfl(1)
 INTEGER, INTENT(IN)                      :: nbea1(1)
 INTEGER, INTENT(IN)                      :: nbea2(1)
 INTEGER, INTENT(IN)                      :: ifla1(1)
 INTEGER, INTENT(IN)                      :: ifla2(1)
 COMPLEX, INTENT(OUT)                     :: dt(1)
 COMPLEX, INTENT(OUT)                     :: dpz(1)
 COMPLEX, INTENT(OUT)                     :: dpy(1)
 INTEGER :: scr1,scr2,scr3,scr4,scr5,ecore,sysbuf
 INTEGER :: z
 DIMENSION NAME(2)
 
 
 
 COMMON /system/ sysbuf
 COMMON /zzzzzz / z(1)
 COMMON /amgmn / mcb(7),nrow,nd,NE,refc,fmach,rfk,tskj(7),isk,nsk
 COMMON /dlbdy/ nj1,nk1,np,nb,ntp,nbz,nby,ntz,nty,nto,ntzs,ntys,  &
     inc,ins,inb,inas,izin,iyin,inbea1,inbea2,insbea,izb,iyb,  &
     iavr,iarb,infl,ixle,ixte,int121,int122,izs,iys,ics,iee,isg,  &
     icg,ixij,ix,idelx,ixic,ixlam,ia0,ixis1,ixis2,ia0p,iria,  &
     inasb,iflax,ifla ,ith1a,ith2a, ecore,next,scr1,scr2,scr3,scr4,scr5,ntbe
 DATA NAME /4HGEND,4HB   /
!   ***   GENERATES THE INFLUENCE COEFFICIENT MATRIX  DT   USING THE
!         FOLLOWING FOUR SUBROUTINES  --  DPPS, DPZY, DZPY,  AND  DYPZ
 nbox = ntp
 lbo  = 1
 lso  = 1
 jbo  = 1
 kb   = 0
 kt   = 0
 DO  i=1,ntbe
   dpy(i) = (0.0,0.0)
   dt(i) = (0.0,0.0)
 END DO
 nbuf = 4
 IF(ntp == 0) nbuf = nbuf - 1
 IF(ntz == 0) nbuf = nbuf - 1
 IF(nty == 0) nbuf = nbuf - 2
 IF(next + nbuf*sysbuf > ecore) CALL mesage(-8,0,NAME)
 ibuf1 = ecore - sysbuf
 ibuf2 = ibuf1 - sysbuf
 nstrip = 0
 j2 = 0
 i2 = 0
 nyflag = 0
 IF(ntp /= 0) CALL gopen(scr1,z(ibuf1),1)
 IF(ntp == 0) GO TO 201
 i1   = 1
 i2   = ntp
 j1   = 1
 j2   = ntp
!  DPP-LOOP
 k    = 1
!  K IS THE PANEL NUMBER ASSOCIATED WITH RECEIVING POINT  I
 ks   = 1
!  KS IS THE STRIP NUMBER ASSOCIATED WITH RECEIVING POINT  I
 nbxr = ncaray(k)
 DO    i=i1,i2
   sgr  = sg(ks)
   cgr  = cg(ks)
   CALL      dppsb(  ks,i,j1,j2,sgr,cgr,           z(iys),z(izs),  &
       nbaray,ncaray,dt,z(1))
   CALL WRITE(scr1,dt,2*ntp,0)
   IF (i == i2)  CYCLE
   IF (i == nbaray(k))  k=k+1
   IF (i == nbxr)  GO TO   50
   CYCLE
   50 CONTINUE
   ks   = ks+1
   nbxr = nbxr+ncaray(k)
 END DO
 CALL WRITE(scr1,0,0,1)
 nstrip = ks
 nzysv= 0
 DO    j=j1,j2
   dt(j)= (0.0,0.0)
 END DO
 nlt1 = 0
 nlt2 = 0
 IF (ntz == 0)  GO TO  180
 IF(nty /= 0) CALL gopen(scr4,z(ibuf2),1)
 i1   = i2+1
 i2   = i2+ntz
!  DPZ-LOOP    **    ALSO USED FOR GENERATING THE DPY-MATRIX  --  SEE
!  COMMENT IN  DPY-LOOP  BELOW
 80 CONTINUE
 kb   = kb+1
!  KB  IS THE BODY NUMBER ASSOCIATED WITH RECEIVING POINT  I
 iz   = 0
 kt   = kt+1
!  KT  IS THE INDEX OF THE ARRAY OF FIRST-AND-LAST-ELEMENTS FOR THETA-1
 icount = 1
 ifl    = nfl(kb)
 nzykb = nbea2(kb)
 ifirst = ifla1(kt)
 ilast = ifla2(kt)
 DO    i=i1,i2
   DO    j=j1,j2
     dpz(j) = (0.0,0.0)
     dpy(j) = (0.0,0.0)
   END DO
   CALL       dpzy(   kb,iz,i,j1,j2,ifirst,ilast,z(iyb),  &
       z(izb),z(iavr),z(iarb),z(ith1a+nlt1),z(ith2a+nlt2),z(int121),  &
       z(int122),nbaray,ncaray,nzykb,dpz,dpy)
   SELECT CASE ( nzykb )
     CASE (    1)
       GO TO 100
     CASE (    2)
       GO TO 100
     CASE (    3)
       GO TO 110
   END SELECT
   100 CONTINUE
   CALL WRITE(scr1,dpz,2*ntp,0)
   IF (nzykb == 1)  GO TO 120
   110 CONTINUE
   CALL WRITE(scr4,dpy,2*ntp,0)
   120 CONTINUE
   IF (iz == nbea1(kb) )  GO TO 130
   IF (iz == ilast.AND.icount < ifl)  GO TO 160
   CYCLE
   130 CONTINUE
   iz     = 0
   IF (nzysv <= 1.AND.nzykb >= 2)  GO TO  140
   GO TO  150
   140 CONTINUE
   lbo  = kb
   lso  = nstrip+lbo
   jbo  = i-nbea1(kb) -nbox+1
   150 CONTINUE
   nzysv = nzykb
   IF(i == i2) EXIT
   kb     = kb+1
   icount = 0
   ifl    = nfl(kb)
   nzykb = nbea2(kb)
   160 CONTINUE
   kt     = kt+1
   icount = icount+1
   ifirst = ifla1(kt)
   ilast = ifla2(kt)
 END DO
 180 CONTINUE
 IF(i2 == ntbe) GO TO 190
!  DPY-LOOP    **    THIS LOOP IS REDUCED TO SETTING THE CORRECT INDICES
!  AND USING THE  DPZ-LOOP  ABOVE
 IF(ntz == 0) CALL gopen(scr4,z(ibuf2),1)
 i1   = i2+1
 i2 = ntbe
 GO TO 80
 190 CALL WRITE(scr1,0,0,1)
 IF(nty /= 0) CALL WRITE(scr4,0,0,1)
 CALL CLOSE(scr1,1)
 CALL CLOSE(scr4,1)
 i1   = 1
 i2   = ntp
 IF (ntz == 0)  GO TO  250
 CALL gopen(scr2,z(ibuf1),1)
!  DZP-LOOP
 k    = 1
!  K  IS THE PANEL NUMBER ASSOCIATED WITH RECEIVING POINT  I
 ks   = 1
!  KS  IS THE STRIP NUMBER ASSOCIATED WITH RECEIVING POINT  I
 nbxr = ncaray(k)
 kb   = 0
!  HERE  KB=0  SERVES AS A FLAG INDICATING THAT THE RECEIVING POINT   I
!  IS ON A PANEL AND NOT   ON A BODY
 j1   = j2+1
 j2   = j2+ntz
 DO    i=i1,i2
   ls   = nstrip+1
   sgr  = sg(ks)
   cgr  = cg(ks)
   CALL       dzpy(kb,ks,ls,   i,j1,j2,nyflag,          sgr,cgr,  &
       fmach,   z(iarb),z(inbea1),dt)
   CALL WRITE(scr2,dt(j1),2*ntz,0)
   IF (i == i2)  CYCLE
   IF (i == nbaray(k))  k =k +1
   IF (i == nbxr)  GO TO  200
   CYCLE
   200 CONTINUE
   ks   = ks+1
   nbxr = nbxr+ncaray(k)
 END DO
 CALL WRITE(scr2,0,0,1)
 201 CONTINUE
 IF(ntz == 0) GO TO 250
 IF(ntp == 0) CALL gopen(scr2,z(ibuf1),1)
 nyflag = 0
!  DZZ-LOOP    **    ALSO USED FOR GENERATING THE  DZY  MATRIX  --  SEE
!  COMMENT IN  DZY-LOOP  BELOW
 kb   = 1
!  KB  IS THE BODY NUMBER ASSOCIATED WITH RECEIVING POINT  I
 ks   = nstrip+1
 iz   = 0
 i1   = i2+1
 i2   = i2+ntz
 sgr  = 0.0
 cgr  = 1.0
 220 CONTINUE
 ls   = nstrip+1
 lsx  = ls
 DO    i=i1,i2
   ls   = lsx
   iz   = iz+1
!  KS IS THE INDEX OF THE Y  AND  Z  COORDINATES OF RECEIVING POINT I
!  IN THE  DZZ-LOOP  KS RUNS FROM  (NSTRIP+1)  THROUGH  (NSTRIP+NBZ)
!  IN THE  DZY-LOOP  KS  RUNS FROM  (NSTRIP+NB-NBY+1) THROUGH  NSTRIP+NB
   CALL       dzpy(kb,ks,ls,   i,j1,j2,nyflag,          sgr,cgr,  &
       fmach,   z(iarb),z(inbea1),dt)
   CALL WRITE(scr2,dt(j1),2*ntz,0)
   IF (iz == nbea1(kb) )  GO TO  230
   CYCLE
   230 CONTINUE
   iz   = 0
   kb   = kb+1
   ks   = ks+1
 END DO
 CALL WRITE(scr2,0,0,1)
 IF(nty == 0) CALL CLOSE(scr2,1)
 IF (nty == 0)  GO TO 320
 IF (nyflag /= 0)  GO TO  250
!  DZY-LOOP    **    THIS LOOP IS REDUCED TO SETTING THE CORRECT INDICES
!  AND USING THE  DZZ-LOOP  ABOVE
 i1 = ntbe-nty+1
 i2 = ntbe
 nyflag = 1
 kb   = lbo
 ks   = lso
 sgr  =-1.0
 cgr  = 0.0
 GO TO  220
 250 CONTINUE
 CALL CLOSE(scr2,1)
 IF (nty == 0)  GO TO 320
 CALL gopen(scr3,z(ibuf1),1)
 i1   = 1
 i2   = ntp
 j1 = ntbe-nty+1
 j2 = ntbe
 IF(ntp == 0) GO TO 275
!  DYP-LOOP
 k    = 1
 ks   = 1
 kb   = 0
 nbxr = ncaray(k)
 sgr  = sg(ks)
 cgr  = cg(ks)
 DO    i=i1,i2
   CALL       dypz(kb,ks,ls,   i,j1,j2,nyflag,          sgr,cgr,  &
       fmach,   z(iarb),z(inbea1),          lbo,lso,jbo,dt)
   CALL WRITE(scr3,dt(j1),2*nty,0)
   IF (i == nbaray(k))  k=k+1
   IF (i == nbxr)  GO TO  260
   CYCLE
   260 CONTINUE
   ks   = ks+1
   nbxr = nbxr+ncaray(k)
   sgr  = sg(ks)
   cgr  = cg(ks)
 END DO
 CALL WRITE(scr3,0,0,1)
 275 CONTINUE
 nyflag = 0
 iz   = 0
 IF (ntz == 0)  GO TO  310
!  DYZ-LOOP    **    ALSO USED FOR GENERATING THE  DYY  MATRIX  --  SEE
!  COMMENT IN  DYY-LOOP  BELOW
 i1   = i2+1
 i2   = i2+ntz
 ks   = nstrip+1
 kb   = 1
 sgr  = 0.0
 cgr  = 1.0
 280 CONTINUE
 DO    i=i1,i2
   ls   = lso
   iz   = iz+1
   CALL       dypz(kb,ks,ls,   i,j1,j2,nyflag,          sgr,cgr,  &
       fmach,   z(iarb),z(inbea1),          lbo,lso,jbo,dt)
   CALL WRITE(scr3,dt(j1),2*nty,0)
   IF (iz == nbea1(kb) )  GO TO  290
   CYCLE
   290 CONTINUE
   iz   = 0
   kb   = kb+1
   ks   = ks+1
 END DO
 CALL WRITE(scr3,0,0,1)
 310 CONTINUE
 IF (nyflag /= 0)  GO TO 320
!  DYY-LOOP    **    THIS LOOP IS REDUCED TO SETTING THE CORRECT INDICES
!  AND USING THE  DYZ-LOOP  ABOVE
 IF(ntp == 0.AND.ntz == 0) CALL gopen(scr3,z(ibuf1),1)
 i1 = ntbe-nty+1
 i2 = ntbe
 nyflag = 1
 kb   = lbo
 ks   = lso
 sgr  =-1.0
 cgr  = 0.0
 GO TO  280
 320 CONTINUE
 CALL CLOSE(scr3,1)
 
!     BUILD SCR5 WITH GEND PART OF A MATRIX
 
 i1   = 1
 i2   = ntp+ntz
 nyflag = 0
 CALL gopen(scr5,z(ibuf1),1)
 ibuf3 = ibuf2 - sysbuf
 ibuf4 = ibuf3 - sysbuf
 IF(ntz /= 0) CALL gopen(scr2,z(ibuf3),0)
 IF(nty /= 0) CALL gopen(scr3,z(ibuf4),0)
 itape = scr1
 IF(i2 == 0) GO TO 365
 330 IF(ntp /= 0) CALL gopen(itape,z(ibuf2),0)
 DO   i=i1,i2
   j1   = 1
   j2   = ntp
   IF(ntp /= 0) CALL fread(itape,dt,2*j2,0)
   IF(i == ntp) CALL fread(itape,0,0,1)
   IF (ntz == 0)  GO TO 340
   j1   = j2+1
   j2   = j2+ntz
   CALL fread(scr2,dt(j1),2*ntz,0)
   IF(i == ntp) CALL fread(scr2,0,0,1)
   340 CONTINUE
   IF (nty == 0)  GO TO 350
   j1   = j2+1
   j2 = j2+nty
   CALL fread(scr3,dt(j1),2*nty,0)
   IF(i == ntp) CALL fread(scr3,0,0,1)
   350 CONTINUE
   CALL WRITE(scr5,dt,2*j2,0)
 END DO
 IF (nty == 0)  GO TO 370
 IF (nyflag /= 0)  GO TO 370
 IF(ntz /= 0.AND.ntp /= 0)CALL fread(scr2,0,0,1)
 IF(nty /= 0.AND.ntp /= 0)CALL fread(scr3,0,0,1)
 CALL CLOSE(itape,1)
 365 CONTINUE
 nyflag = 1
 i1   = i2+1
 i2   = i2+nty
 itape = scr4
 GO TO 330
 370 CONTINUE
 CALL WRITE(scr5,0,0,1)
 CALL CLOSE(scr1,1)
 CALL CLOSE(scr2,1)
 CALL CLOSE(scr3,1)
 CALL CLOSE(scr4,1)
 CALL CLOSE(scr5,1)
 CALL dmpfil(scr5,z(next),ecore-next-100)
 RETURN
END SUBROUTINE gendsb
