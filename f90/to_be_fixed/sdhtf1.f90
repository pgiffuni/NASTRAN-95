SUBROUTINE sdhtf1 (TYPE,reject)
     
!     THIS ROUTINE CONVERTS THE EST DATA FOR ALL THERMAL ELEMENTS TO A
!     COMMON FORMAT. SDHT1B IS CALLED TO PRODUCE THE OUTPUT
 
 
 INTEGER, INTENT(IN)                      :: TYPE
 LOGICAL, INTENT(OUT)                     :: reject
 
 INTEGER :: elid,sub,sil,nesto(100),elem,nest(2),  &
     pointr(8,23),typold,strspt,estwds,nestsc(200),  &
     point1(8,20),point2(8, 3),tube,ftube,chbdy
 DIMENSION       shp(32),dshp(3,32),xjacob(3,3),bxyz(3,32),gpt(32)
 COMMON /condas/ consts(5)
 COMMON /sdr2x4/ dumx(109),strspt,ddrmm,isopl8
 COMMON /sdr2x5/ est(100),elid,sil(32),nq,np,NAME(2),domx(105), dshpb(3,32)
 COMMON /sdr2x6/ sub,imat,af,theta,r(3,32),estscr(200)
 COMMON /gpta1 / nels,last,incr,elem(1)
 COMMON /matin / matid,inflag,eltemp,dum(1),sinth,costh
 EQUIVALENCE     (consts(1),pi),(nestsc(1),estscr(1)),  &
     (nesto(1),sub),(nest(1),est(1)), (point1(1,1),pointr(1,1)),  &
     (point2(1,1),pointr(1,21))
 DATA    typold, numelt, tube, ftube, chbdy / 0,     23,    3,    82,    52 /
!     DATA    HEX   / 16    /
 
!     THE POINTERS TO THE EST DATA ARE
 
!        IM    MAT ID
!        ITH   THETA
!        IA    AREA
!        IG    GRID POINT DATA
!        IS    SIL MINUS 1
!        NP    NO. OF POINTS
!        SUB   SUBROUTINE TYPE
!                       NO.  IS   ITH  IM   IA   IG   NP   SUB
!                      ----  --   ---  --   --   --   --   ----
 DATA   point1 /    1   ,0   ,0   ,4   ,5   ,9   ,2   ,1  &
     ,3   ,0   ,0   ,4   ,5   ,8   ,2   ,1  &
     ,6   ,0   ,5   ,6   ,7   ,15  ,3   ,2  &
     ,9   ,0   ,5   ,6   ,7   ,9   ,3   ,2  &
     ,10  ,0   ,0   ,4   ,5   ,9   ,2   ,1  &
     ,16  ,0   ,6   ,7   ,8   ,10  ,4   ,3  &
     ,17  ,0   ,5   ,6   ,7   ,9   ,3   ,2  &
     ,18  ,0   ,6   ,7   ,8   ,10  ,4   ,3  &
     ,19  ,0   ,6   ,7   ,8   ,16  ,4   ,3  &
     ,34  ,0   ,0   ,16  ,17  ,34  ,2   ,1  &
     ,36  ,0   ,5   ,6   ,0   ,7   ,3   ,4  &
     ,37  ,0   ,6   ,7   ,0   ,8   ,4   ,5  &
     ,39  ,1   ,0   ,2   ,0   ,7   ,4   ,6  &
     ,40  ,1   ,0   ,2   ,0   ,9   ,6   ,7  &
     ,41  ,1   ,0   ,2   ,0   ,11  ,8   ,8  &
     ,42  ,1   ,0   ,2   ,0   ,11  ,8   ,9  &
     ,52  ,1   ,0   ,15  ,16  ,21  ,8   ,10  &
     ,62  ,0   ,6   ,7   ,8   ,10  ,4   ,3  &
     ,63  ,0   ,6   ,7   ,8   ,10  ,4   ,3  &
     ,65  ,0   ,0   ,10  ,0   ,16  ,8   ,16 /
 DATA   point2 /    66  ,0   ,0   ,22  ,0   ,28  ,20  ,16  &
     ,67  ,0   ,0   ,34  ,0   ,40  ,32  ,16  &
     ,76  ,0   ,11  ,12  ,13  ,14  ,8   ,17 /
 
 IF (TYPE ==  ftube) GO TO 115
 IF (TYPE == typold) GO TO 50
 typold = TYPE
 reject = .true.
 DO  i = 1,numelt
   iel = i
   IF (TYPE-pointr(1,i) < 0.0) THEN
     GO TO    30
   ELSE IF (TYPE-pointr(1,i) == 0.0) THEN
     GO TO    40
   ELSE
     GO TO    20
   END IF
 END DO
 30 RETURN
 
 40 reject = .false.
 50 IF ((TYPE >= 65.AND.TYPE <= 67) .AND. strspt == 0) strspt = strspt + 1
 ip = (TYPE-1)*incr
 estwds = elem(ip+12)
 
!     THE LOCATIONS OF DATA FOR EACH PARTICULAR ELEMENT ARE ZEROED OUT
 
 nq = 0
 DO  i = 1,100
   nesto(i) = 0
 END DO
 NAME(1) = elem(ip+1)
 NAME(2) = elem(ip+2)
 elid = nest(1)
 DO  i = 1,32
   nest(i+101) = 0
 END DO
 DO  i = 1,201
   nest(i+137) = 0
 END DO
 IF (TYPE == tube) est(5) = pi*estscr(6)*(estscr(5)-estscr(6))
 IF (TYPE == chbdy .AND. nestsc(2) == 7) est(16) = pi*(estscr(19)+estscr(20))
 is  = pointr(2,iel)
 ith = pointr(3,iel)
 im  = pointr(4,iel)
 ia  = pointr(5,iel)
 ig  = pointr(6,iel)
 sub = pointr(8,iel)
 np  = pointr(7,iel)
 
 IF (sub == 10) sub = sub + nestsc(2) - 1
 inflag = 1
 IF (sub >= 16) inflag = 3
 IF (sub < 2 .OR. sub > 5) GO TO 60
 inflag = 2
 GO TO 70
 60 IF (sub < 6 .OR. sub > 9) GO TO 70
 inflag = 3
 70 CONTINUE
 IF (sub /= 16) GO TO 79
 
!     GET SHAPE FUNCTIONS ETC. FOR STRESS POINT(ALSO DETERMINE THE
!     STRESS POINT, WHICH WILL BE THE GRID POINTS PLUS CENTROID IN
!     ELEMENT COORDINATES
 
 itype = TYPE - 64
 DO  i = 1,np
   gpt(i) = estscr(5*np+7+i)
   DO  j = 1,3
     bxyz(j,i) = estscr(np+4+4*i+j)
   END DO
 END DO
 
!     GET STRESS POINT
 
 y =-1.
 z =-1.
 IF (itype > 1) GO TO 502
 d = 2.
 x = 1.
 GO TO 505
 502 d = 1.
 x = 0.
 505 IF (itype > 1) GO TO 560
 SELECT CASE ( strspt )
   CASE (    1)
     GO TO 510
   CASE (    2)
     GO TO 520
   CASE (    3)
     GO TO 530
   CASE (    4)
     GO TO 510
   CASE (    5)
     GO TO 540
   CASE (    6)
     GO TO 520
   CASE (    7)
     GO TO 530
   CASE (    8)
     GO TO 510
   CASE (    9)
     GO TO 550
 END SELECT
 510 x = x - d
 GO TO 590
 520 x = x + d
 GO TO 590
 530 y = y + d
 GO TO 590
 540 z = z + d
 y = -1.
 GO TO 590
 550 x = 0.
 y = 0.
 z = 0.
 GO TO 590
 560 GO TO (510,520,520,530,530,510,510,570,580,520,  &
     530,510,580,520,520,530,530,510,510,570, 550), strspt
 570 y = y - d
 GO TO 590
 580 z = z + 1.
 y = -1.
 d = 3. - d
 590 CALL ihexss (itype,shp,dshp,xjacob,detj,elid,x,y,z,bxyz)
 
!     GET DERIVATIVES W.R.T.X,Y,Z(REVERSE CALLING SEQUENCE BECAUSE
!     COLUMN-STORED
 
 CALL gmmats (dshp,np,3,0,xjacob,3,3,0,dshpb)
 
 79 CONTINUE
 
 IF (ia > 0) af = estscr(ia)
 matid = nestsc(im)
 IF (matid <= 0) RETURN
 sinth = 0.0
 costh = 1.0
 IF (inflag /= 2) GO TO 80
 theta = estscr(ith)*pi/180.
 IF (theta == 0.0) GO TO 80
 sinth = SIN(theta)
 costh = COS(theta)
 80 itemp = ig + 4*np
 eltemp= estscr(itemp)
 IF (sub /= 16) GO TO 85
 isopl8= 8
 eltemp= 0.
 DO  i = 1,np
   eltemp= eltemp + gpt(i)*shp(i)
 END DO
 85 CONTINUE
 imat  = matid
 CALL hmat (elid)
 
 DO  i = 1,np
   ip = 4*(i-1) + ig
   DO  j = 1,3
     iloc = ip + j
     r(j,i) = estscr(iloc)
   END DO
   isil   = is + i + 1
   sil(i) = nestsc(isil)
 END DO
 
 CALL sdhtff
 GO TO 120
 
!     FTUBE CONVECTION ELEMENT
 
 115 reject =.false.
 i = 0
 nest(i+101) = nestsc(1)
 nest(i+102) = nestsc(2)
 nest(i+103) = nestsc(3)
 est (i+104) = estscr(4)*estscr(5)
 est (i+105) = 0.0
 
 120 RETURN
END SUBROUTINE sdhtf1
