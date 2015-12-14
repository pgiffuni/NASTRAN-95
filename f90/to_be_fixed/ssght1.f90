SUBROUTINE ssght1 (iest,FILE,nequiv  )
!*****
!     THIS ROUTINE CONVERTS THE EST DATA FOR ALL THERMAL ELEMENTS TO A
!     COMMON FORMAT. OPTIONAL TASKS INCLUDE CALCULATING MATOUT DATA AND
!     CONVERTING  SIL VALUES TO UN VALUES.
!*****
 
 INTEGER, INTENT(IN OUT)                  :: iest
 INTEGER, INTENT(IN OUT)                  :: FILE
 INTEGER, INTENT(IN)                      :: nequiv(1)
 INTEGER :: elid, sub, sil, nesto(45), elem, nest(2),zp, bufm
 INTEGER :: TYPE, pointr(8,23), subr(2),flag
 INTEGER :: point1(8,20),point2(8, 3)
 REAL :: est(100)
 LOGICAL :: linear
 
 COMMON/ condas/ consts(5)
 COMMON/ estout/ elid,sub,NAME(2),sil(8),imat,af,theta,r(3,8), mato(6)
 COMMON/  matin/ matid,inflag,eltemp,dum(1),sinth,costh
 COMMON/ hmtout/ bufm(7)
 COMMON/ gpta1 / nelems, last, incr, elem(1)
 COMMON/ hmatdd/ xxx(4), linear
 
 EQUIVALENCE  (consts(1) , pi     )
 EQUIVALENCE  (nesto(1),elid) ,( nest(1), est(1) )
 EQUIVALENCE  (point1(1,1),pointr(1,1)), (point2(1,1),pointr(1,21))
 
 DATA subr / 4HSSGH ,4HT1   /
 DATA numelt / 23 /
!*****
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
     ,67  ,0   ,0   ,34  ,0   ,40  ,32  ,1  &
     ,76  ,0   ,11  ,12  ,13  ,14  ,8   ,17 /
!*****
 CALL delset
 10 CALL READ(*120,*140,iest,TYPE,1,0,flag)
 DO  i =1,numelt
   iel = i
   IF (TYPE   - pointr(1,i) < 0.0) THEN
     GO TO    30
   ELSE IF (TYPE   - pointr(1,i) == 0.0) THEN
     GO TO    40
   ELSE
     GO TO    20
   END IF
 END DO
 30 CALL fwdrec(*150,iest)
 GO TO 10
 
 40 CONTINUE
 zp =   (TYPE-1)*incr
 NAME(1)= elem(zp+1)
 NAME(2)= elem(zp+2)
 nwords = elem(zp+12)
 50 CONTINUE
 CALL READ(*150,*130,iest,est,nwords,0,flag)
 elid = nest(1)
 DO  i = 5,45
   nesto(i) = 0
 END DO
 IF( TYPE == 3) est(5) = pi*est(6)*(est(5)-est(6))
 IF( TYPE == 52 .AND. nest(2) == 7) est(16)=pi*(est(19)+est(20))
 is = pointr(2,iel)
 ith= pointr(3,iel)
 im = pointr(4,iel)
 ia = pointr(5,iel)
 ig = pointr(6,iel)
 sub= pointr(8,iel)
 np = pointr(7,iel)
 
 IF(sub == 10) sub = sub + nest(2)-1
 inflag =1
 IF( sub >= 16) inflag=3
 IF( sub < 2 .OR. sub > 5)  GO TO 60
 inflag =2
 GO TO 70
 60 IF(sub  < 6 .OR. sub > 9)  GO TO 70
 inflag =3
 70 CONTINUE
 IF( ia > 0) af = est(ia)
 matid = nest(im)
 IF(matid <= 0) GO TO 50
 sinth=0.0
 costh=1.0
 IF( inflag /= 2) GO TO 80
 theta= est(ith)*pi/180.0
 IF( theta == 0.0)GO TO 80
 sinth= SIN(theta)
 costh= COS(theta)
 80 itemp = ig + 4*np
 eltemp = est(itemp)
 imat = matid
 linear=.false.
 CALL hmat( elid )
!*****
!     TEST IF NONLINEAR
!*****
 IF( linear ) GO TO 50
 DO  i=1,6
   mato(i)= bufm(i)
 END DO
 DO  i=1,np
   jpoint = 4*(i-1) + ig
   DO  j=1,3
     iloc = jpoint + j
     r(j,i) = est(iloc)
   END DO
   isil= is+i + 1
   ipt    =nest(isil)
   IF( ipt == 0) CYCLE
   sil(i) = nequiv(ipt)
 END DO
!*****
!     WRITE A UNIFORM EST GROUP OF CONVERTED DATA HERE
!*****
 CALL WRITE(FILE,nesto(1), 45, 0 )
!*****
!     RETURN FOR ANOTHER ELEMENT
!******
 GO TO 50
 120 RETURN
!*****
!     DONE WITH THIS ELEMENT TYPE
!*****
 130 IF( flag == 0) GO TO 10
!******
 140 j=-3
 GO TO 160
 150 j = -2
 160 CALL mesage(j,iest,subr)
 RETURN
END SUBROUTINE ssght1
