SUBROUTINE ssg1a (n1,ilist,nedt,ntemp,ncent,casecc,iharm)
     
!     ROUTINE ANALIZES CASECC AND SLT TO BUILD LISTS OF SELECTED
!     LOADS
 
 
 INTEGER, INTENT(OUT)                     :: n1
 INTEGER, INTENT(OUT)                     :: ilist(1)
 INTEGER, INTENT(OUT)                     :: nedt
 INTEGER, INTENT(OUT)                     :: ntemp
 INTEGER, INTENT(OUT)                     :: ncent
 INTEGER, INTENT(IN)                      :: casecc
 INTEGER, INTENT(OUT)                     :: iharm
 INTEGER :: system,slt,bgpdt,cstm,sil,ecpt,mpt,gptt,edt,  &
      core(138),NAME(2),name1(2), idefml(1080),itempl(1080),icomb(1080)
 COMMON /loadx / lc,slt,bgpdt,OLD,cstm,sil,isil,ecpt,mpt,gptt,edt,  &
     n(3),lodc,mass,nobld
 COMMON /BLANK / nrowsp,loadnn
 COMMON /system/ system,nout,dum53(53),itherm
 COMMON /zzzzzz/ icore(1)
 EQUIVALENCE     (icore(1),core(1))
 DATA    NAME  / 4HSSG1,4HA   /
 DATA    name1 / 4HSLT ,4HSSG1/
 
 
!     INITIALIZE.
 
 nedt  = 0
 ntemp = 0
 ncent = 0
 ifound= 0
 n1    = 0
 lc1   = lc - system
 islt  = 0
 CALL OPEN (*20,slt,core(lc1+1),0)
 islt = 1
 CALL READ (*320,*10,slt,ilist(1), -2,0,n1)
 CALL READ (*320,*10,slt,ilist(1),lc1,1,n1)
 
!     ALLOW FOR 360 LOADS
 
 10 IF (n1 <= 360) GO TO 12
 NAME(2) = n1
 CALL mesage (-30,137,NAME)
 12 lc1   = lc1 - system
 llist = n1
 20 CALL OPEN (*350,casecc,core(lc1+1),0)
 ione  = 0
 DO  i = 1,loadnn
   CALL fwdrec (*350,casecc)
 END DO
 ifrst = 0
 40 CALL READ (*110,*350,casecc,core(1),166,1,flag)
 IF (ifrst /= 0) GO TO 41
 ifrst = 1
 ispcn = core(3)
 mpcn  = core(2)
 
!     TEST FOR SYMMETRY BUCKLING, OR DIFFERENTIAL STIFFNESS
 
 41 IF (core(16) /= 0 .OR. core(5) /= 0 .OR. core(138) /= 0) GO TO 40
 IF (core(2) /= mpcn .OR. core(3) /= ispcn) GO TO 110
 iharm  = core(136)
 ione   = 1
 IF (core(6) == 0) GO TO 50
 
!     SEE IF EL DEFORM LOAD ALREADY APPLIED
 
 IF (nedt == 0) GO TO 52
 DO  i = 1,nedt
   IF (idefml(i) == core(6)) GO TO 50
 END DO
 
!     ADD TO LIST
 
 52 CONTINUE
 nedt = nedt + 1
 idefml(nedt) = core(6)
 50 IF (core(7) == 0) GO TO 60
 
!     SEE IF TEMP LOAD ALREADY APPLIED
 
 IF (itherm /= 0) GO TO 60
 IF (ntemp  == 0) GO TO 54
 DO  i = 1,ntemp
   IF (itempl(i) == core(7)) GO TO 60
 END DO
 54 CONTINUE
 ntemp = ntemp + 1
 itempl(ntemp) = core(7)
 60 IF (core(4) == 0) GO TO 40
 IF (islt == 0) CALL mesage (-31,core(4),name1)
 IF (n1   == 0) GO TO 90
 DO  i = 1,n1
   IF (core(4) == IABS(ilist(i))) GO TO 100
 END DO
 
!     MUST LOOK AT LOAD CARDS
 
 90 ifound = ifound + 1
 icomb(ifound) = core(4)
 GO TO 40
 100 ilist (i) = -IABS(ilist(i))
 GO TO 40
 110 CALL CLOSE (casecc,1)
 IF (ione  == 0) GO TO 360
 IF (ntemp == 0) GO TO 130
 DO  i = 1,ntemp
   j = n1 + i
   ilist(j) = itempl (i)
 END DO
 130 IF(nedt == 0) GO TO 150
 DO  i = 1,nedt
   j = n1 + ntemp + i
   ilist(j) = idefml(i)
 END DO
 150 IF (ifound == 0) GO TO 270
 
!     LOOK AT LOAD CARDS
 
 DO  i = 1,n1
   CALL fwdrec (*320,slt)
 END DO
 i = 1
 nogo = 0
 CALL READ (*370,*190,slt,core(1),lc1,1,iflag)
 190 llist = n1 + nedt + ntemp
 IF (llist == 0) GO TO 370
 DO  i = 1,ifound
   j = 1
   200 IF (icomb(i) == core(j)) GO TO 220
   j = j + 6
   210 IF (j-1 > iflag) GO TO 255
   IF (core(j-1) == -1) GO TO 200
   j = j + 2
   GO TO 210
   220 j = j + 3
   230 IF (core(j) == -1) CYCLE
   DO  k = 1,llist
     IF (core(j) /= IABS(ilist(k))) CYCLE
     ilist(k) = -IABS(ilist(k))
     j = j + 2
     GO TO 230
   END DO
   255 CALL mesage (31,icomb(i),name1)
   nogo = 1
 END DO
 IF (nogo /= 0) GO TO 390
 270 IF (islt /= 0) CALL CLOSE (slt,1)
 IF (n1 == 0) GO TO 310
 DO  i = 1,n1
   IF (ilist (i) < 0) THEN
     GO TO   290
   ELSE IF (ilist (i) == 0) THEN
     GO TO   300
   END IF
   280 ilist (i) = 0
   CYCLE
   290 ilist (i) = -ilist(i)
 END DO
 310 RETURN
 
!     ERROR MESSAGES.
 
 320 ip1 = slt
 330 ip2 = -1
 CALL mesage (ip2,ip1,NAME)
 350 ip1 = casecc
 GO TO 330
 360 WRITE  (nout,365)
 365 FORMAT ('0*** MISSING LOAD CARD IN CASE CONTROL')
 CALL mesage (-7,0,NAME)
 370 ip2 = 31
 DO  i = 1,ifound
   ip1 = icomb(i)
   CALL mesage (ip2,ip1,name1)
 END DO
 390 CALL mesage (-61,0,NAME)
 RETURN
 
END SUBROUTINE ssg1a
