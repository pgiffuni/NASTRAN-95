SUBROUTINE mred1b (mode)
     
!     THIS SUBROUTINE PROCESSES THE BDYS AND BDYS1 DATA FOR THE FIXED
!     IDENTIFICATION SET (FIXSET) AND THE BOUNDARY IDENTIFICATION SET
!     (BNDSET) FOR THE MRED1 MODULE.
 
!     INPUT DATA
!     GINO   - GEOM4  - BDYS DATA
!                     - BDYS1 DATA
!     OTHERS - MODE   - SUBROUTINE PROCESSING FLAG
!                     = 1, PROCESS FIXED ID SET
!                     = 2, PROCESS BOUNDARY ID SET
 
!     OUTPUT DATA
!     GINO   - USETX  - S,R,B DEGREES OF FREEDOM
 
!     PARAMETERS
!     INPUT  - NOUS   - FIXED POINTS FLAG
!                       .GE.  0, FIXED POINTS DEFINED
!                       .EQ. -1, NO FIXED POINTS DEFINED
!              GBUF1  - GINO BUFFER
!              KORLEN - CORE LENGTH
!              IO     - OUTPUT OPTION FLAG
!              NAMEBS - BEGINNING ADDRESS OF BASIC SUBSTRUCTURES NAMES
!              EQSIND - BEGINNING ADDRESS OF EQSS GROUP ADDRESSES
!              NSLBGN - BEGINNING ADDRESS OF SIL DATA
!              KBDYC  - BEGINNING ADDRESS OF BDYC DATA
!              USETX  - USETX OUTPUT FILE NUMBER
!              NBDYCC - NUMBER OF BDYC WORDS
!     OUTPUT - DRY    - MODULE OPERATION FLAG
!     OTHERS - LOCUST - BEGINNING ADDRESS OF USET ARRAY
!              IERR   - NO BDYS/BDYS1 DATA ERROR FLAG
!                       .LT. 2, NO ERRORS
!                       .EQ. 2, ERRORS
!              GRPBGN - ABSOLUTE BEGINNING ADDRESS OF EQSS GROUP DATA
!              GRPEND - ABSOLUTE ENDING ADDRESS OF EQSS GROUP DATA
!              GRPIP  - ABSOLUTE ADDRESS OF EQSS DATA GROUP
!              LOCBGN - BEGINNING ADDRESS OF EQSS DATA FOR SUBSTRUCTURE
!              NFOUND - NUMBER OF EQSS DATA ITEMS FOUND FOR SET ID
!              KPNTBD - ARRAY OF BDYC DOF COMPONENTS
!              KPNTSL - ARRAY OF EQSS DOF COMPONENTS
!              INDSIL - ABSOLUTE INDEX INTO SIL DATA
!              NSILUS - ABSOLUTE INDEX INTO USET ARRAY
 
 
 INTEGER, INTENT(IN)                      :: mode
 IMPLICIT INTEGER (a-z)
 EXTERNAL        rshift,andf,orf,complf
 LOGICAL :: bounds,ponly
 REAL :: rz(1)
 DIMENSION       array(3),bdyi(2,2),bdy(2),eqstrl(7),idum(3),  &
     kpntsl(32),ioshft(2),kpntbd(9),modnam(2),usetrl(7)
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm
 COMMON /BLANK / oldnam(2),dry,idum13,nous,idum2(4),gbuf1,  &
     idum14(4),korlen,idum4(5),io,idum5(2),bndset,  &
     fixset,ieig,korbgn,idum12,namebs,eqsind,nslbgn,  &
     idum6,kbdyc,nbdycc,luset,locust,idum3(4),bounds, ponly
 COMMON /zzzzzz/ z(1)
 COMMON /system/ idum7,iprntr,idum8(6),nlpp,idum9(2),line
 COMMON /two   / itwo(32)
 COMMON /bitpos/ idum10(5),ul,ua,uf,us,un,idum11(10),ub,ui
 COMMON /patx  / lcore,nsub(3),fuset
 COMMON /unpakx/ typeu,irowu,nrowu,incru
 EQUIVALENCE     (rz(1),z(1))
 DATA    geom4 , bdyi,usetx/102,1210,12,1310,13,201/
 DATA    modnam/ 4HMRED,4H1B   /
 DATA    ioshft/ 11,2  /
 DATA    item  / 4HUPRT/
 DATA    uprt  , eqst  /301,203/
 
!     TEST FOR FIXED SET INPUT
 
 IF (nous == -1 .AND. mode == 1) GO TO 430
 
!     CHECK FOR LOADS PROCESSING ONLY
 
 IF (ponly) GO TO 345
 
!     PROCESS BDY(S/S1) BULK DATA FOR SPECIFIED BDYC
 
 ishift = ioshft(mode)
 IF (nbdycc == 0) GO TO 335
 korbgn = kbdyc + 4*nbdycc
 IF (korbgn >= korlen) GO TO 300
 IF (andf(rshift(io,ishift),1) == 0) GO TO 10
 CALL page1
 IF (mode == 1) WRITE (iprntr,900)
 IF (mode == 2) WRITE (iprntr,901)
 line  = line + 7
 10 ibits = itwo(ul) + itwo(ua) + itwo(uf)
 IF (mode == 2) ibits = itwo(ui)
 ibits = complf(ibits)
 ierr  = 0
 ibdy  = 0
 ifile = geom4
 
!     SET BULK DATA PROCESSING FLAG AND READ SET ID
!     IBDY .EQ. 1 - BDYS
!     IBDY .EQ. 2 - BDYS1
 
 nxtbdy = 1
 ifound = 0
 CALL preloc (*270,z(gbuf1),geom4)
 20 ibdy = ibdy + 1
 IF (ibdy == 3) GO TO 260
 DO  i = 1,2
   bdy(i) = bdyi(i,ibdy)
 END DO
 CALL locate (*250,z(gbuf1),bdy,iflag)
 GO TO 40
 35 CALL bckrec (geom4)
 nxtbdy = nxtbdy + 1
 IF (nxtbdy > nbdycc) GO TO 20
 CALL READ (*280,*290,geom4,idum,3,0,iflag)
 40 CALL READ (*280,*20,geom4,array,ibdy,0,iflag)
 
!     CHECK REQUEST ID
 
 bdyj = 2
 bdyk = 2
 bdyl = 3
 bdym = 2
 IF (ibdy == 1) GO TO 50
 bdyj = 3
 bdyk = 1
 bdyl = 2
 bdym = 3
 50 iwds = 2 + 4*(nxtbdy-1)
 DO  i = nxtbdy, nbdycc
   IF (z(kbdyc+iwds) == array(1)) GO TO 90
   iwds = iwds + 4
 END DO
 
!     FINISH BDY(S/S1) SET ID READING
 
 60 CALL READ (*280,*290,geom4,array(bdyj),bdyk,0,iflag)
 IF (ibdy - 2 < 0) THEN
   GO TO    70
 ELSE
   GO TO    80
 END IF
 70 IF (array(2) /= -1 .AND. array(3) /= -1) GO TO 60
 GO TO 40
 80 IF (array(3) /= -1) GO TO 60
 GO TO 40
 
!     CONTINUE BDY(S/S1) SET ID PROCESSING
 
 90 CALL READ (*280,*290,geom4,array(bdyj),bdyk,0,iflag)
 IF (ibdy - 2 < 0) THEN
   GO TO   100
 ELSE
   GO TO   110
 END IF
 100 IF (array(2) == -1 .AND. array(3) == -1) GO TO 115
 GO TO 120
 110 IF (array(3) == -1) GO TO 115
 GO TO 120
 
!     CHECK FOR NEXT BDY(S/S1) CARD HAVING SAME SET ID AS CURRENT ID
 
 115 CALL READ (*280,*35,geom4,array,ibdy,0,iflag)
 IF (z(kbdyc+iwds) == array(1)) GO TO 90
 GO TO 35
 
!     LOCATE EQSS DATA FOR SUBSTRUCTURE
 
 120 ifound = 1
 ip     = 2*(z(kbdyc+iwds+1)-1)
 grpbgn = z(eqsind+ip)
 grpend = grpbgn + z(eqsind+ip+1)
 k = z(eqsind+ip+1)/3
 CALL bisloc (*170,array(bdym),z(grpbgn),3,k,locbgn)
 grpip  = grpbgn + locbgn - 1
 loc    = grpip - 3
 130 IF (loc < grpbgn) GO TO 140
 IF (z(loc) < z(grpip)) GO TO 140
 loc    = loc - 3
 GO TO 130
 140 locbgn = loc + 3
 nfound = 1
 loc    = locbgn + 3
 150 IF (loc >= grpend) GO TO 180
 IF (z(locbgn) < z(loc)) GO TO 180
 loc    = loc + 3
 nfound = nfound + 1
 GO TO 150
 
!     CANNOT LOCATE EXTERNAL ID
 
 170  CALL page1
 IF (mode == 1) WRITE (iprntr,902) ufm,array(3),array(2),  &
     array(1),z(namebs+ip),z(namebs+ip+1)
 IF (mode == 2) WRITE (iprntr,903) ufm,array(3),array(2),  &
     array(1),z(namebs+ip),z(namebs+ip+1)
 dry = -2
 GO TO 90
 
!     LOCATE CORRECT IP FOR THIS EXTERNAL ID
 
 180 CALL splt10 (array(bdyl),kpntbd,jwds)
 m = 0
 DO  i = 1, nfound
   j = (3*(i-1)) + 2
   icode = z(locbgn+j)
   CALL decode (icode,kpntsl,kwds)
   DO  k = 1, kwds
     DO  l = 1, jwds
       IF (kpntsl(k) == kpntbd(l)-1) GO TO 200
     END DO
     CYCLE
     
!     CONVERT GRID ID AND COMPONENT TO SIL VALUE
     
     200 IF (andf(rshift(io,ishift),1) == 0) GO TO 220
     IF (line <= nlpp) GO TO 210
     CALL page1
     IF (mode == 1) WRITE (iprntr,900)
     IF (mode == 2) WRITE (iprntr,901)
     line = line + 7
     210 IF (m == 0) WRITE (iprntr,906) array(1),array(bdym),array(bdyl)
     m    = 1
     line = line + 1
     220 indsil = nslbgn + ((2*z(locbgn+j-1))-2)
     nsilus = locust + ((z(indsil)-1)+(k-1))
     kpntbd(l) = 0
     
!     FILL USET ARRAY
!     IF FIXSET - TURN OFF UL, UA, UF BITS AND TURN ON US BIT
!     IF BNDSET - TURN OFF UI BIT AND TURN ON UB BIT
     
     ubors = us
     IF (mode == 2) ubors = ub
     z(nsilus) = andf(z(nsilus),ibits)
     z(nsilus) = orf(z(nsilus),itwo(ubors))
   END DO
 END DO
 
!     CHECK THAT ALL IP FOUND
 
 DO  i = 1,jwds
   IF (kpntbd(i) == 0) CYCLE
   IF (mode == 1) WRITE (iprntr,904) uwm,array(bdym),z(namebs+ip),  &
       z(namebs+ip+1)
   IF (mode == 2) WRITE (iprntr,905) uwm,array(bdym),z(namebs+ip),  &
       z(namebs+ip+1)
   GO TO 90
 END DO
 GO TO 90
 
!     SET NO DATA AVAILABLE FLAG
 
 250 ierr = ierr + 1
 GO TO 20
 
!     END OF ID SET PROCESSING
 
 260 CALL CLOSE (geom4,1)
 IF (ierr   == 2) GO TO 330
 IF (ifound == 0) GO TO 330
 IF (mode   == 1) GO TO 430
 
!     WRITE USETX DATA
 
 CALL gopen (usetx,z(gbuf1),1)
 CALL WRITE (usetx,z(locust),luset,1)
 CALL CLOSE (usetx,1)
 usetrl(1) = usetx
 usetrl(2) = 1
 usetrl(3) = luset
 usetrl(4) = 7
 usetrl(5) = 1
 CALL wrttrl (usetrl)
 
!     VERIFY OLD BOUNDARY UNCHANGED
 
 IF (.NOT.bounds) GO TO 430
 IF (locust+2*luset >= korlen) GO TO 300
 345 CALL softrl (oldnam,item,usetrl)
 IF (usetrl(1) /= 1) GO TO 440
 nrowu = usetrl(3)
 IF (ponly) luset = nrowu
 IF (nrowu /= luset) GO TO 420
 
!     GET OLD UPRT VECTOR
 
 typeu = usetrl(5)
 CALL mtrxi (uprt,oldnam,item,0,itest)
 newust = locust + luset
 IF (ponly) newust = locust
 IF (ponly .AND. newust+nrowu >= korlen) GO TO 300
 irowu = 1
 incru = 1
 CALL gopen (uprt,z(gbuf1),0)
 CALL unpack (*350,uprt,rz(newust))
 GO TO 370
 350 DO  i = 1,luset
   rz(newust+i-1) = 0.0
 END DO
 370 CALL CLOSE (uprt,1)
 IF (ponly) GO TO 405
 
!     GET NEW UPRT VECTOR
 
 lcore = korlen - (newust+luset)
 fuset = usetx
 CALL calcv (uprt,un,ui,ub,z(newust+luset))
 typeu = 1
 nrowu = luset
 CALL gopen (uprt,z(gbuf1),0)
 CALL unpack (*380,uprt,rz(newust+luset))
 GO TO 400
 380 DO  i = 1,luset
   rz(newust+luset+i-1) = 0.0
 END DO
 400 CALL CLOSE (uprt,1)
 
!     CHECK OLD, NEW UPRT VECTORS AND COUNT NUMBER OF ROWS IN 0, 1
!     SUBSETS AND SAVE IN EQST TRAILER FOR USE IN MRED2A
 
 405 isub0 = 0
 isub1 = 0
 DO  i = 1,luset
   IF (rz(newust+i-1) == 0.0) isub0 = isub0 + 1
   IF (rz(newust+i-1) == 1.0) isub1 = isub1 + 1
   IF (ponly) CYCLE
   IF (rz(newust+i-1) /= rz(newust+luset+i-1)) GO TO 420
 END DO
 eqstrl(1) = eqst
 eqstrl(6) = isub0
 eqstrl(7) = isub1
 CALL wrttrl (eqstrl)
 GO TO 430
 
!     BOUNDARY POINTS ARE NOT THE SAME
 
 420 WRITE (iprntr,909) ufm,oldnam
 dry = -2
 430 CONTINUE
 RETURN
 
!     PROCESS SYSTEM FATAL ERRORS
 
 270 imsg = -1
 GO TO 320
 280 imsg = -2
 GO TO 310
 290 imsg = -3
 GO TO 310
 300 imsg = -8
 ifile = 0
 310 CALL CLOSE (geom4,1)
 320 CALL sofcls
 CALL mesage (imsg,ifile,modnam)
 RETURN
 
!     PROCESS MODULE FATAL ERRORS
 
 330 IF (mode == 1) WRITE (iprntr,907) ufm,fixset
 IF (mode == 2) WRITE (iprntr,908) ufm,bndset
 335 dry = -2
 CALL sofcls
 CALL CLOSE (geom4,1)
 RETURN
 440 SELECT CASE ( itest )
   CASE (    1)
     GO TO 450
   CASE (    2)
     GO TO 450
   CASE (    3)
     GO TO 460
   CASE (    4)
     GO TO 470
   CASE (    5)
     GO TO 480
   CASE (    6)
     GO TO 480
 END SELECT
 450 WRITE (iprntr,910) ufm,modnam,item,oldnam
 dry = -2
 RETURN
 460 imsg = -1
 GO TO 490
 470 imsg = -2
 GO TO 490
 480 imsg = -3
 490 CALL smsg (imsg,item,oldnam)
 RETURN
 
 900 FORMAT (//45X,40HTABLE of grid points composing fixed set,  &
     //53X,5HFIXED,/53X,25HSET id   grid point   dof, /53X,  &
     26HNUMBER   id NUMBER    code,/)
 901 FORMAT (1H0,44X,43HTABLE of grid points composing boundary set,  &
     //52X,8HBOUNDARY,/53X,25HSET id   grid point   dof, /53X,  &
     26HNUMBER   id NUMBER    code,/)
 902 FORMAT (a23,' 6624, GRID POINT',i9,' COMPONENT',i9,' SPECIFIED ',  &
     'IN FIXED SET',i9, /5X,'FOR SUBSTRUCTURE ',2A4, ' DOES NOT EXIST.',//////)
 903 FORMAT (a23,' 6611, GRID POINT',i9,' COMPONENT',i9,' SPECIFIED ',  &
     'IN BOUNDARY SET',i9, /5X,'FOR SUBSTRUCTURE ',2A4, ' DOES NOT EXIST.',//////)
 904 FORMAT (a25,' 6625, DEGREES OF FREEDOM AT GRID POINT',i9,  &
     ' COMPONENT SUBSTRUCTURE ',2A4, /32X,'INCLUDED IN A FIXED',  &
 ' SET DO NOT EXIST.  REQUEST WILL BE IGNORED.')
   905 FORMAT (a25,' 6610, DEGREES OF FREEDOM AT GRID POINT',i9,  &
       ' COMPONENT SUBSTRUCTURE ',2A4, /32X,'INCLUDED IN A NON-',  &
       'EXISTING BOUNDARY SET.  REQUEST WILL BE IGNORED.')
   906 FORMAT (52X,2(i8,3X),i6)
   907 FORMAT (a23,' 6626, NO BDYS OR BDYS1 BULK DATA HAS BEEN INPUT TO',  &
       ' DEFINE FIXED SET',i9,1H.)
   908 FORMAT (a23,' 6607, NO BDYS OR BDYS1 BULK DATA HAS BEEN INPUT TO',  &
       ' DEFINE BOUNDARY SET',i9,1H.)
   909 FORMAT (a23,' 6637, OLDBOUND HAS BEEN SPECIFIED BUT THE BOUNDARY',  &
       ' POINTS FOR SUBSTRUCTURE ',2A4,' HAVE BEEN CHANGED.')
   910 FORMAT (a23,' 6215, MODULE ',2A4,8H - item ,a4,' OF SUBSTRUCTURE '  &
       ,       2A4,' PSEUDO-EXISTS ONLY.')
   
 END SUBROUTINE mred1b
