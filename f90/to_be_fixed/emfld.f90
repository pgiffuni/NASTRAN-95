SUBROUTINE emfld
     
!                                            SEE T01191A     ===========
!     COMPUTES TOTAL MAGNETIC FIELD STRENGTH AND INDUCTION FOR
!     EACH ELEMENT IN BASIC COORDINATES BY ADDING HM AND HC
 
!     EMFLD    HOEF1,HEST,CASECC,HCFLD,MPT,DIT,REMFLD,GEOM1,CSTM,HCCEN/
!              HOEH1/V,N,HLUSET $
 
 INTEGER :: hest,hoeh1,hoef1,casecc,estfld,hluset,dit,  &
     FILE,buf1,buf2,buf3,buf4,buf5,sysbuf,otpe,typout,  &
     eltype,subcas,elid,oldcas,oldeid,strspt,  &
     remfl,buf6,idum(2),geom1,cstm,hccen,hcount
 DIMENSION       coord(4),icoord(4),ta(9),temp(3),mcb(7),hmg(3),  &
     hm(3),hc(3),ibuf(150),rbuf(150),nam(2),iz(1),zn(2)
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm,uim,sfm
 COMMON /BLANK / hluset
 COMMON /gpta1 / nelems,last,incr,NE(1)
 COMMON /system/ sysbuf,otpe
 COMMON /unpakx/ typout,ii,nn,incur
 COMMON /zzzzzz/ z(1)
 EQUIVALENCE     (rbuf(1),ibuf(1)),(coord(1),icoord(1)), (z(1),iz(1))
 DATA    hoef1 , hest,casecc,mpt,dit/ 101,102,103, 105,106/
 DATA    remfl , geom1,cstm,hccen   / 107,108,109, 110    /
 DATA    estfld, hoeh1/301,201 /
 DATA    nam   / 4HEMFL,4HD    /, zn/ 4HHOEH,4H1          /
 DATA    hex1  , hex2, hex3    /  4HHEX1,4HHEX2,4HHEX3    /
 
!     CHECK TO SEE IF HOEF1 EXISTS. IF NOT, THEN NO MAG. FIELD REQUESTS
 
 mcb(1) = hoef1
 CALL rdtrl (mcb)
 IF (mcb(1) < 0) GO TO 600
 mcb(1) = hccen
 CALL rdtrl (mcb)
 nn = mcb(3)
 IF (mcb(1) > 0) GO TO 20
 nn = 0
 mcb(1) = remfl
 CALL rdtrl (mcb)
 IF (mcb(1) > 0) GO TO 20
 WRITE  (otpe,10)  uwm
 10 FORMAT (a25,', DATA BLOCKS HCFLD AND REMFL ARE PURGED IN EM ',  &
     'PROBLEM. ALL RESULTS ARE ZERO')
 GO TO 600
 
 20 mcb(1) = hest
 CALL rdtrl (mcb)
 nelx = 3*mcb(2)
 
 typout = 1
 ii     = 1
 incur  = 1
 
!     CREATE ESTFLD WHICH LOOKS LIKE HEST BUT CONTAINS ONLY TYPE, ID,
!     NUMBER OF SILS,SILS,3 X 3 MATERAIL MATRIX,AND 3 X 3 TRANSFORMATION
!     MATRIX FROM LOCAL TO BASIC,BFIELD,AND COORDS OF STRESS POINT FOR
!     NON-RECTANGULAR BFIELD
 
 CALL estmag (hest,estfld,mpt,dit,geom1,iany,kcount)
 
!     KCOUNT SHOULD BE NUMBER OF TERMS IN ROW OF HCCEN
 
 IF (nn == 0) nn = kcount
 IF (nn /= kcount) GO TO 500
 nrows = nn
 
!     NOW FETCH HC AT EACH POINT FROM HCFLD
 
 lcore = korsz(z)
 buf1  = lcore- sysbuf + 1
 buf2  = buf1 - sysbuf
 buf3  = buf2 - sysbuf
 buf4  = buf3 - sysbuf
 buf5  = buf4 - sysbuf
 buf6  = buf5 - sysbuf
 lcore = buf6 - 1
 IF (lcore <= 0) GO TO 550
 
 ncount = 0
 oldcas = 0
 hcount = 0
 
!     COPY HEADER FROM HOEF1 TO HOEH1
 
 FILE = hoeh1
 CALL OPEN (*520,hoeh1,z(buf5),1)
 FILE = hoef1
 CALL OPEN (*520,hoef1,z(buf4),0)
 CALL READ (*530,*30,hoef1,z,lcore,0,iwords)
 GO TO 550
 30 z(1) = zn(1)
 i    = 2
 z(i) = zn(2)
 CALL WRITE (hoeh1,z,iwords,1)
 
!     OPEN CSTM FOR NON-BASIC COORDINATE SYSTEM
 
 ncstm = 0
 IF (iany == 0) GO TO 50
 CALL gopen (cstm,z(buf1),0)
 CALL READ  (*530,*40,cstm,z,lcore,0,ncstm)
 GO TO 550
 40 CALL CLOSE  (cstm,1)
 CALL pretrs (z(1),ncstm)
 
 50 nncr = ncstm + nrows
 nall = nncr  + nelx
 CALL gopen (casecc,z(buf1),0)
 CALL gopen (hccen,z(buf2),0)
 CALL gopen (estfld,z(buf3),0)
 CALL gopen (remfl,z(buf6),0)
 
!     READ ID RECORD FROM HOEF1. COPY TO HOEH1, EXCEPT CHANGE NUMBER OF
!     WORDS FROM +9 TO -9 AS AN INDICATOR FOR TITLES IN OFP. (+9 IS FOR
!     HEAT TRANSFER) ALSO PICK UP SUBCASE NUMBER AND ELEMENT TYPE. IF
!     SAME SUBCASE AS PREVIUUS ONE, USE SAME HCFLD VECTOR. IF NOT,
!     CREATE A NEW ONE
 
 60 CALL READ (*410,*540,hoef1,ibuf,146,1,iwords)
 ibuf(10) = -9
 CALL WRITE (hoeh1,ibuf,146,1)
 eltype = ibuf(3)
 subcas = ibuf(4)
 oldeid = 0
 IF (subcas == oldcas) GO TO 260
 oldcas = subcas
 ncount = 0
 hcount = 0
 CALL REWIND (estfld)
 FILE   = estfld
 CALL fwdrec (*530,estfld)
 
!     IF THIS SUBCASE IS NOT A SUBCOM, UNPACK NEXT COLUMN OF HCFLD. IF
!     IT IS A SUBCOM, BCKREC HCFLD THE SAME NUMBER OF RECORDS AS THERE
!     ARE FACTORS ON THE SUBSEQ AND COMBINE VECTORS TO PRODUCE ONE
!     VECTOR.
 
 IF (lcore < 16) GO TO 550
 70 FILE = casecc
 CALL READ (*530,*540,casecc,z(ncstm+1),16,0,iwords)
 IF (iz(ncstm+1) == subcas) GO TO 80
 CALL fwdrec (*530,casecc)
 FILE = hccen
 CALL fwdrec (*530,hccen)
 FILE = remfl
 CALL fwdrec (*530,remfl)
 GO TO 70
 
!     MATCH ON SUBCASE ID. SEE HOW LONG THE RECORD IS
 
 80 IF (iz(ncstm+16) == 0) GO TO 200
 
!     SUBCOM UNLESS IZ(16).LT.0. IN WHICH CASE IT IS A REPEAT SUBCASE
 
 IF (iz(ncstm+16) > 0) GO TO 90
 CALL bckrec (hccen)
 CALL bckrec (remfl)
 GO TO 200
 
!     SUBCOM. GET NUMBER OF FACTORS AND BCKREC THAT MANY RECORDS ON
!     HCFLD
 
!     OPEN CORE (AFTER NCSTM WORDS OF CSTM)
!             1 - NEXTZ   CASECC
!             NEXTZ+1 - NEXTZ+NROWS   COLUMN OF HCCEN
!             NEXTZ+NROWS+1 - NEXTZ+2*NROWS=NEXTP  HCCEN COMBINATION
!             NEXTP+1 - NEXTP+NELX  COLUMN OF REMFL
!             NEXTP+NELX+1 - NEXTP+2*NELX  REMFL COMBINATION
 
 90 CALL READ (*530,*100,casecc,z(ncstm+17),lcore,0,iwords)
 GO TO 550
 100 lcc  = iz(ncstm+166)
 lsym = iz(ncstm+lcc)
 DO  i = 1,lsym
   CALL bckrec (hccen)
   CALL bckrec (remfl)
 END DO
 nextz  = iwords + 16 + ncstm
 nrows2 = 2*nrows
 nextr  = nextz + nrows
 nelx2  = 2*nelx
 nall2  = nrows2 + nelx2
 nextp  = nextz  + nrows2
 isub   = nextp  + nelx
 IF (nextz+nall2 > lcore) GO TO 550
 
!     SET UP FOR SUBSEQ
 
 DO  i = 1,nall2
   z(nextz+i) = 0.
 END DO
 DO  i = 1,lsym
   coef = z(ncstm+lcc+i)
   IF (coef == 0.) GO TO 160
   nn = nrows
   CALL unpack (*140,hccen,z(nextz+1))
   DO  j = 1,nrows
     z(nextr+j) = z(nextr+j) + coef*z(nextz+j)
   END DO
   140 nn = nelx
   CALL unpack (*170,remfl,z(nextp+1))
   DO  j = 1,nelx
     z(isub+j) = z(isub+j) + coef*z(nextp+j)
   END DO
   CYCLE
   
!     COEF = 0.
   
   160 FILE = hccen
   CALL fwdrec (*530,hccen)
   FILE = remfl
   CALL fwdrec (*530,remfl)
 END DO
 
!     MOVE THE VECTOR IN CORE
 
 DO  i = 1,nrows
   z(ncstm+i) = z(nextr+i)
 END DO
 DO  i = 1,nelx
   z(nncr+i) = z(isub+i)
 END DO
 GO TO 260
 
!     NOT A SUBCOM
!     UNPACK A COLUMN OF HCFLD. FIRST SKIP TO NEXT RECORD ON CASECC
 
 200 FILE = casecc
 CALL fwdrec (*530,casecc)
 nn = nrows
 CALL unpack (*210,hccen,z(ncstm+1))
 GO TO 230
 210 DO  i = 1,nrows
   z(ncstm+i) = 0.
 END DO
 230 nn = nelx
 CALL unpack (*240,remfl,z(nncr+1))
 GO TO 260
 240 DO  i = 1,nelx
   z(nncr+i) = 0.
 END DO
 
!     HCFLD VECTOR IS IN Z(NCSTM+1)-Z(NCSTM+NROWS=NNCR) AND REMFL IS IN
!     Z(NNCR+1)-Z(NNCR+NELX). MATCH ELEMENT TYPE ON HOEF1 WITH ESTFLD
 
 260 FILE = estfld
 270 CALL READ (*530,*540,estfld,iel,1,0,iwords)
 iex = 3
 IF (iel == 66 .OR. iel == 67) iex = 63
 IF (iel == 65) iex = 27
 ipts = iex/3
 
!     SINCE IS2D8 HAS 9 POINTS ON HCCEN BUT ONLY ONE ON HEOF1 AND ESTFLD
!     RESET IPTS
 
 IF (iel == 80) ipts = 9
 IF (iel == eltype) GO TO 290
 
!     NO MATCH. SKIP TO NEXT RECORD, BUT KEEP UP WITH NCOUNT
 
 280 CALL READ (*530,*270,estfld,idum,2,0,iwords)
 CALL fread (estfld,idum,-(idum(2)+19+iex),0)
 ncount = ncount + 1
 hcount = hcount + ipts
 GO TO 280
 
!     MATCH ON ELEMENT TYPE. FIND A MATCH ON ELEMENT ID
 
 290 FILE = hoef1
 CALL READ (*530,*380,hoef1,rbuf,9,0,iwords)
 elid = ibuf(1)/10
 FILE = estfld
 
!     NEXT STATEMENT IS FOR ISOPARAMETRICS WHICH HAVE MULTIPLE POINTS
!     ON HOEF1, BUT ONLY ONE SET OF INFO ON ESTFLD(BUT MULTIPLE COORDS
!     FOR NON-BASIC COORDINATE SYSTEMS). IF MATERIAL IS ALLOWED TO BE
!     TEMPERATURE-DEPENDENT AT SOME LATER DATE IN MAGNETICS PROBLEMS,
!     THEN ESTFLD WILL HAVE MULTIPLE INFO. WRIITEN IN ESTMAG AND THIS
!     STATEMENT CAN BE DELETED
 
 IF (oldeid /= 0) GO TO 310
 
 300 CALL READ (*530,*540,estfld,iz(nall+1),2,0,iwords)
 ncount = ncount + 1
 hcount = hcount + ipts
 ielid  = iz(nall+1)
 ngrids = iz(nall+2)
 nwords = ngrids + 19 + iex
 IF (nall+nwords > lcore) GO TO 550
 CALL READ (*530,*540,estfld,z(nall+1),nwords,0,iwords)
 
 IF (elid == ielid) GO TO 310
 GO TO 300
 
!     MATCH ON ELEMENT ID. PICK UP HM FROM HOEF1(IN ELEMENT COORDS)
!     PICK UP 3 X 3 TRANSFORMATION MATRIX FROM ESTFLD TO CONVERT ELEMENT
!     SYSTEM TO BASIC. THEN MULTIPLY
 
 310 hm(1) = rbuf(4)
 hm(2) = rbuf(5)
 hm(3) = rbuf(6)
!WKBNB 8/94 ALPHA-VMS
 itype = numtyp( hm(2) )
 IF ( itype <= 1 ) hm(2) = 0.
 itype = numtyp( hm(3) )
 IF ( itype <= 1 ) hm(3) = 0.
!WKBNE 8/94 ALPHA-VMS
 CALL gmmats (z(nall+ngrids+10),3,3,0,hm,3,1,0,hmg)
 
!     PICK UP HC FROM HCCEN VECTOR. FOR ALL EXCEPT ISOPARAMETRICS,HCOUNT
!     POINTS TO THE Z COMPONENT OF PROPER HC WHICH STARTS AT Z(NCSTM+1)
 
 IF (rbuf(2) /= hex1 .AND. rbuf(2) /= hex2 .AND. rbuf(2) /= hex3) GO TO 330
 
!     ISOPARAMETRIC SOLIDS
 
 IF (oldeid == elid) GO TO 320
 oldeid = elid
 strspt = 0
 320 strspt = strspt + 1
 IF (strspt >= 21) oldeid = 0
 IF (rbuf(2) == hex1 .AND. strspt >= 9) oldeid = 0
 GO TO 340
 330 strspt = 1
 
!     NEXT LINE IS FOR IS2D8 WHICH HAS 9 POINTS ON HCCEN BUT ONE ON
!     ESTFLD
 
 IF (iel == 80) strspt = 9
 340 isub  = ncstm + 3*(hcount-ipts+strspt-1)
 hc(1) = z(isub+1)
 hc(2) = z(isub+2)
 hc(3) = z(isub+3)
 
 DO  i = 1,3
   rbuf(i+3) = hmg(i) + hc(i)
 END DO
 
!     TO GET INDUCTION B, MULTIPLY H BY MATERIALS
 
 CALL gmmats (z(nall+ngrids+1),3,3,0,rbuf(4),3,1,0,rbuf(7))
 
!     ADD IN REMANENCE Z(NNCR+1)-Z(NNCR+NELX)
 
 isub = nncr + 3*ncount - 3
 rbuf(7) = rbuf(7) + z(isub+1)
 rbuf(8) = rbuf(8) + z(isub+2)
 rbuf(9) = rbuf(9) + z(isub+3)
 
!     CHECK FOR REQUEST FOR NON-BASIC COORD. SYSTEM. TA TRANSFORMS TO
!     BASIC
 
 ifield = iz(nall+ngrids+19)
 IF (ifield == 0) GO TO 370
 icoord(1) = ifield
 
!     NEXT LINE IS FOR IS2D8 WHICH ONLY ONE POINT ON ESTFLD
 
 IF (iel == 80) strspt = 1
 isub = nall + ngrids + 19 + 3*strspt - 3
 coord(2) = z(isub+1)
 coord(3) = z(isub+2)
 coord(4) = z(isub+3)
 CALL transs (coord,ta)
 CALL gmmats (ta,3,3,1,rbuf(7),3,1,0,temp)
 DO  i = 1,3
   rbuf(i+6) = temp(i)
 END DO
 
 370 CONTINUE
 
!     WRITE OUT TO HOEH1
 
 CALL WRITE (hoeh1,rbuf,9,0)
 
!     GET ANOTHER ELEMENT OF THIS TYPE IN THIS SUBCASE
 
 GO TO 290
 
!     END OF ELEMENTS OF PRESENT TYPE AND/OR SUBCASE ON HOEF1
 
 380 CALL WRITE (hoeh1,0,0,1)
 FILE = estfld
 
!     SKIP RECORD BUT KEEP UP WITH NCOUNT
 
 390 CALL READ (*530,*400,estfld,idum,2,0,iwords)
 CALL fread (estfld,idum,-(idum(2)+19+iex),0)
 ncount = ncount + 1
 hcount = hcount + ipts
 GO TO 390
 400 FILE = hoef1
 GO TO 60
 
!     EOF ON HOEF1 - ALL DONE
 
 410 CALL CLOSE (casecc,1)
 CALL CLOSE (hccen,1)
 CALL CLOSE (estfld,1)
 CALL CLOSE (hoef1,1)
 CALL CLOSE (remfl,1)
 CALL CLOSE (hoeh1,1)
 mcb(1) = hoef1
 CALL rdtrl (mcb)
 mcb(1) = hoeh1
 CALL wrttrl (mcb)
 GO TO 600
 
!     FATAL ERROR MESSAGES
 
 500 WRITE  (otpe,510) sfm
 510 FORMAT (a25,', ROW COUNT ON HCCEN IN EMFLD IS NOT CONSISTENT')
 CALL mesage (-61,0,0)
 520 n = -1
 GO TO 560
 530 n = -2
 GO TO 560
 540 n = -3
 GO TO 560
 550 n = -8
 FILE = 0
 560 CALL mesage (n,FILE,nam)
 
 600 RETURN
END SUBROUTINE emfld
