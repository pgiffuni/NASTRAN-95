SUBROUTINE trht1a (casexx,usetd,gptt,trl,ngroup)
     
!     TRHT1A INITIALIZES FOR TRHT MODULE
 
!     ITS TASK IS TO EXTRACT INITIAL CONDITION POINTS FROM CASEXX
!     AND TO PUT INITIAL STUFF ON ICR5
 
 
 INTEGER, INTENT(IN OUT)                  :: casexx
 INTEGER, INTENT(IN)                      :: usetd
 INTEGER, INTENT(IN)                      :: gptt
 INTEGER, INTENT(IN)                      :: trl
 INTEGER, INTENT(OUT)                     :: ngroup
 EXTERNAL        andf
 INTEGER :: sysbuf,iz(160),NAME(2), FILE,andf,two1,mcb(7),ia(1)
 COMMON /bitpos/ isk(11),iue,isk1(3),iud
 COMMON /two   / two1(32)
 COMMON /BLANK / x
 COMMON /system/ sysbuf
 COMMON /trhtx / ik(7),ib(7),icr1,icr2,icr3,icr4,iscr5
 COMMON /trdd1 / nlft1,dit1,nlftp1
 COMMON /zzzzzz/ z(1)
 COMMON /zblpkx/ a(4),ii
 COMMON /packx / it1,it2,ii1,jj1,incr
 EQUIVALENCE     (z(1),iz(1)), (a(1),ia(1))
 DATA    NAME  / 4HTRHT,4H1A  /
 
 
 nz = korsz(z)
 nx = nz
 ibuf1 = nz - sysbuf + 1
 nz = nz - sysbuf
 CALL gopen (casexx,iz(ibuf1),0)
 CALL fread (casexx,iz(1),166,1)
 CALL CLOSE (casexx,1)
 itstep = iz(38)
 nlftp1 = iz(160)
 intmp  = iz(9)
 inltmp = iz(8)
 
!     FIND STUFF ON TRL
 
 FILE = trl
 CALL OPEN (*200,trl,iz(ibuf1),0)
 CALL READ (*220,*10,trl,iz(1),nz,0,iflag)
 GO TO 230
 10 ns =  iz(3)
 CALL skprec (trl,ns)
 30 CALL READ (*240,*40,trl,iz(1),nz,0,iflag)
 GO TO 230
 40 IF (iz(1) /= itstep) GO TO 30
 
!     TSTEP STUFF FOUND
 
 CALL CLOSE (trl,1)
 ngroup = (iflag-1)/3
 
!     MOVE TSETP STUFF TO BOTTOM OF CURE
 
 nz = nx - iflag + 1
 igroup = nz + 1
 DO  i = 2,iflag
   k = igroup + i - 2
   iz(k) = iz(i)
 END DO
 ibuf1 = nz - sysbuf + 1
 ibuf2 = ibuf1 -sysbuf
 nz = ibuf2
 CALL gopen (iscr5,iz(ibuf1),1)
 CALL WRITE (iscr5,iz(igroup),iflag-1,1)
 FILE = usetd
 
!     BRING IN USETD
 
 CALL gopen (usetd,iz(ibuf2),0)
 CALL READ (*220,*60,usetd,iz(1),nz,1,lusetd)
 GO TO 230
 60 CALL CLOSE (usetd,1)
 
!     BUILD SIL TO SILD CONVERTER TABLE
 
 mskue = two1(iue)
 mskud = two1(iud)
 m = 1
 l = 0
 DO  i = 1,lusetd
   IF (andf(iz(i),mskue) /= 0) GO TO 65
   l = l + 1
   IF (andf(iz(i),mskud) == 0) GO TO 67
   iz(l) = m
   65 CONTINUE
   m = m + 1
   CYCLE
   67 iz(l) = 0
 END DO
 
!     FIND STUFF IN GPTT
 
 its = intmp
 CALL makmcb (mcb,iscr5,m-1,2,1)
 ns  = 0
 FILE = gptt
 CALL OPEN (*200,gptt,iz(ibuf2),0)
 
!     POSITION TO HEADER RECORD
 
 ival = nz - 2*l
 CALL READ (*220,*80,gptt,iz(l+1),ival,0,iflag)
 GO TO 230
 
!     PUT OUT TEMPS
 
 80 CONTINUE
 
!     DETERMINE NUMBER OF ELEMENT TEMP RECORDS TO SKIP.
 
 list = l + 3
 k = l + iflag
 82 nsk = iz(k)
 IF (nsk > 0) GO TO 84
 k = k - 3
 IF (k > list) GO TO 82
 
!     SET IPOS TO SKIP ELEMENT TEMP RECORDS AND DUPLICATE HEADER.
 
 84 ipos = -nsk
 mcb(2) = 0
 90 IF (its == 0) GO TO 170
 k = list
 100 IF (iz(k) ==  its) GO TO 110
 k = k + 3
 IF (k > l+iflag) CALL mesage (-31,its,NAME)
 GO TO 100
 
!     FOUND TEMP SET
 
 110 tdflt = 0.0
 IF (iz(k+1) /= -1) tdflt = z(k+1)
 m = l + iflag
 DO  i = 1,l
   j = m + i
   z(j) = tdflt
 END DO
 
!     RECORD NUMBER OF TEMP SET FOUND
 
 ns = iz(k+2)
 IF (ns == 0) GO TO 150
 
!     SKIP TO DESIRED RECORD
 
 132 IF (ns-ipos < 0) THEN
   GO TO   134
 ELSE IF (ns-ipos == 0) THEN
   GO TO   140
 ELSE
   GO TO   136
 END IF
 134 CALL bckrec (gptt)
 ipos = ipos - 1
 GO TO 132
 136 CALL fwdrec (*220,gptt)
 ipos = ipos + 1
 GO TO 132
 140 CALL READ (*220,*145,gptt,a,2,0,iflg)
 IF (ia(1) <= 0) GO TO 140
 j = ia(1) + m
 z (j) = a(2)
 GO TO 140
 145 ipos = ipos + 1
 
!     ALL SET UP OUTPUT
 
 150 inext = m + 1
 DO  i = 1,l
   j = m + i
   ii = iz(i) + m
   IF (ii == m) CYCLE
   IF (ii == inext) GO TO 155
   DO  k = inext,ii
     z(k)  = 0.0
   END DO
   155 z(ii) = z(j)
   inext = ii + 1
 END DO
 j = inext - (m+1)
 CALL WRITE (iscr5,z(m+1),j,0)
 170 CALL WRITE (iscr5,z(1),0,1)
 mcb(2) = mcb(2) + 1
 IF (mcb(2) == 2) GO TO 190
 its = inltmp
 GO TO 90
 
!     ALL DONE
 
 190 CALL CLOSE (iscr5,1)
 CALL CLOSE (gptt,1)
 CALL wrttrl (mcb)
 RETURN
 
!     ERROR MESAGES
 
 200 ip1 = -1
 210 CALL mesage (ip1,FILE,NAME)
 RETURN
 220 ip1 = -2
 GO TO 210
 230 ip1 = -8
 GO TO 210
 240 CALL mesage (-31,itstep,NAME)
 RETURN
END SUBROUTINE trht1a
