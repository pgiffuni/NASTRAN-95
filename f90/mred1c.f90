SUBROUTINE mred1c
     
!     THIS SUBROUTINE CONVERTS THE EQSS DATA AND BGSS DATA TO CORRESPOND
!     TO THE BOUNDARY DEGREES OF FREEDOM (UB) FOR THE MRED1 MODULE
 
!     INPUT DATA
!     SOF - BGSS - BASIC GRID POINT IDENTIFICATION TABLE
 
!     OUTPUT DATA
!     GINO - EQST - TEMPORARY EQSS DATA FILE
 
!     PARAMETERS
!     INPUT - GBUF1  - GINO BUFFER
!             KORLEN - LENGTH OF OPEN CORE
!             NEWNAM - NAME OF NEW SUBSTRUCTURE
!             RGRID  - FREEBODY MODE IDENTIFICATION NUMBERS (SET IN
!                      MRED1)
!                      RGRID(1) .EQ. GRID POINT IDENTIFICATION NUMBER
!                      RGRID(2) .EQ. NUMBER OF CONTRIBUTING SUBSTRUCTURE
!             NCSUBS - NUMBER OF CONTRIBUTING SUBSTRUCTURES
!             NAMEBS - BEGINNING ADDRESS OF BASIC SUBSTRUCTURE NAMES
!             EQSIND - BEGINNING ADDRESS OF EQSS GROUP ADDRESSES
!             NSLBGN - BEGINNING ADDRESS OF SIL DATA
!             NSIL   - NUMBER OF SIL GROUPS
!             LOCUST - BEGINNING ADDRESS OF USET ARRAY
 
 EXTERNAL        andf,orf
 LOGICAL :: bounds
 INTEGER :: oldnam,dry,gbuf1,eqsind,rgrid,z,sildof,ub,estdta,  &
     silind,estwrt,bitpat,eqst,eqstrl,andf,orf
 DIMENSION       bitpat(32),modnam(2),eqstrl(7)
 COMMON /BLANK / oldnam(2),dry,idum1(6),gbuf1,idum2(4),korlen,  &
     newnam(2),idum3(4),rgrid(2),idum4(4),ncsubs,  &
     namebs,eqsind,nslbgn,nsil,idum6(3),locust, idum7(4),bounds
 COMMON /zzzzzz/ z(1)
 COMMON /two   / itwo(32)
 COMMON /bitpos/ idum5(20),ub
 DATA    eqst  , nhbgss,modnam / 203,4HBGSS,4HMRED,4H1C  /
 
!     IF OLDBOUNDS OPTION, GET EQST TRAILER
 
 IF (dry == -2) RETURN
 eqstrl(1) = eqst
 IF (.NOT. bounds) GO TO 5
 CALL rdtrl (eqstrl)
 
!     GET SIL DOF AND DECODE
 
 5 newips = 0
 DO  i = 1,nsil
   sildof = nslbgn + ((2*i) - 1)
   icode  = z(sildof)
   CALL decode (icode,bitpat,nwdsd)
   
!     TEST FOR DOF REMAINING IN BOUNDARY SET
   
   ndof = 0
   kompnt = 0
   DO  j = 1,nwdsd
     k = locust + (z(sildof-1)-1) + (j-1)
     IF (andf(z(k),itwo(ub)) == 0) CYCLE
     k = 32 - bitpat(j)
     kompnt = orf(kompnt,itwo(k))
     ndof = ndof + 1
   END DO
   
!     SAVE NEW SIL DATA
   
   IF (ndof == 0) GO TO 20
   newips = newips + 1
   z(sildof-1) = (8*newips) + ndof
   z(sildof) = kompnt
   CYCLE
   
!     SIL DATA NOT NEEDED
   
   20 z(sildof-1) = -1
 END DO
 
!     WRITE EQSS GROUP 0 DATA ONTO TEMPORARY EQST TABLE
 
 CALL gopen (eqst,z(gbuf1),1)
 CALL WRITE (eqst,newnam,2,0)
 CALL WRITE (eqst,ncsubs,1,0)
 CALL WRITE (eqst,newips,1,0)
 nwds = eqsind - namebs
 CALL WRITE (eqst,z(namebs),nwds,1)
 eqstrl(2) = nwds + 4
 
!     WRITE REMAINING EQSS GROUP DATA ONTO TEMPORARY EQST TABLE
 
 eqstrl(3) = ncsubs
 DO  i = 1,ncsubs
   j = 2*(i-1)
   estdta = z(eqsind+j)
   nwds = z(eqsind+j+1)
   
!     TEST SUBSTRUCTURE COMPONENTS
   
   IF (nwds <= 0) GO TO 60
   DO  j = 1,nwds,3
     silind = nslbgn + (2*(z(estdta+j) - 1))
     IF (rgrid(1) <= 0) GO TO 40
     IF (i /= rgrid(2)) GO TO 40
     IF (rgrid(1) /= z(estdta+j-1)) GO TO 40
     rgrid(1) = z(estdta+j)
     40 IF (z(silind) == -1) CYCLE
     
!     REPLACE IP, SIL NUMBERS AND WRITE DATA
     
     estwrt = estdta + j
     z(estwrt  ) = z(silind)/8
     z(estwrt+1) = z(silind+1)
     CALL WRITE (eqst,z(estwrt-1),3,0)
   END DO
   60 CALL WRITE (eqst,0,0,1)
 END DO
 
!     REDUCE SIL ENTRIES AND STORE NEW SIL DATA AT Z(2*NSIL)
 
 ndof   = 1
 loindx = 0
 newsil = nslbgn + (2*nsil)
 IF ((newsil+(2*nsil)) >= korlen) GO TO 130
 DO  i = 1,nsil
   j = 2*(i-1)
   IF (z(nslbgn+j) == -1) CYCLE
   z(newsil+loindx  ) = ndof
   z(newsil+loindx+1) = z(nslbgn+j+1)
   ndof = ndof + andf(z(nslbgn+j),7)
   loindx = loindx + 2
 END DO
 
!     WRITE SIL DATA ONTO TEMPORARY EQST TABLE
 
 korbgn = namebs
 IF (loindx <= 0) CALL WRITE (eqst,0,0,1)
 IF (loindx > 0) CALL WRITE (eqst,z(newsil),loindx,1)
 eqstrl(4) = loindx
 
!     READ AND WRITE BGSS GROUP 0 DATA
 
 CALL sfetch (oldnam,nhbgss,1,itest)
 IF (itest == 3) GO TO 90
 IF (itest == 4) GO TO 100
 IF (itest == 5) GO TO 110
 CALL suread (z(korbgn),-1,nwdsrd,itest)
 z(korbgn  ) = oldnam(1)
 z(korbgn+1) = oldnam(2)
 nbgss = z(korbgn+2)
 z(korbgn+2) = loindx/2
 CALL WRITE (eqst,z(korbgn),3,1)
 
!     ELIMINATE BGSS DATA NOT REQUIRED
 
 i = 0
 eqstrl(5) = 0
 DO  j = 1,nbgss
   CALL suread (z(korbgn),4,nwdsrd,itest)
   IF (i > (2*nsil)) GO TO 80
   IF (z(nslbgn+i) == -1) GO TO 80
   CALL WRITE (eqst,z(korbgn),4,0)
   eqstrl(5) = eqstrl(5) + 4
   80 i = i + 2
 END DO
 CALL WRITE  (eqst,0,0,1)
 CALL wrttrl (eqstrl)
 
!     CLOSE EQST FILE
 
 CALL CLOSE (eqst,1)
 RETURN
 
!     PROCESS MODULE FATAL ERRORS
 
 90 imsg = -1
 GO TO 120
 100 imsg = -2
 GO TO 120
 110 imsg = -3
 120 CALL smsg (imsg,nhbgss,oldnam)
 RETURN
 
!     PROCESS SYSTEM FATAL ERRORS
 
 130 imsg = -8
 ifile = 0
 CALL sofcls
 CALL mesage (imsg,ifile,modnam)
 
 RETURN
END SUBROUTINE mred1c
