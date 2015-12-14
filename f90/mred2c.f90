SUBROUTINE mred2c (kode)
     
!     THIS SUBROUTINE PROCESSES THE OLDMODES OPTION FLAG FOR THE MRED2
!     MODULE.
 
!     INPUT DATA
!     GINO - LAMAMR   - EIGENVALUE TABLE  FOR SUBSTRUCTURE BEING REDUCED
!            PHISS    - EIGENVCTOR MATRIX FOR SUBSTRUCTURE BEING REDUCED
!     SOF  - LAMS     - EIGENVALUE  TABLE FOR ORIGINAL SUBSTRUCTURE
!            PHIS     - EIGENVCTOR  TABLE FOR ORIGINAL SUBSTRUCTURE
 
!     OUTPUT DATA
!     GINO - LAMAMR   - EIGENVALUE TABLE  FOR SUBSTRUCTURE BEING REDUCED
!            PHISS    - EIGENVCTOR MATRIX FOR SUBSTRUCTURE BEING REDUCED
!     SOF  - LAMS     - EIGENVALUE TABLE  FOR ORIGINAL SUBSTRUCTURE
!            PHIS     - EIGENVCTOR MATRIX FOR ORIGINAL SUBSTRUCTURE
 
!     PARAMETERS
!     INPUT  - GBUF   - GINO BUFFER
!              INFILE - INPUT FILE NUMBERS
!              ISCR   - SCRATCH FILE NUMBERS
!              KORBGN - BEGINNING ADDRESS OF OPEN CORE
!              OLDNAM - NAME OF SUBSTRUCTURE BEING REDUCED
!              MODES  - OLDMODES OPTION FLAG
!              LAMAAP - BEGINNING ADDRESS OF LAMS RECORD TO BE APPENDED
!              NFOUND - NUMBER OF MODAL POINTS USED
!              MODLEN - LENGTH OF MODE USE ARRAY
!     OTHERS - LAMAMR - LAMAMR INPUT FILE NUMBER
!              PHIS   - PHIS INPUT FILE NUMBER
!              LAMS   - LAMS INPUT FILE NUMBER
!              PHISS  - PHISS INPUT FILE NUMBER
 
 
 INTEGER, INTENT(IN OUT)                  :: kode
 LOGICAL :: modes
 INTEGER :: dry,gbuf1,oldnam,z,phis,phiss,rgdfmt
 DIMENSION       modnam(2),itmlst(2)
 COMMON /BLANK / idum1,dry,idum7,gbuf1,idum2(5),infile(12),  &
     idum3(6),iscr(10),korlen,korbgn,oldnam(2),  &
     idum5(9),modes,idum6,lamaap,nfound,modlen
 COMMON /zzzzzz/ z(1)
 COMMON /system/ idum4,iprntr
 EQUIVALENCE     (lamamr,infile(2)),(phis,infile(3)),  &
     (lams,iscr(5)),(phiss,iscr(6))
 DATA    modnam/ 4HMRED,4H2C  /
 DATA    itmlst/ 4HPHIS,4HLAMS/
 DATA    rgdfmt/ 3 /
 
!     TEST OPERATION FLAG
 
 IF (dry == -2) GO TO 200
 IF (kode > 1) GO TO 20
 
!     TEST OLDMODES OPTION FLAG
 
 IF (modes) GO TO 10
 
!     STORE GINO PHIS AS PHIS ON SOF
 
 ifile = phis
 CALL mtrxo (phis,oldnam,itmlst(1),0,itest)
 item  = itmlst(1)
 IF (itest /= 3) GO TO 120
 GO TO 200
 
!     READ SOF PHIS ONTO GINO PHIS SCRATCH FILE
 
 10 CALL mtrxi (phiss,oldnam,itmlst(1),0,itest)
 item = itmlst(1)
 IF (itest /= 1) GO TO 120
 
!     READ SOF LAMS ONTO GINO LAMAMR SCRATCH FILE
 
 CALL sfetch (oldnam,itmlst(2),1,itest)
 item = itmlst(2)
 IF (itest > 1) GO TO 120
 CALL gopen  (lams,z(gbuf1),1)
 CALL suread (z(korbgn),-1,nwdsrd,itest)
 CALL WRITE  (lams,z(korbgn),nwdsrd,1)
 CALL suread (z(korbgn),-1,nwdsrd,itest)
 CALL WRITE  (lams,z(korbgn),nwdsrd,1)
 CALL CLOSE  (lams,1)
 
!     SWITCH FILE NUMBERS
 
 phis   = phiss
 lamamr = lams
 GO TO 200
 
!     STORE LAMAMR (TABLE) AS LAMS ON SOF
 
 20 IF (modes) GO TO 70
 item = itmlst(2)
 CALL DELETE (oldnam,item,itest)
 IF (itest == 2 .OR. itest > 3) GO TO 120
 ifile = lamamr
 CALL gopen (lamamr,z(gbuf1),0)
 CALL fwdrec (*100,lamamr)
 itest = 3
 CALL sfetch (oldnam,itmlst(2),2,itest)
 IF (itest /= 3) GO TO 120
 DO  i = 1,2
   z(korbgn+i-1) = oldnam(i)
 END DO
 z(korbgn+2  ) = rgdfmt
 z(korbgn+3  ) = modlen
 CALL suwrt (z(korbgn),4,2)
 lamwds = modlen - 1
 IF (lamwds < 1) GO TO 55
 DO  i = 1,lamwds
   CALL READ  (*90,*100,lamamr,z(korbgn),7,0,nwds)
   CALL suwrt (z(korbgn),7,1)
 END DO
 55 CALL READ  (*90,*100,lamamr,z(korbgn),7,0,nwds)
 CALL CLOSE (lamamr,1)
 CALL suwrt (z(korbgn),7,2)
 IF (kode == 3) GO TO 60
 CALL suwrt (z(lamaap),modlen,2)
 CALL suwrt (z(lamaap),0,3)
 GO TO 70
 60 DO  i = 1,modlen
   z(korbgn+i-1) = 1
 END DO
 CALL suwrt (z(korbgn),modlen,2)
 CALL suwrt (z(korbgn),0,3)
 70 CONTINUE
 GO TO 200
 
!     PROCESS SYSTEM FATAL ERRORS
 
 90 imsg = -2
 GO TO 110
 100 imsg = -3
 110 CALL sofcls
 CALL mesage (imsg,ifile,modnam)
 GO TO 200
 
!     PROCESS MODULE FATAL ERRORS
 
 120 SELECT CASE ( itest )
   CASE (    1)
     GO TO 130
   CASE (    2)
     GO TO 135
   CASE (    3)
     GO TO 140
   CASE (    4)
     GO TO 150
   CASE (    5)
     GO TO 160
   CASE (    6)
     GO TO 180
 END SELECT
 130 imsg = -9
 GO TO 190
 135 ismg = -11
 GO TO 190
 140 imsg = -1
 GO TO 170
 150 imsg = -2
 GO TO 170
 160 imsg = -3
 170 CALL smsg (imsg,item,oldnam)
 GO TO 200
 180 imsg = -10
 190 dry  = -2
 CALL smsg1 (imsg,item,oldnam,modnam)
 200 RETURN
END SUBROUTINE mred2c
