SUBROUTINE cmrd2b (kode)
     
    !     THIS SUBROUTINE PROCESSES THE OLDMODES OPTION FLAG FOR THE CMRED2
    !     MODULE.
 
    !     INPUT  DATA
    !     GINO - LAMAMR - EIGENVALUE TABLE FOR SUBSTRUCTURE BEING REDUCED
    !            PHISSR - RIGHT HAND EIGENVECTOR MATRIX FOR SUBSTRUCTURE
    !                     BEING REDUCED
    !            PHISSL - LEFT HAND EIGENVECTOR MATRIX FOR SUBSTRUCTURE
    !                     BEING REDUCED
    !     SOF  - LAMS   - EIGENVALUE TABLE FOR ORIGINAL SUBSTRUCTURE
    !            PHIS   - RIGHT HAND EIGENVECTOR TABLE FOR ORIGINAL
    !                     SUBSTRUCTURE
    !            PHIL   - LEFT HAND EIGENVECTOR TABLE FOR ORIGINAL
    !                     SUBSTRUCTURE
 
    !     OUTPUT DATA
    !     GINO - LAMAMR - EIGENVALUE TABLE FOR SUBSTRUCTURE BEING REDUCED
    !            PHISS  - EIGENVECTOR MATRIX FOR SUBSTRUCTURE BEING REDUCED
    !     SOF  - LAMS   - EIGENVALUE TABLE FOR ORIGINAL SUBSTRUCTURE
    !            PHIS   - EIGENVECTOR MATRIX FOR ORIGINAL SUBSTRUCTURE
 
    !     PARAMETERS
    !     INPUT- GBUF   - GINO BUFFER
    !            INFILE - INPUT FILE NUMBERS
    !            ISCR   - SCRATCH FILE NUMBERS
    !            KORBGN - BEGINNING ADDRESS OF OPEN CORE
    !            OLDNAM - NAME OF SUBSTRUCTURE BEING REDUCED
    !            MODES  - OLDMODES OPTION FLAG
    !            NFOUND - NUMBER OF MODAL POINTS USED
    !            LAMAAP - BEGINNING ADDRESS OF LAMS RECORD TO BE APPENDED
    !            MODLEN - LENGTH OF MODE USE ARRAY
    !     OTHERS-LAMAMR - LAMAMR INPUT FILE NUMBER
    !            PHIS   - PHIS INPUT FILE NUMBER
    !            LAMS   - LAMS INPUT FILE NUMBER
    !            PHISS  - PHISS INPUT FILE NUMBER
 
 
    INTEGER, INTENT(IN OUT)                  :: kode
    LOGICAL :: modes
    INTEGER :: dry,gbuf1,oldnam,z,phissr,phissl,phisl,rgdfmt
    DIMENSION       rz(1),modnam(2),itmlst(3)
    CHARACTER (LEN=23) :: ufm
    COMMON /xmssg / ufm
    COMMON /BLANK / idum1,dry,idum7,gbuf1,idum2(5),infile(11),  &
        idum3(6),iscr(11),korlen,korbgn,oldnam(2),  &
        idum5(7),modes,idum6,lamaap,nfound,modlen
    COMMON /zzzzzz/ z(1)
    COMMON /system/ idum4,iprntr
    EQUIVALENCE     (rz(1),z(1)),(lamamr,infile(2)),  &
        (phissr,infile(3)),(phissl,infile(4)), (lams,iscr(5)),(phisl,iscr(6))
    DATA    modnam/ 4HCMRD,4H2B  /
    DATA    itmlst/ 4HPHIS,4HPHIL,4HLAMS/
    DATA    rgdfmt/ 3 /
 
    !     TEST OPERATION FLAG
 
    IF (dry == -2) RETURN
    IF (kode == 3) GO TO 20
 
    !     TEST OLDMODES OPTION FLAG
 
    IF (modes) GO TO 10
 
    !     STORE GINO PHISS(R,L) AS PHI(S,L) ON SOF
 
    ifile = phissr
    IF (kode == 2) ifile = phissl
    item = itmlst(kode)
    CALL mtrxo (ifile,oldnam,item,0,itest)
    IF (itest /= 3) GO TO 120
    RETURN
 
    !     READ SOF PHI(S,L) ONTO GINO PHI(S,L) SCRATCH FILES
 
10  item = itmlst(kode)
    CALL mtrxi (phisl,oldnam,item,0,itest)
    IF (itest /= 1) GO TO 120
 
    !     READ SOF LAMS ONTO GINO LAMS SCRATCH FILE
 
    CALL sfetch (oldnam,itmlst(3),1,itest)
    item = itmlst(3)
    IF (itest > 1) GO TO 120
    CALL gopen  (lams,z(gbuf1),1)
    CALL suread (z(korbgn),-1,nwdsrd,itest)
    CALL WRITE  (lams,z(korbgn),nwdsrd,1)
    CALL suread (z(korbgn),-1,nwdsrd,itest)
    CALL WRITE  (lams,z(korbgn),nwdsrd,1)
    CALL CLOSE  (lams,1)
 
    !     SWITCH FILE NUMBERS
 
    IF (kode == 1) phissr = phisl
    IF (kode == 2) phissl = phisl
    lamamr = lams
    RETURN
 
    !     STORE LAMAMR (TABLE) AS LAMS ON SOF
 
20  IF (modes) GO TO 60
    item = itmlst(3)
    CALL DELETE (oldnam,item,itest)
    IF (itest == 2 .OR. itest > 3) GO TO 120
    ifile = lamamr
    CALL gopen  (lamamr,z(gbuf1),0)
    CALL fwdrec (*100,lamamr)
    itest = 3
    CALL sfetch (oldnam,itmlst(3),2,itest)
    IF (itest /= 3) GO TO 120
    DO  i = 1, 2
        z(korbgn+i-1) = oldnam(i)
    END DO
    z(korbgn+2) = rgdfmt
    z(korbgn+3) = modlen
    CALL suwrt (z(korbgn),4,2)
    lamwds = modlen - 1
    rz(korbgn+6) = 0.0
    DO  i = 1,lamwds
        CALL READ  (*90,*100,lamamr,z(korbgn),6,0,nwds)
        CALL suwrt (z(korbgn),7,1)
    END DO
    CALL READ  (*90,*100,lamamr,z(korbgn),6,0,nwds)
    CALL CLOSE (lamamr,1)
    CALL suwrt (z(korbgn),7,2)
    CALL suwrt (z(lamaap),modlen,2)
    CALL suwrt (z(lamaap),0,3)
60 CONTINUE
   RETURN
 
   !     PROCESS SYSTEM FATAL ERRORS
 
90 imsg = -2
   GO TO 110
100 imsg = -3
110 CALL sofcls
   CALL mesage (imsg,ifile,modnam)
   RETURN
 
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
130 WRITE (iprntr,900) ufm,modnam,item,oldnam
   dry = -2
   RETURN
 
135 WRITE (iprntr,902) ufm,modnam,item,oldnam
   dry = -2
   RETURN
140 imsg = -1
   GO TO 170
150 imsg = -2
   GO TO 170
160 imsg = -3
170 CALL smsg (imsg,item,oldnam)
   RETURN
 
180 WRITE (iprntr,901) ufm,modnam,item,oldnam
   dry = -2
   RETURN
 
900 FORMAT (a23,' 6211, MODULE ',2A4,' - ITEM ',a4,  &
       ' OF SUBSTRUCTURE ',2A4,' HAS ALREADY BEEN WRITTEN.')
901 FORMAT (a23,' 6632, MODULE ',2A4,' - NASTRAN MATRIX FILE FOR I/O',  &
       ' OF SOF ITEM ',a4,', SUBSTRUCTURE ',2A4,', IS PURBED.')
902 FORMAT (a23,' 6215, MODULE ',2A4,' - ITEM ',a4,  &
       ' OF SUBSTRUCTURE ',2A4,' PSEUDO-EXISTS ONLY.')
 
END SUBROUTINE cmrd2b
