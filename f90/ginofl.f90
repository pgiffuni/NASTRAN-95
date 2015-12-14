SUBROUTINE ginofl
     
!     ROUTINE FOR GINOFILE MODULE
 
!     MODULE GINOFILE WILL CAPTURE ONE SCRATCH FILE (301 THRU 309) OF
!     PREVIOUS DAMP MODULE, AND MAKE IT A LEGITIMATE GINO FILE.
!     THE SCRATCH FILE CAN BE A TABLE DATA BLOCK OR A MATRIX DATA BLOCK.
!     USE DMAP ALTER TO PLACE THIS MODULE IMMEDIATELY AFTER ANY NASTRAN
!     EXECUTABLE DMAP MODULE WHOSE SCRATCH FILE IS TO BE CAPTURED
 
!     IT IS USER'S RESPONSIBILITY TO SEE THAT NO FIAT TABLE RE-
!     ARRANGEMENT BY MODULE XSFA BETWEEN THIS GINOFILE MODULE AND THE
!     PREVIOUS INTENDED MODULE
 
!     GINOFILE  /OUTFL/C,N,P1/C,N,P2/C,N,P3  $
 
!     INPUT   FILE = NONE
!     OUTPUT  FILE = OUTFL, ANY UNIQUE NAME
!     SCRATCH FILE = 301
!     PARAMETERS   -
!             P1   = SCRATCH FILE NUMBER, 301,302,303,...,309
!                    (NO DEFAULT)
!             P2   = ADDITIONAL NUMBER OF RECORDS IN P1 FILE TO BE
!                    SKIPPED (NOT INCLUDING HEADER RECORD, WHETHER IT
!                    EXISTS OR NOT, DEFAULT = 0)
!             P3   = NO. OF RECORDS TO BE COPIED TO OUTPUT FILE OUTFL,
!                    STARTING FROM THE P2+1 RECORD, OR UP TO EOF RECORD
!                    (DEFAULT JJ=999999)
 
!     THIS GINOFILE MODULE SHOULD BE MAPPED IN ALL LINKS EXCEPT LINK1
 
!     WRITTEN BY G.CHAN/UNISYS, MAY 1988
!     DEFINITELY THIS IS NOT AN ORDINARY JOB FOR AMATEUR OR SEASONNED
!     PROGRAMMERS, MY FRIENDS
 
!     THE TRICKY PART OF THIS PROGRAM IS THAT GINOFILE MODULE USES ONLY
!     ONE OUTPUT FILE AND ONE SCRATCH FILE, WHICH IS 301
!     THE PROBLEMS HERE ARE (1) HOW TO CAPTURE OTHER SCRATCH FILE OF THE
!     PREVIOUS DMAP MODULE, SAY 303, WHILE ONLY 301 IS AVAILABLE. AND
!     (2) HOW TO CAPTURE SCRATCH1 FILE WHILE THE ORIGINAL 301 GINO DATA,
!     SUCH AS TRAILER, UNIT NUMBER, FILE INDEX ETC. ARE GONE. (THE
!     ORIGINAL 301 GINO DATA HAS BEEN ZEROED OUT TO GIVE ROOM FOR THE
!     NEW SCRATCH1 BEING ASSIGNED TO GINOFILE MODULE).
 
 IMPLICIT INTEGER (a-z)
 LOGICAL :: debug
 CHARACTER (LEN=6) :: tbmx,mxtb(2)
 DIMENSION        trl(7),fn(2),NAME(2),tchi(9)
 DIMENSION        itrl(7)
 
 COMMON /packx /  itypep, jtypep, irowp, jrowp, incrp
 COMMON /unpakx/  itypeu,         irowu, jrowu, incru
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg /  ufm,uwm,uim
 COMMON /system/  sysbuf,nout,skip(21),icfiat
 COMMON /BLANK /  p1,p2,p3
 COMMON /xfiat /  fiat(3)
 COMMON /xfist /  fist(2)
 COMMON /xsortx/  SAVE(6)
 COMMON /zzzzzz/  z(1)
 EQUIVALENCE      (scr,p1)
 DATA    NAME  /  4HGINO,4HFL  /,  scra, tch / 4HSCRA, 4H     /
!WKBR DATA    BLANK /  4H           /,  IZ2 / 2   /,OUTFL  / 201   /
 DATA    iz2 / 2   /,outfl  / 201   /
 DATA    tchi  /  4HTCH1,4HTCH2 ,  4HTCH3,4HTCH4,4HTCH5,4HTCH6,  &
     4HTCH7,4HTCH8 ,  4HTCH9    /
 DATA    mxtb  /  'MATRIX'      ,  'TABLE '  /,debug  /.false./
 
!     CHECK SCRATCH FILE PARAMETER
 
 IF (scr > 300 .AND. scr < 400) GO TO 20
 WRITE  (nout,10) uwm,scr
 10 FORMAT (a25,', SCRATCH FILE PARAMETER ERROR. GINOFILE ABORTED AND'  &
     ,      ' NO OUTPUT GINO FILE CREATED', /5X,'FIRST PARAMETER =',i5)
 GO TO 300
 20 IF (scr < 310) GO TO 40
 WRITE  (nout,30) ufm
 30 FORMAT (a23,', GINOFILE IS PROGRAMMED TO PROCESS ONLY THE FIRST ',  &
     '9 SCRATCH FILES')
 CALL mesage (-61,0,0)
 
!     SETUP CORE, BUFFERS, AND GINO OUTPUT FILE NAME
 
 40 kore  = korsz(z(1))
 ibuf1 = kore - sysbuf - 1
 ibuf2 = ibuf1- sysbuf
 kore  = ibuf2- 1
 CALL fname (outfl,fn)
 
!     RECAPTURE SCRATCH FILE NUMBER, TRAILER, AND INDEX POINTER IN FIAT
!     AND FIST
!     NOTE -
!     IT IS HERE THAT SCRATCH FILE IS LIMITED FROM 301 THRU 309
!     SCRATCH FILES 310 AND HIGHER MAY NOT HAVE UNIQUE 8-LETTER NAMES
!     IN ALL COMPUTERS.
 
 INDEX = 0
 k   = fiat(3)*icfiat - 2
 tch = tchi(scr-300)
 DO  i = 4,k,icfiat
   IF (fiat(i+1) == scra .AND. fiat(i+2) == tch) GO TO 70
 END DO
 WRITE  (nout,60) ufm,scr
 60 FORMAT (a23,', SCRATCH FILE',i4,' DOES NOT EXIST IN FIAT TABLE. ',  &
     'THIS ERROR MAY DUE TO', /5X,  &
     'USER ERROR, OR GINOFILE WAS PRECEDED BY XSFA MODULE')
 CALL mesage (-37,0,NAME)
 70 INDEX = i
 IF (debug) WRITE (6,80) INDEX
 80 FORMAT (5X,'INDEX =',i6)
 IF (scr /= 301) GO TO 90
 
!     IF SCRATCH FILE IS 301, THE TRAILER IN FIAT HAS BEEN INITIALIZED
!     TO ZEROS. MUST RECLAIM THE TRAILER FROM /XSORTX/, SAVED BY WRTTRL
 
!     THE LABEL COMMON /XSORTX/, DEFINED VIA SEMDBD AND AVAILABLE IN ALL
!     LINKS, IS ORIGNALLY USED ONLY BY XSORT2 ROUTINE WHICH WAS EXECUTED
!     IN EARLY LINK1. THUS IT IS SAFE TO SAVE THE SCRATCH 301 TRAILER IN
!     /XSORTX/. NOTE THE OTHER SCARTCH FILES 302 THRU 309 DO NOT HAVE
!     THIS PROBLEM
 
 fiat(INDEX+ 3) = SAVE(1)
 fiat(INDEX+ 4) = SAVE(2)
 fiat(INDEX+ 5) = SAVE(3)
 IF (icfiat == 8) GO TO 90
 fiat(INDEX+ 8) = SAVE(4)
 fiat(INDEX+ 9) = SAVE(5)
 fiat(INDEX+10) = SAVE(6)
 
!     LOCATE 301 IN FIST TABLE AND SWAP FIAT INDEX THAT POINTS TO THE
!     TARGET SCR FILE
 
 90 k = fist(2)*2+2
 DO  i=3,k,2
   IF (fist(i) == 301) GO TO 110
 END DO
 CALL mesage (-37,0,NAME)
 110 fisti  = i
 fisti1 = fist(i+1)
 fist(i+1) = INDEX-1
 IF (debug) WRITE (6,120) i,fist(i),fist(i+1),INDEX
 120 FORMAT (10X,' I,FIST(I),FIST(I+1),INDEX =',4I6)
 
!     NOW, WE CAN READ THE SCRATCH FILE TRAILER
 
 trl(1) = 301
 CALL rdtrl (trl(1))
 trl(1) = scr
 tbmx = mxtb(2)
 IF (trl(7) > 0) tbmx = mxtb(1)
 WRITE  (nout,130) uim,tch,trl,tbmx,fn
 130 FORMAT (a29,' FROM GINOFILE MODULE', /5X,'TRAILER OF SCRA',a4,  &
     ' FILE IN PREVIOUS MODULE = (',i3,1H),5I5,i8,  /5X,a6,  &
     ' CONTENTS OF THIS FILE WILL BE TRANSFERRED TO GINO FILE ', 2A4,/)
 
!     SWAP SCR AND SCRX (301) FILE
!     OPEN SCRS FILE, AND SKIP P2 RECORDS IF REQUESTED BY USER
!     (DEFAULT SKIP 1 HEADER RECORD IF IT EXISTS)
 
 trl2  = trl(2)
 FILE  = scr
 scrx  = 301
 CALL OPEN (*260,scrx,z(ibuf1),0)
 nwds = trl(5)
 IF (nwds == 3) nwds =2
 nwds = trl(3)*nwds
 itypeu = trl(5)
 irowu  = 1
 jrowu  = trl(3)
 incru  = 1
 itypep = itypeu
 jtypep = itypeu
 irowp  = 1
 jrowp  = trl(3)
 incrp  = 1
 itrl(1) = outfl
 itrl(2) = 0
 itrl(3) = trl(3)
 itrl(4) = trl(4)
 itrl(5) = trl(5)
 itrl(6) = 0
 itrl(7) = 0
 CALL rectyp (scrx, irctyp)
 IF (irctyp == 0) GO TO 135
 icrq = nwds -kore
 IF (icrq <= 0) GO TO 145
 CALL mesage (-8, scrx, NAME)
 135 CALL READ (*250,*140,scrx,z,2,1,k)
!WKBR  140 IF (Z(1).NE.SCRA .OR. Z(IZ2).NE.TCH) CALL BCKREC (SCRX,1)
 140 IF (z(1) /= scra .OR. z(iz2) /= tch) CALL bckrec (scrx)
 145 ncol = 0
 IF (p3 <= 0) p3 = 999999
 IF (p2 <= 0) GO TO 160
 DO  ii=1,p2
   CALL fwdrec (*250,scrx)
 END DO
 
!     OPEN OUTPUT GINO FILE AND WRITE A HEADER RECORD
 
 160 FILE = outfl
 CALL OPEN  (*260,outfl,z(ibuf2),1)
 CALL WRITE (outfl,fn,2,1)
 162 CALL rectyp (scrx, irctyp)
 IF (irctyp == 0) GO TO 170
 
!     PROCESS STRING-FORMATED RECORD HERE
 
 CALL unpack (*164, scrx, z)
 GO TO 168
 164 DO  l = 1, nwds
   z(l) = 0
 END DO
 168 CALL pack (z, outfl, itrl)
 GO TO 185
 
!     COPY SCRATCH FILE DATA DIRECTLY TO OUTPUT FILE
 
 170 CALL READ  (*190,*180,scrx,z,kore,0,k)
 CALL WRITE (outfl,z,kore,0)
 GO TO 170
 180 CALL WRITE (outfl,z,k,1)
 185 ncol = ncol + 1
 IF (ncol < p3) GO TO 162
 
!     ALL DONE, CLOSE ALL FILES, WRITE TRAILER, AND ISSUE FRIENDLY
!     MESSAGES
 
 190 CALL CLOSE (scrx ,1)
 CALL CLOSE (outfl,1)
 trl(1) = outfl
 trl(2) = ncol
 IF (ncol > itrl(2)) CALL wrttrl (trl)
 IF (ncol == itrl(2)) CALL wrttrl (itrl)
 WRITE  (nout,200) uim,tch,fn
 200 FORMAT (a29,', DATA TRANSFER FROM PREVIOUS SCRA',a4,' FILE TO ',  &
     2A4,' IS ACCOMPLISHED')
 IF (p2 >      0) WRITE (nout,210) tch,p2
 IF (p3 < 999999) WRITE (nout,220) p3
 210 FORMAT (5X,'FIRST',i5,' RECORDS IN SCRA',a4,' FILE WERE SKIPPED ',  &
     'BEFORE DATA TRANSFER')
 220 FORMAT (5X,'LAST RECORD COPIED WAS RECORD NO.',i5)
 WRITE  (nout,230) fn,trl
 230 FORMAT (5X,'TRAILER OF THE NEW GINO FILE ',2A4,'  = (',i3,1H), 5I5,i8)
 
!     IF SCRATCH FILE CONTAINS MATRIX DATA, CHECK NO. OF COLUMNS
 
 IF (tbmx == mxtb(1) .AND. ncol /= trl2 .AND.  &
     (p2 == 0 .AND. p3 == 999999)) WRITE (nout,240) uim,tch,fn
 240 FORMAT (a29,', POSSIBLE ERROR IN GINOFILE WAS DETECTED', /5X,  &
     'NUMBERS OF COLUMNS IN INPUT FILE SCAR',a4,  &
     ' AND OUTPUT FILE ',2A4,' DISAGREE',//)
 
!     RESET FIST ORIGINAL INDEX FOR SCRATCH FILE 301
 
 fist(fisti+1) = fisti1
 GO TO 300
 
!     ERRORS
 
 250 k =-2
 GO TO 270
 260 k = -1
 270 CALL mesage (k,FILE,NAME)
 
 300 RETURN
END SUBROUTINE ginofl
