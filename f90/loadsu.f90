SUBROUTINE loadsu
     
!     LOADSU SETS UP LOAD INFOTMATION FOR PROLAT FROM NSLT.
!     Z(IST)IS THE STARTING POINT FOR OPEN CORE,Z(MCORE) IS THE LAST
!     AVAILABLE WORD, NTOT IS THE NUMBER OF WORDS PUT INTO OPEN CORE
!     BY THIS ROUTINE. LOAD IS THE LOAD ID.
 
 LOGICAL :: remfl
 INTEGER :: subcas,buf2,scr1,FILE,hest,bgpdt
 DIMENSION       nwords(19),mcb(7),iz(1),l(2),zl(2),nam(2)
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm
 COMMON /zzzzzz/ z(1)
 COMMON /system/ sysbuf,iout
 COMMON /biot  / ng1,ng2,ist,subcas,x1,y1,z1,x2,y2,z2,buf2,remfl,  &
     mcore,load,nslt,scr1,hest,ntot
 EQUIVALENCE     (z(1),iz(1)),(l(1),zl(1))
 DATA    nam   / 4HLOAD,4HSU  /
 DATA    nwords/ 6,6,4,4,6,6,2,5,5,6,6,7,12,10,10,19,38,7,5/
 
 bgpdt  = 103
 mcb(1) = bgpdt
 CALL rdtrl (mcb)
 nrowsp = mcb(2)
 mcb(1) = hest
 CALL rdtrl (mcb)
 nel    = mcb(2)
 nsimp  = 0
 FILE   = nslt
 CALL OPEN (*1001,nslt,z(buf2),0)
 CALL READ (*1002,*10,nslt,z(ist+1),mcore,0,iwords)
 GO TO 1008
 10 nloads = iwords-2
 
!     CHECK LOAD SELECTION AGAINST SIMPLE LOAD ID-S
 
 IF (nloads == 0) GO TO 35
 DO  i = 1,nloads
   IF (iz(ist+2+i) == load) GO TO 80
 END DO
 
!     NOT A SIMPLE LOAD-MUST BEA LOAD COMBINATION. SKIP NLOADS RECORDS
!     AND SEARCH FOR PROPER LOAD ID
 
 DO  i = 1,nloads
   CALL fwdrec (*1002,nslt)
 END DO
 
!     READ 2 WORDS AT A TIME -1,-1 SIGNIFIES END OF LOAD CARD
 
 35 iload = ist + iwords
 40 CALL READ (*1002,*500,nslt,l,2,0,iflag)
 IF (l(1) == load) GO TO 60
 
!     NO MATCH-SKIP TO -1-S
 
 50 CALL fread (nslt,l,2,0)
 IF (l(1) == -1 .AND. l(2) == -1) GO TO 40
 GO TO 50
 
!     MATCH
 
 60 alls  = zl(2)
 70 CALL fread (nslt,l,2,0)
 IF (l(1) == -1 .AND. l(2) == -1) GO TO 90
 nsimp = nsimp + 1
 IF (iload+2*nsimp > mcore) GO TO 1008
 isub  = 2*nsimp - 1
 z(iload+isub) = zl(1)
 iz(iload+isub+1) = l(2)
 GO TO 70
 
!     WE HAVE NSIMP SIMPLE LOADS. FOR ONE LOAD,SET PROPER PARAMETERS
 
 80 nsimp = 1
 alls  = 1.
 iload = ist + iwords
 z(iload+1) = 1.
 iz(iload+2) = load
 
!     FOR EACH SIMPLE LOAD, FIND PROPER LOAD ID AND THEN POSITION TO
!     PROPER LOAD RECORD IN NSLT
 
 90 ntot  = 0
 isimp = iload + 2*nsimp
 DO  ns = 1,nsimp
   
   isub   = iload + 2*ns - 1
   factor = z(isub)
   id     = iz(isub+1)
   ncards = 0
   CALL REWIND (nslt)
   i = 1
   IF (nloads == 0) GO TO 110
   DO  i = 1,nloads
     IF (id == iz(ist+2+i)) GO TO 110
   END DO
   GO TO 499
   
   110 DO  j = 1,i
     CALL fwdrec (*1002,nslt)
   END DO
   
   125 CALL READ  (*1002,*260,nslt,nobld,1,0,iflag)
   CALL fread (nslt,ido,1,0)
   IF (isimp+2 > mcore) GO TO 1008
   iz(isimp+1) = nobld
   iz(isimp+2) = ido
   isimp = isimp + 2
   ntot  = ntot + 2
   
!     SKIP NOBLD=-20. IF NOBLD=24(REMFLUX), STORE ONLY NOBLD AND IDO,
!     BUT SKIP REMFLUX INFO ON NSLT
   
   IF (nobld == -20) GO TO 250
   IF (nobld <=  19) GO TO 245
   ktype = nobld - 19
   SELECT CASE ( ktype )
     CASE (    1)
       GO TO 126
     CASE (    2)
       GO TO 127
     CASE (    3)
       GO TO 128
     CASE (    4)
       GO TO 129
     CASE (    5)
       GO TO 130
   END SELECT
   126 mwords = 3*nrowsp
   GO TO 140
   127 mwords = 12
   GO TO 140
   128 mwords = 48
   GO TO 140
   129 mwords = 9
   GO TO 140
   130 mwords = 3*nel
   mwords = -mwords
   GO TO 141
   
   140 IF(isimp+mwords*ido > mcore) GO TO 1008
   ntot = ntot + mwords*ido
   141 DO  j = 1,ido
     
!     NCARDS TELLS HOW MANY SIMPLE LOAD CARDS HAVE THE PRESENT FACTOR
!     APPLIED TO IT
     
     ncards = ncards + 1
     CALL fread (nslt,z(isimp+1),mwords,0)
     IF (nobld /= 24) isimp = isimp + mwords
   END DO
   
!     DONE WITH CARDS OF PRESENT TYPE-GET ANOTHER TYPE
   
   GO TO 125
   
!     TYPE=-20    SKIP IT
   
   250 CALL fread (nslt,z,-(3*nrowsp),0)
   GO TO 125
   
!     NOT A MAGNETICS TYPE OF LOAD. - SKIP IT
   
   245 WRITE  (iout,246) uwm,load
   246 FORMAT (a25,', IN FUNCTIONAL MODULE PROLATE, LOAD SET',i8, /5X,  &
       'CONTAINS A NONMAGNETIC LOAD TYPE. IT WILL BE IGNORED.')
   DO  i = 1,ido
     CALL fread (nslt,z,-nwords(nobld),0)
   END DO
   
!     EOR ON NSLT-DONE WITH THIS SIMPLE LOAD-GET ANOTHER SIMPLE LOAD
   
!     SUBSTITUTE IN OPEN CORE NCARDS FOR THE SIMPLE LOAD ID. WE NO
!     LONGER NEED THE ID, BUT WE MUST SAVE NCARDS
   
   260 CONTINUE
   iz(isub+1) = ncards
   
 END DO
 
!     DONE
 
!     STORE ALL THIS INFO BACK AT Z(IST) AS FOLLOWS
 
!     ALLS,NSIMP,(LOAD FACTOR,NCARDS) FOR EACH SIMPLE LOAD ID,
!     ALL LOAD INFO FOR EACH SIMPLE LOAD STARTING WITH NOBLD AND IDO
 
 z(ist+1)  = alls
 iz(ist+2) = nsimp
 ns2   = 2*nsimp
 DO  i = 1,ns2
   z(ist+2+i) = z(iload+i)
 END DO
 isub1 = ist + ns2 + 2
 isub2 = iload + 2*nsimp
 DO  i = 1,ntot
   z(isub1+i) = z(isub2+i)
 END DO
 ntot = ntot + 2*nsimp + 2
 CALL CLOSE (nslt,1)
 RETURN
 
 499 load = id
 500 WRITE  (iout,501) ufm,load
 501 FORMAT (a23,', CANNOT FIND LOAD',i8,' ON NSLT IN BIOTSV')
 CALL mesage (-61,0,0)
 
 1001 n =-1
 GO TO 1010
 1002 n =-2
 GO TO 1010
 1008 n =-8
 FILE = 0
 1010 CALL mesage (n,FILE,nam)
 RETURN
END SUBROUTINE loadsu
