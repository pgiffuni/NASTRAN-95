SUBROUTINE xscndm
     
!     THE PURPOSE OF THIS ROUTINE IS TO RETURN TO THE CALLING PROGRAM
!     THE NEXT BCD OR BINARY ENTRY IN DMAP ARRAY.
 
!     IBUFF  = BUFFER AREA WHERE CARD IMAGE IS STORED FOR XRCARD INPUT.
!     IDLMTR = TABLE OF DELIMITER CHARACTERS
!     ITYPE  = TABLE FOR CONVERTING NUMBER TYPE TO WORD LENGTH.
 
!     LAST REVISED BY G.CHAN/UNISYS, 2/90
!     REMOVING LVAX AND .NOT.LVAX AND STANDARDIZED ALL BYTE OPERATIONS
 
 IMPLICIT INTEGER (a-z)
 EXTERNAL        lshift,rshift,orf
 INTEGER :: gnobuf(1),itype(6),dmpcrd(1),idlmtr(8),os(5), oscar(1)
 COMMON /system/ ksystm(100)
 COMMON /xgpic / icold,islsh,iequl,nblank,nxequi,  &
     ndiag,nsol,ndmap,nestm1,nestm2,nexit,  &
     nbegin,nend,njump,ncond,nrept,ntime,nsave,noutpt,  &
     nchkpt,npurge,nequiv,nacpw,nbpc,nawpc, maskhi,masklo,isgnon,nosgn,iallon
 COMMON /xgpie / nscr
 COMMON /zzzzzz/ core(1)
 COMMON /xgpi4 / irturn,insert,iseqn,dmpcnt,  &
     idmpnt,dmppnt,bcdcnt,length,icrdtp,ICHAR,newcrd, modidx,ldmap,isavdw,dmap(1)
 COMMON /xgpi5 / iapp,start,alter(2)
 COMMON /xgpi6 / krud(6),diag14,diag17,diag4,diag25,ifirst, ibuff(20)
 COMMON /passer/ istopf,modnam,kkcomm
 EQUIVALENCE     (ksystm(3),nogo),(core(1),os(1),loscar),  &
     (os(2),osprc),(os(3),osbot),(os(4),ospnt),  &
     (os(5),oscar(1),dmpcrd(1),gnobuf(1))
 DATA    itype / 1,1,2,2,2,4/ ,idlmtr/  &
     4H$   , 4H/   ,4H=   ,4H,   ,4H(   ,4H)   ,4H    ,4H*   /
 DATA    noscr1/ 4HOSCA/,      noscr2/4HR   /
 DATA    npt   / 4HNPTP/,      izero / 0    /
 DATA    nwpc  / 18    /,      ncpw  / 4    /
 
! *** WARNING - NWPC AFFECTS CODE IN XOSGEN SO BEWARE IF YOU CHANGE IT.
 
 kcomma = khrfn1(izero,1,idlmtr(4),1)
 kblank = khrfn1(izero,1,idlmtr(7),1)
 kkcomm = 0
 
!     CHECK FOR OSCAR TABLE OVERFLOW
 
 IF (oscar(osbot)+osbot > icrdtp) GO TO 310
 
!     CHECK FOR CARD READ ERROR
 
 IF (nogo == 2) GO TO 340
 
!     CHECK FOR NEW CARD NEEDED.
 
 IF (newcrd /= 0) GO TO 200
 IF (bcdcnt < 0.0) THEN
   GO TO   330
 ELSE IF (bcdcnt == 0.0) THEN
   GO TO    10
 ELSE
   GO TO   130
 END IF
 
!     BCDCNT = 0, TEST MODE
 
 10 IF (modnam == 0) GO TO 90
 kfl1 = 0
 icom = 0
 DO  kh  = 1,nwpc
   DO  kdh = 1,ncpw
     nchar  = khrfn1(izero,1,ibuff(kh),kdh)
     IF (nchar-kblank == 0) THEN
       GO TO    20
     ELSE
       GO TO    40
     END IF
     20 IF (kfl1 == 0) THEN
       GO TO    70
     END IF
     30 kfl1 = 2
     CYCLE
     40 IF (nchar-kcomma == 0) THEN
       GO TO    60
     END IF
     50 IF (icom == 1 .OR. kfl1 == 2) EXIT
     kfl1 = 1
     CYCLE
     60 kfl1 = 2
     icom = icom + 1
     IF (icom /= 2) CYCLE
     kkcomm = 1
     EXIT
     70 CONTINUE
   END DO
 END DO
 90 IF (dmap(idmpnt) == rshift(iallon,1)) GO TO 180
 IF (dmap(idmpnt) < 0.0) THEN
   GO TO   100
 ELSE IF (dmap(idmpnt) == 0.0) THEN
   GO TO   110
 ELSE
   GO TO   120
 END IF
 
!     BINARY VALUE - TRANSLATE TYPE INTO LENGTH
 
 100 i = IABS(dmap(idmpnt))
 IF (i > 6) GO TO 330
 
!     A MISUNDERSTANDING MAKES THE FOLLOWING STATEMENT NECESSARY.
 
 dmap(idmpnt) = orf(isgnon,i)
 length = itype(i)
 dmppnt = idmpnt
 idmpnt = length + 1 + idmpnt
 irturn = 3
 GO TO 350
 
!     CONTINUE MODE - GET NEXT CARD
 
 110 newcrd = 1
 GO TO 200
 
!     MODE IS BCD, INITIALIZE BCDCNT, DMPPNT, AND CHECK FOR OVERFLOW
 
 120 bcdcnt = dmap(idmpnt)
 idmpnt = idmpnt + 1
 IF (2*bcdcnt+idmpnt > ldmap) GO TO 330
 
!     TEST FOR OPERATOR ENTRY.
 
 130 irturn = 2
 IF (dmap(idmpnt) == iallon) GO TO 150
 140 dmppnt = idmpnt
 idmpnt = idmpnt + 2
 bcdcnt = bcdcnt - 1
 GO TO 350
 
!     DELIMITER FOUND - CHECK FOR COMPLEX NUMBER
 
 150 irturn = 1
 IF (khrfn1(izero,1,dmap(idmpnt+1),1) /=  &
     khrfn1(izero,1,idlmtr(5),1)) GO TO 140
 
!     LEFT PAREN FOUND - SEE IF TWO NUMBERS FOLLOW
 
 IF (dmap(idmpnt+2) == -2 .AND. dmap(idmpnt+4) == -2) GO TO 160
 IF (dmap(idmpnt+2) /= -4 .OR.  dmap(idmpnt+5) /= -4) GO TO 140
 
!     DOUBLE PRECISION COMPLEX NUMBER FOUND - FORM NUMBER CORRECTLY AND
!     SET TYPE CODE.
 
 dmap(idmpnt+5) = dmap(idmpnt+4)
 dmap(idmpnt+4) = dmap(idmpnt+3)
 dmap(idmpnt+3) = -6
 GO TO 170
 
!     SINGLE PRECISION COMPLEX NUMBER FOUND - FORM NUMBER CORRECTLY
 
 160 dmap(idmpnt+4) = dmap(idmpnt+3)
 dmap(idmpnt+3) = -5
 170 bcdcnt = 0
 idmpnt = idmpnt + 3
 GO TO 100
 
!     END OF DMAP INSTRUCTION
 
 180 irturn = 4
 GO TO 350
 
!     GET NEXT CARD IMAGE AND TRANSLATE INTO DMAP ARRAY.
 
 200 ibufct = 1
 ibwrd  = 1
 icall  = 0
 
!     CHECK FOR INSERT TO BE MADE
 
 IF (insert > 0 .OR. insert == -1) GO TO 210
 GO TO 250
 
!     GET NEXT CARD IMAGE FROM ALTER FILE
 
 210 CONTINUE
 CALL READ (*230,*220,npt,ibuff,18,1,l)
 GO TO 260
 
!     NO MORE INSTRUCTIONS TO INSERT FOR THIS ALTER
!     MOVE NEXT ALTER CONTROL TO ALTER CELLS
 
 220 alter(1) = ibuff(1)
 alter(2) = ibuff(2)
 GO TO 240
 
!     END OF ALTER FILE - SET ALTER CELL INFINITE
 
 230 alter(1) = 10000
 240 CONTINUE
 IF (newcrd > 0) GO TO 300
 GO TO 180
 
!     FILL IBUFF WITH CARD IMAGE
 
 250 CALL READ (*320,*260,nscr,ibuff,nwpc,0,lx)
 
!     CHECK INSERT FOR NO PRINT
 
 260 IF (insert < 0) GO TO 270
 
!     PRINTOUT DMAP INSTRUCTION
 
 IF (ifirst == 0) GO TO 270
 IF (diag17 == 0 .AND. (diag14 == 0 .OR. diag14 >= 10)) GO TO 270
 i = 5
 IF (newcrd > 0) i = 6
 CALL xgpimw (i,nwpc,dmpcnt,ibuff)
 
!     CHECK FOR COMMENT CARD
 
 270 IF (khrfn1(izero,1,idlmtr(1),1) == khrfn1(izero,1,ibuff(1),1)) GO TO 200
 
!     CONVERT CARD IMAGE
 
 CALL xrcard (dmap,ldmap,ibuff)
 
!     CHECK FOR BAD CARD FORMAT
 
 IF (dmap(1) == 0) GO TO 180
 
!     TRANSLATE CARD IMAGE INTO DMAP ARRAY
 
 idmpnt = 1
 bcdcnt = 0
 newcrd = 0
 GO TO 10
 
!     DIAGNOSTIC MESSAGES -
 
!     ERROR IN ALTER DECK - CANNOT FIND LOGICAL END OF CARD
 
 300 CALL xgpidg (40,0,0,0)
 GO TO 180
 
!     OSCAR TABLE OVERFLOW
 
 310 CALL xgpidg (14,noscr1,noscr2,dmpcnt)
 CALL xgpidg (-38,2000,0,0)
 
!     THIS DMAP INSTRUCTION NOT FOLLOWED BY END CARD.
 
 320 CALL xgpidg (44,ospnt,0,0)
 GO TO 340
 
!     CANNOT INTERPRET DMAP CARD
 
 330 CALL xgpidg (34,0,dmpcnt,0)
 
!     ABORT - CANNOT CONTINUE COMPILATION
 
 340 nogo   = 2
 irturn = 5
 350 RETURN
END SUBROUTINE xscndm
