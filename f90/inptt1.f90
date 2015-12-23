SUBROUTINE inptt1
     
!     READ DATA BLOCK(S) FROM A NASTRAN USER TAPE WHICH MUST BE SET UP.
 
!     CALL TO THIS MODULE IS
 
!     INPUTT1  /O1,O2,O3,O4,O5/V,N,P1/V,N,P2/V,N,P3/V,N,P4  $
 
!     PARAMETERS P1 AND P2 ARE INTEGER INPUT, P3 AND P4 ARE BCD
 
!               P1= 0, NO ACTION TAKEN BEFORE READ (DEFAULT)
!                 =+N, SKIP FORWARD N DATA BLOCKS BEFORE READ
!                 =-1, USER TAPE IS REWOUND BEFORE READ
!                 =-2, A NEW REEL IS MOUNTED BEFORE READ
!                 =-3, THE NAMES OF ALL DATA BLOCKS ON USER TAPE ARE
!                      PRINTED AND READ OCCURS AT BEGINNING OF TAPE
!                 =-4, AN OUTPUT TAPE IS TO BE DISMOUNTED
!                      AFTER AN END-OF-FILE MARK IS WRITTEN.
!                      A NEW INPUT REEL WILL THEN BE MOUNTED.
!                 =-5, SEARCH USER TAPE FOR FIRST VERSION OF DATA
!                      BLOCKS REQUESTED.
!                      IF ANY ARE NOT FOUND, A FATAL TERMINATION
!                      OCCURS.
!                 =-6, SEARCH USER TAPE FOR FINAL VERSION OF DATA
!                      BLOCKS REQUESTED.
!                      IF ANY ARE NOT FOUND, A FATAL TERMINATION
!                      OCCURS.
!                 =-7, SEARCH USER TAPE FOR FIRST VERSION OF DATA
!                      BLOCKS REQUESTED.
!                      IF ANY ARE NOT FOUND, A WARNING OCCURS.
!                 =-8, SEARCH USER TAPE FOR FINAL VERSION OF DATA
!                      BLOCKS REQUESTED.
!                      IF ANY ARE NOT FOUND, A WARNING OCCURS.
!                 =-9, REWIND AND UNLOAD USER TAPE
 
!               P2= 0, FILE NAME IS INPT
!                 = 1, FILE NAME IS INP1
!                 = 2, FILE NAME IS INP2
!                 = 3, FILE NAME IS INP3
!                 = 4, FILE NAME IS INP4
!                 = 5, FILE NAME IS INP5
!                 = 6, FILE NAME IS INP6
!                 = 7, FILE NAME IS INP7
!                 = 8, FILE NAME IS INP8
!                 = 9, FILE NAME IS INP9
!                 THE MPL DEFAULT VALUE FOR P2 IS 0
 
!               P3=    TAPE ID CODE FOR USER TAPE, AN ALPHANUMERIC
!                      VARIABLE WHOSE VALUE MUST MATCH A CORRESPONDING
!                      VALUE ON THE USER TAPE.
!                      THIS CHECK IS DEPENDENT ON THE VALUE OF
!                      P1 AS FOLLOWS..
!                       *P1*             *TAPE ID CHECKED*
!                        +N                     NO
!                         0                     NO
!                        -1                    YES
!                        -2                    YES (ON NEW REEL)
!                        -3                    YES (WARNING CHECK)
!                        -4                    YES (ON NEW REEL)
!                        -5                    YES
!                        -6                    YES
!                        -7                    YES
!                        -8                    YES
!                        -9                     NO
!                      THE MPL DEFAULT VALUE FOR P3 IS XXXXXXXX
 
 
 EXTERNAL        rshift,andf
 LOGICAL :: tapeup,tapbit
 INTEGER :: oubuf,output,p1,p2,p3,p4,zero,rshift,andf,NONE(2),  &
     trl(7),NAME(2),subnam(2),inn(10),out(5),namex(2),  &
     idhdr(7),idhdrx(7),p3x(2),nt(5,3),dx(3),tapcod(2), bcdbin(4)
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm,uim,sfm
 COMMON /machin/ mach
 COMMON /BLANK / p1,p2,p3(2),p4(2) /system/ ksystm(65)  &
        /zzzzzz/ x(1)
 EQUIVALENCE     (ksystm(1),nb  ), (ksystm( 2),nout),  &
                 (ksystm(9),nlpp), (ksystm(12),line)
 DATA    subnam/ 4HINPT, 4HT1  / , msc      / 4HMSC /
 DATA    out   / 201,202,203,204,205/, mask / 65535 /
 DATA    zero  , mone,mtwo,mtre,mfor/ 0,-1,-2,-3,-4 /,  &
         mfiv  , msix,mete,mnin     /-5,-6,-8,-9    /
 DATA    inn   / 4HINPT,4HINP1,4HINP2,4HINP3,4HINP4 ,  &
                 4HINP5,4HINP6,4HINP7,4HINP8,4HINP9 /
 DATA    idhdr / 4HNAST,4HRAN ,4HUSER,4H TAP,4HE ID,4H CDE,4HE - /
 DATA    bcdbin/ 4HBCD ,4H    ,4HBINA,4HRY          /
 DATA    none  / 4H (NO,4HNE) /, ipt1,ipt4/ 1H1,1H4 /
 
 
 iptx = ipt1
 IF (p4(1) == msc) GO TO 20
 GO TO 100
 
 
 ENTRY input1
!     ============
 
!     INPUT1 HANDELS MSC/OUTPUT1 DATA.
!     INPUT1 IS CALLED FROM INPTT1, WITH P4 = 'MSC', OR IT IS CALLED
!     FROM INPTT4
 
 20 iptx = ipt4
 IF (p3(1) == bcdbin(1) .AND. p3(2) == bcdbin(2)) GO TO 9918
 IF (p3(1) == bcdbin(3) .AND. p3(2) == bcdbin(4)) GO TO 9918
 WRITE  (nout,30) uim
 30 FORMAT (a29,'. INPUTT1 IS REQUESTED TO READ INPUT TAPE GENERATED',  &
     ' IN MSC/OUTPUT1 COMPATIBLE RECORDS')
 
 100 lcor = korsz(x) - 2*nb
 IF (lcor <= 0) CALL mesage (-8,lcor,subnam)
 inbuf = lcor  + 1
 oubuf = inbuf + nb
 tapcod(1) = p3(1)
 tapcod(2) = p3(2)
 IF (p2 < 0 .OR. p2 > 9) GO TO 9907
 in = inn(p2+1)
 IF (iptx == ipt4) WRITE (nout,110) uim,nb,in
 110 FORMAT (a29,', CURRENT NASTRAN BUFFER SIZE IS',i9,' WORDS', /5X,  &
     'SYNCHRONIZED BUFFSIZE IS REQUIRED IN CURRENT NASTRAN AND',  &
     ' THE VERSION THAT WROTE ',a4,' TAPE (OR FILE)', /5X, 3(4H====),/)
 ifile  = in
 IF (mach >= 5) GO TO 120
 tapeup = tapbit(in)
 IF (.NOT.tapeup ) GO TO 9909
 120 IF (p1 < mnin) GO TO 9908
 
 IF (p1 == mnin) GO TO 5000
 IF (p1 < mfor) GO TO 3000
 IF (p1 == mtre) GO TO 2000
 IF (p1 <= zero) GO TO 150
 
 CALL OPEN (*9901,in,x(inbuf),2)
 DO  i = 1,p1
   CALL READ (*9906,*9906,in,namex,2,0,nf)
   CALL skpfil (in,1)
 END DO
 GO TO 250
 
 150 IF (p1 /= mtwo .AND. p1 /= mfor) GO TO 190
 
!     P1 = -2 OR P1 = -4 IS ACCEPTABLE ONLY ON IBM OR UNIVAC
 
 IF (mach /= 2 .AND. mach /= 3) GO TO 9908
 
 iold = -p1/2
 CALL OPEN (*9901,in,x(inbuf),2)
 CALL tpswit (in,iold,1,tapcod)
 
 190 IF (p1 /= mone .AND. p1 /= mtwo .AND. p1 /= mfor) GO TO 230
 
!     OPEN USER TAPE TO READ WITH REWIND AHD TAPE ID CHECK
 
 IF (p1 /= mone .AND. p1 /= mtwo .AND. p1 /= mfor .AND.  &
     iptx == ipt4) GO TO 230
 CALL OPEN (*9901,in,x(inbuf),0)
 CALL READ (*9911,*9912,in,dx,3,0,nf)
 CALL READ (*9911,*9912,in,idhdrx,7,0,nf)
 DO  kf = 1,7
   IF (idhdrx(kf) /= idhdr(kf)) GO TO 9913
 END DO
 CALL READ (*9911,*9912,in,p3x,2,1,nf)
 IF (p3x(1) /= p3(1) .OR. p3x(2) /= p3(2)) GO TO 9910
 CALL skpfil (in,1)
 GO TO 250
 
!     OPEN USER TAPE TO READ WITHOUT REWIND AND NO TAPE ID CHECK
 
 230 CALL OPEN (*9901,in,x(inbuf),2)
 IF (iptx == ipt4) CALL fwdrec (*9912,ix)
 
 250 DO  i = 1,5
   output = out(i)
   trl(1) = output
   CALL rdtrl (trl)
   IF (trl(1) <= 0) CYCLE
   CALL fname (output,NAME)
   IF (NAME(1) == NONE(1) .AND. NAME(2) == NONE(2)) CYCLE
   
!     PASS FILE NAME HEADER RECORD
   
   CALL READ (*9904,*9905,in,namex,2,0,nf)
   
!     READ TRAILER RECORD, SIX WORDS (OR 3 WORDS, IPTX=4 ONLY)
   
   CALL READ (*9904,*300,in,trl(2),6,1,nf)
   GO TO 340
   
!     JUST A NOTE, FROM G.CHAN/UNISYS -
!     LEVEL 17.5 USED 2 RECORDS HERE FOR THE MATRIX NAME (2 BCD WORDS,
!     1ST RECORD) AND 7 TRAILER WORDS (2ND RECORD)
   
   300 IF (iptx /= ipt4 .OR. nf < 3) GO TO 9905
   trl(5) = trl(2)
   trl(6) = trl(3)
   trl(7) = trl(4)
   DO  j = 2,7
     j1 = j/2 + 4
     j2 = MOD(j-1,2)*16
     trl(j) = andf(rshift(trl(j1),j2),mask)
   END DO
   
!     OPEN OUTPUT DATA BLOCK TO WRITE WITH REWIND
   
   340 CALL OPEN (*9902,output,x(oubuf),1)
   
!     COPY CONTENTS OF USER TAPE ONTO OUTPUT DATA BLOCK, INCLUDING
!     FILE NAME IN RECORD 0
   
   CALL cpyfil (in,output,x,lcor,nf)
   
!     CLOSE OUTPUT DATA BLOCK WITH REWIND AND EOF
   
   CALL CLOSE (output,1)
   
!     WRITE TRAILER
   
   trl(1) = output
   CALL wrttrl (trl)
   CALL page2 (-3)
   WRITE  (nout,400) uim,NAME,in,namex
   400 FORMAT (a29,' 4105,     DATA BLOCK ',2A4,' RETRIEVED FROM USER ',  &
       'TAPE',a4, /5X,'NAME OF DATA BLOCK WHEN PLACED ON USER ', 'TAPE WAS ',2A4 )
   
 END DO
 
!     CLOSE NASTRAN USER TAPE WITHOUT REWIND
 
 CALL CLOSE (in,2)
 RETURN
 
!     OBTAIN LIST OF DATA BLOCKS ON USER TAPE.
 
 2000 CALL OPEN (*9901,in,x(inbuf),0)
 CALL READ (*9911,*9912,in,dx,3,0,nf)
 CALL READ (*9911,*9912,in,idhdrx,7,0,nf)
 DO  kf = 1,7
   IF (idhdrx(kf) /= idhdr(kf)) GO TO 9913
 END DO
 CALL READ (*9911,*9912,in,p3x,2,1,nf)
 IF (p3x(1) /= p3(1) .OR. p3x(2) /= p3(2)) GO TO 9914
 2006 CALL skpfil (in,1)
 kf = 0
 2007 CALL page1
 line = line + 5
 WRITE  (nout,2010) in
 2010 FORMAT (1H0,50X,a4,14H FILE contents ,/46X,4HFILE,18X,4HNAME/1H0)
 2020 CALL READ (*2050,*9915,in,namex,2,1,nf)
 CALL skpfil (in,1)
 kf = kf + 1
 line = line + 1
 WRITE (nout,2030) kf,namex
 2030 FORMAT (45X,i5,18X,2A4)
 IF (line - nlpp < 0) THEN
   GO TO  2020
 ELSE
   GO TO  2007
 END IF
 2050 CALL REWIND (in)
 CALL skpfil (in,1)
 GO TO 250
 
 
!     SEARCH MODE
 
 3000 CONTINUE
 
!     EXAMINE OUTPUT REQUESTS AND FILL NAME TABLE
 
 nnt = 0
 DO  i = 1,5
   output = out(i)
   trl(1) = output
   CALL rdtrl (trl)
   IF (trl(1) <= 0) GO TO 3020
   CALL fname (output,NAME)
   IF (iptx == ipt4 .AND. NAME(1) == NONE(1) .AND. NAME(2) == NONE(2)  &
       ) GO TO 3010
   nt(i,1) = 0
   nt(i,2) = NAME(1)
   nt(i,3) = NAME(2)
   nnt = nnt + 1
   CYCLE
   3010 nt(i,2) = NAME(1)
   nt(i,3) = NAME(2)
   3020 nt(i,1) = -1
 END DO
 
 IF (nnt > 0) GO TO 3070
 CALL page2 (-2)
 WRITE  (nout,3060) uwm,iptx
 3060 FORMAT (a25,' 4137,  ALL OUTPUT DATA BLOCKS FOR INPUTT',a1,  &
     ' ARE PURGED.')
 RETURN
 
!     CHECK TAPE ID LABEL.
 
 3070 CALL OPEN (*9901,in,x(inbuf),0)
 CALL READ (*9911,*9912,in,dx,3,0,nf)
 CALL READ (*9911,*9912,in,idhdrx,7,0,nf)
 DO  kf = 1,7
   IF (idhdrx(kf) /= idhdr(kf)) GO TO 9913
 END DO
 CALL READ (*9911,*9912,in,p3x,2,1,nf)
 IF (p3x(1) /= p3(1) .OR. p3x(2) /= p3(2)) GO TO 9910
 CALL skpfil (in,1)
 
 
!     BEGIN SEARCH OF TAPE.
 
 kf = 0
 3110 CALL READ (*3500,*9915,in,namex,2,0,nf)
 kf = kf + 1
 
 DO  i = 1,5
   NAME(1) = nt(i,2)
   NAME(2) = nt(i,3)
   IF (nt(i,1) < 0) CYCLE
   IF (NAME(1) /= namex(1) .OR. NAME(2) /= namex(2)) CYCLE
   nt(i,1) = nt(i,1) + 1
   IF (nt(i,1) == 1 .OR. p1 == msix .OR. p1 == mete) GO TO 3150
   CALL page2 (-3)
   WRITE  (nout,3140) uwm,NAME,kf,in
   3140 FORMAT (a25,' 4138,  DATA BLOCK ',2A4,' (DATA BLOCK COUNT =',i5,  &
       ') HAS PREVIOUSLY BEEN RETRIEVED FROM', /36X ,  &
       'USER TAPE ',a4,' AND WILL BE IGNORED.')
   GO TO 3205
   3150 CALL READ (*9904,*3160,in,trl(2),6,1,nf)
   GO TO 3180
   3160 IF (iptx /= ipt4 .OR. nf < 3) GO TO 9905
   trl(5) = trl(2)
   trl(6) = trl(3)
   trl(7) = trl(4)
   DO  j = 2,7
     j1 = j/2 + 4
     j2 = MOD(j-1,2)*16
     trl(j) = andf(rshift(trl(j1),j2),mask)
   END DO
   3180 output = out(i)
   CALL OPEN (*9902,output,x(oubuf),1)
   CALL cpyfil (in,output,x,lcor,nf)
   CALL CLOSE (output,1)
   trl(1) = output
   CALL wrttrl (trl)
   CALL page2 (-2)
   WRITE  (nout,3185) uim,NAME,in,kf
   3185 FORMAT (a29,' 4139, DATA BLOCK ',2A4,' RETRIEVED FROM USER TAPE ',  &
       a4,' (DATA BLOCK COUNT =',i5,1H))
   IF (nt(i,1) > 1) GO TO 3190
   nnt = nnt - 1
   GO TO 3210
   3190 WRITE  (nout,3195) uwm
   3195 FORMAT (a25,' 4140, SECONDARY VERSION OF DATA BLOCK HAS REPLACED',  &
       ' EARLIER ONE.')
   CALL page2 (-2)
   GO TO 3210
 END DO
 
 3205 CALL skpfil (in,1)
 3210 IF (nnt > 0 .OR. p1 == msix .OR. p1 == mete) GO TO 3110
 GO TO 3900
 
 3500 IF (nnt <= 0) GO TO 3900
 CALL page2 (-7)
 IF (p1 == mfiv .OR. p1 == msix) GO TO 9916
 WRITE  (nout,3510) uwm
 3510 FORMAT (a25,' 4141, ONE OR MORE DATA BLOCKS NOT FOUND ON USER ',  &
     'TAPE.')
 DO  i = 1,5
   IF (nt(i,1) /= 0) CYCLE
   WRITE  (nout,3520) nt(i,2),nt(i,3)
   3520 FORMAT (20X,21HNAME of DATA BLOCK = ,2A4)
 END DO
 IF (p1 == mfiv .OR. p1 == msix) GO TO 9995
 
 3900 CONTINUE
 CALL skpfil (in,-1)
 CALL CLOSE (in,2)
 RETURN
 
 5000 CONTINUE
 CALL unload (in)
 RETURN
 
!     ERRORS
 
 9901 WRITE  (nout,9951) sfm,iptx,in
 9951 FORMAT (a25,' 4107, MODULE INPTT',a1,' UNABLE TO OPEN NASTRAN ',  &
     'FILE ',a4,1H.)
 GO TO 9995
 
 9902 WRITE  (nout,9952) sfm,iptx,output
 9952 FORMAT (a25,' 4108, SUBROUTINE INPTT',a1,' UNABLE TO OPEN OUTPUT',  &
     ' DATA BLOCK',i5)
 GO TO 9995
 
 9904 CALL mesage (-2,ifile,subnam)
 
 9905 CALL mesage (-3,ifile,subnam)
 
 9906 WRITE  (nout,9956) ufm,iptx,p1,in,i
 9956 FORMAT (a22,' 4111, MODULE INPUTT',a1,' IS UNABLE TO SKIP FORWARD'  &
     ,       i10,' DATA BLOCKS ON PERMANENT NASTRAN FILE ',a4,1H., /5X,  &
     'NUMBER OF DATA BLOCKS SKIPPED =',i5)
 line = line + 1
 GO TO 9995
 
 9907 WRITE  (nout,9957) ufm,iptx,p2
 9957 FORMAT (a23,' 4112, MODULE INPUTT',a1,' - ILLEGAL VALUE FOR ',  &
     'SECOND PARAMETER =',i20)
 GO TO 9995
 
 9908 WRITE  (nout,9958) ufm,iptx,p1
 9958 FORMAT (a23,' 4113, MODULE INPUTT',a1,' - ILLEGAL VALUE FOR ',  &
     'FIRST PARAMETER =',i20)
 GO TO 9995
 
 9909 WRITE  (nout,9959) ufm,in
 9959 FORMAT (a23,' 4127, USER TAPE ',a4,' NOT SET UP.')
 GO TO 9995
 
 9910 WRITE  (nout,9960) ufm,p3x,iptx,p3
 9960 FORMAT (a23,' 4136, USER TAPE ID CODE -',2A4,'- DOES NOT MATCH ',  &
     'THIRD INPUTT',a1,' DMAP PARAMETER -',2A4,2H-.)
 GO TO 9995
 
 9911 WRITE  (nout,9961) ufm,iptx,in
 9961 FORMAT (a23,' 4132, MODULE INPUTT',a1,' - END-OF-FILE ENCOUNTERED'  &
     ,   ' WHILE ATTEMPTING TO READ TAPE ID CODE ON USER TAPE ',a4,1H.)
 GO TO 9995
 
 9912 WRITE  (nout,9962) ufm,iptx,in
 9962 FORMAT (a23,' 4133, MODULE INPUTT',a1,' - END-OF-RECORD ',  &
     'ENCOUNTERED WHILE ATTEMPTING TO READ TAPE ID CODE ON ', 'USER TAPE ',a4,1H.)
 GO TO 9995
 
 9913 WRITE  (nout,9963) ufm,iptx,idhdrx
 9963 FORMAT (a23,' 4134, MODULE INPUTT',a1,  &
     ' - ILLEGAL TAPE CODE HEADER = ',7A4)
 GO TO 9995
 
 9914 WRITE  (nout,9964) uwm,p3x,p3
 9964 FORMAT (a25,' 4135, USER TAPE ID CODE -',2A4,'- DOES NOT MATCH ',  &
     'THIRD INPUTT1 DMAP PARAMETER -',2A4,2H-.)
 line = line + 2
 GO TO 2006
 
 9915 WRITE  (nout,9965) sfm,iptx
 9965 FORMAT (a25,' 4106, MODULE INPUTT',a1,' - SHORT RECORD.')
 GO TO 9995
 
 9916 WRITE  (nout,9966) ufm
 9966 FORMAT (a23,' 4142, ONE OR MORE DATA BLOCKS NOT FOUND ON USER ',  &
     'TAPE',/)
 DO  i = 1,5
   IF (nt(i,1) /= 0) CYCLE
   WRITE (nout,9967) nt(i,2),nt(i,3)
   line = line + 1
 END DO
 9967 FORMAT (20X,'NAME OF DATA BLOCK = ',2A4)
 GO TO 9995
 
 9918 WRITE  (nout,9968) ufm,p3
 9968 FORMAT (a23,', ILLEGAL TAPE LABEL NAME -',2A4,'-  POSSIBLY ',  &
     'THE 4TH PARAMETER OF INPTT4 IS IN ERROR')
 
 
 9995 line = line + 2
 CALL mesage (-61,0,0)
 RETURN
 
END SUBROUTINE inptt1
