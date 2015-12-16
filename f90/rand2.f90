SUBROUTINE rand2 (FILE,ilist,load,IF,LEN,llist)
     
!     READS FILE UNTIL IT FINDS DATA RECORD IN LIST - RETURNS LOAD,
!     IF, AND LEN
 
 
 INTEGER, INTENT(IN)                      :: FILE
 INTEGER, INTENT(IN OUT)                  :: ilist(2)
 INTEGER, INTENT(OUT)                     :: load
 INTEGER, INTENT(IN OUT)                  :: IF
 INTEGER, INTENT(OUT)                     :: LEN
 INTEGER, INTENT(IN)                      :: llist
 INTEGER :: idr(10),NAME(2), mid(2,10),itemp(5),  &
     idata(50),ifmt(2,84),ifmtt(11),data1(100),DATA(1),  &
     itb(180),itb1(137),itb2(145),filex
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm
 COMMON /system/ ibuf,nout
 EQUIVALENCE     (itb1(1),itb(1)), (itb2(1),itb(138))
 DATA    ifmtt / 1,11,41,55,61,99,121,155,181,199,237/
 DATA    ioldld/ 0 /
 DATA    ifmt  / 1, 1,    -1,-1,     1, 1,     1, 1,     1, 1,     6, 2,  &
     6, 2,     6, 2,     0, 3,     1, 1,     4, 4,     4, 4,  &
     4, 4,     4, 0,     6, 2,     0, 3,     6, 2,     6, 2,  &
     6, 2,    -1,-1,    -1,-1,    -1,-1,    -1,-1,    -1,-1,  &
     -1,-1,    -1,-1,    -1,-1,    -1,-1,    -1,-1,    -1,-1,  &
     -1,-1,    -1,-1,    -1,-1,     7, 5,    -1,-1,    -1,-1,  &
     -1,-1,    -1,-1,     0, 8,     0, 8,     0, 8,     0, 8,  &
     -1,-1,    -1,-1,    -1,-1,    -1,-1,     0, 9,     0,10,  &
     0,11,     0, 6,     0, 8,    -1,-1,    -1,-1,    -1,-1,  &
     -1,-1,    -1,-1,    -1,-1,    -1,-1,    -1,-1,    -1,-1,  &
     -1,-1,     0, 3,     0, 3,     7, 2,    -1,-1,    -1,-1,  &
     -1,-1,    -1,-1,    -1,-1,    -1,-1,    -1,-1,    -1,-1,  &
     -1,-1,    -1,-1,    -1,-1,    -1,-1,    -1,-1,    -1,-1,  &
     -1,-1,    -1,-1,    -1,-1,    -1,-1,     7, 2,    -1,-1/
 
!     IFMT TABLE (ELEM ID IN GPTABD ORDER) HAS 2 WORDS PER ENTRY
!        WORD1  FORCE  FORMAT POINTER INTO IFMTT TABLE
!        WORD2  STRESS FORMAT POINTER INTO IFMTT TABLE
 
!     IFMTT TABLE  HAS ONE ENTRY PER FORMAT TYPE
!        THE ENTRY IS THE BEGINNING OF THE FORMAT IN THE ITB TABLE
 
 DATA itb1/ 6,   3,   5,   4,   6,  &
     0,   1,   1,   2,   2,  &
     16,   3,   4,  12,   5,  13,   6,  14,   7,   8,  16,   9,  17, 10,  18,  &
     0,   1,   2,   2,   3,   3,   4,   4,   5,   6,   6,   7,   7, 8,   8,  &
     8,   3,   6,   4,   7,   5,   8, 0,   1,   1,   2,   2,   3,   3,  &
     4,   3,   4, 0,   1,   1,  &
     20,   3,   4,   5,   6,   7,  12,  13,  14,  15,  16,   8,   9,  &
     10,  11,  17,  18,  19,  20,  &
     0,   1,   2,   3,   4,   5,   1,   2,   3,   4,   5,   6,   7,  &
     8,   9,   6,   7,   8,   9,  &
     12,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12,  &
     0,   1,   2,   3,   4,   5,   1,   2,   3,   4,   5,  &
     18,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12,  13,  14,  &
     15,  16,  17,  18/
 DATA itb2/ 0,   1,   2,   3,   4,   5,   6,   7,   8,   1,   2,   3,   4,  &
     5,   6,   7,   8,  &
     14,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12,  13,  14,  &
     0,   1,   2,   3,   4,   5,   6,   1,   2,   3,   4,   5,   6,  &
     10,   3,   4,   5,   6,   7,   8,   9,  10,  &
     0,   1,   2,   3,   4,   1,   2,   3,   4,  &
     20,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12,  13,  14,  &
     15,  16,  17,  18,  19,  20,  &
     0,   1,   2,   3,   4,   5,   6,   7,   8,   9,   1,   2,   3,  &
     4,   5,   6,   7,   8,   9,  &
     24,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12,  13,  14,  &
     15,  16,  17,  18,  19,  20,  21,  22,  23,  24,  &
     0,   1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11,   1,  &
     2,   3,   4,   5,   6,   7,   8,   9,  10,  11/
 DATA NAME/ 4HRAND,4H2   /
 DATA mid / 3001  ,4HDISP, 3010  ,4HVELO,  &
     3011  ,4HACCE, 3002  ,4HLOAD,  &
     3003  ,4HSPCF, 3004  ,4HELFO,  &
     3005  ,4HSTRE, 3015  ,4HDISP,  &
     3016  ,4HVELO, 3017  ,4HACCE/
 
 filex = FILE
 
!     POSITION TO + READ ID RECORD
 
 5 CALL fwdrec (*950,FILE)
 CALL READ (*950,*920,FILE,idr,10,1,i)
 idr(5) = idr(5)/10
 
!     IDR(5) = 10*ELEM.ID + DEV.CODE
!     IDR(2) = GINO FILE 3004, 3005 ETC.
!     CONVERT MAJOR ID TO MNEMONIC
 
 DO  i = 1,10
   IF (idr(2) == mid(1,i)) GO TO 20
 END DO
 
!     ILLEGAL FORMAT
 
 GO TO 970
 
!     CHECK FOR  MID
 
 20 IF (ilist(1) /= mid(2,i)) GO TO 5
 ielem = i
 
!     LOOK FOR ID IN LIST
 
 DO  i = 1,llist,5
   IF (idr(5)-ilist(i+1) < 0) THEN
     GO TO     5
   ELSE IF (idr(5)-ilist(i+1) == 0) THEN
     GO TO    50
   ELSE
     GO TO    30
   END IF
   30 CONTINUE
 END DO
 GO TO 5
 
!     ID IS IN  LIST
 
 50 i = i - 1
 IF (i == 0) GO TO 100
 
!     FLIP LIST ORDER
 
 m  = 0
 ll = i
 
!     SAVE CURRENT STUFF AT END OF LIST
 
 54 DO  j = 1,5
   l = i + j
   itemp(j) = ilist(l)
 END DO
 
!     PUSH DOWN LIST
 
 DO  j = 1,ll
   k = i - j + 1
   ilist(k+5) = ilist(k)
 END DO
 
!     RESTORE CURRENT STUFF IN FRONT OF LIST, IN FLIPPED ORDER
 
 DO  j = 1,5
   k = m + j
   ilist(k) = itemp(j)
 END DO
 
!     AGAIN
 
 IF (ilist(i+7) /= itemp(2)) GO TO 100
 m = m + 5
 i = i + 5
 GO TO 54
 
!     FOUND IT
 
!     IDR( 3) = ELEMENT TYPE
!     IDR( 4) = SUBCASE NO.
!     IDR( 9) = FORMAT CODE, 1=REAL, 2=REAL/IMAG, 3=MAG/PHASE
!     IDR(10) = NO. OF WORDS PER ENTRY
!     IELEM   = 6 OR 7 FOR ELFORCE (OEFC2) OR STRESS (OESC2)
 
 100 load = idr(4)
 IF   = 0
 IF (idr(9) == 3) IF = 1
 len1 = idr(10)
 LEN  = len1
 ieltp= idr(3)
 IF (ielem < 6 .OR. ielem > 7) GO TO 150
 
!     EXECUTE THIS PORTION FOR STRESSES AND FORCES
 
!     FIND FORMAT TYPE
 
 IF (ifmt(1,ieltp) == -1) GO TO 930
 
!     PICK UP FORMAT POINTER
 
 ifmtp = ifmt(ielem-5,ieltp)
 IF (ifmtp  == 0) GO TO 970
 j = ifmtt(ifmtp)
 
!     SAVE EXTERNAL DATA LENGTH
 
 LEN  = itb(j)
 
!     SAVE MAP OF ITB
 
 DO  i = 1,len1
   k = j + i - 1
   idata(i) = itb(k)
 END DO
 
!     CONVERT POINTERS TO NEW DATA VALUES
 
 IF (ioldld ==    0) GO TO 131
 IF (ioldld /= load) GO TO 150
 131 ioldld = load
 DO  i = 1,llist,5
   IF (ilist(i) /= mid(2,ielem) .OR. ilist(i+1) /= ilist(2)) EXIT
   k = ilist(i+2)
   IF (k <= len1) GO TO 141
   
!     POINTER OUT OF RANGE
   
   CALL mesage (52,ilist(i),ilist(i+1))
   k = len1
   141 k = j + k - 1 + len1
   ilist(i+2) = itb(k)
 END DO
 150 ichk  = 1234321
 lenx  = LEN
 
!     FILE AND LEN WERE SAVED LOCALLY IN FILEX AND LENX, SO THAT THEY
!     CAN BE USED IN RAND2A
 
 RETURN
 
 
 ENTRY rand2a (DATA)
!     ===================
 
!     WILL OBTAIN DATA AND REFORMAT IF NECESSARY
 
!     READ DATA
 
 IF (ichk /= 1234321) CALL mesage (-37,0,NAME)
 CALL READ (*910,*920,filex,DATA(1),len1,0,iflag)
 IF (ielem < 6) RETURN
 
!     APPLY DATA MAP  I.E. REARRANGE DATA ACCORDING TO DATA MAP
 
 DO  i = 1,lenx
   data1(i) = 0
 END DO
 data1(1) = DATA(1)
 DO  i = 2,len1
   j = idata(i)
   data1(j) = DATA(i)
 END DO
!WKBR 9/93 DO 190 I = 1,LEN
 DO  i = 1,lenx
   DATA(i)  = data1(i)
 END DO
 RETURN
 
!     FILE ERRORS
 
 910 ip1 = -2
 911 CALL mesage (ip1,filex,NAME)
 920 ip1 = -3
 GO TO 911
 930 WRITE  (nout,940) uwm,ieltp
 940 FORMAT (a25,' 2185, CURRENTLY RAND2 ROUTINE DOES NOT PROCESS ',  &
     'ELEMENT TYPE',i5)
 GO TO 5
! 950 LOAD = 0
 950 CALL REWIND (filex)
!WKBI 9/93
 WRITE(nout,9901)
 9901  FORMAT(' THE FOLLOWING I/O ERROR OCCURRED MOST LIKELY BECAUSE'  &
     ,/,' THERE WAS A PLOT REQUEST FOR A POINT THAT DOES NOT EXIST.')
 CALL mesage (-2,FILE,NAME)
 GO TO 150
 970 ip1 = -7
 GO TO 911
END SUBROUTINE rand2
