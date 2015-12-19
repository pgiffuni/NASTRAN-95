SUBROUTINE sorti (inpfl,outfl,nwds,keywrd,l,nx)
     
!     WITH ENTRY POINT SORTI2 TO SORT TABLE BY 2 KEY WORDS
 
!     THIS SORTING ROUTINE WAS CALLED SORT BEFORE, AND IS NOW RENAMED
!     SORTI. IT IS CAPABLE FOR IN-CORE SORTING AND FILE SORT.
 
!     THE NEW SUBROUTINE SORT IS A TRUNCATED VERSION OF THIS ROUTINE
!     ONLY FOR IN-CORE SORTING. IT CAN HANDLE INTEGER, REAL, BCD(A4),
!     BCD(A8), BCD(A7), AND 2-KEY SORTINGS.
 
!     (95 PERCENT OF NASTRAN ROUTINES ACTUALLY CALL SORT. THE REMAINING
!       5 PERCENT CALL SORTI)
 
!     IF INPFL AND OUTFL ARE ZERO, CALLING ROUTINE SHOULD CALL SORT
!     FOR EFFICIENCY
 
!     THE OLD SHUTTLE EXCHANGE, WHICH WAS VERY SLOW, IS NOW REPLACED BY
!     A SUPER FAST SORTER, A MODIFIED SHELL SORT.
 
!     THIS MODIFIED VERSION ALSO SORTS TABLE OF ANY LENGTH (PREVIOUSLY N
!     OF WORDS PER ENTRY, NWDS, WAS LIMITED TO 20)
 
 
 INTEGER, INTENT(IN)                      :: inpfl
 INTEGER, INTENT(IN OUT)                  :: outfl
 INTEGER, INTENT(IN)                      :: nwds
 INTEGER, INTENT(IN OUT)                  :: keywrd
 INTEGER, INTENT(IN OUT)                  :: l(nwds,2)
 INTEGER, INTENT(IN)                      :: nx
 INTEGER :: scra,scrb,scrc,dist1,dist2,dummy,total,out,  &
     subr(2), temp,FILE,r,bufa,bufb,bufc, sysbuf,bufin,two,two31
 COMMON /setup / nfile(6),bufin
 COMMON /system/ sysbuf,dum38(38),nbpw
 COMMON /two   / two(16)
 EQUIVALENCE     (nfile(1),scrb),(nfile(2),scrc),(nfile(3),scra)
 DATA    subr  / 4HSORT, 4HI    /
 
 key2 = 1
 GO TO 10
 
 
 ENTRY sorti2 (inpfl,outfl,nwds,keywrd,l,nx)
!     ==========================================
 
 key2 = 2
 
!     IF INPFL EQ 0, CORE BLOCK L OF LENGTH NX IS TO BE SORTED
!     IF INPFL NE 0, INPFL IS TO BE SORTED USING BLOCK L
 
 10 keywd = IABS(keywrd)
 nnn   = nx
 IF (nnn < nwds) GO TO 350
 j = 30
 IF (nbpw >= 60) j = 62
 two31 = 2**j
 IF (inpfl == 0) GO TO 30
 bufa  = nx - sysbuf + 1
 
!     MINIMUM CORE REQUIREMENT = 2 X NUMBER OF WORDS PER ENTRY
 
 nz  = bufa - 1
 IF (nz < nwds+nwds) GO TO 360
 CALL OPEN (*370,scra,l(bufa,1),1)
 nn  = (nz/nwds)*nwds
 nnn = nn
 out = scra
 nrec= 0
 20 CALL READ (*430,*170,inpfl,l,nn,0,nnn)
 
!     SORT PHASE --
 
 30 LEN = nnn/nwds
 IF (LEN*nwds /= nnn) GO TO 365
 m = LEN
 IF (keywrd >= 0) GO TO 40
 
!                     - INTEGER SORT ONLY -
!     IF ORIGINAL ORDER IS TO BE MAINTAINED WHERE DUPLICATE KEYWORDS MAY
!     OCCUR, ADD INDICES TO THE KEYWORDS (GOOD FOR BOTH POSITIVE AND
!     NEGATIVE RANGES, AND BE SURE THAT KEYWORDS ARE NOT OVERFLOWED),
!     SORT THE DATA, AND REMOVE THE INDICES LATER
 
!     IF ANY KEYWORD OVERFLOWS, SWITCH TO SHUTTLE EXCHANGE METHOD
!     LIMIT IS THE MAX VALUE BEFORE INTEGER OVERFLOW
 
 IF (LEN >= two(16) .AND. nbpw <= 32) GO TO 130
 limit = (two31-LEN)/LEN
 DO  i = 1,LEN
   j = l(keywd,i)
   IF (IABS(j) > limit) GO TO 124
   j = j*LEN + i
   k = -1
   IF (j < 0) k = -LEN
   l(keywd,i) = j + k
 END DO
 IF (key2 == 1) GO TO 40
 DO  i = 1,LEN
   j = l(keywd+1,i)
   IF (IABS(j) > limit) GO TO 120
   j = j*LEN + i
   k = -1
   IF (j < 0) k = -LEN
   l(keywd+1,i) = j + k
 END DO
 
!     SORT BY
!     MODIFIED SHELL METHOD, A SUPER FAST SORTER
 
 40 m = m/2
 IF (m == 0) GO TO 110
 j = 1
 k = LEN - m
 45 i = j
 50 n = i + m
 IF (l(keywd,i)-l(keywd,n) < 0) THEN
   GO TO   105
 ELSE IF (l(keywd,i)-l(keywd,n) == 0) THEN
   GO TO    60
 ELSE
   GO TO    95
 END IF
 60 IF (key2 == 1) GO TO 105
 IF (l(keywd+1,i)-l(keywd+1,n) > 0) THEN
   GO TO    95
 ELSE
   GO TO   105
 END IF
 95 DO  r = 1,nwds
   temp   = l(r,i)
   l(r,i) = l(r,n)
   l(r,n) = temp
 END DO
 i = i - m
 IF (i >= 1) GO TO 50
 105 j = j + 1
 IF (j-k > 0) THEN
   GO TO    40
 ELSE
   GO TO    45
 END IF
 110 IF (keywrd >= 0) GO TO 160
 DO  i = 1,LEN
   l(keywd,i) = l(keywd,i)/LEN
   IF (key2 == 2) l(keywd+1,i) = l(keywd+1,i)/LEN
 END DO
 GO TO 160
 
!     SORT BY
!     SHUTTLE EXCHANGE METHOD, A SLOW SORTER
!     (THIS WAS NASTRAN ORIGINAL SORTER, MODIFIED FOR 2-D ARRAY
!     OPERATION WITH 20-COLUMN LIMITATION REMOVED)
 
 120 IF (i <= 1) GO TO 123
 j = i - 1
 DO  i = 1,j
   l(keywd+1,i) = l(keywd+1,i)/LEN
 END DO
 123 i = LEN
 124 IF (i <= 1) GO TO 130
 j = i - 1
 DO  i = 1,j
   l(keywd,i) = l(keywd,i)/LEN
 END DO
 
 130 DO  ii = 2,LEN
   jj = ii - 1
   IF (l(keywd,ii)-l(keywd,jj) < 0) THEN
     GO TO   135
   ELSE IF (l(keywd,ii)-l(keywd,jj) == 0) THEN
     GO TO   133
   ELSE
     GO TO   155
   END IF
   133 IF (key2 == 1) CYCLE
   IF (l(keywd+1,ii) >= l(keywd+1,jj)) CYCLE
   135 jj = jj - 1
   IF (jj <= 0) GO TO 140
   IF (l(keywd,ii)-l(keywd,jj) < 0) THEN
     GO TO   135
   ELSE IF (l(keywd,ii)-l(keywd,jj) == 0) THEN
     GO TO   137
   ELSE
     GO TO   140
   END IF
   137 IF (key2 == 2) IF (l(keywd+1,ii)-l(keywd+1,jj)) 135,140,140
   140 jj = jj + 2
   DO  i = 1,nwds
     temp = l(i,ii)
     m = ii
     DO  j = jj,ii
       l(i,m) = l(i,m-1)
       m = m - 1
     END DO
     l(i,jj-1) = temp
   END DO
   155 CONTINUE
 END DO
 
!     IF CORE SORT, SORT IS COMPLETED. IF FILE SORT, WRITE BLOCK ON
!     SCRATCH FILE TO BE MERGED LATER.
 
 160 IF (inpfl == 0) GO TO 350
 165 CALL WRITE (scra,l,nnn,1)
 nrec = nrec + 1
 IF (nnn-nn == 0) THEN
   GO TO    20
 ELSE
   GO TO   180
 END IF
 170 IF (nnn > 0) THEN
   GO TO   175
 ELSE
   GO TO   180
 END IF
 175 IF (nnn-nwds-nwds < 0) THEN
   GO TO   165
 ELSE
   GO TO    30
 END IF
 180 CALL CLOSE (scra,1)
 
!     IF ONLY ONE RECORD, BYPASS MERGE
 
 IF (nrec == 1) GO TO 320
 
!     COMPUTE OPTIMUM DISTRIBUTION OF SORTED RECORDS ON TWO SCRATCH
!     FILES FOR MERGE PHASE USING FIBONACCI SEQUENCE
 
 level = 0
 dist1 = 1
 dist2 = 0
 total = 1
 190 dummy = total - nrec
 IF (dummy >= 0) GO TO 195
 dist1 = dist1 + dist2
 dist2 = dist1 - dist2
 total = dist1 + dist2
 level = level + 1
 GO TO 190
 195 bufb  = bufa - sysbuf
 bufc  = bufb - sysbuf
 IF (bufc < 1) GO TO 360
 nn = bufb - 1
 
!     COPY N SORTED RECORDS ONTO SECOND SCRATCH FILE
 
 CALL OPEN (*370,scra,l(bufa,1),0)
 CALL OPEN (*380,scrb,l(bufb,1),1)
 n = dist2 - dummy
 DO  i = 1,n
   200 CALL READ  (*440,*205,scra,l,nn,0,nflag)
   CALL WRITE (scrb,l,nn,0)
   GO TO 200
   205 CALL WRITE (scrb,l,nflag,1)
 END DO
 CALL CLOSE (scrb,1)
 CALL CLOSE (scra,2)
 nfile(4) = scrb
 nfile(5) = scrc
 k = 4
 
!     MERGE PHASE ---
!     INPUT FILE WITH GREATER NUMBER IF RECORDS = IN1
!     INPUT FILE WITH LESSER  NUMBER OF RECORDS = IN2
!     EACH PASS MERGES ALL RECORDS FROM IN2 WITH LIKE NUMBER OF RECORDS
!     (INCLUDING DUMMY RECORDS) FROM IN1 ONTO OUT. FOR NEXT PASS IN1
!     BECOMES IN2, IN2 BECOMES OUT, AND OUT BECOMES IN1.
 
 DO  i = 1,level
   k = k - 1
   IF (k == 0) k = 3
   in1 = nfile(k)
   in2 = nfile(k+1)
   out = nfile(k+2)
   last= 2
   CALL OPEN (*390,in1,l(bufa,1),2)
   CALL OPEN (*400,in2,l(bufb,1),2)
   CALL OPEN (*410,out,l(bufc,1),1)
   DO  j = 1,dist2
     if1 = nwds
     if2 = nwds
     CALL READ (*450,*275,in1,l,nwds,0,if1)
     IF (dummy > 0.0) THEN
       GO TO   280
     END IF
     210 CALL READ (*460,*290,in2,l(1,2),nwds,0,if2)
     220 IF (l(keywd,1)-l(keywd,2) < 0) THEN
       GO TO   260
     ELSE IF (l(keywd,1)-l(keywd,2) == 0) THEN
       GO TO   230
     ELSE
       GO TO   270
     END IF
     230 IF (key2 == 2) IF (l(keywd+1,1)-l(keywd+1,2)) 260,260,270
     260 CALL WRITE (out,l,nwds,0)
     CALL READ  (*450,*275,in1,l,nwds,0,if1)
     IF (if2 > 0) THEN
       GO TO   220
     ELSE
       GO TO   260
     END IF
     270 CALL WRITE (out,l(1,2),nwds,0)
     CALL READ  (*460,*290,in2,l(1,2),nwds,0,if2)
     IF (if1 > 0) THEN
       GO TO   220
     ELSE
       GO TO   270
     END IF
     275 IF (if2 > 0) THEN
       GO TO   270
     ELSE
       GO TO   300
     END IF
     280 dummy = dummy - 1
     if2 = 0
     290 IF (if1 > 0) THEN
       GO TO   260
     ELSE
       GO TO   300
     END IF
     300 CALL WRITE (out,0,0,1)
   END DO
   dist2 = dist1 - dist2
   dist1 = dist1 - dist2
   IF (dist2 == 0) last = 1
   CALL CLOSE (in1,last)
   CALL CLOSE (in2,1)
   CALL CLOSE (out,1)
 END DO
 
!     COPY PHASE ---
!     IF OUTPUT FILE IS NOT SPECIFIED, NFILE(6) WILL CONTAIN NAME OF
!     SCRATCH FILE CONTAINING OUTPUT
 
 320 nfile(6) = out
 IF (outfl == 0) GO TO 350
 CALL OPEN (*410,out,l(bufa,1),0)
 IF (inpfl /= outfl) GO TO 330
 CALL CLOSE (inpfl,1)
 CALL OPEN  (*420,inpfl,l(bufin,1),1)
 330 CALL READ  (*470,*340,out,l,nz,0,nflag)
 CALL WRITE (outfl,l,nz,0)
 GO TO 330
 340 CALL WRITE (outfl,l,nflag,1)
 CALL CLOSE (out,1)
 350 RETURN
 
!     ERRORS
 
 360 j = -8
 FILE = 0
 GO TO  500
 365 j = -37
 GO TO  500
 370 FILE = scra
 GO TO  480
 380 FILE = scrb
 GO TO  480
 390 FILE = in1
 GO TO  480
 400 FILE = in2
 GO TO  480
 410 FILE = out
 GO TO  480
 420 FILE = inpfl
 GO TO  480
 430 FILE = inpfl
 GO TO  490
 440 FILE = scra
 GO TO  490
 450 FILE = in1
 GO TO  490
 460 FILE = in2
 GO TO  490
 470 FILE = out
 GO TO  490
 480 j = -1
 GO TO  500
 490 j = -2
 500 CALL mesage (j,FILE,subr)
 RETURN
END SUBROUTINE sorti
