SUBROUTINE setinp
     
 IMPLICIT INTEGER (a-z)
 EXTERNAL         rshift,complf
 REAL :: fwrd
 DOUBLE PRECISION :: dwrd
 DIMENSION        NAME(2),el(1),gp(1),card(65),typ(100),awrd(2),  &
     pcard(20),pocard(200)
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg /  ufm,uwm,uim
 COMMON /BLANK /  skp11,nsets,skp12(8),pcdb,skp2(9),  &
     merr,plot,msetid,skp3(7),mset,ipcdb
 COMMON /system/  bufsiz,nout,nogo,nin,nsk(81),intr
 COMMON /gpta1 /  ntypes,last,incr,NE(1)
 COMMON /zzzzzz/  x(1)
 EQUIVALENCE      (x(1),el(1),gp(1))
 EQUIVALENCE      (word,awrd(1),iwrd,fwrd,dwrd)
 DATA    inprew,  outrew,rew,norew,eor / 0, 1, 1, 3, 1000000      /
 DATA    blnk  ,  STOP,GO,NAME /4H    ,4HSTOP,4HGO  ,4H set,3HINP /
 DATA    set   ,  incl,  excl,  elem,  grid,  poin,  exce,  TO    /  &
     3HSET ,  4HINCL,4HEXCL,4HELEM,4HGRID,4HPOIN,4HEXCE,2HTO  /
 DATA    thru  ,  all,    ilxx /4HTHRU,3HALL, 2HXX /
 
 CALL delset
 b1 = korsz(x) - 5*bufsiz + 1
 b2 = b1 + bufsiz
 b3 = b2 + bufsiz
 b4 = b3 + bufsiz
 nogo  = 0
 org   = 0
 porg  = -1
 allon = complf(0)
 pocard(200) = rshift(allon,1)
 encard = pocard(200)
 
!     OPEN ALL NECESSARY FILES
 
 iorew = inprew
 IF (intr <= 0) GO TO 10
 pcdb  = ipcdb
 iorew = outrew
 10 CALL OPEN (*210,pcdb,x(b1),iorew)
 IF (intr <= 0) GO TO 50
 
 WRITE (nout,270)
 20 DO  j = 1,20
   pcard(j) = blnk
 END DO
 DO  j = 1,199
   pocard(j) = blnk
 END DO
 CALL xread (*28,pcard)
 IF (pcard(1) == STOP) GO TO 220
 IF (pcard(1) ==   GO) GO TO 40
 CALL xrcard (pocard,199,pcard)
 CALL ifp1pc (1,icont,pocard,org,porg)
 IF (nogo == 0) GO TO 30
 nogo = 0
 28 WRITE (nout,300)
 GO TO 20
 30 WRITE (1,290) pcard
 ie = 1
 DO  j = 1,199
   IF (pocard(j) /= 0) GO TO 32
   DO  jc = 1,5
     IF (pocard(j+jc) /= blnk) GO TO 32
   END DO
   nw = j
   GO TO 34
   32 IF (pocard(j) /= encard) CYCLE
   nw = j
   GO TO 34
 END DO
 nw = 80
 34 CALL WRITE (pcdb,pocard,nw,ie)
 GO TO 20
 40 CALL CLOSE (pcdb,rew)
 IF (intr > 10) nout = 1
 CALL OPEN (*210,pcdb,x(b1),inprew)
 50 IF (intr <= 0) CALL fread (pcdb,0,-2,1)
 CALL gopen (plot,x(b2),outrew)
 CALL gopen (mset,x(b3),outrew)
 CALL gopen (msetid,x(b4),outrew)
 CALL rdmodx (pcdb,mode,word)
 
!     READ MODE FLAG.  SHOULD BE ALPHABETIC
 
 100 CALL READ (*200,*200,pcdb,mode,1,0,i)
 IF (mode < 0) THEN
   GO TO   101
 ELSE IF (mode == 0) THEN
   GO TO   100
 ELSE
   GO TO   102
 END IF
 101 i = 1
 IF (mode == -4) i = 2
 CALL fread (pcdb,0,-i,0)
 GO TO 100
 102 IF (mode < eor) GO TO 103
 CALL fread (pcdb,0,0,1)
 GO TO 100
 103 mode = mode + 1
 CALL rdword (mode,word)
 CALL rdword (mode,word)
 IF (word == set .AND. mode == 0) GO TO 115
 
!     THIS CARD IS A PLOT CONTROL CARD
 
 105 CALL bckrec (pcdb)
 106 CALL READ (*200,*110,pcdb,card,65,1,i)
 WRITE  (nout,108)
 108 FORMAT ('  ARRAY CARD OF 65 TOO SAMLL')
 CALL mesage (-37,0,NAME)
 110 CALL WRITE (plot,card,i,1)
 IF (card(i) == 0.0) THEN
   GO TO   106
 ELSE
   GO TO   100
 END IF
 
!     THIS CARD DEFINES A NEW SET
 
 115 ASSIGN 116 TO tra
 CALL rdmode (*250,*105,*100,mode,word)
 116 setid = iwrd
 nelx  = 0
 ngpx  = b1
 nt    = 0
 xx    = 1
 elgp  = 0
 
 IF (mode <= 0) CALL rdmode (*136,*121,*175,mode,word)
 121 CALL rdword (mode,word)
 
!     CHECK FOR AN -INCLUDE- OR -EXCLUDE- CARD
 
 125 IF (word /= incl .AND. word /= excl .AND. word /= exce) GO TO 128
 126 IF (word == incl) xx = 1
 IF (word == excl) xx =-1
 IF (word == exce) xx =-xx
 IF (mode == 0) CALL rdmode (*136,*127,*175,mode,word)
 127 CALL rdword (mode,word)
 128 IF (word == grid) GO TO 131
 IF (word /= elem) GO TO 147
 
!     ELEMENTS ARE TO BE INCLUDED OR EXCLUDED (BY ID OR TYPE)
 
 elgp = 0
 IF (mode < 0) THEN
   GO TO   136
 ELSE IF (mode == 0) THEN
   GO TO   135
 ELSE
   GO TO   121
 END IF
 
!     A LIST OF GRID POINTS IS TO BE INCLUDED OR EXCLUDED (PERTAIN ONLY
!     TO DEFORMED PLOTS)
 
 131 IF (mode <= 0) CALL rdmode (*131,*132,*175,mode,word)
 132 CALL rdword (mode,word)
 IF (word /= poin .OR. mode /= 0) GO TO 125
 elgp = 1
 
!     A LIST OF ELEMENT OR GRID POINT ID-S CAN BE EXPLICITLY LISTED, OR
!     PREFERABLY A RANGE CAN BE SPECIFIED (SEPARATED BY THE WORD -TO-
!     OR -THRU-)
 
 135 CALL rdmode (*136,*121,*175,mode,word)
 136 ASSIGN 137 TO tra
 GO TO 250
 137 IF (nelx+1 >= ngpx) CALL mesage (-8,0,NAME)
 IF (elgp /= 0) GO TO 138
 nelx = nelx + 1
 el(nelx) = ISIGN(iwrd,xx)
 GO TO 140
 138 ngpx = ngpx - 1
 gp(ngpx) = ISIGN(iwrd,xx)
 
 140 CALL rdmode (*250,*141,*175,mode,word)
 141 CALL rdword (mode,word)
 IF (word /= TO .AND. word /= thru) GO TO 125
 IF (mode /= 0) GO TO 125
 ASSIGN 142 TO tra
 CALL rdmode (*250,*125,*175,mode,word)
 142 IF (nelx+2 >= ngpx) CALL mesage (-8,0,NAME)
 IF (elgp /= 0) GO TO 143
 el(nelx+1) = TO
 el(nelx+2) = iwrd
 nelx = nelx + 2
 GO TO 135
 143 gp(ngpx-1) = TO
 gp(ngpx-2) = ISIGN(iwrd,xx)
 ngpx = ngpx - 2
 GO TO 135
 
!     AN ELEMENT TYPE CAN BE INCLUDED OR EXCLUDED
 
 145 IF (mode <= 0) CALL rdmode (*136,*146,*175,mode,word)
 146 CALL rdword (mode,word)
 147 IF (word == incl .OR. word == excl .OR. word == exce) GO TO 126
 IF (word == grid .OR. word == elem) GO TO 128
 IF (word /= all) GO TO 150
 i  = ntypes + 1
 149 nt = nt + 2
 
!     SECOND WORD FOR EACH TYP LOCATES ELEMENT INCLUDE/EXCLUDE SEARCH
!     POINTER.  ELEMENT ID-S GIVEN PRIOR TO NELX ARE SKIPPED
 
 typ(nt-1) = ISIGN(i,xx)
 typ(nt  ) = nelx + 1
 elgp = 0
 GO TO 145
 
 150 DO  i = 1,ntypes
   idx = (i-1)*incr
   
!     SKIP ELEMENTS WITH
!       1 GRID
!       SCALAR CONNECTIONS POSSIBLE
!       SPECIAL PLOTTER MNEMONIC OF -XX-
   
   IF (NE(idx+10) <= 1 .OR. NE(idx+11) /= 0) CYCLE
   IF (NE(idx+16) == ilxx) CYCLE
   IF (awrd(1) == NE(idx+1) .AND. awrd(2) == NE(idx+2)) GO TO 149
 END DO
 WRITE  (nout,155) ufm,awrd
 155 FORMAT (a23,' 699,',2A4,' ELEMENT IS INVALID')
 nogo = 1
 elgp = 0
 GO TO 145
 
!     A SET HAS BEEN COMPLETELY DEFINED.  FIRST, WRITE THE SET ID
 
 175 IF (nelx == 0 .AND. nt == 0) GO TO 100
 CALL WRITE (msetid,setid,1,0)
 CALL WRITE (mset,setid,1,0)
 
!     WRITE THE SET OF EXPICIT ELEMENT ID-S
 
 CALL WRITE (mset,nelx,1,0)
 CALL WRITE (mset,el,nelx,0)
 
!     DELETE ALL ELEMENT TYPE DUPLICATES + WRITE REMAINING ONES
 
 n = 0
 IF (nt == 0) GO TO 178
 DO  j = 1,nt,2
   xx = typ(j)
   IF (xx == 0) CYCLE
   DO  i = j,nt,2
     IF (i == j .OR. IABS(xx) /= IABS(typ(i))) CYCLE
     
!     DELETE BOTH IF NEGATIVE OF OTHER
     
     IF (xx == -typ(i)) typ(j) = 0
     typ(i) = 0
   END DO
   IF (typ(j) == 0) CYCLE
   n = n + 2
   typ(n-1) = xx
   typ(n  ) = typ(j+1)
 END DO
 178 CALL WRITE (mset,n,1,0)
 CALL WRITE (mset,typ,n,0)
 
!     WRITE THE SET OF EXPLICIT GRID POINT ID-S
 
 n = b1 - ngpx
 CALL WRITE (mset,n,1,0)
 CALL WRITE (mset,gp(ngpx),n,1)
 nsets = nsets + 1
 GO TO 100
 
!     END OF -PCDB-
 
 200 CALL clstab (mset,rew)
 CALL clstab (plot,rew)
 CALL clstab (msetid,norew)
 CALL CLOSE  (pcdb,rew)
 IF (nsets == 0) WRITE (nout,205) uim
 205 FORMAT (a29,', NO SETS EXIST IN PLOT PACKAGE')
 IF (nogo /= 0) CALL mesage (-61,0,0)
 210 RETURN
 220 nogo = 1
 RETURN
 
!     READ AN INTEGER
 
 250 IF (mode == -1) GO TO 260
 IF (mode == -4) iwrd = dwrd
 IF (mode /= -4) iwrd = fwrd
 260 GO TO tra, (116,137,142)
 
 270 FORMAT (' ENTER PLOT DEFINITION OR ''GO'' IF DONE.')
 290 FORMAT (20A4)
 300 FORMAT (' BAD CARD TRY AGIAN')
END SUBROUTINE setinp
