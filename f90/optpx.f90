SUBROUTINE optpx (dtyp)
     
!     PROCESS PLIMIT CARDS INTO ELEMENT SECTIONS THAT MAY BE READ BY
!     OPTP1D
!     MPT ASSUMED PREPOSITIONED TO PLIMIT CARDS.
 
 
 INTEGER, INTENT(IN)                      :: dtyp(1)
 INTEGER :: count,ycor,b1p1,npow,ept,NAME(2),sysbuf,outtap,  &
      etp(21),any,all,stor(21),BLANK,eject, scrth1,ENTRY,x(7),iy(1)
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /BLANK / skp1(2),count,skp2(2),ycor,b1p1,npow,skp3(4),nklw,  &
     mpt,ept,skp5(4),scrth1,neltyp,ENTRY(21)
 COMMON /optpw1/ kcor,k(10)
 COMMON /zzzzzz/ core(1)
 COMMON /names / nrd,noeor,nwrt,nweor
 COMMON /system/ sysbuf,outtap
 COMMON /gpta1 / ntypes,last,incr,NE(1)
 EQUIVALENCE     (stor(1),k(10)),(core(1),x(1)),(x(7),iy(1))
 DATA    etp   / 21*0 /,  all / 4HALL /,  BLANK / 1H  /,  &
     NAME  / 4H opt,  4HPX        /
 
 maxw  = 0
 iall  = 0
 any   = 0
 nocor = 0
 nogo  = 0
 nx    = 1
 ASSIGN 10 TO iret
 
!     MAKE PRELIMINARY PASS
 
 10 imhere = 10
 CALL READ (*310,*110,mpt,k,9,0,nwds)
 IF (k(1) == all) GO TO 30
 DO  i = 1,ntypes
   IF (dtyp(i) == 0) CYCLE
   idx = incr*(i-1) + 1
   IF (NE(idx  ) /= k(1)) CYCLE
   IF (NE(idx+1) == k(2)) GO TO 40
 END DO
 GO TO 50
 
!     ALL SPECIFIED
 
 30 iall = iall + 1
 GO TO 10
 
!     LEGAL ELEMENT TYPE
 
 40 i = dtyp(i)
 etp(i) = etp(i) + 1
 any = any + 1
 GO TO 10
 
!     ILLEGAL ELEMENT TYPE
 
 50 nogo = nogo + 1
 IF (nogo > 1) GO TO  70
 CALL page2 (-4)
 WRITE  (outtap,60) ufm
 60 FORMAT (a23,' 2290, THE FOLLOWING ILLEGAL ELEMENT TYPES FOUND ON',  &
     ' PLIMIT CARD')
 70 stor(nx  ) = k(1)
 stor(nx+1) = k(2)
 nx = nx + 2
 IF (nx < 20) GO TO 10
 80 i = eject(2)
 IF (i == 0) GO TO 90
 CALL page2 (-2)
 WRITE  (outtap,60) ufm
 90 WRITE  (outtap,100) stor
 100 FORMAT (1H0,9X,10(2A4,1X))
 nx = 1
 GO TO iret, (10,130)
 
!     LAST PLIMIT
 
 110 IF (nx <= 1) GO TO 130
 ASSIGN 130 TO iret
 DO  i = nx,20
   stor(i) = BLANK
 END DO
 GO TO 80
 
!     CONTINUE PROCESSING LEGAL CARDS UNLESS ANY = 0
 
 130 IF (any == 0 .AND. iall == 0) GO TO 300
 CALL bckrec (mpt)
 imhere = 130
 CALL READ (*310,*320,mpt,stor(1),3,noeor,nwds)
 
 loc1 = 1
 
!     START OF OUTPUT LOOP
 
 DO  n = 1,ntypes
   ide = dtyp(n)
   IF (ide <= 0) CYCLE
   idx = ENTRY(ide)
   idx = incr*(idx-1)
   nen = 0
   nde = etp(ide)
   IF (nde <= 0) GO TO 160
   nwds = 0
   
   imhere = 140
   DO  m = 1,nde
     140 CALL READ (*310,*320,mpt,stor(1),9,noeor,nwds)
     IF (stor(1) /= NE(idx+1)) GO TO 140
     IF (stor(2) /= NE(idx+2)) GO TO 140
     CALL optpx1 (*260,stor,nogo,nen,loc1)
   END DO
   CALL bckrec (mpt)
   imhere = 150
   CALL READ (*310,*320,mpt,stor(1),3,noeor,nwds)
   
!     CHECK IF ALL SPECIFIED
   
   160 IF (iall <= 0) GO TO 190
   imhere = 170
   DO  m = 1,iall
     170 CALL READ (*310,*320,mpt,stor(1),9,noeor,nwds)
     IF (stor(1) /= all) GO TO 170
     CALL optpx1 (*260,stor,nogo,nen,loc1)
   END DO
   CALL bckrec (mpt)
   imhere = 180
   CALL READ (*310,*320,mpt,stor(1),3,noeor,nwds)
   
!     CONTINUE PROCESSING LEGAL CARDS - SORT ON SECOND WORD
   
   190 IF (nen == 0) CYCLE
   CALL sort (0,0,4,2,iy(loc1),nen)
   
!     CHECK SECOND WORD
   
   i1  = iy(loc1  )
   i2  = iy(loc1+1)
   i3  = iy(loc1+2)
   i4  = iy(loc1+3)
   loc2= loc1 + nen
   l   = loc2
   IF (l+4 > ycor) nwds = 1
   nx = nen - 3
   IF (nx < 5) GO TO 250
   DO  m = 5,nx,4
     j  = loc1 + m - 1
     j1 = iy(j  )
     j2 = iy(j+1)
     
     IF (i1 >= j1) GO TO 220
     IF (i2 >= j1) GO TO 220
     
!     CHECK FOR EXPANDING THE THRU
     
     IF (i2 /=    j1-1) GO TO 200
     IF (i3 /= iy(j+2)) GO TO 200
     IF (i4 /= iy(j+3)) GO TO 200
     i2 = j2
     IF (m /= nx) CYCLE
     iy(nx) = i1
     GO TO 250
     
!     OUTPUT PLIMIT DATA IN SETS OF 4
     
     200 IF (nogo > 0 .OR. nwds > 0) GO TO 210
     iy(l  ) = i1
     iy(l+1) = i2
     iy(l+2) = i3
     iy(l+3) = i4
     210 l = l + 4
     IF (l+3 > ycor) nwds = nwds + 4
     i1 = j1
     i2 = j2
     i3 = iy(j+2)
     i4 = iy(j+3)
     CYCLE
     
!     OVERLAPPING RANGE ERROR CONDITION
     
     220 CALL page2 (-2)
     WRITE  (outtap,230) ufm,i1,i2,j1,j2
     230 FORMAT (a23,' 2291, PLIMIT RANGE INCORRECT FOR',i8,' THRU',i8,  &
         ' AND',i8,' THRU',i8,'.')
     i1 = j1
     i2 = j2
     nogo = nogo + 1
   END DO
   
!     AFTER ELEMENTS THAT MAY BE OPTIMIZED, FLUSH BUFFER.
   
   250 IF (l+3 > ycor) GO TO 260
   iy(l  ) = iy(nx  )
   iy(l+1) = iy(nx+1)
   iy(l+2) = iy(nx+2)
   iy(l+3) = iy(nx+3)
   l = l + 3
   GO TO 280
   
!     INSUFFICIENT CORE FOR ELEMENTS OF THIS TYPE
   
   260 CALL page2 (-2)
   nocor = 1
   nwds  = nwds + 3
   WRITE  (outtap,270) ufm,NE(idx+1),NE(idx+2),nwds
   270 FORMAT (a23,' 2292, INSUFFICIENT CORE FOR PLIMIT DATA, ELEMENT ',  &
       2A4,i5,' WORDS SKIPPED.')
   nogo = nogo + 1
   
!     WRITE ONTO SCRATCH FILE
   
   280 IF (nogo > 0) CYCLE
   maxw = MAX0(l,maxw)
   stor(1) = ide
   stor(2) = (l-loc2+1)/4
   CALL WRITE (scrth1,stor(1),2,noeor)
   
!     AFTER ELEMENT TYPE, NUMBER WORDS - WRITE DATA
   
   CALL WRITE (scrth1,iy(loc2),l-loc2+1,nweor)
   
 END DO
 
!     END OF OUTPUT LOOP
 
 CALL eof (scrth1)
 300 IF (nogo  == 0) nklw  = maxw
 IF (nogo  > 0) count = -1
 IF (nocor /= 0) nklw  = -64
 RETURN
 
!     ILLEGAL EOF (310), EOR (320)
 
 310 j = -2
 nwds = -222
 GO TO 330
 320 j = -3
 330 WRITE  (outtap,340) imhere,nwds
 340 FORMAT ('  ERROR IN OPTPX.  IMHERE=',i4,',  NWDS=',i6)
 CALL mesage (j,mpt,NAME)
 GO TO 300
END SUBROUTINE optpx
