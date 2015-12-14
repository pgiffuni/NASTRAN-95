SUBROUTINE wrtmsg (filex)
     
 
 INTEGER, INTENT(IN)                      :: filex
 EXTERNAL        lshift,rshift,andf,orf,complf
 INTEGER :: FILE, title,ttlsav(32,6),count,lst(50),  &
     for(100),ret,eject,rew,BLANK,formax,mask1(5),  &
     mask2(5),pos,andf,orf,rshift,complf,sysx
!WKBI
 CHARACTER (LEN=1) :: formt(400)
 COMMON /output/ title(32,6)
 COMMON /machin/ mach
 COMMON /system/ sysx(41)
!WKBI
 EQUIVALENCE     (formt, for)
 EQUIVALENCE  (xlst,lst)
 EQUIVALENCE     (sysx( 2),mo   ), (sysx( 9),maxlin),  &
     (sysx(12),count), (sysx(39),nbpc  ), (sysx(40),nbpw ), (sysx(41),ncpw  )
 DATA    lstmax, rew,formax,BLANK/ 50,1,100,4H     /
 
 n2cpw  = ncpw/2
 n2cpw1 = n2cpw - 1
 nbpc2  = 2*nbpc
 mask1(1) = rshift(complf(0),nbpc2)
 mask2(1) = complf(mask1(1))
 DO  i  = 2,n2cpw
   mask1(i) = orf(mask2(1),rshift(mask1(i-1),nbpc2))
   mask2(i) = complf(mask1(i))
 END DO
 FILE = filex
 
 DO  j = 1,6
   DO  i = 1,32
     ttlsav(i,j) = title(i,j)
   END DO
 END DO
 
 30 count = maxlin
 40 CALL READ (*500,*30,FILE,n,1,0,nf)
 IF (n < 0) THEN
   GO TO   100
 ELSE IF (n == 0) THEN
   GO TO   130
 ELSE
   GO TO   110
 END IF
 
!     A TITLE OR SUBTITLE FOLLOWS.
 
 100 n = -n
 IF (n <= 6) CALL fread (FILE,title(1,n),32,0)
 IF (n > 6) CALL fread (FILE,0,-32,0)
 GO TO 30
 
!     A MESSAGE FOLLOWS...N = NUMBER OF LIST ITEMS.
 
 110 IF (n <= lstmax) GO TO 120
 CALL fread (FILE,0,-n,0)
 GO TO 130
 120 IF (n /= 0) CALL fread (FILE,lst,n,0)
 
!     READ THE CORRESPONDING FORMAT...NF = SIZE OF THE FORMAT.
 
 130 CALL fread (FILE,nf,1,0)
 IF (nf < 0) THEN
   GO TO   140
 ELSE IF (nf == 0) THEN
   GO TO   150
 ELSE
   GO TO   160
 END IF
 140 count = count - nf
 GO TO 130
 150 count = maxlin
 GO TO 130
 160 IF (nf <= formax) GO TO 170
 CALL fread (FILE,0,-nf,0)
 GO TO 30
 170 CALL fread (FILE,for,nf,0)
 
!     CONDENSE FOR ARRAY TO ACQUIRE CONTIGUOUS HOLLERITH STRINGS.
 
 IF (ncpw == 4) GO TO 300
 DO  i = 2,nf
   k1 = 1
   pos= 2*i - 1
   j  = (pos+n2cpw1)/n2cpw
   k2 = pos - n2cpw*(j-1)
   ASSIGN 200 TO ret
   GO TO 240
   200 CONTINUE
   k1 = 2
   IF (k2+1 <= n2cpw) GO TO 210
   k2 = 1
   j  = j + 1
   GO TO 220
   210 k2 = k2 + 1
   220 CONTINUE
   ASSIGN 230 TO ret
   GO TO 240
   230 CONTINUE
   CYCLE
   240 IF (k2-k1 < 0) THEN
     GO TO   250
   ELSE IF (k2-k1 == 0) THEN
     GO TO   260
   ELSE
     GO TO   270
   END IF
   250 for(j) = orf(andf(for(j),mask1(k2)),  &
       lshift(andf(for(i),mask2(k1)),(nbpc2*(k1-k2))))
   GO TO 280
   260 for(j) = orf(andf(for(j),mask1(k2)),andf(for(i),mask2(k1)))
   GO TO 280
   270 for(j) = orf(andf(for(j),mask1(k2)),  &
       rshift(andf(for(i),mask2(k1)),(nbpc2*(k2-k1))))
   GO TO 280
   280 CONTINUE
   GO TO ret, (200,230)
 END DO
 300 CONTINUE
 
!     PRINT THE LINE
 
 IF (eject(1) == 0) GO TO 450
 DO  j = 4,6
   DO  i = 1,32
     IF (title(i,j) /= BLANK) GO TO 420
   END DO
   count = count - 1
   CYCLE
   420 WRITE  (mo,430) (title(i,j),i=1,32)
   430 FORMAT (2X,32A4)
 END DO
 WRITE  (mo,430)
 count = count + 1
 
 450 IF(n == 0 .AND. (mach == 5 .OR. mach == 12) )GO TO 470
 IF (mach == 5 .OR. mach == 12 ) GO TO 460
 CALL forwrt ( formt, lst, n )
 GO TO 40
 460 WRITE (mo,for,ERR=465) (lst(j),j=1,n)
 465 CONTINUE
 GO TO 40
 470 WRITE (mo,for)
 GO TO 40
 
!     END OF MESSAGE FILE
 
 500 CALL CLOSE (FILE,rew)
 DO  j = 1,6
   DO  i = 1,32
     title(i,j) = ttlsav(i,j)
   END DO
 END DO
 RETURN
END SUBROUTINE wrtmsg
