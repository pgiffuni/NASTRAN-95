SUBROUTINE ruler (rule,icp,zrct,onct,list,n,buff,iopt)
     
!     DETERMINES STRING OF ZEROS AND ONES IN LIST BY APPLYING RULE TO
!     CP.
 
 
 INTEGER, INTENT(IN)                      :: rule
 INTEGER, INTENT(IN)                      :: icp(1)
 INTEGER, INTENT(OUT)                     :: zrct
 INTEGER, INTENT(OUT)                     :: onct
 INTEGER, INTENT(OUT)                     :: list(1)
 INTEGER, INTENT(IN)                      :: n
 REAL, INTENT(IN OUT)                     :: buff(1)
 INTEGER, INTENT(IN OUT)                  :: iopt
 EXTERNAL        orf
 INTEGER :: oct, eol,orf,zct
 
 COMMON /two   / two1(32)
 COMMON /zntpkx/ a1(4),l,eol
 
!     PICK UP PARAMETERS
 
 eol = 0
 r   = rule
 namcp = icp(1)
 zct = 0
 oct = 0
 ASSIGN 150 TO is
 IF (r >= 0.0) ASSIGN 140 TO is
 r   = ABS(r)
 l   = 0
 j1  = 0
 m   = 0
 n1  = n
 IF (namcp == 0) GO TO 50
 CALL gopen (namcp,buff,0)
 CALL intpk (*50,namcp,0,1,0)
 GO TO 60
 50 m   = n1
 eol = 1
 60 DO  i = 1,n1
   j = (i+31)/32
   IF (m   >= i) GO TO 90
   IF (eol == 0) GO TO 80
   l = n1
   a1(1) = 0.0
   GO TO 90
   80 CALL zntpki
   90 IF (l == i) GO TO 110
   m = l
   a = 0.0
   GO TO 120
   110 a = a1(1)
   120 IF (iopt == 1 .OR. j <= j1) GO TO 130
   j1 = j
   list(j) = 0
   130 GO TO is, (140,150)
   140 IF (a-r == 0.0) THEN
     GO TO   190
   ELSE
     GO TO   160
   END IF
   150 IF (a-r < 0.0) THEN
     GO TO   160
   ELSE IF (a-r == 0.0) THEN
     GO TO   190
   ELSE
     GO TO   200
   END IF
   160 oct = oct + 1
   IF (iopt == 1) GO TO 180
   k = i - ((i-1)/32)*32
   list(j) = orf(list(j),two1(k))
   CYCLE
   180 list(i) = oct
   CYCLE
   190 zct = zct + 1
   IF (iopt /= 0) list(i) = -zct
 END DO
 zrct = zct
 onct = oct
 IF (namcp /= 0) CALL CLOSE (namcp,1)
 RETURN
END SUBROUTINE ruler
