SUBROUTINE cdtfbs (dx,dy,iobuf,fileu,nrow)
     
!     CDTFBS IS A SPECIAL VERSION OF GFBS USED BY COMPLEX DETERMINATE
!     METHOD
 
!     DEFINITION OF INPUT PARAMETERS
 
!     FILEU = MATRIX CONTROL BLOCK FOR THE UPPER TRIANGLE U
!     DX    = THE LOAD VECTOR B
!     DY    = THE SOLUTION VECTOR X
!     IOBUF = THE INPUT BUFFER
 
 
 DOUBLE PRECISION, INTENT(IN OUT)         :: dx(1)
 DOUBLE PRECISION, INTENT(IN)             :: dy(1)
 INTEGER, INTENT(IN OUT)                  :: iobuf(1)
 INTEGER, INTENT(IN)                      :: fileu(7)
 INTEGER, INTENT(IN)                      :: nrow
 INTEGER :: typear   ,cdp      ,parm(4) 
 DOUBLE PRECISION :: dtemp
 COMMON /names /  rd        ,rdrew    ,wrt      ,wrtrew   ,  &
     rew       ,norew    ,eofnrw   ,rsp      ,  &
     rdp       ,csp      ,cdp      ,sqr      ,  &
     rect      ,diag     ,lowtri   ,uprtri   , sym       ,row      ,identy
 COMMON /unpakx/  itypex    ,iunpak   ,junpak   ,incrx
 DATA    parm(3), parm(4)   /4HDETF,4HBS  /
 
 incrx  = 1
 itypex = cdp
 typear = cdp
 
!     BEGIN BACKWARD PASS
 
 ioff = fileu(7) - 1
 parm(2) = fileu(1)
 CALL OPEN (*80,fileu,iobuf,rd)
 DO  i = 1,nrow
   iunpak = 0
   j  = nrow - i + 1
   jj = j + j
   CALL bckrec (fileu)
   CALL unpack (*90,fileu,dy)
   CALL bckrec (fileu)
   ising = 0
   k = (junpak-iunpak+1)*2
   ju = junpak + junpak
   GO TO 30
   10 ising = 1
   
!     DIVIDE BY THE DIAGONAL
   
   dtemp    = (dx(jj)*dy(k-1)-dx(jj-1)*dy(k))/(dy(k)**2+dy(k-1)**2)
   dx(jj-1) = (dx(jj-1)*dy(k-1)+dx(jj)*dy(k))/(dy(k)**2+dy(k-1)**2)
   dx(jj  ) = dtemp
   20 k  = k  - 2
   ju = ju - 2
   junpak = junpak - 1
   IF (k == 0) GO TO 60
   IF (dy(k) == 0.d0 .AND. dy(k-1) == 0.d0) GO TO 20
   30 IF (junpak-j < 0) THEN
     GO TO    50
   ELSE IF (junpak-j == 0) THEN
     GO TO    10
   END IF
   40 jk = (j-ioff)*2
   dx(jk-1) = dx(jk-1) - dx(ju-1)*dy(k-1) + dx(ju  )*dy(k)
   dx(jk  ) = dx(jk  ) - dx(ju  )*dy(k-1) - dx(ju-1)*dy(k)
   GO TO 20
   50 CONTINUE
   dx(ju-1) = dx(ju-1) - dx(jj-1)*dy(k-1) + dx(jj  )*dy(k)
   dx(ju  ) = dx(ju  ) - dx(jj  )*dy(k-1) - dx(jj-1)*dy(k)
   GO TO 20
   60 IF (ising == 0) GO TO 90
 END DO
 CALL CLOSE (fileu,rew)
 RETURN
 
 80 parm(1) = -1
 GO TO 100
 90 parm(1) = -5
 100 CALL mesage (parm(1),parm(2),parm(3))
 RETURN
END SUBROUTINE cdtfbs
