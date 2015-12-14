SUBROUTINE gptlbl (gplst,x,u,deform,buf)
     
 
 INTEGER, INTENT(IN)                      :: gplst(1)
 REAL, INTENT(IN)                         :: x(3,1)
 REAL, INTENT(IN)                         :: u(2,1)
 INTEGER, INTENT(IN OUT)                  :: deform
 INTEGER, INTENT(IN OUT)                  :: buf
 INTEGER :: exgpid,rew,gp,gpt,gpx
 
 COMMON /BLANK / ngp,skp1(9),skp2(5),exgpid
 COMMON /pltdat/ skpplt(20),skpa(3),cntx
 DATA    inprew, rew / 0,1 /
 
 CALL gopen (exgpid,gplst(buf),inprew)
 CALL typint (0,0,0,0,0,-1)
 DO  gp = 1,ngp
   CALL fread (exgpid,gpt,1,0)
   CALL fread (exgpid,gpx,1,0)
   gpx = gplst(gpx)
   
!     IF THE GRID POINT INDEX IS 0 (NOT IN SET) OR NEGATIVE (EXCLUDED),
!     NEVER PUT A LABEL AT THAT GRID POINT.
   
   IF (gpx <= 0) CYCLE
   
!     TYPE THE GRID POINT ID
   
   IF (deform /= 0) GO TO 111
   xx = x(2,gpx)
   yy = x(3,gpx)
   GO TO 112
   111 xx = u(1,gpx)
   yy = u(2,gpx)
   112 CALL typint (xx+cntx,yy,1,gpt,0,0)
 END DO
 
 CALL CLOSE (exgpid,rew)
 CALL typint (0,0,0,0,0,1)
 RETURN
END SUBROUTINE gptlbl
