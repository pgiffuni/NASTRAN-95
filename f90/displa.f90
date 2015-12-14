SUBROUTINE displa (gplst,x,s,u,pen,deform,label,pt,b1)
     
 
 INTEGER, INTENT(IN OUT)                  :: gplst(1)
 REAL, INTENT(IN)                         :: x(3,1)
 REAL, INTENT(IN)                         :: s(2,1)
 REAL, INTENT(IN OUT)                     :: u(3,1)
 INTEGER, INTENT(OUT)                     :: pen
 INTEGER, INTENT(IN OUT)                  :: deform
 INTEGER, INTENT(OUT)                     :: label(50)
 REAL, INTENT(OUT)                        :: pt(8)
 INTEGER, INTENT(IN OUT)                  :: b1
 INTEGER :: scr1,ect2,axis,daxis,elid, gpts(12),gp,color,offset
 REAL :: maxdef
 DIMENSION       SIGN(3),a(4), xx(4),yy(4), mvect(3), msg(13)
 COMMON /BLANK / skip(5),ngpset,sk(6),ect2,skp(7),merr,ski(6),scr1
 COMMON /xxparm/ skip1(39),maxdef,defmax,axis(3),daxis(3),  &
     skip2(110),ncntr,cntr(50),icntvl,skppar(6), sk18(18),color
 COMMON /pltdat/ skip3(2),xmin
 DATA    nmsg  / 13   /,  &
     msg   / 4H(33X, 4H,41H, 4H*** , 4HINCO, 4HMPLE, 4HTE p,  &
     4HLOT , 4HDUE , 4HTO i, 4HNPUT, 4H OR , 4HFILE, 4H.)  /
 DATA    mvect / 3*0   /,        kbar,kt3,kq4 / 2HBR,2HT3,2HQ4 /
 
 
 CALL gopen (scr1,gplst(b1),1)
 IF (ABS(defmax) > 1.e-8) GO TO 5
 CALL wrtprt (merr,mvect,msg,nmsg)
 GO TO 120
 
 5 DO  i = 1,3
   SIGN(i) = daxis(i)/axis(i)
 END DO
 DO  gp = 1,ngpset
   DO  i  = 1,3
     j  = axis(i)
     ij = IABS(j)
     a(ij) = SIGN(ij)*u(i,gp)
   END DO
   dmax = maxdef
   IF (dmax < .00001) dmax = 1.0
   DO  i = 1,3
     u(i,gp) = a(i)*(defmax/dmax)
   END DO
 END DO
 INDEX = icntvl - 9
 ncntr = MIN0(ncntr,50)
 IF (cntr(1) /= cntr(2)) GO TO 40
 IF (INDEX <= 3) conmin = u(INDEX,1)
 IF (INDEX > 3) conmin = SQRT(u(1,1)**2 + u(2,1)**2 + u(3,1)**2)
 conmax = conmin
 DO  gp = 1,ngpset
   IF (INDEX > 3) GO TO 25
   conmin = AMIN1(conmin,u(INDEX,gp))
   conmax = AMAX1(conmax,u(INDEX,gp))
   CYCLE
   25 d = SQRT(u(1,gp)**2 + u(2,gp)**2 + u(3,gp)**2)
   conmin = AMIN1(conmin,d)
   conmax = AMAX1(conmax,d)
 END DO
 delta = (conmax-conmin)/FLOAT(ncntr-1)
 cntr(1) = conmin
 j = ncntr - 1
 DO  i = 2,j
   cntr(i) = cntr(i-1) + delta
 END DO
 cntr(ncntr) = conmax
 40 CALL line (0.,0.,0.,0.,pen,+1)
 DO  i = 1,ncntr
   label(i) = 3
 END DO
 50 CALL READ  (*100,*100,ect2,itype,1,0,m)
 offset = 0
 IF (itype == kbar) offset = 6
 IF (itype == kt3 .OR. itype == kq4) offset = 1
 CALL fread (ect2,ngppe,1,0)
 55 CALL fread (ect2,elid,1,0)
 IF (elid == 0) GO TO 50
 CALL fread (ect2,0,-1,0)
 CALL fread (ect2,gpts,ngppe,0)
 IF (offset /= 0) CALL fread (ect2,0,-offset,0)
 IF (ngppe <= 2) GO TO 55
 ij = 1
 ik = 3
 60 j  = 0
 DO  i = ij,ik
   j  = j + 1
   ig = gpts(i)
   ig = IABS(gplst(ig))
   IF (INDEX <= 3) a(j) = u(INDEX,ig)
   IF (INDEX > 3) a(j) = SQRT(u(1,ig)**2 +u(2,ig)**2 +u(3,ig)**2)
   IF (deform /= 0) GO TO 65
   pt(2*j-1) = x(2,ig)
   pt(2*j  ) = x(3,ig)
   CYCLE
   65 pt(2*j-1) = s(1,ig)
   pt(2*j  ) = s(2,ig)
 END DO
 pt(7) = pt(1)
 pt(8) = pt(2)
 a(4)  = a(1)
 DO  i = 1,ncntr
   IF (color == 0) GO TO 75
   j = IABS(color)
   IF (ncntr <= j) pen = i*j/ncntr
   IF (ncntr > j) pen = 1 + i/(ncntr/j)
   IF (pen   > j) pen = j
   75 CONTINUE
   DO  j = 1,3
     xx(j) = xmin - 1.0
     d = a(j) - a(j+1)
     IF (ABS(a(j  )-cntr(i)) > ABS(d) .OR.  &
         ABS(a(j+1)-cntr(i)) > ABS(d)) CYCLE
     IF (d == 0.0) d = 1.0
     xx(j) = pt(2*j-1) + (pt(2*j+1)-pt(2*j-1))*(a(j)-cntr(i))/d
     yy(j) = pt(2*j  ) + (pt(2*j+2)-pt(2*j  ))*(a(j)-cntr(i))/d
   END DO
   xx(4) = xx(1)
   yy(4) = yy(1)
   DO  j = 1,3
     IF (xx(j) < xmin .OR. xx(j+1) < xmin) CYCLE
     CALL line (xx(j),yy(j),xx(j+1),yy(j+1),pen,0)
     label(i) = label(i) + 1
     IF (label(i) /= 4) CYCLE
     label(i) = 0
     CALL WRITE (scr1,i,1,0)
     CALL WRITE (scr1,xx(j),1,0)
     CALL WRITE (scr1,yy(j),1,0)
   END DO
 END DO
 IF (ngppe == 3 .OR. ij == 3) GO TO 55
 gpts(ngppe+1) = gpts(1)
 ij = 3
 ik = 5
 GO TO 60
 100 CALL bckrec (ect2)
 DO  gp = 1,ngpset
   DO  i  = 1,3
     j  = axis(i)
     ij = IABS(j)
     a(i) = SIGN(i)*u(ij,gp)
   END DO
   DO  i = 1,3
     u(i,gp) = a(i)*(dmax/defmax)
   END DO
 END DO
 120 CALL CLOSE (scr1,1)
 RETURN
END SUBROUTINE displa
