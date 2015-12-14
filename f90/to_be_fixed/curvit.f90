SUBROUTINE curvit (indep,ni,dep,nd,ifile,z,iz,lz,mclose,toler,  &
        mcsid,xscale,yscale)
     
!     PERFORMS LOCAL INTERPOLATION
 
!     INDEP  = X,Y COORDINATES OF INDEPENDENT ELEMENT CENTERS (2 X NI)
!     DEP    = X,Y COORDINATES OF DEPENDENT GRID POINTS  (2 X ND)
!     IFILE  = FILE TO WRITE SPECIAL FORM ROWS OF G-MATRIX
!     Z      = REAL AREA OF CORE, LENGTH = LZ.
!     IZ     = EQUIVALENT INTEGER AREA OF CORE, LENGTH = LZ.
!     MCLOSE = NUMBER OF CLOSEST INDEPENDENT POINTS TO USE
!     TOLER  = PERCENT OF DISTANCE FROM A DEPENDENT POINT TO
!              INDEPENDENT POINT NUMBER -NCLOSE- POINTS FURTHER OUT ARE
!              ALLOWED TO BE SUCH AS TO BE INCLUDED IN A LOCAL
!              INTERPOLATION.
 
 
 REAL, INTENT(IN)                         :: indep(2,1)
 INTEGER, INTENT(IN)                      :: ni
 REAL, INTENT(IN)                         :: dep(2,1)
 INTEGER, INTENT(IN)                      :: nd
 INTEGER, INTENT(IN OUT)                  :: ifile
 REAL, INTENT(OUT)                        :: z(1)
 INTEGER, INTENT(OUT)                     :: iz(1)
 INTEGER, INTENT(IN)                      :: lz
 INTEGER, INTENT(IN OUT)                  :: mclose
 REAL, INTENT(IN)                         :: toler
 INTEGER, INTENT(IN OUT)                  :: mcsid
 REAL, INTENT(IN OUT)                     :: xscale
 REAL, INTENT(IN OUT)                     :: yscale
 INTEGER :: sysbuf, subr(2), itemp(2), rd, rdrew, wrt,  &
     wrtrew, clsrew, cls, eor
 
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm, uwm
 COMMON /system/ sysbuf, ioutpt
 COMMON /names / rd, rdrew, wrt, wrtrew, clsrew, cls
 DATA    subr  / 4HCURV ,4HIT  /, eor, noeor / 1, 0 /
 
 nclose = MIN0(mclose,ni)
 IF (nclose <= 2) nclose = ni
 
!     COMPUTE TOLERANCE MULTIPLIER WITH RESPECT TO SQUARES.
!     TOLERANCE IS IN PERCENT OF DISTANCE TO POINT NUMBER -NCLOSE- IN
!     FINAL LIST
 
 tol = (1.0 + toler/100.0)**2
 
!     THUS IF DISTANCE FROM THE DEPENDENT POINT TO INDEPENDENT POINT
!     NUMBER -NCLOSE- = LSQ, ADDITIONAL INDEPENDENT POINTS WILL BE
!     INCLUDED IF THE SQUARE OF THEIR DISTANCE TO THE DEPENDENT POINT
!     IS .LE. TOL TIMES LSQ.
 
 
!     ALLOCATE BUFFER FOR -IFILE- AND OPEN -IFILE-.
 
 ibuf = lz - sysbuf
 jz   = ibuf - 1
 icrq = -jz
 IF (jz <= 0) GO TO 900
 CALL gopen (ifile,iz(ibuf),1)
 
!     EACH ROW OF G-MATRIX WILL BE WRITTEN AS A LOGICAL RECORD
!     WITH PAIRS OF
!                 1- INDEPENDENT POINT INDEX
!                 2- G VALUE
 
 
!     SHORT CUT WILL BE TAKEN IF ALL INDEPENDENT POINTS ARE TO BE USED
!     FOR INTERPOLATION AT EACH DEPENDENT POINT.
 
 IF (nclose == ni) GO TO 550
 
!     MASTER LOOP ON DEPENDENT POINTS. EACH DEPENDENT POINT RESULTS IN
!     A VARIABLE LENGTH ROW OF G-MATRIX DEPENDING ON HOW MANY
!     INDEPENDENT POINTS ARE SELECTED FOR USE. (AT LEAST 3 MUST BE USED)
 
 80 DO  i = 1,nd
   
!     LIST OF DISTANCE SQUARES OF ALL INDEPENDENT POINTS TO
!     CURRENT DEPENDENT POINT IS FORMED.
   
!     SELECTION OF THE -NCLOSE- SMALLEST VALUES IS THEN MADE.
   
!     THEN ANY OTHER INDEPENDENT POINTS WITHIN TOLERANCE RANGE OF
!     POINT NUMBER -NCLOSE- IN LIST ARE ADDED.
   
   fmax = 0.0
   x    = dep(1,i)
   y    = dep(2,i)
   icrq = ni - jz
   IF (ni > jz) GO TO 900
   DO  j = 1,ni
     z(j) = (xscale*(indep(1,j)-x))**2 + (yscale*(indep(2,j)-y))**2
     IF (z(j) <= fmax) CYCLE
     fmax = z(j)
   END DO
   fmax = 2.0*fmax + 1.0
   
!     ALLOCATE FOR LIST OF INDEXES TO THE MINIMUMS.
   
   ilist = ni + 1
   nlist = ni
   
!     FIND -NCLOSE- SMALLEST VALUES.
   
   DO  j = 1,nclose
     fmin = fmax
     
     DO  k = 1,ni
       IF (fmin-z(k) > 0.0) THEN
         GO TO   150
       ELSE
         GO TO   160
       END IF
       150 fmin = z(k)
       idx  = k
     END DO
     
!     ADD INDEX TO THIS MINIMUM TO THE LIST
     
     icrq = nlist + 1 - jz
     IF (icrq > 0) GO TO 900
     iz(nlist+1) = idx
     nlist = nlist + 1
     
!     RESET THIS VALUE SO IT CAN NOT BE USED AGAIN
     
     z(idx) = fmax
   END DO
   
!     ADD ANY ADDITIONAL INDEPENDENT POINTS WITHIN TOLERANCE RANGE OF
!     LAST ONE SELECTED ABOVE.
   
   fmax = tol*fmin
   DO  j = 1,ni
     IF (z(j) > fmax) CYCLE
     icrq = nlist + 1 - jz
     IF (icrq > 0) GO TO 900
     iz(nlist+1) = j
     nlist = nlist + 1
   END DO
   
!     LIST IS COMPLETE THUS MOVE IT TO THE BEGINNING OF THE CORE BLOCK.
   
   j = 0
   DO  k = ilist,nlist
     j = j + 1
     iz(j) = iz(k)
   END DO
   ilist = 1
   nlist = j
   ipts  = j
   
!     HERE AND IZ(ILIST) TO IZ(NLIST) CONTAINS LIST OF
!     POSITION INDEXES OF INDEPENDENT POINT COORDINATES TO BE USED.
   
!     NOW SET UP LIST OF XY-CCORDINATES OF THESE INDEPENDENT POINTS
!     FOR THE SSPLIN CALL.
   
   ixy  = nlist + 1
   nxy  = nlist + 2*ipts
   icrq = nxy - jz
   IF (nxy > jz) GO TO 900
   jxy  = nlist
   DO  j = ilist,nlist
     k    = iz(j)
     z(jxy+1) = indep(1,k)
     z(jxy+2) = indep(2,k)
     jxy  = jxy + 2
   END DO
   
!     NOW READY FOR SSPLIN ROUTINE CALL.
   
   CALL ssplin (ipts,z(ixy),1,dep(1,i),0,0,0,1,0,z(jxy+1),jz-jxy, ising)
   IF (ising /= 2) GO TO 300
   
!     ILL-CONDITION FOR THIS DEPENDENT POINT - NO SOLUTION POSSIBLE.
   
   CALL page2 (4)
   WRITE  (ioutpt,250) uwm,i,mcsid
   250 FORMAT (a25,' 2252. (CURVIT-1) LOCAL INTERPOLATION USING INDE',  &
       'PENDENT VALUES WITHIN RANGE OF THE', /5X,i7,'-TH SORTED ',  &
       'ORDER GRID ID INVOLVED WITH RESPECT TO MATERIAL COORDIN',  &
       'ATE SYSTEM ID',i9, /5X,'CAN NOT BE COMPLETED.  ILL-CONDI',  &
       'TION MAY HAVE RESULTED FROM ALIGNMENT OF INDEPENDENT ',  &
       'VALUE COORDINATES.', /5X,  &
       'OUTPUT FOR THE GRID ID IN QUESTION WILL NOT APPEAR.')
   ipts = 0
   GO TO 340
   
!     REPLACE INDEPENDENT POINT XY PAIRS WITH SPECIAL FORM DEPENDENT
!     POINT G-MATRIX OUTPUT ROW.
   
   300 k1 = ilist
   k2 = jxy + 1
   DO  j = ixy,nxy,2
     iz(j ) = iz(k1)
     z(j+1) = z(k2)
     k1 = k1 + 1
     k2 = k2 + 1
   END DO
   
   340 CALL WRITE (ifile,iz(ixy),2*ipts,eor)
   
!  GO PROCESS NEXT DEPENDENT POINT.
   
 END DO
 GO TO 800
 
!     CHECK FOR SUFFICIENT CORE FOR SHORT CUT.
 
 550 n = ni + 3
 n = n**2 + 3*n + ni*nd + n*nd
 IF (n > jz) GO TO 80
 
!     CALL SSPLIN AND GET G-MATRIX STORED BY ROWS.
 
 CALL ssplin (ni,indep(1,1),nd,dep(1,1),0,0,0,1,0,z(1),jz,ising)
 IF (ising /= 2) GO TO 650
 n = 0
 WRITE (ioutpt,250) uwm,n,mcsid
 
!     OUTPUT NULL ROW FOR EACH DEPENDENT POINT.
 
 DO  i = 1,nd
   CALL WRITE (ifile,0,0,eor)
 END DO
 GO TO 800
 
!     OUTPUT ROWS OF G-MATRIX WITH INDEXES.
 
 650 k = 0
 DO  i = 1,nd
   DO  j = 1,ni
     k = k + 1
     itemp(1) = j
     itemp(2) = iz(k)
     CALL WRITE (ifile,itemp(1),2,noeor)
   END DO
   CALL WRITE (ifile,0,0,eor)
 END DO
 
!     ALL G-MATRIX ROWS COMPLETE. (ROWS SINGULAR ARE EMPTY LOGICAL
!     RECORDS IN -IFILE- )
 
 800 CALL CLOSE (ifile,clsrew)
 RETURN
 
 900 CALL mesage (-8,icrq,subr)
 RETURN
END SUBROUTINE curvit
