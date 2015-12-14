SUBROUTINE xyfind (*,*,*,majid,idz)
     
 
 , INTENT(IN OUT)                         :: *
 , INTENT(IN OUT)                         :: *
 , INTENT(IN OUT)                         :: *
 INTEGER, INTENT(IN OUT)                  :: majid(11)
 INTEGER, INTENT(IN OUT)                  :: idz
 LOGICAL :: random    ,retry
 INTEGER :: FILE      ,vector    ,vecid   ,  &
     z         ,eor       ,flag      ,subc
 COMMON /zzzzzz/ z(1)
 COMMON /xywork/ FILE      ,tcurve(32),ntops     ,PRINT   ,  &
     ifile     ,xaxis(32) ,nbots     ,plot    ,  &
     vector    ,yaxis(32) ,vecid(5)  ,punch   ,  &
     major     ,ytaxis(32),subc(5)   ,center  ,  &
     random    ,ybaxis(32),idin(153) ,buf(100),  &
     ivalue(60),iat       ,idout(300),outopn  ,  &
     steps     ,nat       ,paplot    ,knt
 DATA    eor   / 1 /
 
!     THIS SUBROUTINE LOCATES THE ID RECORD FOR A PARTICULAR ELEMENT OR
!     POINT ID AND IF THIS IS A RANDOM PLOT IT CONSIDERS THE COMPONENT
 
 k = 1
 retry = .false.
 itemp = idz
 IF (subc(FILE) < 0.0) THEN
   GO TO    15
 END IF
 1 CONTINUE
 IF (knt < 0) THEN
   GO TO     3
 ELSE IF (knt == 0) THEN
   GO TO    15
 ELSE
   GO TO     7
 END IF
 3 CONTINUE
 isav = idin(4)
 5 CALL READ (*80,*110,ifile,idin(1),146,1,flag)
 IF (isav == idin(4)) GO TO 21
 CALL fwdrec (*100,ifile)
 GO TO 5
 7 CONTINUE
 isav = idin(4)
 GO TO 11
 9 CALL fwdrec (*100,ifile)
 11 CALL READ (*80,*110,ifile,idin(1),146,1,flag)
 IF (idin(4) == isav) GO TO 9
 GO TO 21
 15 CALL REWIND (ifile)
 CALL fwdrec (*100,ifile)
 vecid(FILE) = 0
 20 CALL READ (*80,*110,ifile,idin(1),146,eor,flag)
 21 CONTINUE
 IF (major /= idin(2)) GO TO 25
 IF (subc(FILE) ==  0) GO TO 30
 IF (subc(FILE) == idin(4)) GO TO 30
 25 CONTINUE
 CALL fwdrec (*100,ifile)
 k = k + 1
 GO TO 20
 
!     MATCH ON MAJOR ID MADE
 
 30 vecid(FILE) = vector
 40 IF (idin(5)/10 == z(idz)) GO TO 90
 itemp = -1
 50 CALL fwdrec (*100,ifile)
 CALL READ (*80,*110,ifile,idin(1),146,eor,flag)
 IF (major == idin(2)) GO TO 40
 
!     ELEMENT DATA ARE NOT IN ASCENDING SORT LIKE GRID DATA, BUT ARE
!     SORTED BY ELEMENT NAME, THEN BY ELEMENT NUMBER.
!     SINCE IT IS POSSIBLE FOR THE DESIRED ELEMENT TO BE AHEAD OF THE
!     CURRENT POSITION OF FILE, REWIND AND TRY AGAIN TO FIND MISSING
!     ELEMENT DATA FOR FORCES AND STRESSES.
 
 80 IF (knt == 0 .OR. retry .OR. subc(FILE) == 0) GO TO 82
 retry = .true.
 GO TO 15
 82 IF (subc(FILE) /= 0) GO TO 85
 subc(FILE) = -1
 RETURN
 
 85 CONTINUE
 vecid(FILE) = 0
 idz = itemp
 CALL REWIND (ifile)
 CALL fwdrec (*100,ifile)
 RETURN 3
 
!     IF RANDOM CHECK COMPONENT FOR MATCH
 
 90 IF (z(idz+1) /= idin(6) .AND. random) GO TO 50
 IF (subc(FILE) == 0) RETURN
 IF (subc(FILE) /= idin(4)) GO TO 50
 RETURN
 
!     EOF HIT WHEN AN EOF SHOULD NOT HAVE BEEN HIT
 
 100 RETURN 1
 
!     EOR HIT WHEN AN EOR SHOULD NOT HAVE BEEN HIT
 
 110 RETURN 2
 
END SUBROUTINE xyfind
