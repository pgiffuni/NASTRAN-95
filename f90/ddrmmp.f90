SUBROUTINE ddrmmp(*,z,ncore,lused,ixytyp,icase,buff,anyxy)
!*****
!  BUILD LIST OF POINTS IN SORT FOR WHICH XYCDB OUTPUT REQUESTS EXIST
!  OF FILE TYPE -IXYTYP- AND OF SUBCASE 0 AND SUBCASE -ICASE-.
!*****
 
 INTEGER, INTENT(OUT)                     :: z(1)
 INTEGER, INTENT(IN)                      :: ncore
 INTEGER, INTENT(OUT)                     :: lused
 INTEGER, INTENT(IN OUT)                  :: ixytyp
 INTEGER, INTENT(IN OUT)                  :: icase
 INTEGER, INTENT(IN OUT)                  :: buff(1)
 LOGICAL, INTENT(OUT)                     :: anyxy
 INTEGER :: loc(6), xycdb
 
 
!/////
 COMMON/system/ sysbuf, iout
 COMMON/names / rd, rdrew, wrt, wrtrew, clsrew, cls
 COMMON/ddrmc1/ dummy(362), ierror
!/////
 
 DATA xycdb/ 108 /, noeor / 0 /
 
 lused = 0
 anyxy = .false.
 CALL OPEN(*100,xycdb,buff,rdrew)
 CALL fwdrec(*300,xycdb)
 CALL fwdrec(*300,xycdb)
 
!     FIND ENTRIES IN SUBCASE 0 OF THIS TYPE IF ANY.
 
 5 CALL READ(*300,*300,xycdb,loc,6,noeor,nwds)
 IF( loc(1)  > 0) THEN
   GO TO    20
 END IF
 10 IF( loc(2) /= ixytyp ) GO TO 5
 
!     SAVE ID IN TABLE
 
 IF( lused  > 0) THEN
   GO TO    12
 ELSE
   GO TO    11
 END IF
 
!      ADD TO LIST IF NOT A REPEAT ID
 
 12 IF( loc(3) == z(lused) ) GO TO 5
 11 lused = lused + 1
 IF( lused > ncore ) GO TO 1000
 z(lused) = loc(3)
 GO TO 5
 
!     FIND ENTRIES IN SUBCASE -ICASE- OF THIS TYPE IF ANY EXIST.
 
 15 CALL READ(*300,*300,xycdb,loc,6,noeor,nwds)
 20 IF( loc(1) - icase  < 0) THEN
   GO TO    15
 ELSE IF ( loc(1) - icase  == 0) THEN
   GO TO    30
 ELSE
   GO TO   300
 END IF
 30 IF( loc(2) - ixytyp  < 0) THEN
   GO TO    15
 ELSE IF ( loc(2) - ixytyp  == 0) THEN
   GO TO    40
 ELSE
   GO TO   300
 END IF
 40 lused = lused + 1
 IF( lused > ncore ) GO TO 1000
 z(lused) = loc(3)
 GO TO 15
 
!     LIST IS NOW COMPLETE THUS SORT IT, AND REMOVE REPEATED IDS.
 
 300 CALL CLOSE( xycdb, clsrew )
 IF( lused  > 0) THEN
   GO TO   301
 ELSE
   GO TO   100
 END IF
 301 CALL sort( 0, 0, 1, 1, z(1), lused )
 anyxy = .true.
 
 j = 1
 IF( lused == 1 ) GO TO 305
 DO  i = 2,lused
   IF( z(i) == z(j) ) CYCLE
   j = j + 1
   z(j) = z(i)
 END DO
 
 305 lused = j
 100 RETURN
 
!     INSUFFICIENT CORE ALTERNATE RETURN.
 
 1000 ierror = 859
 RETURN 1
END SUBROUTINE ddrmmp
