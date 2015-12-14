SUBROUTINE rename (name1,name2,z,nz,itest)
     
!     THIS ROUTINE RENAMES SUBSTRUCTURE NAME1 TO NAME2. SOF ITEMS EQSS,
!     BGSS, CSTM, LODS, LOAP AND PLTS ARE REWRITTEN TO REFLECT THE NEW
!     NAME.  THESE ITEMS ARE CHANGED FOR NAME1 AND ANY HIGHER LEVEL
!     SUBSTRUCTURE FOR WHICH NAME1 IS A COMPONENT.  NO CHANGES ARE MADE
!     TO SECONDARY SUBSTRUCTURES OF NAME1 WHICH RETAIN THEIR ORIGINAL
!     NAMES.  ALSO NO CHANGES ARE MADE TO THE SOLUTION DATA (ITEM SOLN)
!     FOR SUBSTRUCTURE NAME1 OR ANY HIGHER LEVEL SUBSTRUCTURES.
 
!     VALUES RETURNED IN ITEST ARE
!        1 - NORMAL RETURN
!        4 - SUBSTRUCTURE NAME1 DOES NOT EXIST
!       10 - SUBSTRUCTURE NAME2 ALREADY EXISTS
 
 
 INTEGER, INTENT(IN)                      :: name1(2)
 INTEGER, INTENT(IN)                      :: name2(2)
 INTEGER, INTENT(OUT)                     :: z(2)
 INTEGER, INTENT(IN)                      :: nz
 INTEGER, INTENT(OUT)                     :: itest
 EXTERNAL        andf
 LOGICAL :: ditup    ,mdiup    ,higher
 INTEGER :: NAME(2)  ,nameh(2)  , eog      ,BLANK    ,ps        ,andf    ,  &
     namsub(2)
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm      ,uwm      ,uim
 COMMON /itemdt/ nitem    ,items(7,1)
 COMMON /sof   / ditdum(6),iodum(8) ,mdidum(4),nxtdum(15),ditup   , mdiup
 COMMON /system/ sysbuf, nout
 COMMON /zzzzzz/ buf(1)
 DATA    ps    / 1 /
 DATA    eog   / 4H$eog/, BLANK /4H    /
 DATA    namsub/ 4HRENA,4HME  /
 
 
!     CHECK IF NAME2 ALREADY EXISTS
 
 CALL fdsub (name2,ind)
 IF (ind /= -1) GO TO 1000
 
!     CHANGE DIT ENTRY FOR SUBSTRUCTURE NAME1
 
 CALL fdsub (name1,ind)
 IF (ind < 0) GO TO 1100
 CALL fdit (ind,idit)
 IF (name1(1) /= buf(idit) .OR. name1(2) /= buf(idit+1)) GO TO 1200
 buf(idit  ) = name2(1)
 buf(idit+1) = name2(2)
 ditup   = .true.
 
 NAME(1) = name2(1)
 NAME(2) = name2(2)
 higher  = .false.
 
!     CHANGE TABLE ITEMS WHICH CONTAIN SUBSTRUCTRUE NAME
!     SUBSTRUCTURE NAME.
!     HIGHER = .FALSE. - WE ARE WORKING WITH SUBSTRUCTURE NAME1
!     HIGHER = .TRUE.  - WE ARE WORKING WITH A SUBSTRUCTURE FOR
!                        WHICH NAME1 IS A COMPONENT
 
 10 CALL fdsub (NAME,ind)
 CALL fmdi (ind,imdi)
 ips = andf(buf(imdi+ps),1023)
 DO  itm = 1,nitem
   IF (items(2,itm) > 0) CYCLE
   item = items(1,itm)
   inum = items(3,itm)/1000000
   iloc = (items(3,itm) - inum*1000000)/1000
   incr =  items(3,itm) - inum*1000000 - iloc*1000
   
!     PROCESS THE FOLLOWING ITEMS
   
!     SUBSTRUCTRUE NAME1
!     DONT PROCESS THE ITEM IF THIS IS A SECONDARY SUBSTRUCTURE
!     AND THE ACTUAL ITEM IS STORED FOR THE PRIMARY (I.E. BGSS,CSTM(
!     HIGHER LEVEL SUBSTRUCTURE
!     DONT PROCESS THE ITEM IF IT IS ACTUALLY STORED FOR THE
!     PRIMARY SUBSTRUCTURE (I.E. BGSS,CSTM)
   
   IF (.NOT.higher .AND. iloc == 0 .AND. ips /= 0) CYCLE
   IF (higher .AND. iloc == 0) CYCLE
   
!     READ ITEM INTO OPEN CORE
   
   irw = 1
   CALL sfetch (NAME,item,irw,itest)
   IF (itest /= 1) CYCLE
   ncore = nz
   icore = 1
   20 CALL suread (z(icore),ncore,nwds,itest)
   IF (itest == 3) GO TO 30
   IF (itest == 1) GO TO 1300
   z(icore+nwds) = eog
   icore = icore + nwds + 1
   ncore = ncore - nwds - 1
   GO TO 20
   30 nwds  = icore + nwds - 1
   
!     CHANGE ANY OCCURANCE OF NAME1 WITH NAME2
   
   IF (higher) GO TO 40
   
!     SUBSTRUCTURE NAME1 - NAME SHOULD BE IN WORDS 1 AND 2 OF GROUP 0
   
   IF (z(1) /= name1(1) .OR. z(2) /= name1(2)) CYCLE
   z(1) = name2(1)
   z(2) = name2(2)
   
!     SEARCH THE LIST OF COMPONENT SUBSTRUCTURES FOR NAME1
   
   40 ncomp = z(inum)
   iloc2 = iloc + incr*ncomp - 1
   DO  i = iloc,iloc2,incr
     IF (z(i) /= name1(1) .OR. z(i+1) /= name1(2)) CYCLE
     z(i  ) = name2(1)
     z(i+1) = name2(2)
     GO TO 60
   END DO
   
!     DELETE OLD ITEM
   
   60 CALL DELETE (NAME,item,itest)
   
!     WRITE NEW ITEM TO SOF
   
   itest = 3
   irw   = 2
   CALL sfetch (NAME,item,irw,itest)
   itest = 3
   CALL suwrt (z(1),nwds,itest)
   
 END DO
 
!     GET NEXT HIGHER LEVEL SUBSTRUCTURE FOR WHICH NAME1 IS A
!     COMPONENT AND PERFORM SAME PROCEDURE
 
 CALL fndnxl (NAME,nameh)
 IF (nameh(1) == BLANK .OR. nameh(1) == NAME(1) .AND.  &
     nameh(2) == NAME(2)) GO TO 110
 NAME(1) = nameh(1)
 NAME(2) = nameh(2)
 higher = .true.
 GO TO 10
 
!     NO HIGHER LEVEL SUBSTRUCTURES LEFT - PRINT INFORMATION MESSAGE
!     AND RETURN
 
 110 WRITE  (nout,120) uim,name1,name2
 120 FORMAT (a29,' 6229, SUBSTRUCTURE ',2A4,' HAS BEEN RENAMED TO ', 2A4)
 itest = 1
 RETURN
 
!     ERROR RETURNS
 
 
!     SUBSTRUCTURE NAME2 ALREADY EXIST ON THE SOF
 
 1000 WRITE  (nout,1010) uwm,name1,name2
 1010 FORMAT (a25,' 6230, SUBSTRUCTURE ',2A4,' HAS NOT BEEN RENAMED ',  &
     'BECAUSE ',2A4,' ALREADY EXISTS ON THE SOF.')
 itest = 10
 RETURN
 
!     SUBSTRUCTURURE NAME1 DOES NOT EXIST
 
 1100 itest = 4
 RETURN
 
!     DIT FORMAT ERROR
 
 1200 CALL errmkn (21,5)
 
!     INSUFFICIENT CORE TO HOLD ITEM
 
 1300 CALL sofcls
 CALL mesage (-8,0,namsub)
 RETURN
END SUBROUTINE rename
