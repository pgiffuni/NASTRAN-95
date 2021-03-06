SUBROUTINE errmkn (n,ierr)
     
!     SENDS ERROR MESSAGES.  N IS THE INDEX OF THE SUBROUTINE CALLING
!     ERROR, AND IERR IS AN ERROR CODE.
 
 
 INTEGER, INTENT(IN OUT)                  :: n
 INTEGER, INTENT(IN OUT)                  :: ierr
 DIMENSION       isubr(26)
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm,uim,sfm
 COMMON /system/ nbuff,nout
 DATA    isubr / 4HCRSU,4HB   ,4HDSTR,4HOY  ,4HFDIT,4H    ,  &
     4HFMDI,4H    ,4HFNXT,4H    ,4HGETB,4HLK  ,  &
     4HRETB,4HLK  ,4HSETE,4HQ   ,4HSJUM,4HP   ,  &
     4HSURE,4HAD  ,4HRENA,4HME  ,4HEXO2,4H    , 4HEXIO,4H1   /
 
 WRITE (nout,1000) sfm,isubr(n),isubr(n+1)
 
 CALL sofcls
 SELECT CASE ( ierr )
   CASE (    1)
     GO TO 10
   CASE (    2)
     GO TO 20
   CASE (    3)
     GO TO 30
   CASE (    4)
     GO TO 40
   CASE (    5)
     GO TO 50
   CASE (    6)
     GO TO 60
   CASE (    7)
     GO TO 70
   CASE (    8)
     GO TO 80
   CASE (    9)
     GO TO 90
   CASE (   10)
     GO TO 100
 END SELECT
 10 WRITE (nout,1010)
 GO TO 900
 20 WRITE (nout,1020)
 GO TO 900
 30 WRITE (nout,1030)
 GO TO 900
 40 WRITE (nout,1040)
 GO TO 900
 50 WRITE (nout,1050)
 GO TO 900
 60 WRITE (nout,1060)
 GO TO 900
 70 WRITE (nout,1070)
 GO TO 900
 80 WRITE (nout,1080)
 GO TO 900
 90 WRITE (nout,1090)
 GO TO 900
 100 WRITE (nout,1100)
 GO TO 900
 900 CALL mesage (-61,0,0)
 RETURN
 
 1000 FORMAT (a25,' 6224, SOF UTILITY SUBROUTINE ',2A4)
 1010 FORMAT (5X,'I IS TOO LARGE OR NXTTSZ HAS NOT BEEN PROPERLY ', 'UPDATED')
 1020 FORMAT (5X,'ILLEGAL BLOCK NUMBER')
 1030 FORMAT (5X,'ERROR IN SETTING UP THE LIST IMORE')
 1040 FORMAT (5X,'NXTCUR IS TOO LARGE')
 1050 FORMAT (5X,'ERROR IN UPDATING DIT')
 1060 FORMAT (5X,'ERROR IN UPDATING MDI')
 1070 FORMAT (5X,'ERROR IN LINKING BLOCKS OF DIT')
 1080 FORMAT (5X,'LINK THROUGH COMBINED SUBSTRUCTURES IS NOT CIRCULAR')
 1090 FORMAT (5X,'ERROR IN LINKING SOF BLOCKS')
 1100 FORMAT (5X,'INTERNAL ARRAY DIMENSION EXCEEDED')
END SUBROUTINE errmkn
