SUBROUTINE preloc (*,buf,FILE)
     
!     PRELOC OPENS AND POSITIONS REQUESTED FILE TO FIRST DATA RECORD.
!     LOCATE POSITIONS FILE TO REQUESTED DATA RECORD WITHIN FILE.
 
 
 , INTENT(IN OUT)                         :: *
 INTEGER, INTENT(OUT)                     :: buf(2)
 INTEGER, INTENT(IN)                      :: FILE
 EXTERNAL        andf
 INTEGER :: nam(2),trl(7),nm1(2),flag  ,andf  ,  &
     two   ,ret   ,bff(1),id(2) ,flg
 COMMON /two   / two(32)
 DATA    nam   ,        nm1          / 4HPREL, 4HOC  ,4HLOCA,4HTE  /
 
!     OPEN FILE. IF PURGED, GIVE ALTERNATE RETURN.
!     OTHERWISE SKIP HEADER RECORD
 
 trl(1) = FILE
 CALL rdtrl (trl)
 IF (trl(1) < 0) GO TO 10
 IF (trl(2)+trl(3)+trl(4)+trl(5)+trl(6)+trl(7) == 0) GO TO 10
 CALL OPEN (*10,FILE,buf(2),0)
 CALL fwdrec (*2,FILE)
 buf(1) = FILE
 icheck = 123456789
 RETURN
 10 RETURN 1
 
!     FATAL FILE ERRORS
 
 2 CALL mesage (-2,FILE,nam)
 3 CALL mesage (-3,trl,nm1)
 
 
 ENTRY locate (*,bff,id,flg)
!     ===========================
 
!     ENTRY TO POSITION DATA RECORD.
 
!     READ TRAILER FOR FILE. IF BIT NOT ON OR FILE PURGED,
!     GIVE ALTERNATE RETURN.
 
!WKBD IF (ICHECK .NE. 123456789) CALL ERRTRC ('LOCATE  ',10)
 trl(1) = bff(1)
 CALL rdtrl (trl)
 IF (trl(1) < 0) RETURN 1
 k = (id(2)-1)/16
 l =  id(2)- 16*k
 IF (andf(trl(k+2),two(l+16)) == 0) RETURN 1
 
!     READ THREE ID WORDS FROM DATA RECORD.
!     IF END-OF-FILE, REPOSITION FILE TO FIRST DATA RECORD AND RETRY.
!     IF ID WORD MATCHES USER, RETURN.
 
 last = 0
 ASSIGN 20 TO ret
 20 CALL READ (*50,*20,trl(1),trl(2),3,0,flag)
 IF (trl(2) /= id(1)) GO TO 22
 flg = trl(4)
 RETURN
 
!     SKIP RECORD. READ ID WORDS FROM NEXT RECORD. IF MATCH,RETURN.
!     IF END-OF FILE, REPOSITION TO FIRST DATA RECORD AND RETRY.
!     IF NO MATCH, TEST FOR RETURN TO ORIGINAL FILE POSITION. IF SO,
!     QUEUE MESSAGE AND GIVE ALTERNATE RETURN. IF NOT, CONTINUE SEARCH.
 
 22 ASSIGN 30 TO ret
 25 CALL fwdrec (*2,trl(1))
 30 CALL READ (*50,*3,trl(1),trl(5),3,0,flag)
 IF (trl(5) /= id(1)) GO TO 32
 flg = trl(7)
 RETURN
 
 32 IF (trl(5) /= trl(2)) GO TO 25
 35 CALL sswtch (40,j)
!WKBD IF (J .NE. 0) CALL ERRTRC ('LOCATE  ',35)
 CALL mesage (30,72,id)
 CALL fwdrec (*2,trl(1))
 RETURN 1
 
!     CODE TO POSITION FILE TO FIRST DATA RECORD.
 
 50 CALL REWIND (trl(1))
 IF (last /= 0) GO TO 35
 last = 1
 CALL fwdrec (*2,trl(1))
 GO TO ret, (20,30)
END SUBROUTINE preloc
