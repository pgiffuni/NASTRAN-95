SUBROUTINE cmmcon (nce)
     
!     THIS SUBROUTINE DETERMINES WHETHER MORE THAN ONE CONNECTION ENTRY
!     HAS BEEN SPECIFIED FOR A GIVEN IP NUMBER.
 
 
 INTEGER, INTENT(OUT)                     :: nce
 LOGICAL :: mcon
 INTEGER :: scconn,buf1,z,score,scmcon,buf3,aaa(2)
 COMMON /cmb001/ scr1,scr2,scbdat,scsfil,scconn,scmcon
 COMMON /cmb002/ buf1,buf2,buf3,junk(2),score,lcore,inpt,outt
 COMMON /cmb003/ junk2(38),npsub,junk3(2),mcon
 COMMON /zzzzzz/ z(1)
 DATA    aaa   / 4HCMMC,4HON   /
 
!     READ CONNECTION ENTRIES INTO OPEN CORE
 
 nwd   = 2 + npsub
 mcon  = .true.
 ifile = scconn
 CALL OPEN (*700,scconn,z(buf1),0)
 j   = 0
 nce = 0
 90 CALL READ (*200,*100,scconn,z(score+j),10,1,nnn)
 100 nce = nce + 1
 z(score+j) = nce
 j = j + nwd
 GO TO 90
 200 CALL CLOSE (scconn,1)
 
!     SWEEP THROUGH CONNECTION ENTRIES AND DETERMINE THOSE THAT
!     REPRESENT MULTIPLE CONNECTIONS.
 
 mcon  = .false.
 ncem1 = nce - 1
 
 DO  k = 1,ncem1
   DO  i = 1,npsub
     ist = score + i + (k-1)*nwd + 1
     IF (z(ist) == 0) CYCLE
     DO  j = 1,nce
       IF (k == j) CYCLE
       isub = score + 1 + i + (j-1)*nwd
       IF (z(ist) /= z(isub)) CYCLE
       iloc = i + 1
       z(ist -iloc) = -1*IABS(z(ist -iloc))
       z(isub-iloc) = -1*IABS(z(isub-iloc))
       mcon = .true.
     END DO
   END DO
 END DO
 
 IF (.NOT.mcon) RETURN
 
!     GENERATE OUTPUT FILE OF CONNECTION ENTRY IDS
 
 ifile = scmcon
 CALL OPEN (*700,scmcon,z(buf1),1)
 DO  i = 1,nce
   loc = score + (i-1)*nwd
   IF (z(loc) < 0) CALL WRITE (scmcon,IABS(z(loc)),1,0)
 END DO
 CALL WRITE (scmcon,0,0,1)
 CALL CLOSE (scmcon,1)
 RETURN
 
 700 CALL mesage (-1,ifile,aaa)
 RETURN
END SUBROUTINE cmmcon
