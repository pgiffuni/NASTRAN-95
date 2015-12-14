SUBROUTINE procom (procos,procof,casecc,ncoefs,ngrids)
     
!     PROCOM COMBINES PROCOF CASES FOR SUBCOM-S AND REPCASES
 
 
 INTEGER, INTENT(IN)                      :: procos
 INTEGER, INTENT(IN)                      :: procof
 INTEGER, INTENT(IN)                      :: casecc
 INTEGER, INTENT(IN)                      :: ncoefs
 INTEGER, INTENT(IN OUT)                  :: ngrids
 INTEGER :: buf1,buf2,buf3,FILE, info(7), iz(1),nam(2)
 COMMON /system/ ibuf
 COMMON /zzzzzz/ z(1)
 EQUIVALENCE     (z(1),iz(1))
 DATA    i166  , i16  ,nam   / 166, 16, 4HPROC,4HOM  /
 
 lcore = korsz(z)
 buf1  = lcore - ibuf + 1
 buf2  = buf1  - ibuf
 buf3  = buf2  - ibuf
 lcore = buf3  - 1
 IF (lcore < ncoefs .OR. lcore < ngrids) GO TO 108
 CALL gopen (procos,z(buf1),0)
 CALL gopen (procof,z(buf2),1)
 
!     CHECK EACH SUBCASE FOR REPCASE OR SUBCOM-IF NONE(JUST COPY SET OF
!     5 RECORDS FROM PROCOS TO PROCOF
 
 FILE = casecc
 CALL gopen (casecc,z(buf3),0)
 10 FILE = casecc
 CALL READ (*90,*20,casecc,z(1),lcore,0,iwords)
 GO TO 108
 20 IF (iz(i16) /= 0) GO TO 30
 
!     NOT A SUBCOM - MIGHT BE REPCASE
 
 25 FILE = procos
 CALL fread (procos,z,103,1)
 CALL WRITE (procof,z,103,1)
 CALL fread (procos,z,ncoefs,1)
 CALL WRITE (procof,z,ncoefs,1)
 CALL fread (procos,z,ncoefs,1)
 CALL WRITE (procof,z,ncoefs,1)
 CALL fread (procos,z,ngrids,1)
 CALL WRITE (procof,z,ngrids,1)
 CALL fread (procos,z,ngrids,1)
 CALL WRITE (procof,z,ngrids,1)
 
!     GO BACK FOR ANOTHER CASE CONTROL RECORD
 
 GO TO 10
 
!     REPCASE OR SUBCOM
 
 30 IF (iz(i16) > 0) GO TO 45
 
!     REPCASE
 
 DO  i = 1,5
   CALL bckrec (procos)
 END DO
 GO TO 25
 
!     SUBCOM
 
 45 lcc  = iz(i166)
 lsym = iz(lcc)
 DO  i = 1,lsym
   DO  j = 1,5
     CALL bckrec (procos)
   END DO
 END DO
 ntot = 2*(ncoefs+ngrids)
 IF (iwords+2*ntot > lcore) GO TO 108
 inew = iwords + ntot
 DO  i = 1,ntot
   z(inew+i) = 0.
 END DO
 DO  i = 1,lsym
   coef = z(lcc+i)
   IF (coef == 0.) GO TO 75
   CALL fread (procos,info,103,1)
   CALL fread (procos,z(iwords+1),ncoefs,1)
   CALL fread (procos,z(iwords+ncoefs+1),ncoefs,1)
   CALL fread (procos,z(iwords+2*ncoefs+1),ngrids,1)
   CALL fread (procos,z(iwords+2*ncoefs+ngrids+1),ngrids,1)
   DO  j = 1,ntot
     z(inew+j) = z(inew+j) + coef*z(iwords+j)
   END DO
   CYCLE
   75 DO  k = 1,5
     CALL fwdrec (*102,procos)
   END DO
   
 END DO
 
!     WRITE TO PROCOF- 1ST BE SURE THAT ISYM IS 0 TO ACCOUNT FOR
!     POSSIBLE SYMMETRY-ANTISYMMETRY COMBINATION
 
 info(6) = 0
 CALL WRITE (procof,info,103,1)
 CALL WRITE (procof,z(inew+1),ncoefs,1)
 CALL WRITE (procof,z(inew+ncoefs+1),ncoefs,1)
 CALL WRITE (procof,z(inew+2*ncoefs+1),ngrids,1)
 CALL WRITE (procof,z(inew+2*ncoefs+ngrids+1),ngrids,1)
 
!     GO BACK FOR ANOTHER SUBCASE
 
 GO TO 10
 
!     DONE
 
 90 CALL CLOSE (casecc,1)
 CALL CLOSE (procos,1)
 CALL CLOSE (procof,1)
 info(1) = procos
 CALL rdtrl (info)
 info(1) = procof
 CALL wrttrl (info)
 RETURN
 
 102 CALL mesage (-2,0,nam)
 108 CALL mesage (-8,0,nam)
 RETURN
END SUBROUTINE procom
