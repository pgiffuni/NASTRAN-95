SUBROUTINE cmckcd
     
!     THIS SUBROUTINE DETERMINES WHETHER MANUALLY SPECIFIED CONNECTION
!     ENTRIES ARE ALLOWABLE BASED ON THE PRESCRIBED GEOMETRIC TOLERANCE.
 
 INTEGER :: scsfil,combo,score,ist(7),scconn,ce(9),aaa(2), buf2,outt
 DIMENSION       ipnum(7),coord(7,3),diff2(3)
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm,uim
 COMMON /cmb001/ scr1,scr2,scbdat,scsfil,scconn,scmcon,sctoc,  &
     geom4,casecc,sccstm,scr3
 COMMON /cmb002/ buf1,buf2,buf3,buf4,buf5,score,lcore,intp,outt
 COMMON /cmb003/ combo(7,5),conset,iauto,toler,npsub,conect,tran,  &
     mcon,restct(7,7),isort,origin(7,3),iprint
 COMMON /zzzzzz/ z(1)
 COMMON /BLANK / step,idry
 DATA    aaa   / 4HCMCK,4HCD   /
 
!     READ ALL BGSS INTO OPEN CORE
 
 it    = 2
 ierr  = 0
 llco  = lcore
 j     = 0
 ifile = scsfil
 CALL OPEN (*200,scsfil,z(buf2),0)
 DO  i = 1,npsub
   nrec  = combo(i,5) + 1
   DO  jj = 1,nrec
     CALL fwdrec (*210,scsfil)
   END DO
   CALL READ (*210,*20,scsfil,z(score+j),llco,1,nnn)
   GO TO 220
   20 ist(i) = score + j
   j = j + nnn
   llco  = llco - nnn
   CALL skpfil (scsfil,1)
 END DO
 CALL CLOSE (scsfil,1)
 
!     READ CONNECTION ENTRIES AND LOAD INTO COORD ARRAY
 
 ifile = scconn
 CALL OPEN (*200,scconn,z(buf2),0)
 40 CALL READ (*180,*50,scconn,ce,10,1,nnn)
 
!     LOAD COORD ARRAY
!     CE(3)... UP TO CE(9) ARE INTERNAL POINT NO.
!     IZ(IADD) IS THE COORD (CSTM) ID OF THE INTERNAL PTS.
!     Z(IADD+1,+2,+3) ARE THE COORD. ORIGINS
 
 50 npt  = 0
 DO  i = 1,npsub
   IF (ce(i+2) > 0.0) THEN
     GO TO    60
   ELSE
     GO TO    80
   END IF
   60 npt  = npt + 1
   iadd = 4*(ce(i+2)-1) + ist(i)
   ipnum(npt) = ce(i+2)
   DO  j = 1,3
     coord(npt,j) = z(iadd+j)
   END DO
 END DO
 
!     COMPARE ALL PAIRS OF COORDINATES AGAINST TOLER.
 
 nptm1 = npt - 1
 DO  i = 1,nptm1
   it = it - 1
   jj = i  + 1
   DO  j = jj,npt
     DO  kk = 1,3
       diff2(kk) = (coord(j,kk)-coord(i,kk))**2
     END DO
     sum  = 0.0
     DO  kk = 1,3
       sum  = sum + diff2(kk)
     END DO
     dist = SQRT(sum)
     IF (dist <= toler) CYCLE
     IF (it > 1) GO TO 120
     WRITE  (outt,110) ufm
     110 FORMAT (a23,' 6514, ERRORS HAVE BEEN FOUND IN MANUALLY SPECIFIED',  &
         ' CONNECTION ENTRIES. SUMMARY FOLLOWS')
     ierr = 1
     idry =-2
     it   = 2
     120 IF (it > 2) GO TO 140
     WRITE  (outt,130) (ce(kdh),kdh=1,nnn)
     130 FORMAT ('0*** GEOMETRIC ERRORS HAVE BEEN FOUND IN THE FOLLOWING',  &
         ' CONNECTION ENTRY', /5X,9I10)
     it   = 3
     140 WRITE (outt,150) ipnum(i),(coord(i,mm),mm=1,3),  &
         ipnum(j),(coord(j,mm),mm=1,3)
     150 FORMAT ('0*** IP NUMBER',i10,13H coordinates  ,3E16.6,4H AND, /,  &
         '     IP NUMBER',i10,13H coordinates  ,3E16.6,  &
         ' ARE NOT WITHIN TOLER UNITS.')
   END DO
 END DO
 GO TO 40
 
 180 IF (ierr == 0) WRITE (outt,190) uim
 190 FORMAT (a29,' 6516, ALL MANUAL CONNECTIONS SPECIFIED ARE ',  &
     'ALLOWABLE WITH RESPECT TO TOLERANCE')
 CALL CLOSE (scconn,1)
 GO TO 250
 
 200 imsg = -1
 GO TO 230
 210 imsg = -2
 GO TO 230
 220 imsg = -8
 230 CALL mesage (imsg,ifile,aaa)
 
 250 RETURN
END SUBROUTINE cmckcd
