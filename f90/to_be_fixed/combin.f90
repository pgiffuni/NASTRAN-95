SUBROUTINE combin (pg,ilist,nlist)
     
 
 INTEGER, INTENT(IN OUT)                  :: pg
 INTEGER, INTENT(IN)                      :: ilist(1)
 INTEGER, INTENT(IN)                      :: nlist
 INTEGER :: sysbuf, NAME(2),hcflds,hcfld,hccens,hccen,otpe,  &
     remfls,remfl,mcb(7)
 DIMENSION       ary(1), alpha(360),loadn(360),loadnn(360),  &
     iary(1),alpha1(360),lodc1(7),head(2)
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /loadx / lc,n(13),lodc,mass
 COMMON /BLANK / nrowsp
 COMMON /system/ sysbuf,otpe,dum52(52),iprec
 COMMON /zzzzzz/ core(1)
 COMMON /loads / nload,iptr
 COMMON /zntpkx/ a(4),ll,ieol,IEOR
 COMMON /packx / ita,itb,ii,jj,incur
 EQUIVALENCE     (core(1),iary(1),ary(1))
 
!     ALSO COMBINE HCFLD AND REMFL IN MAGNETOSTATIC PROBLEMS
 
 DATA    hcflds, hcfld /304,202/
 DATA    remfls, remfl /305,203/
 DATA    hccens, hccen /307,204/
 DATA    NAME  / 4HCOMB,4HIN   /
 
 
 ita = 1
 itb = iprec
 ii  = 1
 
!     PERFORM CHECKS IN E AND M PROBLEM
!     IN E AND M PROBLEM, REMFLS AND HCFLDS MUST HAVE THE SAME NUMBER
!     OF COLUMNS AS PG
 
 mcb(1) = remfls
 CALL rdtrl (mcb)
 nperms = 0
 IF (mcb(1) <= 0) GO TO 1
 nperms = mcb(2)
 1 mcb(1) = hcflds
 CALL rdtrl (mcb)
 nhc = 0
 IF (mcb(1) <= 0) GO TO 2
 nhc = mcb(2)
 2 IF (nhc /= nperms) GO TO 300
 IF (nhc == 0) GO TO 5
 mcb(1) = pg
 CALL rdtrl (mcb)
 IF (nhc /= mcb(2)) GO TO 300
 5 CONTINUE
 mcb(1) = hccens
 CALL rdtrl (mcb)
 ns = 0
 IF (mcb(1) <= 0) GO TO 6
 ns = mcb(2)
 6 IF (ns /= nhc) GO TO 300
 jj    = nrowsp
 incur = 1
 lcore = lc
 ibuf1 = lcore
 lcore = lcore - sysbuf
 CALL OPEN (*200,lodc,core(lcore+1),1)
 CALL fname (lodc,head)
 CALL WRITE (lodc,head,2,1)
 lcore = lcore - sysbuf
 CALL OPEN (*190,pg,core(lcore+1),0)
 CALL makmcb (lodc1,lodc,nrowsp,2,iprec)
 nlj = iptr
 nl1 = 0
 DO  i = 1,nload
   DO   j = 1,nrowsp
     core(j) = 0.0
   END DO
   nlj = nlj + nl1*2 + 1
   nl1 = iary(nlj)
   DO  k = 1,nl1
     kk  = nlj + (k-1)*2 + 1
     loadn(k) = iary(kk)
     IF (loadn(k) < 0) GO TO 150
     alpha(k) = ary(kk+1)
   END DO
   kk = 1
   kl = 0
   DO  k = 1,nlist
     IF (ilist(k) == 0) THEN
       GO TO    60
     END IF
     30 kl = kl + 1
     DO  j = 1,nl1
       IF (loadn(j)-ilist(k) == 0) THEN
         GO TO    50
       ELSE
         GO TO    40
       END IF
     END DO
     CYCLE
     50 loadnn(kk) = kl
     alpha1(kk) = alpha(j)
     kk = kk + 1
   END DO
   kk = 1
   DO  j = 1,nl1
     inull = 0
     IF (j /= 1) GO TO 70
     CALL skprec (pg,1)
     70 CALL intpk (*120,pg,0,1,0)
     80 IF (loadnn(j)-kk == 0) THEN
       GO TO   100
     END IF
     90 IF (inull == 1) GO TO 91
     IF (IEOR  == 0) CALL skprec (pg,1)
     91 CONTINUE
     kk = kk + 1
     inull = 0
     GO TO 70
     100 IF (inull == 1) GO TO 130
     IF (ieol == 0) THEN
       GO TO   110
     ELSE
       GO TO   130
     END IF
     110 CALL zntpki
     core(ll) = core(ll) + a(1)*alpha1(j)
     GO TO 100
     120 inull = 1
     GO TO 80
     130 kk = kk + 1
   END DO
   150 CALL pack (core,lodc1(1),lodc1)
   CALL REWIND (pg)
 END DO
 CALL wrttrl (lodc1(1))
 CALL CLOSE  (lodc1(1),1)
 CALL CLOSE  (pg,1)
 IF (pg == hcflds) GO TO 170
 IF (pg == remfls) GO TO 180
 IF (pg == hccens) RETURN
 
!     DO MAGNETOSTATIC FIELDS FOR USE IN EMFLD
 
 lodc1(1) = hcflds
 CALL rdtrl (lodc1)
 
!     IF HCFLD IS PURGED, SO MUST REMFLS
 
 IF (lodc1(2) <= 0) RETURN
 pg   = hcflds
 lodc = hcfld
 nrowsp = 3*nrowsp
 GO TO 5
 
!     DO REMFLS
 
 170 lodc1(1) = remfls
 CALL rdtrl (lodc1)
 IF (lodc1(2) <= 0) RETURN
 pg   = remfls
 lodc = remfl
 nrowsp = lodc1(3)
 GO TO 5
 
!     HCCENS
 
 180 lodc1(1) = hccens
 CALL rdtrl (lodc1)
 IF (lodc1(2) <= 0) RETURN
 pg   = hccens
 lodc = hccen
 nrowsp = lodc1(3)
 GO TO 5
 190 ip1 = pg
 195 CALL mesage (-1,ip1,NAME)
 200 IF (lodc == hcfld) RETURN
 ip1 = lodc
 GO TO 195
 300 WRITE  (otpe,350) ufm
 350 FORMAT (a23,', IN AN E AND M PROBLEM, SCRATCH DATA BLOCKS HCFLDS',  &
     ' AND REMFLS HAVE DIFFERENT NUMBERS OF COLUMNS.', /10X,  &
     ' THIS MAY RESULT FROM SPCFLD AND REMFLU CARDS HAVING THE',  &
     ' SAME LOAD SET ID')
 CALL mesage (-61,0,0)
 RETURN
END SUBROUTINE combin
