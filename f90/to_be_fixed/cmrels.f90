SUBROUTINE cmrels
     
!     THIS SUBROUTINE ENFORCES THE RELES DATA SPECIFIED FOR THE
!     COMB1 MODULE.
 
 EXTERNAL        andf
 LOGICAL :: first
 INTEGER :: ix(7,3),scbdat,z,score,buf1,buf2,scconn,ps1,ps2
 INTEGER :: list(32),andf,stce,aaa(2)
 COMMON /cmb001/ scr1,scr2,scbdat,scsfil,scconn
 COMMON /cmb002/ buf1,buf2,buf3,buf4,buf5,score,lcore,inpt,outt
 COMMON /cmb003/ junk(38),npsub
 COMMON /zzzzzz/ z(1)
 DATA    aaa   / 4HCMRE,4HLS   /
 
 ifile = scbdat
 kj = 0
 DO  i = 1,7
   DO  j = 1,3
     ix(i,j) = 0
   END DO
 END DO
 DO  i = 1,npsub
   first = .true.
   CALL OPEN (*150,scbdat,z(buf1),0)
   CALL skpfil (scbdat,3)
   30 CALL READ (*60,*170,scbdat,id,1,0,n)
   IF (id == i) GO TO 40
   CALL fwdrec (*60,scbdat)
   GO TO 30
   40 CALL READ (*160,*50,scbdat,z(score+kj),lcore,1,nw)
   GO TO 180
   50 IF (first) ix(i,2) = score + kj
   first = .false.
   ix(i,3) = ix(i,3) + nw/2
   kj = kj + nw
   lcore = lcore - nw
   ix(i,1) = 1
   GO TO 30
   60 CALL CLOSE (scbdat,1)
 END DO
 DO  i = 1,npsub
   IF (ix(i,1) == 0) CYCLE
   ist = ix(i,2)
   nw  = ix(i,3)*2
   CALL sort (0,0,2,1,z(ist),nw)
 END DO
 ifile = scconn
 CALL OPEN (*150,scconn,z(buf2),0)
 nwrd = 2 + npsub
 nce  = 0
 stce = score + kj
 90 CALL READ (*110,*100,scconn,z(score+kj),lcore,1,nnn)
 GO TO 180
 100 kj  = kj + nwrd
 nce = nce + 1
 GO TO 90
 110 CALL CLOSE (scconn,1)
 nce = nwrd*nce
 DO  i = 1,nce,nwrd
   ii = i - 1
   icode = z(stce+ii+1)
   CALL decode (icode,list,nc)
   IF (nc /= 2) CYCLE
   ps1  = list(1) + 1
   ps2  = list(2) + 1
   ist1 = ix(ps1,2)
   ist2 = ix(ps2,2)
   nw1  = ix(ps1,3)
   nw2  = ix(ps2,3)
   IF (ix(ps1,1) == 0) GO TO 120
   kid  = z(stce+ii+1+ps1)
   CALL bisloc (*120,kid,z(ist1),2,nw1,iw)
   z(stce+ii) = z(stce+ii) - andf(z(stce+ii),z(ist1+iw))
   120 IF (ix(ps2,1) == 0) CYCLE
   kid = z(stce+ii+1+ps2)
   CALL bisloc (*130,kid,z(ist2),2,nw2,iw)
   z(stce+ii) = z(stce+ii) - andf(z(stce+ii),z(ist2+iw))
 END DO
 CALL OPEN (*150,scconn,z(buf1),1)
 DO  i = 1,nce,nwrd
   ii = i - 1
   IF (z(stce+ii) /= 0) CALL WRITE (scconn,z(stce+ii),nwrd,1)
 END DO
 CALL eof (scconn)
 CALL CLOSE (scconn,1)
 RETURN
 
 150 imsg = -1
 GO TO 190
 160 imsg = -2
 GO TO 190
 170 imsg = -3
 GO TO 190
 180 imsg = -8
 190 CALL mesage (imsg,ifile,aaa)
 RETURN
END SUBROUTINE cmrels
