SUBROUTINE cmdisc
     
!     THIS SUBROUTINE DETERMINES THE DISCONNECTED DEGREES OF FREEDOM
!     AND GENERATES  DISCONNECTION ENTRIES  WHICH ARE MERGED WITH THE
!     CONNECTION ENTRIES
 
 EXTERNAL        orf
 INTEGER :: scsfil,z,score,combo,iptr(7),scconn,buf3,ce(9),  &
     orf,scdisc,de(9),aaa(2),scr1,scr2,buf2,outt
 COMMON /cmb001/ scr1,scr2,scbdat,scsfil,scconn,scmcon,sctoc, geom4,casecc
 COMMON /cmb002/ buf1,buf2,buf3,buf4,buf5,score,lcore,inpt,outt
 COMMON /cmb003/ combo(7,5),conset,iauto,toler,npsub
 COMMON /cmb004/ tdat(6),nipnew
 COMMON /zzzzzz/ z(1)
 DATA    aaa   / 4HCMDI,4HSC   /
 
 
 nwd = npsub+2
 isvcor = score
 itot = 0
 ilen = lcore
 nn = 0
 kk = score
 CALL OPEN (*200,scsfil,z(buf3),0)
 
!     LOOP ON THE NUMBER OF PSEUDO STRUCTURES READING THE SIL,C TABLE
!     INTO CORE FOR EACH.  THE ARRAY IPTR(I) POINTS TO THE START OF
!     THE I-TH TABLE IN CORE
 
 DO  i = 1,npsub
   ncsub = combo(i,5)
   
!     FIND SIL, C TABLE
   
   DO  j = 1,ncsub
     CALL fwdrec (*210,scsfil)
   END DO
   kk = kk + nn
   iptr(i) = kk
   CALL READ (*210,*20,scsfil,z(kk),lcore,1,nn)
   GO TO 220
   
!     ZERO OUT SIL VALUES, LOCATION WILL STORE CNEW
   
   20 DO  j = 1,nn,2
     z(kk+j-1) = 0
   END DO
   lcore = lcore - nn
   itot  = itot + nn
   CALL skpfil (scsfil,1)
 END DO
 CALL CLOSE (scsfil,1)
 
!     ALL EQSS HAVE BEEN PROCESSED, NOW SCAN THE CONNECTION ENTRIES
!     AND GET CNEW VALUES.
 
 CALL OPEN (*200,scconn,z(buf3),0)
 
!     READ AND PROCESS CONNECTION ENTRIES ONE AT A TIME
 
 50 CALL READ (*80,*60,scconn,ce,10,1,nn)
 60 DO  i = 1,npsub
   IF (ce(2+i) == 0) CYCLE
   
!     TRANSLATE CODED IP TO ACTUAL IP, COMPUTE LOCATION IN OPEN CORE
!     AND UPDATE CNEW
   
   ip = ce(2+i) - 1000000*(ce(2+i)/1000000)
   loc = iptr(i) + 2*ip - 2
   z(loc) = orf(z(loc),ce(1))
 END DO
 GO TO 50
 
!     ALL CONNECTIONS HAVE BEEN ACCOUNTED FOR,NOW DETERMINE DISCONN.
 
 80 CONTINUE
 scdisc = scr1
 IF (scr1 == scconn) scdisc = scr2
 CALL OPEN (*200,scdisc,z(buf2),1)
 DO  i = 1,npsub
   IF (i < npsub) LEN = iptr(i+1) - iptr(i)
   IF (i == npsub) LEN = itot - iptr(i)
   istrt = iptr(i)
   DO  j = 1,LEN,2
     DO  kdh = 1,9
       de(kdh) = 0
     END DO
     ip  = j/2 + 1
     loc = istrt + j - 1
     
!     POINT IS TOTALLY DISCONNECTED
     
     IF (z(loc) == z(loc+1)) CYCLE
     IF (z(loc) /= 0) GO TO 100
     
!     POINT IS TOTALLY CONNECTED
     
     de(1) = z(loc+1)
     de(2) = 2**i
     de(2+i) = ip
     GO TO 110
     
!     POINT IS PARTIALLY DISCONNECTED
     
     100 de(1) = z(loc+1) - z(loc)
     de(2) = 2**i
     de(2+i) = ip
     110 CALL WRITE (scdisc,de,nwd,1)
   END DO
 END DO
 CALL eof (scdisc)
 CALL CLOSE (scdisc,1)
 kk = score
 lcore = ilen
 CALL OPEN (*200,scdisc,z(buf2),0)
 CALL REWIND (scconn)
 id = 1
 140 CALL READ (*150,*160,scdisc,z(kk),lcore,1,nnn)
 GO TO 220
 150 id = 2
 CALL READ (*170,*160,scconn,z(kk),lcore,1,nnn)
 GO TO 220
 160 kk = kk + nwd
 lcore = lcore - nwd
 IF (lcore < nwd) GO TO 220
 IF (id == 1) GO TO 140
 GO TO 150
 170 CALL CLOSE (scconn,1)
 CALL CLOSE (scdisc,1)
 CALL OPEN (*200,scconn,z(buf3),1)
 LEN = kk - score
 nipnew = LEN/nwd
 DO  i = 1,LEN,nwd
   z(score+i) = IABS(z(score+i))
 END DO
 CALL sort (0,0,nwd,2,z(score),LEN)
 DO  i = 1,LEN,nwd
   CALL WRITE (scconn,z(score+i-1),nwd,1)
 END DO
 CALL eof (scconn)
 CALL CLOSE (scconn,1)
 CALL CLOSE (scdisc,1)
 score = isvcor
 lcore = ilen
 RETURN
 
 200 imsg = -1
 GO TO 230
 210 imsg = -2
 GO TO 230
 220 imsg = -8
 230 CALL mesage (imsg,ifile,aaa)
 RETURN
END SUBROUTINE cmdisc
