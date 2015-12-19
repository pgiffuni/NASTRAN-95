SUBROUTINE cmcomb (nps,nent,ndof,ic)
     
!     THIS SUBROUTINE COMBINES CONNECTION ENTRIES THAT HAVE BEEN SPECIFI
!     ON SEVERAL CONCT OR CONCT1 CARDS.
 
 
 INTEGER, INTENT(IN)                      :: nps
 INTEGER, INTENT(IN OUT)                  :: nent
 INTEGER, INTENT(IN)                      :: ndof
 INTEGER, INTENT(IN OUT)                  :: ic(nent,nps,ndof)
 EXTERNAL        orf
 LOGICAL :: match
 INTEGER :: ce(9),ceid,scconn,scmcon,buf1,buf2,savce,orf,z,  &
     scr2,buf3,score,comset,io(10),saconn,aaa(2)
 DIMENSION  list(32),krow(6),iertab(2000)
 COMMON /cmb001/ scr1,scr2,junk(2),scconn,scmcon
 COMMON /cmb002/ buf1,buf2,buf3,junk1(2),score,lcore
 COMMON /zzzzzz/ z(1)
 COMMON /BLANK / istep,idry
 DATA    aaa   / 4HCMCO,4HMB   /
 
!     CE IS THE CONNECTION ENTRY
!     KROW(I) IS THE NO. OF ROWS IN THE ITH DOF MATRIX
 
 iersub = 0
 itomny = 0
 ifile  = scconn
 CALL OPEN (*400,scconn,z(buf1),0)
 ifile = scmcon
 CALL OPEN (*400,scmcon,z(buf2),0)
 nrec  = -1
 npss  = nps - 1
 nword = nps + 1
 ient  = 0
 DO  i = 1,6
   krow(i) = 0
 END DO
 savce = 0
 20 CALL READ (*410,*190,scmcon,ceid,1,0,nnn)
 nrec  = ceid - savce - 1
 savce = ceid
 
!     GO FIND ENTRY NO. CEID
 
 ifile = scconn
 IF (nrec == 0) GO TO 40
 DO  i = 1,nrec
   CALL fwdrec (*420,scconn)
 END DO
 
!     READ IN CONNECTION ENTRY
 
 40 CALL READ (*410,*50,scconn,ce,10,1,nnn)
 
!     FIND WHICH DOF ARE PRESENT IN CONNECTION ENTRY
 
 50 CALL decode (ce(1),list,ncomp)
 loop180:  DO  i = 1,ncomp
   icomp = list(i) + 1
   IF (krow(icomp) == 0) GO TO 170
   
!     FIND FIRST NON-ZERO ENTRY IN CURRENT CE
   
   DO  j = 1,npss
     IF (ce(j+2) == 0) CYCLE
     isub = j
     GO TO 70
   END DO
   
!     NOW HAVE FOUND FIRST NON-ZERO, SEARCH FOR POSSIBLE
!     MATCHING ENTRIES IN MATRIX
   
   70 nloop = krow(icomp)
   DO  j = 1,nloop
     match  = .false.
     nersub = 0
     DO  jj = isub,npss
       IF (ic(j,jj,icomp) == 0 .OR. ce(jj+2) == 0) CYCLE
       IF (ic(j,jj,icomp)-ce(jj+2) == 0) THEN
         GO TO   100
       END IF
       80 IF (iersub+nersub > 2000) itomny = 1
       IF (iersub+nersub > 2000) GO TO 90
       iertab(iersub+nersub+1) = icomp
       iertab(iersub+nersub+2) = jj
       iertab(iersub+nersub+3) = ic(j,jj,icomp)
       iertab(iersub+nersub+4) = ce(jj+2)
       nersub = nersub + 4
       90 CONTINUE
       CYCLE
       100 match = .true.
     END DO
     IF (match) iersub = iersub + nersub
     IF (.NOT.match) CYCLE
     DO  jj = isub,npss
       IF (ce(jj+2) /= 0 .AND. ic(j,jj,icomp) /= 0) CYCLE
       ic(j,jj,icomp) = ic(j,jj,icomp) + ce(jj+2)
     END DO
     ic(j,npss+1,icomp) = orf(ic(j,npss+1,icomp),ce(2))
     CYCLE loop180
   END DO
   150 DO  jj = 1,npss
     ic(nloop+1,jj,icomp) = ce(jj+2)
   END DO
   ic(nloop+1,npss+1,icomp) = ce(2)
   krow(icomp) = krow(icomp) + 1
   CYCLE loop180
   170 nloop = 0
   GO TO 150
 END DO loop180
 GO TO 20
 190 CONTINUE
 IF (iersub == 0) GO TO 200
 
!     GENERATE ERROR TABLE AND TERMINATE
 
 CALL CLOSE (scconn,1)
 CALL CLOSE (scmcon,1)
 CALL cmtrce (iertab,iersub,itomny)
 idry = -2
 RETURN
 
 200 CONTINUE
 CALL CLOSE (scconn,1)
 ifile = scr2
 CALL OPEN (*400,scr2,z(buf3),1)
 DO  k = 1,ndof
   irow = krow(k)
   IF (irow > 0) THEN
     GO TO   210
   ELSE
     GO TO   240
   END IF
   210 DO  i = 1,irow
     io(1) = k
     io(2) = ic(i,nps,k)
     DO  j = 1,npss
       io(j+2) = ic(i,j,k)
     END DO
     CALL WRITE (scr2,io(1),nps+1,0)
   END DO
   240 CONTINUE
 END DO
 CALL WRITE (scr2,io(1),0,1)
 CALL CLOSE (scr2,1)
 CALL OPEN (*400,scr2,z(buf3),0)
 CALL READ (*410,*250,scr2,z(score),lcore,1,nwd)
 GO TO 430
 250 CALL sort (0,0,nps+1,2,z(score),nwd)
 CALL CLOSE (scr2,1)
 CALL OPEN (*400,scr2,z(buf3),1)
 ifin = score + nwd - 1
 iinc = nps + 1
 DO  i = score,ifin,iinc
   IF (z(i) == 0.0) THEN
     GO TO   310
   END IF
   260 comset = z(i)
   ibeg = i + iinc
   loop280:  DO  j = ibeg,ifin,iinc
     IF (z(j) == 0) CYCLE loop280
     IF (z(j+1) > z(i+1)) EXIT loop280
     DO  k = 1,npss
       IF (z(i+k+1) /= z(j+k+1)) CYCLE loop280
     END DO
     comset = 10*comset+z(j)
     z(j) = 0
   END DO loop280
   290 CALL encode (comset)
   io(1) = comset
   DO  kk = 1,nps
     io(1+kk) = z(i+kk)
   END DO
   CALL WRITE (scr2,io,nps+1,1)
   310 CONTINUE
 END DO
 CALL REWIND (scmcon)
 ifile = scmcon
 CALL READ (*410,*320,scmcon,z(score),lcore,1,nmcon)
 320 nce = 0
 saconn = scconn
 CALL OPEN (*400,scconn,z(buf1),0)
 330 CALL READ (*360,*340,scconn,ce,10,1,nnn)
 340 nce = nce + 1
 DO  i = 1,nmcon
   IF (nce == z(score+i-1)) GO TO 330
 END DO
 CALL WRITE (scr2,ce,nps+1,1)
 GO TO 330
 360 CALL CLOSE (scmcon,1)
 CALL CLOSE (scconn,1)
 CALL CLOSE (scr2,1)
 scconn = scr2
 scr2   = saconn
 RETURN
 
 400 imsg = -1
 GO TO 440
 410 imsg = -2
 GO TO 440
 420 imsg = -3
 GO TO 440
 430 imsg = -8
 440 CALL mesage (imsg,ifile,aaa)
 RETURN
END SUBROUTINE cmcomb
