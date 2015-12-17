SUBROUTINE cmauto
     
!     THIS SUBROUTINE PROCESSES THE AUTOMATIC CONNECTION OF
!     SUBSTRUCTURES IN THE COMB1 MODULE
 
 EXTERNAL        rshift,andf
 LOGICAL :: PRINT,found,tdat,back,iauto
 INTEGER :: scsfil,scconn,buf1,buf2,snext(8),st,nwd(8),score,  &
     z,spk,snk,ce(9),svkk,andf,aaa(2),ssil(8),nsil(8),  &
     sts,combo,restct,outt,NAME(14),rshift,ihd(12), ibits(2),jbits(2)
 DIMENSION       rz(1),a(3),b(3)
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /cmb001/ scr1,scr2,scbdat,scsfil,scconn,scmcon,sctoc, geom4,casecc
 COMMON /cmb002/ buf1,buf2,buf3,buf4,buf5,score,lcore,inpt,outt
 COMMON /cmb003/ combo(7,5),conset,iauto,toler,npsub,conect,tran,  &
     mcon,restct(7,7),isort,origin(7,3),iprint
 COMMON /cmb004/ tdat(6)
 COMMON /zzzzzz/ z(1)
 COMMON /BLANK / step,idry
 COMMON /output/ ititl(96),ihead(96)
 COMMON /system/ xxx,iot,junk(6),nlpp,junk1(2),line,junk2(2),  &
     idat(3),junk7(7),isw
 EQUIVALENCE     (z(1),rz(1))
 DATA    aaa   / 4HCMAU, 2HTO /, iblnk / 4H    /
 DATA    ihd   / 4H sum, 4HMARY, 4H of , 4H aut, 4HOMAT, 4HICAL,  &
     4HLY g, 4HENER, 4HATED, 4H con, 4HNECT, 4HIONS/
 
 nlin  = 1000
 found = .false.
 PRINT = .false.
 IF (andf(rshift(iprint,10),1) == 1) PRINT = .true.
 np2 = 2*npsub
 DO  i = 1,np2,2
   j = i/2 + 1
   NAME(i  ) = combo(j,1)
   NAME(i+1) = combo(j,2)
 END DO
 DO  i = 1,96
   ihead(i) = iblnk
 END DO
 j = 1
 DO  i = 75,86
   ihead(i) = ihd(j)
   j = j + 1
 END DO
 isavs = score
 isavl = lcore
 ifile = scconn
 CALL OPEN (*310,scconn,z(buf2),3)
 IF (iauto) GO TO 40
 CALL CLOSE (scconn,1)
 RETURN
 
 40 ifile = scsfil
 CALL OPEN (*310,scsfil,z(buf1),0)
 ssil(1) = score
 nout = npsub + 2
 idir = isort + 1
 DO  i = 1,npsub
   sts = ssil(i)
   ncsub = combo(i,5)
   DO  j = 1,ncsub
     CALL fwdrec (*320,scsfil)
   END DO
   
!     READ SIL,C FOR THE I-TH PSEUDOSTRUCTURE
   
   CALL READ (*320,*60,scsfil,z(sts),lcore,1,nsil(i))
   GO TO 330
   60 lcore = lcore - nsil(i)
   snext(i) = score + nsil(i)
   score = score + nsil(i)
   st = snext(i)
   
!     READ BGSS FOR THE I-TH PSEUDOSTRUCTURE
   
   CALL READ (*320,*70,scsfil,z(st),lcore,1,nwd(i))
   GO TO 330
   70 snext(i+1) = snext(i) + nwd(i)
   ssil(i+1)  = snext(i) + nwd(i)
   score = score + nwd(i)
   CALL skpfil (scsfil,1)
   lcore = lcore - nwd(i)
   ni = nwd(i) + st
   
!     WRITE THE IP NUMBER OVER THE CID IN THE BGSS
!     WILL BE USED AFTER SORTING
   
   DO  j = st,ni,4
     jj = (j-st+4)/4
     IF (z(j)+1 == 0.0) THEN
       GO TO    90
     END IF
     80 z(j) = jj
     CYCLE
     90 z(j) = -jj
   END DO
 END DO
 
!     SORT EACH BGSS IN THE SPECIFIED COORDINATE DIRECTION
 
 DO  i = 1,npsub
   st = snext(i)
   CALL sortf (0,0,4,idir,rz(st),nwd(i))
 END DO
 i    = 1
 130 k    = 0
 kk   = 0
 back = .false.
 svkk = 0
 ic1  = ssil(i)
 nipi = nwd(i)/4
 j    = i + 1
 IF (restct(i,j) /= 1) GO TO 280
 140 ic2  = ssil(j)
 nipj = nwd(j)/4
 150 spk  = snext(i) + k + 1
 IF (z(spk-1) < 0) GO TO 260
 a(1) = rz(spk  )
 a(2) = rz(spk+1)
 a(3) = rz(spk+2)
 160 snk  = snext(j) + kk + 1
 IF (z(snk-1) < 0) GO TO 270
 b(1) = rz(snk  )
 b(2) = rz(snk+1)
 b(3) = rz(snk+2)
 IF (a(isort) < b(isort)-toler) GO TO 250
 IF (b(isort) < a(isort)-toler) GO TO 270
 IF (back) GO TO 170
 back = .true.
 svkk = kk
 170 CONTINUE
 asej = a(isort)
 bsej = b(isort)
 xsej = asej - bsej
 DO  mm = 1,3
   IF (mm == isort) CYCLE
   asej = a(mm)
   bsej = b(mm)
   xsej = a(mm) - b(mm)
   IF (ABS(xsej) > toler) GO TO 270
 END DO
 
!     GENERATE THE NEW CONNECTION ENTRY
 
 DO  kdh = 1,9
   ce(kdh) = 0
 END DO
 ce(2)   = 2**(i-1) + 2**(j-1)
 ce(2+i) = IABS(z(spk-1))
 ce(2+j) = IABS(z(snk-1))
 m1 = IABS(z(spk-1))
 m2 = IABS(z(snk-1))
 ce(1) = andf(z(ic1+2*m1-1),z(ic2+2*m2-1))
 found = .true.
 
!     WRITE THE CONNECTION ENTRY ON SCCONN
 
 IF (ce(1) /= 0) CALL WRITE (scconn,ce,nout,1)
 IF (  .NOT.PRINT) GO TO 240
 IF (ce(1) == 0) GO TO 240
 IF (nlin < nlpp) GO TO 220
 200 nlin = 0
 CALL page
 WRITE  (outt,210) (NAME(kdh),kdh=1,np2)
 210 FORMAT (/14X,22HCONNECTED   connection,29X,22HPSEUDOSTRUCTURE  nam&
     &es, /17X,3HDOF,9X,4HCODE,3X,7(3X,2A4)//)
 nlin = nlin + 10
 220 CALL bitpat (ce(1),ibits)
 CALL bitpat (ce(2),jbits)
 nlin = nlin + 1
 IF (nlin > nlpp) GO TO 200
 WRITE (outt,230) ibits(1),ibits(2),jbits(1),jbits(2), (ce(kdh+2),kdh=1,npsub)
 230 FORMAT (16X,a4,a2,5X,a4,a3,2X,7(3X,i8))
 240 CONTINUE
 GO TO 270
 250 kk   = svkk
 back = .false.
 260 k    = k + 4
 IF (k/4 < nipi) GO TO 150
 k    = 0
 kk   = 0
 svkk = 0
 back = .false.
 GO TO 280
 270 kk = kk + 4
 IF (kk/4 < nipj) GO TO 160
 GO TO 250
 280 j = j + 1
 IF (j <= npsub) GO TO 140
 i = i + 1
 j = i
 IF (i < npsub) GO TO 130
 WRITE  (outt,290)
 290 FORMAT (//40X,'NOTE - GRID POINTS IN PSEUDOSTRUCTURE INTERNAL',  &
     ' GRID NUMBERS')
 CALL CLOSE (scconn,1)
 CALL CLOSE (scsfil,1)
 score = isavs
 lcore = isavl
 IF (found .OR. tdat(1).OR.tdat(2)) RETURN
 
 WRITE  (outt,300) ufm
 300 FORMAT (a23,' 6531, NO CONNECTIONS HAVE BEEN FOUND DURING ',  &
     'AUTOMATIC CONNECTION PROCEDURE.')
 idry = -2
 RETURN
 
 310 imsg = -1
 GO TO 350
 320 imsg = -2
 GO TO 350
 330 imsg = -8
 350 CALL mesage (imsg,ifile,aaa)
 RETURN
END SUBROUTINE cmauto
