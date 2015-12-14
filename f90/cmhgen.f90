SUBROUTINE cmhgen
     
!     THIS SUBROUTINE GENERATES THE (H) TRANSFORMATION MATRICES FOR
!     COMPONENT SUBSTRUCTURES IN A COMBINE OPERATION AND WRITES THEM
!     ON THE SOF
 
 LOGICAL :: frsfil
 INTEGER :: bufex,scr1,mcb(7),ihead(2),nam(2),scr3,buf1,  &
     cnam,combo,z,ssil,lcore,score,scsfil,buf2,buf3,  &
     scbdat,scconn,listo(32),listn(32),aaa(2),buf4, ce(10)
 DIMENSION       t(6,6),tp(6,6),tpp(6,6),colout(6),tid(6,6), ttran(6,6)
 COMMON /cmb001/ scr1,scr2,scbdat,scsfil,scconn,scmcon,sctoc,  &
     geom4,casecc,sccstm,scr3
 COMMON /cmb002/ buf1,buf2,buf3,buf4,buf5,score,lcore,inpt,outt
 COMMON /cmb003/ combo(7,5),conset,iauto,toler,npsub
 COMMON /cmb004/ tdat(6),nipnew,cnam(2)
 COMMON /zzzzzz/ z(1)
 COMMON /packx / iin,iout,iiii,nnnn,incr
 DATA    zero  / 0.0 /,aaa/ 4HCMHG,4HEN   /,ihead/ 4HHORG,4H    /
 DATA    tid   / 1.,0.,0.,0.,0.,0., 0.,1.,0.,0.,0.,0.,  &
     0.,0.,1.,0.,0.,0., 0.,0.,0.,1.,0.,0., 0.,0.,0.,0.,1.,0., 0.,0.,0.,0.,0.,1. /
 DATA    nheqss/ 4HEQSS   /
 
!     READ SIL,C FROM SOF FOR COMBINED STRUCTURE
 
 incr  = 1
 bufex = lcore - buf2 + buf3
 lcore = bufex - 1
 IF (lcore < 0) GO TO 320
 ioefil = 310
 CALL OPEN (*300,ioefil,z(bufex),0)
 mcb(1) = scr1
 mcb(4) = 2
 mcb(5) = 1
 iin    = 1
 iout   = 1
 CALL sfetch (cnam,nheqss,1,itest)
 nsub = 0
 DO  i = 1,npsub
   nsub = nsub + combo(i,5)
 END DO
 CALL sjump (nsub+1)
 CALL suread (z(score),-1,nsilnw,itest)
 
!     LOOP ON NUMBER OF PSEUDO-STRUCTURES BEING COMBINED
 
 ssil  = score + nsilnw
 lcore = lcore - nsilnw
 ifile = scr3
 CALL OPEN (*300,scr3,z(buf1),0)
 ifile = scsfil
 CALL OPEN (*300,scsfil,z(buf2),0)
 
 DO  i = 1,npsub
   frsfil = .true.
   mcb(2) = 0
   mcb(6) = 0
   mcb(7) = 0
   
!     READ SIL,C FOR COMPONENT SUBSTRUCTURE
   
   ncs = combo(i,5) + 2
   DO  j = 1,ncs
     CALL fwdrec (*310,scsfil)
   END DO
   ifile = ioefil
   CALL READ (*310,*30,ioefil,z(ssil),lcore,1,nslold)
   GO TO 320
   30 ishptr = ssil + nslold
   ifile  = scsfil
   CALL READ (*310,*40,scsfil,z(ishptr),lcore,1,lhptr)
   GO TO 320
   40 CALL skpfil (scsfil,1)
   
!     COMPUTE NUMBER OF ROWS IN MATRIX
   
   icode = z(ssil+nslold-1)
   CALL decode (icode, listo, ncom)
   mcb(3) = z(ssil+nslold-2) + ncom - 1
   
!     READ CONNECTION ENTRIES
   
!     READ TRANSFORMATION MATRIX FOR PSEUDOSTRUCTURE
   
   ifile = scr3
   CALL READ (*310,*50,scr3,ttran,37,1,nnn)
   50 CONTINUE
   CALL skpfil (scr3,-1)
   IF (i /= 1) CALL skpfil (scr3,1)
   ifile = scconn
   CALL OPEN (*300,scconn,z(buf3),0)
   ifile = scr1
   CALL OPEN (*300,scr1,z(buf4),1)
   CALL WRITE (scr1,ihead,2,1)
   ipnew = 0
   60 CALL READ (*250,*70,scconn,ce,10,1,nnn)
   70 ipnew  = ipnew + 1
   locipn = score + 2*(ipnew-1) + 1
   IF (ce(i+2) == 0) GO TO 230
   ipold  = ce(i+2)
   locipo = ssil + 2*(ipold-1) + 1
   icode  = z(locipn)
   CALL decode (icode,listn,ncn)
   icode  = z(locipo)
   CALL decode (icode,listo,nco)
   
   iaddh = ishptr + ipold - 1
   idh   = z(iaddh)
   IF (idh-1 < 0) THEN
     GO TO    80
   ELSE IF (idh-1 == 0) THEN
     GO TO   100
   ELSE
     GO TO   120
   END IF
   
!     IDENTITY MATRIX
   
   80 CONTINUE
   DO  i1 = 1,6
     DO  i2 = 1,6
       t(i1,i2) = tid(i1,i2)
     END DO
   END DO
   GO TO 160
   
!     TRANS MATRIX
   
   100 CONTINUE
   DO  i1 = 1,6
     DO  i2 = 1,6
       t(i1,i2) = ttran(i1,i2)
     END DO
   END DO
   GO TO 160
   
!     MATRIX DUE TO GTRAN
   
   120 CONTINUE
   idhm1 = idh - 1
   DO  i1 = 1,idhm1
     CALL fwdrec (*310,scr3)
   END DO
   CALL READ (*310,*140,scr3,t,37,1,nnn)
   140 DO  i1 = 1,idh
     CALL bckrec (scr3)
   END DO
   160 CONTINUE
   
!     DELETE ROWS OF (T) FOR EACH COLD EQUAL TO ZERO
   
   DO  j1 = 1,nco
     ir = listo(j1) + 1
     DO  j2 = 1,6
       tp(j1,j2) = t(ir,j2)
     END DO
   END DO
   nrow = nco
   
!     DELETE COLUMNS OF (T) FOR EACH CNEW EQUAL TO ZERO
   
   DO  j1 = 1,ncn
     ic = listn(j1) + 1
     DO  j2 = 1,nrow
       tpp(j2,j1) = tp(j2,ic)
     END DO
   END DO
   ncol = ncn
   DO  i1 = 1,ncol
     DO  i2 = 1,nrow
       colout(i2) = tpp(i2,i1)
     END DO
     iiii = z(locipo-1)
     nnnn = iiii + nrow - 1
     CALL pack (colout,scr1,mcb)
   END DO
   GO TO 60
   230 iiii = 1
   nnnn = 1
   icode = z(locipn)
   CALL decode (icode,listn,ncn)
   DO  i1 = 1,ncn
     CALL pack (zero,scr1,mcb)
   END DO
   GO TO 60
   250 CONTINUE
   CALL CLOSE (scconn,1)
   CALL wrttrl (mcb)
   CALL CLOSE (scr1,1)
   nam(1) = combo(i,1)
   nam(2) = combo(i,2)
   CALL mtrxo (scr1,nam,ihead(1),z(buf4),itest)
   CALL skpfil (scr3,1)
 END DO
 
 CALL CLOSE (scsfil,1)
 CALL CLOSE (scr3,1)
 CALL CLOSE (ioefil,1)
 lcore = bufex + buf2 - buf3
 RETURN
 
 300 imsg = -1
 GO TO 330
 310 imsg = -2
 GO TO 330
 320 imsg = -8
 330 CALL mesage (imsg,ifile,aaa)
 RETURN
END SUBROUTINE cmhgen
