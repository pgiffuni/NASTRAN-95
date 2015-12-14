SUBROUTINE cmsfil
     
    !     THIS SUBROUTINE GENERATES THE WORKING SUBSTRUCTURE FILE AND
    !     APPLIES ALL TRANSFORMATIONS
 
    EXTERNAL        rshift,andf
    LOGICAL :: tdat,xtran,xcstm,iauto
    !WKBI ALPHA-OSF 9/94
    INTEGER :: scr3, scsfil
    INTEGER :: aaa(2),andf,rshift,bufex,outt,  &
        scbdat,buf4,score,trn,sym,combo,nam(2),z,scmcon,  &
        sccstm,cstmid,cgid,buf2,buf3,buf1,ecpt1,twojm1
    DIMENSION       rz(1),ecpt(4),dofn(6),list(32),tg(3,3),tg6(6,6),  &
        tsave(6,6),tmat(6,6),tc(3,3),tt(3,3),xx(3)
    COMMON /cmb001/ scr1,scr2,scbdat,scsfil,scconn,scmcon,sctoc,  &
        geom4,casecc,sccstm,scr3
    COMMON /cmb002/ buf1,buf2,buf3,buf4,buf5,score,lcore,intp,outt
    COMMON /cmb003/ combo(7,5),conset,iauto,toler,npsub,conect,tran,  &
        mcon,restct(7,7),isort,origin(7,3),iprint
    COMMON /cmb004/ tdat(6)
    COMMON /gtmatx/ loc1,kltran,trn,tt6(6,6),tc6(6,6)
    COMMON /zzzzzz/ z(1)
    EQUIVALENCE     (ecpt1,ecpt(1)), (rz(1),z(1))
    DATA    aaa   / 4H cms, 4HFIL  /
    DATA    nhbgss, nhcstm,nheqss  / 4HBGSS,4HCSTM,4HEQSS /
 
    bufex  = lcore - buf1 + buf2
    lcore  = bufex - 1
    IF (lcore < 0) GO TO 620
    llco   = lcore
    ioefil = 310
    ifile  = scbdat
    CALL OPEN (*600,scbdat,z(buf4),0)
 
    !     READ GTRAN DATA INTO OPEN CORE
 
    IF (tdat(3) .OR. tdat(6)) CALL skpfil (scbdat,1)
    IF (.NOT.tdat(6)) GO TO 20
    ksgtrn = score
    CALL READ (*610,*10,scbdat,z(ksgtrn),llco,1,klgtrn)
    GO TO 620
10  llco   = llco - klgtrn
    kfgtrn = ksgtrn + klgtrn - 1
    GO TO 30
 
    !     READ TRANS DATA INTO OPEN CORE
 
20  kstran = score
    GO TO 40
30  kstran = kfgtrn + 1
40  IF (tdat(3) .OR. tdat(6)) CALL skpfil (scbdat,1)
    IF (.NOT. tdat(3)) GO TO 60
    CALL READ (*610,*50,scbdat,z(kstran),llco,1,kltran)
    GO TO 620
50  loc1   = kstran
    llco   = llco - kltran
    xtran  = .true.
    kftran = kstran + kltran - 1
    GO TO 70
60  loc1   = 0
    xtran  = .false.
70  CALL CLOSE (scbdat,1)
    ifile  = scsfil
    CALL OPEN (*600,scsfil,z(buf4),1)
    ifile  = sccstm
    CALL OPEN (*600,sccstm,z(buf3),1)
    kkc    = 0
 
    !     LOOP ON EACH PSEUDOSTRUCTURE
 
    ifile  = scr3
    CALL OPEN (*600,scr3,z(buf1),1)
    ifile  = ioefil
    CALL OPEN (*600,ioefil,z(bufex),1)
    llcold = llco
 
    DO  i = 1,npsub
        llco   = llcold
        nam(1) = combo(i,1)
        nam(2) = combo(i,2)
        trn    = combo(i,3)
        sym    = combo(i,4)
        ncomp  = combo(i,5)
   
        !     READ BGSS FOR I-TH PSEUDOSTRUCTURE
   
        ksbgss = kstran
        IF (xtran) ksbgss = kstran + kltran
        CALL sfetch (nam,nhbgss,1,itest)
        ngrp   = 1
        CALL sjump (ngrp)
        CALL suread (z(ksbgss),llco,klbgss,itest)
        IF (klbgss == llco .AND. itest /= 2) GO TO 620
        llco   = llco - klbgss
        kfbgss = ksbgss + klbgss - 1
        nip    = klbgss/4
   
        !     READ CSTM FOR THIS PSEUDOSTRUCTURE
   
        CALL sfetch (nam,nhcstm,1,itest)
        xcstm  = .false.
        loc2   = 0
        IF (itest == 3) GO TO 80
        kscstm = kfbgss + 1
        ngrp   = 1
        CALL sjump (ngrp)
        CALL suread (z(kscstm),llco,klcstm,itest)
        IF (klcstm == llco .AND. itest /= 2) GO TO 620
        llco   = llco - klcstm
        loc2   = kscstm
        xcstm  = .true.
        kfcstm = kscstm + klcstm - 1
80  CONTINUE
   
    !     DEFINE OPEN CORE ARRAYS FOR NEW CID AND HPTR
   
    ksncid = kfbgss + 1
    IF (loc2 /= 0) ksncid = kfcstm + 1
    klncid = nip
    llco   = llco - klncid
    kfncid = ksncid + klncid - 1
   
    kshptr = kfncid + 1
    klhptr = nip
    llco   = llco - klhptr
    kfhptr = kshptr + klhptr - 1
   
    !     SET ARRAYS TO ZERO
   
    DO  j = ksncid,kfncid
        z(j) = 0
    END DO
    DO  j = kshptr,kfhptr
        z(j) = 0
    END DO
   
    !     GET THE TRANS AND SYMT MATRIX FOR THIS PSEUDOSTRUCTURE
   
    CALL gtmat1 (sym,tt)
   
    !     TRANSFORM THE COORDINATES IN THE BGSS, NOTE THAT THE ORIGINS
    !     FOR TRANSLATION ARE STORED IN ARRAY ORIGIN.
   
    IF (trn+sym == 0) GO TO 130
    DO  j = ksbgss,kfbgss,4
        IF (z(j) == -1) CYCLE
        CALL gmmats (tt,3,3,0 ,rz(j+1),3,1,0 ,xx)
        DO  jj = 1,3
            rz(j+jj) = xx(jj) + origin(i,jj)
        END DO
    END DO
130 CONTINUE
   
    !     TRANSFORM DEGREES OF FREEDOM FOR EACH EQSS CONTAINED
    !     IN THE PSEUDOSTRUCTURE.
   
    CALL WRITE (scr3,tt6,36,1)
    nhmat = 1
    CALL sfetch (nam,nheqss,1,itest)
    kneqss = kfhptr + 1
    CALL suread (z(kneqss),4,kleqss,itest)
    CALL suread (z(kneqss),-1,kleqss,itest)
    llco  = llco - 2*ncomp
    ifile = scmcon
    CALL OPEN (*600,scmcon,z(buf2),1)
    DO  j = 1,ncomp
        kseqss = kfhptr + 2*ncomp
        CALL suread (z(kseqss),llco,kleqss,itest)
        IF (kleqss == 0) GO TO 340
        IF (kleqss == llco .AND. itest /= 2) GO TO 620
        kfeqss = kseqss + kleqss - 1
     
        !     LOOP ON EACH IP IN THE EQSS AND GENERATE TRANSFORMATION MATRIX
     
        DO  jj = kseqss,kfeqss,3
            ip    = z(jj+1)
            icomp = z(jj+2)
       
            !     GET CSTM FOR THIS IP
       
            cstmid = z(ksbgss+4*ip-4)
            ecpt1  = cstmid
            IF (cstmid < 0) ecpt1 = 0
            DO  jdh = 1,3
                ecpt(jdh+1) = rz(ksbgss+4*ip-4+jdh)
            END DO
            CALL gtmat2 (loc2,klcstm,ecpt,tc)
       
            !     TEST FOR POSSIBLE GTRAN
       
            igtran = 0
            IF (.NOT.tdat(6)) GO TO 170
            cgid = 1000000*j + z(jj)
            DO  k = ksgtrn,kfgtrn,5
                IF (z(k+3) == cgid .AND. z(k) == i .AND. z(k+1) == j) GO TO 160
            END DO
       
            !     NO GTRAN
       
            GO TO 170
160         CALL gtmat3 (z(k+4),tg,tg6,ikind)
            igtran = 1
            GO TO 180
170         CALL gtmat3 (-1,tg,tg6,ikind)
       
            !     ALL TRANSFORMATIONS HAVE BEEN FOUND, COMPUTE THE FINAL MATRIX TMAT
       
180         CALL gmmats (tg6  ,6,6,1,tt6,6,6,0,tsave)
            CALL gmmats (tsave,6,6,0,tc6,6,6,0,tmat )
       
            !     DECODE DEGREES OF FREEDOM AND FORM VECTOR
       
            CALL decode (icomp,list,ndof)
       
            !     FIND NEW DEGREES OF FREEDOM AND UPDATE EQSS
       
            IF (cstmid /= 0 .AND. igtran == 0) GO TO 220
            loop200:  DO  i1 = 1,6
                dofn(i1) = 0.0
                DO  i2 = 1,ndof
                    l = list(i2) + 1
                    IF (ABS(tmat(l,i1)) < 1.0E-4) CYCLE
                    dofn(i1) = 1.0
                    CYCLE loop200
                END DO
            END DO loop200
            icode = 0
            DO  i1 = 1,6
                icode = icode + dofn(i1)*2**(i1-1)
            END DO
            GO TO 230
220         icode   = icomp
230         z(jj+2) = icode
       
            !     WRITE IP,C ON SCRATCH TO COMPUTE NEW SIL,C
       
            CALL WRITE (scmcon,ip,1,0)
            CALL WRITE (scmcon,icode,1,0)
       
            !     UPDATE CID NUMBERS
       
            iadd   = ksbgss + 4*ip - 4
            iadd1  = ksncid + ip - 1
            ikkind = ikind  + 1
            GO TO (240,240,250,250,240,240,260,260,290,290,290,290,270,270,  &
                280,280,290,290,290,290,290,290,290,290,290,290,290,290,  &
                280,280,280,280,280,280,280), ikkind
240         z(iadd1) = z(iadd)
            IF (z(iadd) == -1) z(iadd1) = -100000000
            IF (z(iadd) == -2) z(iadd1) = -200000000
            GO TO 290
       
            !     COMMENTS FROM G.CHAN/UNISYS  9/92
            !     250 AND 240 ARE IDENTICAL HERE. IS IT POSSIBLY AN ERROR HERE?
       
250         z(iadd1) = z(iadd)
            IF (z(iadd) == -1) z(iadd1) = -100000000
            IF (z(iadd) == -2) z(iadd1) = -200000000
            GO TO 290
260         z(iadd1) = 0
            GO TO 290
270         z(iadd1) = -trn
            GO TO 290
280         z(iadd1) = -z(k+4)
290     CONTINUE
       
        !     SET POINTERS FOR H MATRIX
       
        itis = 0
        iadd2 = kshptr + ip - 1
        IF (cstmid   < 0) GO TO 300
        IF (z(iadd2) > 2) CYCLE
        IF (ikkind == 3 .OR. ikkind == 4 .OR. ikkind == 13 .OR.  &
            ikkind == 14) itis = 1
        IF (ikkind == 1 .OR. ikkind == 2 .OR. ikkind == 5 .OR.  &
            ikkind == 6) itis = 2
        IF (itis - 1 < 0) THEN
            GO TO   320
        ELSE IF (itis - 1 == 0) THEN
            GO TO   300
        ELSE
            GO TO   310
        END IF
300     z(iadd2) = 0
        CYCLE
310     z(iadd2) = 1
        CYCLE
320 CONTINUE
    nhmat = nhmat + 1
    z(iadd2) = nhmat
    CALL WRITE (scr3,tmat,36,1)
END DO
     
!     INSERT MULTIPLE IP CODE
     
IF (ncomp /= 1) CALL eqscod (kseqss,kleqss,z(1))
     
!     WRITE EQSS ON FILE SCSFIL
     
340 CALL WRITE (scsfil,z(kseqss),kleqss,1)
twojm1 = 2*(j-1)
IF (andf(rshift(iprint,19),1) == 1)  &
    CALL cmiwrt (1,nam,z(kneqss+twojm1),kseqss,kleqss,z,z)
END DO
CALL eof (scr3)
CALL CLOSE (scmcon,1)
   
!     GENERATE NEW SIL,C LIST
   
ifile = scmcon
CALL OPEN (*600,scmcon,z(buf2),0)
CALL READ (*360,*360,scmcon,z(kseqss),llco,1,nnn)
GO TO 620
360 CALL sort (0,0,2,1,z(kseqss),nnn)
ksej = kseqss + nnn - 1
i1 = kseqss
i2 = kseqss + 2
370 IF (i2-kseqss >= nnn) GO TO 400
IF (z(i1)   == z(i2)) GO TO 380
i1 = i1 + 2
i2 = i2 + 2
GO TO 370
380 DO  j = i2,ksej
    z(j-2) = z(j)
END DO
ksej = ksej - 2
nnn  = nnn  - 2
GO TO 370
400 CONTINUE
    z(kseqss) = 1
    DO  j = 3,nnn,2
        jj = j - 1
        icode = z(kseqss+jj-1)
        CALL decode (icode,list,ndof)
        z(kseqss+jj) = z(kseqss+jj-2) + ndof
    END DO
    CALL WRITE  (scsfil,z(kseqss),nnn,1)
    CALL suread (z(kseqss),llco,kleqss,itest)
    CALL WRITE  (ioefil,z(kseqss),kleqss,1)
    CALL CLOSE  (scmcon,1)
   
    !     PRINT EQSS SIL LIST IF REQUESTED
   
    IF (andf(rshift(iprint,19),1) == 1)  &
        CALL cmiwrt (8,nam,0,kseqss,kleqss,z,z)
   
    !     UPDATE CSTM NUMBERING SYSTEM
    !     KKC IS TRANSFORMED SYSTEM COORD. ID
   
    ip  = 0
    DO  i6 = ksncid,kfncid
        ip  = ip + 1
        loc = ksbgss + 4*(ip-1)
        IF (z(i6) == 100000000) CYCLE
        IF (z(i6) < 0.0) THEN
            GO TO   430
        ELSE IF (z(i6) == 0.0) THEN
            GO TO   520
        END IF
420     i1  = kscstm
        i2  = kfcstm
        GO TO 440
430     IF (z(i6) == -100000000) GO TO 530
        i1  = kstran
        i2  = kftran
440     kkc = kkc + 1
        IF (iauto) GO TO 480
        IF (kkc > 1) GO TO 460
        CALL page1
        CALL page2 (5)
        WRITE  (outt,450)
450     FORMAT (//45X,'SUMMARY OF OVERALL SYSTEM COORDINATES', //36X,  &
            'PSEUDO STRUCTURE ID.   SYSTEM COORD.ID    USER COORD.ID',/)
460     CALL page2 (1)
        WRITE  (outt,470) i,kkc,z(loc)
470     FORMAT (43X,i6,14X,i6,11X,i6)
480     look4 = z(i6)
        DO   j6 = i1,i2,14
            IF (IABS(z(i6)) /= z(j6)) CYCLE
            IF (z(i6) <= 0) GO TO 490
            CALL gmmats (tt,3,3,0,z(j6+5),3,3,0,tc)
            CALL WRITE (sccstm,kkc,1,0)
            CALL WRITE (sccstm,z(j6+1),4,0)
            CALL WRITE (sccstm,tc,9,0)
            CYCLE
490     CONTINUE
        CALL WRITE (sccstm,kkc,1,0)
        CALL WRITE (sccstm,z(j6+1),13,0)
    END DO
       
    !     FIND OTHER CIDS THAT ARE THE SAME
       
    iip = 0
    DO  j6 = ksncid,kfncid
        iip = iip + 1
        IF (z(j6) /= look4) CYCLE
        loc = ksbgss + 4*(iip-1)
        z(loc) = kkc
        z(j6) = 100000000
    END DO
    CYCLE
520 z(loc) = 0
    CYCLE
530 z(loc) = -1
    IF (z(i6) == -200000000) z(loc) = -2
END DO
     
!     WRITE PROCESSED BGSS
     
CALL WRITE (scsfil,z(ksbgss),klbgss,1)
IF (andf(rshift(iprint,18),1) ==  1)  &
    CALL cmiwrt (2,nam,nam,ksbgss,klbgss,z(1),z(1) )
     
!     WRITE ARRAY OF H POINTERS
     
CALL WRITE (scsfil,z(kshptr),klhptr,1)
CALL eof (scsfil)
END DO
   
CALL CLOSE (scr3,1)
CALL CLOSE (scsfil,1)
CALL WRITE (sccstm,tmat,0,1)
CALL CLOSE (sccstm,1)
CALL CLOSE (ioefil,1)
lcore = bufex + buf1 - buf2
RETURN
   
600 imsg = -1
GO TO 630
610 imsg = -2
GO TO 630
620 imsg = -8
630 CALL mesage (imsg,ifile,aaa)
RETURN
END SUBROUTINE cmsfil
