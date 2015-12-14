SUBROUTINE cmsofo
     
    !     THIS ROUTINE GENERATES THE NEW SOF DATA FOR A COMBINATION.
 
    EXTERNAL        rshift,andf
    LOGICAL :: tdat,tocopn,lonly
    INTEGER :: pora,papp,lods,loap,  &
        buf1,buf3,ce(9),cnam,combo,scsfil,scr1,buf2,score,  &
        scconn,sbgss,sctoc,z,aaa(2),getip,ent(5),  &
        namold(14),eog,scbdat,sccstm,outt,rshift,andf
    DIMENSION       sav1(3),sav2(9),rz(1), tmat(9),ecpt(4),rent(3),list(32)
    CHARACTER (LEN=23) :: ufm
    COMMON /xmssg / ufm
    COMMON /cmb001/ scr1,scr2,scbdat,scsfil,scconn,scmcon,sctoc,  &
        geom4,casecc,sccstm
    COMMON /cmb002/ buf1,buf2,buf3,buf4,buf5,score,lcore,inpt,outt
    COMMON /cmb003/ combo(7,5),conset,iauto,toler,npsub,conect,tran,  &
        mcon,restct(7,7),isort,origin(7,3),iprint,tocopn
    COMMON /cmb004/ tdat(6),nipnew,cnam(2),lonly
    COMMON /zzzzzz/ z(1)
    COMMON /BLANK / step,dry,pora
    EQUIVALENCE     (rz(1),z(1)),(itran,ecpt(1))
    DATA    aaa   / 4HCMSO,4HFO   / , eog / 4H$eog /
    DATA    papp  , lods,loap             / 4HPAPP,4HLODS,4HLOAP  /
    DATA    nheqss, nhbgss,nhcstm,nhplts  / 4HEQSS, 4HBGSS,4HCSTM,4HPLTS  /
 
    !     GET NAMES OF BASIC COMPONENTS FROM THE TABLE OF CONTENTS
 
    ifile = sctoc
    IF (.NOT.tocopn) CALL OPEN (*720,sctoc,z(buf2),0)
    CALL REWIND (sctoc)
    k = 0
    j = 0
10  j = j + 1
 
    CALL READ (*30,*730,sctoc,0,-3,0,nnn)
    CALL READ (*30,*20,sctoc,z(score+k),lcore,1,nnn)
20  k = k + nnn
    IF (j < npsub) GO TO 10
30  IF (.NOT.tocopn) CALL CLOSE (sctoc,1)
    istnm = score
    score = score + k
    nnames = k
 
    nsub = 0
    DO  j = 1,npsub
        nsub = nsub + combo(j,5)
    END DO
 
    !     WRITE THE FIRST GROUP OF THE EQSS
 
    IF (lonly) GO TO 330
 
    np2 = 2*npsub
    DO  i = 1,np2,2
        j = i/2 + 1
        namold(i  ) = combo(j,1)
        namold(i+1) = combo(j,2)
    END DO
    npp = npsub
    CALL setlvl (cnam,npp,namold,itest,29)
    IF (itest == 8) GO TO 700
    itest = 3
    CALL sfetch (cnam,nheqss,2,itest)
    itest = 1
    CALL suwrt (cnam,2,itest)
    CALL suwrt (nsub,1,itest)
    CALL suwrt (nipnew,1,itest)
    itest = 2
    CALL suwrt (z(istnm),nnames,itest)
 
    ifile = scconn
    CALL OPEN (*720,scconn,z(buf1),0)
    ifile = scsfil
    CALL OPEN (*720,scsfil,z(buf2),0)
    ifile = scr1
    CALL OPEN (*720,scr1,z(buf3),1)
    DO  i = 1,npsub
        ipn = 0
        kk  = 0
60      ipn = ipn + 1
        CALL READ (*80,*70,scconn,ce,10,1,nnn)
70      IF (ce(i+2) == 0) GO TO 60
        z(score+kk  ) = ce(  1)
        z(score+kk+1) = ce(i+2)
        z(score+kk+2) = ipn
        kk = kk + 3
        GO TO 60
80      noipn = (kk)/3
        CALL sort (0,0,3,2,z(score),kk)
   
        !     READ BGSS FROM SUBFIL
   
        ifile = scsfil
        npsp1 = combo(i,5) + 1
        DO  j = 1,npsp1
            CALL fwdrec (*730,scsfil)
        END DO
        sbgss = score + kk
        lcore = lcore - kk
        CALL READ (*730,*100,scsfil,z(sbgss),lcore,1,lbgss)
        GO TO 740
100 CONTINUE
    DO  j = 1,kk,3
        jj     = j - 1
        getip  = z(score+jj+1)
        ent(1) = z(sbgss+4*getip-4  )
        ent(2) = z(sbgss+4*getip-4+1)
        ent(3) = z(sbgss+4*getip-4+2)
        ent(4) = z(sbgss+4*getip-4+3)
        ent(5) = z(score+jj+2)
        CALL WRITE (scr1,ent,5,0)
    END DO
    CALL WRITE (scr1,ent,0,1)
    IF (i == 1) CALL REWIND (scsfil)
    IF (i /= 1) CALL skpfil (scsfil,-1)
    IF (i /= 1) CALL skpfil (scsfil, 1)
    ncomp = combo(i,5)
    DO  j = 1,ncomp
        CALL READ (*720,*120,scsfil,z(sbgss),lcore,1,nnn)
        GO TO 740
120     IF (nnn == 0) GO TO  160
        DO  jj = 1,nnn,3
            kid = z(sbgss+jj)
            CALL bisloc (*150,kid,z(score+1),3,noipn,nwd)
130         IF (z(score+nwd) /= z(score+nwd-3)) GO TO 140
            IF (nwd <= 1) GO TO 140
            nwd = nwd - 3
            GO TO 130
140         ent(1) = z(sbgss+jj -1)
            ent(2) = z(score+nwd+1)
            ent(3) = z(score+nwd-1)
            CALL WRITE (scr1,ent,3,0)
            IF (z(score+nwd) /= z(score+nwd+3)) CYCLE
            IF (nwd+3 >= noipn*3) CYCLE
            nwd = nwd + 3
            GO TO 140
150     CONTINUE
        END DO
160     CALL WRITE (scr1,0,0,1)
    END DO
    CALL skpfil (scsfil,1)
    CALL REWIND (scconn)
END DO
CALL CLOSE (scsfil,1)
CALL CLOSE (scr1,1)
 
!     WRITE OUT EQSS ONTO SOF
 
ifile = scr1
CALL OPEN (*720,scr1,z(buf3),0)
DO  i = 1,npsub
    ncomp = combo(i,5)
    CALL fwdrec (*730,scr1)
    DO  j = 1,ncomp
        CALL READ (*730,*190,scr1,z(sbgss),lcore,1,nnn)
        GO TO 740
190     itest = 2
        CALL suwrt (z(sbgss),nnn,itest)
    END DO
END DO
 
!     WRITE OUT MASTER SIL,C LIST FOR NEW STRUCTURE
 
CALL REWIND (scconn)
kk   = 0
isil = 1
220 CALL READ (*240,*230,scconn,ce,10,1,nnn)
230 CALL decode (ce(1),z(buf2),ndof)
z(sbgss+kk  ) = isil
z(sbgss+kk+1) = ce(1)
isil = isil + ndof
kk   = kk + 2
GO TO 220
240 CONTINUE
    itest = 2
    CALL suwrt (z(sbgss),kk,itest)
    itest = 3
    CALL suwrt (ent,0,itest)
 
    !     WRITE BGSS ONTO SOF
 
    lcc = lcore
    kk  = 0
    CALL REWIND (scr1)
    DO  i = 1,npsub
        ncomp = combo(i,5)
        CALL READ (*730,*250,scr1,z(sbgss+kk),lcc,1,nw)
        GO TO 740
250     kk  = kk  + nw
        lcc = lcc - nw
        DO  j = 1,ncomp
            CALL fwdrec (*280,scr1)
        END DO
    END DO
280 CALL sort (0,0,5,5,z(sbgss),kk)
    itest = 3
    CALL sfetch (cnam,nhbgss,2,itest)
    itest = 1
    CALL suwrt (cnam,2,itest)
    itest = 2
    CALL suwrt (nipnew,1,itest)
    kkk  = 0
290 isub = sbgss + kkk
    ent(1) = z(isub  )
    ent(2) = z(isub+1)
    ent(3) = z(isub+2)
    ent(4) = z(isub+3)
    IF (kkk > 0 .AND. z(isub+4) == z(isub-1)) GO TO 300
    itest = 1
    CALL suwrt (ent,4,itest)
300 kkk = kkk + 5
    IF (kkk < kk) GO TO 290
    itest = 2
    CALL suwrt (ent,0,itest)
    itest = 3
    CALL suwrt (ent,0,itest)
    CALL CLOSE (scr1,1)
    CALL CLOSE (scconn,1)
 
    !     PROCESS CSTM ITEM
 
    CALL OPEN (*720,sccstm,z(buf3),0)
    CALL READ (*310,*310,sccstm,z(score),lcore,1,nnn)
    GO TO 740
310 IF (nnn == 0) GO TO 320
    itest = 3
    CALL sfetch (cnam,nhcstm,2,itest)
    itest = 2
    CALL suwrt (cnam,2,itest)
    itest = 2
    CALL suwrt (z(score),nnn,itest)
    itest = 3
    CALL suwrt (0,0,itest)
320 CALL CLOSE (sccstm,1)
 
    !     PROCESS LODS ITEM
 
330 nlv  = 0
    ncs  = nnames/2
    j    = 0
    litm = lods
    IF (pora == papp) litm = loap
    DO  i = 1,npsub
        namold(1) = combo(i,1)
        namold(2) = combo(i,2)
        CALL sfetch (namold,litm,1,itest)
        IF (itest == 3) GO TO 340
        CALL suread (ce,4,nout,itest)
        nlv = nlv + ce(3)
        jdh = 1
        CALL sjump (jdh)
        CALL suread (z(score+j),-2,nout,itest)
        j = j + nout
        CYCLE
340     z(score+j  ) = 0
        z(score+j+1) = eog
        j = j + 2
    END DO
    itest = 3
    CALL sfetch (cnam,litm,2,itest)
    itest = 1
    CALL suwrt (cnam,2,itest)
    CALL suwrt (nlv ,1,itest)
    CALL suwrt (ncs ,1,itest)
    itest = 2
    CALL suwrt (z(istnm),nnames,itest)
    itest = 3
    CALL suwrt (z(score),j,itest)
    IF (lonly) GO TO 580
 
    !     PROCESS PLTS ITEM
 
 
    !     FIND OLD PLTS TRANSFORMATIONS
 
    nout  = 0
    j     = 0
    nnsub = 0
    DO  i = 1,npsub
        namold(1) = combo(i,1)
        namold(2) = combo(i,2)
        CALL sfetch (namold,nhplts,1,itest)
        IF (itest == 3) CYCLE
        CALL suread (z(score+j),3,nout,itest)
        nnsub = nnsub + z(score+j+2)
        CALL suread (z(score+j),-1,nout,itest)
        j = j + nout
    END DO
    npwd  = j
    istrn = score + npwd
    llco  = lcore - npwd
    itest = 3
    CALL sfetch (cnam,nhplts,2,itest)
    itest = 1
    CALL suwrt (cnam ,2,itest)
    CALL suwrt (nnsub,1,itest)
    nt = 0
    IF (.NOT.tdat(3)) GO TO 390
    CALL OPEN (*720,scbdat,z(buf1),0)
    CALL skpfil (scbdat,2)
    CALL READ (*730,*380,scbdat,z(istrn),lcore,1,nt)
    GO TO 740
380 CALL pretrs (z(istrn),nt)
390 IF (.NOT.tocopn) CALL OPEN (*720,sctoc,z(buf2),1)
    CALL REWIND (sctoc)
    ist  = istrn + nt
    llco = llco  - nt
    j = 0
400 j = j + 1
    itran = combo(j,3)
    DO  i = 1,9
        tmat(i) = 0.0
    END DO
    CALL READ (*530,*420,sctoc,z(ist),llco,1,nnn)
    GO TO 740
420 IF (itran == 0) GO TO 440
    DO  i = 2,4
        ecpt(i) = 0.0
    END DO
    CALL transs (ecpt,tmat)
    GO TO 450
440 tmat(1) = 1.0
    tmat(5) = 1.0
    tmat(9) = 1.0
 
    !     DETERMINE SYMMETRY
 
450 IF (combo(j,4) == 0) GO TO 470
    CALL decode (combo(j,4),list,ndir)
    DO  i = 1,ndir
        idir = list(i) + 1
        idir = 4 - idir
        tmat(idir  ) = -tmat(idir  )
        tmat(idir+3) = -tmat(idir+3)
        tmat(idir+6) = -tmat(idir+6)
    END DO
    470 DO  i = 1,3
        rent(i) = origin(j,i)
    END DO
    nnn = nnn - 1
    DO  i = 3,nnn,2
   
        !     PROCESS OLD TRANSFORMATIONS
   
        DO  kdh = 1,npwd,14
            IF (z(ist+i) == z(score+kdh-1) .AND. z(ist+i+1) == z(score+kdh))  &
                GO TO 500
        END DO
        CYCLE
500     CALL gmmats (tmat,3,3,0, rz(score+kdh+1),3,1,0, sav1)
        DO   ii = 1,3
            sav1(ii) = sav1(ii)+rent(ii)
        END DO
        CALL gmmats (tmat,3,3,0, rz(score+kdh+4),3,3,0, sav2)
        itest = 1
        CALL suwrt (z(ist+i),2,itest)
        CALL suwrt (sav1(1) ,3,itest)
        CALL suwrt (sav2(1) ,9,itest)
        CYCLE
    END DO
    IF (j < npsub) GO TO 400
530 IF (.NOT.tocopn) CALL CLOSE (sctoc,1)
    itest = 2
    CALL suwrt (0,0,itest)
    itest = 3
    CALL suwrt (0,0,itest)
    CALL CLOSE (scbdat,1)
    CALL eqsout
 
    !     PROCESS OUTPUT REQUESTS
 
    IF (andf(rshift(iprint,12),1) /= 1) GO TO 550
 
    !     WRITE EQSS FOR NEW STRUCTURE
 
    CALL sfetch (cnam,nheqss,1,itest)
    CALL suread (z(score),4,nout,itest)
    CALL suread (z(score),-1,nout,itest)
    ist = score + nout
    DO  i = 1,nsub
        CALL suread (z(ist),-1,nout,itest)
        iadd = score + 2*(i-1)
        CALL cmiwrt (1,cnam,z(iadd),ist,nout,z,z)
    END DO
    CALL suread (z(ist),-1,nout,itest)
    CALL cmiwrt (8,cnam,0,ist,nout,z,z)
550 IF (andf(rshift(iprint,13),1) /= 1) GO TO 560
 
    !     WRITE BGSS FOR NEW STRUCTURE
 
    CALL sfetch (cnam,nhbgss,1,itest)
    ngrp = 1
    CALL sjump (ngrp)
    ist = score
    CALL suread (z(ist),-1,nout,itest)
    CALL cmiwrt (2,cnam,cnam,ist,nout,z,z)
560 IF (andf(rshift(iprint,14),1) /= 1) GO TO 570
 
    !     WRITE CSTM ITEM
 
    CALL sfetch (cnam,nhcstm,1,itest)
    IF (itest == 3) GO TO 570
    ngrp = 1
    CALL sjump (ngrp)
    ist = score
    CALL suread (z(ist),-1,nout,itest)
    CALL cmiwrt (3,cnam,cnam,ist,nout,z,z)
570 IF (andf(rshift(iprint,15),1) /= 1) GO TO 580
 
    !     WRITE PLTS ITEM
 
    CALL sfetch (cnam,nhplts,1,itest)
    ist = score
    CALL suread (z(ist), 3,nout,itest)
    CALL suread (z(ist),-1,nout,itest)
    CALL cmiwrt (4,cnam,cnam,ist,nout,z,z)
580 IF (andf(rshift(iprint,16),1) /= 1) GO TO 600
 
    !     WRITE LODS ITEM
 
    CALL sfetch (cnam,lods,1,itest)
    IF (itest == 3) GO TO 600
    CALL suread (z(score), 4,nout,itest)
    CALL suread (z(score),-1,nout,itest)
    ist   = score + nout
    itype = 5
    IF (litm == loap) itype = 7
    DO  i = 1,nsub
        iadd = score + 2*(i-1)
        CALL suread (z(ist),-1,nout,itest)
        CALL cmiwrt (itype,cnam,z(iadd),ist,nout,z,z)
        itype = 6
    END DO
600 CONTINUE
    RETURN
 
700 WRITE  (outt,710) ufm
710 FORMAT (a23,' 6518, ONE OF THE COMPONENT SUBSTRUCTURES HAS BEEN ',  &
        'USED IN A PREVIOUS COMBINE OR REDUCE.')
    imsg = -37
    GO TO 750
720 imsg = -1
    GO TO 750
730 imsg = -2
    GO TO 750
740 imsg = -8
750 CALL sofcls
    CALL mesage (imsg,ifile,aaa)

    RETURN
END SUBROUTINE cmsofo
