SUBROUTINE mpyado (zz      ,z      ,zd      )
     
!     MPYAD PERFORMS THE MATRIX OPERATION
!       (+/-)A    * B (+/-)C = D   OR
!       (+/-)A(T) * B (+/-)C = D
 
!     LAST REVISED  1/92 BY G.CHAN/UNISYS
!     . NEW METHOD 4T WAS ADDED WHICH IS FASTER THAN METHOD 2T UNDER
!       CERTAIN CONDITIONS.
!     . NEW ROUTINE FOR DIAGONAL, IDENTITY, AND ROW VECTOR MATRICES
!     . USER CAN REVERT TO ORIGINAL MPYAD ALL METHODS, BY DIAG 41
 
 
!     LEGEND:
 
!     +---+ + +  IS A MATRIX       +-            \    AN ELEMENT OF
!     |   | | |  BY COLUMNS        |  IS A        \   A, B, OR C IN
!     |   | | |  IN MULTIPLE       |  COLUMN       \     AND
!     |   | | |  PASSES            |             OR   /  AN ELEMENT
!     +---+ + +                    +-                /   OF D OUT
 
!     +-------+  IS A MATRIX
!     |       |  BY ROWS           => OR   INDICATES MATRICES C AND D
!     +-------+  IN MULTIPLE       <=      ARE USING SAME CORE SPACE
!     +-------+  PASSES
!     +-------+
 
!     UPPER CASE LETTER INDICATES UNPACKED MATRIX OR COLUMN
!     LOWER CASE LETTER INDICATES MATRIX OR COLUMN IN STRINGS FORM
 
 
!     METHOD 1NT AND 1T       METHOD 2NT               METHOD 2T
!                                        B                        +-
!         +----+ + +                    /                         |
!         |    | | |                   /                          |B
!     a   | B  | | |          + + +---+ +-                        |
!      \  |    | | |          | | |   | |            +----------+ +-
!       \ +----+ + +          | | | a | |C           |    a     |
!         +----+              | | |   | |            +----------+
!         |    |              + + +---+ +- +-        +----------+ \
!         | D  | <= C                      |         +----------+  \
!         |    |                        => |D                     \ C
!         |    |                           |                       \
!         +----+                           +-                       D
 
!      METHOD 3T     +-                  METHOD 4T     +-
!                    |                                 |
!                    |b                                |b(BANDED)
!                    |                                 |
!                    +-                                +-
!         + + +----+ +-    +-              +---------+ +-
!         | | |    | |     |               |         | |
!         | | | A  | |D  + |C              |    a    | |C
!         | | |    | |     |               |         | |
!         + + +----+ +-    +-              +---------+ +-  +-
!                        ADD ON            +---------+     |
!                        LAST              +---------+     |D(FULL)
!                        PASS                           => |
!                                                          +-
 
 
 INTEGER, INTENT(OUT)                     :: zz(6)
 REAL, INTENT(OUT)                        :: z(1)
 DOUBLE PRECISION, INTENT(OUT)            :: zd(1)
 LOGICAL :: last  ,null
 EXTERNAL        andf  ,orf   ,lshift
!WKBI 9/93
 INTEGER :: prntyp(4), namea(2), nameb(2), namec(2), NAMED(2), prca
 INTEGER :: p     ,q     ,r     ,t     ,op    ,opa   ,  &
     opb   ,opc   ,opbc  ,op2   ,one1  ,one2  ,p1    ,  &
     pp1   ,pp2   ,prc   ,prec  ,prec1 ,bcd   ,rcb   ,  &
     rcd   ,rc    ,rd    ,rdrew ,wrt   ,wrtrew,cls   ,  &
     clsrew,andf  ,orf   ,eol   ,eor   ,acol  ,acol1 ,  &
     acoln ,acore ,apoint,point ,bcol  ,buf1  ,buf2  ,  &
     buf3  ,buf4  ,bufi  ,blk   ,BLOCK ,row   ,rowa  ,  &
     arow  ,arow1 ,arown ,crow  ,drow  ,TYPE  ,typea ,  &
     typeb ,typec ,typed ,typebd,typd  ,typd1 ,FILE  ,  &
     filea ,fileb ,filec ,filed ,cfile ,dfile ,efile ,  &
     scrtch,signab,signc ,firstl,sysbuf,FORM  ,flag  , densc
 DOUBLE PRECISION :: ad(2) ,bd(2) ,dd(2) , xnd
 DIMENSION       b(4)  , mpy(3),bcd(2),zero(4)      ,xns(1),  &
     NAME(2)      ,blk(15)      ,method(6)
!NVXNB
 COMMON /logout/ lout
!NVXNE
 COMMON /machin/ mach  ,ihalf ,jhalf /mpyqt4/ qt(2) ,ll4   ,jmp(2)
 COMMON /mpyadx/ filea(7)     ,fileb(7)     ,filec(7)     ,  &
     filed(7)     ,nz    ,t     ,signab,signc ,prec1 , scrtch,time  &
     /system/ ksystm(152) /TYPE  / prc(2),nwds(4)      ,rc(4)  &
     /names / rd    ,rdrew ,wrt   ,wrtrew,clsrew,cls /zblpkx/ d(4)  ,drow  &
     /zntpkx/ a(4)  ,ip    ,eol   ,eor  &
     /packx / typed ,typd1 ,one1  ,pp1   ,incr1  &
     /unpakx/ typebd,one2  ,pp2   ,incr2
 COMMON /ntime / nitems,tmio  ,tmbpak,tmipak,tmpak ,tmupak,tmgstr,  &
     tmpstr,tmt(4),tml(4)  &
     /mpyadz/ rcb   ,rcd   ,ll    ,lll   ,jb    ,nbx   ,ndx   ,  &
     jmax1x,acol  ,acol1 ,acoln ,acore ,apoint,bcol  ,  &
     crow  ,firstl,na    ,nb    ,nd    ,nwda  ,nwdb  ,  &
     nwdd  ,prec  ,jmax  ,incra ,BLOCK(20) /zzzzzz/ xnd(8500)
 EQUIVALENCE     (ksystm( 1),sysbuf) , (ksystm( 2),mout ) ,  &
     (ksystm(58),ksys58) , (ksystm(40),nbpw )
!WKBI 10/93  &
 ,(ksystm(55),iprec )
 EQUIVALENCE     (a(1)    ,ad(1)   ) , (b(1)    ,bd(1)  ) ,  &
     (d(1)    ,dd(1)   ) , (filea(2),m      ) ,  &
     (filea(3),n,rowa  ) , (filea(5),typea  ) ,  &
     (fileb(2),q       ) , (fileb(3),r      ) ,  &
     (fileb(5),typeb   ) , (filec(5),typec  ) ,  &
     (filed(5),typd    ) , (nzz     ,buf1   ) ,  &
     (acoln   ,arown   ) , (filec(7),densc  )
 EQUIVALENCE     (BLOCK(2),TYPE    ) , (BLOCK(3),FORM   ) ,  &
     (BLOCK(4),row     ) , (BLOCK(5),point  ) , (BLOCK(6),nbrstr  ) ,  &
     (BLOCK(8),flag    ) , (xnd(1)  ,xns(1) ) ,  &
     (acol1   ,arow1   ) , (acol    ,arow   ) , (mpy(1)  ,NAME(1) )
 
 DATA    NAME  / 4HMPYA, 4HD   /, jbegn /  4HBEGN/, jend  / 3HEND/  &
     time1 / 0. /  , time2 /  0.    /, zero  /  4*0   /,  &
     method/ 4H1 nt, 4H1 t ,  4H2 nt,  4H2 t  , 4H3 t , 3H4 t/
!WKBI 9/93
 DATA prntyp / 2HRS, 2HRD, 2HCS, 2HCD /
!NVXNB
 IF (typea == 0) typea = iprec
 IF (typeb == 0) typeb = iprec
 IF (typec == 0) typec = iprec
!NVXNE
!WKBNB 7/94 SPR94008
 itypea = typea
 itypeb = typeb
 itypec = typec
!WKBNE 7/94 SPR94008
 
!     CHECK TO SEE IF THE INPUT MATRICES ARE CONFORMABLE
 
 CALL sswtch (19,l19)
 CALL sswtch (41,l41)
 nogo = 0
 FILE = 0
 noab = 0
 IF (filea(6) == 0 .OR. fileb(6) == 0) noab = 1
 irowb = filea(2)
 irowc = filea(3)
 IF (t /= 0) t = 1
 IF (t == 0) GO TO 30
 irowb = filea(3)
 irowc = filea(2)
 30 IF (noab == 1) GO TO 50
 IF (fileb(3) /= irowb) nogo = 1
 IF (filec(1) <= 0) GO TO 40
 IF (filec(2) /= fileb(2) .OR. filec(3) /= irowc) nogo = 1
 40 IF (nogo == 1) GO TO 560
 
!     PERFORM GENERAL INITIALIZATION
 
 50 mpy(3) = jbegn
 IF (filed(1) > 0) CALL conmsg (mpy,3,0)
 nout  = lout
 
!  -- USE SINGLE PRECISION ON MACHINES WITH 60 OR 64 BITS PER WORD
 
 IF (nbpw >= 60) prec1 = 1
 opb   = rdrew
 opc   = rdrew
 op2   = wrtrew
 op    = cls
 cfile = filec(1)
 IF (cfile == 0) typec = 1
 b(2)  = 0.
 b(3)  = 0.
 b(4)  = 0.
 typd1 = typd
 one1  = 1
 one2  = 1
 p     = n
 IF (t /= 0) p = m
 pp1   = p
 incr1 = 1
 IF (cfile == 0  .OR. filec(6) == 0) cfile = 0
 IF (fileb(6) == 0 .AND. cfile == 0) pp1   = 1
 incr2    = 1
 filed(2) = 0
 filed(6) = 0
 filed(7) = 0
 mpass3   = 0
 time3    = 1.0E+10
 prec     = prec1
 IF (prec /= 2) prec = 1
 IF (prec1 == 0 .AND. (prc(typea) == 2 .OR. prc(typeb) == 2 .OR.  &
     prc(typec) == 2)) prec = 2
 
!     ELIMINATE METHOD THREE FROM SELECTION FOR THIS BAD CASE
!     (I.E. TRANSPOSE AND MIXED MATRIX PRECISION)
 
 it = t
 IF (it /= 0 .AND. prec == 1 .AND. prc(typeb) == 2) it = 0
 IF (it /= t .AND. l19 /= 0) WRITE (nout,60) typea,typeb,typec
 60 FORMAT ('0METHOD 3T IS ELIMINATED FROM SELECTION/MPYAD@60',/1X,  &
     'MATRIX TYPES A,B,C =',3I3)
 
!     COMPUTE TYPE AND PRECISION OF D MATRIX
!     RCD    = 1 FOR REAL,   2 FOR COMPLEX
!     PREC   = 1 FOR SINGLE, 2 FOR DOUBLE
!     TYPED  = 1 FOR RSP, 2 FOR RDP, 3 FOR CSP, AND 4 FOR CDP
!     PRC(1) = 1 FOR S.P.   PRC(2) = 2 FOR D.P.
 
 rcd = 0
 IF (prec == 2) GO TO 70
 IF (andf(typea,1) == 0) typea = typea - 1
 IF (andf(typeb,1) == 0) typeb = typeb - 1
 IF (andf(typec,1) == 0) typec = typec - 1
 70 IF (typea > 2 .OR. typeb > 2 .OR. typec > 2) rcd = 2
 typed = rcd + prec
 IF (rcd == 0) rcd = 1
 
!     RCA/B/D   = 1 IF A/B/D IS REAL, = 2 IF A/B/D IS COMPLEX
!     NWDA/B/D  = NUMBER OF WORDS PER ELEMENT OF A/B/D
!     NBX/DX    = NUMBER OF ELEMENTS PER COLUMN OF B/C ORD
!     NB/D      = NUMBER OF WORDS PER COLUMN OF B/C OR D
!     NZZ       = BUF1 = POINTER TO FIRST GINO BUFFER
!     BUF2/3    = POINTER TO SECOND AND THIRD GINO BUFFERS
!     JJ        = MAX. NO. OF COLNS OF B AND D THAT MAY BE HELD IN CORE
!     MPASS1/2/3 = NUMBER OF PASSES REQUIRED FOR METHOD ONE/TWO/THREE
!     JZB/JZDB  = POINTER TO FIRST ELEMENT OF B FOR SP/DP REFERENCE
!     JB        = POINTER TO FIRST ELEMENT OF B FOR PRECISION OF PROBLEM
!     ACORE     = POINTER TO FIRST WORD FOR STORAGE OF PACKED COLUMNS
!                 OF A MATRIX FOR METHOD TWO
!     KSYS58    = SYSTEM(58), METHOD REQUESTED BY USER IF IT IS NON-ZERO
 
 
!     TURN TRANSPOSE FLAG OFF IF INPUT MATRIX A IS SYMMETRIC, AND SURELY
!     THAT COLUMNS EQUEL ROWS, AND DIAG 41 IS OFF.
 
!     IF INPUT A OR B IS DIAGONAL, ROW VECTOR, OR IDENTITY MATRICES,
!     MATRICES ARE NOT IN MIXED PRECISTION TYPES, AND DIAG 41 FLAG IS
!     OFF AND SYSTEM(94) IS NOT 1, BRANCH OFF TO SPECIAL SUBROUTINE
!     MPY-D-R-I
 
 k = filea(4)
 IF (k == 6 .AND. m == n .AND. l41 == 0) t = 0
!        SYMMETRIC    COLN=ROW     DIAG41 OFF
 IF (l41 == 1 .OR. MOD(ksystm(94),10) == 1) GO TO 80
 j = fileb(4)
 IF (k /= 3 .AND. k /= 7 .AND. k /= 8 .AND.  &
     j /= 3 .AND. j /= 7 .AND. j /= 8) GO TO 80
!         DIAGONAL     ROW VCTR     IDENTITY
 
 IF (typea /= typeb .OR. typea /= typd) GO TO 80
 k = MAX0(m,n,q,r)
 j = k*2 + 1
 k = k + 1
 CALL mpydri (z,z,z(k),z(k),z(j),z(j))
 GO TO 380
 
 80 rcb   = rc(typeb)
 nbx   = r*rcb
 nwdb  = nwds(typeb)
 nwdb1 = nwdb + 1
 nb    = r*nwdb
 ndx   = p*rcd
 nd    = p*nwds(typed)
 nzz   = IABS(nz) - sysbuf + 1
 buf2  = buf1 - sysbuf
 buf3  = buf2 - sysbuf
 buf4  = buf3 - sysbuf
 jj    = (nzz-1)/(nb+nd)
 icrq  = nb + nd - nzz + 1
 IF (icrq > 0) GO TO 530
 mpass1= (q-1)/jj + 1
 jzb   = jj*nd  + 1
 jzdb  = jj*ndx + 1
 jb    = jzb
 IF (prc(typeb) == 2) jb = jzdb
 nwda  = nwds(typea)
 prca  = prc(typea)
 na    = nwda*n
 nwda1 = nwda + 1
 nwdd  = nwds(typed)
 acore = nd + 1
 IF (t /= 0) acore = nb + 1
 acore = ((acore+1)/2)*2 + 1
 IF (signab == 0 .AND. prec == 1 .AND. (prc(typea) == 2 .OR.  &
     prc(typeb) == 2)) typed = rcd + 1
 IF (noab == 1   .OR. signab == 0) GO TO 1100
 IF (signab == 1 .OR. signab == -1) GO TO 100
 WRITE  (mout,90)
 90 FORMAT ('0*** USER FATAL MESSAGE 2398, MPYAD REQUIRES SIGN OF ',  &
     'A*B TO BE -1, 0, OR +1')
 GO TO 540
 100 CALL mpyq (z)
 
!     CALCULATE ESTIMATED EXECUTION TIMES AND SELECT METHOD.
 
 ncore = buf3 - acore
 icrq  = -ncore
 IF (icrq > 0) GO TO 530
 core  = FLOAT(ncore/nwda)
 fn    = filea(2)
 fm    = filea(3)
 fp    = fileb(2)
 rhoa  = AMIN1(1.e-4*FLOAT(filea(7)),1.0)
 rhob  = AMIN1(1.e-4*FLOAT(fileb(7)),1.0)
 rhoc  = AMIN1(1.e-4*FLOAT(filec(7)),1.0)
 rhod  = AMAX1(rhoa,rhob)
 arith = fm*fn*(tmt(typed) + (1.0-rhoa)*tml(typed))
 aterm = (fm*rhoa+5.0)*fn*tmipak
 bterm = FLOAT(r)*fp*0.5*(1.0+rhob)*tmupak
 dterm = fm*fp*0.5*(1.0+rhod)*tmpak
 cterm = 0
 IF (cfile /= 0) cterm = fm*fp*0.5*(1.0+rhoc)*tmupak
 time1 = (fm*fn*fp*rhoa*tmt(typed) + FLOAT(mpass1)*aterm + bterm  &
     +  dterm + cterm)*1.0E-6
 
 mpass2= (2.0-rhoa)*fm*fn*rhoa/core + 1.0
 fr    = mpass2
 IF (t /= 0) GO TO 110
 time2 = (fp*rhoa*rhob*arith + aterm  &
     + (fr+1.0)/2.0*(fn*rhob+10.0)*fp*tmipak + fr*dterm  &
     + (fr-1.0)*0.5*fm*fp*(1.0+rhod)*tmupak + cterm)*1.0E-6
 GO TO 120
 
 110 fnt   = fn*fm*rhob
 p1    = AMIN1((fnt/FLOAT(fileb(6))+fp)/2.0,fnt,fp)
 fp1   = p1
 cterm2= 0.
 IF (cfile /= 0) cterm2 = (fn*rhoc+5.0)*fp*tmipak
 bterm = fm*fp*0.5*(1.0+rhob)*tmupak
 dterm2= (fn*rhod+5.0)*fp
 time2 = (fp1*rhoa*arith + (fm*rhoa+5.0)*fn*tmipak + fr*bterm  &
     + (fr+1.0)/2.0*dterm2*tmbpak + (fr-1.0)/2.0*dterm2*tmipak + cterm2)*1.0E-6
 
 bufi  = buf4
 IF (filec(1) == 0) bufi = buf3
 nbrrow= MIN0((bufi-orf(nd+1,1))/na,m)
 mpass3= (m-1)/nbrrow + 1
 fr    = mpass3
 time3 = (fm*fn*fp*rhob*tmt(typed) + fm*fn*0.5*(1.0+rhoa)*tmupak  &
     + fr*fp*(fn*rhob+5.0 )*tmipak + (fr+1.0)/4.0*fn*fp*(1.0+rhod)*tmpak  &
     + (fr-1.0)/4.0*fn*fp*(1.0+rhod)*tmupak + cterm2)*1.e-6
 120 CALL tmtogo (itimgo)
 IF (core <= 0.0) time2 = AMAX1(time1,time3) + 1.0
 time  = AMIN1(time1,time2,time3)
 itime = time + 1
 IF (itimgo <= itime .AND. filed(1) > 0) GO TO 550
 
!     PRINT TIMING MESSAGE AND IF OUTPUT FILE IS PURGED RETURN
 
 ielems = fn*fm*rhoa + 0.5
 jelems = FLOAT(r)*fp*rhob
!WKBNB 9/93
 IF(l19 == 0) GO TO 137
 CALL fname ( filea, namea )
 CALL fname ( fileb, nameb )
 CALL fname ( filec, namec )
 CALL fname ( filed, NAMED )
 WRITE( nout,136, IOSTAT=ierr )
!WKBR 7/94/SPR 94008 *         NAMEA, N, M, IELEMS, RHOA, PRNTYP( TYPEA )
!KWBR 7/94 SPR 94008 *,        NAMEB, R, Q, JELEMS, RHOB, PRNTYP( TYPEB )  &
 namea, n, m, ielems, rhoa, prntyp( itypea )  &
     ,        nameb, r, q, jelems, rhob, prntyp( itypeb )
 136 FORMAT(  &
     '  /-----------------------------------------------------------/' ,/  &
     ,'  /     MATRIX      ROWS   COLS     TERMS  DENS    TYPE       /' ,/  &
     ,'  /-----------------------------------------------------------/' ,/  &
     ,'  /  A- ',2A4,i8,i7,i10,f7.4, 5X, a2 ,/  &
     ,'  /  B- ',2A4,i8,i7,i10,f7.4, 5X, a2 )
 ielems = fn*fm*rhoc + .5
 IF (cfile == 0) GO TO 11140
 WRITE( nout,11136, IOSTAT=ierr )  &
     namec, filec(3), filec(2), ielems, rhoc, prntyp(itypec)
 11136 FORMAT( '  /  C- ',2A4,i8,i7,i10, f7.4, 5X, a2 )
 11140 WRITE( nout, 11137 ) NAMED, prntyp(typed)
 11137 FORMAT('  /  D- ',2A4,8X, 7X, 10X, 7X,   5X, a2 )
 WRITE( nout, 11138 ) signab, signc, t, core, mpass1,mpass2,  &
     mpass3, time1, time2, time3
 11138 FORMAT('  /  SIGNAB =',i4,'  SIGNC =',i4,'  TIME EST=',i9  &
     ,      ' MEMORY =',f8.0  &
     ,/,    '  /  MPASS1 =',i4, '  MPASS2=',i4, '  MPASS3=',i4  &
     ,/,    '  /  TIME1  =',e9.2,' TIME2=',e9.2,' TIME3=',e9.2,/  &
     ,'  /-----------------------------------------------------------/' )
 137 CONTINUE
!WKBNE 9/93
 
 180 IF (filed(1) < 0) GO TO 1600
 
 j = ksys58
 IF (j < 0 .OR. j > 3 .OR. (j == 3 .AND. it == 0)) j = 0
 IF (j  /= 0) THEN
    SELECT CASE ( j )
     CASE (    1)
       GO TO 200
     CASE (    2)
       GO TO 600
     CASE (    3)
       GO TO 1300
   END SELECT
 END IF
 IF (it /= 0) GO TO 190
!WKBNB 2/95 NCL93004
 IF ( mpass1 < mpass2 ) GO TO 200
 IF ( mpass2 < mpass1 ) GO TO 600
!WKBNE 2/95 NCL93004
 IF (time1 < time2) GO TO 200
 GO TO 600
!WKBD 2/95 NCL93004 190 IF (TIME1.LT.TIME2 .AND. TIME1.LT.TIME3) GO TO 200
!WKBNB 2/95 NCL93004
 190 CONTINUE
 IF ( mpass1 < mpass2 .AND.  &
     ( time1 < time3 .OR. mpass1 < mpass3 ) ) GO TO 200
 IF ( mpass2 < mpass1 .AND.  &
     ( time2 < time3 .OR. mpass2 < mpass3 ) ) GO TO 200
 IF ( time1 < time2 .AND. time1 < time3 ) GO TO 200
 IF (time2 < time3) GO TO 600
!WKBNE 2/95 NCL93004
 GO TO 1300
 
!               *********************
!               *                   *
!               *    METHOD  ONE    *
!               *    MPY1NT $ 1T    *
!               *                   *
!               *********************
 
!     BUILD MATRIX PRODUCT JMAX COLUMNS PER PASS OF A MATRIX
!     WHERE JMAX=JJ EXCEPT ON FINAL PASS
 
 200 jcol = 1
 230 WRITE  (nout,240) method(t+1),mpass1,time1
 240 FORMAT ('    METHOD TO BE USED:',a4,', NBR PASSES =',i4,  &
     ',  EST. TIME =',f9.1)
 250 jmax  = MIN0(jcol+jj-1,q)
 IF (jmax == q) op = clsrew
 jmax1 = jmax  - jcol
 jmax  = jmax1 + 1
 jmax1x= jmax1*ndx
 IF (fileb(6) == 0) GO TO 270
 
!     READ AND UNPACK JMAX COLUMNS OF THE B MATRIX
 
 FILE = fileb(1)
 jz   = jzb
 typebd = typeb*signab
 nbd  = nb
 opbc = opb
 pp2  = r
 ASSIGN 270 TO mm
 GO TO 400
 
!     READ AND UNPACK JMAX COLUMNS OF THE C MATRIX
 
 270 FILE = filec(1)
 jz   = 1
 typebd = typed*signc
 nbd  = nd
 opbc = opc
 pp2  = p
 ASSIGN 280 TO mm
 GO TO 400
 
!     OPEN AND POSITION A MATRIX TO FIRST COLUMN
 
 280 IF (fileb(6) == 0) GO TO 340
 FILE = filea(1)
 CALL OPEN (*500,filea,z(nzz),rdrew)
 290 CALL fwdrec (*510,filea)
 
!     SET POINTERS
!     L   = COLUMN NUMBER
!     LL  = POINTER TO LTH ROW OF B MATRIX
!     LLL = POINTER TO LTH ROW OF D MATRIX
 
 l   = 1
 ll  = jb
 lll = 1
 
!     CALL INTPK TO INITIATE READING THE LTH COLUMN OF THE A MATRIX
!     IF COLUMN IS NULL, BYPASS ARITHMETIC
 
 310 CALL intpk (*320,filea,0,typed,0)
 
!     FORM EITHER  A(I,L)*B(L,J) + D(I,J)
!              OR  A(L,I)*B(I,J) + D(L,J)
!           WHERE  J RUNS ACROSS COLUMNS OF B AND D NOW IN CORE
 
 CALL mpy1v (zz,z,zd)
 
!     POSITION POINTERS FOR NEXT COLUMN OF A
 
 320 ll  = ll + rcb
 lll = lll+ rcd
 l   = l  + 1
 IF (l <= m) GO TO 310
 
!     CLOSE AND REWIND FILE CONTAINING A MATRIX
 
 CALL CLOSE (filea,clsrew)
 
!     OPEN FILE CONTAINING D MATRIX TO WRITE
 
 340 FILE = filed(1)
 CALL OPEN (*500,filed,z(nzz),op2)
 
!     IF FIRST COLUMNS OF D, WRITE HEADER
 
 IF (op2 == wrt) GO TO 360
 CALL fname (filed,bcd)
 CALL WRITE (filed,bcd,2,1)
 
!     PACK AND WRITE JMAX COLUMNS OF THE D MATRIX
 
 360 jz = 1
 DO  j = 1,jmax
   CALL pack (z(jz),filed,filed)
   jz = jz + nd
 END DO
 
!     TEST FOR END OF MULTIPLICATION
!     CLOSE FILE CONTAINING D MATRIX
 
 CALL CLOSE (filed,op)
 
!     SET OP FLAGS FOR OPEN CALLS FOR NEXT PASS
 
 opb = rd
 opc = rd
 op2 = wrt
 
 jcol = jcol + jj
 IF (jcol <= q) GO TO 250
 380 mpy(3) = jend
 CALL conmsg (mpy,3,0)
 GO TO 1600
 
!     INTERNAL SUBROUTINE TO READ JMAX COLUMNS OF THE B OR C MATRICES
!     ELEMENTS ARE SET TO ZERO IF COLUMN IS NULL OR MATRIX ABSENT
 
!     OPEN AND POSITION FILE IF MATRIX IS PRESENT
 
 400 IF (FILE) 410,420,410
 410 CALL OPEN (*500,FILE,z(nzz),opbc)
 IF (jcol /= 1) GO TO 420
 CALL fwdrec (*510,FILE)
 
!     LOOP THROUGH JMAX COLUMNS OF MATRIX
 
 420 DO  j = 1,jmax
   
!     UNPACK THE JTH COLUMN IF MATRIX IS PRESENT
   
   IF (FILE) 440,450,440
   440 CALL unpack (*450,FILE,z(jz))
   GO TO 470
   
!     ZERO COLUMN
   
   450 k2 = jz + nbd - 1
   DO  k = jz,k2
     z(k) = 0.
   END DO
   
!     POSITION POINTERS TO NEXT COLUMN OF MATRIX
   
   470 jz = jz + nbd
 END DO
 
!     CLOSE FILE IF MATRIX IS PRESENT
 
 IF (FILE) 480,490,480
 480 CALL CLOSE (FILE,op)
 
!     RETURN
 
 490 GO TO mm, (270,280)
 
 
!     ERROR CONDITIONS
 
 500 mm = -1
 GO TO 520
 510 mm = -2
 520 CALL mesage (mm,FILE,NAME)
 GO TO 1600
 
 530 mm = -8
 FILE = icrq
 GO TO 520
 540 mm = -37
 GO TO 520
 550 mm = -50
 FILE = itime
 GO TO 520
 560 CALL fname (filea,zz(1))
 CALL fname (fileb,zz(3))
 CALL fname (filec,zz(5))
 IF (filec(2) /= fileb(2) .OR. filec(3) /= irowc) nogo = 1
 WRITE (mout,570) zz(1),zz(2),filea(2),filea(3),zz(3),zz(4),  &
     fileb(2),fileb(3),zz(5),zz(6),filec(2),irowc
 570 FORMAT (3(4X,2A4,2I7))
 mm = -55
 GO TO 520
 
 
!               *********************
!               *                   *
!               *    METHOD  TWO    *
!               *    MPY2NT $ 2T    *
!               *    AND   MPY4T    *
!               *                   *
!               *********************
 
 600 CONTINUE
 
!     INITIALIZE FOR METHODS 2NT, 2T AND 4T.
!     METHOD 4T DOES NOT HANDLE COMPLEX MATRIX-D FROM REAL MATRICES A
!     AND B (LL4 = 6).
 
 mt4   = 0
 mt2   = t
 IF (MOD(ksystm(94),100)/10 == 1) GO TO 620
 IF (t == 0 .OR. l41 == 1 .OR. ll4 == 6) GO TO 620
 IF (mpass2 <= 2 .OR. cfile == 0 .OR. densc < 700) GO TO 620
 mt2   = 0
 mt4   = 2
 acore = nb + nd + 1
 acore = ((acore+1)/2)*2 + 1
 jzb   = nd + 1
 jb    = nd/prec1 + 1
 620 dfile = filed(1)
 efile = scrtch
 BLOCK(1) = filea(1)
 cfile = filec(1)
 opa   = rdrew
 typec = typed*signc
 firstl= buf3 - 1
 640 WRITE (nout,240) method(t+3+mt4),mpass2,time2
 
!     BEGIN PASS
 
!     OPEN DFILE TO WRITE.
!     READ AS MANY COLUMNS (OR ROWS) OF A AS CAN BE HELD
!     IN CORE IN PACKED FORM ON THIS PASS.
 
 acol1 = 1
 650 FILE = dfile
 CALL OPEN  (*500,dfile,z(buf3),wrtrew)
 CALL fname (filed(1),bcd)
 CALL WRITE (dfile,bcd,2,1)
 filed(2) = 0
 filed(6) = 0
 filed(7) = 0
 FILE = filea(1)
 CALL gopen (filea,z(buf2),opa)
 apoint = acore
 l    = firstl
 acol = acol1
!WKBR 9/94   660 IF ( (APOINT+NA+2) .GE. L-2) GO TO 530
! ABOVE CHECK WAS OVER-ZEALOUS IN CHECKING FOR AVAILABLE MEMORY
! BECAUSE OF THE CHECK TWO LINES AFTER STATEMENT 670
 660 IF ( (apoint+2) >= l-2) GO TO 750
 zz(l  ) = 0
 zz(l-1) = 0
 BLOCK(8)=-1
 CALL getstr (*730,BLOCK)
 incra = 1
 IF (prc(TYPE) == 2 .AND. prc(typea) == 1) incra = 2
 zz(l) = apoint
 670 kr1 = apoint + 2
 krn = kr1 + nbrstr*nwda - 1
 IF (krn >= l-2) GO TO 740
 
!     MOVE STRING FROM BUFFER TO CORE AND COMPLETE STRING DEFINITION
!     WORDS
 
 IF (prc(TYPE) /= 2 .OR. prc(typea) /= 1) GO TO 690
 
!  -- THIS CODE NECESSARY FOR UNIVAC DOUBLE PRECISION TO SINGLE PRC.
 
 inc   = 1
 incra = 1
 IF (TYPE == 4) inc = 2
 krn = kr1 + nbrstr*inc - 1
 DO  ii = kr1,krn
   z(ii) = xnd(point)
   point = point + incra
 END DO
 GO TO 710
 690 IF (prc(TYPE) == 2) point = point*2 - 1
 DO  ii = kr1,krn
   z(ii) = xns(point)
   point = point + incra
 END DO
 710 zz(apoint  ) = row
 zz(apoint+1) = nbrstr
 zz(l-1) = zz(l-1) + 1
 apoint  = krn + 1
 
!     GET NEXT STRING DEFINITION
 
 CALL endget (BLOCK)
 CALL getstr (*730,BLOCK)
 GO TO 670
 
!     END-OF-COLUMN -
!     SAVE LAST NON-ZERO TERM POSTION FOR MTHOD 4T, THEN
!     TEST FOR ALL COLUMNS
 
!     SINCE GINO IS MDS, MAKE SURE THAT THE LAST VALUES IN NBRSTR AND
!     ROW HERE ARE STILL VALID. OTHERWISE THEY MUST BE SAVED FIRST (AT
!     710) AND USED ON NEXT LINE.
 
 730 IF (mt4 == 2) zz(l-1) = orf(zz(l-1),lshift(row+nbrstr-1,ihalf))
!                                    NBR + LAST NON-ZERO TERM COLUMN NO.
 
 l = l - 2
 acol = acol + 1
 IF (acol <= m) GO TO 660
 
!     ALL COLUMNS OF A ARE IN - THIS IS THE LAST PASS
 
 acoln = m
 CALL CLOSE (filea(1),clsrew)
 GO TO 760
 
!     ALL COLUMNS OF A WILL NOT FIT ON THIS PASS.
 
 740 CALL bckrec (filea(1))
 750 CALL CLOSE  (filea(1),cls)
 acoln = acol - 1
 
!     IF CFILE IS PRESENT, OPEN IT.
!     IF THIS IS THE FIRST PASS, SKIP HEADER RECORD.
!     OPEN BFILE AND SKIP HEADER RECORD.
!     INITIALIZE COLUMN (OR ORW) COUNTER, BCOL, TO 1, AND BRANCH ON T.
 
 760 IF (cfile == 0) GO TO 770
 FILE = cfile
 CALL OPEN (*500,cfile,z(buf1),rdrew)
 CALL fwdrec (*510,cfile)
 770 FILE = fileb(1)
 CALL OPEN (*500,fileb(1),z(buf2),rdrew)
 CALL fwdrec (*510,fileb(1))
 bcol = 1
 780 IF (mt2 == 1) GO TO 900
 
!     UNPACK A COLUMN OF C.
 
 IF (cfile == 0) GO TO 810
 typebd = typec
 IF (mt4 /= 0) one2 = 1
 pp2    = p
 CALL unpack (*810,cfile,z)
 GO TO 830
 810 DO  ii = 1,nd
   z(ii) = 0.
 END DO
 830 IF (mt4 /= 0) GO TO 850
 
!     NON-TRANSPOSE CASE, METHOD 2NT
!     ==============================
 
!     INITIATE INTERPRETATION OF A COLUMN OF B.
 
!     ITYPSG = TYPED*SIGNAB
 itypsg = typeb*signab
 CALL intpk (*860,fileb(1),0,itypsg,0)
 
!     FOR EACH NON-ZERO ELEMENT B(I) IN THE CURRENT COLMN OF B SUCH
!     THAT FOR I.GE.ACOL1 .AND I.LE.ACOLN, FORM ALL PRODUCTS OF
!     D(K,I) = A(K,I)*B(I) + C(K,I)
 
 CALL mpy2nv (zz,z,zd)
 GO TO 860
 
!     TRNASPOSE CASE, METHOD 4T
!     =========================
 
!     UNPACK A BANDED COLUMN OF MATRIX B, RANGING FROM ONE2 THRU PP2.
!     FOR THE RANGE  MAX0(ONE2,ACOL1) THRU MIN0(PP2,ACOLN), FORM ALL
!     PRODUCTS
!     D(I,K) = A(I,J)*B(J,K) + C(I,K)
 
 850 typebd = typeb*signab
 one2   = 0
 CALL unpack (*860,fileb,z(jzb))
 
!     WE HAVE HERE -
!     ACLO1, ACOLN = COLUMNS OF MATRIX A IN CORE
!     BCOL = CURRENTLY WE ARE WORKING ON THE BCOL COLUMN OF MATRIX B,
!            WHICH IS ALSO THE WORKING COLUMNS OF MATRIX D AND MATRIX C
!     Z(JZB) THRU Z(ACORE-1) CONTAIN THE BCOL COLUMN OF MATRIX B
 
 CALL mpy4t (z,z,z)
 
!     PACK CURRENT COLUMN ONTO DFILE FOR BOTH 2NT AND 4T METHOD, AND
!     GO TO TEST FOR END OF PASS.
 
 860 CALL pack (z,dfile,filed)
 GO TO 980
 
!     TRANSPOSE CASE, METHOD 2T
!     =========================
 
!     INITIATE BUILDING OF A PACKED COLUMN OF D.
!     UNPACK A COLUMN OF B IN CORE. IF NULL, COPY COLUMN FROM C TO D.
!     INITIATE INTERPRETATION OF A COLUMN OF C.
 
 900 CALL bldpk (typed,typd,dfile,0,0)
 typebd = typeb*signab
 pp2  = r
 CALL unpack (*910,fileb(1),z)
 eol  = 1
 crow = 16777215
!            16777215 = 2**24 - 1
 
 IF (cfile == 0) GO TO 960
 CALL intpk (*960,cfile,0,typec,0)
 crow = 0
 GO TO 960
 910 IF (cfile == 0) GO TO 970
 CALL intpk (*970,cfile,0,typec,0)
 920 CALL zntpki
 crow = ip
 930 DO  ii = 1,nwdd
   d(ii) = a(ii)
 END DO
 drow = crow
 CALL zblpki
 950 IF (eol == 0) GO TO 920
 GO TO 970
 
!     FOR ALL NON-NULL ROWS OF A IN CORE, FORM A(I,J)*B(J) + C(I)
 
 960 CALL mpy2tv (zz,z,zd)
 IF (arown == m .OR. crow == 16777215) GO TO 970
 IF (crow > arown) GO TO 930
 GO TO 950
 
!     TERMINATE CURRENT COLUMN OF D.
 
 970 CALL bldpkn (dfile,0,filed)
 
!     BOTH TRANSPOSE (2T AND 4T) AND NON-TRANSPOSE (2NT) CASES
 
!     TEST FOR COMPLETION OF PASS. IF COMPLETE, TEST ALL PASSES.
 
 980 bcol = bcol + 1
 IF (bcol <= q) GO TO 780
 CALL CLOSE (fileb,clsrew)
 IF (cfile /= 0) CALL CLOSE (cfile,clsrew)
 CALL CLOSE (dfile,clsrew)
 IF (acoln == m) GO TO 1010
 
!     NOT LAST PASS - SWITCH C AND D FILES AND CONTINUE
 
 opa   = rd
 typec = typed
 IF (acol1 == 1) GO TO 990
 k     = cfile
 cfile = dfile
 dfile = k
 GO TO 1000
 990 cfile = dfile
 dfile = efile
 1000 acol1 = acoln + 1
 GO TO 650
 
!     LAST PASS -
!     MAKE SURE D MATRIX IS ON PROPER FILE.
!     IF NOT, SWITCH FIST AND FIAT UNIT NBRS IN /XFIAT/
 
 1010 IF (dfile /= filed(1)) CALL filswi (dfile,filed)
 GO TO 380
 
!     A MATRIX OR B MATRIX IS NULL - COPY C MATRIX TO D MATRIX
 
 1100 time = 0.0
 IF (filed(1) < 0) GO TO 1600
 IF (q <= 0) q = filec(2)
 filed(2) = 0
 filed(6) = 0
 filed(7) = 0
 WRITE  (nout,1140)
 1140 FORMAT ('             MPYAD - NULL MATRIX PRODUCT')
 CALL gopen (filed,z(buf1),wrtrew)
 IF (cfile == 0) GO TO 1150
 IF (typec == signc*typd) GO TO 1170
 GO TO 1200
 
!     PACK NULL COLUMNS BECAUSE C MATRIX IS NULL
 
 1150 pp1 = 1
 DO  acol = 1,q
   CALL pack (zero,filed,filed)
 END DO
 GO TO 1190
 
!     USE CPYSTR TO COPY C TO D
 
 1170 BLOCK(1) = cfile
 blk(1)   = filed(1)
 CALL gopen (cfile,z(buf2),rdrew)
 DO  ii = 1,q
   CALL cpystr (BLOCK,blk,0,0)
 END DO
 CALL CLOSE (cfile,clsrew)
 filed(2) = q
 filed(5) = filec(5)
 filed(6) = filec(6)
 filed(7) = filec(7)
 1190 IF (filec(1) > 0) filed(4) = filec(4)
 CALL CLOSE (filed,clsrew)
 GO TO 380
 
!     USE INTPK/BLDPK TO COPY C TO D BECAUSE TYPES CONFLICT
 
 1200 CALL gopen (cfile,z(buf2),rdrew)
 DO  ii = 1,q
   CALL bldpk (typd,typd,filed,BLOCK,1)
   itypsg = signc*typd
   CALL intpk (*1220,filec,0,itypsg,0)
   1210 CALL zntpki
   CALL bldpki (a,ip,filed,BLOCK)
   IF (eol == 0) GO TO 1210
   1220 CALL bldpkn (filed,BLOCK,filed)
 END DO
 CALL CLOSE (cfile,clsrew)
 GO TO 1190
 
 
!               *********************
!               *                   *
!               *    METHOD THREE   *
!               *       MPY3T       *
!               *                   *
!               *********************
 
!     TRANSPOSE CASE ONLY, METHOD 3T
!     ==============================
 
 1300 CONTINUE
 WRITE (nout,240) method(5),mpass3,time3
 BLOCK(1) = fileb(1)
 acore = orf(nd+1,1)
 cfile = scrtch
 dfile = filed(1)
 IF (MOD(mpass3,2) /= 0) GO TO 1340
 cfile = filed(1)
 dfile = scrtch
 1340 arow1 = 1
 last  = .false.
 opa   = rdrew
 
!     BEGIN PASS BY FILLING CORE WITH UNPACKED COLUMNS OF A
 
 1350 arown = MIN0(arow1+nbrrow-1,m)
 IF (arown == m) last = .true.
 CALL gopen (filea,z(buf1),opa)
 typebd = typea*signab
 pp2    = n
 apoint = acore
 DO  arow = arow1,arown
   CALL unpack (*1360,filea,z(apoint))
   GO TO 1380
   1360 k2 = apoint + na - 1
   DO  ii = apoint,k2
     z(ii) = 0.
   END DO
   1380 apoint = apoint + na
 END DO
 ii = cls
 IF (last) ii = clsrew
 CALL CLOSE (filea,ii)
 incra = (arown-arow1)*na
 
!     PREPARE TO PASS B MATRIX AND C MATRIX FROM LAST PASS
 
 IF (arow1 /= 1) CALL gopen (cfile,z(buf2),rdrew)
 CALL gopen (dfile,z(buf3),wrtrew)
 CALL gopen (fileb,z(buf1),rdrew )
 IF (last .AND. filec(1) /= 0) CALL gopen (filec,z(buf4),rdrew)
 filed(2) = 0
 filed(6) = 0
 filed(7) = 0
 typebd   = typed
 pp2 = arown
 k2  = arown*nwdd
 
 DO  bcol = 1,q
   IF (arow1 /= 1) GO TO 1420
   
!     FIRST PASS OR NULL COLUMN ON CFILE - SET COLUMN OF D TO ZERO
   
   1400 DO  ii = 1,k2
     z(ii) = 0.
   END DO
   null  = .true.
   IF (last) GO TO 1430
   GO TO 1500
   
!     INTERMEDIATE PASS OR LAST PASS - UNPACK COLUMN FROM PREVIOUS PASS
   
   1420 CALL unpack (*1400,cfile,z)
   null = .false.
   IF (.NOT.last) GO TO 1500
   
!     LAST PASS - ADD COLUMN FROM C MATRIX (IF PRESENT)
   
   1430 IF (filec(1) == 0) GO TO 1500
   itypsg = typed*signc
   CALL intpk (*1500,filec,0,itypsg,0)
   null = .false.
   1440 CALL zntpki
   SELECT CASE ( typed )
     CASE (    1)
       GO TO 1450
     CASE (    2)
       GO TO 1460
     CASE (    3)
       GO TO 1470
     CASE (    4)
       GO TO 1480
   END SELECT
   1450 z(ip) = z(ip) + a(1)
   GO TO 1490
   1460 zd(ip) = zd(ip) + ad(1)
   GO TO 1490
   1470 z(2*ip-1) = z(2*ip-1) + a(1)
   z(2*ip  ) = z(2*ip  ) + a(2)
   GO TO 1490
   1480 zd(2*ip-1) = zd(2*ip-1) + ad(1)
   zd(2*ip  ) = zd(2*ip  ) + ad(2)
   1490 IF (eol == 0) GO TO 1440
   
!     FOR EACH NON-ZERO TERM B(J) IN THE CURRENT COLUMN OF B FORM
!     D(I,K) = D(I,K) + A(I,J)*B(J,K)
   
   1500 CALL mpy3t (*1510,z(acore),z(acore),z(1),z(1))
   GO TO 1520
   
!     PACK NULL COLUMN
   
   1510 IF (.NOT.null) GO TO 1520
   pp1 = 1
   CALL pack (zero,dfile,filed)
   CYCLE
   
!     PACK NON-NULL COLUMN
   
   1520 pp1 = arown
   CALL pack (z,dfile,filed)
   
!     TEST FOR END OF CURRENT PASS
   
 END DO
 
 IF (arow1 /= 1) CALL CLOSE (cfile,clsrew)
 CALL CLOSE (dfile,clsrew)
 CALL CLOSE (fileb,clsrew)
 IF (last) GO TO 1540
 
!     NOT LAST PASS - SWITCH FILES AND CONTINUE
 
 ii    = cfile
 cfile = dfile
 dfile = ii
 arow1 = arown + 1
 opa   = rd
 GO TO 1350
 
!     LAST PASS - SIGNAL END AND RETURN
 
 1540 IF (filec(1) /= 0) CALL CLOSE (filec,clsrew)
 GO TO 380
 
 1600 RETURN
END SUBROUTINE mpyado
