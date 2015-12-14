SUBROUTINE criggp (n23)
     
    !     ******************************************************************
    !                                                                      *
    !     THIS SUBROUTINE GENERATES COEFFICIENTS FOR RIGID ELEMENTS        *
    !     FOR USE BY SUBROUTINE GP4.  THE DATA SO GENERATED IS             *
    !     COMPATIBLE WITH MPC SET DATA.                                    *
    !                                                                      *
    !     (MODIFIED BY G.CHAN/SPERRY TO REDUCE EXCESSIVE OPENINGS,         *
    !     CLOSINGS, AND READINGS OF THE BGPDT FILE (IN 2ND METHOD).        *
    !     WITHOUT THIS MODIFICATION, A PROBLEM OF 2000 RIGID ELEMENTS,     *
    !     FOR EXAMPLE, WOULD REQUIRE MORE THAN 10,000 OPENS AND 10,000     *
    !     CLOSES AND OVER 10 MILLION CALLS TO SUBROUTINE READ   10/86)     *
    !                                                                      *
    !     (MODIFIED AGAIN BY G.CHAN/UNISYS TO INCLUDE CRROD, CRBAR, CRBE1, *
    !     CRBE2, CRBE3, CRTRPLT, AND CRSPLINE RIGID ELEMENTS    11/88)     *
    !                                                                      *
    !     ******************************************************************
 
    !     EXTERNAL          ORF    ,LSHIFT
 
    INTEGER, INTENT(OUT)                     :: n23

    LOGICAL :: again  ,genre  ,l38    ,debug

    INTEGER :: geomp  ,bgpdt  ,cstm   ,rgt    ,scr1       ,  &
        buf(20),mask16 ,gpoint ,z      ,flag       ,  &
        file   ,ret    ,ret1   ,ic(1)  ,mcode(2)   , buf1   ,buf2   ,buf3   ,buf4
    INTEGER :: crigdr(2),crigd1(2)    ,crigd2(2),crigd3(2),  &
        crtrpt(2),crspli(2)    ,crrod(2) ,crbar(2) , crbe1(2) ,crbe2(2)     ,crbe3(2)
    INTEGER :: rdrew

    DOUBLE PRECISION :: dz(1)

    DIMENSION         rz(1)  ,NAME(2),indcmp(6)
    DIMENSION         a(36)  ,b(6)   ,ib(6)  ,c(18)

    CHARACTER (LEN=29) :: uim
    CHARACTER (LEN=25) :: uwm
    CHARACTER (LEN=23) :: ufm

    COMMON /xmssg /   ufm    ,uwm    ,uim
    COMMON /machin/   mach   ,ihalf  ,jhalf
    COMMON /zzzzzz/   z(1)
    COMMON /gp4fil/   geomp  ,bgpdt  ,cstm   ,rgt    ,scr1
    COMMON /gp4prm/   buf    ,buf1   ,buf2   ,buf3   ,buf4   ,knkl1  ,  &
        mask16 ,nogo   ,gpoint ,kn
    COMMON /names /   rd     ,rdrew  ,wrt    ,wrtrew ,clsrew
    COMMON /system/   ksystm(55)

    EQUIVALENCE      (z(1),rz(1),dz(1))
    EQUIVALENCE      (ic(1)     ,c(1) )
    EQUIVALENCE      (ksystm(2) , nout)
    EQUIVALENCE      (ksystm(55),iprec)

    DATA    crigd1/  5310,   53/, crigd2/  5410,   54/,  &
            crigd3/  8310,   83/, crigdr/  8210,   82/,  &
            crrod /  6510,   65/, crbar /  6610,   66/,  &
            crtrpt/  6710,   67/, crbe1 /  6810,   68/,  &
            crbe2 /  6910,   69/, crbe3 /  7010,   70/,  &
            crspli/  7110,   71/
    DATA    NAME  /  4HCRIG,2HGP/
    DATA    mset  /  4HMSET    /
    DATA    a     /  36*0.     /
    DATA    debug /  .false.   /
    DATA    l38   /  .false.   /
 
 
    !     ****************************************************************
    !      OPEN CORE -
    !                                        ALLOCATED BY
    !     !<--- ALLOCATED BY GP4 --->!<------- CRIGGP -   --------->!
    !     +--------+--------+------+-+----+-----+------   ----+-----+-----+
    !     ! EQEXIN ! EQEXIN !SORTED!U!CSTM!BGPDT!   ...       ! DEP ! GINO!
    !     !1ST REC !2ND REC ! SIL  !S!    !     !             ! SIL !BFFRS!
    !     +--------+--------+------+-+----+-----+------   ----+-----+-----+
    !      1      KN       KM        /    /    / \_KNKL1      /     /
    !                            KNKL2  KNKL3 KNKL4         MU   BUF4
 
    !      OPEN CORE FORMAT, STARTS WITH Z(KNKL2)
    !      (KNKL2 = INITIAL VALUE OF KNKL1)
 
    !      NUMBER OF WORDS                 CONTENTS
 
    !           NCSTM    ***     COORDINATE SYSTEM TRANSFORMATION TABLE
    !           NBGPDT   ***     BASIC GRID POINT DEFINITION TABLE
    !                    *       (ONLY IF ENOUGH OPEN CORE SPACE AVAILABLE)
    !                    ***     SIL 1   ***
    !                    *       SIL 2     *
    !             6      *       SIL 3     *  INDEPENDENT GRID POINT
    !                    *       SIL 4     *
    !                    *       SIL 5     *
    !                    ***     SIL 6   ***
    !                    ***     SIL 1          ***             ***
    !                    *       DEGREE OF FREEDOM*               *
    !                    *       INTERNAL INDEX   *               *
    !                    *       SIL 2            *               *
    !                    *       DEGREE OF FREEDOM*  FIRST DEPEND.*  ALL
    !                    *       INTERNAL INDEX   *  GRID POINT   *  DEPEND.
    !           3*MDEP   *           .            *               *  GRID
    !                    *           .            *               *  POINTS
    !                    *       SIL 6            *               *
    !                    *       DEGREE OF FREEDOM*               *
    !                    ***     INTERNAL INDEX ***             ***
    !                    ***     INDEPENDENT GRID POINT BGPDT TABLE
    !             4      *          WORD 1     COORDINATE SYSTEM ID-INTEGER
    !                    ***        WORD 2-4 = X, Y, Z, IN BASIC SYSTEM-REAL
    !          4*MDBGP   ***     DEPENDENT GRID POINT BGPDT TABLE
    !                    ***
    !         36*MDBGP   *       ROW STORED GG MATRIX (SINGLE PRECISION)
    !                    ***     36 ELEMENTS * NO. DEPEND. GRID PT.
    !          9*IPREC   ***     INDEPEND. GRID PT TRANSFORMATION MAT.-REAL
    !          9*IPREC   ***     DEPEND. GRID PT TRANSFORMATION MATRIX-REAL
    !         36*IPREC   ***     GG MATRIX  -  REAL  36 ELEMENTS
    !                    ***          .        ***
    !                    *            .          *
    !                    *  AVAILABLE OPEN CORE  *
    !                    *            .          *
    !                    *            .          *
    !                    ***          .        ***
    !                    ***
    !           MDEP     *       DEPENDENT SILS
    !                    ***
    !                    ***
    !          BUFFERS   *       GINO BUFFERS
    !                    ***
 
    !     *************************************************************
    !     NOTE  IPREC = 1   SINGLE PRECISION
    !           IPREC = 2   DOUBLE PRECISION
    !           MDEP  =     NUMBER DEPENDENT SILS
    !           MDBGP =     NUMBER DEPENDENT GRID POINTS
    !     *************************************************************
 
    mask15= jhalf/2
    kn2   = kn/2
    ncstm = 0
    knkl2 = knkl1
    kiold = 0
    ibuf1 =-99
    again = .false.
    CALL sswtch (20,j)
    IF (j == 1) debug =.true.
    CALL sswtch (38,j)
    IF (j == 1) l38 =.true.
    CALL page2(-4)
    WRITE  (nout,200) uim
200 FORMAT (a29,' 3113, RIGID ELEMENTS ARE BEING PROCESSED IN GP4',/)
 
    !     OPEN CSTM AND READ INTO CORE, FROM Z(KNKL2) THRU Z(KNKL3)
 
    left = buf4 - knkl1
    FILE = cstm
    CALL OPEN (*300,cstm,z(buf2),rdrew)
    CALL skprec (cstm,1)
    CALL READ (*1230,*270,cstm,z(knkl2),left,1,ncstm)
    GO TO 1280
 
    !     IF CORE WAS FILLED WITHOUT HITTING AN EOR, CALL MESAGE
 
270 IF (iprec == 1) CALL pretrs (z(knkl1),ncstm)
    IF (iprec == 2) CALL pretrd (z(knkl1),ncstm)
    CALL CLOSE (cstm,clsrew)
    GO TO 300
 
    !     IF THERE IS ENOUGH CORE AVAILABLE, OPEN AND READ BGPDT INTO OPEN
    !     CORE, FROM Z(KNKL3+1) THRU Z(KNKL4), CLOSE BGPDT FILE, AND RESET
    !     VARIOUS POINTERS FOR BUILDING UP RGT DATA. (AGAIN=.FALSE.)
    !     THIS METHOD USES ONLY ONE OPEN, ONE CLOSE, AND ONE READ.
 
    !     HOWEVER, IF THERE IS NOT ENOUGH CORE FOR BGPDT DATA AND THE NEEDED
    !     SPACE FOR BUILDING UP RGT DATA, SET AGAIN TO .TRUE., AND REPEAT
    !     DATA PROCESSING BY READING DATA DIRECTLY OFF THE BGPDT FILE EACH
    !     TIME WHEN THE BGPDT DATA IS NEEDED.   THIS SECOND METHOD USES ONLY
    !     ONE OPEN, ONE CLOSE, AND MULTIPLE READS.
 
    !     IN THE SECOND METHOD, TWO POINTERS, KIOLD AND KINEW, ARE USED TO
    !     COMPUTE PRECISELY WHERE TO READ DATA OFF THE BGPDT FILE
 
290 again = .true.
    CALL WRITE  (rgt,0,0,1)
    CALL bckrec (rgt)
    knkl3 = 0
    knkl1 = knkl2
    nbgpdt= knkl1 + ncstm
    CALL CLOSE (bgpdt,clsrew)
300 FILE  = bgpdt
    CALL OPEN (*1210,bgpdt,z(buf2),rdrew)
    CALL fwdrec (*1240,bgpdt)
    kiold = 0
 
    !     CALCULATE STARTING POINT
    !     AND READ BGPDT INTO OPEN CORE
 
    knkl1 = knkl1 + ncstm
    IF (again) GO TO 310
    knkl3 = knkl1
    CALL READ (*1230,*310,bgpdt,z(knkl3+1),buf4-knkl3,1,nbgpdt)
    imhere = 305
    IF (debug) WRITE (nout,1255) imhere
    knkl3 = 0
    nbgpdt= knkl1
    again = .true.
    CALL bckrec (bgpdt)
310 IF (.NOT.again) CALL CLOSE (bgpdt,clsrew)
    knkl4 = knkl3 + nbgpdt
    knkl1 = knkl4 + 1
    mu    = buf4  - 1
    irdg  = 0
    itype = 0
    genre = .false.
 
    !     *************************************************************
 
    !     CRIGD1, CRIDG2, AND CRBE2 RIGID ELEMENTS ARE PROCESSED HERE
 
    !     *************************************************************
 
    !     LOCATE CRIGD1 DATA IN THE INPUT FILE
 
    FILE = geomp
    CALL locate (*500,z(buf1),crigd1,flag)
    irdg = 1
    GO TO 1000
 
    !     LOCATE CRIGD2 DATA ON INPUT FILE
 
500 FILE = geomp
    CALL locate (*600,z(buf1),crigd2,flag)
    irdg = 2
    imhere = 500
    IF (debug) WRITE (nout,4400) imhere
    GO TO 1000
 
    !     LOCATE CRBE2 DATA ON INPUT FILE
 
600 FILE = geomp
    CALL locate (*4000,z(buf1),crbe2,flag)
    irdg = 3
    imhere = 600
    IF (debug) WRITE (nout,4400) imhere
 
1000 CONTINUE
     IF (debug) WRITE (nout,1005) irdg
1005 FORMAT ('0 IRDG/CRIGGP =',i6)
 
     !     READ ELEMENT ID AND INDEPENDENT GRID POINT NUMBER
 
1730 ifile = geomp
     nwds  = 2
     GO TO 1734
1732 ifile = scr1
     nwds  = 9
1734 FILE  = ifile
     CALL READ (*1230,*1240,ifile,buf,nwds,0,flag)
     IF ((debug.OR.l38) .AND. buf(1) /= ibuf1) WRITE (nout,1735) buf(1)
1735 FORMAT (5X,'ELEMENT',i8,' IS BEING PROCESSED')
     IF (.NOT.genre) GO TO 1739
     ibuf1 = buf(1)
 
     !     SET UP INDEPENDENT D.O.F. FOR THE GENERAL RIGID ELEMENTS,
     !     CRIGID3 AND CRBE1, AND ALSO THE CRBAR AND CRTRPLT ELEMENTS
     !     WHICH WERE CONVERTED TO CRIGID3 FORMAT BY IFS3P
 
     DO  i = 1,6
         indcmp(i) = buf(i+2)
     END DO
     itype = buf(9)
     IF (itype /= 0) GO TO 1739
     DO  i = 1,36
         a(i) = 0.0
     END DO
     INDEX = 0
     ilast = 0
     DO  i = 1,6
         IF (indcmp(i) /= i) CYCLE
         j    = 6*ilast + i
         a(j) = 1.0
         ilast= ilast + 1
     END DO
     nind = ilast
 
1739 ASSIGN 1740 TO ret
     ASSIGN 1743 TO ret1
     idr   = buf(1)
     gpoint= buf(2)
     ntype = 1
     GO TO 7060
 
     !     STORE SIL FOR INDEPENDENT DEGREES OF FREEDOM
 
     1740 DO  i=1,6
         z(knkl1+i-1) = gpoint + i - 1
     END DO
1743 kinew = k - 2*kn
     ASSIGN 1750 TO ret
     ASSIGN 1745 TO ret1
 
     !     READ DEPENDENT GRID POINTS
 
     j = knkl1 + 3
     mdbgp = 0
     mdep  = 0
1745 CALL READ (*1230,*1240,ifile,buf,7,0,flag)
     IF (buf(1) == -1) GO TO 1760
     mdbgp = mdbgp + 1
     gpoint= buf(1)
     ntype = 2
     GO TO 7060
1750 CONTINUE
     IF (nogo /= 0) GO TO 1745
 
     !     STORE DEPENDENT GRID POINT SIL, DOF, AND INTERNAL INDEX
 
     DO  i = 1,6
         IF (buf(i+1) == 0) CYCLE
         j = j + 3
         l = j
         z(l) = gpoint + i - 1
         z(l+1) = i
         z(l+2) = k - 2*kn
         mdep   = mdep + 1
     END DO
     GO TO 1745
 
     !     HERE WHEN ALL DEPENDENT GRID POINTS FOR AN ELEMENT HAVE BEEN READ
 
1760 more = 0
     i = knkl1 + 6 + 3*mdep + 4 + 4*mdbgp + (9+9+36*mdbgp+36)*iprec
 
     !     CHECK FOR OPEN CORE AVAILABILITY
 
     imhere = 176
     IF (i      >= mu) GO TO 1250
     IF (buf(2) ==  0) more = 1
     IF (nogo   /=  0) GO TO 3645
 
     !     LOCATE DATA IN BGPDT FOR INDEPENDENT GRID POINT
 
     iopen = knkl1 + 6 + 3*mdep
     IF (again) GO TO 1761
     ki4 = knkl3 + kinew*4
     IF (ki4 > knkl4) GO TO 1290
     z(iopen    ) = z(ki4 -3)
     z(iopen + 1) = z(ki4 -2)
     z(iopen + 2) = z(ki4 -1)
     z(iopen + 3) = z(ki4   )
     GO TO 1763
1761 FILE = bgpdt
     IF (kinew > kiold) GO TO 1762
     CALL bckrec (bgpdt)
     kiold = 0
1762 ki4 = (kinew-kiold-1) * 4
     IF (ki4 > 0) CALL READ (*1230,*1240,bgpdt,buf,-ki4,0,flag)
     CALL READ (*1230,*1240,bgpdt,buf,4,0,flag)
     z(iopen    ) = buf(1)
     z(iopen + 1) = buf(2)
     z(iopen + 2) = buf(3)
     z(iopen + 3) = buf(4)
1763 kiold = kinew
 
     !     SORT DEPENDENT DEGREE OF FREEDOM LIST ON BGPDT REFERENCE NUMBER
 
     i = mdep*3
     CALL sort (0,0,3,3,z(knkl1+6),i)
 
     j = 0
     m = 0
     indx  = knkl1 + 5
     indxx = knkl1 + 6 + 3*mdep + 4
     DO  i = 1,mdep
         k = indx + 3*i
         kinew = z(k)
         IF (kiold == kinew) GO TO 1767
         j = j + 1
   
         !     READ GRID POINT INFORMATION
   
         m = m + 1
         n = indxx + (m-1)*4
         IF (again) GO TO 1764
         ki4 = knkl3 + kinew*4
         IF (ki4 > knkl4) GO TO 1290
         z(n    ) = z(ki4 -3)
         z(n + 1) = z(ki4 -2)
         z(n + 2) = z(ki4 -1)
         z(n + 3) = z(ki4   )
         GO TO 1766
1764     FILE = bgpdt
         IF (kinew > kiold) GO TO 1765
         CALL bckrec (bgpdt)
         kiold = 0
1765     ki4 = (kinew-kiold-1)*4
         IF (ki4 > 0) CALL READ (*1230,*1240,bgpdt,buf,-ki4,0,flag)
         CALL READ (*1230,*1240,bgpdt,buf,4,0,flag)
         z(n    ) = buf(1)
         z(n + 1) = buf(2)
         z(n + 2) = buf(3)
         z(n + 3) = buf(4)
1766     kiold = kinew
1767     z(k)  = j
     END DO
 
     IF (iprec == 2) GO TO  3200
 
     !     FORM REFERENCE GRID POINT TRANSFORMATION MATRIX
 
     iba = knkl1 + 6 + 3*mdep
     ita = iba + 4 + 4*mdbgp + 36*mdbgp
     IF (z(iba) /= 0) CALL transs (rz(iba),rz(ita))
 
     !     PREPARE POINTERS USED TO FORM THE G MATRIX
 
     itb = ita + 9
     itc = itb - 1
 
     !     SET INDEXES FOR TRANSFORMATION MATRIXES AND GG MATRIXES TO
     !     FIRST ELEMENT - 1 FOR SUBROUTINE FORMGG
 
     ita = ita - 1
     ig  = indxx + 4*mdbgp - 1
     igg = ig + (36*mdbgp) + 9 + 9
     indx= knkl1 + 3
     m   = -1
 
     !     BEGIN LOOP TO FORM THE G MATRIX
 
     DO  i = 1,mdep
         k  = indx + i*3
         mm = z(k+2)
         IF (mm == m) GO TO 3030
         ibb = indxx + (mm-1)*4
   
         !     FORM DEPENDENT DEGREE OF FREEDOM TRANSFORMATION MATRIX
   
         IF (z(ibb) /= 0) CALL transs (rz(ibb),rz(itb))
   
         !     FORM THE GG MATRIX
   
         CALL formgg (igg,ita,itc,iba,ibb)
3030 CONTINUE
   
     !     SELECT PROPER ROW BASED ON COMPONENT NUMBER AND STORE IN G
     !     ACCORDING TO PARTITIONING VECTOR OF REFERENCE GRID POINT.
   
     m  = mm
     mm = z(k+1)
     DO  ij = 1,6
         indxxx = igg + (mm-1)*6 + ij
         rz(ig+ij) = rz(indxxx)
     END DO
     ig = ig + 6
 END DO
 GO TO 3300
 
 !     FORM REFERENCE GRID POINT TRANSFORMATION MATRIX (DOUBLE PREC.)
 
3200 ibase = (knkl1 + 6 + 3*mdep + 4 + 4*mdbgp + 36*mdbgp) / 2 + 1
 iba = knkl1 + 6 + 3*mdep
 ita = ibase
 IF (z(iba) /= 0) CALL transd (rz(iba),dz(ita))
 
 !     PREPARE POINTERS USED TO FORM THE G MATRIX
 
 itb = ita + 9
 itc = itb - 1
 
 !     SET INDEXES FOR TRANSFORMATION MATRIXES AND GG MATRIXES TO
 !     FIRST ELEMENT - 1 FOR SUBROUTINE FORMGG
 
 ita = ita - 1
 ig  = indxx + 4*mdbgp - 1
 igg = ibase + 9 + 9 - 1
 indx= knkl1 + 3
 m   = -1
 
 !     BEGIN LOOP TO FORM THE G MATRIX
 
 DO  i = 1, mdep
     k  = indx  + i*3
     mm = z(k+2)
     IF (mm == m) GO TO 3230
     ibb = indxx + (mm-1)*4
   
     !     FORM DEPENDENT DEGREE OF FREEDOM TRANSFORMATION MATRIX
   
     IF (z(ibb) /= 0) CALL transd (rz(ibb),dz(itb))
   
     !     FORM THE GG MATRIX
   
     CALL formg2 (igg,ita,itc,iba,ibb)
3230 CONTINUE
   
     !     SELECT PROPER ROW BASED ON COMPONENT NUMBER AND STORE IN G
   
     m  = mm
     mm = z(k+1)
     DO  ij = 1,6
         indxxx = igg + (mm-1)*6 + ij
         rz(ig+ij) = dz (indxxx)
     END DO
     ig = ig + 6
 END DO
3300 ig = indxx + 4*mdbgp - 1
 
 !     WRITE THE CODED COLUMN-ROW NUMBERS AND ELEMENTS OF THE GM
 !     MATRIX ON RGT FILE SO AS TO MAKE RIGID ELEMENT DATA
 !     COMPATIBLE WITH MPC SET DATA
 !     (REVISED 7/86, CODED COLUMN-ROW NUMBERS ARE NOT USED HERE.
 !     THEY WILL BE RE-CODED IN GP4 IF NEEDED)
 
 k  = 0
 IF (genre .AND. itype == 0) GO TO 3380
 mu = mu - mdep
 
 !     TEST FOR OPEN CORE AVAILABILITY
 
 imhere = 3380
 IF (iopen >= mu) GO TO 1250
3380 CONTINUE
     indx = knkl1 + 3
     DO  i = 1, mdep
         IF (genre .AND. itype == 0) GO TO 3390
         z(mu+i) = z(indx + i*3)
3390     krow    = z(indx + i*3)
         mcode(2)= krow
         IF (krow > mask15) n23 = 3
         DO  j = 1,6
             k    = k+1
             kcol = z(knkl1+j-1)
             mcode(1) = kcol
             IF (kcol > mask15) n23 = 3
             IF (genre .AND. itype == 0) GO TO 3440
             rz(ig+k) = -rz(ig+k)
             IF (genre .AND. itype == 1) GO TO 3400
             CALL WRITE (rgt,mcode,2,0)
             CALL WRITE (rgt,rz(ig+k),1,0)
             CYCLE
3400         ic(j) = ib(j)
             IF (ic(j) > mask15) n23 = 3
             CYCLE
3440         IF (INDEX  >= nind) GO TO 3460
             IF (indcmp(j) /= j) GO TO 3460
             INDEX = INDEX + 1
             ib(INDEX) = kcol
3460         a(6*ilast+j) = rz(ig+k)
         END DO
         IF (   .NOT.genre) GO TO 3635
         IF (itype == -1) GO TO 3635
         IF (itype ==  1) GO TO 3625
         INDEX = INDEX + 1
         ib(INDEX) = krow
         ilast = ilast + 1
         CYCLE
3625     CALL gmmats (rz(ig+k-5),1,6,0,a,6,6,0,b)
         DO  j = 1, 6
             CALL WRITE (rgt,ic(j),1,0)
             CALL WRITE (rgt,krow ,1,0)
             CALL WRITE (rgt,b(j) ,1,0)
         END DO
3635     mcode(1) = krow
         coeff = 1.0
         CALL WRITE (rgt,mcode,2,0)
         CALL WRITE (rgt,coeff,1,0)
     END DO
     IF (.NOT.genre .OR. itype /= 0) GO TO 3645
     !     NO NEED TO COMPUTE DETERMINANT SINCE IT IS NOT USED SUBSEQUENTLY.
     ising = -1
     CALL invers (6,a,6,b,0,det,ising,c)
 
     !     CHECK TO SEE IF GENERAL RIGID ELEMENTS (CRIGD3, CRBE1, CRBAR, AND
     !     CRTRPLT) ARE PROPERLY DEFINED
 
     IF (ising /= 2) GO TO 3645
     WRITE  (nout,6130) ufm,idr
     nogo = 1
3645 IF (more == 0) GO TO 3650
     IF (.NOT.genre) GO TO 1730
     GO TO 1732
 
3650 IF (genre) CALL CLOSE (scr1,1)
     IF (irdg <= 8) THEN
         SELECT CASE ( irdg )
             CASE (    1)
                 GO TO 500
             CASE (    2)
                 GO TO 600
             CASE (    3)
                 GO TO 4000
             CASE (    4)
                 GO TO 4100
             CASE (    5)
                 GO TO 4200
             CASE (    6)
                 GO TO 4300
             CASE (    7)
                 GO TO 5200
         END SELECT
     END IF
     CALL errtrc ('CRIGGP  ',3655)
 
     !     ******************************************************************
 
     !     CRBAR, CRTRPLT, CRIGD3, AND CRBE1 ELEMENTS ARE PROCESSED HERE.
     !     THE CRBAR AND CRTRPLT HAVE THE SAME DATA FORMAT AS THAT OF THE
     !     GENERAL RIGID ELEMENT CRIGD3.
     !     CRBE1 WAS MADE EXACTLY SAME AS CRIGD3 IN IFS3P ROUTINE.
 
     !     ******************************************************************
 
     !     LOCATE CRBAR DATA ON INPUT FILE
 
4000 FILE = geomp
     imhere = 4000
     IF (debug) WRITE (nout,4400) imhere
     CALL locate (*4100,z(buf1),crbar,flag)
     irdg = 4
     GO TO 5000
 
     !     LOCATE CRTRPLT DATA ON INPUT FILE
 
4100 FILE = geomp
     imhere = 4100
     IF (debug) WRITE (nout,4400) imhere
     CALL locate (*4200,z(buf1),crtrpt,flag)
     irdg = 5
     GO TO 5000
 
     !     LOCATE CRIGD3 DATA ON INPUT FILE
 
4200 CALL locate (*4300,z(buf1),crigd3,flag)
     imhere = 4200
     IF (debug) WRITE (nout,4400) imhere
     irdg = 6
     GO TO 5000
 
     !     LOCATE CRBE1 DATA ON INPUT FILE
 
4300 CALL locate (*5200,z(buf1),crbe1,flag)
     imhere = 4300
     IF (debug) WRITE (nout,4400) imhere
4400 FORMAT ('0  I AM HERE/CRIGGP =',i6)
     irdg  = 7
 
5000 genre = .true.
     more  = 1
     IF (debug) WRITE (nout,1005) irdg
 
     !     OPEN SCR1 FILE TO WRITE
 
     CALL OPEN (*1210,scr1,z(buf4),1)
 
     !     READ ELEMENT ID
 
5010 FILE = geomp
     CALL READ (*1230,*1240,geomp,buf,1,0,flag)
     idr = buf(1)
 
     !     READ INDEPENDENT GRID POINTS AND THEIR COMPONENT NUMBERS
 
     n = 0
     j = knkl1
5020 CALL READ (*1230,*1240,geomp,buf,1,0,flag)
     IF (buf(1) == mset) GO TO 5040
     n = n + 7
     z(j) = buf(1)
     CALL READ (*1230,*1240,geomp,z(j+1),6,0,flag)
     j = j + 7
     GO TO 5020
5040 nind = n/7
 
     !     CHECK TO SEE IF THE NUMBER OF INDEPENDENT GRID POINTS
     !     IS MORE THAN ONE AND SET TYPE FLAG
 
     itype = -1
     IF (nind == 1) GO TO 5050
     itype = 0
     j = knkl1
 
     !     WRITE THE INDEPENDENT GRID POINTS AS A PSEUDO CRIGD2 ELEMENT
 
 
     !     WRITE THE ELEMENT ID
 
     CALL WRITE (scr1,idr,1,0)
 
     !     WRITE THE FIRST INDEPENDENT GRID POINT AND ITS COMPONENT NUMBERS
 
     CALL WRITE (scr1,z(j),7,0)
 
     !     WRITE THE TYPE FLAG
 
     CALL WRITE (scr1,itype,1,0)
 
     !     WRITE THE REMAINING INDEPENDENT GRID POINTS AND THEIR
     !     COMPONENT NUMBERS
 
     j = j + 7
     n = n - 7
     CALL WRITE (scr1,z(j),n,0)
     DO  l =1,7
         buf(l) = -1
     END DO
     buf(2) =  0
     CALL WRITE (scr1,buf,7,0)
     itype = 1
 
     !     WRITE THE FIRST INDEPENDENT GRID POINT AND ALL THE
     !     DEPENDENT GRID POINTS AS A PSEUDO CRIGD2 ELEMENT
 
5050 j = knkl1
 
     !     WRITE THE ELEMENT ID
 
     CALL WRITE (scr1,idr,1,0)
 
     !     WRITE THE FIRST INDEPENDENT GRID POINT AND ITS COMPONENT NUMBERS
 
     CALL WRITE (scr1,z(j),7,0)
 
     !     WRITE THE TYPE FLAG
 
     CALL WRITE (scr1,itype,1,0)
 
     !     PROCESS THE DEPENDENT GRID POINTS AND THEIR COMPONENT NUMBERS
 
5060 CALL READ (*1230,*1240,geomp,buf,7,0,flag)
     IF (buf(1) == -1) GO TO 5070
     CALL WRITE (scr1,buf,7,0)
     GO TO 5060
5070 IF (buf(2) == -1) more = 0
     DO  l = 1,7
         buf(l) = -1
     END DO
     buf(2) =  0
     IF (more == 0) buf(2) = -1
     CALL WRITE (scr1,buf,7,0)
     IF (more == 1) GO TO 5010
     CALL WRITE (scr1,0,0,1)
 
     !     CLOSE SCR1, AND OPEN IT FOR READ
 
     CALL CLOSE (scr1,1)
     CALL OPEN (*1210,scr1,z(buf4),0)
     imhere = 5085
     IF (debug) WRITE (nout,4400) imhere
     GO TO 1732
 
     !     *********************************************************
 
     !     CRBE3 AND CRSPLINE ELEMENTS ARE PROCESSED HERE
 
     !     *********************************************************
 
     !     LOCATE CRBE3 DATA ON INPUT FILE
 
5200 FILE = geomp
     irdg = 8
     CALL locate (*5300,z(buf1),crbe3,flag)
     imhere = 5200
     IF (debug) WRITE (nout,4400) imhere
     GO TO 5400
 
     !     LOCATE CRSPLINE DATA ON INPUT FILE
 
5300 FILE = geomp
     irdg = 9
     CALL locate (*5800,z(buf1),crspli,flag)
     imhere = 530
     IF (debug) WRITE (nout,4400) imhere
 
5400 j = irdg-7
     IF (debug) WRITE (nout,1005) irdg
     IF (iprec == 1) CALL crspls (*5600,j,mu,knkl3+1,z(knkl1),again, n23)
     IF (iprec == 2) CALL crspld (*5600,j,mu,knkl3+1,z(knkl1),again, n23)
     SELECT CASE ( j )
         CASE (    1)
             GO TO 5300
         CASE (    2)
             GO TO 5800
     END SELECT
5600 WRITE  (nout,5610) ufm
5610 FORMAT (a23,' 8, INSUFFICIENT CORE FOR CRBE3 OR CRSPLINE RIGID ',  &
         'ELEMENT COMPUTATION')
     nogo = 1
 
     !     *********************************************************
 
     !     CRIGDR AND CRROD (RIGID ROD ELEMENTS) ARE PROCESSED HERE
     !     (CRROD DATA FORMAT WAS CONVERTED TO CRIGDR FORMAT IN IFS3P)
 
     !     *********************************************************
 
     !     LOCATE CRIGDR AND CRROD DATA ON INPUT FILE
 
5800 genre= .false.
     nwds = 4
     FILE = geomp
     CALL locate (*5900,z(buf1),crigdr,flag)
     irdg = 10
     imhere = 5800
     IF (debug) WRITE (nout,4400) imhere
     GO TO 6000
5900 FILE = geomp
     irdg = 11
     CALL locate (*7000,z(buf1),crrod,flag)
     imhere = 5900
     IF (debug) WRITE (nout,4400) imhere
 
 !     ***************************************************************
 
 !                  OPEN CORE FORMAT FOR RIGID ROD
 
 !      NUMBER OF WORDS                 CONTENTS
 
 !           NCSTM    ***     COORDINATE SYSTEM TRANSFORMATION TABLE
 !           NBGPDT   ***     BASIC GRIP POINT DEFINITION TABLE
 !                    ***     SIL 1 ***
 !             3      *       SIL 2   * INDEPENDENT GRID POINT
 !                    ***     SIL 3 ***
 !                    ***     INDEPENDENT GRID POINT BGPDT TABLE
 !             4      *          WORD 1     COORDINATE SYSTEM ID-INTEGER
 !                    ***        WORD 2-4   X, Y, Z, IN BASIC SYSTEM-REAL
 !                    ***     SIL 1 ***
 !             3      *       SIL 2   * DEPENDENT GRID POINT
 !                    ***     SIL 3 ***
 !             4      ***     DEPENDENT GRID POINT BGPDT TABLE
 !                    ***                   ***
 !                    *                       *
 !                    *  AVAILABLE OPEN CORE  *
 !                    *                       *
 !                    *                       *
 !                    ***                   ***
 !                    ***
 !           MDEP     *       DEPENDENT SILS
 !                    ***
 !                    ***
 !          BUFFERS   *       GINO BUFFERS
 !                    ***
 
 !     **************************************************************
 
 
 !     CHECK AVAILABILITY OF CORE
 
6000 CONTINUE
     IF (debug) WRITE (nout,1005) irdg
     itest  = knkl1 + 14 + 27*iprec + 2
     IF (itest >= mu) GO TO 1250
 
     !     READ ELEMENT DATA
 
6010 CALL READ (*1230,*7000,geomp,buf,nwds,0,flag)
     idr   = buf(1)
     idepgp= buf(3)
     icomp = buf(4)
 
     !     PROCESS THE INDEPENDENT GRID POINT
 
     FILE = bgpdt
     j = knkl1
     gpoint = buf(2)
     ASSIGN 6020 TO ret
     ASSIGN 6050 TO ret1
     GO TO 7060
 
     !     STORE SIL VALUES
 
6020 IF (nogo  == 0) GO TO 6030
     IF (j == knkl1) GO TO 6050
     GO TO 6010
6030 z(j  ) = gpoint
     z(j+1) = gpoint + 1
     z(j+2) = gpoint + 2
     kinew  = k - 2*kn
 
     !     LOCATE DATA IN BGPDT
 
     IF (again) GO TO 6035
     ki4 = knkl3 + kinew*4
     IF (ki4 > knkl4) GO TO 1290
     z(j+3) = z(ki4-3)
     z(j+4) = z(ki4-2)
     z(j+5) = z(ki4-1)
     z(j+6) = z(ki4  )
     GO TO 6045
6035 IF (kinew > kiold) GO TO 6040
     CALL bckrec (bgpdt)
     kiold = 0
6040 ki4 = (kinew-kiold-1) * 4
     IF (ki4 > 0) CALL READ (*1230,*1240,bgpdt,buf,-ki4,0,flag)
     CALL READ (*1230,*1240,bgpdt,buf,4,0,flag)
 
     !     STORE BASIC GRID POINT DATA
 
     z(j+3) = buf(1)
     z(j+4) = buf(2)
     z(j+5) = buf(3)
     z(j+6) = buf(4)
6045 kiold  = kinew
     IF (j /= knkl1) GO TO 6060
 
     !     PROCESS THE DEPENDENT GRID POINT
 
6050 j = j + 7
     gpoint = idepgp
     ASSIGN 6010 TO ret1
     GO TO 7060
6060 IF (iprec == 1) CALL crdrd  (*6065,*6125,mu,icomp,n23)
     IF (iprec == 2) CALL crdrd2 (*6065,*6125,mu,icomp,n23)
     GO TO 6010
6065 WRITE  (nout,6070) ufm,idr
6070 FORMAT (a23,' 3133, RIGID ELEMENT',i9,' HAS ZERO LENGTH')
     nogo = 1
     GO TO 6010
6125 WRITE  (nout,6130) ufm,idr
6130 FORMAT (a23,' 3134, RIGID ELEMENT',i9,' IS NOT PROPERLY DEFINED')
     nogo = 1
     GO TO 6010
 
7000 IF (irdg == 10) GO TO 5900
 
     IF (again) CALL CLOSE (bgpdt,clsrew)
     IF (nogo /= 0) CALL mesage (-61,0,NAME)
     CALL WRITE (rgt,0,0,1)
 
     !     WRITE A LIST OF DEPENDENT SIL VALUES FOR RIGID ELEMENTS ONTO THE
     !     RGT IN SORTED FORM
 
     jrigid = mu + 1
     m = buf4 - jrigid
     CALL sort  (0,0,1,1,z(jrigid),m)
     CALL WRITE (rgt,z(jrigid),m,1)
     j = buf4-1
     IF (debug) WRITE (nout,7010) (z(i),i=jrigid,j)
7010 FORMAT (/,'  CRIGGP/@7010  DEPEND.SIL LIST:',/,(5X,10I7))
     knkl1 = knkl2
 
     !     CLOSE RGT FILE AND RETURN
 
     CALL CLOSE (rgt,clsrew)
     RETURN
 
     !     **********************************************************
 
     !     INTERNAL SUBROUTINE TO PERFORM BINARY SEARCH IN EQEXIN
     !     AND CONVERT THE EXTERNAL NUMBER TO A SIL VALUE
 
7060 klo = 0
     khi = kn2
     lastk = 0
7070 k= (klo+khi+1)/2
     IF (lastk == k) GO TO 1350
     lastk = k
     IF (gpoint-z(2*k-1) < 0.0) THEN
         GO TO  7090
     ELSE IF (gpoint-z(2*k-1) == 0.0) THEN
         GO TO  7150
     ELSE
         GO TO  7100
     END IF
7090 khi= k
     GO TO 7070
7100 klo= k
     GO TO 7070
7150 k = z(2*k) + 2*kn
     gpoint= z(k)
     GO TO ret, (1740,1750,6020)
 
     !     **********************************************************
 
     !     FATAL ERROR MESSAGES
 
1210 j= -1
     GO TO 1260
1230 j= -2
     GO TO 1260
1240 j= -3
     GO TO 1260
1250 IF (again) GO TO 1280
     CALL CLOSE (scr1,clsrew)
     WRITE  (nout,1255) imhere
1255 FORMAT (///,' *** CRIGGP/GP4 NEEDS MORE OPEN CORE.',  &
         /5X,' CRIGGP REVERTED TO USE SLOW METHOD',i9,//)
     GO TO 290
1260 CALL mesage (j,FILE,NAME)
1280 j= -8
     GO TO 1260
1290 WRITE  (nout,1300) knkl1,knkl3,knkl4,ki4
1300 FORMAT (//,' *** SYSTEM FATAL ERROR IN CRIGGP',4I10)
     j =-61
     GO TO 1260
1330 nogo= 1
     CALL mesage (30,n,buf)
     GO TO ret1, (1743,1745,6010,6050)
1350 IF (genre .AND. itype == 1 .AND. ntype == 1) GO TO 1743
     buf(1) = gpoint
     buf(2) = irdg*100000000 + idr
     n = 151
     GO TO 1330

 END SUBROUTINE criggp
