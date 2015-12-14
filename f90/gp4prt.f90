SUBROUTINE gp4prt (ibuf)
     
!     1. PRINTS  DOF VS. DISP. SETS IF DIAG 21 ON.
!     2. PRINTS  DISP. SETS VS. DOF IF DIAG 22 ON.
!     3. CREATES SUBSTRUCTURE COUPLING DATA TABLE.
 
!     IF IBUF IS .LT. 0, SOME FILES MAY NOT BE CLOSED PROPERLY WHEN
!     THIS ROUTINE IS CALLED
 
 
 INTEGER, INTENT(IN OUT)                  :: ibuf
 EXTERNAL        andf,orf
 INTEGER :: title(2,8),z,sysbuf,andf,eqexin,orf,FILE,NAME(3),  &
     d21,d22,msk(12),iflg(8),scr1,zdum(10),zcom(10),  &
     erec1,buf,two,um,uo,ur,usg,usb,ul,ua,uf,us,un,ug,  &
     dash,upbit(12),BLANK,tdb204,nam204(2),trl(7), exflag,extype,sbit(12),  &
     ifrmat(32),iifrmt(2),iafrmt(2)
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm,uim
 COMMON /system/ sysbuf,nout,junk(6),nlpp,mtemp,npage,line
 COMMON /two   / two(32)
 COMMON /bitpos/ um,uo,ur,usg,usb,ul,ua,uf,us,un,ug
 COMMON /zzzzzz/ z(1)
 COMMON /BLANK / luset,mpcf1,mpcf2,single,omit1,react,nskip,  &
     repeat,nosets,nol,noa,idsub,iautsp
 DATA    scr1  / 301/  , eqexin / 103   /, tdb204 / 204 /
 DATA    NAME  / 4HGP4P, 4HRT   , 4H    /
DATA    ibegn , iend  / 4HBEGN , 4HEND /
DATA    title / 4H    , 4H mpc , 4H    , 4H spc ,  &
    4H    , 4HOMIT , 4HANAL, 4HYSIS ,  &
    4H  su, 4HPORT , 4HPERM, 4H spc ,  &
    4HBDRY, 4H spc , 4HAUTO, 4H spc /
DATA    BLANK / 1H  /
DATA    dash  / 1H- /
DATA    ifrmat/ 4H(13X, 4H,i6,, 4H3X,i, 4H8,1X, 4H,a1,, 4HI2,1,  &
    4HX   , 4H,1X,, 4H  i6, 4H,1X,, 4H  i6, 4H,1X,,  &
    4H  i6, 4H,1X,, 4H  i6, 4H,1X,, 4H  i6, 4H,1X,,  &
    4H  i6, 4H,1X,, 4H  i6, 4H,1X,, 4H  i6, 4H,1X,, 4H  i6, 4H,1X,,  &
    4H  i6, 4H,1X,, 4H  i6, 4H,1X,, 4H  i6, 4H)   /
DATA    iifrmt/ 4H,1X,, 4H  i6/
DATA    iafrmt/ 4H,3X,, 4H  a4/
DATA    iprint/ 0   /


CALL sswtch (21,d21)
CALL sswtch (22,d22)
IF (d21 /= 1 .AND. d22 /= 1 .AND. idsub <= 0) RETURN
IF (iprint == 1) RETURN
iprint  = 1
NAME(3) = ibegn
CALL conmsg (NAME,3,0)
buf  = IABS(ibuf)
FILE = eqexin
IF (ibuf < 0) CALL CLOSE (eqexin,1)
CALL OPEN (*1220,eqexin,z(buf),0)
CALL fwdrec (*1230,eqexin)
CALL fwdrec (*1230,eqexin)
erec1 = luset +1
CALL READ (*1230,*2900,eqexin,z(erec1),buf-erec1,1,kn)
GO TO 9001
2900 CALL CLOSE (eqexin,1)
CALL sort (0,0,2,2,z(erec1),kn)

IF (d21 /= 1) GO TO 3000
ku = 1
msk(ku+1 ) = two(usb)
msk(ku+2 ) = two(usg)
msk(ku+3 ) = two( ul)
msk(ku+4 ) = two( ua)
msk(ku+5 ) = two( uf)
msk(ku+6 ) = two( un)
msk(ku+7 ) = two( ug)
msk(ku+8 ) = two( ur)
msk(ku+9 ) = two( uo)
msk(ku+10) = two( us)
msk(ku+11) = two( um)
DO  ku = 1,12
  sbit(ku) = 0
END DO
CALL page1
line = line + 2
WRITE (nout,1900) uim
line = line + 4
WRITE (nout,1902)
i  = erec1
kl = 0
DO  k = 1,kn,2
  itm = z(k+i)/10
  itm = z(k+i) - 10*itm
  l   = 6
  IF (itm == 2) l = 1
  DO  kk = 1,l
    kl = kl + 1
    iu = z(kl)
    ip = z(i+k-1)
    idof = kk
    IF (andf(msk(11),iu) == 0) GO TO 2914
    IF (andf(msk(2),iu) /= 0 .OR. andf(msk(3),iu) /= 0) GO TO 2914
    sbit(1)   = sbit(1) + 1
    upbit(1)  = sbit(1)
    ifrmat(8) = iifrmt(1)
    ifrmat(9) = iifrmt(2)
    GO TO 2916
    2914 upbit(1)  = BLANK
    ifrmat(8) = iafrmt(1)
    ifrmat(9) = iafrmt(2)
    2916 DO  ku = 2,12
      INDEX = 2*(ku-1) + 8
      IF (andf(msk(ku),iu) == msk(ku)) GO TO 2920
      upbit(ku) = BLANK
      ifrmat(INDEX  ) = iafrmt(1)
      ifrmat(INDEX+1) = iafrmt(2)
      CYCLE
      2920 sbit (ku) = sbit(ku) + 1
      upbit(ku) = sbit(ku)
      ifrmat(INDEX  ) = iifrmt(1)
      ifrmat(INDEX+1) = iifrmt(2)
    END DO
    IF (l == 1) idof = 0
    line = line + 1
    IF (line <= nlpp) GO TO 2945
    CALL page1
    WRITE (nout,1902)
    line = line + 5
    2945 WRITE (nout,ifrmat) kl,ip,dash,idof,upbit
  END DO
END DO
WRITE (nout,1901) sbit
line = line + 2

3000 IF (d22 /= 1 .AND. idsub <= 0) RETURN
msk(1)  = two(um )
msk(2)  = two(us )
msk(3)  = two(uo )
msk(4)  = two(ua )
msk(5)  = two(ur )
msk(6)  = two(usg)
msk(7)  = two(usb)
exflag  = 0
extype  = 0
IF (d22 /= 1) GO TO 3010
CALL page1
line = line + 2
WRITE (nout,1907) uim
line = line + 4
3010 FILE = scr1
IF (ibuf < 0) CALL CLOSE (scr1,1)
CALL OPEN (*1220,scr1,z(buf),1)
DO  imk = 1,8
  iflg(imk) = 0
  i  = erec1
  ip = 0
  kl = 0
  DO  k = 1,kn,2
    itm = z(k+i)/10
    itm = z(k+i) - 10*itm
    l   = 6
    IF (itm == 2) l = 1
    DO  kk = 1,l
      kl = kl + 1
      iu = z(kl)
      IF (z(i+k-1) < ip) exflag = 1
      ip = z(i+k-1)
      IF (l == 1) GO TO 3920
      idof = kk
      extype = orf(extype,2)
      GO TO 3930
      3920 idof = 0
      extype = orf(extype,1)
      3930 IF (imk /= 8) GO TO 3940
      IF (andf(iu,msk(2)) == 0) CYCLE
      IF (andf(iu,msk(6)) /= 0 .OR. andf(iu,msk(7)) /= 0) CYCLE
      GO TO 3945
      3940 IF (andf(iu,msk(imk)) /= msk(imk)) CYCLE
      3945 CALL WRITE (scr1,10*ip+idof,1,0)
      iflg(imk) = 1
    END DO
  END DO
  IF (iflg(imk) /= 1) CYCLE
  CALL WRITE (scr1,0,0,1)
END DO
CALL WRITE (scr1,z(1),luset,1)
CALL CLOSE (scr1,1)
CALL OPEN  (*1220,scr1,z(buf),0)
iflag = 0
DO  i = 1,8
  IF (iflg(i) /= 1) CYCLE
  iflag = iflag + 1
  CALL READ (*1230,*4010,scr1,z(1),buf,1,kn)
  CALL page2 (-4)
  WRITE (nout,9501) FILE
  GO TO 4600
  4010 CONTINUE
  IF (idsub <= 0 .OR. i /= 4) GO TO 4040
  CALL CLOSE (scr1,2)
  FILE = tdb204
  CALL OPEN  (*1220,tdb204,z(buf),1)
  CALL fname (tdb204,nam204)
  CALL WRITE (tdb204,nam204,2,1)
  CALL WRITE (tdb204,z(1),kn,1)
  CALL CLOSE (tdb204,1)
  trl(1) = tdb204
  trl(2) = 0
  trl(3) = kn
  trl(4) = 0
  trl(5) = idsub
  trl(6) = exflag
  trl(7) = extype
  CALL wrttrl (trl)
  CALL OPEN (*1220,scr1,z(buf),2)
  4040 CONTINUE
  IF (d22 /= 1) CYCLE
  ipas = kn/10
  irem = kn - 10*ipas
  IF (iflag > 1) line = nlpp
  id1  =-9
  inos = 0
  IF (ipas < 1) GO TO 4105
  DO  k = 1,ipas
    DO  j = 1,10
      inos = inos + 1
      zdum(j) = z(inos)/10
      zcom(j) = z(inos) - 10*zdum(j)
    END DO
    line = line + 1
    IF (iflag == 1 .AND. k == 1) GO TO 4060
    IF (line <= nlpp) GO TO 4090
    CALL page1
    4060 WRITE (nout,1910) title(1,i),title(2,i)
    line = line + 5
    4090 CONTINUE
    id1 = id1 + 10
    WRITE (nout,1913) id1,(zdum(kk),zcom(kk),kk=1,10)
  END DO
  4105 IF (irem == 0) CYCLE
  DO  j = 1,irem
    inos = inos + 1
    zdum(j) = z(inos)/10
    zcom(j) = z(inos) - zdum(j)*10
  END DO
  line = line + 1
  IF (iflag == 1 .AND. ipas == 0) GO TO 4120
  IF (line <= nlpp) GO TO 4400
  CALL page1
  4120 WRITE (nout,1910) title(1,i),title(2,i)
  line = line + 5
  4400 CONTINUE
  id1 = id1 + 10
  WRITE (nout,1913) id1,(zdum(kk),zcom(kk),kk=1,irem)
END DO

!     RE-ESTABLISH USET IN OPEN CORE.

4600 CALL READ (*1230,*9001,scr1,z(1),luset,1,kn)
CALL CLOSE (scr1,1)
NAME(3) = iend
CALL conmsg (NAME,3,0)

!     TERMINATE RUN IF DIAG 21 OR 22, AND DIAG 20 ARE REQUESTED BY UESER
!     SIMLUTANEOUSLY

CALL sswtch (20,j)
IF (j == 0 .OR. d21+d22 == 0) RETURN
WRITE  (nout,4700)
4700 FORMAT (10X,25HJOB terminated by diag 20)
CALL pexit

1220 j = -1
GO TO 1260
1230 j = -2
1260 CALL mesage (j,FILE,NAME)
RETURN

1900 FORMAT (a29,' 2118, SUBROUTINE GP4PRT - DIAG 21 SET-DOF VS. DISP',  &
    ' SETS FOLLOWS.')
1901 FORMAT (1H0, 34H--- c o l u m n   t o t a l s --- , 12I7)
1902 FORMAT (1H0,14X,5H(sil), /14X,  &
    48HINT dof   ext gp. dof  sauto     sb     sg      ,  &
    49HL      a      f      n      g      r      o      ,  &
    8HS      m, /1H , 131(1H-))
1907 FORMAT (a29,' 2119, SUBROUTINE GP4PRT - DIAG 22 SET DISP SETS VS',  &
    '. DOF FOLLOWS')
1910 FORMAT (1H0,52X,2A4,17H displacement set ,/  &
    1H0,15X,3H-1-,8X,3H-2-,8X,3H-3-,8X,3H-4-,8X,3H-5-,  &
    8X,3H-6-,8X,3H-7-,8X,3H-8-,8X,3H-9-,7X,4H-10- ,/1H )
1913 FORMAT (1H ,i6,1H=,10(1X,i8,1H-,i1))

!     ERRORS

9001 CALL page2(-4)
WRITE  (nout,9501) uwm,FILE
9501 FORMAT (a25,' 2110, INSUFFICIENT CORE TO HOLD CONTENTS OF GINO ',  &
    'FILE',i4, //5X, 'FURTHER PROCESSING OF THIS DATA BLOCK IS ABANDONED.')

CALL CLOSE (FILE,1)
RETURN

END SUBROUTINE gp4prt
