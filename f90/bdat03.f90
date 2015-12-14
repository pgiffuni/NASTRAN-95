SUBROUTINE bdat03
     
    !     THIS SUBROUTINE PROCESSES TRANS BULK DATA, GENERATES THE
    !     TRANSFORMATION MATRIX, AND WRITES TO SCBDAT.
 
    EXTERNAL        rshift,andf
    LOGICAL :: tdat
    INTEGER :: buf1,trans(2),geom4,combo,aaa(2),outt,buf2,buf4,  &
        andf,rshift,ihd(10),z
    DIMENSION       temp(9),xax(3),yax(3),zax(3),v2(3),out(9)
    CHARACTER (LEN=23) :: ufm
    COMMON /xmssg / ufm
    COMMON /cmb001/ scr1,scr2,scbdat,scsfil,scconn,scmcon,sctoc, geom4,casecc
    COMMON /zzzzzz/ z(1)
    COMMON /cmb002/ buf1,buf2,buf3,buf4,buf5,score,lcore,intp,outt
    COMMON /cmb003/ combo(7,5),conset,iauto,toler,npsub,conect,tran,  &
        mcon,restct(7,7),isort,origin(7,3),iprint
    COMMON /output/ ititl(96),ihead(96)
    COMMON /system/ xxx,iot,junk(6),nlpp,junk1(2),line,junk2(2), idat(3)
    COMMON /cmb004/ tdat(6)
    COMMON /BLANK / step,idry
    DATA    ihd   / 4H  su , 4HMMAR , 4HY of , 4H pro  , 4HCESS    ,  &
        4HED t , 4HRANS , 4H bul , 4HK da  , 4HTA      /
    DATA    trans / 310,3 /,    aaa / 4HBDAT,4H03   /  , iblnk / 4H    /
 
    ngtrn = z(buf4)
    inum  = 1
    ierr  = 0
    DO  i = 1,7
        DO  j = 1,3
            origin(i,j) = 0.0
        END DO
    END DO
    DO  i = 1,96
        ihead(i) = iblnk
    END DO
    j = 1
    DO  i = 76,85
        ihead(i) = ihd(j)
        j = j + 1
    END DO
    CALL locate (*220,z(buf1),trans(1),flag)
    ifile = geom4
40  CALL READ (*300,*130,geom4,id,1,0,n)
    DO  i = 1,npsub
        it = 1
        IF (id == combo(i,3)) GO TO 80
    END DO
    IF (ngtrn == 0) GO TO 70
    DO  i = 1,ngtrn
        it = 2
        IF (id == z(buf4+i)) GO TO 80
    END DO
70 CONTINUE
   CALL READ (*300,*310,geom4,temp,-9,0,nnn)
   GO TO 40
80 tdat(3) = .true.
   IF (it == 1) combo(i,3) = -combo(i,3)
   IF (it == 2)  z(buf4+i) = -z(buf4+i)
   CALL READ (*300,*310,geom4,temp,9,0,nnn)
   IF (it /= 1) GO TO 100
   DO  ll = 1,3
       origin(i,ll) = temp(ll)
   END DO
100 CONTINUE
 
    !     DEFINE Z-AXIS
 
    zax(1) = temp(4) - temp(1)
    zax(2) = temp(5) - temp(2)
    zax(3) = temp(6) - temp(3)
 
    !     DEFINE Y-AXIS
 
    v2(1)  = temp(7) - temp(1)
    v2(2)  = temp(8) - temp(2)
    v2(3)  = temp(9) - temp(3)
    yax(1) = zax(2)*v2(3) - zax(3)*v2(2)
    yax(2) = zax(3)*v2(1) - zax(1)*v2(3)
    yax(3) = zax(1)*v2(2) - zax(2)*v2(1)
 
    !     DEFINE X-AXIS
 
    xax(1) = yax(2)*zax(3) - zax(2)*yax(3)
    xax(2) = yax(3)*zax(1) - zax(3)*yax(1)
    xax(3) = yax(1)*zax(2) - zax(1)*yax(2)
 
    !     CHANGE TO UNIT VECTORS
 
    zmag = SQRT(zax(1)**2 + zax(2)**2 + zax(3)**2)
    ymag = SQRT(yax(1)**2 + yax(2)**2 + yax(3)**2)
    xmag = SQRT(xax(1)**2 + xax(2)**2 + xax(3)**2)
    DO  i = 1,3
        zax(i) = zax(i)/zmag
        yax(i) = yax(i)/ymag
        xax(i) = xax(i)/xmag
    END DO
    CALL WRITE (scbdat,id,1,0)
    CALL WRITE (scbdat, 1,1,0)
    CALL WRITE (scbdat,temp(1),3,0)
    out(1) = xax(1)
    out(2) = yax(1)
    out(3) = zax(1)
    out(4) = xax(2)
    out(5) = yax(2)
    out(6) = zax(2)
    out(7) = xax(3)
    out(8) = yax(3)
    out(9) = zax(3)
    CALL WRITE (scbdat,out,9,0)
    IF (andf(rshift(iprint,6),1) /= 1) GO TO 120
    inum = inum + 1
    IF (MOD(inum,2) == 0) CALL page
    WRITE (outt,430) id
    WRITE (outt,440) (temp(i),i=1,3)
    WRITE (outt,420) ( out(i),i=1,9)
120 CONTINUE
    GO TO 40
130 CONTINUE
 
    !     PROCESS REPEATED GTRAN IDS
 
    IF (ngtrn < 2) GO TO 160
    ngtrn1 = ngtrn - 1
    DO  i = 1,ngtrn1
        IF (z(buf4+i) >= 0) CYCLE
        kk = i + 1
        DO  j = kk,ngtrn
            IF (IABS(z(buf4+i)) == z(buf4+j)) z(buf4+j) = -z(buf4+j)
        END DO
    END DO
160 npm1 = npsub - 1
    DO  i = 1,npm1
        IF (combo(i,3) >= 0) CYCLE
        kk = i + 1
        DO  j = kk,npsub
            IF (IABS(combo(i,3)) /= combo(j,3)) CYCLE
            combo(j,3) = -combo(j,3)
            DO  jdh = 1,3
                origin(j,jdh) = origin(i,jdh)
            END DO
        END DO
    END DO
 
    !     TEST TO SEE THAT ALL TRANS HAVE BEEN FOUND
 
    DO  i = 1,npsub
        IF (combo(i,3) <= 0) CYCLE
        ierr = 1
        WRITE (outt,400) ufm,combo(i,3)
    END DO
    IF (ngtrn == 0) GO TO 220
    DO  i = 1,ngtrn
        IF (z(buf4+i) <= 0) CYCLE
        ierr = 1
        WRITE (outt,410) ufm,z(buf4+i)
    END DO
220 CALL eof (scbdat)
    CALL WRITE (scbdat,id,1,1)
    CALL CLOSE (scbdat,1)
    DO  i = 1,npsub
        combo(i,3) = IABS(combo(i,3))
    END DO
    IF (ierr == 1) idry = -2
    RETURN
 
300 imsg = -2
    GO TO 320
310 imsg = -3
320 CALL mesage (imsg,ifile,aaa)
    RETURN
 
400 FORMAT (a23,' 6511, THE REQUESTED TRANS SET ID',i9,  &
        ' HAS NOT BEEN DEFINED BY BULK DATA.')
410 FORMAT (a23,' 6513, THE TRANS SET ID',i9,' REQUESTED BY A GTRAN ',  &
        'BULK DATA CARD HAS NOT BEEN DEFINED.')
420 FORMAT (43X,5H*****,42X,5H*****, /3(43X,1H*,50X,1H*, /43X,1H*,1X,  &
        3E15.6,4X,1H*,/),43X,1H*,50X,1H*, /43X,5H*****,42X,5H*****)
430 FORMAT (//48X,34HTRANS set identification NUMBER = ,i8)
440 FORMAT ( /50X,37HCOORDINATES of origin in basic system ,  &
        /45X,3E15.6, //58X,21HTRANSFORMATION matrix/)

END SUBROUTINE bdat03
