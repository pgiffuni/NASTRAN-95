SUBROUTINE biotsv (xx,yy,zz,hcx,hcy,hcz)
     
    !     THIS ROUTINE COMPUTES THE MAGNETIC FIELD AT A POINT (XX,YY,ZZ)
    !     DUE TO MAGNETIC SOIRCES. THE ROUTINE IS USED BY PROLATE IN
    !     COMPUTING HC POTENTIALS USING LINE INTEGRALS. AT Z(IST) IS STORED
    !     LOAD INFO. NEEDED FOR THIS SUBCASE (WHICH COULD BE A LOAD
    !     COMBINATION) AS STORED BY ROUTINE LOADSU. THE INFO. IS STORED AS
    !     FOLLOWS -
 
    !     OVERALL SCALE FACTOR - ALLS
    !     NUMBER OF SIMPLE LOADS - NSIMP
    !     SCALE FACTOR FOR 1ST SIMPLE LOAD
    !     NUMBER OF LOAD CARDS FOR 1ST SIMPLE LOAD
    !     SCALE FACTOR FOR 2ND SIMPLE LOAD
    !     NUMBER OF LOAD CARDS FOR 2ND SIMPLE LOAD
    !      .
    !     ETC.
    !      .
    !     TYPE(NOBLD) OF 1ST CARD FOR 1ST SIMPLE LOAD
    !     NUMBER OF CARDS FOR THIS TYPE - IDO
    !     LOAD INFO FOR THIS TYPE FOR 1ST SIMPLE LOAD
    !     ANOTHER TYPE FOR 1ST SIMPLE LOAD
    !      .
    !     ETC
    !      .
    !     LOAD CARDS FOR SUBSEQUENT SIMPLE LOADS FOR THIS SUBCASE
 
 
    REAL, INTENT(IN OUT)                     :: xx
    REAL, INTENT(IN OUT)                     :: yy
    REAL, INTENT(IN OUT)                     :: zz
    REAL, INTENT(OUT)                        :: hcx
    REAL, INTENT(OUT)                        :: hcy
    REAL, INTENT(OUT)                        :: hcz
    INTEGER :: hest,bgpdt,scr1,FILE,buf2,subcas
    DIMENSION       nam(2),iz(1),buf(50),ibuf(50),mcb(7),hc(3),hc1(3), hc2(3)
    COMMON /biot  / ng1,ng2,ist,subcas,x1,y1,z1,x2,y2,z2,buf2,remfl,  &
        mcore,load,nslt,scr1,hest,ntot
    COMMON /system/ sysbuf,iout
    COMMON /zzzzzz/ z(1)
    EQUIVALENCE     (z(1),iz(1)),(buf(1),ibuf(1))
    DATA    nam   / 4HBIOT,4HSV  /
 
    hcx    = 0.
    hcy    = 0.
    hcz    = 0.
    scr1   = 301
    bgpdt  = 103
    mcb(1) = bgpdt
    CALL rdtrl (mcb)
    nrowsp = mcb(2)
    mcb(1) = scr1
    CALL rdtrl (mcb)
    n3     = mcb(3)
    ngrids = n3/3
 
    alls   = z(ist+1)
    nsimp  = iz(ist+2)
    isimp  = ist + 2*nsimp + 2
 
    !     LOOP ON NUMBER OF SIMPLE LOADS
 
    DO  ns = 1,nsimp
        nc     = 0
        hc(1)  = 0.
        hc(2)  = 0.
        hc(3)  = 0.
   
        factor = z(ist+2*ns+1)
        ncards = iz(ist+2*ns+2)
15      nobld  = iz(isimp+1)
        ido    = iz(isimp+2)
        isimp  = isimp + 2
   
   
        ktype  = nobld - 19
        SELECT CASE ( ktype )
            CASE (    1)
                GO TO 20
            CASE (    2)
                GO TO 30
            CASE (    3)
                GO TO 40
            CASE (    4)
                GO TO 50
            CASE (    5)
                GO TO 60
        END SELECT
20      mwords = 3*nrowsp
        GO TO 70
30      mwords = 12
        GO TO 70
40      mwords = 48
        GO TO 70
50      mwords = 9
        GO TO 70
60      mwords = 0
   
        70 DO  j = 1,ido
     
            SELECT CASE ( ktype )
                CASE (    1)
                    GO TO 145
                CASE (    2)
                    GO TO 150
                CASE (    3)
                    GO TO 150
                CASE (    4)
                    GO TO 150
                CASE (    5)
                    GO TO 180
            END SELECT
     
        !     SPCFLD DATA STARTS AT Z(ISIMP+1)
     
145     CONTINUE
     
     
        !     NG1 AND NG2 ARE THE SIL NUMBERS OF THE END POINTS OF THE LINE
        !     INTEGRL WITH (X1,Y1,Z1) AND (X2,Y2,Z2) BEING THE COORDINATES.
        !     LINEARLY INTERPOLATE TO (XX,YY,ZZ). THE SILS ARE POINTERS INTO
        !     THE SPCFLD DATA
     
        isub   = isimp + 3*ng1
        hc1(1) = z(isub-2)
        hc1(2) = z(isub-1)
        hc1(3) = z(isub)
        isub   = isimp + 3*ng2
        hc2(1) = z(isub-2)
        hc2(2) = z(isub-1)
        hc2(3) = z(isub)
148     tlen   = SQRT((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)
        xlen   = SQRT((xx-x1)**2 + (yy-y1)**2 + (zz-z1)**2)
        ratio  = xlen/tlen
        hc(1)  = hc(1) + (1.-ratio)*hc1(1) + ratio*hc2(1)
        hc(2)  = hc(2) + (1.-ratio)*hc1(2) + ratio*hc2(2)
        hc(3)  = hc(3) + (1.-ratio)*hc1(3) + ratio*hc2(3)
        GO TO 240
     
        !     CEMLOOP,GEMLOOP,MDIPOLE
     
        150 DO  k = 1,mwords
            buf(k) = z(isimp+k)
        END DO
        ltype  = ktype - 1
        SELECT CASE ( ltype )
            CASE (    1)
                GO TO 155
            CASE (    2)
                GO TO 160
            CASE (    3)
                GO TO 165
        END SELECT
155     CALL axloop (buf,ibuf,xx,yy,zz,hca,hcb,hcc)
        GO TO 170
160     CALL geloop (buf,ibuf,xx,yy,zz,hca,hcb,hcc)
        GO TO 170
165     CALL dipole (buf,ibuf,xx,yy,zz,hca,hcb,hcc)
     
170     hc(1)  = hc(1) + hca
        hc(2)  = hc(2) + hcb
        hc(3)  = hc(3) + hcc
        GO TO 240
     
        !     REMFLUX - BRING IN VALUES FROM SCR1 AFTER POSITIONING TO PROPER
        !     CASE
     
180     CALL gopen (scr1,z(buf2),0)
        ic = subcas - 1
        IF (ic == 0) GO TO 200
        DO  i = 1,ic
            CALL fwdrec (*520,scr1)
        END DO
     
200     isimp1 = 6*ngrids + ntot
        CALL fread (scr1,z(isimp1+1),n3,1)
     
        CALL CLOSE (scr1,1)
     
        !     MUST MATCH NG1 AND NG2 TO SIL-S IN CORE TO LOCATE REMFLUX INFO ON
        !    SCR1
     
        ing1 = 0
        ing2 = 0
        DO  i = 1,ngrids
            IF (ng1 == iz(i)) GO TO 205
            IF (ng2 == iz(i)) GO TO 210
            CYCLE
205         ing1 = i
            IF (ing2 == 0) CYCLE
            GO TO 230
210         ing2 = i
            IF (ing1 == 0) CYCLE
            GO TO 230
        END DO
        GO TO 510
230     isub   = 3*ing1 + isimp1
        hc1(1) = z(isub-2)
        hc1(2) = z(isub-1)
        hc1(3) = z(isub)
        isub   = 3*ing2 + isimp1
        hc2(1) = z(isub-2)
        hc2(2) = z(isub-1)
        hc2(3) = z(isub)
     
        !     INTERPOLATE AS WITH SPCFLD
     
        GO TO 148
     
        !     DONE FOR ONE CARD OF PRESENT TYPE  - GET ANOTHER
     
240     isimp = isimp + mwords
        nc = nc + 1
     
    END DO
   
    !     CHECK TO SEE IF WE ARE DONE WITH THIS LOAD FACTOR
   
    IF (nc < ncards) GO TO 15
   
    !     DONE WITH THIS SIMPLE LOAD. APPLY INDIVIDUAL AND OVERALL SCALE
    !     FACTORS THEN GET ANOTHER SIMPLE LOAD
   
    fac = factor*alls
    hcx = hcx + fac*hc(1)
    hcy = hcy + fac*hc(2)
    hcz = hcz + fac*hc(3)
   
END DO
 
!     DONE
 
RETURN
 
510 WRITE  (iout,511) ng1,ng2
511 FORMAT ('0*** LOGIC ERROR, SILS',2I8,  &
    ' CANNOT BE FOUND IN PROLATE LIST IN BIOTSV')
CALL mesage (-61,0,0)
 
520 CALL mesage (-2,FILE,nam)

RETURN
END SUBROUTINE biotsv
