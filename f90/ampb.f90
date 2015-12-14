SUBROUTINE ampb(phidh,gtka,d1jk,d2jk,d1je,d2je,useta,  &
    djh1,djh2,gki,scr1,scr2,scr3)
     
    !     THE PURPOSE OF THIS SUBROUTINE IS TO SOLVE FOR THE DJH MATRICES.
    !      IT ALSO COMPUTES GKI FOR LATER USE.
    !      THE STEPS ARE,
 
    !     1. PHIDH GOES TO   1       1      1
    !                        1 PHIA  1      1
    !                        1 ----- 1 ---- 1
    !                        1       1      1
    !                        1       1      1
 
    !     2. GKI =GTKA$PHIA
 
    !     3. DJI1=D1JK*GKI
    !     4. DJI2=D2JK*GKI
    !     5.
    !     6. DJH1= 1 DJI1 1 D1JE 1
    !              1      1      1
    !     7. DJH2= 1 DJI2 1 D2JE 1
 
 
 
 
    INTEGER, INTENT(IN)                      :: phidh
    INTEGER, INTENT(IN OUT)                  :: gtka
    INTEGER, INTENT(IN OUT)                  :: d1jk
    INTEGER, INTENT(IN OUT)                  :: d2jk
    INTEGER, INTENT(IN OUT)                  :: d1je
    INTEGER, INTENT(IN OUT)                  :: d2je
    INTEGER, INTENT(IN)                      :: useta
    INTEGER, INTENT(IN)                      :: djh1
    INTEGER, INTENT(IN)                      :: djh2
    INTEGER, INTENT(IN OUT)                  :: gki
    INTEGER, INTENT(IN OUT)                  :: scr1
    INTEGER, INTENT(IN OUT)                  :: scr2
    INTEGER, INTENT(IN)                      :: scr3
    INTEGER :: phia,dji1,dji2,mcb(7),ud,ua,ue
 
    COMMON /BLANK/noue
    COMMON /patx/lc,ns0,ns1,ns2,iuset
    COMMON /bitpos/um,uo,ur,usg,usb,ul,ua,uf,us,un,ug,ue,up,une,ufe,ud  &
        ,ups,usa,uk,upa
    COMMON/zzzzzz/z(1)
    COMMON/system/sysbuf,nout,skip(52),iprec
    COMMON /ampcom/ ncolj
 
    !-----------------------------------------------------------------------
 
    mcb(1)=phidh
    CALL rdtrl(mcb)
    noh=mcb(2)
 
    !     DETERMINE IF PHIDH MUST BE MODIFIED
 
    IF(noue == -1)GO TO 20
 
    !     BUILD PARTITIONING VECTORS
 
    iuset = useta
    lc=korsz(z)
    CALL calcv(scr1,ud,ua,ue,z)
    CALL ampb1(scr2,noh-noue,noue)
 
    !     PERFORM PARTITION
    !                       RP   CP
    CALL ampb2(phidh,scr3,0,0,0,scr2,scr1,0,0)
    phia=scr3
    GO TO 30
 
    !     NO MOD REQUIRED
 
20  phia=phidh
30  CONTINUE
 
   !     COMPUTE GKI
 
    CALL ssg2b(gtka,phia,0,gki,1,iprec,1,scr1)
 
   !     START COMPUTATION OF DJH MATRICES
 
    dji1=scr3
    dji2=scr3
    IF(noue > 0)GO TO 40
    dji1=djh1
    dji2=djh2
40  CONTINUE
    CALL ssg2b(d1jk,gki,0,dji1,1,iprec,1,scr1)
    IF(noue == -1)GO TO 50
    CALL merged(dji1,d1je,0,0,djh1,scr2,0,0,ncolj)
50  CONTINUE
    CALL ssg2b(d2jk,gki,0,dji2,1,iprec,1,scr1)
    IF(noue == -1)GO TO 60
    CALL merged(dji2,d2je,0,0,djh2,scr2,0,0,ncolj)
60  CONTINUE

    RETURN
END SUBROUTINE ampb
