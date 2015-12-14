SUBROUTINE cnstrc (gp,ele,buf,MAX)
     
    !     THIS SUBROUTINE BUILDS THE ELSETS FILE
    !     THIS SUBROUTINE IS CALLED ONLY BY DPLTST, WHICH IS THE DRIVER OF
    !     DMAP MODULE PLTSET
    !     THE SUBROUITNE PLTSET OF THE PLOT MODULE HAS NOTHING TO DO WITH
    !     THIS SUBROUTINE
 
    !     REVISED 10/1990 BY G.CHAN/UNISYS TO INCLUDE OFFSET FOR BAR, TRIA3
    !     AND QUAD4 ELEMENTS
 
 
    INTEGER, INTENT(OUT)                     :: gp(1)
    INTEGER, INTENT(OUT)                     :: ele(1)
    INTEGER, INTENT(IN OUT)                  :: buf(1)
    INTEGER, INTENT(IN OUT)                  :: MAX
    INTEGER :: ae       ,b1       ,b2       ,b3       ,  &
        bufsiz   , elid     ,ERR(2)   ,etype   ,  &
        exgp     , gpt      ,gpts(32) ,outnor  ,  &
        outrew   ,rew      ,setid    ,setnum   ,SIGN    ,  &
        TYPE(50) ,NAME(2)  ,msg1(14) ,eid(2)   ,offset  ,  &
        br       ,t3       ,q4       ,off(6)
    COMMON /BLANK / ngp      ,nsets    ,skp1(8)  ,skp2     ,exgpid  ,  &
        skp3(8)  ,merr     ,skp4     ,grid     ,ect2    , skp5(6)  ,mset     ,ect1
    COMMON /system/ bufsiz
    COMMON /names / rd       ,inprew   ,outnor   ,outrew   ,rew     , norew
    COMMON /gpta1 / ntyps    ,last     ,incr     ,NE(1)
    EQUIVALENCE     (eid(1)  ,elid)
    DATA    NAME  / 4H cns   ,4HTRC  / ,ae   /    72       /
    DATA    nmsg1 / 14       /
    DATA    msg1  / 4H(33X   ,4H,44H   ,4HNO p   ,4HLOTA   ,4HBLE   ,  &
        4HSTRU   ,4HCTUR   ,4HAL e   ,4HLEME   ,4HNTS   ,  &
        4HEXIS   ,4HT in   ,4H set   ,4H,i8)   /
    DATA    br, t3, q4       /2HBR     ,2HT3     ,2HQ4     /
 
    b1 = 1
    b2 = b1 + bufsiz
    b3 = b2 + bufsiz
    CALL gopen (mset,buf(b3),inprew)
    CALL gopen (ect2,buf(b2),outrew)
 
    DO  setnum = 1,nsets
        CALL fread (mset,setid,1,0)
        DO  i = 1,ngp
            gp(i) = 0
        END DO
   
        !     READ THE EXPLICIT ELEMENT NUMBERS IN THIS SET.
   
        CALL fread (mset,nel,1,0)
        IF (nel >= MAX) CALL mesage (-8,0,NAME)
        ele(nel+1) = 0
        CALL fread (mset,ele,nel,0)
   
        !     READ THE ELEMENT TYPES TO BE INCLUDED OR EXCLUDED IN THIS SET.
   
        CALL fread (mset,ntypes,1,0)
        CALL fread (mset,TYPE,ntypes,0)
   
        !     GENERATE AN ECT FOR THE ELEMENTS INCLUDED IN THIS SET.
   
        CALL gopen (ect1,buf(b1),inprew)
20      CALL READ (*300,*300,ect1,etype,1,0,i)
   
        !     CHECK WHETHER OR NOT THIS ELEMENT TYPE IS TO BE EXCLUDED.
   
        mtype = -1
        labgp =  1
        IF (etype == ae) labgp = -2
        IF (ntypes == 0) GO TO 50
        DO  i = 1,ntypes,2
            IF (-etype == TYPE(i)) GO TO 40
        END DO
        GO TO 50
40      mtype = TYPE(i+1)
   
        !     THIS ELEMENT TYPE MAY BE INCLUDED AS A TYPE AND/OR SOME OF THEM
        !     MAY BE INCLUDED SPECIFICALLY. READ -NGPPE- = NUMBER OF GRID
        !     POINTS PER ELEMENT FOR THIS TYPE.
   
50      CALL fread (ect1,ngppe,1,0)
        IF (ngppe <= 0) GO TO 290
        idx    = (etype-1)*incr
        neltyp = 0
        ne16   = NE(idx+16)
        offset = 0
        IF (ne16 == br) offset = 6
        IF (ne16 == t3 .OR. ne16 == q4) offset = 1
   
        !     CHECK WHETHER OR NOT THIS ELEMENT TYPE IS TO BE INCLUDED.
   
        IF (ntypes == 0 .OR. mtype >= 0) GO TO 70
        DO  i = 1,ntypes,2
            IF (etype == TYPE(i) .OR. TYPE(i) == ntyps+1) GO TO 200
        END DO
   
        !     NOW CHECK WHETHER OR NOT ANY OF THE ELEMENTS OF THIS TYPE ARE
        !     EXPLICITLY INCLUDED. PUT ALL SUCH ON THE NEW ECT (ECT2).
   
70      CALL READ (*280,*280,ect1,eid,2,0,i)
        CALL fread (ect1,gpts,ngppe,0)
        IF (offset /= 0) CALL fread (ect1,off,offset,0)
        IF (nel <= 0) GO TO 70
        m = 0
        n = 1
   
        !     FOR TYPES DELETED ONLY SEARCH LIST AFTER TYPE WAS KNOWN TO BE
        !     DELETED (2ND WORD OF TYPE)
   
        IF (mtype > 0) n = mtype
        IF (n   > nel) GO TO 110
80      CALL intlst (ele,n,SIGN,n1,n2)
        IF (SIGN < 0) GO TO 90
        IF (elid >= n1 .AND. elid <= n2) m = 1
        GO TO 100
90      IF (elid >= n1 .AND. elid <= n2) m = 0
100     IF (n <= nel) GO TO 80
110 CONTINUE
    IF (m      == 0) GO TO 70
    IF (neltyp /= 0) GO TO 120
    CALL WRITE (ect2,NE(idx+16),1,0)
    CALL WRITE (ect2,ngppe,1,0)
120 CALL WRITE (ect2,eid,2,0)
    CALL WRITE (ect2,gpts,ngppe,0)
    IF (offset /= 0) CALL WRITE (ect2,off,offset,0)
    neltyp = neltyp + 1
    DO  i = 1,ngppe
        j = gpts(i)
        gp(j) = labgp
    END DO
   
    !     AERO ELEMENT - CENTER ONLY LABELED
   
    IF (etype == ae) gp(j) = 1
    GO TO 70
   
    !     THIS ELEMENT TYPE IS TO BE INCLUDED, EXCEPT THE ONES EXPLICITLY
    !     EXCLUDED
   
    !     ONLY SEARCH LIST AFTER TYPE WAS INCLUDED
   
200 mtype = TYPE(i+1)
210 CALL READ (*280,*280,ect1,eid,2,0,i)
    CALL fread (ect1,gpts,ngppe,0)
    IF (offset /= 0) CALL fread (ect1,off,offset,0)
    IF (nel <= 0) GO TO 250
    m = 1
    n = 1
    IF (mtype > 0) n = mtype
    IF (n   > nel) GO TO 250
220 CALL intlst (ele,n,SIGN,n1,n2)
    IF (SIGN > 0) GO TO 230
    IF (elid >= n1 .AND. elid <= n2) m = 0
    GO TO 240
230 IF (elid >= n1 .AND. elid <= n2) m = 1
240 IF (n    <= nel) GO TO 220
    IF (m    ==   0) GO TO 210
250 IF (neltyp /= 0) GO TO 260
    CALL WRITE (ect2,NE(idx+16),1,0)
    CALL WRITE (ect2,ngppe,1,0)
260 CALL WRITE (ect2,eid,2,0)
    CALL WRITE (ect2,gpts,ngppe,0)
    IF (offset /= 0) CALL WRITE (ect2,off,offset,0)
    DO  i = 1,ngppe
        j = gpts(i)
        gp(j) = labgp
    END DO
   
    !     AERO ELEMENT - CENTER ONLY LABELED
   
    IF (etype == ae) gp(j) = 1
    neltyp = neltyp + 1
    GO TO 210
   
    !     END OF NEW ECT FOR THIS ELEMENT TYPE
   
280 IF (neltyp > 0) CALL WRITE (ect2,0,1,0)
    GO TO 20
   
    !     SKIP THIS ELEMENT TYPE (NON-EXISTENT)
   
290 CALL fread (ect1,0,0,1)
    GO TO 20
   
    !     END OF ECT FOR THIS ELEMENT SET
   
300 CALL CLOSE (ect1,rew)
    CALL WRITE (ect2,0,0,1)
   
    !     FLAG ALL GRID POINTS TO BE EXCLUDED FROM A DEFORMED SHAPE.
   
    CALL fread (mset,ngpts,1,0)
    IF (ngpts >= MAX) CALL mesage (-8,0,NAME)
    ele(ngpts+1) = 0
    CALL fread (mset,ele,ngpts,1)
    IF (ngpts <= 0) GO TO 400
    CALL gopen (exgpid,buf(b1),inprew)
    DO  gpt = 1,ngp
        CALL fread (exgpid,exgp,1,0)
        CALL fread (exgpid,ingp,1,0)
        m = 0
        n = 1
310     CALL intlst (ele,n,SIGN,n1,n2)
        IF (SIGN > 0) GO TO 320
        IF (exgp >= n1 .AND. exgp <= n2) m = ingp
        GO TO 330
320     IF (exgp >= n1 .AND. exgp <= n2) m = 0
330     IF (n <= ngpts) GO TO 310
        IF (m ==     0) CYCLE
        IF (gp(m) /= -2) gp(m) = -gp(m)
    END DO
    CALL CLOSE (exgpid,rew)
   
    !     GENERATE A GRID POINT LIST FOR THIS SET (CONVERT THE INTERNAL
    !     GRID POINT NUMBERS TO POINTERS TO THE GRID POINTS PECULIAR TO
    !     THIS SET)
   
400 CALL gopen (grid,buf(b1),outnor)
    ngpts = 0
    DO  i = 1,ngp
        IF (gp(i) == 0) CYCLE
        ngpts = ngpts+1
        gp(i) = ISIGN(ngpts,gp(i))
    END DO
    IF (ngpts /= 0) GO TO 420
    ERR(1) = 1
    ERR(2) = setid
    CALL wrtprt (merr,ERR,msg1,nmsg1)
   
420 CALL WRITE (grid,ngpts,1,0)
    CALL WRITE (grid,gp,ngp,0)
    IF (setnum /= nsets) CALL CLOSE (grid,norew)
END DO
 
!     ALL DONE. THE SET DEFINITION FILE (MSET) + THE SHORT ECT FILE
!     (ECT1) WILL NOT BE NEEDED AGAIN.
 
CALL clstab (grid,rew)
CALL clstab (ect2,rew)
CALL CLOSE  (mset,rew)

RETURN
END SUBROUTINE cnstrc
