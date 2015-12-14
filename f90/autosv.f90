SUBROUTINE autosv
     
    !     THIS ROUTINE GENERATES OSCAR ENTRIES FOR PARAMTERS
    !     THAT ARE TO BE SAVED IMPLICITLY
 
    EXTERNAL         lshift,andf
    INTEGER :: savnam,osprc,osbot,ospnt,oscar(1),os(5),xsav(2), dmpcnt,andf,vps
    COMMON  /xgpic / junk(25),maskhi,junk1(2),nosgn
    COMMON  /autosm/ nwords,savnam(100)
    COMMON  /zzzzzz/ core(1)
    COMMON  /xvps  / vps(1)
    COMMON  /xgpi4 / junk4(2),iseqn,dmpcnt
    COMMON  /autohd/ ihead
    EQUIVALENCE      (core(1),os(1),loscar), (os(2),osprc),  &
        (os(3),osbot), (os(4),ospnt), (os(5),oscar(1))
    DATA     xsav  / 4HXSAV,4HE    /
 
    !     UPDATE OSCAR PARAMETERS
 
    ihead = 1
    osprc = osbot
    osbot = oscar(osbot) + osbot
    ospnt = osbot
    iseqn = oscar(osprc+1) + 1
 
    !     LOAD HEADER
 
    oscar(ospnt  ) = 6
    oscar(ospnt+1) = iseqn
    oscar(ospnt+2) = 4 + lshift(8,16)
    oscar(ospnt+3) = xsav(1)
    oscar(ospnt+4) = xsav(2)
    oscar(ospnt+5) = dmpcnt
    CALL xlnkhd
 
    !     HAVING THE VPS POINTERS FOR EACH PARAMETER, FIND THE
    !     DISPLACEMENT IN COMMON
 
    j = osprc + 6 + 3*oscar(osprc+6) + 1
    IF (andf(oscar(osprc+2),maskhi) == 1) j = j+1+3*oscar(j)
    j  = j + 1
    n3 = j+1
    n1 = oscar(j)
    n2 = 1
    oscar(ospnt+6) = nwords
    oscar(ospnt  ) = oscar(ospnt) + 1
    ipt = 1
    ist = n3
140 IF (oscar(ist) > 0) GO TO 110
 
    !     SEEE IF PARAMETER IS IN SAVE LIST
 
    !     LL = ANDF(OSCAR(IST),NOSGN )  REPLACED BY NEXT CARD, OCT. 1983
    ll = andf(oscar(ist),maskhi)
    l  = andf(vps(ll-1) ,maskhi)
    DO  i = 1,nwords
        IF (andf(oscar(ist),nosgn) == savnam(i)) GO TO 120
    END DO
 
    !     NOT TO BE SAVED, GO TO NEXT PARAMETER
 
    ist = ist + 1
    n2  = n2  + l
    GO TO 140
 
    !     CONSTANT PARAMETER, SKIP IT
 
110 nwd = oscar(ist)
    ist = ist+nwd+1
    n2  = n2 + nwd
    GO TO 140
 
    !     PARAMETER TO BE SAVED, PUT IN OSCAR
 
120 oscar(ospnt+6+2*i-1) = savnam(ipt)
    oscar(ospnt+6+2*i  ) = n2
    oscar(ospnt) = oscar(ospnt) + 2
    ipt = ipt + 1
    ist = ist + 1
    n2  = n2  + l
    IF (ipt <= nwords) GO TO  140
    ihead = 0
    RETURN
END SUBROUTINE autosv
