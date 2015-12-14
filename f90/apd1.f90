SUBROUTINE apd1 (fst,ns,fct,nc,ls,lc)
 
    REAL, INTENT(IN OUT)                     :: fst(1)
    INTEGER, INTENT(IN)                      :: ns
    REAL, INTENT(IN OUT)                     :: fct(1)
    INTEGER, INTENT(IN)                      :: nc
    LOGICAL, INTENT(IN OUT)                  :: ls
    LOGICAL, INTENT(IN OUT)                  :: lc
 
    INTEGER :: iz(1),NAME(2),FILE,ncary(2),silc,necta(6),  &
        cp,acsid,eid,eidb,cid(5),cidbx,auset(6,2),silb,  &
        rdrew,clsrew,ays(5),key(5),sildx(4),acsix(4),back,  &
        scr1,scr2,scr3,scr4,scr5,ecta,bgpa,gpla,useta,  &
        sila,cstma,acpt,buf10,buf11,buf12,wtrew,acsib,pid
    REAL :: rb1(3),acpl(3,3),vx1(3),vx2(3),axic(3),  xb(5)
    COMMON /BLANK / nk,nj,luseta
    COMMON /system/ sysbuf,NOT
    COMMON /apd1c / eid,pid,cp,nspan,nchord,lspan,lchord,igid,  &
        x1,y1,z1,x12,x4,y4,z4,x43,xop,x1p,alzo,mcstm,  &
        ncst1,ncst2,cidbx,acsid,iacs,silb,ncrd,  &
        scr1,scr2,scr3,scr4,scr5,ecta,bgpa,gpla,useta,  &
        sila,cstma,acpt,buf10,buf11,buf12,next,left,isiln
    COMMON /apd1d / icpl(14),yp4,s1,c1,xp2,xp3,xp4,ra1(3)
    COMMON /apd12c/ key,auset,usa,uk,ncam2,nasb,ippc
    COMMON /zzzzzz/ z(1)
    EQUIVALENCE     (icpl(3),rb1(1)),(icpl(6),acpl(1,1)),  &
        (necta(2),cid(1)),(key(2),np),(key(3),nstrip),  &
        (key(4),ntp),(key(5),f),(ays(1),ys),(ays(2),zs),  &
        (ays(3),ee),(ays(4),sg),(ays(5),cg),(axic(1),xic),  &
        (axic(2),delx),(axic(3),xlam),(sildx(1),icid),  &
        (sildx(3),silc),(acsix(1),acsib),(z(1),iz(1)),  &
        (acsix(2),vx2(1)),(necta(1),eidb)
    DATA    rdrew , clsrew,wtrew / 0,1,1 /
    DATA    NAME  / 4HAPD1,1H    /
 
    key(1) = 1
    silc   = silb
 
    !     IF NEW IGRID SET INITIALIZE
 
    IF (.NOT.ls) GO TO 30
    np   = 0
    ntp  = 0
    nbox = 0
    nasb = 0
    nstrip = 0
    CALL gopen (scr3,z(buf10),wtrew)
    CALL gopen (scr4,z(buf11),wtrew)
    CALL gopen (scr5,z(buf12),wtrew)
 
    !     MAKE COORD SYSTEM AND GET POINTS IN PROPER SYSTEM
 
30  CALL apdcs
    sg = s1
    cg = c1
    acsib = mcstm
 
    !     CHECK FOR ASSOCIATED BODIES
 
    DO  j = 1,6
        IF (iz(ippc+j) == 0) GO TO 45
        nasb = nasb + 1
    END DO
45 CONTINUE
 
   !     GENERATE BOXES
 
   ncrdp= ncrd
   np   = np + 1
   fsj1 = apdf(fst,1,nspan)
   yj1  = fsj1*yp4
   dj1  = fsj1*xp4
   cj1  = (1.0-fsj1)*xp2 + fsj1*(xp3-xp4)
   eidb = eid - 1
   DO  j = 1,ns
       yj   = yj1
       dj   = dj1
       cj   = cj1
       fsj1 = apdf(fst,j+1,nspan)
       yj1  = fsj1*yp4
       dj1  = fsj1*xp4
       cj1  = (1.0-fsj1)*xp2 + fsj1*(xp3-xp4)
       ee   = .5*(yj1-yj)
       ysp  = yj + ee
       nstrip = nstrip + 1
       fci1 = apdf(fct,1,nchord)
       xi1j = dj + fci1*cj
       xi1j1= dj1+ fci1*cj1
       ds   = 1.0/(yj1-yj)
       ys   = ysp*cg + ra1(2)
       zs   = ysp*sg + ra1(3)
       CALL WRITE (scr3,ays(1),5,0)
       DO  i = 1,nc
           ntp  = ntp + 1
           xij  = xi1j
           xij1 = xi1j1
           fci1 = apdf(fct,i+1,nchord)
           xi1j = dj + fci1*cj
           xi1j1= dj1+ fci1*cj1
           aij  = (1.0-xop)*xij  + xop*xi1j
           aij1 = (1.0-xop)*xij1 + xop*xi1j1
           xic  = .5*(aij+aij1)  + ra1(1)
           xlam = (aij1-aij)*ds
           delx = .50*(-xij+xi1j - xij1+xi1j1)
           CALL WRITE (scr4,axic(1),3,0)
           xic  = xic - ra1(1)
           eidb = eidb + 1
           nbox = nbox + 1
           cid(1) = cidbx + i + (nc+1)*(j-1)
           cid(2) = cid(1) + 1
           cid(3) = cid(1) + nc + 1
           cid(4) = cid(3) + 1
           cid(5) = eidb
           ncid = cid(4)
           nj   = nj + 1
           nk   = nk + 2
           vx1(3) = 0
           IF (j /= 1) GO TO 310
           IF (i /= 1) GO TO 300
           ASSIGN 300 TO back
           icid = cid(1)
           vx1(1) = xij
           vx1(2) = yj
           kk = 1
           GO TO 340
300        ASSIGN 310 TO back
           icid = cid(2)
           vx1(1) = xi1j
           vx1(2) = yj
           kk = 1
           GO TO 340
310        IF (i /= 1) GO TO 320
           ASSIGN 320 TO back
           icid = cid(3)
           vx1(1) = xij1
           vx1(2) = yj1
           kk = 1
           GO TO 340
320        ASSIGN 330 TO back
           icid = cid(4)
           vx1(1) = xi1j1
           vx1(2) = yj1
           kk = 1
           GO TO 340
330        ASSIGN 360 TO back
           icid = cid(5)
           vx1(1) = xic + .25*delx
           vx1(2) = ysp
           kk = 2
340        CALL gmmats (acpl,3,3,0, vx1,3,1,0, vx2)
           DO  k = 1,3
               vx2(k) = vx2(k) + rb1(k)
           END DO
           CALL WRITE (bgpa,acsix,4,0)
           CALL WRITE (gpla,icid,1,0)
           CALL WRITE (useta,auset(1,kk),6,0)
           ncrd = ncrd + 1
           silc = silc + 6
           isiln= isiln + 6
           sildx(4) = isiln
           luseta = silc
           sildx(2) = 10*silc + 1
           CALL WRITE (sila,silc,1,0)
           CALL WRITE (scr2,isiln,1,0)
           CALL WRITE (scr2,silc,1,0)
           CALL WRITE (scr1,icid,2,0)
           GO TO back, (300,310,320,330,360)
360        cid(1) = iapd(i  ,j  ,nc,ncrdp)
           cid(2) = iapd(i+1,j  ,nc,ncrdp)
           cid(4) = iapd(i  ,j+1,nc,ncrdp)
           cid(3) = iapd(i+1,j+1,nc,ncrdp)
           cid(5) = cid(3) + 1
           CALL WRITE (ecta,necta(1),6,0)
       END DO
   END DO
   cidbx = ncid
   ncary(1) = nc
   ncary(2) = nbox
   CALL WRITE (scr5,ncary,2,0)
 
   !     ADD PROPERITY CARD POINTERS FOR APD2
 
   CALL WRITE (scr5,ippc,1,0)
   silb = silc
   IF(.NOT.lc) RETURN
 
   !     WRITE ACPT TABLE
 
   f = x1p - xop
   CALL WRITE (acpt,key,5,0)
 
   !     COPY STUFF FROM SCRATCH FILES TO ACPT
 
   FILE = scr5
   k = 3
   ASSIGN 410 TO iret
   GO TO 375
410 ASSIGN 420 TO iret
   FILE = scr3
   k = 5
   GO TO 375
420 ASSIGN 430 TO iret
   FILE = scr4
   k = 3
375 CALL WRITE (FILE,0,0,1)
   CALL CLOSE (FILE,clsrew)
   CALL gopen (FILE,z(buf12),rdrew)
   DO  i = 1,k
380    CALL READ (*480,*390,FILE,xb(1),k,0,j)
   
       !     SKIP PROPERTY CARD POINTERS
   
       IF (i == 3 .AND. FILE == scr5) GO TO 380
       CALL WRITE (acpt,xb(i),1,0)
       GO TO 380
390    CALL REWIND (FILE)
       CALL skprec (FILE,1)
   END DO
   CALL CLOSE (FILE,clsrew)
   GO TO iret, (410,420,430)
430 CALL WRITE (acpt,0,0,1)
   RETURN
 
   !     ERROR MESAGES
 
480 ip1 = -2
   CALL mesage (ip1,FILE,NAME)

   RETURN
END SUBROUTINE apd1
