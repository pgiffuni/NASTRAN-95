SUBROUTINE contor (gplst,x,u,dd,z,iz,ppen,deform,b1,opcor)
     
 
    INTEGER, INTENT(IN OUT)                  :: gplst(1)
    REAL, INTENT(IN)                         :: x(3,1)
    REAL, INTENT(IN)                         :: u(2,1)
    REAL, INTENT(IN OUT)                     :: dd(3,1)
    REAL, INTENT(OUT)                        :: z(1)
    INTEGER, INTENT(OUT)                     :: iz(1)
    INTEGER, INTENT(IN)                      :: ppen
    INTEGER, INTENT(IN OUT)                  :: deform
    INTEGER, INTENT(IN)                      :: b1
    INTEGER, INTENT(IN)                      :: opcor
    INTEGER :: pen, parm,est, stress,sort,scr1,bufsiz, b2,b3,elid,ERR,sub(2),  &
        color,elmid,gpts(12),pedge,esym,offset
    REAL :: xb(8), center,rcntrl,rcolor
    DIMENSION       labl(50), ibegin(2),pt(8)
    CHARACTER (LEN=29) :: uim
    CHARACTER (LEN=25) :: uwm
    CHARACTER (LEN=23) :: ufm
    COMMON /xmssg / ufm,uwm,uim
    COMMON /BLANK / skip(10),parm,skp4,est,skp8(11),stress,sort, newoes,scr1
    COMMON /xxparm/ skip2(157),ncntr,cntr(50),icntvl,skip6(5),iset,  &
        sk18(18),color,layer
    COMMON /pltdat/ skip3(2),xmin
    COMMON /system/ bufsiz,nout
    COMMON /drwdat/ jset,skip7(14),pedge
    EQUIVALENCE     (ibegin(1),lines),(ibegin(2),igdpt)
    DATA    kbar  , kt3,kq4/ 2HBR,2HT3,2HQ4 /
    DATA    sub   / 4HCONT , 4HOR           /
 
    b2  = b1 - 2*bufsiz
    b3  = b2 - bufsiz
    lopcor  = opcor/5
    id  = 1
    ERR = 0
    irr = 0
    ncntr = IABS(ncntr)
 
    !     COLOR = 0 IS NO COLOR CONTOUR,
    !     COLOR = 1 TO  31 IS DRAW CONTOUR LINES IN COLOR.
    !     COLOR =-1 TO -31 ID COLOR FILL ELEMENTS BY STRESS
 
    !     THIS IS THE CODE FOR THE COLOR BAR SCALE
    !     AT THE TOP OF THE PLOT
 
    IF (color == 0) GO TO 40
    CALL line (0.0,0.0,0.0,0.0,32,-1)
    icolor = IABS(color)
    rcolor = 535.0/icolor
    !              535.0 IS BASED ON 1000X1000 PLOT SCREEN COORDINATE FRAME
    DO  ib = 1,icolor
        pen   = ib + 31
        xb(1) = 368.14 + rcolor*(ib-1)
        xb(2) = 959.67
        xb(3) = xb(1) + rcolor
        xb(4) = xb(2)
        xb(5) = xb(3)
        xb(6) = 977.51
        xb(7) = xb(1)
        xb(8) = xb(6)
        ik    = 0
        ik1   = 1
        ik2   = 2
        ik3   = 3
        ik4   = 4
        DO  ii = 1,4
            ik1 = ik1 + ik
            ik2 = ik2 + ik
            ik3 = ik3 + ik
            ik4 = ik4 + ik
            ik  = 2
            IF (ii == 4) pen = 0
            IF (ii /= 4) GO TO 10
            ik1 = 7
            ik2 = 8
            ik3 = 1
            ik4 = 2
10          CALL line (xb(ik1),xb(ik2),xb(ik3),xb(ik4),pen,0)
        END DO
    END DO
    CALL line (0.0,0.0,0.0,0.0,pen,+1)
40  isav = lopcor + id
    ival = lopcor + isav
    icen = lopcor + ival
    lopcor= lopcor- 1
    pen = ppen
    IF (icntvl > 9 .AND. icntvl < 14 .AND. pedge /= 1) GO TO 50
    CALL CLOSE (parm,2)
    IF (icntvl <= 9 .OR. icntvl > 13) CALL create (gplst,x,u,deform,  &
        conmin,conmax,z(ival),z(icen),lopcor,b1,b2)
    IF (iset /= jset) CALL order (gplst,z(ival),z(isav),z(icen),  &
        z(id),lopcor,b1,b2,b3)
    iset = jset
50  IF (icntvl > 9 .AND. icntvl < 14)  &
        CALL displa (gplst,x,u,dd,pen,deform,labl,pt,b1)
    IF (icntvl > 9 .AND. icntvl < 14) GO TO 420
    IF (conmin == conmax) GO TO 470
    IF (color < 0) CALL CLOSE (est,2)
    CALL gopen (scr1,gplst(b1),1)
    CALL CLOSE (scr1,2)
    CALL gopen (sort,gplst(b2),2)
    CALL gopen (stress,gplst(b1),0)
 
    !     BUFFERS ASSIGNEMENT HERE -
    !     B1 IS USED BY STRESS (SCRATCH1/301), AND BY SCR1 (SCRATCH4/304)
    !     FOR SHORT PERIODS OF TIME ONLY
    !     B2 IS USED BY SORT (SCRATCH2/302)
    !     B3 IS USED BY EST (ELEST/103) AND BY SCR1 (SCRATCH4/304)
 
    ncntr = MIN0(ncntr,50)
    IF (cntr(1) /= cntr(2)) GO TO 90
 
    !     IF INTERVALS SPECIFIED, DEFINE CONTOUR VALUES
 
    delta = (conmax-conmin)/FLOAT(ncntr-1)
    cntr(1) = conmin
    j = ncntr - 1
    DO  i = 2,j
        cntr(i) = cntr(i-1) + delta
    END DO
    cntr(ncntr) = conmax
90  CALL line (0.,0.,0.,0.,pen,-1)
    DO  i = 1,ncntr
        labl(i) = 3
    END DO
 
    !     READ AND STORE CONTOUR VALUES AND CENTROIDS
 
    elid   = 0
    lopcox = lopcor + 1
    DO  i = 1,lopcox
        is = isav + i - 1
        iz(is) = 0
    END DO
    IF (color >= 0) GO TO 130
    CALL gopen (est,gplst(b3),2)
    CALL bckrec (est)
    imhere = 120
120 CALL READ (*280,*270,est,esym,1,0,m)
    irr = 0
    CALL fread (est,ngppe,1,0)
130 CALL fwdrec (*415,sort)
140 CALL READ (*415,*415,sort,iflag,1,0,m)
    IF (iflag ==  0) GO TO 415
    IF (iflag == -2) GO TO 130
    CALL fread (sort,ibegin,2,0)
    CALL READ  (*415,*150,sort,iz(id),lines,1,i)
150 iread = 0
    nel = 0
    loop170:  DO  i = 1,lines
        ic = icen + 2*(i-1)
        iv = ival + i - 1
        id1= id   + i - 1
        is = isav + i - 1
        DO  j = 1,lopcor
            js = isav + j - 1
            jv = ival + j - 1
            jc = icen + 2*(j-1)
            IF (iz(js)  == 0)  CYCLE loop170
            IF (iz(id1) /= iz(js)) CYCLE
            z(iv  ) = z(jv)
            z(ic  ) = z(jc)
            z(ic+1) = z(jc+1)
            iz(is ) = iz(id1)
            nel = nel + 1
            CYCLE loop170
        END DO
    END DO loop170
    IF (elid > 0) GO TO 190
180 CALL READ (*290,*290,stress,essym,1,0,m)
190 CALL READ (*290,*180,stress,elid,1,0,m)
    IF (elid == 0) GO TO 180
    CALL fread (stress,v,1,0)
    CALL fread (stress,pt,2,0)
    DO  i = 1,lines
        id1 = id  + i - 1
        is = isav + i - 1
        iv = ival + i - 1
        ic = icen + 2*i - 2
        IF (iz(id1) /= elid) CYCLE
        IF (iz(is)  == elid) CYCLE
        z(iv  ) = v
        z(ic  ) = pt(1)
        z(ic+1) = pt(2)
        iz(is ) = elid
        IF (color >= 0) GO TO 240
        imhere = 200
        ASSIGN 206 TO irtn
195     jrr = 0
200     offset = 0
        IF (esym == kbar) offset = 6
        IF (esym == kt3 .OR. esym == kq4) offset = 1
201     CALL READ (*204,*205,est,elmid,1,0,m)
        IF (elmid == 0) GO TO 203
        CALL fread (est,0,-1,0)
        CALL fread (est,gpts,ngppe,0)
        IF (offset /= 0) CALL fread (est,0,-offset,0)
        IF (elmid == elid) GO TO 210
        GO TO 201
   
203     jrr = jrr + 1
        IF (jrr <= 1) CALL bckrec (est)
        imhere = 203
        ASSIGN 120 TO irtn
        CALL READ (*204,*204,est,esym,1,0,m)
        CALL fread (est,ngppe,1,0)
        GO TO 200
   
204     ERR = ERR + 1
        IF (ERR > 3) GO TO 285
        CALL REWIND (est)
        CALL skprec (est,1)
        GO TO irtn, (120,205,206)
   
205     imhere = 205
        ASSIGN 205 TO irtn
206     CALL READ (*204,*205,est,esym,1,0,m)
        CALL fread (est,ngppe,1,0)
        GO TO 195
   
        !     START TO CONTOUR FILL HERE
   
210     rcolor = icolor
        rcntrl = ncntr
        DO  ik = 1,ncntr
            pen = 32 + (1.0-(rcntrl-ik+1)/rcntrl)*rcolor
            ik1 = ik + 1
            IF (ik == ncntr) ik1 = ik
            IF (v < cntr(ik) .OR. v > cntr(ik1)) CYCLE
            DO  j = 1,ngppe
                k  = j + 1
                ig = gpts(j)
                ig = IABS(gplst(ig))
                IF (j == ngppe) k = 1
                ig1 = gpts(k)
                ig1 = IABS(gplst(ig1))
                IF (j == ngppe) pen = 0
                CALL line (x(2,ig),x(3,ig),x(2,ig1),x(3,ig1),pen,0)
            END DO
            EXIT
        END DO
240     nel = nel + 1
        EXIT
    END DO
260 IF (nel >= lines) GO TO 300
    GO TO 190
 
270 CALL bckrec (est)
    irr = irr + 1
    IF (irr < 3) GO TO 120
 
    !     END OF FILE ON EST
 
280 ERR = ERR + 1
    IF (ERR > 3) GO TO 285
    CALL REWIND (est)
    CALL skprec (est,1)
    GO TO 120
285 WRITE  (nout,286) uim,elid,imhere,ERR,irr,ngppe
286 FORMAT (a29,', CONTOUR FAILED TO LOCATE ELMENT ID =',i8, /5X,  &
        'IMHERE =',i5, 5X,'ERR,IRR,NGPPE =',3I8)
    GO TO 190
 
    !     END OF FILE ON STRESS
 
290 CALL REWIND (stress)
    CALL fwdrec (*415,stress)
    IF (iread == 1) GO TO 140
    iread = 1
    GO TO 180
 
    !     END DATA SEARCH
 
300 l  = lines
    is = lines + isav
    iz(is) = 0
    IF (lines > 3) GO TO 310
    xmid   = z(icen+4)
    ymid   = z(icen+5)
    cenval = z(ival+2)
    l = 1
    GO TO 350
310 ig = IABS(gplst(igdpt))
    IF (deform /= 0) GO TO 320
    xmid = x(2,ig)
    ymid = x(3,ig)
    GO TO 330
320 xmid = u(1,ig)
    ymid = u(2,ig)
330 sum1 = 0.0
    sum2 = 0.0
    DO  i = 1,lines
        iv = ival + i - 1
        ic = icen + 2*i - 2
        s  = SQRT((xmid-z(ic))**2 + (ymid-z(ic+1))**2)
        sum1 = sum1 + z(iv) * s
        sum2 = sum2 + s
    END DO
    cenval = sum1/sum2
350 iv = ival + lines
    ic = icen + 2*lines
    z(iv  ) = z(ival)
    z(ic  ) = z(icen)
    z(ic+1) = z(icen+1)
 
    !     PLOT CONTOURS.
 
    IF (color < 0) GO TO 140
    rcolor = icolor
    rcntrl = ncntr
 
    CALL CLOSE (est,2)
    CALL gopen (scr1,gplst(b3),3)
 
    DO  i = 1,ncntr
        IF (color /= 0) pen = 1 + (1.0-(rcntrl-i+1)/rcntrl)*rcolor
        DO  j = 1,l
            pt(1) = xmin - 1.0
            pt(3) = pt(1)
            pt(5) = pt(1)
            jc    = icen + 2*j - 2
            jv    = ival + j - 1
            d     = (z(jv) - z(jv+1))
            IF (ABS(z(jv)-cntr(i)) > ABS(d) .OR. ABS(z(jv+1)-cntr(i))  &
                > ABS(d)) GO TO 360
            IF (d == 0.0) d = 1.0
            pt(1) = z(jc  ) + (z(jc+2)-z(jc  ))*(z(jv)-cntr(i))/d
            pt(2) = z(jc+1) + (z(jc+3)-z(jc+1))*(z(jv)-cntr(i))/d
360         d = z(jv+1) - cenval
            IF (ABS(z(jv+1)-cntr(i)) > ABS(d) .OR. ABS(cenval-cntr(i))  &
                > ABS(d)) GO TO 370
            IF (d == 0.0) d = 1.0
            pt(3) = z(jc+2) + (xmid-z(jc+2))*(z(jv+1)-cntr(i))/d
            pt(4) = z(jc+3) + (ymid-z(jc+3))*(z(jv+1)-cntr(i))/d
370         d = cenval - z(jv)
            IF (ABS(cenval-cntr(i)) > ABS(d) .OR.  &
                ABS(z(jv)-cntr(i)) > ABS(d)) GO TO 380
            IF (d == 0.0) d = 1.0
            pt(5) = xmid + (z(jc  )-xmid)*(cenval-cntr(i))/d
            pt(6) = ymid + (z(jc+1)-ymid)*(cenval-cntr(i))/d
380         pt(7) = pt(1)
            pt(8) = pt(2)
            DO  k = 1,5,2
                IF (pt(k) < xmin .OR. pt(k+2) < xmin) CYCLE
                CALL line (pt(k),pt(k+1),pt(k+2),pt(k+3),pen,0)
                labl(i) = labl(i) + 1
                IF (labl(i) /= 4) CYCLE
                labl(i) = 0
                CALL WRITE (scr1,i,1,0)
                CALL WRITE (scr1,pt(k),2,0)
            END DO
        END DO
    END DO
 
    CALL CLOSE (scr1,2)
    CALL gopen (est,gplst(b3),2)
    GO TO 140
 
415 CALL CLOSE (sort,1)
    CALL CLOSE (stress,1)
    CALL CLOSE (scr1,1)
    !     IF (COLOR .LT. 0) CALL CLOSE (EST,1)
    !     IF (COLOR .GE. 0) CALL GOPEN (EST,GPLST(B3),2)
420 CALL line (0.,0.,0.,0.,pen,+1)
    IF (color == 0) GO TO 430
    CALL typflt (0.0,0.0,0,0,0,-1)
    CALL typflt (368.14,990.0,1,cntr(1),-8,0)
    center = (cntr(1)+cntr(ncntr))/2.0
    CALL typflt (585.90,990.0,1,center,-8,0)
    CALL typflt (796.3,990.0,1,cntr(ncntr),-8,0)
    CALL typflt (0.0,0.0,0,0,0,+1)
    IF (color < 0) GO TO 460
430 CALL gopen (scr1,gplst(b1),0)
    IF (color == 0) CALL typint (0.,0.,0,0,0,-1)
440 CALL READ (*450,*450,scr1,i,1,0,m)
    CALL fread (scr1,pt,2,0)
    IF (color == 0) CALL typint (pt(1),pt(2),1,i,1,0)
    GO TO 440
450 IF (color == 0) CALL typint (0.,0.,0,0,0,+1)
    CALL CLOSE (scr1,1)
460 CALL pltopr
470 IF ((icntvl > 9 .AND. icntvl < 14) .AND. pedge /= 1) RETURN
    CALL gopen (parm,gplst(b2),2)

    RETURN
END SUBROUTINE contor
