SUBROUTINE rcovds
     
!     THIS ROUTINE GENERATES THE DYNAMIC SOLUTION ITEM FOR RIGID
!     FORMATS 8 AND 9
 
 INTEGER :: dry        ,step       ,fss        ,rfno        ,  &
            rd         ,rdrew      ,wrt        ,wrtrew      ,  &
            rew        ,sysbuf     ,rc         ,eqss        ,  &
            soln       ,srd        ,swrt       ,eoi         ,  &
            eog        ,iz(5)      ,upv        ,trl(7)      ,  &
            buf1       ,dload      ,dlt        ,casess      ,  &
            geom4      ,loadc(2)   ,tolppf     ,NAME(2)     ,  &
            FILE       ,dit        ,tabloc(13) ,casecc(2)
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm        ,uwm
 COMMON /BLANK / dry        ,loop       ,step       ,fss(2)      ,  &
                 rfno       ,neigv      ,lui        ,uinms(2,5)  ,  &
                 nosort     ,uthres     ,pthres     ,qthres
 COMMON /rcovcr/ icore      ,lcore      ,buf1       ,buf2        ,  &
                 buf3       ,buf4       ,sof1       ,sof2        , sof3
 COMMON /rcovcm/ mrecvr     ,ua         ,pa         ,qa          ,  &
                 iopt       ,rss(2)     ,energy     ,uimpro      ,  &
                 range(2)   ,ireq       ,lreq       ,lbasic
 COMMON /names / rd         ,rdrew      ,wrt        ,wrtrew      ,  &
                 rew        ,norew      ,eofnrw
 COMMON /zzzzzz/ z(1)
 COMMON /system/ sysbuf     ,nout
 COMMON /condas/ pi         ,twopi      ,radeg      ,degra
 EQUIVALENCE     (z(1),iz(1)) ,  (iscale,scale)
 DATA    NAME  / 4HRCOV,4HDS  /
 DATA    eqss  , soln,lods    / 4HEQSS,4HSOLN,4HLODS /
 DATA    srd   , swrt,eog,eoi / 1,2,2,3      /
 DATA    upv   , dlt,casess,geom4,tolppf,dit /  &
         106   , 108,101   ,102  ,111   ,107 /
 DATA    loadc / 500,  5      /
 DATA    tabloc/ 4,1105,11,1,1205,12,2,1305,13,3,1405,14,4 /
 DATA    casecc/ 4HCASE,4HCC  /
 
!     CREATE SOLN FOR RIGID FORMAT 8 OR 9
 
!     GET NUMBER OF BASIC SUBSTRUCTURES (NS) FROM EQSS AND CREATE
!     GROUP 0 OF SOLN AT TOP OF OPEN CORE
 
 lcore = buf1 - 1
 CALL sfetch (fss,eqss,srd,rc)
 IF (rc == 1) GO TO 110
 CALL smsg (rc-2,eqss,fss)
 GO TO 810
 110 CALL suread (z,2,nwds,rc)
 CALL suread (ns,1,nwds,rc)
 IF (lcore < 2*ns+5) GO TO 9008
 CALL suread (z,1,nwds,rc)
 iz(1) = fss(1)
 iz(2) = fss(2)
 iz(3) = rfno
 iz(4) = ns
 
!     GET THE BASIC SUBSTRUCTURE NAMES FROM EQSS
 
 DO  i = 1,ns
   CALL suread (z(3*i+3),2,nwds,rc)
 END DO
 
!     GET THE NUMBER OF LOAD VECTORS FOR EACH SUBSTRUCTURE FORM LODS
 
 CALL sfetch (fss,lods,srd,rc)
 IF (rc == 1) GO TO 160
 CALL smsg (rc-2,lods,fss)
 GO TO 9200
 160 j = 1
 CALL sjump (j)
 DO  i = 1,ns
   CALL suread (z(3*i+5),1,nwds,rc)
   CALL sjump  (j)
 END DO
 
!     GET THE NUMBER OF TIME OR FREQUENCY STEPS FROM UPV OR UPVC
 
 trl(1) = upv
 CALL rdtrl (trl)
 nstep = trl(2)
 IF (rfno == 9) nstep = nstep/3
 iz(5) = nstep
 
!     GET THE REQUESTED DLOAD SET FROM CASE CONTROL
 
 FILE = casess
 CALL gopen (casess,z(buf1),rdrew)
 180 CALL fread (casess,trl,2,1)
 IF (trl(1) /= casecc(1) .OR. trl(2) /= casecc(2)) GO TO 180
 CALL fread (casess,0,-12,0)
 CALL fread (casess,dload,1,0)
 CALL CLOSE (casess,rew)
 
!     CHECK IF DLOAD SET POINTS TO A DLOAD COMBINATION CARD OR A
!     SIMPLE LOAD CARD BY LOOKING AT SET IDS IN HEADER RECORD OF DLT
 
 i = 3*ns + 6
 FILE = dlt
 CALL OPEN (*9001,dlt,z(buf1),rdrew)
 CALL READ (*9002,*200,dlt,z(i),lcore-i,1,nwds)
 GO TO 9008
 200 idlset = i + 3
 ldlset = idlset + iz(i+2) - 1
 ildset = ldlset + 1
 lldset = i + nwds - 1
 idload = lldset + 1
 IF (idlset > ldlset) GO TO 215
 DO  i = idlset,ldlset
   IF (iz(i) == dload) GO TO 220
 END DO
 
!     NO DLOAD MATCH - MUST BE SIMPLE RLOAD OR TLOAD
 
 215 z(idload) = 1.0
 iz(idload+1) = dload
 ldload = idload + 1
 GO TO 270
 
!     DLOAD MATCH FOUND - READ DLOAD DATA FORM DLT RECORD 1
 
 220 CALL fread (dlt,trl,2,0)
 IF (trl(1) == dload) GO TO 240
 230 CALL fread (dlt,trl,2,0)
 IF (trl(1) /= -1) GO TO 230
 GO TO 220
 240 i = idload
 iscale = trl(2)
 250 CALL fread (dlt,z(i),2,0)
 IF (iz(i) == -1) GO TO 260
 z(i) = z(i)*scale
 i = i + 2
 IF (i > lcore) GO TO 9008
 GO TO 250
 260 ldload = i - 1
 
!     READ THE RLOAD AND TLOAD DATA FORM DLT AND SAVE REQUESTED CARDS
 
 270 iload = ldload + 1
 l = iload
 IF (idlset <= ldlset) CALL fwdrec (*9002,dlt)
 DO  i = ildset,lldset
   DO  j = idload,ldload,2
     IF (iz(j+1) == iz(i)) GO TO 290
   END DO
   CALL fwdrec (*9002,dlt)
   CYCLE
   
!     SAVE RLOAD DATA IF RIGID FORMAT 8
!     SAVE TLOAD DATA IF RIGID FORMAT 9
   
   290 CALL fread (dlt,itype,1,0)
   IF (itype <= 2 .AND. rfno == 8) GO TO 300
   IF (itype >= 3 .AND. rfno == 9) GO TO 300
   CALL fwdrec (*9002,dlt)
   CYCLE
   
   300 iz(l) = itype
   CALL fread (dlt,iz(l+1),7,1)
   iz(j+1) = -l
   l = l + 8
   IF (l > lcore) GO TO 9008
 END DO
 
 lload = l - 1
 CALL CLOSE (dlt,rew)
 
!     READ THE LOADC DATA FROM GEOM4 AND SAVE ANY THAT WAS REQUESTED
!     ON TLOAD OR RLOAD CARDS
 
!     NOTE - UNTIL A MODULE FRLG IS WRITTEN NO RLOAD CARD MAY REQUEST A
!            SCALAR LOAD
 
 nsload = 0
 iloadc = lload  + 1
 lloadc = iloadc - 1
 isload = iloadc
 lsload = isload - 1
 
 IF (rfno == 8) GO TO 500
 
 CALL preloc (*500,z(buf1),geom4)
 CALL locate (*500,z(buf1),loadc,i)
 iold = 0
 i1 = iloadc
 i2 = i1
 320 CALL READ (*9002,*370,geom4,trl(1),2,0,nwds)
 iscale = trl(2)
 IF (iold == trl(1)) GO TO 360
 ihit = 0
 DO  i = iload,lload,8
   IF (trl(1) /= iz(i+1)) CYCLE
   iz(i+1) = -i1
   ihit = ihit + 1
 END DO
 IF (ihit > 0) GO TO 350
 340 CALL fread (geom4,trl(1),4,0)
 IF (trl(3) /= -1) GO TO 340
 GO TO 320
 
!     THIS LOADC DATA WAS REQUESTED - SAVE THE DATA AND A POINTER TO IT
 
 350 iold = trl(1)
 i1 = i2
 iz(i1) = 0
 i2 = i1 + 1
 360 CALL fread (geom4,z(i2),4,0)
 IF (iz(i2+2) == -1) GO TO 320
 iz(i1) = iz(i1) + 1
 z(i2+3) = z(i2+3)*scale
 i2 = i2 + 4
 IF (i2 > lcore) GO TO 9008
 GO TO 360
 
!     CONVERT LOADC LOAD SETS TO INTERNAL LOAD IDS BY USING THE LODS
!     ITEM
 
 370 lloadc = i2 - 1
 IF (iloadc > lloadc) GO TO 500
 CALL sfetch (fss,lods,srd,rc)
 i = 1
 CALL sjump (i)
 ilod = 1
 idat0= lloadc + 1
 idat = idat0  + 1
 ndat = lcore  - lloadc
 isub = 6
 lsub = 3*ns + 5
 
!     FOR EACH BASIC READ THE LODS DATA INTO CORE
 
 DO  i = isub,lsub,3
   CALL suread (z(idat0),ndat,nwds,rc)
   IF (rc /= 2) GO TO 9008
   j = iloadc
   380 i1 = j + 1
   i2 = j + iz(j)*4
   DO  k = i1,i2,4
     IF (iz(k) /= iz(i) .OR. iz(k+1) /= iz(i+1)) CYCLE
     
!     FOUND LOADC DATA FOR THIS BASIC - CONVERT LOAD SET ID
     
     iz(k  ) = 0
     iz(k+1) = 0
     nwds = idat0 + nwds - 1
     DO  l = idat,nwds
       IF (iz(l) == iz(k+2)) GO TO 395
     END DO
     WRITE (nout,6316) uwm,iz(k+2),z(i),z(i+1),fss
     iz(k+2) = -1
     CYCLE
     
     395 iz(k+2) = ilod + l - idat
     
   END DO
   j = i2 + 1
   IF (j < lloadc) GO TO 380
   
   ilod = ilod + iz(idat0)
 END DO
 
!     CREATE A LIST OF INTERNAL LOAD VECTORS REQUESTED - ALSO CHECK IF
!     ANY BASIC NAMES WERE NOT FOUND
 
 isload = lloadc + 1
 lsload = isload - 1
 nsload = 0
 j  = iloadc
 420 i1 = j + 1
 i2 = j + iz(j)*4
 outer_loop:  DO  k = i1,i2,4
   IF (iz(k) == 0) GO TO 430
   WRITE (nout,6315) uwm,z(k),z(k+1),fss,iz(k+2),fss
   iz(k+2) = -1
   CYCLE outer_loop
   430 IF (iz(k+2) < 0) CYCLE outer_loop
   IF (nsload  == 0) GO TO 455
   DO  i = isload,lsload
     IF (iz(i) == iz(k+2)) CYCLE outer_loop
   END DO
   455 nsload = nsload + 1
   lsload = lsload + 1
   IF (lsload > lcore) GO TO 9008
   iz(lsload) = iz(k+2)
 END DO outer_loop
 j = i2 + 1
 IF (j < lloadc) GO TO 420
 
!     SORT LIST OF IDS
 
 CALL sort (0,0,1,1,z(isload),nsload)
 
!     MAKE ONE MORE PASS THROUGH THE LOAC DATA CONVERTING THE
!     INTERNAL LOAD IDS TO A RELATIVE POSITION IN THE LOAD LIST
!     STARTING AT ISLOAD
 
 j  = iloadc
 470 i1 = j + 1
 i2 = j + iz(j)*4
 DO  k = i1,i2,4
   IF (iz(k+2) < 0) CYCLE
   DO  l = isload,lsload
     IF (iz(k+2) == iz(l)) GO TO 490
   END DO
   CYCLE
   490 iz(k+2) = l - isload
 END DO
 j = i2 + 1
 IF (j < lloadc) GO TO 470
 
!     OK - NOW WE CAN WRITE OUT GROUP 0 OF THE SOLN ITEM
 
 500 CALL CLOSE (geom4,rew)
 rc = 3
 CALL sfetch (fss,soln,swrt,rc)
 CALL suwrt (z(1),3*ns+5,1)
 CALL suwrt (nsload,1,1)
 IF (nsload > 0) CALL suwrt (z(isload),nsload,1)
 CALL suwrt (0,0,eog)
 
!     COPY THE FREQUENCY STEPS FROM PPF OR THE TIME STEPS FROM TOL
!     FOR GROUP 1 OF THE SOLN ITEM
 
 istep = isload
 lstep = istep + nstep - 1
 IF (lstep > lcore) GO TO 9008
 FILE = tolppf
 CALL OPEN (*9001,tolppf,z(buf1),rdrew)
 CALL fread (tolppf,trl,2,0)
 CALL fread (tolppf,z(istep),nstep,0)
 CALL CLOSE (tolppf,rew)
 
 CALL suwrt (z(istep),nstep,eog)
 
!     IF ANY SCALAR LOADS EXIST CALCULATE THE SCALE FACTORS FOR EACH
!     LOAD AND WRITE THEM TO THE SOF - 1 GROUP PER TIME OR FREQUENCY
!     STEP
 
 IF (nsload == 0) GO TO 800
 ivec = lstep + 1
 lvec = ivec + nsload - 1
 IF (lvec > lcore) GO TO 9008
 
!     CALL PRETAB TO READ IN THE REQUIRED TABLE DATA - FIRST MAKE A
!     LIST OF REQUESTED TABLE IDS
 
 itab0 = lvec + 1
 iz(itab0) = 0
 itab = itab0 + 1
 ltab = itab  - 1
 DO  j = iload,lload,8
   IF (iz(j+1) >= 0) CYCLE
   itype = iz(j)
   SELECT CASE ( itype )
     CASE (    1)
       GO TO 510
     CASE (    2)
       GO TO 510
     CASE (    3)
       GO TO 520
     CASE (    4)
       GO TO 570
   END SELECT
   510 i1 = j + 2
   i2 = j + 3
   GO TO 530
   520 i1 = j + 2
   i2 = j + 2
530  outer_loop2: DO  k = i1,i2
     IF (iz(k) ==   0) CYCLE outer_loop2
     IF (ltab < itab) GO TO 550
     DO  l = itab,ltab
       IF (iz(l) == iz(k)) CYCLE outer_loop2
     END DO
     550 ltab = ltab + 1
     IF (ltab > lcore) GO TO 9008
     iz(ltab ) = iz(k)
     iz(itab0) = iz(itab0) + 1
   END DO outer_loop2
   570 CONTINUE
 END DO
 
 IF (iz(itab0) == 0) GO TO 585
 itabd = ltab + 1
 CALL pretab (dit,z(itabd),iz(itabd),z(buf1),lcore-itabd,ltabd,  &
     z(itab0),tabloc)
 ltabd = itabd + ltabd - 1
 585 CONTINUE
 
!     LOOP OVER EACH TIME OR FREQUENCY STEP
 
 DO  i = istep,lstep
   
!     ZERO A VECTOR IN CORE FOR THE SCALE FACTORS
   
   DO  j = ivec,lvec
     iz(j) = 0
   END DO
   
!     PASS THROUGH THE DLOAD DATA
   
   DO  j = idload,ldload,2
     IF (iz(j+1) >= 0) CYCLE
     
!     PROCESS THE TLOAD OR RLOAD DATA THIS DLOAD ENTRY POINTS TO
     
     ild = -iz(j+1)
     IF (iz(ild+1) >= 0) CYCLE
     itype = iz(ild  )
     ildc  =-iz(ild+1)
     
!     CALCULATE THE SCALE FACTOR FOR THE CARD FOR THIS TIME OR FREQUENCY
!     STEP
     
     SELECT CASE ( itype )
       CASE (    1)
         GO TO 600
       CASE (    2)
         GO TO 640
       CASE (    3)
         GO TO 680
       CASE (    4)
         GO TO 720
     END SELECT
     
!     RLOAD1 DATA
     
     600 scale = 0.0
     GO TO 760
     
!     RLOAD2 DATA
     
     640 scale = 0.0
     GO TO 760
     
!     TLOAD1 DATA
     
     680 CALL tab (iz(ild+2),z(i),scale)
     GO TO 760
     
!     TLOAD2 DATA
     
     720 scale = 0.0
     tt    = z(i) - z(ild+2)
     IF (tt == 0.0) GO TO 730
     IF (tt < 0.0 .OR. tt > z(ild+3)) GO TO 760
     scale = tt**z(ild+7)*EXP(z(ild+6)*tt)*COS(twopi*z(ild+4)*tt  &
         + z(ild+5)*degra)
     GO TO 760
     730 IF (z(ild+7) /= 0.0) GO TO 760
     scale = COS(z(ild+5))
     
!     NOW APPLY THIS SCALE FACTOR TO EACH LOADC ENTRY.
!     TOTAL SCALE FACTOR = T(R)LOAD FACTOR*DLOAD FACTOR*LOADC FACTOR
     
     760 CONTINUE
     IF (scale == 0.0) CYCLE
     i1 = ildc + 1
     i2 = ildc + iz(ildc)*4
     DO  k = i1,i2,4
       IF (iz(k+2) < 0) CYCLE
       ifac = ivec + iz(k+2)
       z(ifac) = z(ifac) + scale*z(j)*z(k+3)
     END DO
     
   END DO
   
!     WRITE OUT THESE FACTORS TO THE NEXT GROUP OF THE SOF
   
   CALL suwrt (z(ivec),nsload,eog)
 END DO
 
!     FINISHED
 
 800 CALL suwrt (0,0,eoi)
 810 CALL sofcls
 RETURN
 
!     DIAGNOSTICS
 
 6315 FORMAT (a25,' 6315, RCOVR MODULE - SUBSTRUCTURE ',2A4,' IS NOT A',  &
                  ' COMPONENT OF ',2A4, /32X,'LOAD SET',i9,' FOR THAT ',  &
                  'SUBSTRUCTURE WILL BE IGNORED IN CREATING', /32X,  &
                  'THE SOLN ITEM FOR FINAL SOLUTION STRUCTURE ',2A4)
 6316 FORMAT (a25,' 6316, RCOVR MODULE IS UNABLE TO FIND LOAD SET ',i8,  &
                  ' FOR SUBSTRUCTURE ',2A4, /32X,'AMONG THOSE ON LODS.  ',  &
                  'IT WILL BE IGNORED IN CREATING THE SOLN ITEM FOR FINAL',  &
             /32X,'SOLUTION STRUCTURE ',2A4)
             
 9001 n = 1
 GO TO 9100
 9002 n = 2
 GO TO 9100
 9008 n = 8
 9100 CALL mesage (n,FILE,NAME)
 9200 CALL sofcls
 iopt = -1
 CALL CLOSE (casess,rew)
 CALL CLOSE (dlt,rew)
 CALL CLOSE (geom4,rew)
 CALL CLOSE (tolppf,rew)
 CALL CLOSE (dit,rew)
 
 RETURN
END SUBROUTINE rcovds
