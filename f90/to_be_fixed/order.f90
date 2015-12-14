SUBROUTINE order (gplst,id,rest,grids,idtab,lcor,b1,b2,b3)
     
 
 INTEGER, INTENT(IN OUT)                  :: gplst(1)
 INTEGER, INTENT(OUT)                     :: id(1)
 INTEGER, INTENT(OUT)                     :: rest(2)
 INTEGER, INTENT(IN)                      :: grids(1)
 INTEGER, INTENT(OUT)                     :: idtab(2)
 INTEGER, INTENT(IN)                      :: lcor
 INTEGER, INTENT(IN OUT)                  :: b1
 INTEGER, INTENT(IN OUT)                  :: b2
 INTEGER, INTENT(IN OUT)                  :: b3
 LOGICAL :: spill
 INTEGER :: tp, igrd(2),scr4, isym(14),itype(14),hold(3),elid,sils(34),  &
     est,sil,scr2,ecpt, three(3),offset
 COMMON /BLANK / ngp,skp(11),est,skip1(3),sil,skip2(5),ecpt,oes1,  &
     scr1,scr2,newoes,scr4
 EQUIVALENCE     (three(1),iflag),(three(2),nelmt),(three(3),igdpt)
 EQUIVALENCE     (kq4,isym(13)),(kt3,isym(14))
 DATA    isym  / 2HSH,2HT1,2HTB,2HTP,2HTM,2HQP,2HQM,2HT2,2HQ2,2HQ1,  &
     2HM1,2HM2,2HQ4,2HT3/    ,kbar/2HBR/
 DATA    itype / 4,6,7,8,9,15,16,17,18,19,62,63,64,83/
 DATA    ntype / 14 /
 
!     BUILD A TABLE FOR GPECT POINTERS TO ELID AND ITS ORDERED GRID PTS
 
 spill = .false.
 j = 1
 i = 3
 idtab(1) = 0
 jspill = 1
 lcorx  = lcor
 newin  = scr4
 newout = scr2
 2    CALL READ (*130,*12,est,tp,1,0,m)
 offset = 0
 IF (tp == kbar) offset = 6
 IF (tp == kt3 .OR. tp == kq4) offset = 1
 3    CALL fread (est,ngppe,1,0)
 idtab(i-1) = ngppe
 
!     SKIP PAST THE NON-CONTOUR ELEMENTS
 
 DO  k = 1,ntype
   IF (tp == isym(k)) GO TO 8
 END DO
 6    CALL fread (est,elid,1,0)
 IF (elid == 0) GO TO 2
 j = 1 + ngppe + offset
 CALL fread (est,0,-j,0)
 GO TO 6
 
!     CONSTRUCT IDTAB  1. 0, 2.NGPPE, 3.ELID, 4.ELIDPTR, 5.REPEAT.  3,4
!     FOR ALL ELEMENTS OF THIS TYPE, 6.REPEAT 1-5 FOR ALL ELEMENTS IN
!     THE SET. CONSTRUCT GRIDS  1-NGPPE. GRIDS FOR 1ST ELEMENT, NEXT.
!     REPEAT 1ST FOR ALL ELEMENTS IN THE IDTAB
 
 8    CALL READ (*12,*12,est,idtab(i),2,0,m)
 i = i + 2
 IF (idtab(i-2) /= 0) GO TO 10
 
!     END OF ELEMENTS OF THIS TYPE
 
 tp = idtab(i-1)
 GO TO 3
 
 10   CALL fread (est,grids(j),ngppe,0)
 IF (offset /= 0) CALL fread (est,0,-offset,0)
 j = j + ngppe
 IF (i >= lcorx) GO TO 14
 GO TO 8
 
!     TABLE FIT INTO CORE
 
 12   CALL bckrec (est)
 GO TO 16
 
!     SPILL OCCURS - TABLE DID NOT FIT
 
 14   spill = .true.
 
!     END OF TABLE
 
 16   lidtab = i - 1
 IF (lidtab <= 2) THEN
    SELECT CASE ( jspill )
     CASE (    1)
       GO TO 130
     CASE (    2)
       GO TO 140
   END SELECT
 END IF
 lgrids = j - 1
 lastng = ngppe
 SELECT CASE ( jspill )
   CASE (    1)
     GO TO 18
   CASE (    2)
     GO TO 140
 END SELECT
 18   CALL OPEN (*130,ecpt,gplst(b2),0)
 CALL gopen (scr2,gplst(b1),1)
 CALL fwdrec (*120,ecpt)
 igdpt = 0
 20   igdpt = igdpt + 1
 IEOR  = 0
 IF (gplst(igdpt) /= 0) GO TO 25
 CALL fwdrec (*120,ecpt)
 GO TO 20
 25   nelmt = 0
 iflag =-1
 
!      ECPT--1. PIVOT POINT, 2. DEG.FREEDOM, 3. -LENGTH, 4. ELID POINTER
!      5. ELTYPE, 6.SILS (THERE ARE (LENGTH-2) OF THEM), 7. REPEAT ITEMS
!      (3-6) FOR ALL ELEMENTS ATTACHED TO PIVOT, 8. EOR, 9. REPEAT ITEMS
!      (1-8) FOR ALL GRIDS IN THE PROBLEM.
 
 CALL READ (*120,*120,ecpt,igrd,2,0,m)
 30   CALL READ (*120,*75,ecpt,length,1,0,m)
 CALL fread (ecpt,sils,-length,0)
 tp = sils(2)
 DO  i = 1,ntype
   IF (tp == itype(i)) GO TO 33
 END DO
 GO TO 30
 
!     MATCH ELIDPTR WITH ITS ELID AND GRID POINTS IF POSSIBLE
 
 33   j = 1
 DO  i = 1,lidtab,2
   IF (idtab(i) == 0) THEN
     GO TO    35
   ELSE
     GO TO    40
   END IF
   35   ngppe = idtab(i+1)
   CYCLE
   40   IF (idtab(i+1) == sils(1)) GO TO 55
   j = j + ngppe
 END DO
 
!     IF NOT IN THE TABLE, IS THERE SPILL(IE IS TABLE NOT COMPLETE).
!     NO SPILL, SKIP HIM.  YES SPILL, FLAG HIM.
 
 IF (.NOT.spill) THEN
    SELECT CASE ( jspill )
     CASE (    1)
       GO TO 30
     CASE (    2)
       GO TO 145
   END SELECT
 END IF
 elid  =-sils(1)
 nelmt = nelmt + 1
 GO TO 70
 
!     FOUND ELEMENT IN THE TABLE
 
 55   elid = idtab(i)
 DO  i = 1,ngppe
   k = j + i - 1
   IF (igdpt == grids(k)) EXIT
 END DO
 65   iafter = i - (i/ngppe)*ngppe + j
 ibefor = j + i - 2
 IF (i == 1) ibefor = ibefor + ngppe
 nelmt = nelmt + 1
 rest(2*nelmt-1) = grids(iafter)
 rest(2*nelmt  ) = grids(ibefor)
 70   id(nelmt) = elid
 IF (nelmt < lcor/2) THEN
    SELECT CASE ( jspill )
     CASE (    1)
       GO TO 30
     CASE (    2)
       GO TO 145
   END SELECT
 END IF
 GO TO 80
 75   IF (nelmt == 0) SELECT CASE ( jspill )
   CASE (    1)
     GO TO 20
   CASE (    2)
     GO TO 140
 END SELECT
 IEOR = 1
 
!     ORDER ELEMENTS IF WE HAVE REACHED END OF EST FILE
 
 80   IF (spill) GO TO 112
 IF (nelmt <= 2)  GO TO 110
 INDEX = 3
 iall  = 2*nelmt
 ione  = rest(1)
 itwo  = rest(2)
 85   IF (ione == itwo) GO TO 105
 DO  i = INDEX,iall,2
   IF (itwo == rest(i)) GO TO 95
 END DO
 GO TO 110
 95   IF (i == INDEX) GO TO 100
 j = (INDEX+1)/2
 k = (i+1)/2
 hold(1) = id(j)
 id(j)   = id(k)
 id(k)   = hold(1)
 hold(2) = rest(INDEX  )
 hold(3) = rest(INDEX+1)
 rest(INDEX  ) = rest(i  )
 rest(INDEX+1) = rest(i+1)
 rest(i  ) = hold(2)
 rest(i+1) = hold(3)
 100  INDEX = INDEX + 2
 itwo  = rest(INDEX-1)
 IF (INDEX < iall) GO TO 85
 IF (ione  /= itwo) GO TO 110
 
!     INTERIOR ELEMENTS
 
 105  CALL WRITE (newout,three,3,0)
 CALL WRITE (newout,id,nelmt,1)
 IF (igdpt < ngp) THEN
    SELECT CASE ( jspill )
     CASE (    1)
       GO TO 20
     CASE (    2)
       GO TO 140
   END SELECT
 END IF
 GO TO 120
 
!     BORDER ELEMENTS
 
 110  iflag = -2
 112  CALL WRITE (newout,three,3,0)
 j = -1
 DO  i = 1,nelmt
   j = j + 2
   CALL WRITE (newout,id(i),1,0)
   CALL WRITE (newout,rest(j),2,0)
 END DO
 iq = 2*nelmt
 CALL WRITE (newout,0,0,1)
 SELECT CASE ( jspill )
   CASE (    1)
     GO TO 118
   CASE (    2)
     GO TO 140
 END SELECT
 118  IF (IEOR) 25,25,119
 119  IF (igdpt < ngp) GO TO 20
 120  CALL CLOSE (ecpt,1)
 125  CALL WRITE (newout,0,1,1)
 CALL CLOSE (newout,1)
 
!     IF NO SPILL -  RETURN
 
 130  IF (.NOT.spill) GO TO 170
 
!    COME HERE IF WE HAVE SPILL
 
 i = newout
 newout = newin
 newin  = i
 CALL gopen (newin,gplst(b1),0)
 CALL gopen (newout,gplst(b2),1)
 jspill = 2
 ngppe  = lastng
 idtab(1) = 0
 idtab(2) = lastng
 i = 3
 j = 1
 spill = .false.
 GO TO 8
 
!     TABLE CONSTRUCTED SO RETURN HERE
 
 140  CALL READ (*160,*160,newin,three,3,0,m)
 nelmt = 0
 145  CALL READ (*160,*75,newin,sils(1),3,0,m)
 IF (sils(1) > 0.0) THEN
   GO TO   155
 END IF
 150  sils(1) = -sils(1)
 GO TO 33
 155  elid  = sils(1)
 nelmt = nelmt + 1
 rest(2*nelmt-1) = sils(2)
 rest(2*nelmt  ) = sils(3)
 GO TO 70
 
!     END OF FILE
 
 160  CALL CLOSE (newin,1)
 GO TO 125
 
!     OUTPUT FILE MUST BE SCRATCH 2
 
 170  IF (newout == scr2) RETURN
 CALL gopen  (newout,gplst(b1),0)
 CALL gopen  (scr2,gplst(b2),1)
 CALL cpyfil (newout,scr2,rest,lcor,m)
 CALL CLOSE  (scr2,1)
 CALL CLOSE  (newout,1)
 RETURN
END SUBROUTINE order
