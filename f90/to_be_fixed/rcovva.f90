SUBROUTINE rcovva (in,intyp,outt,outu,outv,outa,ssnm,rz,dz,cz)
     
!     THIS SUBROUTINE COMPUTES THE VELOCITIES AND ACCELERATIONS FOR
!     FOR A GIVEN DISPLACEMENT VECTOR
 
!     INTYP = 0   IN CONTAINS U ONLY AND V AND A ARE CALCULATED
!     INTYP = 1   U CONTAINS U, V AND A SO THEY ARE SPLIT ONTO OUTU,
!                 OUTV AND OUTA
!     INTYP =-1   OUTU, OUTV AND OUTA ARE MERGED ONTO OUTT
 
 
 INTEGER, INTENT(IN OUT)                  :: in
 INTEGER, INTENT(IN)                      :: intyp
 INTEGER, INTENT(IN)                      :: outt
 INTEGER, INTENT(IN)                      :: outu
 INTEGER, INTENT(IN)                      :: outv
 INTEGER, INTENT(IN)                      :: outa
 INTEGER, INTENT(IN OUT)                  :: ssnm(2)
 REAL, INTENT(IN OUT)                     :: rz(4)
 DOUBLE PRECISION, INTENT(OUT)            :: dz(1)
 COMPLEX, INTENT(IN OUT)                  :: cz(2)
 INTEGER :: dry        ,rss        ,ua         ,NAME(2)    ,  &
     rc         ,rfno       ,mcbu(7)    ,mcbv(7)    ,  &
     srd        ,soln       ,buf1       ,buf2       ,  &
     buf3       ,buf4       ,inblk3(15) ,outblk(15) , sysbuf     ,  &
      mcb(7)     ,mcba(7)    ,inblk(15)  ,  &
     oblk1(15)  ,oblk2(15)  ,oblk3(15)  ,temp(4)    ,  &
     FILE       , inblk1(15) ,inblk2(15)
 REAL :: freq       ,rscale     ,iscale
 DOUBLE PRECISION :: rval       ,ival
 COMPLEX :: scale
 CHARACTER (LEN=1) ::
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm swm*27
 COMMON /xmssg /  ufm        ,uwm        ,uim        ,sfm        , swm
 COMMON /BLANK /  dry        ,loop       ,step       ,fss(2)     ,  &
     rfno       ,neigv      ,lui        ,uinms(2,5) ,  &
     nosort     ,uthres     ,pthres     ,qthres
 COMMON /rcovcr/  icore      ,lcore      ,buf1       ,buf2       ,  &
     buf3       ,buf4       ,sof1       ,sof2       , sof3
 COMMON /rcovcm/  mrecvr     ,ua         ,pa         ,qa         ,  &
     iopt       ,rss(2)     ,energy     ,uimpro     ,  &
     range(2)   ,ireq       ,lreq       ,lbasic
 COMMON /system/  sysbuf     ,nout
 COMMON /names /  rd         ,rdrew      ,wrt        ,wrtrew     ,  &
     rew        ,norew      ,eofnrw     ,rsp        ,  &
     rdp        ,csp        ,cdp        ,square     ,  &
     rect       ,diag       ,upper      ,lower      , sym
 COMMON /condas/  phi        ,twophi
 COMMON /packx /  itinp      ,itoutp     ,irp        ,nrp        , incrp
 COMMON /unpakx/  itinu      ,iru        ,nru        ,incru
 COMMON /TYPE  /  pr(2)      ,nwords(4)
 EQUIVALENCE      (temp(1),scale,rscale) ,(temp(2),iscale)       ,  &
     (inblk1(1),oblk1(1))   ,(inblk2(1),oblk2(1))   ,  &
     (inblk3(1),oblk3(1))   ,(outblk(1),inblk(1))
 DATA    NAME  /  4HRCOV,4HVA   /
 DATA    srd   /  1             /
 DATA    soln  /  4HSOLN        /
 
!     GET DISPLACEMENT TRAILER AND DETERMINE TYPE
 
 IF (outt /= 0 .AND. outu+outv+outa /= 0 .AND. intyp >= 0) GO TO 9007
 
 FILE = in
 IF (intyp < 0) FILE = outu
 mcb(1) = FILE
 CALL rdtrl (mcb)
 IF (mcb(1) < 0) GO TO 9001
 ncol  = mcb(2)
 nrow  = mcb(3)
 iprec = mcb(5)
 nword = nwords(iprec)
 nwcol = nrow*nword
 
!     SET UP PACK UNPACK COMMONS
 
 itinu = iprec
 iru   = 1
 nru   = nrow
 incru = 1
 itinp = iprec
 itoutp= iprec
 irp   = 1
 nrp   = nrow
 incrp = 1
 
!     BRANCH ON TYPE OF DISPLACEMENTS OR RIGID FORMAT
 
 IF (intyp > 0) GO TO 400
 IF (intyp < 0) GO TO 500
 
 IF (rfno > 9) GO TO 9007
 SELECT CASE ( rfno )
   CASE (    1)
     GO TO 600
   CASE (    2)
     GO TO 600
   CASE (    3)
     GO TO 100
   CASE (    4)
     GO TO 9007
   CASE (    5)
     GO TO 9007
   CASE (    6)
     GO TO 9007
   CASE (    7)
     GO TO 9007
   CASE (    8)
     GO TO 200
   CASE (    9)
     GO TO 400
 END SELECT
 
!     NORMAL MODES
 
!     CHECK IF VECTORS ARE COMPLEX
 
 100 IF (iprec >= 3) GO TO 200
 
!     REAL NORMAL MODES
 
!     V =  U*OMEGA
!     A = -V*OMEGA
 
 IF (lcore < nwcol) GO TO 6313
 item = soln
 CALL sfetch (ssnm,soln,srd,rc)
 IF (rc /= 1) GO TO 6000
 n = 1
 CALL sjump (n)
 IF (n < 0) GO TO 6100
 
 CALL gopen (in,rz(buf1),rdrew)
 IF (outt /= 0) CALL gopen (outt,rz(buf2),wrtrew)
 IF (outu /= 0) CALL gopen (outu,rz(buf2),wrtrew)
 IF (outv /= 0) CALL gopen (outv,rz(buf3),wrtrew)
 IF (outa /= 0) CALL gopen (outa,rz(buf4),wrtrew)
 CALL makmcb (mcb,outt,nrow,rect,iprec)
 CALL makmcb (mcbu,outu,nrow,rect,iprec)
 CALL makmcb (mcbv,outv,nrow,rect,iprec)
 CALL makmcb (mcba,outa,nrow,rect,iprec)
 
!     LOOP THROUGH EACH COLUMN
 
 DO  i = 1,ncol
   
!     GET SCALE FACTOR FOR THIS COLUMN
   
   CALL suread (rz,7,nwds,rc)
   IF (rc /= 1) GO TO 6200
   rscale = rz(4)
   
   CALL unpack (*110,in,rz)
   GO TO 120
   110 DO  j = 1,nwcol
     rz(j) = 0.0
   END DO
   120 IF (outt /= 0) CALL pack (rz(1),outt,mcb)
   IF (outu /= 0) CALL pack (rz(1),outu,mcbu)
   
   DO  j = 1,2
     IF (iprec == 2) GO TO 140
     
     DO  k = 1,nrow
       rz(k) = rscale*rz(k)
     END DO
     GO TO 160
     
     140 DO  k = 1,nrow
       dz(k) = rscale*dz(k)
     END DO
     
     160 IF (outt /= 0) CALL pack (rz(1),outt,mcb)
     IF (outv /= 0 .AND. j == 1) CALL pack (rz(1),outv,mcbv)
     IF (outa /= 0 .AND. j == 2) CALL pack (rz(1),outa,mcba)
     
     rscale = -rscale
   END DO
 END DO
 
 
 CALL CLOSE (in,rew)
 IF (outt /= 0) CALL CLOSE (outt,rew)
 IF (outu /= 0) CALL CLOSE (outu,rew)
 IF (outv /= 0) CALL CLOSE (outv,rew)
 IF (outa /= 0) CALL CLOSE (outa,rew)
 IF (outt /= 0) CALL wrttrl (mcb)
 IF (outu /= 0) CALL wrttrl (mcbu)
 IF (outv /= 0) CALL wrttrl (mcbv)
 IF (outa /= 0) CALL wrttrl (mcba)
 GO TO 600
 
!     COMPLEX NORMAL MODES
 
!     V = U*POLE
!     A = V*POLE
 
!     FREQUENCY RESPONSE
 
!     V = U*TWOPHI*FREQ*I
!     A = V*TWOPHI*FREQ*I
 
 200 IF (lcore < nwcol) GO TO 6313
 item = soln
 CALL sfetch (ssnm,soln,srd,rc)
 IF (rc /= 1) GO TO 6000
 n = 1
 CALL sjump (n)
 IF (n < 0) GO TO 6100
 
 CALL gopen (in,rz(buf1),rdrew)
 IF (outt /= 0) CALL gopen (outt,rz(buf2),wrtrew)
 IF (outu /= 0) CALL gopen (outu,rz(buf2),wrtrew)
 IF (outv /= 0) CALL gopen (outv,rz(buf3),wrtrew)
 IF (outa /= 0) CALL gopen (outa,rz(buf4),wrtrew)
 CALL makmcb (mcb ,outt,nrow,rect,iprec)
 CALL makmcb (mcbu,outu,nrow,rect,iprec)
 CALL makmcb (mcbv,outv,nrow,rect,iprec)
 CALL makmcb (mcba,outa,nrow,rect,iprec)
 
!     LOOP THROUGH EACH COLUMN
 
 DO  i = 1,ncol
   
!     GET SCALE FACTOR FOR THIS COLUMN
   
   IF (rfno == 8) GO TO 204
   CALL suread (cz(1),7,nwds,rc)
   IF (rc /= 1) GO TO 6200
   scale = cz(2)
   GO TO 206
   204 CALL suread (freq,1,nwds,rc)
   IF (rc /= 1) GO TO 6200
   scale = twophi*freq*(0.0,1.0)
   
   206 CALL unpack (*210,in,cz(1))
   GO TO 230
   210 DO  j = 1,nwcol
     rz(j) = 0.0
   END DO
   230 IF (outt /= 0) CALL pack (cz(1),outt,mcb)
   IF (outu /= 0) CALL pack (cz(1),outu,mcbu)
   
   DO  j = 1,2
     IF (iprec > 3) GO TO 250
     
     DO  k = 1,nrow
       cz(k) = scale*cz(k)
     END DO
     GO TO 270
     
     250 nt = nrow*2
     DO  k = 1,nt,2
       rval = dz(k  )
       ival = dz(k+1)
       dz(k  ) = rscale*rval - iscale*ival
       dz(k+1) = rscale*ival + iscale*rval
     END DO
     
     270 IF (outt /= 0) CALL pack (cz(1),outt,mcb)
     IF (outv /= 0 .AND. j == 1) CALL pack (cz(1),outv,mcbv)
     IF (outa /= 0 .AND. j == 2) CALL pack (cz(1),outa,mcba)
     
   END DO
   
 END DO
 
 CALL CLOSE (in,rew)
 IF (outt /= 0) CALL CLOSE (outt,rew)
 IF (outu /= 0) CALL CLOSE (outu,rew)
 IF (outv /= 0) CALL CLOSE (outv,rew)
 IF (outa /= 0) CALL CLOSE (outa,rew)
 IF (outt /= 0) CALL wrttrl (mcb)
 IF (outu /= 0) CALL wrttrl (mcbu)
 IF (outv /= 0) CALL wrttrl (mcbv)
 IF (outa /= 0) CALL wrttrl (mcba)
 GO TO 600
 
!     THE DISPLACEMENT FILE ALREADY CONTAINS THE VELOCITIES AND
!     ACCELERATIONS SO WE JUST SANT TO SPLIT THEM UP
 
 400 IF (lcore < 0) GO TO 6313
 CALL gopen (in,rz(buf1),rdrew)
 IF (outu /= 0) CALL gopen (outu,rz(buf2),wrtrew)
 IF (outv /= 0) CALL gopen (outv,rz(buf3),wrtrew)
 IF (outa /= 0) CALL gopen (outa,rz(buf4),wrtrew)
 
 inblk(1) = in
 oblk1(1) = outu
 oblk2(1) = outv
 oblk3(1) = outa
 FILE = in
 ncol = ncol/3
 
 DO  i = 1,ncol
   IF (outu /= 0) CALL cpystr (inblk,oblk1,0,i)
   IF (outu == 0) CALL fwdrec (*9002,in)
   IF (outv /= 0) CALL cpystr (inblk,oblk2,0,i)
   IF (outv == 0) CALL fwdrec (*9002,in)
   IF (outa /= 0) CALL cpystr (inblk,oblk3,0,i)
   IF (outa == 0) CALL fwdrec (*9002,in)
 END DO
 
 CALL CLOSE (in,rew)
 IF (outu /= 0) CALL CLOSE (outu,rew)
 IF (outv /= 0) CALL CLOSE (outv,rew)
 IF (outa /= 0) CALL CLOSE (outa,rew)
 mcb(2) = ncol
 mcb(1) = outu
 IF (outu /= 0) CALL wrttrl (mcb)
 mcb(1) = outv
 IF (outv /= 0) CALL wrttrl (mcb)
 mcb(1) = outa
 IF (outa /= 0) CALL wrttrl (mcb)
 GO TO 600
 
!     THE DISPLACEMENTS, VELOCITIES AND ACCLERATIONS ALREADY EXIST AND
!     ARE TO BE MERGED TOGETHER
 
 500 IF (lcore < 0) GO TO 6313
 CALL gopen (outu,rz(buf1),rdrew)
 CALL gopen (outv,rz(buf2),rdrew)
 CALL gopen (outa,rz(buf3),rdrew)
 CALL gopen (outt,rz(buf4),wrtrew)
 
 inblk1(1) = outu
 inblk2(1) = outv
 inblk3(1) = outa
 outblk(1) = outt
 
 j = 1
 DO  i = 1,ncol
   CALL cpystr (inblk1,outblk,0,j)
   j = j + 1
   CALL cpystr (inblk2,outblk,0,j)
   j = j + 1
   CALL cpystr (inblk3,outblk,0,j)
   j = j + 1
 END DO
 
 CALL CLOSE (outu,rew)
 CALL CLOSE (outv,rew)
 CALL CLOSE (outa,rew)
 CALL CLOSE (outt,rew)
 mcb(1) = outt
 mcb(2) = ncol*3
 CALL wrttrl (mcb)
 
!     NORMAL RETURN
 
 600 RETURN
 
!     ERRORS
 
 6000 IF (rc == 6) GO TO 9100
 CALL smsg (rc-2,item,ssnm)
 GO TO 9200
 6100 CALL smsg (7,item,ssnm)
 GO TO 9200
 6200 CALL smsg (rc+4,item,ssnm)
 GO TO 9200
 6313 WRITE  (nout,6314) swm,rss
 6314 FORMAT (a25,' 6313, INSUFFICIENT CORE FOR RCOVR MODULE WHILE ',  &
     'TRYING TO PROCESS', /34X,'PRINTOUT DATA BLOCKS FOR ', 'SUBSTRUCTURE ',2A4)
 GO TO 9200
 9001 n = 1
 GO TO 9100
 9002 n = 2
 GO TO 9100
 9007 n = 7
 9100 CALL mesage (n,FILE,NAME)
 9200 in = 0
 CALL CLOSE (in,rew)
 IF (outt /= 0) CALL CLOSE (outt,rew)
 IF (outu /= 0) CALL CLOSE (outu,rew)
 IF (outv /= 0) CALL CLOSE (outv,rew)
 IF (outa /= 0) CALL CLOSE (outa,rew)
 
 RETURN
END SUBROUTINE rcovva
