SUBROUTINE reig
     
!     READ   KAA,MAA,MR,DM,EED,USET,CASECC/LAMA,PHIA,MI,OEIGS/C,N,IPROB
!            /V,N,NUMMOD/C,N,ICASE/C,N,XLAMDA $
 
 INTEGER :: sysbuf   ,eigr(4)  ,icore(12),casecc   ,FILE     ,  &
     error(3) ,feerx    ,sdet     ,udet     ,inv      ,  &
     sinv     ,uinv     ,givi     ,sturm    ,givn(7)  ,  &
     dm       ,eed      ,lama     ,phia     ,mi       ,  &
     uset     ,pout     ,scr1     ,scr2     ,scr3     ,  &
     scr4     ,scr5     ,scr6     ,scr7     ,option   , subnam(2),ix(7)    ,optn2
 REAL :: lmin     ,lmax     ,lfreq    ,mb(1)
 DOUBLE PRECISION :: dcore(1)
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm      ,uwm      ,uim
 COMMON /BLANK / iprob(2) ,nummod   ,icase    ,xlamda  &
     /invpwx/ ifilk(7) ,ifilm(7) ,ifillm(7),ifilvc(7),iscr1    ,  &
     iscr2    ,iscr3    ,iscr4    ,iscr5    ,iscr6    ,  &
     iscr7    ,iscr8    ,idump    ,lmin     ,lmax     ,  &
     noest    ,ndplus   ,ndmnus   ,eps      ,novect  &
     /givn  / dum(100) ,n        ,lfreq    ,order    ,d1       ,  &
     hfreq    ,d2       ,nv       ,nprt     ,d4       , nfr
 COMMON /feercx/ ifkaa(7) ,ifmaa(7) ,iflelm(7),iflvec(7),sr1fle   ,  &
     sr2fle   ,sr3fle   ,sr4fle   ,sr5fle   ,sr6fle   ,  &
     sr7fle   ,sr8fle   ,dmpfle   ,nord     ,xlmbda   ,  &
     neig     ,mord     ,ibk      ,critf    ,northo   ,  &
     iflrva   ,iflrvc   ,iepx /ntime / lntime   ,tcons(15)  &
     /regean/ im(7)    ,ik(7)    ,iev(7)   ,scr1     ,scr2     ,  &
     scr3     ,scr4     ,scr5     ,lcore    ,rmax     ,  &
     rmin     ,mz       ,nev      ,epsi     ,rminr    ,  &
     NE       ,nit      ,nevm     ,scr6     ,scr7     ,  &
     nfound   ,lamda    ,ibuck    ,nsym
 COMMON /sturmx/ sturm    ,shftpt   ,KEEP     ,ptshft   ,nr       , shftzo  &
     /reigkr/ option   ,optn2  &
     /system/ sysbuf   ,nout     ,nogo     ,ksys(51) ,jprec
 COMMON /condas/ pi       ,twopi    ,radeg    ,degra    , fps
 COMMON /packx / itp1     ,itp2     ,iip      ,jjp      ,incrp
 COMMON /unpakx/ itu      ,iiu      ,jju      ,incru /zzzzzz/ core(1)
 EQUIVALENCE     (givn(1) ,core(1))
 EQUIVALENCE     (tcons(4),apc    ) ,(tcons(5),apu      ) ,  &
     (tcons(8),mb(1)  ) ,(error(2),subnam(1)) , (dcore(1) ,core(1) ,icore (1))
 DATA            eigr                                   ,casecc/  &
     307      ,3        ,107      ,1        ,107   /,  &
     sdet     ,udet     ,inv      ,sinv     ,i0    /  &
     4HSDET   ,4HUDET   ,4HINV    ,4HSINV   ,0     /,  &
     uinv     ,givi     ,kaa      ,maa      ,mr    /  &
     4HUINV   ,4HGIV    ,101      ,102      ,103   /,  &
     dm       ,eed      ,uset     ,lama     ,phia  /  &
     104      ,105      ,106      ,201      ,202   /,  &
     mi       ,pout     ,icr1     ,icr2     ,mode  /  &
     203      ,204      ,301      ,302      ,4HMODE/,  &
     error                        ,feerx    ,mgiv  /  &
     4HEED    ,4HREIG   ,4H       ,4HFEER   ,4HMGIV/
 
 
 ibuck  = 1
 lcore  = korsz(core) - sysbuf - 3
 llcore = lcore - sysbuf
 CALL gopen (lamda,core(lcore+1),1)
 CALL CLOSE (lamda,2)
 IF (iprob(1) /= mode) ibuck = 3
 sturm  = -1
 KEEP   = 0
 shftpt = 0.0
 ptshft = 0.0
 nr     = 0
 shftzo = 0.0
 CALL OPEN (*10,casecc,core(lcore+1),0)
 CALL skprec (casecc,icase)
 CALL fread (casecc,icore,166,1)
 CALL CLOSE (casecc,1)
 method = icore(5)
 GO TO 20
 10 method = -1
 20 FILE   = eed
 CALL preloc (*170,core(lcore+1),eed)
 CALL locate (*40,core(lcore+1),eigr(ibuck),iflag)
 30 CALL READ (*40,*40,eed,core(1),18,0,iflag)
 IF (method == icore(1) .OR. method == -1) GO TO 50
 GO TO 30
 
!     NO SET NUMBER FOUND
 
 40 CALL mesage (-32,method,error)
 
!     FOUND DATA CARD
 
 50 norm = icore(10)
 CALL CLOSE (eed,1)
 
!     TEST THE SIZE OF THE K AND M MATRICES VIA THEIR TRAILERS
 
 CALL rdtrl (ik(1))
 CALL rdtrl (im(1))
 IF (im(2) == ik(2) .AND. im(3) == ik(3)) GO TO 51
 
!     K AND M MATRICES ARE NOT OF THE SAME SIZE
 
 WRITE  (nout,200) ufm
 200 FORMAT (a23,' 3131, INPUT STIFFNESS AND MASS MATRICES ARE NOT ',  &
     'COMPATIBLE.')
 CALL mesage (-37,0,error(2))
 
!     K AND M MATRICES ARE COMPATIBLE
 
 51 CONTINUE
 
!     CHECK TO SEE IF THE INPUT STIFFNESS AND/OR MASS MATRIX IS NULL
 
 IF (ik(6) == 0 .OR. im(6) == 0) CALL mesage (-60,0,0)
 
!     SET FLAG FOR THE METHOD OF ANALYSIS AND THE PROPER
!     TYPE OF DECOMPOSITION
 
 option = icore(2)
 optn2  = icore(3)
 IF (option == givi .OR. option == udet .OR. option == uinv) GO TO 53
 IF (option == feerx .OR. option == mgiv) GO TO 53
 IF (option == sdet  .OR. option == sinv) GO TO 52
 option = udet
 IF (icore(2) == inv) option = uinv
 IF (im(4) /= 6 .OR. ik(4) /= 6) GO TO 53
 option = sdet
 IF (icore(2) == inv) option = sinv
 GO TO 53
 52 IF (im(4) == 6 .AND. ik(4) == 6) GO TO 53
 WRITE (nout,2100) uwm
 option = udet
 IF (icore(2) == sinv) option = uinv
 53 isil  = icore(12)
 i     = 9
 epsii = core(i)
 IF (ibuck == 3) GO TO 60
 
!     CONVERT FREQUENCY TO LAMDA
 
 IF ((icore(2) == givi .OR. icore(2) == mgiv) .AND. icore(7) > 0) GO TO 55
 IF (core(i0+4) >= 0.0) GO TO 55
 WRITE (nout,2000) uwm
 core(i0+4) = 0.0
 55 core(i0+4) = fps*core(i0+4)*core(i0+4)
 IF (icore(2) /= feerx) core(i0+5) = fps*core(i0+5)*core(i0+5)
 60 CONTINUE
 core4  = core(i0+4)
 core5  = core(i0+5)
 icore6 = icore(6)
 icore7 = icore(7)
 icore8 = icore(8)
 IF (icore(2) == givi .OR. icore(2) == mgiv) GO TO 70
 IF (icore(2) == feerx) GO TO 65
 IF (icore(7) == 0) icore(7) = 3*icore(6)
 icore7 = icore(7)
 
!     FEER, INVERSE POWER AND DETERMINANT METHODS
 
!     CHECK IF IT IS A NORMAL MODES PROBLEM OR A BUCKLING PROBLEM
 
 65 IF (ibuck == 3) GO TO 80
 
!     NORMAL MODES PROBLEM
 
!     CHECK FOR APPEND
 
 IF (nummod <= 0) GO TO 70
 ix(1) = phia
 CALL rdtrl (ix)
 IF (ix(1) <= 0 .OR. ix(2) <= 0) GO TO 70
 
!     NEW EIGENVALUES AND EIGENVECTORS WILL BE APPENDED TO THOSE
!     PREVIOUSLY CHECKPOINTED
 
 nr = ix(2)
 IF (nummod < nr) nr = nummod
 WRITE (nout,2200) uim,nr
 
!     RETRIEVE EIGENVALUES AND EIGENVECTORS PREVIOUSLY CHECKPOINTED.
 
!     COPY OLD EIGENVALUES FROM LAMA FILE TO ICR1 FILE.
 
!     COPY OLD EIGENVECTORS FROM PHIA FILE TO ICR2 FILE.
 
 CALL read7 (nr,lama,phia,icr1,icr2)
 GO TO 80
 
!     NO APPEND
 
!     CHECK IF RIGID BODY MODES ARE TO BE COMPUTED SEPARATELY
 
 70 ix(1) = mr
 CALL rdtrl (ix)
 IF (ix(1) < 0) GO TO 75
 
!     COMPUTE RIGID BODY MODES
 
 CALL read1 (dm,mr,scr4,scr5,scr3,icr2,uset,nr,icr1,scr6)
 
!     RIGID BODY EIGENVALUES ARE ON ICR1
 
!     RIGID BODY EIGENVECTORS ARE ON ICR2
 
 75 IF (option == givi .OR. option == mgiv) GO TO 100
 80 IF (option == feerx) GO TO  95
 IF (option == sdet) GO TO 109
 IF (option == udet) GO TO 110
 
 
!     INVERSE POWER METHOD
!     ********************
 
 lmin   = core4
 lmax   = core5
 noest  = icore6
 ndplus = icore7
 ndmnus = 0
 IF (ibuck == 3) ndmnus = icore8
 eps = epsii
 IF (eps <=      0.) eps = .0001
 IF (eps < .000001) eps = .000001
 CALL rdtrl (ifilk(1))
 CALL rdtrl (ifilm(1))
 novect = nr
 CALL invpwr
 method = 2
 nummod = novect
 GO TO 140
 
 
!     FEER METHOD
!     ***********
 
 95 iflrva = icr1
 iflrvc = icr2
 xlmbda = core4
 neig   = icore7
 iepx   = icore8
 IF (ibuck == 3) neig = icore6
 northo = nr
 critf  = core5
 ix(1)  = kaa
 CALL rdtrl (ix)
 n = ix(2)
 IF (critf == 0.) critf = .001/n
 CALL feer
 method = 2
 nummod = mord + nr
 CALL sswtch (26,l26)
 IF (nummod > neig .AND. l26 /= 0) nummod = neig
 ifilk(2) = nord
 GO TO 140
 
 
!     GIVENS METHOD
!     *************
 
 100 lfreq = core4
 hfreq = core5
 method= 3
 nfr   = nr
 nprt  = icore6
 nv    = icore7
 givn(   1) = kaa
 givn(i0+2) = maa
 givn(i0+3) = phia
 DO  i = 1,4
   givn(i+3) = eigr(i)
 END DO
 CALL givens
 nnv = givn(1)
 nummod = n
 GO TO 145
 
 
!     DETERMINANT METHOD
!     ******************
 
 109 nsym  = 1
 110 method= 4
 rmin  = core4
 rmax  = core5
 IF (rmin == 0.0 ) rmin = rmax*1.0E-4
 rminr = -.01*rmin
 nev   = icore6
 IF (ibuck == 3 .AND. epsii /= 0.0) epsi = epsii
 nevm  = icore7
 CALL rdtrl (im(1))
 iev(3) = ik(3)
 IF (nevm > ik(3)) nevm = ik(3)
 mz = nr
 
!     PICK UP UNREMOVED FREE BODY MODES
 
 IF (icore8 > nr) mz = -icore8
 iev(2) = nr
 CALL detm
 nummod = nfound + nr
 ifilk(2) = iev(3)
 
!     SORT EIGENVECTORS AND VALUES
 
 140 IF (nummod == 0) GO TO 160
 CALL read3 (nummod,ifilk(2),lamda,iev,phia,lama)
 145 IF (method == 2 .OR. nummod == 1) GO TO 150
 
!     CHECK ORTHOGONALITY
 
 ifilvc(1) = phia
 CALL rdtrl (ifilvc(1))
 CALL read4 (lama,ifilvc(1),scr1,epsii,maa)
 150 CONTINUE
 
!     SET FLAG FOR GIVENS METHOD FOR USE IN READ2 ROUTINE
 
 dum(1) = 0.0
 IF (method == 3) dum(1) = 1.0
 nv = nnv
 
!     FORM MODAL MASS, NORMALIZE AND FORM SUMMARY FILE.
 
 CALL read2 (maa,phia,scr1,norm,isil,xxx,mi,lama,pout,icr2, epsii,scr6)
 GO TO 165
 160 nummod = -1
 CALL read5 (pout)
 165 IF (nogo == 14) WRITE (nout,166)
 166 FORMAT ('0*** THIS NASTRAN JOB WILL BE TERMINATED')
 RETURN
 
 170 ip1 = -1
 180 CALL mesage (ip1,FILE,subnam)
 GO TO 180
 
!     ERROR MESSAGES
 
 2000 FORMAT (a25,' 2367, FREQUENCY F1 (FIELD 4) ON THE EIGR BULK DATA',  &
     ' CARD IS NEGATIVE', /5X,  &
     'IT IS ASSUMED TO BE ZERO FOR CALCULATION PURPOSES.',/)
 2100 FORMAT (a25,' 2368, SYMMETRIC DECOMPOSITION IS SPECIFIED ON THE ',  &
     'EIGR BULK DATA CARD, BUT', /5X,  &
     'UNSYMMETRIC DECOMPOSITION WILL BE USED AS THIS IS THE ',  &
     'PROPER TYPE OF DECOMPOSITION FOR THIS PROBLEM.')
 2200 FORMAT (a29,' 3143, THE EIGENVALUES AND EIGENVECTORS FOUND IN ',  &
     'THIS ANALYSIS WILL BE APPENDED', /5X,'TO THE',i8,  &
     ' EIGENVALUES AND EIGENVECTORS COMPUTED EARLIER.')
 
END SUBROUTINE reig
