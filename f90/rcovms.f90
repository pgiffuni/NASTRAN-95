SUBROUTINE rcovms
     
!     THIS ROUTINE GENERATES THE MODAL SOLUTION ITEM FOR RIGID FORMAT 3
 
 LOGICAL :: mrecvr
 INTEGER :: dry        ,step       ,fss        ,rfno       ,  &
     soln       ,rc         ,swrt       ,srd        ,  &
     eog        ,eoi        ,buf1       ,z          ,  &
     rd         ,rdrew      ,wrt        ,wrtrew     ,  &
     rew        ,NAME(2)    ,FILE       ,sof3
 COMMON /BLANK / dry        ,loop       ,step       ,fss(2)     ,  &
     rfno       ,neigv      ,lui        ,uinms(2,5) ,  &
     nosort     ,uthres     ,pthres     ,qthres
 COMMON /rcovcr/ icore      ,lcore      ,buf1       ,buf2       ,  &
     buf3       ,buf4       ,sof1       ,sof2       , sof3
 COMMON /rcovcm/ mrecvr     ,ua         ,pa         ,qa         ,  &
     iopt       ,rss(2)     ,energy     ,uimpro     ,  &
     range(2)   ,ireq       ,lreq       ,lbasic
 COMMON /names / rd         ,rdrew      ,wrt        ,wrtrew     ,  &
     rew        ,norew      ,eofnrw
 COMMON /zzzzzz/ z(1)
 DATA    lams  , soln /   4HLAMS,4HSOLN   /
 DATA    srd   , swrt,eog,eoi   / 1,2,2,3 /
 DATA    lama  / 102  /,  i7    / 7 /
 DATA    NAME  / 4HRCOV,  4HMS      /
 
 
!     CREATE SOLN FOR RIGID FORMAT 3
 
 IF (mrecvr) GO TO 500
 
!     WRITE GROUP 0
 
 rc = 3
 CALL sfetch (fss,soln,swrt,rc)
 CALL suwrt (fss,2,1)
 CALL suwrt (rfno,1,1)
 CALL suwrt (neigv,1,eog)
 
!     IF NO EIGENVALUES, GO HOME
 
 IF (neigv <= 0) GO TO 430
 
!     COPY RECORD 2 OF LAMA OR CLAMA TO GROUP 1 OF SOLN
 
 FILE = lama
 CALL OPEN (*9001,lama,z(buf1),rdrew)
 CALL fwdrec (*9002,lama)
 CALL fread (lama,itype,1,1)
 nw = 7
 IF (itype == 90) nw = 6
 z(i7) = 0
 i = 1
 410 CALL READ (*9002,*420,lama,z(1),nw,0,nwds)
 CALL suwrt (z,7,i)
 GO TO 410
 420 CALL suwrt (0,0,eog)
 CALL CLOSE (lama,rew)
 
!     FINISH
 
 430 CALL suwrt (0,0,eoi)
 RETURN
 
!     FOR MODAL RECOVER COPY THE LAMS ITEM TO SOLN
 
 500 CALL sfetch (fss,lams,srd,rc)
 IF (rc /= 1) GO TO 6000
 CALL suread (z(1),-2,n,rc)
 IF (n > sof3) GO TO 9008
 CALL sfetch (fss,soln,swrt,rc)
 CALL suwrt (z(1),n,eoi)
 RETURN
 
!     ERROR RETURNS
 
 6000 CALL smsg (rc-2,lams,fss)
 GO TO 9200
 9001 n = 1
 GO TO 9100
 9002 n = 2
 GO TO 9100
 9008 n = 8
 GO TO 9100
 9100 CALL mesage (n,FILE,NAME)
 9200 CALL sofcls
 iopt = -1
 CALL CLOSE (lama,rew)
 RETURN
END SUBROUTINE rcovms
