SUBROUTINE hccom(itype,lcore,icore,nextz,kcount)
     
! COMBINES HC CENTROID INFO ON SCR6 TO HCCENS
 
 
 INTEGER, INTENT(IN OUT)                  :: itype
 INTEGER, INTENT(IN OUT)                  :: lcore
 INTEGER, INTENT(IN OUT)                  :: icore
 INTEGER, INTENT(IN)                      :: nextz
 INTEGER, INTENT(IN)                      :: kcount
 INTEGER :: scr6,hccens,iz(1),id(2),nam(2),mcbh(7)
 LOGICAL :: incore,bldp,eor
 DIMENSION hc(63)
 COMMON/system/idum,iout
 COMMON/zzzzzz/z(1)
 COMMON/packx/ita,itb,ii,jj,incr
 COMMON/zblpkx/a(4),irow
 EQUIVALENCE (z(1),iz(1))
 DATA mcbh/307,0,0,2,1,0,0/
 DATA scr6,hccens/306,307/
 DATA nam/4HHCCO,4HM   /
 
 ita=1
 itb=1
 ii=1
 incr=1
 icount=0
 nskip=0
 nskip1=0
 in=0
 eor=.false.
 bldp=.false.
 incore=.true.
 
! IF TYPE IS 24 JUST PACK ZEROS ON HCCENS
 
 IF(itype == 24)GO TO 60
 
 CALL gopen(scr6,z(lcore+1),0)
 
! SCR6 HAS 3 ENTRIES PER ELEMENT-ID,NUMBER OF POINTS AT WHICH HC IS
! COMPUTED=N, AND 3*N HC VALUES--THERE IS ONE RECORD PER CARD TYPE
! ON SCR6 FOR THIS SUBCASE
! IF  .NOT. INCORE, THEN WE ARE BACK HERE DUE TO SPILL LOGIC AND  ARE
! TRYING TO FINISH THE FIRST RECORD. SO WE MUST SKIP THE PART OF THE
! RECORD PREVIOUSLY READ.
 
 5 IF(.NOT.incore)CALL fread(scr6,id,-nskip,0)
 inword=0
 10 CALL READ (*1002,*20,scr6,id,2,0,nwds)
 nskip=nskip+2
 inext=nextz+inword
 nwords=3*id(2)
 IF(inext+nwords > icore)GO TO 80
 CALL fread(scr6,z(inext),nwords,0)
 
! INWORD IS THE NUMBER OF WORDS READ INTO CORE ON THIS READ
! NSKIP IS THE TOTAL NUMBER OF WORDS READ FROM SCR6 FROM THIS RECORD
! ICOUNT IS THE TOTAL NUMBER OF WORDS SAVED IN CORE FROM THIS RECORD
 
 inword=inword+nwords
 nskip=nskip+nwords
 icount=icount+nwords
 GO TO 10
 
 20 eor=.true.
 IF(.NOT.incore)GO TO 95
 
! CHECK ON COUNT CONSISTENCY
 
 IF(icount /= kcount)GO TO 500
 
! EOR ON SCR6, I.E. END OF HC FOR A GIVEN CARD TYPE IN THIS SUBCASE.
! IF OTHER CARD TYPES EXIST IN THIS SUBCASE, THEY ARE IN SUBSEQUENT
! RECORDS. ADD RESULTS TO PREVIOUS ONES
 
 30 jcount=0
 35 CALL READ (*50,*30,scr6,id,2,0,nwds)
 nwords=3*id(2)
 CALL fread(scr6,hc,nwords,0)
 
! ADD TO PREVIOUS HC FOR THIS ELEMENT- ALL ELEMENTS SHOULD BE  ON SCR6
! IN SAME ORDER IN EVERY RECORD
 
 nj=nextz+jcount-1
 DO  i=1,nwords
   z(nj+i)=z(nj+i)+hc(i)
 END DO
 jcount=jcount+nwords
 IF((.NOT.incore).AND.jcount == inword)GO TO 90
 GO TO 35
 
! INFO WILL NOT FIT IN CORE - SPILL LOGIC
 
 80 incore=.false.
 90 CALL fwdrec (*1002,scr6)
 
! SKIP APPROPRIATE NUMBER OF WORDS IN THIS RECORD TO ACCOUNT FOR
! THE PORTION OF THIS RECORD PREVIOUSLY READ
 
 95 CALL READ (*50,*1003,scr6,id,-nskip1,0,nwds)
 GO TO 30
 
 
! DONE FOR THIS SUBCASE. PACK RESULTS. CLOSE SCR6 AND REOPEN TO WRITE
! NEXT SUBCASE (IF ALL DATA CAN FIT INTO CORE)
 
 50 IF(incore)GO TO 57
 
! SPILL LOGIC-PACK OUT INWORD WORDS. THEN REWIND  SCRL AND SKIP DOWN
! AS NECESSARY
 
 IF(.NOT.bldp)CALL bldpk(1,1,scr6,0,0)
 bldp=.true.
 DO  k=1,inword
   a(1)=z(nextz+k-1)
   irow=in+k
   CALL zblpki
 END DO
 IF(eor)GO TO 58
 
 in=in+inword
 CALL REWIND(scr6)
 CALL fwdrec (*1002,scr6)
 nskip=nskip-2
 nskip1=nskip
 GO TO 5
 
 57 CALL CLOSE(scr6,1)
 jj=icount
 mcbh(3)=jj
 CALL pack(z(nextz),hccens,mcbh)
 GO TO 70
 
! DONE FOR THIS SUBCASE (SPILL LOGIC)
 
 58 CALL CLOSE(scr6,1)
 mcbh(3)=icount
 CALL bldpkn(scr6,0,mcbh)
 GO TO 70
 
 
! PACK A COLUMN OF ZEROS CORRESPONDING TO REMFLUX
 
 60 mcbh(3)=kcount
 CALL bldpk(1,1,hccens,0,0)
 CALL bldpkn(hccens,0,mcbh)
 
 70 CALL wrttrl(mcbh)
 IF(itype == 24)GO TO 75
 
! CHECK ON COUNT CONSISTENCY
 
 IF(incore)GO TO 75
 IF(icount /= kcount)GO TO 500
 
 75 CALL gopen(scr6,z(lcore+1),1)
 RETURN
 
 500 WRITE(iout,501)
 501 FORMAT(58H0***system fatal error,logic error,counts are off in hccom)
 CALL mesage(-61,0,0)
 
 1002 CALL mesage(-2,scr6,nam)
 1003 CALL mesage(-3,scr6,nam)
 RETURN
END SUBROUTINE hccom
