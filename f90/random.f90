SUBROUTINE random
     
!     RANDOM ANALYSIS MODULE
 
!     INPUTS   CASECC,XYCB,DIT,DISP,SPCF,LOAD,STRESS,FORCE,PSDL  (9)
 
!     OUTPUTS  PSDF,AUTO  (2)
 
!     SCRATCHES (0)
 
!     PARAMETERS 1 INTEGER
 INTEGER :: casecc,xycb,dit,ifile(5),psdl,psdf,auto
 COMMON /BLANK/ icoup
 DATA xycb,dit,psdl,ifile,casecc/101,102,103,104,105,106,107,108, 109/
 DATA psdf,auto /201,202/
 DATA nfile              /5/
 
!     INITIALIZE + SET UP
 
 CALL rand7(ifile,nfile,psdl,dit,icoup,nfreq,npsdl,ntau,ltab, casecc,xycb)
 IF( icoup < 0) THEN
   GO TO    10
 ELSE IF ( icoup == 0) THEN
   GO TO    20
 ELSE
   GO TO    30
 END IF
 10 RETURN
 
!     UNCOUPLED
 
 20 CALL rand5(nfreq,npsdl,ntau,xycb,ltab,ifile,psdf,auto,nfile)
 GO TO 10
 
!     COUPLED
 
 30 CALL rand8(nfreq,npsdl,ntau,xycb,ltab,ifile,psdf,auto,nfile)
 GO TO 10
END SUBROUTINE random
