BLOCK DATA vdrbd
!VDRBD
! BLOCK DATA FOR THE VECTOR DATA RECOVERY MODULE (VDR).
!*****
 INTEGER :: usetd,casecc,eqdyn ,oeigs ,pp    ,xycdb ,pnl   ,outfle  &
     ,opnl1,scr1  ,scr2  ,buf1  ,buf2  ,buf3  ,cei   ,frq ,trn  ,DIRECT,xset0 ,buf
 
 DIMENSION nam(2)    ,buf(50)      ,masks(6)  &
     ,cei(2)    ,frq( 2)      ,trn(2)       ,modal(2) ,DIRECT(2)
 
 COMMON/vdrcom/vdrcom,idisp ,ivel  ,iacc  ,ispcf ,iloads,istr  &
     ,ielf  ,iadisp,iavel ,iaacc ,ipnl  ,ittl  ,ilsym  &
     ,ifrout,idload,casecc,eqdyn ,usetd ,infile,oeigs  &
     ,pp    ,xycdb ,pnl   ,outfle,opnl1 ,scr1  ,scr2  &
     ,buf1  ,buf2  ,buf3  ,nam   ,buf   ,masks ,cei  &
     ,frq   ,trn   ,DIRECT,xset0 ,vdrreq,modal
 
! DATA DEFINING POSITION OF PARAMETERS IN CASE CONTROL RECORD.
 
 DATA   idisp  / 20/ ,ivel   / 32/ ,iacc   / 29/ ,ispcf  / 35/  &
     ,iloads / 17/ ,istr   / 23/ ,ielf   / 26/ ,iadisp /151/  &
     ,iavel  /154/ ,iaacc  /157/ ,ipnl   / 10/ ,ittl   / 39/  &
     ,ilsym  /200/ ,ifrout /145/ ,idload / 13/
 
! DATA DEFINING GINO FILE NAMES
 
 DATA   casecc /101/ ,eqdyn  /102/ ,usetd  /103/ ,infile /104/  &
     ,oeigs  /105/ ,pp     /105/ ,xycdb  /106/ ,pnl    /107/  &
     ,outfle /201/ ,opnl1  /202/ ,scr1   /301/ ,scr2   /302/
 
! MISC DATA
 
 DATA   buf   /50*0/ ,nam    /4HVDR ,4H    /  &
     ,masks /4,8,16,32,64,128/,xset0/100000000/
 
! DATA DEFINING RIGID FORMATS AND PROBLEM TYPES
 
 DATA   cei    /4HCEIG,4HEN  /,frq   /4HFREQ,4HRESP/  &
     ,trn    /4HTRAN,4HRESP/,modal /4HMODA,4HL   / ,DIRECT /4HDIRE,4HCT  /

END 

