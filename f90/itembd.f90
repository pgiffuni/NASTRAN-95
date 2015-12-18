BLOCK DATA itembd
!ITEMBD
!     ITEMBD BLOCK DATA
 
!     TO ADD NEW ITEMS TO THE SOF THE FOLLOWING CODE CHANGES MUST BE
!     MADE.
 
!        1) INCREASE THE DIMENSION OF ITEM IN THE ITEMDT COMMON BLOCK.
!        2) INCREASE THE DIMENSION OF THE ITEMXX ARRAY AND ADD ANY
!           ADDITIONAL ARRAYS AND EQUIVALENCES IF NECESSARY.
!        3) INCREASE THE VALUE OF NITEM IN THE DATA STATEMENT.
!        4) ADD THE NEW DATA DESCRIBING THE NEW ITEMS.
!        5) SUBROUTINE EXO2 MUST BE CHANGED IF THE NEW ITEM IS A TABLE.
!           THIS ROUTINE PROCESSES THE SOFOUT(EXTERNAL) STATEMENT.
!        6) SUBROUTINE SOFTOC MUST BE CHANGED IF THE NEW ITEMS WILL
!           INCREASE THE NUMBER IF ITEMS TO MORE THEN 27.
 
!     NOTE... IF THE NUMBER OF ITEMS IS DECREASED THE LENGTH OF THE
!             ITEMDT COMMON BLOCK SHOULD NOT BE DECREASED.  SOFS CREATED
!             ON THE OLDER SYSTEM WILL REQUIRE THE EXTRA SPACE WHEN
!             RESTORING THE ITEM STRUCTURE FOR THAT SOF.
 
 INTEGER :: item01(7,10)  ,item02(7,10)  ,item03(7,5)
 
 COMMON / itemdt /       nitem    ,item(7,25)
 
 EQUIVALENCE   ( item01(1,1) , item(1, 1) ) ,( item02(1,1) , item(1,11) )  &
     ,( item03(1,1) , item(1,21) )
 
!        NITEM     = NUMBER OF ITEMS
!        ITEM(1,I) = ITEM NAME
!        ITEM(2,I) = ITEM TYPE
!                         LE 0 - TABLE ITEM
!                         GE 1 - MATRIX ITEM
!        ITEM(3,I) = EQUIV DATA FOR GROUP 0 OF ITEMS TO BE COPYIED
!                         X 1000 - WORD WITH NUMBER OF NAMES
!                         X 100  - WORD WITH FIRST NAME
!                         X 1    - NUMBER OF WORDS FOR EACH NAME
!        ITEM(4,I) = IMAGE SUBSTRUCTURE DATA
!                         0 - ITEM IS ONLY A POINTER TO PRIMARY DATA
!                         1 - UNIQUE DATA, RETURN TO FREE BLOCK LIST
!        ITEM(5,I) = SECONDARY SUBSTRUCTURE DATA
!                         0 - ITEM IS ONLY A POINTER TO PRIMARY DATA
!                         1 - UNIQUE DATA, RETURN TO FREE BLOCK LIST
!        ITEM(6,I) = HIGHER LEVEL SUBSTRUCTURE DATA
!                         0 - ITEM DOES NOT PERTAIN TO HIGER LEVEL
!                         1 - ITEM DESCRIBES HIGHER LEVEL
!        ITEM(7,1) = EDIT DATA.  EACH BIT IS SET IF THAT ITEM IS IN
!                    THE COORESPONDING EDIT GROUP.  EXAMPLE - A VALUE
!                    OF 36 WOULD CAUSE THE ITEM TO BE DELETED BY
!                    EDIT(32) OR EDIT(4)
 
!***********************************************************************
 
 DATA nitem / 25 /
 
 
!          NAME   TYPE     EQUIV     IMAGE    SECONDARY   HIGHER   EDIT
 
 DATA item01 / 4HEQSS   ,0        ,3005002  ,1        ,1        ,0      ,32  &
     ,4HBGSS   ,0        ,0        ,0        ,0        ,0      ,32  &
     ,4HCSTM   ,0        ,0        ,0        ,0        ,0      ,32  &
     ,4HLODS   ,0        ,4005002  ,1        ,1        ,0      ,36  &
     ,4HPLTS   ,0        ,3004014  ,1        ,1        ,0      ,32  &
     ,4HKMTX   ,1        ,0        ,0        ,0        ,0      ,33  &
     ,4HMMTX   ,1        ,0        ,0        ,0        ,0      ,34  &
     ,4HPVEC   ,1        ,0        ,0        ,0        ,0      ,36  &
     ,4HPOVE   ,1        ,0        ,0        ,1        ,1      ,48  &
     ,4HUPRT   ,1        ,0        ,0        ,1        ,1      ,48 /
 DATA item02 /  &
     4HHORG   ,1        ,0        ,0        ,1        ,1      ,560  &
     ,4HUVEC   ,1        ,0        ,1        ,1        ,0      ,40  &
     ,4HQVEC   ,1        ,0        ,1        ,1        ,0      ,40  &
     ,4HSOLN   ,0        ,0        ,1        ,1        ,0      ,40  &
     ,4HPAPP   ,1        ,0        ,0        ,0        ,0      ,100  &
     ,4HPOAP   ,1        ,0        ,0        ,1        ,1      ,112  &
     ,4HLOAP   ,0        ,4005002  ,1        ,1        ,0      ,100  &
     ,4HLMTX   ,1        ,0        ,0        ,1        ,1      ,48  &
     ,4HGIMS   ,1        ,0        ,0        ,1        ,1      ,48  &
     ,4HPHIS   ,1        ,0        ,0        ,1        ,1      ,288 /
 DATA item03 /  &
     4HLAMS   ,0        ,0        ,0        ,1        ,1      ,288  &
     ,4HK4MX   ,1        ,0        ,0        ,0        ,0      ,160  &
     ,4HBMTX   ,1        ,0        ,0        ,0        ,0      ,160  &
     ,4HPHIL   ,1        ,0        ,0        ,1        ,1      ,288  &
     ,4HHLFT   ,1        ,0        ,0        ,1        ,1      ,560 /

END BLOCK DATA
