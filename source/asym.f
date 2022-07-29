C     *****************************************************************         

      FUNCTION ASYM(I)
      CHARACTER*2 ASYM
      CHARACTER*2 ELEMNT(86), AS
      DATA ELEMNT/'H ','He',
     1 'Li','Be','B ','C ','N ','O ','F ','Ne',
     2 'Na','Mg','Al','Si','P ','S ','Cl','Ar',
     3 'K ','Ca','Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu',
     4 'Zn','Ga','Ge','As','Se','Br','Kr',
     5 'Rb','Sr','Y ','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag',
     6 'Cd','In','Sn','Sb','Te','I ','Xe',
     7 'Cs','Ba','La','ce','pr','nd','pm','sm','eu','gd','tb','dy',
     8 'ho','er','tm','yb','lu','Hf','Ta','W ','Re','Os','Ir','Pt',
     9 'Au','Hg','Tl','Pb','Bi','Po','At','Rn'/
      ASYM=ELEMNT(I)
      RETURN
      END

      SUBROUTINE UPPER(AS)
      CHARACTER*2 AS
      NSP=ICHAR(' ')
      ND=ICHAR('A')-ICHAR('a')
      DO 10 I=1,2
         J=ICHAR(AS(I:I))
         IF(J.NE.NSP)J=J+ND
         AS(I:I)=CHAR(J)
  10  CONTINUE
      RETURN
      END

      SUBROUTINE UPPER10(AS)
      CHARACTER*(*) AS
      NSP=ICHAR(' ')
      ND=ICHAR('A')-ICHAR('a')
      DO 10 I=1,10
         J=ICHAR(AS(I:I))
         IF(J.NE.NSP)J=J+ND
         AS(I:I)=CHAR(J)
  10  CONTINUE
      RETURN
      END
