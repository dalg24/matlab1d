function  [CHLSAT]=hlsat( PRE, TSAT )
% C=======================================================================
% C
% C     Calcul de l'enthalpie liquide de saturation
% C
% C=======================================================================
% *>PRE     Pression en Pa (1bar=1.0E5 Pa)
% *>TSAT    Temperature de saturation issue de la fonction CTSAT
% C=======================================================================

C(1)= 0.18637E-02;
C(2)=-0.50352E+05;
C(3)= 0.36540E-12;
C(4)=-0.30413E-05;
C(5)= 0.40047E+04;
C(6)=-0.95261E-08;
C(7)=-0.25785E+00;
C(8)= 0.20641E+08;
C(9)= 399.980E+00;

P(1)= C(1)*PRE+C(2);
P(2)=(C(3)*PRE+C(4))*PRE +C(5);
P(3)=(C(6)*PRE+C(7))*PRE +C(8);
P(4)=1./(C(9)-TSAT);
% C
% C=======================================================================
% C

CHLSAT=P(1)+P(2)*TSAT+P(3)*P(4);