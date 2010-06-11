function [CTLIQ]=tliq( HL, PRE )
% C=======================================================================
% C
% C     Calcul de la temperature liquide en fonction de P et Hl
% C
% 
% C=======================================================================
% *>HL      Enthalpie liquide 
% *>PRE     Pression en Pa
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
P(2)=(C(3)*PRE+C(4))*PRE+C(5);
P(3)=(C(6)*PRE+C(7))*PRE+C(8);

% C=======================================================================
TEMPO=(P(2)*C(9)-P(1)+HL)^2 + 4*P(2)*((P(1)-HL)*C(9)+P(3));

CTLIQ=0.5E+0/P(2)*(P(2)*C(9)-P(1)+HL-sqrt(TEMPO));

