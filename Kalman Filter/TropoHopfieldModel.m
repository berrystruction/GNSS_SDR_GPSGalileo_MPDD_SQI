function [TropoDelay,Sigma2Tropo]=TropoHopfieldModel(Elevation,H)

global Const

Zhdy     				= 0;

Dhdy     				= 0;
Dwet     				= 0;
M_El     				= 0;



Num     = 1e-6*Const.K1*Const.Rd*Const.P0_45;
Den     = Const.Gm;
Zhdy    = Num/Den;

Num     = 1e-6*Const.K2*Const.Rd;
Den     = (Const.Gm*(Const.Lambda0_45+1)) - (Const.Beta0_45*Const.Rd);
Zwet    = (Num/Den)*(Const.e0_45/Const.T0_45);

Num     = 1-((Const.Beta0_45*H)/Const.T0_45);
Exp1    = Const.G/(Const.Rd*Const.Beta0_45);
Exp2    = (((Const.Lambda0_45+1)*Const.G)/(Const.Rd*Const.Beta0_45))-1;

if(Num<=0)
    fprintf('\nNum<=0\n')
    return
end

Dhdy = Zhdy*power(Num,Exp1);
Dwet = Zwet*power(Num,Exp2);

M_El = 1.001/sqrt( (0.002001 + power(sin(Elevation),2)));

Sigma2Tropo = power((M_El*0.12),2);

TropoDelay = (-1.0*(Dhdy+Dwet)*M_El);
