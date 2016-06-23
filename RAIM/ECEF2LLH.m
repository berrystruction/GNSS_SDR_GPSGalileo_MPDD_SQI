function [LLH]=ECEF2LLH(ECEF)
	
global Const


a=6378137;
b=6356752.3142;
e=sqrt(1-(power(b,2)/power(a,2)));
e1=e*(a/b);

num=ECEF(3)*a;
den=b*(sqrt(power(ECEF(1),2)+power(ECEF(2),2)));

if(den==0)
    return
end


Teta=atan2(num,den);

%Compute of Latitude

num=ECEF(3)+(power(e1,2)*b*power(sin(Teta),3));
den=(sqrt(power(ECEF(1),2)+power(ECEF(2),2))) -(power(e,2)*a*power(cos(Teta),3));
LLH(1) = atan2(num,den)*Const.RADIANTI_GRADI;

%Compute of longitude.

LLH(2) = atan2(ECEF(2),ECEF(1))*Const.RADIANTI_GRADI;


%Compute of altitude.

LLH(3)=(sqrt(power(ECEF(1),2)+power(ECEF(2),2))/ cos(LLH(1)*Const.GRADI_RADIANTI))-(a/(sqrt(1-(power(e,2)*power(sin(LLH(1)*Const.GRADI_RADIANTI),2)))));

