function [Azimuth,Elevation]=ComputeAzEl(SvPos,Solution)

global Const

P = 0; % user distance from earth centre projected on the equatorial plane
R = 0; % user distance from earth centre
S = 0;

E = zeros(3,3);
D = zeros(1,3);

P = norm(Solution(1:2));
R = norm(Solution(1:3));

if(R==0 || P==0)
    return
end

E(1,1) = -1.0*(Solution(2)/P);
E(1,2) = +1.0*(Solution(1)/P);
E(1,3) = 0.0;
E(2,1) = -1.0*(Solution(1)*Solution(3)/(P*R));
E(2,2) = -1.0*(Solution(2)*Solution(3)/(P*R));
E(2,3) = P/R;
E(3,1) = Solution(1)/R;
E(3,2) = Solution(2)/R;
E(3,3) = Solution(3)/R;

for(k=1:3)
    for(t=1:3)
        D(k) = D(k) + (SvPos(t)-Solution(t)) * E(k,t);
    end
end

S = D(3)/norm(D);

if(S==1)
    Elevation = 0.5*Const.PI;
else
    Elevation = atan(S/sqrt(1-(S*S)));
end

if(D(2)==0 && D(1)>0)

    Azimuth = 0.5*Const.PI;

else

    if(D(2)==0 && D(1)<0)

        Azimuth = 1.5*Const.PI;

    else

        Azimuth = atan(D(1)/D(2));
        if(D(2)<0.0)
            Azimuth = Azimuth + Const.PI;
        else
            if(D(2)>0 && D(1)<0)
                Azimuth = Azimuth + (2*Const.PI);
            end
        end
    end
end
