function corrRho=Ionospheric_Klobuchar_correction(xyz,satPositions,alpha_coeff,beta_coeff,transmitTime,c)
Nsat=length(satPositions);
corrRho=zeros(Nsat,1);
for SVindex= 1:Nsat
    [Delta_I]=Error_Ionospheric_Klobuchar(xyz,satPositions(SVindex,:),alpha_coeff,beta_coeff,transmitTime);
    corrRho(SVindex)=Delta_I*c;
end

%xyzdt(1:3),satPositions(:,activeChnList),eph(trackResults(activeChnList(1)).PRN).alpha_coeff,eph(trackResults(activeChnList(1)).PRN).beta_coeff,


end