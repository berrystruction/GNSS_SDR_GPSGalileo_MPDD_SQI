function [lat, lon, h] = xyz2llh(x, y, z)
%XYZ2LLH		[lat, lon, h]=xyz2llh(x, y, z)
%
%Calculates latitude, longitude, and height relative to the
%WGS-84 reference ellipsoid, given ECEF cartesian coordinates.
%Inputs can be vectors.
%
% INPUT:
%       x:      vector of X positions
%       y:      vector of Y positions
%       z:      vector of Z positions
%
% OUTPUT:
%       lat:    vector of latitude
%       lon:    vector of longitude
%       h:      vector of heights
%
% NOTE: lat and lon are returned in radians


%Set constants
	%Note: either define a and b
 	
		a=6378137.0;
		b=6356752.3142;
		f=(a-b)/a;

	%or a and f:

		%a=6378137.0;
		%f=1/298.257223563;		
       	%b=a-f*a;
				
	e2=2*f - f^2;
	E2=(a^2-b^2)/(b^2);
	p=sqrt(x.^2 + y.^2);
   
%Calculate longitude

	lon = atan2(y,x);

%Calculate latitude
   
    theta=atan((z*a)./(p*b));
	lat=atan((z+E2*b*sin(theta).^3)./(p-e2*a*cos(theta).^3) );
   
%Calculate height
   
	N=a./sqrt(1-e2*sin(lat).^2);
	h=p./cos(lat)-N;