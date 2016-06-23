function [LAF]=LAF_SQM(u,d,M,version)

LOS=1;
%covarianceMethod=true;
%v2=true;
%version=3;
%version=2;

switch version
    %if version==2
    case 2
        [N,U,w,y,ynew,D,e,powers,error_vec,Particular_ratio]=LAFmethod_v2(u,d,M,LOS);
        %else
        % volevo fare lo shift come nel DLL col metodo early late
        %[N,U,w,y,ynew,D,e,powers,error_vec,Particular_ratio]=LAFmethod_v3(u',d',M,LOS,1);
    case 3
        % implemento il metodo delle autocorrelazioni al posto delle
        % covarianze con CLS
        [N,U,w,y,ynew,D,e,powers,error_vec,Particular_ratio]=LAFmethod_autocorrelation_CLS(u',d',M,LOS,1);
    case 4
        % implemento il metodo delle autocorrelazioni al posto delle
        % covarianze con LS
        [N,U,w,y,ynew,D,e,powers,error_vec,Particular_ratio]=LAFmethod_autocorrelation_LS(u',d',M,LOS,1);
    case 5
        [N,U,w1,y1,ynew,D1,e,powers,error_vec,Particular_ratio]=LAFmethod_autocorrelation_IHT(u',d',M,LOS,1,1);
        errLAFref=sum((D1-y1).^2);
        w=w1; y=y1; D=D1; LAF.k=1; LAF.errLAF=errLAFref;
    case 6
        % implemento il metodo delle autocorrelazioni al posto delle
        % covarianze con LS con equality constraint
        [N,U,w,y,ynew,D,e,powers,error_vec,Particular_ratio]=LAFmethod_autocorrelation_LS_EC(u',d',M,LOS,1);
   
        % nuovo metodo LAF allungato
        %[N,U,w,y,ynew,D,e,powers,error_vec,Particular_ratio]=LAF_New_method(u,d,M,LOS,tott);
        %end
end

LAF.u=u;
LAF.d=d;
LAF.N=N;
LAF.U=U;
LAF.w=w;
LAF.y=y;
LAF.ynew=ynew;
LAF.D=D;
LAF.e=e;
%LAF.powers=powers;
%LAF.error_vec=error_vec;
LAF.Particular_ratio=Particular_ratio;

end