function [code_m] =codematrixgen(code,maxlags,offset)

code_m=zeros(size(code))';

% Genero matrice codici shiftati di 1*Ts x volta (ATTENTO AL SEGNO DELL'INDICE ii)
for ii = -(maxlags+offset/2):maxlags-offset/2
    code_m(:,ii+maxlags+1+offset/2)= circshift(code',ii);
    
end

end
