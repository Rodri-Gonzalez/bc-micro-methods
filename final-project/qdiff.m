% Solves seller differential equation to find optimal quality, given mg
% cost draw from D.                  

function dq = qdiff(c,q,cdfD,pdfD,cdfH,pdfH,F_0,gamm,alph,delt,p,R_ub)

dq      = - (2*gamm*delt*(p./c - q) .* cdfH(q)*alph*(1-F_0).*pdfD(c))./...
              ((delt + alph.*(1-F_0).*cdfD(c))^3 + gamm.*delt.*cdfH(q).*(delt + alph.*(1-F_0).*cdfD(c)) - ...
               gamm.*delt.*(p./c - q).*pdfH(q).*(delt + alph.*(1-F_0).*cdfD(c)));

           
dq = dq(:);