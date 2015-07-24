function [connmat, gapmat] = connectPV2PV(nPV)
% randomly (until newer data) connect PV 2 PV with a chance of 77%.
% as in: Gibson, Connors, 1999
%Eh, o Gibson einai layers 4 and 6, eno o Galarreta, Hestrin, 1999 (same
%issue) einai smoatosensory/visual L5:
% 66% electrical coupling <=80? appart.
% 18% GABAergic connections
% 4.5% Reciprocal GABAergic connections

connmat = zeros(nPV);
gapmat = zeros(nPV);
start = 2;
for i=1:nPV
    for j=start:nPV
        rn = rand(1);
        if(rn < 0.18)
            if (rand(1) < 0.25)
                connmat(i,j) = 1;
                connmat(j,i) = 1;
                gapmat(i,j) = 1;
            else
                if(rand(1) < 0.5)
                    connmat(i,j) = 1;
                else
                    connmat(j,i) = 1;
                end
                gapmat(i,j) = 1;
            end
        else
            if(rn < 0.66)
                gapmat(i,j) = 1;
            end
        end
    end
    start = start + 1;
end
connmat(logical(eye(nPV))) = rand(1,nPV)<0.85 ;

gabaPercentage = sum(sum((connmat | connmat') .* triu(ones(nPV),1))) / ((numel(connmat)-nPV)/2);
gabaRecipPercentage = sum(sum((connmat & connmat') .* triu(ones(nPV),1))) / ((numel(connmat)-nPV)/2);
gapPercentage = sum(sum(gapmat .* triu(ones(nPV),1))) / ((numel(gapmat)-nPV)/2);

fprintf('Percentage of GABA connected pairs: %f(%%)\n', gabaPercentage*100)
fprintf('Percentage of reciprocally GABA connected pairs: %f(%%)\n', gabaRecipPercentage*100)
fprintf('Percentage of gap connected pairs: %f(%%)\n', gapPercentage*100)



end