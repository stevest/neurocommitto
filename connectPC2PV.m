function connmat = connectPC2PV(dist)
    % randomly (until newer data) connect PC 2 PV with a chance of 23%.
    % as in:
    % Kawaguchi, 2009
    
    connmat = zeros(size(dist));
    for i=1:size(dist,1) % PC
        for j=1:size(dist,2) % PV
            if(rand(1) <= 0.23)
                connmat(i,j) = 1;
            end
        end
    end
    

end