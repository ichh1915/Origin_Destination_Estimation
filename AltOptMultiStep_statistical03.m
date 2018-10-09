function [Xhat,Phat,E_Y] = AltOptMultiStep_statistical03( Y0,Xpri,Pinit,PC,itN,e_tol )
% Note the dimension of Y0
% Y0 is the observed data with convolution head and tail [24 by 63]
[n_l,n_n,tau_max] = size(Pinit); % [24 by 9 by 4]
[Xm,tau_max_x] = size(Xpri);    % [9 by 60]
if Xm~=n_n
    error('Xm should = n_n');
end

Xhat = Xpri;
Phat = Pinit;

tic;
E_Y = zeros(itN,1);
for itn = 1:itN
    %% fix X^t and optimise P^t
    % Xcircle needs to be cut the head and tail
    X0 = reshape( Xhat,1,n_n,tau_max_x);
    X0 = MToep( X0,tau_max );
    X0 = X0( tau_max:tau_max_x,: );     
    P0 = reshape(Phat,n_l*n_n*tau_max,1);
    
    cvx_begin 
        cvx_precision low
        variables Pcircle(n_l,n_n*tau_max) %Psym =reshape(Pcircle,n_l*n_n*tau_max,1)*reshape(Pcircle,1,n_l*n_n*tau_max)        
        maximize sum(sum((Y0').*log(X0*(Pcircle')) - (X0*(Pcircle'))))       %{Poisson}
                  %matrix_frac((reshape(Pcircle,n_l*n_n*tau_max,1,1)-P0),(sum_square_abs(reshape(Pcircle,n_l*n_n*tau_max,1,1)')-P0*P0'))+...
                  %mvnpdf(reshape(Pcircle,n_l*n_n*tau_max,1,1),P0,COV) 
        subject to
            % probability and speed constraints
            Pcircle(~PC.NZIndex) == 0
            0 <= Pcircle(PC.NZIndex) <= 1
            
            % observability constraint
            for oNn = 1:PC.OListLen %size(PC.StartingEdge,2)
                sum( Pcircle(PC.StartingEdge(:,oNn),oNn) ) == 1
            end
            
            % flow constraint
            for jn = 1:PC.OListLen
                sum( Pcircle(PC.InEdge(:,jn),PC.oNode_InIndex(:,jn)),1 ) - ...
                    sum( Pcircle(PC.OutEdge(:,jn),PC.oNode_OutIndex(:,jn)),1 ) >= 0
            end
    cvx_end
    Pcircle = full(Pcircle);
    Phat = reshape(Pcircle,n_l,n_n,tau_max);
    
    %% fix P^t and optimise X^t
    P0 = MToep( Phat,tau_max_x );
    nRow_Ptilde = size(P0,1);
    P0 = P0( n_l*(tau_max-1)+1 : nRow_Ptilde-n_l*(tau_max-1) , :);    
    X0 = reshape(Xhat,Xm*tau_max_x,1);
    
    cvx_begin 
        cvx_precision low
        variables Xcircle(Xm*tau_max_x,1) %[540 by 1]
        maximize sum(sum(X0.*log(Xcircle) - Xcircle)) +...       % {Poisson}
                 sum(sum(Y0(:).*log(P0*Xcircle) - (P0*Xcircle))) % {Poisson}
             
        subject to
            Xcircle >= 0
            
    cvx_end
    Xhat = reshape(Xcircle,Xm,tau_max_x); %[63 by 36]
    
    %% Find Yhat
    Yhat = MVconv( Phat,Xhat );
    Yhat = Yhat(:,tau_max:tau_max_x);
    E_Y(itn) = norm(Y0-Yhat,'fro')/norm(Y0,'fro');
    
    %% display etc.
    if E_Y(itn) <= e_tol
        break;
    end
    if mod(itn,2) == 0
        fprintf('  itn=%d, E_Y=%e, Eratio=%e, time=%f\n',...
            itn,E_Y(itn),(E_Y(itn-1)-E_Y(itn))/E_Y(itn-1),toc);
        tic;
    end
    
end

E_Y = E_Y(1:itn);

end
