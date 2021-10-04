function [Ac1,Ac2,Bc,Cc,Dc,Aq,Bq] = ABCD_continuous(Ts_x,nl,nlb,Cw,Cb,UAwb,UAba)
% Advected finite volume method model of a form
%
%   dTdt = (Ac1 + Ac2*p)T + Bc*u
%   dpdt = nonlinear function of m_dot, mw
%   dQdt = Aq*T + Bq*u,
%
% where T = [Tw;Tb] and u = Tai

betaL= UAwb/Cw;
betaA = UAba/Cb;
betaB = 1/nlb*UAwb/Cb;

nb = nl/nlb;
N = nl+1;   % there is nl+1 elements in the pipe in the continuous model!
M = nb; 

Im = speye(M);
Z1M = sparse(1,M);

% A1 matrix Tw-to-Tw (= blkdiag(0,-speye(N-1)*betaL) );
idx = (1:N)';
A1 = sparse(idx(2:end),idx(2:end),-betaL); %ok
% A2 matrix Tw-to-Tw (= blkdiag(-betaL, sparse(N-2,N-2), betaL);    
A2 = sparse(idx,idx,[-betaL,zeros(1,N-2),betaL]); %ok
    
if  nlb == 1  
    % Concatenation is too slow. Reformulated by sparse indexing.
%     B1 = [Z1M; Im*betaL];
    B1_ = sparse(idx(2:end),idx(1:end-1),1,N,M); 
    B1  = B1_ * betaL; %ok
    
%     B2 = [Im*betaL;Z1M] - B1; 
    B2_ = sparse(idx(1:end-1),idx(1:end-1),1,N,M) - B1_; %ok    
    B2  = B2_ * betaL; %ok
    
%     C1 = [Z1M' Im*betaB];
    C1 = B1_' * betaB; %ok
    
%     C2 = [Im*betaB Z1M'] - C1;
    C2 = B2_' * betaB; %ok
    
    D = Im*(-betaB - betaA);
else    
    % Kron formulation is too slow. Reformulated by sparse indexing.
    
%     B1 = [kron(Im,[0; ones(nlb-1,1)*betaL]); Z1M] + [Z1M; kron(Im,[Z1nlb'; betaL])];
    iRow = (2:N)'; iCol = reshape(repmat(1:M,nlb,1),nlb*M,1);
    B1_ = sparse(iRow, iCol, 1); %ok
    B1  = B1_*betaL;
    
%     B2 = [kron(Im,[betaL; Z1nlb_m1']); Z1M] + [Z1M; kron(Im,[Z1nlb_m1'; -betaL])];
    iRow = (1:nlb:N-nlb)';   iCol = (1:M)';
    kr1 = sparse(iRow, iCol, 1,N,M); 
    iRow = (nlb+1:nlb:N)';   iCol = (1:M)';
    kr2 = sparse(iRow, iCol, -1,N,M); 
    B2_ = kr1 + kr2; 
    B2  = B2_ * betaL; %ok
        
%     C1 = [kron(Im,[0 ones(1,nlb-1)*betaB]) Z1M'] + [Z1M' kron(Im,[Z1nlb_m1 betaB])];
    C1 = B1_' * betaB; %ok
    
%     C2 = [kron(Im,[betaB Z1nlb_m1]) Z1M'] + [Z1M' kron(Im,[Z1nlb_m1 -betaB])];
    C2 = B2_' * betaB; %ok
    
    D = Im*(-nlb*betaB - betaA);
end

% Temperature 
Ac1 = [A1 B1; C1 D];
Ac2 = [A2 B2; C2 sparse(M,M)];
Bc = [sparse(N,1); ones(M,1)*betaA];
Cc = [sparse(1,N-1), 1, sparse(1,M)];     % Two = Tw(end)     
Dc = 0;    

% % Position
% dpdt = m_dot/(mw/nl);  % = 1/Ts_x;

% Heat flow
Aq = [sparse(1,N), 1/Ts_x*UAba/nb*ones(1,M)];   % dQdt = sum(UAba/nb*(Tb - Tai));
Bq = -1/Ts_x*UAba;

end
