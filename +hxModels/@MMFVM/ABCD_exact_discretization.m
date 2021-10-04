function [A, B, C, D] = ABCD_exact_discretization(obj, u, Ts_x)

[im_dot, iV_dot, iTwi, iTai] = deal(1,2,3,4); % Input indexes
Nw_ = obj.Nw;
Nwb_ = obj.Nwb;
Nb_ = obj.Nb;

% --- Calculation type ----------------------------------------
% When UA models are constant, then the computation of
% expm(Ts(t)*Ac) can be decomposed by eigen-decomposition,
% where the P and P^-1 matrices remain constant and only the
% diagonal matrix varies with flow.
%
% For nonconstant UA models there is no better option then
% using 'expm', so far.

Ad_calculation_type = 'expm';
%             Ad_calculation_type = 'eig_decomp';
%             Bd_calculation_type = 'by ss values';
Bd_calculation_type = 'matrix';

u_old = obj.u_persist;
UAba = obj.UAba_persist;
Ad = obj.Ad_persist;
B = obj.B_persist;

if (norm(u - u_old) > 1e-8)
    
    % Recompute Ad only when there is a change in either m_dot or V_dot (UA change)
    if (norm([u(im_dot);u(iV_dot)] - [u_old(im_dot);u_old(iV_dot)]) > 1e-8)
        
        % --- Heat transfer model
        UAwb = obj.UAwb_model(u(im_dot)*3600);
        UAba = obj.UAba_model(u(iV_dot)*3600);
        
        %% --- State transition "A" matrix calculation -----------------------------
        
        % The state transition matrix is
        %    Ad = phi_p(1,0) = expm(Ts*(Ac1 + 1/2Ac2)),
        % which can be efficiently computed as expm(A) = expm(P*J*P_inv)=
        %    P*expm(J)*P_inv = P*diag(exp(ji))*P_inv
        % Then, Ad = expm(Ts*(Ac1 + 1/2Ac2))= P*diag(exp(Ts*ji))*P_inv
        % Note however that when UA models are not constatn,
        % then the factors change too and the decomposition
        % looses its meaning.
        
        %%% !!!!!!!!!!!!!!!!!!!!!!!!!! TODO  Add Aq to Ac1
        [Ac1,Ac2,~,~,~,Aq] = obj.ABCD_continuous(Ts_x,Nw_,Nwb_,obj.Cw,obj.Cb,UAwb,UAba);
        
        % State transition matrix (for p = 1)
        % Ac = Ac1 + Ac2/2
        Ac = [Ac1, zeros(size(Ac1,1),1); ...
            Aq , 0] ...
            + ...
            [Ac2/2, zeros(size(Ac1,1),1);
            zeros(1,size(Ac2,2))    ,0];
        
        switch Ad_calculation_type
            case 'eig_decomp'
                [P_A,J_A] = eig(full(Ac));
                dJ_A = diag(J_A);
                %                             P_A_inv= inv(P_A);
                
                % if UAs are constant, then the part above is
                % calculated once and this part below will be
                % calculated after change in water flow.
                expJ = exp(Ts_x*dJ_A);
                Ad = P_A*diag(expJ)/P_A;
                %                 Ad = sparse(Ad.*(Ad > 1e-15)); %cut the noise
            case 'expm'
                Ad = expm(Ts_x*Ac);
            otherwise
                error('Unknown Ad matrics computation method: "%s"',Ad_calculation_type);
        end        
        
        obj.Ad_persist = Ad;
        obj.UAba_persist = UAba;
    end
    
    %% --- Input "B" matrix calculation ----------------------------------------
    
    switch Bd_calculation_type
        case 'by ss values'
            
            %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! dava jine vysledky nez MATRIX pro Nwb!=1     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!
            % Formulation by steady state profiles
            [yss,xss] = ss_pde(obj, u);
            Twss = xss(1:Nw_);
            Tbss = xss(Nw_+1:Nw_+Nb_);
            
            Tss = [Twss;
                yss(1);             % add Two
                Tbss;
                yss(2)];            % Qss
            
            Tp0 = [u(iTwi);             % add Twi
                Tss(1:Nw_);
                Tss(Nw_+2:Nw_+1+Nb_);
                0];                 % Q(p=0) = 0
            
            Bd = [zeros(size(Tss)) ,(Tss - Ad*Tp0)*1/u(iTai)];
        case 'matrix'
            % Formulation into Bd*[Twi;Tai] matrix model
            
            [~,~,bw,bb,c,ek] = ss_pde(obj,u);
            
            B1 = [bw ones(size(bw))-bw;
                ek 1-ek;
                bb/(c+1) ones(size(bb))-bb/(c+1)];
            b1q = UAba*(ek-1)/(log(ek)*(c+1))* [1 -1];
            B2 = [1  0
                bw ones(size(bw))-bw;
                bb/(c+1) ones(size(bb))-bb/(c+1)];
            b2q = [0, 0];
            
            Bd = [B1;b1q] - Ad*[B2;b2q];
        otherwise
            error('Unknown Bd computation method: "%s"',Bd_calculation_type);
    end
    
    B = [1,0;Bd];
    obj.B_persist = B;
    obj.u_persist = u;
end
%%
S = spdiags(ones(Nw_+Nb_+2,1),-1,Nw_+Nb_+3,Nw_+Nb_+2);
iRow = 1:Nw_+Nb_+1; iCol = [1:Nw_+1,Nw_+3:Nw_+Nb_+2];
K = sparse(iRow,iCol,1,Nw_+Nb_+2,Nw_+Nb_+3);
A = S*Ad*K;
C = sparse([1,2],[Nw_+Nb_+3, Nw_+2],1,2,Nw_+Nb_+3); % [Q;Two]
D = sparse(2,2);
end

function [ v ] = cf_matrixexponential( A,b )
% cf_matrixexponential
% computes v = exp(A)*b where A is a real, symmetric and negative semidefinite
% The vector b is real
% This function is much faster than expm(A)*b as it avoids the explicit
% computation of exp(A).
% Thomas Schmelzer, December 2008
v = cf_matrixfunction(@exp,A,b,14);
end

function [v] = cf_matrixfunction(fct,A,b,n,varargin)
% compute the product v = f(A)*b without the explicit computation of f(A)
% The matrix A has to be negative semidefinite (that is, all eigenvalues
% should be negative)

% First the method computes a Caratheodory-Fejer approximation, that is, a
% a uniform rational approximation of the function f of type (n,n) on
% the negative real axis.

% The matrix A has to be real and n should be even
% The vector b has to be real, too.

[zi,ci,r_inf] = cf_realaxis(fct,n,varargin{:});

% identity matrix
I = speye(size(A));
v = real(r_inf*b);
for i=1:n/2
    v = v + 2*real(ci(i)*((A-zi(i)*I)\b));
end
end

% Caratheodory-Fejer method for constructing a uniform rational approximation
% of the function f of type (n,n) on the negative real axis.
% This is a modified version of a code published in
% Trefethen, Weideman, Schmelzer
% Talbot Quadratures and Rational Approximations
% BIT, 2006

% zi      -- poles of the rational function
% ci      -- residues associated with those poles
% r_inf   -- r at infinity
% decay   -- SVD of the Hankel matrix associated with the problem
% s       -- error introduced by the (n,n) approximation.
function [zi,ci,r_inf,s,decay] = cf_realaxis(fct,n,varargin)

K = 75;                                 % no of Cheb coeffs
nf = 1024;                              % no of pts for FFT
w = exp(2i*pi*(0:nf-1)/nf);             % roots of unity
t = real(w);                            % Cheb pts (twice over)
scl = 9;                                % scale factor for stability
F = zeros(size(t));
g = (t~=-1);
F(g) = feval(fct,scl*(t(g)-1)./(t(g)+1),varargin{:});
c = real(fft(real(F)))/nf;              % Cheb coeffs of F
f = polyval(c(K+1:-1:1),w);             % analytic part f of F
[U,S,V] = svd(hankel(c(2:K+1)));        % SVD of Hankel matrix
s = S(n+1,n+1);                         % singular value
decay = diag(S);                        % decay of the error
u = U(K:-1:1,n+1)'; v = V(:,n+1)';      % singular vectors
zz = zeros(1,nf-K);                     % zeros for padding
b = fft([u zz])./fft([v zz]);           % finite Blaschke product
rt = f-s*w.^K.*b;                       % extended function r-tilde
zr = roots(v); qj = zr(abs(zr)>1);      % poles
qc = poly(qj);                          % coeffs of denominator
pt = rt.*polyval(qc,w);                 % numerator
ptc = real(fft(pt)/nf);                 % coeffs of numerator
ptc = ptc(n+1:-1:1); ci = 0*qj;
for k = 1:n                             % calculate residues
    q = qj(k); q2 = poly(qj(qj~=q));
    ci(k) = polyval(ptc,q)/polyval(q2,q);
end
zi = scl*(qj-1).^2./(qj+1).^2;          % poles in z-plane
ci = 4*ci.*zi./(qj.^2-1);               % residues in z-plane
r_inf = 0.5*(feval(fct,0,varargin{:})+sum(ci./zi));    % r at infinity

[m, order]=sort(imag(zi));              % sort poles
zi=zi(order);                           % reorder poles
ci=ci(order);                           % reorder residues
end
