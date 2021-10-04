function dxdt = dxdt_rhs(obj,t,x,u,Ts_u,z,iLim)

if obj.isDisplayProgress
    %             if tn < floor(t/10)
    %                 disp(num2str(t));
    %                 tn = floor(t/10);
    %             end
end

[im_dot, iV_dot, iTwi, iTai] = deal(1,2,3,4); % Input indexes

Cw_ = obj.Cw;
Cb_ = obj.Cb;
mw_ = obj.mw;

n = numel(z);
x1 = x(1:n,1);
x2 = x(n+1:2*n,1);

1;
% --- Index of the current input
iInp = floor(t/Ts_u) + 1;

%         % impose boundary condition (first two elements on the left)
%         if strcmp(method,'KT')
%             x1(1,1) = Twi(iInput); % Tw(1) = Twi
%             x1(n,1) = x1(n-1,1);   % Tw(n) = Tw(n-1)
%             x2(1,1) = x2(2,1);     % Tb(1) = Tb(2)
%             x2(n,1) = x2(n-1,1);   % Tb(n) = Tb(n-1)
%         end

% --- Spatial derivatives -----------------------------------------
% --- Flux term
flux_fun = @(t,x) x;
dx1dz = flux(n,z,t,x1,flux_fun,u(iInp,iTwi),iLim);

%
% --- Get actual heat exchange coefficients -------------------
h_wb = obj.UAwb_model(u(iInp,im_dot)*3600);
h_ba = obj.UAba_model(u(iInp,iV_dot)*3600);

dx1dt = -u(iInp,im_dot)/mw_ * dx1dz - h_wb/Cw_*(x1 - x2);
dx2dt = 1/(Cb_) * (+h_wb*(x1-x2) - h_ba*(x2 - u(iInp,iTai)));

dxdt = [dx1dt;dx2dt];

    function [fz] = flux(n,z,t,y,flux_fun,yB,iLimiter)        
        %{
        %...  function flux returns the first derivative, fz, of a
        %...  flux function f(x) over the spatial domain z0 < z < zL from Van Leer's
        %...  slope limiter approximations.
        %...
        %...  The basic principle of slope limiters is to approximate the solution
        %...  x by combining a first-order, oscillation-free, monotone scheme with
        %...  a higher-order scheme in regions where the solution is sufficiently
        %...  smooth. In region of high-solution gradients, the first-order scheme
        %...  is prefered.
        %...
        %...  The starting point of this procedure is the finite volume method in
        %...  which the spatial domain is subdivided into cells
        %...
        %...    [z(i)-dz(i)/2  z(i)+dz(i)/2]
        %...
        %...  In a first-order approximation x(i) is supposed to be constant on
        %...  the local domain.
        %...
        %...           |           |                 |              |
        %...           |           |                 |              |
        %...           |           |                 |              |
        %...           |           |                 |     x(i+1)   |
        %...           |           |                 |--------------|
        %...           |           |                 |              |
        %...           |           |                 |              |
        %...           |           |                 |              |
        %...           |           |      x(i)       |              |
        %...           |           |-----------------|              |
        %...           |           |                 |              |
        %...           |           |                 |              |
        %...           |           |                 |              |
        %...           |           |                 |              |
        %...           |           |                 |              |
        %...           |  x(i-1)   |                 |              |
        %...           |-----------|                 |              |
        %...        ___|___________|_________________|______________|_______
        %...                   z(i)-dz(i)/2       z(i)+dz(i)/2
        %...
        %...   The PDE
        %...
        %...    x_t = - f_z(x,t)                                                                  (1)
        %...
        %...  is integrated by the method of lines following :
        %...
        %...    dx(i)     [f(x(i+1/2)) - f(x(i-1/2))]
        %...    ----- = - ---------------------------
        %...     dt                 dz(i)
        %...
        %...  where
        %...
        %...    x(i+1/2) = x(z(i) + dz(i)/2)
        %...
        %...    x(i-1/2) = x(z(i) - dz(i)/2)
        %...
        %...
        %...  However, the approximation of the solution to the piecewise constant
        %...  values x(i) is only first order, and in order to construct a higher
        %...  (second) order approximation, the solution will be approximated by
        %...  piecewise linear functions. The slopes of those functions will be
        %...  limited by the slope limiter phi given in the vanleer case by
        %...
        %...  phi = (r + |r|)/(1 + r)
        %...
        %...  where r is the ratio of two consecutive solution derivatives :
        %...
        %...  if the solution flows from left to right (f_x > 0), we have
        %...
        %...             x(i+1)-x(i)
        %...             -----------
        %...             z(i+1)-z(i)
        %...    r(i) = ---------------
        %...             x(i)-x(i-1)
        %...             -----------
        %...             z(i)-z(i-1)
        %...
        %...  It can be shown that x(i-1/2) and x(i+1/2) (called xL and xR in
        %...  the code) are given by
        %...
        %...  if f_x > 0 :
        %...
        %...                         phi(i-1)*dz(i-1)
        %...    x(i-1/2) = x(i-1) +  ----------------*(x(i-1)-x(i-2))
        %...                         dz(i-1)+dz(i-2)
        %...
        %...                       phi(i)*dz(i)
        %...    x(i+1/2) = x(i) +  -------------*(x(i)-x(i-1))
        %...                       dz(i)+dz(i-1)
        %...
        %...
        %...  argument list
        %...
        %...     z          independent variable (input)
        %...
        %...     n          number of grid points in the z domain including the
        %...                boundary points (input)
        %...
        %...     t          time (input)
        %...
        %...     x          dependent variable (input)
        %...
        %...     flux       matlab function which computes f(x)
        %...                call: f = fun(t,x) where where flux = 'fun'     (input)
        %...
        %}
        
        dz = zeros(n,1);
        fz = zeros(n,1);                
        yR = zeros(n,1);
        yL = zeros(n,1);
        
        % --- Grid differences
        %                                      0  z(1) dz(1) z(2)
        dz(1) = ((z(1)-0) + (z(2)-z(1))/2); %  |   *    |     *
        
        %                                     z(i-1) dz(i)/2 z(i) dz(i)/2   z(i+1)
        dz(2:n-1) = (z(3:n)-z(1:n-2))/2;    %   *    |        *          |    *
        
        %                                        z(n-1)   dz(n)/2   z(n) dz(n)/2 1
        dz(n) = ((z(n)-z(n-1))/2 + (1-z(n)));  %   *     |            *          |
        
        assert(all((dz(1)- dz)<eps));
        
        % --- Computation of ratio of succesive gradients r(i) and
        %     a slope limiter function phi(i) for an elements right face
        zn = [z(1)-dz(1);z];
        yn = [yB;y];
        phi = hxModels.hxModelCommon.flux_limiter(zn,yn,iLimiter);
        
        % --- Flux calculation for all elements
        % Flux ~ first order approximation + in smooth regions (phi~=0) second order
        %        approximation
        for i = 1:n
            % --- Computation of Right face value
            % https://en.wikipedia.org/wiki/Flux_limiter
            % yR = yR_1order + phi*(yR_2order - yR_1order)
            % yR_1order = y(i) %upwind constant extrapolation
            % yR_2order = 3/2*y(i) - 1/2*y(i-1) % upwind linear extrapolation
            % => yR = y(i) + phi(i) * (3/2*y(i)-1/2*y(i-1) - y(i)) =
            %       = y(i) + phi(i) * 1/2 * (y(i) - y(i-1))
            % where phi(i) = limiter_func(r(i)) and
            %       r(i) = (y(i+1)-y(i))/(y(i)-y(i-1)
            %                 yR(i) = y(i) + phi(i) * (y(i) - y(i-1)) * dz(i) / (dz(i)+dz(i-1));
            
            yR(i) = yn(i+1) + phi(i+1)/2 * (yn(i+1) - yn(i));
            
            if i == 1
                yL(i) = yB; % Add boundary value here
                fxL = feval(flux_fun,t,yL(i));
            else
                yL(i) = yR(i-1);
            end
            % --- Calculate flux for an element
            % https://en.wikipedia.org/wiki/MUSCL_scheme
            % Godunov scheme uses piecewise constant first order approximation
            % Higher order approximation try to reconstruct averaged state
            %  f(z)/dz ~ (f(xR) - f(xL)) / dz
            fxR = feval(flux_fun,t,yR(i));
            fz(i) = (fxR - fxL) / dz(i);
            % Save f(xL) for next iteration
            fxL = fxR;
        end
    end

    function [fz] = KT_centered_limiter_fz(ne,n,z,t,x,flux,dflux_dx)
        %...
        %...  function KT_centered_limiter returns the first derivative, fz,
        %...  of a flux-function vector f(x) over the spatial domain z0 < z < zL.
        %...
        %...  This function implements a centered scheme for convection-diffusion PDEs
        %...  proposed by Kurganov and Tadmor (2000), which can be used as a general
        %...  finite-volume method, independently of the eigenstructure of the problem.
        %...  The only required information is an estimation of the velocity of the wave.
        %...
        %...  argument list
        %...
        %...     ne         number of equations (input)
        %...
        %...     n          number of grid points in the z domain including the
        %...                boundary points (input)
        %...
        %...     z          independent variable (input) : z(n)
        %...
        %...     t          time (input)
        %...
        %...     x          dependent variable (input) : x(n,ne)
        %...
        %...     flux       matlab function which computes f(x) : f(n,ne)
        %...                call: f = fun(t,x) where flux = 'fun'     (input)
        %...
        %...     dflux_dx   matlab function which computes the derivative of f with
        %...                respect to x : dfdx(n,ne,ne)
        %...                call: f_x = fun(t,x) where dflux_dx = 'fun'
        %...
        %...
        for i=1:ne
            xtmp(1:n) = x(1:n,i);
            
            xtmpz(1) = (xtmp(2)-xtmp(1))/(z(2)-z(1));
            
            xtmpz1(1:n-2) = 2*(xtmp(3:n)-xtmp(2:n-1))./(z(3:n)-z(2:n-1))';
            xtmpz2(1:n-2) = (xtmp(3:n)-xtmp(1:n-2))./(z(3:n)-z(1:n-2))';
            xtmpz3(1:n-2) = 2*(xtmp(2:n-1)-xtmp(1:n-2))./(z(2:n-1)-z(1:n-2))';
            testsgn       = sign(xtmpz1)+sign(xtmpz2)+sign(xtmpz3);
            for j=1:n-2
                if (testsgn(j) == 3) || (testsgn(j) == -3)
                    xtmpz(j+1)=sign(testsgn(j))*min([abs(xtmpz1(j)) abs(xtmpz2(j)) abs(xtmpz3(j))]);
                else
                    xtmpz(j+1)=0;
                end
            end
            
            xtmpz(n)    = (xtmp(n)-xtmp(n-1))/(z(n)-z(n-1));
            
            xtp(1:n-1)  = xtmp(2:n)-xtmpz(2:n).*(z(2:n)-z(1:n-1))'/2;
            xtm(1:n-1)  = xtmp(1:n-1)+xtmpz(1:n-1).*(z(2:n)-z(1:n-1))'/2;
            
            xp(1:n-1,i) = xtp(1:n-1);
            xm(1:n-1,i) = xtm(1:n-1);
        end
        
        dfludxp = feval(dflux_dx,t,xp,ne);
        dfludxm = feval(dflux_dx,t,xm,ne);
        
        
        for i=1:n-1
            size(dfludxp);
            ajac(1:ne,1:ne) = dfludxp(i,1:ne,1:ne);
            bajac           = balance(ajac);
            vp(1:ne)        = eig(bajac);
            ajac(1:ne,1:ne) = dfludxm(i,1:ne,1:ne);
            bajac           = balance(ajac);
            vp(ne+1:2*ne)   = eig(bajac);
            aspeed(i,1)     = max(abs(vp));
        end
        
        flzp = feval(flux,t,xp,ne);
        flzm = feval(flux,t,xm,ne);
        
        for j=1:ne
            fz(1,j)=0;
            fz(2:n-1,j)=(flzp(2:n-1,j)-flzp(1:n-2,j)+flzm(2:n-1,j)-flzm(1:n-2,j)...
                -aspeed(2:n-1,1).*(xp(2:n-1,j)-xm(2:n-1,j))...
                +aspeed(1:n-2,1).*(xp(1:n-2,j)-xm(1:n-2,j)))./(z(3:n)-z(1:n-2));
            fz(n,j)=0;
        end
    end
end