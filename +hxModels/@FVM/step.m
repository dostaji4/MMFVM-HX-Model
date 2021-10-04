function [y, t] = step(obj,t , u, Ts)

x_ = obj.x;
z_ = obj.z;
Ts_u_ = obj.Ts_u;
iLim = obj.iLimiter;

[im_dot, iV_dot, ~, iTai] = deal(1,2,3,4); % Input indexes

% --- CFL condition
Ts_max = obj.mw_e / u(im_dot);
assert(Ts <= Ts_max,'CFL condition failed at time %.2f. Sampling time is %f, but the maximum sampling time according to CFL is %f.',t,Ts, Ts_max);

% --- Time integration
switch obj.iIntMethod
    case 1 % RK 1st order1
        k1 = dxdt_rhs(obj,0,x_,u,Ts_u_,z_,iLim);
        x_ = x_ + Ts * k1;
        
    case 2 % RK 2nd order - midpoint (slope at midpoint)
        % http://lpsa.swarthmore.edu/NumInt/NumIntSecond.html
        k1 = dxdt_rhs(obj,0,x_,u,Ts_u_,z_,iLim);
        x1 = x_ + k1*Ts/2;
        k2 = dxdt_rhs(obj,0,x1,u,Ts_u_,z_,iLim); % assumes no dependency on time (because of bad implementation of continuous implementation)
        x_ = x_ + k2*Ts;
        
    case 3 % RK 2nd order - endpoint (endpoint slopes average)
        % http://lpsa.swarthmore.edu/NumInt/NumIntSecond.html
        k1 = dxdt_rhs(obj,0,x_,u,Ts_u_,z_,iLim);
        x1 = x_ + k1*Ts;
        k2 = dxdt_rhs(obj,0,x1,u,Ts_u_,z_,iLim); % assumes no dependency on time (because of bad implementation of continuous implementation)
        x_ = x_ + (k1+k2)/2 * Ts;
        
    case 4 % RK 4th order
        % http://lpsa.swarthmore.edu/NumInt/NumIntFourth.html
        k1 = dxdt_rhs(obj,0,x_,u,Ts_u_,z_,iLim);
        x1 = x_ + k1*Ts/2;
        k2 = dxdt_rhs(obj,0,x1,u,Ts_u_,z_,iLim);
        x2 = x_ + k2*Ts/2;
        k3 = dxdt_rhs(obj,0,x2,u,Ts_u_,z_,iLim);
        x3 = x_ + k3*Ts;
        k4 = dxdt_rhs(obj,0,x3,u,Ts_u_,z_,iLim);
        x_ = x_ + (k1 + 2*k2 + 2*k3 + k4)/6 * Ts;
    otherwise
        error('Undefined integration method "%d".',obj.iIntMethod);
end

% --- Step time
t = t + Ts;

Nw_ = obj.Nw;
Nb_ = obj.Nb;
Tb = x_(Nw_+1:Nw_+Nb_);

% --- Calculate outputs

% Q = 1/Nb_ * feval(@obj.UAba_model,u(iV_dot)*3600).*sum((Tb-repmat(u(iTai),Nb_,1)), 1);
Q = 1/Nb_ * obj.UAba_model(u(iV_dot)*3600).*sum((Tb-repmat(u(iTai),Nb_,1)), 1);
Two = x_(Nw_);

y = [Q; Two];
obj.x = x_;

end