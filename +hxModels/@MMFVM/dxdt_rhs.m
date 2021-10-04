function dxdt = dxdt_rhs(obj,t,x,u,Ts_u,tau)
% dxdt = f(x,u) of y = h(x,u)
[im_dot, iV_dot, iTwi, iTai] = deal(1,2,3,4); % Input indexes
idx_Tw_ = obj.idx_Tw;
Nw_ = obj.Nw;

% --- Get temperature state and shift ---------------------
T = x([1:idx_Tw_(end),obj.idx_Tb]); 
p = x(obj.idx_p);

% --- Calculate inputs -----------------------------------------
k_v = floor(t/Ts_u) + 1; % index for velocity
% k_b = floor((t-tau)/Ts_u) + 1; % index for boundary conditions
k_b = k_v;

% --- Heat exchange coefficients --------------------------
UAwb = obj.UAwb_model(u(k_b,im_dot)*3600);
UAba = obj.UAba_model(u(k_b,iV_dot)*3600);

% %% HACK - assumes constant flow during initial
% step_number = u(k,im_dot)*t/ ((obj.mw_e*Nw_)/u(k,im_dot)) ;
% if step_number <= 1
%     Nw_s = 2*Nw_;
% elseif step_number == Nw_
%     Nw_s = 2*Nw_;
% else
%     Nw_s = Nw_;
% end

Ts_y = ((obj.mw_e*Nw_)/Nw_)/u(k_v,im_dot);  % = 1/dpdt

% --- Derivatives -----------------------------------------
[Ac1,Ac2,Bc,~,~,Aq,Bq] = obj.ABCD_continuous(Ts_y,Nw_,obj.Nwb,obj.Cw,obj.Cb,UAwb,UAba);

% dT = (Nw_+1)/Nw_ *((Ac1 + Ac2*p)*T + Bc*u(k_b,iTai));
% dQdt = (Nw_+1)/Nw_ * (Aq*T          + Bq*u(k_b,iTai));
% dpdt = (Nw_+1)/Nw_ * 1/Ts_y; % (m_dot/mw)*Nw
% dT = diag([(Nw_+1)/Nw_ *ones(Nw_+1,1);ones(Nw_,1)]) *((Ac1 + Ac2*p)*T + Bc*u(k_b,iTai));
% dQdt = (Nw_+1)/Nw_ * (Aq*T          + Bq*u(k_b,iTai));
% dpdt = (Nw_+1)/Nw_ * 1/Ts_y; % (m_dot/mw)*Nw
dT = ((Ac1 + Ac2*p)*T + Bc*u(k_b,iTai));
dQdt = (Aq*T          + Bq*u(k_b,iTai));
dpdt = 1/Ts_y; % (m_dot/mw)*Nw

% T = [Twi,Tw,Tb]
dTwidt = dT(1);
dTwdt = dT(idx_Tw_);
dTwodt = 0; 
dTbdt = dT(idx_Tw_(end)+1:end);

% x = [Twi; Tw; Two; Tb; Q; p]
dxdt = [dTwidt; dTwdt; dTwodt; dTbdt; dQdt; dpdt];
end