function y = step(obj,u,Ts)

% --- Input checking and initialization -----------------------
if ~obj.isInitialized
    error('Run initialization first.');
end
assert(numel(u) == 4,'Input dimension must be 4 ([u(m_dot),u(V_dot),u(Twi),u(Tai)]).');
[im_dot, iV_dot, iTwi, iTai] = deal(1,2,3,4); % Input indexes

% --- Get actual heat exchange coefficients -------------------
UAwb = obj.UAwb_model(u(im_dot)*3600);
UAba = obj.UAba_model(u(iV_dot)*3600);

% --- Get model matrices --------------------------------------
[A,B,C,D] = obj.ABCD(u(im_dot),...
    obj.Nw, obj.Nwb, obj.mw, obj.Cw, obj.Cb, UAwb, UAba, Ts);

% --- Update --------------------------------------------------
xkp1 = A * obj.x + B * [u(iTwi); u(iTai)];
y    = C * obj.x + D * [u(iTwi); u(iTai)]; %[Two; Q; Qmax]
obj.x = xkp1;

y = [y(2), y(1)]; %[Q; Two]

% Two = y(1);
% Q = y(2);
% Qmax = y(3);
% 
% Tao = u(Tai) + Q / (u(V_dot)*obj.cv_std_air);
% if Qmax > 0
%     Q_Qmax = Q/Qmax;
% else, Q_Qmax = 0;
% end

end