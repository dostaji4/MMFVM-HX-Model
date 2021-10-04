function [Q, Two, t_y] = run(obj, m_dot, V_dot, Twi, Tai, Ts_u)

if strcmp(obj.simulation_mode,'discrete')
    % --- Call the common run method
    [Q, Two, t_y] = run@hxModels.hxModelCommon(obj, m_dot, V_dot, Twi, Tai, Ts_u);    
    return
end

% === CONTINUOUS MODE HERE AFTER ========================================

% --- Input checks and routines
narginchk(6,6);
tmp = {m_dot, V_dot, Twi, Tai};
nData = max(cellfun(@numel,tmp));
for iInput = 1:numel(tmp)
    if numel(tmp{iInput}) == 1
        % Replicate nData times
        tmp{iInput} = tmp{iInput} * ones(nData,1);
    elseif numel(tmp{iInput}) == nData
        % Ensure column vector
        tmp{iInput} = tmp{iInput}(:);
    else
        error('Incorrect dimension of input #%d "%s". Specify either a scalar value or a vector of length %d.',iInput,inputname(iInput),nData);
    end
end
u = cell2mat(tmp); % size = [time, 4]

% --- Initialization
obj.initialize;
Nwb_ = obj.Nwb;
Nb_ = obj.Nb;
assert(Nwb_ == 1,'FVM defined only for Nwb=1.');

obj.Ts_u = Ts_u;
if isempty(obj.Ts_y)
    obj.Ts_y = obj.Ts_y_max;
end

% --- Grid with boundary points (no boundary points except 'KT', handled directly)
z_ = [  obj.z(1:end)];
x0 = [interp1(obj.z,obj.Tw0,z_,'pchip','extrap');interp1(obj.z,obj.Tb0,z_,'pchip','extrap')];
n = numel(z_);
didx = numel(z_) - numel(obj.z);
idxTw = (1+didx):n;
idxTb = n+(1+didx):2*n;
idxTwo = n;
    
if strcmp(obj.flux_limiter,'KT')
    % KT
    % --- Grid with boundary points (only at the left boundary, flow assumed always positive)
    dz_b = dz_o/10; %boundary size (do not exceed 'dz');
    
    z_ = [dz_b/4; z_o ; 1-dz_b/4]; %add point before and after
    x0 = [interp1(z_o,obj.Tw0,z_,'pchip','extrap');interp1(z_o,obj.Tb0,z_,'pchip','extrap')];
    n = numel(z_);
    idxTw = 2:n-1;
    idxTb = n+2:2*n-1;
    idxTwo = n-1;
end

% --- SIMULATE ------------------------------------------------------------
if obj.isDisplayProgress
    tic
    fprintf('## Started simulation of "%s" model \n',obj.model_type);
end

t_end = Ts_u * (nData - 1);
time = 0:obj.Ts_y:t_end;
iLim = obj.iLimiter;

% --- Run the simulation --------------------------------------    
opt = odeset('RelTol',1e-6,'NormControl','on');
[t_y,xout] = ode45(@(t,x)dxdt_rhs(obj,t,x,u,Ts_u,z_,iLim),time,x0,opt);


if obj.isDisplayProgress
    t_elapsed = toc;
    fprintf('\tDone. Time elapsed %.2f s.\n',t_elapsed);
end

if obj.Ts_y ~= Ts_u
    t_u = 0:Ts_u:(nData-1)*Ts_u;
    Tai_r = interp1(t_u,Tai,t_y);
    V_dot_r = interp1(t_u,V_dot,t_y);
else
    Tai_r = Tai(:);
    V_dot_r = V_dot(:);
end

Tb = xout(:,idxTb); % omit the boundary points
Tw = xout(:,idxTw);
Two = xout(:,idxTwo);
Q = 1/Nb_ * obj.UAba_model(V_dot_r*3600).*sum((Tb-repmat(Tai_r,1,Nb_)), 2);

obj.saved.x.addArray([Tw,Tb]);
obj.saved.y.addArray([Q,Two]); %Nw_ because of KT method
obj.saved.t.addArray(t_y);    
end

