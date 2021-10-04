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
Nw_ = obj.Nw;
mw_e_ = obj.mw_e;

obj.Ts_u = Ts_u;
if isempty(obj.Ts_y)
    obj.Ts_y = obj.Ts_y_max;
end
t_end = Ts_u * (nData - 1);
time_equidistant = 0:obj.Ts_y:t_end;

% --- Prepare integrator
integrator = @ode23; %should always use 23, 45 not accurate
RelTol = 1e-8;
AbsTol = 1e-3;
optionsMMFVM = odeset('RelTol',RelTol,...
    'AbsTol',AbsTol,...
    'NormControl','on',...    'OutputFcn',@outFcn,... %manually called
    'Events',@(t,x) zeroEvents(t,x));
optionsMMFVM.MaxStep = m_dot(1)/(obj.mw/Nw_);
optionsMMFVM.InitialStep = optionsMMFVM.MaxStep/2;

optionsFVM = optionsMMFVM;
optionsFVM.Events = [];

obj.FVMmodel.Ts_u = Ts_u;

% --- SIMULATE ------------------------------------------------------------
if obj.isDisplayProgress
    tic
    fprintf('## Started simulation of "%s" model \n',obj.model_type);
end

% --- Initial state
% x = [Twi; Tw; Two; Tb; Q; p]
x_ = [Twi(1);...      Twi
    obj.Tw0;...       Tw
    obj.Tw0(end);...  Two    
    obj.Tb0;...       Tb
    0;...             Q
    0];%              p = 0
obj.x = x_;

% --- Saving
alloc_block = min(1000,nData);
xs = zeros(obj.nx,alloc_block);
ys = zeros(obj.ny,alloc_block);
ts = zeros(1, alloc_block);
alloc_size = alloc_block; % actual allocated size
nSaved = 0;

t_sim = 0;
tstart = t_sim;
tau = 0;
t_half = 0;
tau_buffer = zeros(2*(Nw_),1);
% Save initial
outFcn(0,x_,'init');

while t_sim < t_end
    
    % Input index
    k = floor(t_sim/obj.Ts_u) + 1;
    
    % --- Select sampling
    % Prepare equidistant time
    idx_t = find(t_sim < time_equidistant,1);
    time_eq = [tstart, time_equidistant(idx_t:end)];
    
    switch obj.sampling_mode
        case 'equidistant'
            time = time_eq;
        case 'output'
            time = [tstart t_end]; % will get chopped by zeroEvents
    end
    if k == 199 && strcmp(obj.sampling_mode,'output')
        1;        
    end
    Ts_c = mw_e_/m_dot(k); 
    if (Ts_c > obj.Ts_y_max) 
        
        obj.iteratingFVM = 1;
        
        if obj.switch_to_FVM == 1
            obj.switch_to_FVM = 0;
            % --- Regrid water elements
            [zn] = hxModels.hxModelCommon.equidistant_grid(obj.FVMmodel.Nw); %positions of             
            
            %for grid scaling uncomment the following 2 lines
%             zn_ = zn * (1+obj.dz/2);  
%             zn = zn_;

            [Twn_e] = hxModels.MMFVM.grid_value_interpolation(x_(1:obj.idx_Two),...
                x_(obj.idx_p),obj.iInterpolation_scheme,...
                obj.z, t_sim, zn);
            
            % --- Regrid body elements
            Tb = x_(obj.idx_Tb);
            Tbn_e = interp1(obj.z,Tb,zn,'makima','extrap');
            
            % x_FVM = [Tw;Tb]
            obj.FVMmodel.x = [Twn_e;Tbn_e];
                        
            %{
            % refactor back and make energy bilance  
%             Twn = hxModels.hxModelCommon.refactor(Twn_e',...
%                 zn_,zn,'average')';
%             Tbn = hxModels.hxModelCommon.refactor(Tbn_e',...
%                 zn_,zn,'average')';
            
            Tw_state = x_(1:obj.idx_Tw(end));
            Tw = hxModels.MMFVM.grid_value_interpolation(x_(1:obj.idx_Two),...
                x_(obj.idx_p),obj.iInterpolation_scheme,...
                obj.z, t_sim);
                
            [im_dot, iV_dot, ~, iTai] = deal(1,2,3,4); % Input indexes
                      
            Qwb = 1/numel(Tw) * obj.UAwb_model(u(k,im_dot)*3600).*sum(Tw-Tb,1);                        
            Qwbn_m = 1/numel(Twn_e) * obj.UAwb_model(u(k,im_dot)*3600).*sum(Twn_e-Tbn_e,1);
            fprintf('Qwb: Old - interpolated = %.2f W. \n',Qwb-Qwbn_m);
            % Twn_ = averaged back from the interpolated Twn to obj.z
            % positions
            Twn_a = hxModels.hxModelCommon.refactor(Twn_e',...            
                    zn,obj.z,'average');
            Qwbn_a = 1/numel(Twn_a) * obj.UAwb_model(u(k,im_dot)*3600).*sum(Twn_a'-Tb,1);
            fprintf('Qwb: Old - averaged interpolated = %.2f W. \n',Qwb-Qwbn_a);
                        
            Qba = 1/numel(Tb) * obj.UAba_model(u(k,iV_dot)*3600).*sum(Tb-repmat(u(k,iTai),numel(Tb),1),1);
            Qban_m = 1/numel(Tbn_e) * obj.UAba_model(u(k,iV_dot)*3600).*sum(Tbn_e-repmat(u(k,iTai),numel(Tbn_e),1),1);            
            fprintf('Qba: Old - interpolated = %.2f W. \n',Qba-Qban_m);
            % Tbn_a = averaged back from the interpolated Tbn to obj.z
            Tbn_a = hxModels.hxModelCommon.refactor(Tbn_e',...            
                    zn,obj.z,'average');
            Qban_a = 1/numel(Tbn_a) * obj.UAba_model(u(k,iV_dot)*3600).*sum(Tbn_a'-repmat(u(k,iTai),numel(Tbn_a),1),1);            
            fprintf('Qba: Old - averaged interpolated = %.2f W. \n',Qba-Qban_a);
            
            figure(11); cla
            plot(obj.z,Tb,'^'); hold on
            plot(zn,Tbn_e,'+-');            
            plot(obj.z,Tbn_a,'o');
            
            plot([-obj.dz/2;obj.z] + x_(obj.idx_p)*obj.dz,Tw_state,'^');
%             plot(obj.z,Tw,'*-'); hold on
            plot(zn,Twn_e,'+-');            
            plot(obj.z,Twn_a,'o');
            %}          
        end

        % --- Extended FVM grid
        [zn] = hxModels.hxModelCommon.equidistant_grid(obj.FVMmodel.Nw); %positions of         
        
        % --- Step the internal FVM model - always equidistant        
        Ts_CFL_max = obj.FVMmodel.mw_e/m_dot(k); % CFL limit given time step
        Ts = min(obj.Ts_y,Ts_CFL_max); % CFL states maximal sampling to be Ts_c        
        [y_,t_sim] = step(obj.FVMmodel, t_sim, u(k,:),Ts);   
                
        % Plot profile
        hxModels.MMFVM.grid_value_interpolation(obj.FVMmodel.x(obj.FVMmodel.idx_Tw),...
            0,obj.iInterpolation_scheme,...
            obj.FVMmodel.z, t_sim);
        
        % --- Get values at MMFVM.z positions             
        TwFVM = obj.FVMmodel.x(obj.FVMmodel.idx_Tw);
        TbFVM = obj.FVMmodel.x(obj.FVMmodel.idx_Tb);
        Tw_ = hxModels.hxModelCommon.refactor(TwFVM',...            
            zn,obj.z,obj.plot_water_state_method);
        Tb_ = hxModels.hxModelCommon.refactor(TbFVM',...
            zn,obj.z,'average');
            
        % x = [Twi; Tw; Two; Tb; Q; p]            
        [~, iV_dot, iTwi, iTai] = deal(1,2,3,4);
        x_ = [u(k,iTwi);...
            Tw_';...
            y_(2);... 
            Tb_';...
            y_(1);... %Q out of FVM not valid
            0];         %p = 0
                      
        % --- Save
        outFcn(t_sim, x_,'');
                
        % --- Prepare new run
        k = floor(t_sim/obj.Ts_u) + 1;
        Ts_c = mw_e_/m_dot(k); 
        if Ts_c <= obj.Ts_y_max
           x_(obj.idx_Q) = 0; % reset heat output (integrates)
        end
        tstart = t_sim;
    else        
        assert(Ts_c <= obj.Ts_y_max);
        obj.iteratingFVM = 0;
        
        % --- Step the model
        f = @(t,x)dxdt_rhs(obj,t,x,u,Ts_u,tau);
        [ts_,~,te,xe] = integrator(f,time,x_,optionsMMFVM);        
        
        if numel(te) < 2
            break %time has reached an end
        else
            t_sim = te(2);
            x_ = xe(2,:)';
            t_half = te(1);
            tau = t_sim - time(1);   
            
           % Save
           outFcn(t_sim,x_,'');
        end
        
        if obj.switch_to_FVM == 0
            assert(abs(x_(obj.idx_p)-1) < AbsTol); % p != 1 after every integration
            
            % A good guess of a valid first timestep is the length of the last valid
            % timestep, so use it for faster computation.  'refine' is 4 by default.
            % Note: do not use odeset, it is incredibly slow!
            if (ts_(end)-ts_(end-1))>0
                optionsMMFVM.InitialStep = ts_(end)-ts_(end-1);
            end
            optionsMMFVM.MaxStep = ts_(end)-ts_(1);
           
            % --- Shift and reset states of the MMFVM model
            % x = [Twi; Tw; Two; Tb; Q; p]        
            k = floor(t_sim/obj.Ts_u) + 1;
            x_ = [Twi(k);...                  Twi
                x_(1:obj.idx_Tw(end)-1);... Tw - shift water elements
                x_(obj.idx_Tw(end));...     Two - last Tw element
                x_(obj.idx_Tb);...          Tb - no change
                0;...                       Q = 0 - reset
                0];%                        p = 0 - reset
        else            
            % x = [Twi; Tw; Two; Tb; Q; p]        
            x_(obj.idx_Q) = 0;
            obj.x = x_;

        end
        
        % --- Prepare new run
        tstart = t_sim;
    end
   
end

if obj.isDisplayProgress
    t_elapsed = toc;
    fprintf('\tDone. Time elapsed %.2f s.\n',t_elapsed);
end

obj.saved.x.addArray(xs(:,1:nSaved)');
obj.saved.y.addArray(ys(:,1:nSaved)');
obj.saved.t.addArray(ts(:,1:nSaved)'); 
% obj.saved.isMMFVM.addArray

% --- Output
Q = obj.saved.y.get(1);
Two = obj.saved.y.get(2);
t_y = obj.saved.t.get;

    function [val,isterminal,dir] = zeroEvents(t, x)
        
        isterminal = [1;0];
        dir = [1;1];         
        % --- Leave to FVM when m_dot < m_dot_min
        k = floor(t/obj.Ts_u) + 1;
        Ts_c = obj.mw_e/m_dot(k);
        if Ts_c > obj.Ts_y_max
            val = [0;0];
            obj.switch_to_FVM = 1;
        else            
            % --- Catch the exact time when p = 1
            p = x(obj.idx_p);
            val = [p-1; % Stop at p = 1
                   p-0.5]; % Pass by 0.5
        end
    end

    function status = outFcn(t,x,flag)
        status = 0;
        isSave = 0;
        switch flag
            case 'init'
                if t(1) == 0
                    isSave = 1;
                    t = t(1);
                    p = x(obj.idx_p,:);
                end
            case []
                p = x(obj.idx_p,:);
                switch obj.sampling_mode
                    case 'equidistant'
                        isSave = abs(p - 1) > AbsTol;
                    case 'output'
                        isSave = abs(p - 1) < AbsTol;
                end
            case 'done'
        end
        % --- Save state
        if any(isSave) || obj.iteratingFVM
            n = numel(t);
            for i = 1:n
                if obj.iteratingFVM
                    x_interp = x;
                    
%                     [~, iV_dot, ~, iTai] = deal(1,2,3,4); % Input indexes
%                     k = floor(t(i)/obj.Ts_u) + 1;                    
%                     Q = 1/Nb_ * obj.UAba_model(u(k,iV_dot)*3600).*sum((x(obj.idx_Tb,i)-repmat(u(k,iTai),Nb_,1)), 1);
                    Q = x_interp(obj.idx_Q,i);
                    
                    W = obj.FVMmodel.x(obj.FVMmodel.idx_Tw(end-2));
                    P = obj.FVMmodel.x(obj.FVMmodel.idx_Tw(end-1));
                    E = obj.FVMmodel.x(obj.FVMmodel.idx_Tw(end));                    
                    
                    phi = hxModels.hxModelCommon.flux_limiter([-1,0,1],[W,P,E],11);                    
                    
                    %FOuw
                    m0 = [0, 0, 1];
                    %SOuw
                    m1 = 1/2*[0, -1, +3];
                    %TOuw
                    m2 = 1/8*[3, -10, +15];                    
                    
                    Two = m0*[W;P;E];
%                     Two = (m0 + phi(2)*(m1-m0))*[W;P;E];

%                     Two_switch = x_interp(obj.idx_Two);
%                     Two_ = min(Two,Two_switch);
                    y = [Q;Two];                
                elseif strcmp(obj.sampling_mode,'equidistant')
                    % --- Grid value interpolation
                    
                    if mod(t(i),obj.Ts_y) ~= 0 && ~obj.iteratingFVM %abs(p(i) - 1) < AbsTol
                        % Skip output when p=1
                        continue;
                    end
                    
                    [Tw,Two] = hxModels.MMFVM.grid_value_interpolation(...
                        x(1:obj.idx_Two,i),p(i),obj.iInterpolation_scheme,...
                        obj.z,t(i));
                    
                    % x = [Twi; Tw; Two; Tb; Q; p]
                    x_interp = [NaN; Tw; Two; x(obj.idx_Tb(1):end,i)];
                    [~, iV_dot, ~, iTai] = deal(1,2,3,4); % Input indexes
                    k = floor(t(i)/obj.Ts_u) + 1;                    
                    Q = 1/Nb_ * obj.UAba_model(u(k,iV_dot)*3600).*sum((x(obj.idx_Tb,i)-repmat(u(k,iTai),Nb_,1)), 1);
                    
                    % Q = x(obj.idx_Q,i)/p;
                    y = [Q;Two];
                    
                else % 'Ts_y'
%                     % [Twi; Tw; Two; Tb; Q; p]
                    if t(i) ~= 0
                        x_interp = [NaN; ...
                            x(1:obj.idx_Tw(end-1),i);... shifted state
                            x(obj.idx_Tw(end),i); ...
                            x(obj.idx_Tb(1):end,i)];
                    else
                        x_interp = x;
                    end
                    [im_dot, iV_dot, iTwi, iTai] = deal(1,2,3,4); % Input indexes
%                     k = floor(t(i)/obj.Ts_u) + 1;
%                     Q = 1/Nb_ * obj.UAba_model(u(k,iV_dot)*3600).*sum((x(obj.idx_Tb,i)-repmat(u(k,iTai),Nb_,1)), 1);                    
                    Q = x_interp(obj.idx_Q,i);
                    
                    W = x_interp(obj.idx_Tw(end-1)); %states already shifted
                    P = x_interp(obj.idx_Tw(end));
                    E = x_interp(obj.idx_Two);
                    

                    
                    %FO                    
                    m0 = [0 0 1];
                    
%                     phi = hxModels.hxModelCommon.flux_limiter([-1,0,1],[W,P,E],1);                    
%                     %SO
%                     m1 = 1/2*[0 1 1];
%                     %TO
%                     m2 = 1/8*[-1 6 3];
                    
                    Two = m0*[W;P;E]; %Two value only valid at complete ejection.
%                     Two = m1*[W;P;E];
%                     Two = (m0 + phi(2)*(m1-m0))*[W;P;E]; %creates bias in
    %                     steady state for systems with decay.
%                     Two = (m1 + phi(2)*(m2-m1))*[W;P;E];
                    
%                       % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%                       % HACK HACK HACK HACK (assumes constant velocity between consecutive MMFVM samplings)
%                       step_number = nSaved; 
%                       if step_number == 0
%                       elseif step_number <= Nw_  
% %                           tau_buffer(2*step_number+[-1:0]) = [t_half - tstart; t(i)- t_half];
%                           t(i) = t_half;
%                       else % step_number > Nw_ + 1  
%                           
% %                           last_half = t(i) - t_half;
% %                           t(i) = t(i) - (sum([tau_buffer;last_half]) - sum(tau_buffer(2:end)));                        
% %                           % leave first add last
% %                           tau_buffer = [tau_buffer(2:end); last_half];                         
%                           t(i) = t(i) - tau;
%                       end
%                         t(i) = t_half;                    

                    y = [Q;Two];
                end
                
                obj.x = x_interp;
                
                % --- Saving
                nSaved = nSaved + 1;
                if nSaved >= alloc_size
                    if obj.isStateSaving
                        xs = [xs zeros(obj.nx, alloc_block)]; 
                    end
                    ys = [ys zeros(obj.ny, alloc_block)];
                    ts = [ts zeros(1, alloc_block)];
                    alloc_size = alloc_size + alloc_block;
                end
                    if obj.isStateSaving
                        xs(:,nSaved) = x_interp;
                    end
                ys(:,nSaved) = y;
                ts(nSaved) = t(i);
            end
        end
    end

end

