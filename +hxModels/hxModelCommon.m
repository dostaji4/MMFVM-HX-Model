classdef hxModelCommon < handle
    % Abstract class for various heat exchanger simulation methods
    
    properties (Abstract, Constant)
        model_type        
    end
    
    properties (Abstract)
        simulation_mode
    end
    
    properties
        Nw
        Nwb
        isStateSaving = 1;
        isSaving = 1;
        isDisplayProgress = 1;
        isPlotSteadyState = 1;
        Tw0
        Tb0
        UAba_model
        UAwb_model
        Cb
        mw
        Ts_y_max = 5; %[s]
        Ts_y         % Prefered sampling of state and output
        Ts_u         % Input sampling
    end   
    
    properties(Hidden)
        plot_water_state_method = 'average'; %'middle', 'average'
    end
    
    properties (Access = protected)
        saved = struct('t',0,...
            'y',0,...
            'x',0,...
            'tss',0,...
            'yss',0,...
            'xss',0);  % Structure of saveTimeSnapshot data for saving
        Nb
        Cw
        mw_e
        isInitialized = 0;        
        nx           % State dimension
        x            % State
        ny           % Output dimension
        idx_Tw       % Indexes of water elements in the state vector
        idx_Tb       % Indexes of body element in the state vector
        z            % State grid positions
        dz           % Distance between grid nodes
    end
    
    properties (Constant, Access = protected)
        cp_w = 4180;        % [J/kg/K]
        cv_std_air = 1195;  % [J/m3/K] @ 25C/50%/101.235
        m_dot_max = 300/3600;    % Maximal water flow after which it gets clamped [kg/s]
        
        k_w = 0.5980; %[W/(m*K)] thermal conductivity of water
        rho_w = 998.1618; %[kg/m3] density
        nu_w = 1.0035e-06; %[m2/s] kinematic viscosity
        Pr_w = 7.0092; % Prandlt number of water        
    end
    
    properties (Access = {?hxModel, ?hxModels.hxModelCommon}, Constant)        
        simulation_modeList = {'discrete','continuous'};
    end
    
    methods (Access = {?hxModel,?hxModels.hxModelCommon})
        [y, t] = step(obj,t, u, Ts, varargin)
    end
    
    methods (Abstract, Access = protected)
        initialize_hx(obj)
    end
    
    methods
        % Default run function for all models
        function [Q, Two, t_y] = run(obj, m_dot, V_dot, Twi, Tai, Ts_u)
            
            % --- Input checks and routines
            narginchk(6,6);
            %             assert(nargin == 6, 'Not enough input arguments.');
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
            u = cell2mat(tmp);
            obj.Ts_u = Ts_u;
            
            if isempty(obj.Ts_y)
                obj.Ts_y = obj.Ts_y_max;
            end
            
            % --- Initialization
            obj.initialize;
            t_sim = 0;
            t_end = Ts_u * (nData - 1);
            alloc_block = 1000;
            xs = zeros(obj.nx,alloc_block);
            ys = zeros(obj.ny,alloc_block);
            ts = zeros(1, alloc_block);
            alloc_size = alloc_block; % actual allocated size
            k = 1;
            
            if obj.isDisplayProgress
                tic
                fprintf('## Started simulation of "%s" model \n',obj.model_type);
            end
            
            % --- Run the simulation --------------------------------------
            while t_sim <= t_end
                
                % --- Index of the current input
                iInput = floor(t_sim/obj.Ts_u) + 1;
                
                % --- Save state x(k)
                if obj.isStateSaving
                    %                     obj.saved.x.add(obj.x);
                    xs(:,k) = obj.x;
                end
                ts(:,k) = t_sim;
                
                % --- Step model
                % Performace note: invoke methods by f(obj, ..) rather then obj.f(..)
                [y, t_sim] = step(obj,t_sim, u(iInput,:),obj.Ts_y);
                
                % --- Save output y(k) and t(k)
                ys(:,k) = y;
                k = k+1;
                
                % Allocate more space for saving
                if k == alloc_size
                    xs = [xs zeros(obj.nx, alloc_block)]; %#ok<AGROW>
                    ys = [ys zeros(obj.ny, alloc_block)]; %#ok<AGROW>
                    ts = [ts zeros(1, alloc_block)];      %#ok<AGROW>
                    alloc_size = alloc_size + alloc_block;
                end
                
                % IMPLEMENT PROGRESS DISPLAY HERE
            end
            
            if obj.isDisplayProgress
                t_elapsed = toc;
                fprintf('\tDone. Time elapsed %.2f s.\n',t_elapsed);
            end
            
            obj.saved.x.addArray(xs(:,1:k-1)');
            obj.saved.y.addArray(ys(:,1:k-1)');
            obj.saved.t.addArray(ts(:,1:k-1)');
            % --- Output
            Q = obj.saved.y.get(1);
            Two = obj.saved.y.get(2);
            t_y = obj.saved.t.get;
        end
        
        function [Qss, Two_ss, t_y, Twss, Tbss] = run_ss(obj, m_dot, V_dot, Twi, Tai, Ts_u)
            
            % --- Input checks and routines
            narginchk(6,6);
            %             assert(nargin == 6, 'Not enough input arguments.');
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
            u = cell2mat(tmp);
            obj.Ts_u = Ts_u;
            
            if isempty(obj.Ts_y)
                obj.Ts_y = obj.Ts_y_max;
            end
            
            % Initialization
            obj.initialize;
            t_sim = 0;
            t_end = Ts_u * (nData-1);
            alloc_block = 1000;
            nxss = obj.Nw + obj.Nb;
            xs = zeros(nxss,alloc_block);
            ys = zeros(obj.ny,alloc_block);
            ts = zeros(1, alloc_block);
            alloc_size = alloc_block; % actual allocated size
            k = 1;
            
            % --- Run the ss simulation -----------------------------------
            while t_sim <= t_end
                
                % --- Index of the current input
                iInput = floor(t_sim/obj.Ts_u) + 1;
                
                % --- Step model
                % Performace note: invoke methods by f(obj, ..) rather then obj.f(..)
                [yss, xss] = ss_pde(obj,u(iInput,:));
                
                % --- Save
                xs(:,k) = xss;
                ys(:,k) = yss;
                ts(:,k) = t_sim;
                k = k+1;
                % Allocate more space for saving
                if k == alloc_size
                    xs = [xs zeros(nxss, alloc_block)]; %#ok<AGROW>
                    ys = [ys zeros(obj.ny, alloc_block)]; %#ok<AGROW>
                    ts = [ts zeros(1, alloc_block)];      %#ok<AGROW>
                    alloc_size = alloc_size + alloc_block;
                end
                
                % --- Advance time
                t_sim = t_sim + obj.Ts_y;
            end
            
            obj.saved.xss.addArray(xs(:,1:k-1)');
            obj.saved.yss.addArray(ys(:,1:k-1)');
            obj.saved.tss.addArray(ts(:,1:k-1)');
            
            % --- Output
            Qss = obj.saved.yss.get(1);
            Two_ss = obj.saved.yss.get(2);
            t_y = obj.saved.tss.get;
            Twss = obj.saved.xss.get(1:obj.Nw);
            Tbss = obj.saved.xss.get(obj.Nw+1:obj.Nw+obj.Nb);
        end
        
        function h = plot_outputs(obj,varargin)
            t = obj.saved.t.get;
            assert(~isempty(t), 'There are no data to plot. Run the "run" method first.');
            ax = zeros(obj.ny,1);
            tit = {'Heat output','Water outlet'};
            label = {'Heat flow [W]','Temperature [$^\circ$C]'};
            for i = 1:obj.ny
                ax(i) = subplot(obj.ny,1,i);
                h(i) = plot(t,obj.saved.y.get(i),varargin{:}); hold on
                title(tit{i},'FontWeight','normal','Interpreter','Latex');
                ylabel(label{i},'Interpreter','Latex');
            end
            xlabel('Time [s]','Interpreter','latex');
            linkaxes(ax,'x');
            h = h(1);
        end
        
        function plot_ss_outputs(obj, varargin)
            tss = obj.saved.tss.get;
            assert(~isempty(tss), 'There are no data to plot. Run the "run_ss" method first.');
            ax = zeros(obj.ny,1);
            tit = {'Heat [W]','Two [degC]'};
            for i = 1:obj.ny
                ax(i) = subplot(obj.ny,1,i);
                plot(tss,obj.saved.yss.get(i),varargin{:}); hold on
                title(tit{i});
            end
            linkaxes(ax,'x');
        end
        
        function h = plot_states(obj,varargin)
            t = obj.saved.t.get;
            assert(~isempty(t), 'There are no data to plot. Run the "run" method first.');
            
            Tw = obj.saved.x.get(obj.idx_Tw);
            Tb = obj.saved.x.get(obj.idx_Tb);
            
            n_plot = 4; % best is to have  Nw = (2k+1)*n_plot, k from natural numbers
                       
%             --- Water elements - get interpolated spatial temperatures
%             if strcmp(obj.model_type,'MMFVM') && strcmp(obj.plot_water_state_method,'middle')
%                 x_plot_mid = (obj.Nw+1/2)/obj.Nw * obj.equidistant_grid(n_plot);
%                 x_mid = (obj.Nw+1)/obj.Nw * obj.equidistant_grid(obj.Nw+1);
%                 Two = obj.saved.x.get(obj.idx_Tw(end)+1);
%                 Tw_plot = obj.refactor([Tw,Two],x_mid,x_plot_mid,obj.plot_water_state_method);
%             else
                x_plot_mid = obj.equidistant_grid(n_plot);% middle points of plotting positions
                x_mid = obj.equidistant_grid(obj.Nw);
                Tw_plot = obj.refactor(Tw,x_mid,x_plot_mid,obj.plot_water_state_method);
%             end
            
            
            % --- Body elements - get interpolated spatial temperatures
            x_plot_mid = obj.equidistant_grid(n_plot);% middle points of plotting positions
            x_mid = obj.equidistant_grid(obj.Nb);
            Tb_plot = obj.refactor(Tb,x_mid,x_plot_mid,'average');
            
            ax(1) = subplot(2,1,1);
            h = plot(t,Tb_plot,varargin{:}); hold on
            title('Body volumes','FontWeight','normal','Interpreter','Latex');
            ylabel('Temperature [$^\circ$C]','Interpreter','Latex');
            ax(2) = subplot(2,1,2);
            plot(t,Tw_plot,varargin{:}); hold on
            title('Water volumes','FontWeight','normal','Interpreter','Latex');
            ylabel('Temperature [$^\circ$C]','Interpreter','Latex');
            xlabel('Time [s]','Interpreter','latex');
            linkaxes(ax,'x');
            h = h(1);
        end
        
        function plot_ss_states(obj, varargin)
            tss = obj.saved.tss.get;
            assert(~isempty(tss), 'There are no data to plot. Run the "run_ss" method first.');
            
            Tw = obj.saved.xss.get(1:obj.Nw);
            Tb = obj.saved.xss.get(obj.Nw+1:obj.Nw+obj.Nb);
            
            n_plot = 5; % best is to have  Nw = (2k+1)*n_plot, k from natural numbers
            
            x_plot_mid = obj.equidistant_grid(n_plot);% middle points of plotting positions
            
            % --- Water elements - get interpolated spatial temperatures
            x_mid = obj.equidistant_grid(obj.Nw);
            Tw_plot = obj.refactor(Tw,x_mid,x_plot_mid,obj.plot_water_state_method);
            
            % --- Body elements - get interpolated spatial temperatures
            x_mid = obj.equidistant_grid(obj.Nb);
            Tb_plot = obj.refactor(Tb,x_mid,x_plot_mid);
            
            ax(1) = subplot(2,1,1);
            plot(tss,Tb_plot,varargin{:}); hold on
            title('Body elements [degC]');
            ax(2) = subplot(2,1,2);
            plot(tss,Tw_plot,varargin{:}); hold on
            title('Water elements [degC]');
            linkaxes(ax,'x');
        end
        
        % Default initialize function
        function initialize(obj)
            
            % --- Assert common settings
            props = {'Nw','Nwb','mw','Cb','Tw0','Tb0','UAba_model','UAwb_model'};
            for i = 1:numel(props)
                assert(~isempty(obj.(props{i})),'The property "%s" has not been set.',props{i});
            end
            assert(mod(obj.Nw,obj.Nwb) == 0, 'The ratio Nw/Nwb must be an integer.');
            obj.Nb = obj.Nw / obj.Nwb;
            
            % --- Initialize HX
            obj.initialize_hx;
            
            % --- Saving
            obj.initialize_saving;
            
            obj.isInitialized = 1;
        end
    end
    
    methods(Access = protected)
        
        function initialize_saving(obj)
            assert(~isempty(obj.nx) && ~isempty(obj.ny),'Initialize state and output sizes, nx and ny, first.');
            
            % Allocate time saving
            obj.saved.t = hxModels.saveTimeSnapshot(1);
            obj.saved.tss = hxModels.saveTimeSnapshot(1);
            
            % Allocate output saving
            obj.saved.y = hxModels.saveTimeSnapshot(obj.ny);
            obj.saved.yss = hxModels.saveTimeSnapshot(obj.ny);
            
            % Allocate state saving
            if obj.isStateSaving
                obj.saved.x = hxModels.saveTimeSnapshot(obj.nx);
                nxss = obj.Nw + obj.Nb;
                obj.saved.xss = hxModels.saveTimeSnapshot(nxss);
            end
        end
        
        % Default steady state calculation
        function [yss, xss] = ss(obj,uss)
            [yss, xss] = ss_pde(obj,uss);
        end
        % Calculate exact heat exchanger steady state profile defined by PDE
        function [yss, xss, varargout] = ss_pde(obj,uss)
            
            % --- Input checking and initialization -----------------------
            assert(numel(uss) == 4,'Input dimension must be 4 ([u(m_dot),u(V_dot),u(Twi),u(Tai)]).');
            [im_dot, iV_dot, iTwi, iTai] = deal(1,2,3,4); % Input indexes
            
            UAwb = obj.UAwb_model(uss(im_dot)*3600);
            UAba = obj.UAba_model(uss(iV_dot)*3600);
            
            %coefficients
            a = UAwb/(uss(im_dot)*obj.cp_w);
            c = UAba/UAwb;
            k = -a*c/(c+1);
            ek = exp(k);
            
            % analytic ss temperature profile
            Tbss_f = @(x) uss(iTai) + 1/(c+1)*(uss(iTwi)-uss(iTai))*ek.^x;     %ok
            Twss_f = @(x) uss(iTai) +         (uss(iTwi)-uss(iTai))*ek.^x;     %ok
            
            xw_edge = linspace(0,1,obj.Nw+1); %there is N+1 edges
            xw_step = xw_edge(2) - xw_edge(1);
            xw_mid = xw_edge(1:end-1) + xw_step/2;
            
            xb_edge = linspace(0,1,obj.Nb+1); %there is N+1 edges
            xb_step = xb_edge(2) - xb_edge(1);
            xb_mid = xb_edge(1:end-1) + xb_step/2;
            
            % analytic element ss temperature (integral over element)
            s = xw_step;  % spatial step
            Tbe_ai_f = @(x0) uss(iTai) + (uss(iTwi) - uss(iTai))/(s*k* (c+1) ) *  ...
                ek.^(x0) * (ek^s - 1);                 %ok
            Twe_ai_f = @(x0) uss(iTai) + (uss(iTwi) - uss(iTai))/(s*k        ) *  ...
                ek.^(x0) * (ek^s - 1);                 %ok
            
            if nargout > 0
                % approximate the integral over element by midpoint function
                % value
                %                 Tbss = Tbss_f(xb_mid)';
                %                 Twss = Twss_f(xw_mid)';
                
                % true value of integral over element
                Tbss = Tbe_ai_f(xb_edge(1:end-1))';
                Twss = Twe_ai_f(xw_edge(1:end-1))';
                
                % b vectors needed for formulation of B matrix
                idx = (0:1:obj.Nw-1)';
                s = xw_step;
                bw = ek.^(idx*s)*(ek^s - 1)/(s*k);
                
                idx = (0:1:obj.Nb-1)';
                s = xb_step;  % spatial step
                bb = ek.^(idx*s)*(ek^s - 1)/(s*k);
                
                Twoss = Twss_f(1);
                % analytic heat output
                Qss = UAba*(uss(iTwi)-uss(iTai))*1/(c+1)*1/k*(ek-1);
                yss = [Qss,Twoss];
                xss = [Twss;Tbss];
                
                if nargout > 2
                    varargout = {bw, bb, c, ek};
                end
            else
                fprintf('## Steady state:\n');
                hold on
                
                % true profile
                n_plot = 100;
                x_pl = linspace(0,1,n_plot);
                p = plot(...
                    x_pl, Tbss_f(x_pl),...
                    x_pl, Twss_f(x_pl),...
                    'LineWidth',1);
                p(1).DisplayName = 'B pde';
                p(2).DisplayName = 'W pde';
                
                % midpoint approximation
                Tw_m = Twss_f(xw_mid);
                Tb_m = Tbss_f(xw_mid);
                p = plot(...
                    xw_mid,Tb_m,'+',...
                    xw_mid,Tw_m,'+',...
                    'LineWidth',1);
                p(1).DisplayName = 'B midpoint approx.';
                p(2).DisplayName = 'W midpoint approx.';
                
                fprintf('\tNorm_1(True-Midpoint approximation) = %e degC\n',...
                    norm([Tw_m;Tb_m] - ...
                    [Twe_ai_f(xw_edge(1:end-1));Tbe_ai_f(xw_edge(1:end-1))],1));
                
                % analytical integral
                p = plot(...
                    xw_mid,Tbe_ai_f(xw_edge(1:end-1)),'o',...
                    xw_mid,Twe_ai_f(xw_edge(1:end-1)),'o',...
                    'LineWidth',1);
                p(1).DisplayName = 'B a. integral';
                p(2).DisplayName = 'W a. integral';
                
                % mean values
                for i = 1:obj.Nw
                    Twe_ni(i) = integral(@(x) 1/xw_step*Twss_f(x),xw_edge(i),xw_edge(i+1),'RelTol',0,'AbsTol',1e-12);
                    Tbe_ni(i) = integral(@(x) 1/xw_step*Tbss_f(x),xw_edge(i),xw_edge(i+1),'RelTol',0,'AbsTol',1e-12);
                end
                p = plot(...
                    xw_mid,Tbe_ni,'x',...
                    xw_mid,Twe_ni,'x',...
                    'LineWidth',1);
                p(1).DisplayName = 'B n. integral';
                p(2).DisplayName = 'W n. integral';
                
                fprintf('\tNorm_1(True-Numerical integral) = %e degC\n',...
                    norm([Twe_ni;Tbe_ni] - ...
                    [Twe_ai_f(xw_edge(1:end-1));Tbe_ai_f(xw_edge(1:end-1))],1));
                
                ylabel('Temperature [degC]');
                xlabel('Spatial position x [0-1]');
                legend('off');
                legend('show');
            end
        end
        
        function copy(obj, srcObj)
            mc = ?hxModels.hxModelCommon;
            for i = 1:numel(mc.PropertyList)
                p = mc.PropertyList(i);
                if  any(strcmpi(p.GetAccess,{'public','protected'})) && ...
                        p.Abstract == 0 && ...
                        p.Constant == 0
                    obj.(p.Name) = srcObj.(p.Name);
                end
            end
        end
        
    end
    
    methods(Static, Access = protected)
        function phi = flux_limiter(z,y,iLimiter)
            
            %...     z          independent variable (input)
            %...
            %...     y          dependent variable (input)
            %...
            %...     iLimiter   limiter function number (see the code)
            
            
            delta=1.0e-05; % defines gradient under which first order approximation is always used
            n = numel(z);
            assert(numel(y) == n);
            r = NaN(n,1);
            phi = NaN(n,1);
            
            % --- Flux limiter calculation for all elements
            % Flux ~ first order approximation (phi=0) or in smooth regions (phi>0) second order
            %        approximation
            for i = 1:n
                
                if i == 1 || i == n
                    phi(i) = 0; % default low order scheme at boundary
                    continue;
                end
                
                % --- Computation of ratio of succesive gradients r(i) and
                %     a slope limiter function phi(i) for an elements' East face
                
                if abs(y(i)-y(i-1)) < delta
                    % Gradient too shallow, use low order scheme for better numerical
                    % accuracy
                    phi(i) = 0;
                else
                    % Ratio of succesive gradients
                    r(i)=((y(i+1)-y(i)) / (z(i+1)-z(i))) / ...
                        ((y(i)-y(i-1)) / (z(i)-z(i-1)));
                    
                    % --- Limiter function
                    %             {'Van_Leer', 'superbee', 'smart', 'mc', 'minmod', 'Koren', 'CHARM', 'Osher', 'UMIST','ospre','HCUS','HQUICK','van_Albada_1','van_Albada_1','mg','Sweby','fo'};
                    % Wouver: Simulating PDE in Matlab/Octave
                    % https://en.wikipedia.org/wiki/Flux_limiter
                    switch iLimiter
                        case 1 % Van_Leer
                            % Van Leer limiter
                            phi(i) = ( r(i)+abs(r(i))) / (1+r(i));
                        case 2 % superbee
                            phi(i) = max([0 min(2*r(i),1) min(r(i),2)]);
                        case 3 % smart
                            phi(i) = max([0, min([2*r(i) .75*r(i)+.25 4])]); %wiki
                            %                     phi(i)=max(0, min([4*r(i) .75*r(i)+.25 2])); %wouver
                        case 4 % mc (Monotized Central)
                            phi(i)=max([0, min([2*r(i) (1+r(i))/2 2])]);
                        case 5 % minmod
                            phi(i) = max([0, min([1 r(i)])]);
                        case 6 % Koren
                            phi(i) = max([0, min([2*r(i) (1+2*r(i))/3 2])]);
                        case 7 % CHARM
                            phi(i) = (r(i) > 0)* (r(i)*(3*r(i)+1))/(r(i)+1)^2;
                        case 8 % Osher
                            beta = 1.5; % beta in [1, 2]
                            phi(i) = max([0, min(r(i),beta)]);
                        case 9 % UMIST
                            phi(i) = max([0,min([2*r(i), (.25 + .75*r(i)), (.75 + .25*r(i)), 2])]);
                        case 10 % ospre
                            phi(i) = 1.5*(r(i)^2 + r(i))/(r(i)^2 + r(i) + 1);
                        case 11 % HCUS
                            phi(i) = 1.5*(r(i) + abs(r(i))) / (r(i) + 2);
                        case 12 % HQUICK
                            phi(i) = 2*(r(i)+abs(r(i)))/(r(i) + 3);
                        case 13 % van_Albada_1
                            phi(i) = r(i)*(r(i)+1) / (r(i)^2 + 1);
                        case 14 % van_Albada_2
                            phi(i) = 2*r(i) / (r(i)^2 + 1);
                        case 15 % mg (Minmod Generalized)
                            beta = 1.5; % beta in [1,2]
                            phi(i) = max([0, min([beta*r(i), (1+r(i))/2, beta])]);
                        case 16 % Sweby
                            beta = 1.5; %beta in [1 2]
                            phi(i) = max([0,min([beta*r(i),1]),min([r(i) beta])]);
                        case 17 %fo (First order)
                            phi(i) = 0;
                        otherwise
                            error('Unknown limiter "%d".',iLimiter);
                    end
                end
            end
        end
        
        function [z,dz] = equidistant_grid(n)
            dz = 1/n;
            z = (1:n)'*dz - dz/2;
        end
        
        function ynew = refactor(yold,zold,znew,method)
            % method = 'middle' or 'average'
            nold = numel(zold);
            nnew = numel(znew);
            if nargin < 4
                method = 'average';
            end
            dzo = diff(zold);
            assert(all((dzo(1) - dzo)<eps));
            dzo = dzo(1);
            dzn = diff(znew);
            assert(all((dzn(1) - dzn)<eps));
            dzn = dzn(1);
            assert(size(yold,2) == nold,'Make sure the value matrix "yold" has the same number of collumns as there is elements in "zold".');
            if nold >= nnew
                done = 0;
                
                % !!! ASSUMES THE SAME Z_MAX FOR ZOLD AND ZNEW !!!
                % Select or interpolate in between
                % Is the right number of elements?  - simply choose the
                % right indexes (Nw != (2k+1)*n_plot => k=Nw/10-0.5, mod(k)?=0
                zold_end = zold(end) + dzo/2;
                znew_end = znew(end) + dzn/2;
                k = nold/(2*nnew)-0.5;
                if mod(k,1) == 0 && (zold_end - znew_end) < eps 
                    switch method
                        case 'middle' % Select middle position
                            idx = (1:nnew)*(2*k+1) - floor(k+0.5);
                            if all((zold(idx) - znew) < eps)
                                ynew = yold(:,idx);
                                done = 1;
                            end
                        case 'average' % Average over big cell
                            ynew = squeeze(mean(reshape(yold,size(yold,1),2*k+1,[]),2));
                            if size(yold,1) ~= size(ynew,1)
                                ynew = ynew';
                            end
                            done = 1;
                        otherwise
                            error('Unknown method "%s".',method);
                    end
                end
                
                if ~done                    
                    % otherwise interpolate                    
                    switch method 
                        case 'middle'
                            nr = size(yold,1);
                            ynew = zeros(nr,nnew);
                            for i = 1:nr
                                ynew(i,:) = interp1(zold,yold(i,:),znew,'pchip','extrap');
                            end
                        case 'average'
                            nr = size(yold,1);
                            ynew = zeros(nr,nnew);                            
                            for i = 1:nr
                                f = @(p) interp1(zold,yold(i,:),znew+p,'pchip','extrap')';                                
                                ynew(i,:) = 1/dzn * integral(f,-dzn/2,+dzn/2,'ArrayValued',true);
                            end
                            
                        otherwise
                            error('Unknown method "%s".',method);
                    end
                end
            else
                % Approximate profile, extrapolate
                
            end
            
        end
    end
    
    methods(Static)
        function out = get_all_model_types
            s = what('hxModels');
            out = s.classes;
        end
    end
    
end
