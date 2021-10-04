classdef MMFVM < hxModels.hxModelCommon
    
    properties (Constant)
        model_type = 'MMFVM';
    end
    
    properties
        simulation_mode = 'discrete'; %'continuous';
        sampling_mode = 'output'; %'equidistant';
        interpolation_scheme = 'quadratic'; 
    end
    
    properties (Access = {?hxModel, ?hxModelCommon}, Constant)        
        sampling_modeList = {'equidistant','output'};
        interpolation_schemeList = {'constant','linear','quadratic','TVD'};
    end
    
    properties (Access = private)
        u_sum = zeros(1,4);
        t_y_last_hit = 0;
        y_persist
        u_persist = [Inf;Inf;Inf;Inf];
        Ad_persist = NaN;
        B_persist = NaN;
        UAba_persist = NaN;        
        idx_Q   % index of Q in the state vector
        idx_Two % index of Two in the state vector
        idx_p
        
        iInterpolation_scheme
        FVMmodel % object for low flow simulations
        iteratingFVM
        switch_to_FVM = 0;
        FVM_Nw_multiplication_factor = 5; %1,5,15,25,etc.        
    end
    
    methods (Static, Access = private)
        [Ac1,Ac2,Bc,Cc,Dc,Aq,Bq] = ABCD_continuous(Ts_x,nl,nlb,Cw,Cb,UAwb,UAba);
    end
    
    methods (Static, Access = protected)
        [Tw, Two] = grid_value_interpolation(TwTwo,p,z,scheme,varargin);
    end
    
    methods (Access = {?hxModel,?hxModels.hxModelCommon})
        % y = [Q, Two]
        % u = [m_dot, V_dot, Twi, Tai]
        [y, t] = step(obj,t, u, Ts);
        dxdt = dxdt_rhs(obj,t,x,u,Ts_u,tau);
    end
    
    methods
        % y = [Q, Two]
        % u = [m_dot, V_dot, Twi, Tai]
        [Q, Two, t_y] = run(obj, m_dot, V_dot, Twi, Tai, Ts_u)
    end
    
    methods (Access = protected)
        [A, B, C, D] = ABCD_exact_discretization(obj, u, Ts_x);
        
        function initialize_hx(obj)
            
            % --- Method specific assertions
            props = {'Ts_y_max','simulation_mode','sampling_mode','interpolation_scheme'};
            hasMemberList = [0,1,1,1];
            for i = 1:numel(props)
                if hasMemberList(i)
                    % Check non-emptiness
                    assert(~isempty(obj.(props{i})),...
                        'Property "%s" not specified. Set the property with one of the following options: ''%s''.',...
                        props{i}, strjoin(obj.([props{i} 'List']),''', '''));
                    % Check validity
                    assert(ismember(obj.(props{i}),obj.([props{i} 'List'])),...
                        'Unknown setting "%s" for property "%s". Use one the the following options: ''%s''.',...
                        obj.(props{i}),props{i},strjoin(obj.([props{i} 'List']),''', '''));
                else
                    % Check non-emptiness
                    assert(~isempty(obj.(props{i})), 'Property "%s" not set. Please specify it first.',...
                        props{i});
                end
            end
            obj.iInterpolation_scheme = -1 + find(ismember(obj.interpolation_schemeList,obj.interpolation_scheme));
            
            % --- Initialize precomputed constants
            obj.Cw = obj.mw * obj.cp_w;  %[J/K] heat capacity of water inside HX
            obj.mw_e = obj.mw / obj.Nw;
%             if isempty(obj.Ts_y_max)
%                 obj.Ts_y_max = obj.mw_e / obj.m_dot_min; % Sampling interval used for energy flow model [s]
%             end
            
            % --- Initial water and body temperatures
            tmp = {obj.Tw0,obj.Tb0};
            nData = [obj.Nw, obj.Nb];
            for i = 1:numel(tmp)
                assert(~isempty(tmp{i}),'Initialize water and body temperatures first. Allowed values are: scalar for constant profile or a vector of dimension Nw=%d (Nb=%d).',obj.Nw,obj.Nb);
                if numel(tmp{i}) == 1
                    % Replicate nData times
                    tmp{i} = tmp{i} * ones(nData(i),1);
                elseif numel(tmp{i}) == nData(i)
                    % Ensure column vector
                    tmp{i} = tmp{i}(:);
                else
                    error('Incorrect dimension of input #%d "%s". Specify either a scalar value or a vector of length %d.',i,inputname(i),nData);
                end
            end
            [obj.Tw0,obj.Tb0] = tmp{:};
            
            % --- Initial outputs
            obj.y_persist = [0; obj.Tw0(end)];
            
            % --- State space model
            switch obj.simulation_mode
                case 'discrete'
                    % x = [Twi; Tw; Two; Tb; Q]
                    obj.x = [NaN; obj.Tw0; NaN; obj.Tb0; NaN]; % state properly initialized in the first step
                    obj.nx = 1 + obj.Nw + 1 + obj.Nb + 1; 
                    obj.ny = 2;
                    
                    obj.idx_Tw = 1 + (1:obj.Nw);                    
                    obj.idx_Tb = 1 + obj.Nw + 1 + (1:obj.Nb);                    
                    obj.idx_Two = 1 + obj.Nw + 1;
                    obj.idx_Q = 1 + obj.Nw + 1 + obj.Nb + 1;
                case 'continuous'
                    % x = [Twi; Tw; Two; Tb; Q; p]
                    obj.x = [NaN; obj.Tw0; NaN; obj.Tb0; NaN; 0]; % state properly initialized in the first step
                    obj.nx = 1 + obj.Nw + 1 + obj.Nb + 1 + 1; 
                    obj.ny = 2;
                    
                    obj.idx_Tw = 1 + (1:obj.Nw);
                    obj.idx_Tb = 1 + obj.Nw + 1 + (1:obj.Nb);
                    obj.idx_Two = 1 + obj.Nw + 1;
                    obj.idx_Q = 1 + obj.Nw + 1 + obj.Nb + 1;
                    obj.idx_p = 1 + obj.Nw + 1 + obj.Nb + 1 + 1;
            end
            
            % --- Define grid (grid without boundary points)
            [obj.z,obj.dz] = hxModels.hxModelCommon.equidistant_grid(obj.Nw);
            
            % --- Create object for FVM low flow simulations
            obj.FVMmodel = hxModels.FVM;
            obj.FVMmodel.copy(obj);             
            obj.FVMmodel.isSaving = 0; % disable double saving
            obj.FVMmodel.Nw = obj.FVM_Nw_multiplication_factor*obj.Nw;            
            obj.FVMmodel.Nwb = 1;
            obj.FVMmodel.Tw0 = obj.Tw0(1);
            obj.FVMmodel.Tb0 = obj.Tb0(1);
            obj.FVMmodel.initialize;
        end
        
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
            
            % Alocate mode saving
            obj.saved.isMMFVM = hxModels.saveTimeSnapshot(1);
        end
    end
    
end

