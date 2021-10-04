classdef saveTimeSnapshot < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here    
    properties        
        nChannels = 1;
        nDataSaved = 0;
    end
    
    properties (Access = private)
       allocationBlock = 1000;         
       array;
    end
    
    methods
        function obj = saveTimeSnapshot(varargin)        
            switch nargin
                case 0
                    nChannels = 1;
                case 1
                    nChannels = varargin{1};
                otherwise
                    error('Too many input arguments.');
            end
            obj.array = zeros(obj.allocationBlock,nChannels);            
            obj.nChannels = nChannels;
            obj.nDataSaved = 0;
        end
        
        function add(obj, value)
            assert(numel(value) == obj.nChannels,...
                'Wrong dimension. The data has %d elements, but the object was initiated for %d channels.',numel(value),obj.nChannels);            
            if obj.nDataSaved > size(obj.array,1)
                obj.array = [obj.array; zeros(obj.allocationBlock,obj.nChannels)];                
            end
            obj.nDataSaved = obj.nDataSaved + 1;
            obj.array(obj.nDataSaved,:) = value(:)';                        
        end
        
        function addArray(obj, arrayIn)            
            s = size(arrayIn);
            assert(s(2) == obj.nChannels,'Insert an array, where number of columns equals the number of channels (=%d)',obj.nChannels); %
            nn = obj.nDataSaved + s(1);
            no = size(obj.array,1);
            if nn > no
               m = ceil((nn-no)/obj.allocationBlock);
               obj.array = [obj.array; zeros(m * obj.allocationBlock,obj.nChannels)];                
            end
            obj.array(obj.nDataSaved+1 : obj.nDataSaved + s(1),:) = arrayIn;
            obj.nDataSaved = obj.nDataSaved + s(1);
        end
        
        function out = get(obj,varargin)            
            switch nargin
                case 1
                    out = obj.array(1:obj.nDataSaved,:);                        
                case 2 
                    % Give channel
                    assert(all(varargin{1} <= obj.nChannels),'Index out of channel number.')
                    out = obj.array(1:obj.nDataSaved,varargin{1});                        
                otherwise
                    error('Too many input arguments.');
            end
        end

% THIS IS VERY SLOW!!
%         function varargout = subsref(obj,S)            
%                 switch S(1).type
%                     case '()'       
%                         % Return array
%                         array_ = obj.get;
%                         [varargout{1:nargout}] = builtin('subsref',array_,S);
%                         return
%                     case '.'
%                         if strcmp(S(1).type,'.') && ...
%                             strcmp(S(1).subs,'get') && ...
%                             numel(S) > 1                            
%                                 % Return array
%                                 array_ = obj.get;
%                                 [varargout{1:nargout}] = builtin('subsref',array_,S(2));
%                             return
%                         end
%                         
%                         % call regular object methods
%                         [varargout{1:nargout}] = builtin('subsref',obj,S);
%                         return
%                     otherwise
%                         error('Indexing "%s" not allowed. Use "()" indexing.',S.type);
%                 end            
%         end
        
        function n = numArgumentsFromSubscript(obj, s, ic)
            n = builtin('numArgumentsFromSubscript', obj, s, ic);
        end

        function plot(obj)
            plot(repmat([1:obj.nDataSaved]',1,obj.nChannels),obj.array(1:obj.nDataSaved,:));
        end
    end
end

