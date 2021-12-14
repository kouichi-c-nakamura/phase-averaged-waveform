classdef Chan
    % Chan class is an abstract class to store time series data
    % Data with an uniform Time vector.
    %
    % Written by Kouichi C. Nakamura Ph.D.
    % MRC Brain Network Dynamics Unit
    % University of Oxford
    % kouichi.c.nakamura@gmail.com
    % 15-Aug-2017 15:26:17
    %
    % See Also 
    % Record, MetaEventChan, WaveformChan, EventChan,
    % MarkerChan, ChanInfo
    
    properties (Dependent = true)
        ChanTitle
        DataUnit % The measurement unit for Data
        Start % Time at the beginning of the file in Seconds.
        SRate % Sampling Rate in Hertz.
        Header % empty or structure inherited from Spike2 exported MATLAB MAT file excluding the 'values' field. ChanComment$ is stored in Header.comment.
        Path % string for file path of the corresponding .mat file
        ChanNumber % positive integer. Channel number in the original Spike2 file
        SInterval % Sampling interval in Second.
        Length % Length of Data and Time.
    end
    
    properties (Dependent = true, SetAccess = protected)
        % Time % In Second. Time has the same length as Data. Time vector is generated at run-time.
        MaxTime % Time at the end of the file in Seconds.
        TimeUnit % Second (constant in ChanInfo)
    end
    
    properties (Abstract, Dependent = true)
        Data % Data has the same length as Time.
    end
    
    properties (Hidden, SetAccess = protected, GetAccess = public)
        ChanInfo_(1,1)  = ChanInfo();% ChanInfo class object    %TODO need to be changed to ChanInfo_(1,1) ChanInfo % since R2019b
    end
    
    methods
        %% Property get methods; dependent on ChanInfo_ (Hidden)
        function ChanTitle = get.ChanTitle(obj)
            ChanTitle = obj.ChanInfo_.ChanTitle;
        end
        
        function Start = get.Start(obj)
            Start = obj.ChanInfo_.Start;
        end
        
        function SRate = get.SRate(obj)
            SRate = obj.ChanInfo_.SRate;
        end
        
        function DataUnit = get.DataUnit(obj)
            DataUnit = obj.ChanInfo_.DataUnit;
        end
        
        function Header = get.Header(obj)
            Header = obj.ChanInfo_.Header;
        end
        
        function MaxTime = get.MaxTime(obj)
            MaxTime = obj.ChanInfo_.MaxTime;
        end
        
        function ChanNumber = get.ChanNumber(obj)
            ChanNumber = obj.ChanInfo_.ChanNumber;
        end
        
        function Path = get.Path(obj)
            Path = obj.ChanInfo_.Path;
        end
 
        function timeunit = get.TimeUnit(obj)
            timeunit = obj.ChanInfo_.TimeUnit;
        end
        
        function Length = get.Length(obj)
            % if not refer to ChanInfo ... for MarkerChan %TODO
            Length = obj.ChanInfo_.Length; 
            %TODO repeated access of this property slows things down.
            % Potentially, tweak set.Data() method to reflect the updated
            % Data value and hold the property value within the object???
            
        end
        
        function SInterval = get.SInterval(obj)
            %SInterval = get.SInterval(obj)
            %SInterval = 1/obj.SRate;
            
            SInterval = obj.ChanInfo_.SInterval;
        end
        
        function chaninfo = get.ChanInfo_(obj)
            % obj.ChanInfo_.Length = obj.Length; % overwrite %TODO
            chaninfo = obj.ChanInfo_;
            
            
        end
        
        
        %% Property set methods
        
        function obj = set.Start(obj, newstart)
            
            %% parse
            narginchk(2,2);
            
            p = inputParser;
            
%             vf_newstart = @(x) validateattributes(...
%                 {'double'},...
%                 {'scalar', 'real', 'finite', 'nonnan'},...
%                 mfilename, 'newstart', 2);% accept minus value
            
            vf_newstart = @(x) isa(x, 'double') &&...
                isscalar(x) &&...
                isreal(x) &&...
                isfinite(x) &&...
                ~isnan(x) ; % accept minus value
            
            addRequired(p, 'newstart', vf_newstart);
            parse(p, newstart);
            
            %% set
            
            obj.ChanInfo_.Start = newstart; %TODO
            
        end
        
        function obj = set.SRate(obj, newsrate)
            
            %% parse
            narginchk(2,2);
            
            p = inputParser;
            
            vf_srate = @(x) isa(x, 'double') &&...
                isscalar(x) &&...
                isreal(x) &&...
                x > 0;
            
            addRequired(p, 'newsrate', vf_srate);
            parse(p, newsrate);
            
            %% set
            
            obj.ChanInfo_.SRate = newsrate;
            
        end
        
        function obj = set.SInterval(obj, newinterval)
            
            %% parse
            narginchk(2,2);
            
            p = inputParser;
            
            vf_sinterv = @(x) isa(x, 'double') &&...
                isscalar(x) &&...
                isreal(x) &&...
                x > 0;
            
            addRequired(p, 'newinterval', vf_sinterv);
            parse(p, newinterval);
            
            %% set
            
            obj.ChanInfo_.SRate = 1/newinterval;
            
        end
        
        
        function obj = set.ChanTitle(obj, chantitle)
            
            %% parse
            narginchk(2,2);
            
            p = inputParser;
            
            vf_name = @(x) ischar(x) &&...
                isrow(x) || isempty(x);
            
            addRequired(p, 'chantitle', vf_name);
            parse(p, chantitle);
            
            %% set
            
            obj.ChanInfo_.ChanTitle = chantitle;
            
        end
        
        function obj = set.DataUnit(obj, newdataunit)
                                    
            assert( ischar(newdataunit) && isrow(newdataunit) || isempty(newdataunit),...
                'K:Chan:DataUnit:set:newdataunit:invalid',...
            'newdataunit must be a string');
            
            obj.ChanInfo_.DataUnit = newdataunit;
            
        end
        
        function obj = set.ChanNumber(obj,channumber)
            obj.ChanInfo_.ChanNumber = channumber;
        end
        
        function obj = set.Header(obj, newheader)
            
            obj.ChanInfo_.Header = newheader;
            
        end
        
        %%
        
        function time = time(obj)
            % TODO you could implement subsref to this method
            % ex)
            % time(obj, 1..3)
            
            if obj.Length >= 1
                
                %NOTE slightly faster with variables
                start = obj.Start;
                sInt = obj.SInterval;
                len = obj.Length;
                
                time = (start:sInt:(start+sInt*(len - 1)))';%SLOW
                
                % time = (obj.Start:obj.SInterval:(obj.Start+obj.SInterval*(obj.Length - 1)))';
            elseif obj.Length == 0
                time = [];
            end
        end

        
        %         function MaxTime = get.MaxTime(obj)
        %             if obj.Length >= 1
        %                 MaxTime = obj.Start + obj.SInterval*(obj.Length - 1);
        %             elseif obj.Length == 0
        %                 MaxTime = [];
        %             end
        %         end
        
        
        %         function obj1 = minus(obj1, obj2)
        %             %TODO
        %         end
        %
        %         function obj1 = times(obj1, obj2)
        %             %TODO
        %         end
        %
        %         function obj1 = mtimes(obj1, obj2)
        %             %TODO
        %         end
        %
        %         function obj1 = ldivide(obj1, obj2)
        %             %TODO
        %         end
        %
        %         function obj1 = rdivide(obj1, obj2)
        %             %TODO
        %         end
        %
        %         function obj1 = mldivide(obj1, obj2)
        %             %TODO
        %         end
        %
        %         function obj1 =mrdivide(obj1, obj2)
        %             %TODO
        %         end
        
        function outnum = mean(this)
            outnum = mean(this.Data);
        end
        
        function outnum = std(this, flag, dim)
            
            narginchk(1,3);
            
            switch nargin
                case 1
                    outnum = std(this.Data);
                case 2
                    outnum = std(this.Data, flag);
                case 3
                    outnum = std(this.Data, flag, dim);
            end
        end
        
        function outnum = var(this, w)
            
            narginchk(1,2)
            
            switch nargin
                case 1
                    outnum = var(this.Data);
                case 2
                    outnum = var(this.Data, w);
            end
        end
        
        function outnum = mode(this)
            outnum = mode(this.Data);
        end
        
        function outnum = median(this)
            outnum = median(this.Data);
        end
        
        function outnum = max(this)
            outnum = max(this.Data);
        end
        
        function outnum = min(this)
            outnum = min(this.Data);
        end
        
        objout = vertcat(obj1,obj2,varargin) %OK
        
        function tsout = chan2ts(chan) %OK
            tsout = timeseries(chan.Data, chan.time, 'Name', chan.ChanTitle);
            tsout = setuniformtime(tsout,'StartTime',chan.Start,'Interval',chan.SInterval);
            tsout.UserData = chan.Header;
        end
        
        
        function chaninfo = getChanInfo(obj)
            chaninfo = obj.ChanInfo_; %TODO test
        end
        
        function tf = testProperties(obj)
            tf = K_testProperties(obj);
        end
        
    end
    
    methods (Abstract)
        plot(obj)
        
        %       chan2struct(obj) % convert chan to Spike2 structure format
    end
    
    
    methods (Static, Abstract)
        chanout = ts2chan(ts)
        % chanout = ts2chan(ts) convert a timeseries object to a Chan object
        
        %     p = inputParser;
        %     vf = @(x) isa(x, 'timeseries');
        %     addRequired(p, 'ts', vf);
        %     chanout = Chan(ts.Data, ts.Time, ts.Name)
        
        
    end
    
    methods (Static)
                
        function varargout = constructChan(varargin)
            % obj = constructChan(datastruct)
            % [obj1, obj2, obj3, ...]  = constructChan(datastruct1, datastruct2, datastruct3, ...)
            % 
            % obj = constructChan(_______, 'ref', refstruct)
            %
            % INPUT VARIABLES
            % datastruct      structure holding Spike2-derived recording
            %                 data 
            %
            % OPTIONAL PARAM/VAL PAIR
            %
            % 'ref'           ref is structure in event or waveform
            %                 format that is derived from Spike2
            
            
            %% Parse
            narginchk(1,inf)
            
            datastruct = varargin;
            
            ref = [];
            if nargin > 2 && ischar(varargin{nargin-1}) && ...
                    strcmpi(varargin{nargin-1}, 'ref') 
                
                ref = varargin{nargin};
                
                assert( isstruct(ref),...
                    'K:Chan:constructChan:ref:notstruct',...
                    'The reference is not structure');
                
                
                assert( ChanInfo.vf_structEvent(ref) ||  ChanInfo.vf_structWaveform(ref),...
                    'K:Chan:constructChan:ref:noteventorwaveform',...
                    'The reference is not in event or waveform data type');
                
                datastruct = varargin(1:(nargin-2));

            end
            
            
            assert( all(cellfun(@isstruct, datastruct)),...
                'K:Chan:constructChan:datastruct:notstruct',...
                'At least one of the variable is not structure');
            
            
            isevent = cellfun(@(x) ChanInfo.vf_structEvent(x), datastruct);
            iswaveform = cellfun(@(x) ChanInfo.vf_structWaveform(x), datastruct);
            ismarker = cellfun(@(x) ChanInfo.vf_structMarker(x), datastruct);
            
            
            if any(~isevent & ~iswaveform & ~ismarker)
                error( 'K:Chan:constructChan:datastruct:unexpectedtype',...
                    'The %dth input variable is NOT event, waveform or marker type.\n',...
                    find(~isevent & ~iswaveform & ~ismarker));
            end
            
            if any(ismarker) && isempty(ref)
                  error( 'K:Chan:constructChan:datastruct:marker:noref',...
                    'Input includes marker data type but reference ref was not provided.');              
            end
                       
            
            %% Job
            
            n = length(datastruct);
            if nargout < n
                if nargout == 0
                    n = 1;
                else
                    n = nargout;
                end
            end
            

            varargout = cell(1, n);
            for i = 1:n
                if isevent(i)
                    varargout{i} = EventChan(datastruct{i});
                    
                elseif iswaveform(i)
                    varargout{i} = WaveformChan(datastruct{i});
                    
                elseif ismarker(i)
                    varargout{i} = MarkerChan(datastruct{i}, ref);
                    
                end
            end
            
        end
        
    end
    
    
end

