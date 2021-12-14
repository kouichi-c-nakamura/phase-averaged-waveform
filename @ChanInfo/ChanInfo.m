classdef ChanInfo
    % ChanInfo object stores meta data for Chan object
    %
    % Written by Kouichi C. Nakamura Ph.D.
    % MRC Brain Network Dynamics Unit
    % University of Oxford
    % kouichi.c.nakamura@gmail.com
    % 15-Aug-2017 15:23:44
    %
    % See also
    % ChanInfo_test, Chan, WaveformChan, MetaEventChan, EventChan, MarkerChan
    % RecordInfo
    
    properties
        ChanTitle = '';            % corresponds to "title" ffield of Spike2 struct (not the name of struct variable)
        DataUnit = '';             % The measurement unit for Data
        Start = 0;                 % Time at the beginning of the file in Seconds.
        SRate = 1;                 % Sampling Rate in Hertz.
        Header = struct;           % vertcat, resample, extractTime
        Path = '';                 % folder path for source data file
        ChanNumber = 0;            % channel number in original Spike2 file
        Length = 0;                 % Length of Data and Time.
    end
        
    properties (Dependent = true, SetAccess = private)
        % Time                     % In Second. Time has the same length as Data. Time vector is generated at run-time.
        SInterval                  % Sampling interval in Second.
        MaxTime                    % Time at the end of the file in Seconds.
    end
    
    properties (Constant)
        TimeUnit = 'second'; % Second (constant)
    end
    
    properties
        ChanStructVarName = '';    % variable name of struct in .mat file (not necessarily the same as ChanTitle)
        PathRef = '';              % only for MarkerChan. Path of the reference binned data.
        ChanTitleRef = '';         % only for MarkerChan. Path of the reference binned data.
        ChanStructVarNameRef = ''; % only for MarkerChan.
    end
    
    properties (SetAccess = protected)
        Bytes = [];                % file size in bytes.
    end

    properties (Dependent, SetAccess = protected)
        LastModDate  
    end
    
    properties (Hidden, SetAccess = protected)
        DateNum_ = [];
    end
    
    
    methods
        %% Constructor
        function obj = ChanInfo(varargin)
            % obj = ChanInfo(matfilepath, chantitle)
            % obj = ChanInfo(matfilepathMk, chantitleMk, matfilepathBin, chantitleBin)
            %
            % obj = ChanInfo(struct)
            % obj = ChanInfo(structMk, structBin)
            %
            % obj = ChanInfo(Chan) .... tautology?
            %                         Chan.getChanInfo() instead?
            %
            %% Call with file path
            % obj = ChanInfo(matfilepath, chantitle)
            % obtain infomation from a saved .mat file
            %   matfilepath       string. A file path for the .mat file that
            %                     contains Spike2 data.
            %   chantitle             string. the name of the struct in the .mat file
            %                     Note: Not necessarily identical to
            %                     resultant ChanTitle property,
            %
            % obj = ChanInfo(matfilepathMk, chantitleMk, matfilepathBin, chantitleBin)
            % for MarkerChan you need another channel (binned data type) as
            % well
            %    matfilepathMk    string. The file path for the Marker
            %                     channel
            %    chantitleMk          string. the name of the struct
            %                     (MarkerChan, not binned) in the .mat file
            %    matfilepathBin   string. The file path for the binned
            %                     channel as reference
            %    chantitleBin         string. the name of the struct for binned
            %                     data
            %
            %
            %% Call with workspace variables
            % obj = ChanInfo(struct)
            % obtain information from a struct (for WaveformChan or
            % EventChan) in workspace
            % struct             struct in Spike2 format (binned)
            %
            % obj = ChanInfo(structMk, structBin)
            % obtain information from a struct (for MarkerChan or
            % EventChan) in workspace
            % structMk           struct for MarkerChan in Spike2 format (not binned)
            % structBin          struct for another binned channel from the
            %                    same file for reference
            %
            %
            %% Call with Chan object in workspace
            % obj = ChanInfo(chan)
            %   obtain information from a Chan object in workspace
            %
            
            narginchk(0, 4);
            
            p = inputParser;
            
            switch nargin
                case 0
                    % object with defualt values
                    
                
                case 1
                    if isa(varargin{1}, 'struct')
                        %% obj = ChanInfo(struct)
                        % constructor with workspace struct variables
                        %
                        % obtain information from a struct in workspace
                        % WaveformChn or EventChan
                        
                        S = varargin{1};
                        obj = struct2ChanInfo(obj, S);
                        obj.Path = '';
                        
                    elseif isa(varargin{1}, 'Chan')
                        %% obj = ChanInfo(chan)
                        % constructor with a Chan object in workspace
                        
                        %TODO not yet implemented
                        chan = varargin{1};
                        obj = chan.getChanInfo(); %TODO implement the method first
                        
                    else
                        error('K:ChanInfo:ChanInfo:argin:invalid',...
                            'Single input argument must be either struct or Chan class');
                    end
                    
                case 2
                    if isstruct(varargin{1}) && isstruct(varargin{2})
                        %% obj = ChanInfo(structMk, structBin)
                        % constructor with workspace struct variables
                        %
                        % special for MarkerChan
                        % You cannot specify Start of Marker channel on its
                        % own.
                        % It requires another channel for reference.
                        % MarkerChan.Length is differently defined from otheres
                        
                        S = varargin{1};
                        Sbin = varargin{2};
                        
                        obj = struct2ChanInfoMk(obj, S, Sbin);
                        obj.Path = '';
                        
                    else
                        %% obj = ChanInfo(matfilepath, chantitle)
                        % constructor with .mat file that contain data
                        %
                        % for WaveformChan or EventChan
                        
                        matfilepath = varargin{1};
                        chantitle = varargin{2};
                        
                        %% parse
                        vf_matfilepath = @(x) ~isempty(x) && isrow(x) && ischar(x);
                        addRequired(p, 'matfilepath', vf_matfilepath);
                        
                        vf_chantitle =@(x) ~isempty(x) && isrow(x) && ischar(x);
                        addRequired(p, 'chantitle', vf_chantitle);
                        
                        parse(p, matfilepath, chantitle);
                        
                        assert(all(~isempty(regexp(matfilepath, '\.mat$', 'once'))), ...
                            'K:ChanInfo:MarkerChan:matfilepath:invalid',...
                            'matfilepath %s doesn''t contain .mat suffix.',...
                            matfilepath);
                        
                        %% job
                        
                        varname = ChanInfo.matchStructSpk2Title(matfilepath, chantitle);
                       
                        s = load(matfilepath, varname);
                                                
                        S = s.(varname);
                        clear s
                        
                        obj = struct2ChanInfo(obj, S, varname);
                        obj.Path = matfilepath; 
                        
                        listing = dir(matfilepath);
                        obj.DateNum_ = listing.datenum;
                        obj.Bytes = listing.bytes;
                        
                        
                    end
                    
                case 3
                    error('K:ChanInfo:MarkerChan:inputarg:invalid',...
                        'Number of input arguments must be 0, 1, 2 or 4.');
                    
                case 4
                    %% obj = ChanInfo(matfilepathMk, chantitleMk, matfilepathBin, chantitleBin)
                    %
                    % for MarkerChan
                    % You cannot specify Start of Marker channel on its
                    % own.
                    % It requires another channel for reference.
                    % MarkerChan.Length is differently defined from otheres
                    
                    matfilepathMk = varargin{1};
                    chantitleMk = varargin{2};
                    matfilepathBin = varargin{3};
                    chantitleBin = varargin{4};
                    
                    %% parse
                    
                    vf_chantitle =@(x) ~isempty(x) && isrow(x) && ischar(x);
                    addRequired(p, 'chantitleMk', vf_chantitle);
                    
                    vf_chantitle =@(x) ~isempty(x) && isrow(x) && ischar(x);
                    addRequired(p, 'chantitleBin', vf_chantitle);
                    
                    vf_matfilepath = @(x) ~isempty(x) && isrow(x) && ischar(x);
                    
                    addRequired(p, 'matfilepathMk', vf_matfilepath);
                    addRequired(p, 'matfilepathBin', vf_matfilepath);
                    
                    parse(p, matfilepathMk, chantitleMk, matfilepathBin, chantitleBin);
                    
                    assert(all(~isempty(regexp(matfilepathMk, '\.mat$', 'once'))), ...
                        'K:ChanInfo:MarkerChan:matfilepath:invalid',...
                        'matfilepath %s doesn''t contain .mat suffix.',...
                        matfilepathMk);
                    
                    assert(all(~isempty(regexp(matfilepathBin, '\.mat$', 'once'))), ...
                        'K:ChanInfo:MarkerChan:matfilepath:invalid',...
                        'matfilepath %s doesn''t contain .mat suffix.',...
                        matfilepathBin);
                    
                    %% job
                    
                    varnameMk = ChanInfo.matchStructSpk2Title(matfilepathMk, chantitleMk);
                    s = load(matfilepathMk, varnameMk);
                    
                    varnameBin = ChanInfo.matchStructSpk2Title(matfilepathBin, chantitleBin);
                    sbin = load(matfilepathBin, chantitleBin);
                       
                    S = s.(chantitleMk);
                    Sbin = sbin.(chantitleBin);
                    clear s sbin
                    
                    obj = struct2ChanInfoMk(obj, S, Sbin, varnameMk, varnameBin);
                    obj.Path = matfilepathMk;
                    
                    obj.ChanTitleRef = chantitleBin;
                    obj.PathRef = matfilepathBin;
                    
                    listing = dir(matfilepathMk);
                    obj.DateNum_ = listing.datenum;
                    obj.Bytes = listing.bytes;

            end
        end
        
        %% property get methods
        function sinterval = get.SInterval(obj)
            %sinterval = get.SInterval(obj)
            sinterval = 1/obj.SRate;
        end
        
        function maxtime = get.MaxTime(obj)
            if obj.Length >= 1
                maxtime = obj.Start + obj.SInterval*(obj.Length - 1);
            elseif obj.Length == 0
                maxtime = [];
            end
        end
        
        function lastmoddate = get.LastModDate(obj)
            if ~isempty(obj.DateNum_)
                lastmoddate = datestr(obj.DateNum_, 'yyyy-mm-dd HH:MM:SS.FFF');
            else
                lastmoddate = '';
            end
        end
        
        %% property set methods
        function obj = set.ChanTitle(obj, chantitle)
            
            %% parse
            narginchk(2,2);
            
            p = inputParser;
            
            vf_name = @(x) ischar(x) &&...
                isrow(x) || isempty(x);
            
            addRequired(p, 'chantitle', vf_name);
            parse(p, chantitle);
            
            %% set
            
            obj.ChanTitle = chantitle;
            
        end
        
        function obj = set.DateNum_(obj, datenum)
            assert(isempty(datenum) || isscalar(datenum) && isa(datenum, 'double'),...
                'K:RecordInfo:setDateNum_:datenum:invalid',...
                'datenum must be scalar double or empty');
            if isempty(datenum)
                datenum = [];
            end
            obj.DateNum_ = datenum; 
        end

        function obj = set.Bytes(obj, bytes)
            assert(isempty(bytes) || isscalar(bytes) && isa(bytes, 'double'),...
                'K:RecordInfo:setBytes:bytes:invalid',...
                'bytes must be scalar double or empty');
            
            if isempty(bytes)
                bytes = [];
            end
            obj.Bytes = bytes; 
        end
        
        %% common methods
        function time = time(obj)
            % TODO you could implement subsref to this method
            % ex)
            % time(obj, 1..3)
            
            if obj.Length >= 1
                time = (obj.Start:obj.SInterval:(obj.Start+obj.SInterval*(obj.Length - 1)))';
            elseif obj.Length == 0
                time = [];
            end
        end
        
        
        function chan = loadChan(obj)
            % chan = loadChan(obj)
            %
            % load .mat file specified by obj.Path and get an object of
            % Chan class
            %
            % TODO this call must be available for subclass of ChanInfo as
            % well
            % 
            
            assert(~isempty(obj.Path), 'K:ChanInfo:loadChan:Path:empty',...
                'Path property of ChanInfo object %s is empty.', inputname(1));
            
            s = load(obj.Path, obj.ChanStructVarName);
            S = s.(obj.ChanStructVarName);
            clear s;
            
            if ChanInfo.vf_structWaveform(S)
                chan = WaveformChan(S);
                
            elseif ChanInfo.vf_structEvent(S)
                chan = EventChan(obj);
                
                %TODO this doesn't allow you to create EventChan with
                %subclass of ChanInfo
                
                % obj = EventChan(ChanInfo)?
                

            elseif ChanInfo.vf_structMarker(S)
                %% load Sbin for reference
                
                assert(~isempty(obj.PathRef), 'K:ChanInfo:loadChan:Path:empty',...
                'PathRef property of ChanInfo object %s is empty.', inputname(1));
                
                sbin = load(obj.PathRef, obj.ChanStructVarNameRef);
                Sbin = sbin.(obj.ChanStructVarNameRef);
                clear sbin;
                
                chan = MarkerChan(S, Sbin);
                
            else
                error('K:ChanInfo:loadData:S',...
                    'The specified variable name doesn''t seem to be Spike2 struct format');
            end
        end
        
        function tf = testProperties(obj)
            tf = K_testProperties(obj);
        end
        
        function props = get(obj)
            % props = get(obj)
            % props    struct with field names corresponding to property
            %          names of ChanInfo
            
            propnames = properties(obj);
            
            for i = 1:length(propnames)
                props.(propnames{i}) = obj.(propnames{i});
            end

        end
        
        function eqstate = eq(obj1, obj2)
            
            if ~strcmp(class(obj1), class(obj2)) % different class
                eqstate = false;
            else
                
                propnames = properties(obj1);
                eqstate = true;
                
                for i = 1:length(propnames)
                    if  ~isequaln(obj1.(propnames{i}), obj2.(propnames{i}))
                        eqstate = false;
                        
                        %disp(propnames{i});
                        break
                    end
                end
            end
            
        end
       
        function s = saveobj(obj)
            
            s.Bytes = obj.Bytes;
            s.ChanNumber = obj.ChanNumber;
            s.ChanStructVarName = obj.ChanStructVarName;
            s.ChanStructVarNameRef = obj.ChanStructVarNameRef;
            s.ChanTitle = obj.ChanTitle;
            s.ChanTitleRef = obj.ChanTitleRef;
            s.DataUnit = obj.DataUnit;
            s.DateNum_ = obj.DateNum_;
            s.Header = obj.Header;
            s.Length = obj.Length;
            s.Path = obj.Path;
            s.PathRef = obj.PathRef;
            s.SRate = obj.SRate;
            s.Start = obj.Start;
            % s.TimeUnit = obj.TimeUnit;

        end
        
    end
    
    methods (Access = protected)
        
        obj = struct2ChanInfo(obj, S, varname)
            % obj = Struct2ChanInfo(obj, S, varname)
            %
            % protected function
            % used for the constructor ChanInfo()
            % construct obj from struct input S
            % varname is the variable name in workspace
            % EventChan and WaveformChan
            
       
        
        obj = struct2ChanInfoMk(obj, S, Sbin, varname1, varname2)
            % obj = struct2ChanInfoMk(obj, S, Sbin, varname1, varname2)
            %
            % protected function
            % used for the constructor ChanInfo()
            % construct obj from struct input S
            % Sbin is struct from a binned channel for reference
            
    end
    methods (Access = protected, Static)

        
        function varname = matchStructSpk2Title(matfilepath, chantitle)
            %varname = matchStructSpk2Title(matfilepath, chantitle)
            % matfilepath   string. including .mat extension
            % chantitle     string. The title filed of Spike2 struct
            %
            % varname       variable name of the struct that contains
            %               Spike2 data for a channel
            %
            % Check if the variable name or the chantitle field of structSpk2
            % in matfilepath matches chantitle
            
            
            matObj = matfile(matfilepath);
            if nnz(strcmp(chantitle, who(matObj))) == 1 % match with varname
                varname = chantitle;
                
                
            else % no match with varnames
                
                S = load(matfilepath);
                finames = fieldnames(S);
                chantitles = cell(length(finames), 1);
                for i = 1:length(finames)
                    chantitles{i} = S.(finames{i}).title;
                end
                tf = strcmp(chantitle, chantitles);
                clear S
                
                if nnz(tf) == 0
                    error('K:ChanInfo:MarkerChan:chanstructvarname:nomatch',...
                        ['The inout arg %s doesn''t match any of variable ',...
                        'names of struct objects for channels or title fields of them ',...
                        'in the specified .mat file %s'], chantitle, matfilepath);
                elseif nnz(tf) >= 2
                    error('K:ChanInfo:MarkerChan:chanstructvarname:notunique',...
                        ['The inout arg %s matches more than one ',...
                        'title fields of struct objects for channels ',...
                        'in the specified .mat file %s'], chantitle, matfilepath);
                else
                    varname = finames{tf};
                end
            end
        end
    end
    
    methods (Static)
        
        tf = vf_structWaveform(x)
        tf = vf_structEvent(x)
        tf = vf_structMarker(x)
        tf = vf_structFile(x)
        
        function obj = loadobj(s)
            
            obj = ChanInfo;
            obj.Bytes = s.Bytes;
            obj.ChanNumber = s.ChanNumber;
            obj.ChanStructVarName = s.ChanStructVarName;
            obj.ChanStructVarNameRef = s.ChanStructVarNameRef;
            obj.ChanTitle = s.ChanTitle;
            obj.ChanTitleRef = s.ChanTitleRef;
            obj.DataUnit = s.DataUnit;
            obj.DateNum_ = s.DateNum_;
            obj.Header = s.Header;
            obj.Length = s.Length;
            obj.Path = s.Path;
            obj.PathRef = s.PathRef;
            obj.SRate = s.SRate;
            obj.Start = s.Start;
            % obj.TimeUnit = s.TimeUnit;
        end
        
    end
end

