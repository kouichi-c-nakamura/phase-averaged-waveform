classdef WaveformChan < Chan
    % WaveformChan class is a subclass of Chan class. WaveformChan
    % class can store time series data Data in waveform format with an
    % uniform Time vector.
    %
    % Data
    % To support NaNs, Data is now stored as double (64bit) at the expense of
    % memory usage.
    %
    % When you export data from Spike2, choose "Align and bin all data at
    % ..." option and choose either int16, single or double format for
    % waveform channels. Then, the constuctor of WaveformChan class can
    % convert it into int16 for storage.
    %
    %
    % Time
    % Time is a uniform vector of time (double) in second. Time vector is
    % generated at run-time with Start, SRate, and Length properties to
    % reduce disk space requirement. Data and Time have the same Length.
    %
    % Written by Kouichi C. Nakamura Ph.D.
    % MRC Brain Network Dynamics Unit
    % University of Oxford
    % kouichi.c.nakamura@gmail.com
    % 15-Aug-2017 15:25:18
    %
    % See Also Chan, Record, MetaEventChan, EventChan,
    % MarkerChan
    
    
    properties (GetAccess = private, SetAccess = private, Hidden)
        Data_ % hidden property to store data as double
    end
    
    properties (Dependent = true)
        Data % Waveform data (double). Data has the same Length as Time.
    end
    
    properties %(SetAccess = private)
        Scale = 1/double(intmax('int16')); % used for conversion of Data between int16 and double (but not compatible with ChanScale())
        Offset = 0; % used for conversion of Data between int16 and double (but not compatible with ChanScale())
        % 
        % When Spike2 exports waveform data into a MAT file, the structure
        % contains scale and offset fields, and they are used to convert
        % int16 values back into double. However, confusingly, they are NOT
        % compatible with those obtained with ChanScale() and ChanOffset()
        % Spike2 functions.
        %        
        %     y axis value = (16-bit ADC value)*scale + offset
        %
        % On the other hand, ChanSclae() and ChanOffset() are compatible
        % with the following equation:
        %
        %     y axis value = (16-bit value)*chanscale/6553.6 + chanoffset
        %
        % See also
        % int16TOdouble, doubleTOint16
        % int16TOdoubleSpike2, doubleTOint16Spike2, getScaleOffset
    end
    
    methods
        function obj = WaveformChan(varargin)
            % obj = WaveformChan(data, start, srate)
            % obj = WaveformChan(dataInt16, start, srate, scale, offset)
            % obj = WaveformChan(____, chantitle)
            % obj = WaveformChan(S)
            %
            % Constructs an EventChan object. All input arguments are
            % optional.
            %
            % data      a column vector of numeric values (single or double).
            %           Default is [].
            %
            % dataInt16 a column vector of int16 data (signed 16 bit integer).
            %
            % start     a scalar number equal to or larger than 0. Default
            %           is 0.
            %
            % srate     Sampling rate [Hz]. A positive scalar number. Default is
            %           1 Hz.
            %
            % scale, offset
            %           scalar number            
            %           They are used to convert int16 values to doble, but
            %           note that these values are different from those
            %           obtained with ChanScale() or ChanOffset() Spike2
            %           functions.
            %
            %           When Spike2 exports waveform data into a MAT file,
            %           the structure contains scale and offset fields, and
            %           they are used to convert int16 values back into
            %           double. However, confusingly, they are NOT
            %           compatible with those obtained with ChanScale() and
            %           ChanOffset() Spike2 functions.
            %
            %           y axis value = (16-bit ADC value)*scale + offset
            %
            %           On the other hand, ChanSclae() and ChanOffset() are
            %           compatible with the following equation:
            %
            %           y axis value = (16-bit value)*chanscale/6553.6 + chanoffset
            %
            % chantitle      
            %           chantitle of EventChan. Could be different from the
            %           variable name. Default is ''.
            %
            % S         structure
            %           Instead of specifying each parameters, you can
            %           simply use Struct exported by Spike2 as a sole input.
            %           Struct must contain 'values' field holding
            %           double, single or int16 values of a waveform channel.
            %
            
            %% parse inputs
            narginchk(0,6);
            
            info = ChanInfo();
            
            info.Start = 0;
            info.SRate = 1;
            info.ChanTitle = '';
            info.DataUnit = 'mV'; % default
            
            obj.Scale = [];
            obj.Offset = [];
            
            p = inputParser;
            aO = @addOptional;
            aR = @addRequired;
            
            vf_data = @(x) isnumeric(x) &&...
                iscolumn(x) &&...
                isa(x, 'double') || isa(x, 'single');
            
            vf_start = @(x) isa(x, 'double') &&...
                isscalar(x) &&...
                isreal(x) &&...
                ~isnan(x) ; % accept minus value
            
            vf_srate = @(x) isa(x, 'double') &&...
                isscalar(x) &&...
                isreal(x) &&...
                x > 0;
            
            vf_name = @(x) ischar(x) &&...
                isrow(x) || isempty(x);
            
            if nargin >= 0 && nargin <=4
                if  nargin ==0 || ~isa(varargin{1}, 'struct')
                    % obj = WaveformChan(data, start, srate, chantitle)
                    
                    aO(p, 'data', [], vf_data);
                    aO(p, 'start', 0, vf_start);
                    aO(p, 'srate', 1, vf_srate);
                    aO(p, 'chantitle', '', vf_name);
                    
                    parse(p, varargin{:});
                    
                    data = double(p.Results.data);
                    
                    obj.Offset = mean([max(data), min(data)]);
                    Mi = double(intmax('int16'));
                    Mx = nanmax(data - obj.Offset);
                    obj.Scale = Mx/Mi;
                    
                    obj.Data_ = data;
                    info.Start = p.Results.start;
                    info.SRate = p.Results.srate;
                    info.ChanTitle = p.Results.chantitle;
                    info.Length = size(obj.Data_,1);
                    
                    obj.ChanInfo_ = info; %TODO
                    
                else % if isa(varargin{1}, 'struct')
                    % obj = WaveformChan(Struct)
                    
                    narginchk(1,1);
                    
                    S = varargin{1};
                    
                    if ~ChanInfo.vf_structWaveform(S)
                        error('K:Chan:WaveformChan:WaveformChan:struct:invalid',...
                            'struct doesn''t seem Spike2 waveform format');
                    end
                    
                    
                    if iscolumn(S.values) && ~isa(S.values, 'int16') &&...
                            isa(S.values, 'double') || isa(S.values, 'single')
                        % S.values is double or single
                        
                        data = double(S.values);
                        
                        % obj.Offset = mean([max(data), min(data)]);
                        % Mi = double(intmax('int16'));
                        % Mx = max(data - obj.Offset);
                        % obj.Scale = Mx/Mi;
                        
                        % no conversion needed
                        obj.Offset = S.offset;
                        obj.Scale = S.scale;
                        
                        obj.Data_ = data;
                        info.Start = S.start;
                        info.SRate = 1/S.interval;
                        info.ChanTitle = S.title;
                        info.DataUnit = S.units;
                        info.Header = rmfield(S, 'values');
                        info.Length = size(obj.Data_,1);

                        obj.ChanInfo_ = info;
                        
                        
                    elseif iscolumn(S.values) && isa(S.values, 'int16')
                        % S.values is int16
                        
                        obj.Scale = S.scale;
                        obj.Offset = S.offset;
                        
                        obj.Data_ = WaveformChan.int16TOdouble(S.values, obj.Scale, obj.Offset);

                        
                        info.Start = S.start;
                        info.SRate = 1/S.interval;
                        info.ChanTitle = S.title;
                        info.DataUnit = S.units;
                        info.Header = rmfield(S, 'values');
                        info.Length = size(obj.Data_,1);
                        
                        obj.ChanInfo_ = info; %TODO
                        
                    else
                        error('K:WaveformChan:data:invalidtype','S.values must be int16, single, or double type');
                    end
                    
                end
                
                
            elseif nargin >= 5
                switch class(varargin{1})
%                     case {'double','single'}
%                         % obj = WaveformChan(data, start, srate, scale, offset, chantitle)
% 
%                         
%                         aR(p, 'data',  vf_data);
%                         aR(p, 'start', vf_start);
%                         aR(p, 'srate', vf_srate);
%                         
%                         vf_scaleoffset = @(x) ~isempty(x) &&...
%                             isnumeric(x) &&...
%                             isscalar(x) &&...
%                             isa(x, 'double') ;
%                         aR(p, 'scale', vf_scaleoffset);
%                         aR(p, 'offset', vf_scaleoffset);
%                         
%                         aO(p, 'chantitle', '', vf_name);
%                         
%                         parse(p, varargin{:});
%                         
%                         data = double(p.Results.data);
%                         
%                         obj.Offset = p.Results.offset;
%                         obj.Scale = p.Results.scale;
%                         
%                         obj.Data_ = data;
%                         info.Start = p.Results.start;
%                         info.SRate = p.Results.srate;
%                         info.ChanTitle = p.Results.chantitle;
%                         
%                         obj.ChanInfo_ = info;
                                                
                        
                    case 'int16'                        
                        % obj = WaveformChan(dataInt16, start, srate, scale, offset, chantitle)
                        
                        vf_data2 = @(x) ~isempty(x) &&...
                            isnumeric(x) &&...
                            iscolumn(x) &&...
                            isa(x, 'int16') ;
                        aR(p, 'data', vf_data2);
                        
                        vf_start2 = @(x)  ~isempty(x) &&...
                            isnumeric(x) &&...
                            isscalar(x) &&...
                            x >= 0 ;
                        aR(p, 'start', vf_start2);
                        
                        vf_srate2 = @(x) ~isempty(x) &&...
                            isnumeric(x) &&...
                            isscalar(x) &&...
                            x > 0;
                        aR(p, 'srate', vf_srate2);
                        
                        vf_scaleoffset = @(x) ~isempty(x) &&...
                            isnumeric(x) &&...
                            isscalar(x) &&...
                            isa(x, 'double') ;
                        aR(p, 'scale', vf_scaleoffset);
                        aR(p, 'offset', vf_scaleoffset);
                        
                        aO(p, 'chantitle', '', vf_name);
                        
                        parse(p, varargin{:});
                        
                        obj.Scale = p.Results.scale;
                        obj.Offset = p.Results.offset;
                        
                        obj.Data_ = WaveformChan.int16TOdouble(p.Results.data, obj.Scale, obj.Offset);
                        
                        info.Start = p.Results.start;
                        info.SRate = p.Results.srate;
                        info.ChanTitle = p.Results.chantitle;
                        info.Length = size(obj.Data_,1);

                        
                        obj.ChanInfo_ = info; %TODO
                        
                    otherwise
                        error('Unexpected syntax.')
                end
                
              
                
            else
                error('K:WaveformChan:nargin', 'Number of input arguments allowed is 0, 4 or 6');
            end
            
        end
        
        function Data = get.Data(obj)
            
            Data = obj.Data_;
            
        end
        
        function obj = set.Data_(obj, newdata_)
           
            obj.Data_ = newdata_;
            
            info = obj.ChanInfo_;
            info.Length = size(obj.Data_,1);
            obj.ChanInfo_ = info;
            
        end
        
        
        function obj = set.Data(obj, newdata)
            % obj = set.Data(obj, newdata)
            %
            % newdata     numeric column vector in double or single class
            
            %% parse inputs
            narginchk(2,2);
            
            p = inputParser;
            vf_newdata = @(x) isnumeric(x) &&...
                iscolumn(x) &&...
                all(isa(x, 'double')) || ...
                all(isa(x, 'single')) ;
            addRequired(p, 'newdata', vf_newdata);
            parse(p, newdata);
            
            %% job
            
            % update obj.Offset and obj.Scale if necessary
            %             if max(newdata) > max(obj.Data) || ...
            %                     min(newdata) < min(obj.Data)
            %
            %                 obj.Offset = mean([max(max(obj.Data), max(newdata)),...
            %                     min(min(obj.Data), min(newdata))]);
            %                 Mi = double(intmax('int16'));
            %                 Mx = max(max(obj.Data - obj.Offset), (max(newdata) - obj.Offset));
            %                 obj.Scale = Mx/Mi;
            %             end
            %
            %             WtoInt16 = @(x) int16((x-obj.Offset)./obj.Scale);
            %             obj.DataInt16 = WtoInt16(newdata);
            
            obj.Data_ = double(newdata); % store as double
            
        end
        
        function obj2 = resample(obj, newRate,varargin)
            % obj2 = resample(obj, newRate)
            % Resample Data with a new sampling rate newRate
            %
            % See Also WaveformChan
            
            newdata = resample(obj.Data, newRate, obj.SRate);
            obj2 = WaveformChan(newdata, obj.Start, newRate, obj.ChanTitle);
            obj2.DataUnit = obj.DataUnit;
            obj2.Header = obj.Header;
        end
        
        function obj2 = setDataInt16(~, dataInt16, scale, offset)
            % obj2 = setDataInt16(~, dataInt16, scale, offset)
            %
            % Allows you to directly set Data with data in int16 format,
            % with scale, and offset parameters.
            %
            % See Also WaveformChan
            
            %% parse
            narginchk(4, 4);
            
            p = inputParser;
            
            vf_dataInt16 = @(x) ~isempty(x) &&...
                isnumeric(x) &&...
                iscolumn(x) &&...
                isa(x, 'int16') ;
            
            addRequired(p, 'dataInt16', vf_dataInt16);
            
            vf_double = @(x) ~isemptry(x) &&...
                isnumeric(x) &&...
                isscalar(x) &&...
                isa(x, 'double') ;
            
            addRequired(p, 'scale', vf_double);
            addRequired(p, 'offset', vf_double);
            
            parse(p, dataInt16, scale, offset);
            
            %% job
            obj2.Scale = scale;
            obj2.Offset = offset;
            
            Int16toW = @(x) double(x).*obj.Scale + obj.Offset;
            
            obj2.Data_ = Int16toW(dataInt16); %TODO
            
        end
        
        function struct = chan2struct(obj, varargin)
            % struct = chan2struct(obj, format)
            %
            % format    optional string
            %           'int16' (default), 'single' or 'double'
            
            %% parse
            narginchk(1,2);
            
            p = inputParser;
            
            vf_format = @(x) ischar(x) && ...
                isrow(x) && ...
                any(strcmpi(x, {'int16','single','double'}));
            
            addOptional(p, 'format', 'int16', vf_format);
            
            parse(p, varargin{:});
            
            format = p.Results.format;
            
            %% job
            
            header = obj.Header;
            
            struct.title = obj.ChanTitle;
            if isempty(header)
                struct.comment = 'No comment';
            else
                struct.comment = obj.Header.comment;
            end
            struct.interval = obj.SInterval;
            struct.scale = obj.Scale;
            struct.offset = obj.Offset;
            struct.units = obj.DataUnit;
            struct.start = obj.Start;
            struct.length = obj.Length;
            
            if ~isempty(header)
                if isfield(header, 'title')
                    header = rmfield(header, 'title');
                end
                if isfield(header, 'comment');
                    header = rmfield(header, 'comment');
                end
                if isfield(header, 'interval');
                    header = rmfield(header, 'interval');
                end
                if isfield(header, 'scale')
                    header = rmfield(header, 'scale');
                end
                if isfield(header, 'offset')
                    header = rmfield(header, 'offset');
                end
                if isfield(header, 'units')
                    header = rmfield(header, 'units');
                end
                if isfield(header, 'start')
                    header = rmfield(header, 'start');
                end
                if isfield(header, 'length')
                    header = rmfield(header, 'length');
                end
            end
            
            switch lower(format)
                case 'int16'
                    if ~isempty(obj.Scale) && ~isempty(obj.Offset)
                        data = obj.Data;
                        indnan = isnan (data);
                        if any(indnan)
                            warning('K:Chan:WaveformChan:chan2struct:NaN',...
                                ['obj.Data contains NaNs while int16 is chosen for output format.\n',...
                                'NaNs are replaced with zeros.']);
                            zr = zeros(size(data));
                            data(indnan) = zr(indnan);
                            clear zr
                        end
                        
                        % WtoInt16 = @(x) int16((x-obj.Offset)./obj.Scale);
                        % struct.values = WtoInt16(data);
                        struct.values = WaveformChan.doubleTOint16(data, obj.Scale, obj.Offset);
                    else
                        if isempty(obj.Scale)
                            error('K:Chan:WaveformChan:chan2struct:Scale',...
                                'obj.Scale is empty');
                        elseif isempty(obj.Offset)
                            error('K:Chan:WaveformChan:chan2struct:Offset',...
                                'obj.Offset is empty');
                        end
                    end
                    
                case 'single'
                    struct.values = single(obj.Data);
                case 'double'
                    struct.values = obj.Data;
            end
            
            %% save the other fields of Header property
            if ~isempty(header)
                finames = fieldnames(header);
                
                for i = 1:length(finames)
                    struct.(finames{i}) = header.(finames{i});
                end
            end
            
        end
        
        
        function s = saveobj(obj)
            
            s.ChanInfo_ = obj.ChanInfo_;
            s.Data_ = obj.Data_;
            s.Offset = obj.Offset;
            s.Scale = obj.Scale;
            
        end
        
    end
    
    %----------------------------------------------------------------------
    
    methods (Static)
        function chanout = ts2chan(ts)
            % chanout = WaveformChan.ts2chan(ts)
            %
            % Converts a timeseris object to WaveformChan object.
            %
            % See Also WaveformChan
            
            p = inputParser;
            vf = @(x) isa(x, 'timeseries');
            addRequired(p, 'ts', vf);
            parse(p, ts);
            
            sInterval = ts.TimeInfo.Increment;
            if isempty(sInterval)
                sInterval = (ts.Time(end) - ts.Time(0))/(ts.Length - 1);
            end
            
            chanout = WaveformChan(ts.Data, ts.Time(1), 1/sInterval, ts.Name);
            chanout.Header = ts.UserData;

        end
        
        function w = int16TOdouble(dataInt16, scale, offset)
            % w = WaveformChan.int16TOdouble(dataInt16, scale, offset)
            %
            % "Mat file data format" > "Waveform RealWave and binned
            % RealMark data drawn as a waveform"
            %
            % The scale and offset values convert 16-bit ADC data to real
            % values using the equation:
            %
            % real = (16-bit ADC value)*scale + offset
            %
            %
            % NOTE
            % The scale and offset being used here are incompatible with
            % those provided by ChanScale() and ChanOffset()
            %
            % 
            % See also
            % WaveformChan.doubleTOint16
            % WaveformChan.int16TOdoubleSpike2
            % ChanScale() [Spike2 doc]
            % int16
            
            arguments
               dataInt16 (:,1) int16 %TODO does this accept double input without an error?
               scale (1,1) double {mustBeReal}
               offset (1,1) double {mustBeReal}
                
            end
            
            w = double(dataInt16).*scale + offset;

        end
        
        function w = int16TOdoubleSpike2(Spike2Int16, chanscale, chanoffset)
            % w = WaveformChan.int16TOdouble(Spike2Int16, chanscale, chanoffset)
            %
            % 
            % y axis value = (16-bit value)*scale/6553.6 + offset
            %
            % INPUT ARGUMENTS
            % Spike2Int16    Signed int16 used internally by Spike2 for
            %                waveform data
            %
            % chanscale    output of Spike2 ChanScale()  
            %
            % chanoffset   output of Spike2 Spike2offset()
            %
            % NOTE
            % The scale and offset being used here are incompatible with
            % those provided by ChanScale() and ChanOffset()
            %
            % 
            % See also
            % WaveformChan.int16TOdouble
            % ChanScale(), ChanOffset() [Spike2 doc]
            % int16
            
            arguments
                Spike2Int16 int16 %TODO does this accept double input without an error?
                chanscale (1,1) double {mustBeReal}
                chanoffset (1,1) double {mustBeReal}
            end
            
            w = double(Spike2Int16).*chanscale./6553.6 + chanoffset;

        end
        
        
        function i = doubleTOint16(dataDouble, scale, offset)
            % i= WaveformChan.doubleTOint16(dataDouble, scale, offset)
            %
            % real = (16-bit ADC value)*scale + offset
            %
            % See also
            % WaveformChan.int16TOdouble
            % WaveformChan.doubleTOint16Spike2
            % int16, double
            
            arguments
                dataDouble double
                scale (1,1) double {mustBeReal}
                offset (1,1) double {mustBeReal}
            end
            
            i = int16((dataDouble - offset)./scale);

        end
        
        function i = doubleTOint16Spike2(Spike2Double, Spike2scale, Spike2offset)
            % i= WaveformChan.doubleTOint16(dataDouble, scale, offset)
            %
            % y axis value = (16-bit value)*scale/6553.6 + offset
            %
            % See also
            % WaveformChan.doubleTOint16
            % int16, double
            
            arguments
                Spike2Double double
                Spike2scale (1,1) double {mustBeReal}
                Spike2offset (1,1) double {mustBeReal}
            end
                        
            i = int16((Spike2Double - Spike2offset).*6553.6./Spike2scale);

        end
        
        function [scale,offset] = getScaleOffset(doubledata,fmt)
            % getScaleOffset is a static method of WaveformChan class and
            % returns scale and offset for the given doubledata to
            % accomodate it into int16 format.
            %
            % SYNTAX
            % [scale,offset] = WaveformChan.getScaleOffset(doubledata)
            % [scale,offset] = WaveformChan.getScaleOffset(doubledata,'int16')
            % [chanscale,chanoffset] = WaveformChan.getScaleOffset(doubledata,'spike2')
            %
            %
            % INPUT ARGUMENTS
            % data        a column vector of double | array
            %             If array, array is internally converted to a
            %             vector by data(:)
            %
            %
            % fmt         'spike2' | 'int16' (default)
            %             (Optional)
            %             If 'spike2', the following equation is used:
            %
            %               y axis value = (16-bit ADC value)*scale/6553.6 + offset
            %
            %             If 'int16', the following equation is used:
            %
            %               real = (16-bit ADC value)*scale + offset
            %
            %
            % OUTPUT ARGUMENTS
            % scale, offset
            %             double scalar 
            %
            % chanscale,chanoffset
            %             double scalar
            %             Note that they are not supposed to be used as
            %             Scale and Offset properties of WaveformChan.
            %
            % your signature comes here
            %
            % See also
            % WaveformChan.doubleTOint16, WaveformChan.doubleTOint16Spike2,
            % WaveformChan.int16TOdouble, WaveformChan.int16TOdoubleSpike2
            
            arguments
                doubledata (:,:) double
                fmt (1,:) char {mustBeMember(fmt,{'spike2','int16'})} = 'int16'
            end
            
            maxdata = max(doubledata(:));
            mindata = min(doubledata(:));
            
            switch fmt
                case 'spike2'
                    
%                     max(data) = double(intmax('int16'))*scale/6553.6 + offset;
%                     
%                     min(data) = double(intmin('int16'))*scale/6553.6 + offset;
                    
%                     double(int16(5))
%                     
%                     double(intmax('int16'))
%                     
%                     double(intmax('int16') - intmin('int16'))*scale/6553.6  = (max(data) - min(data))
                    
           
                    %NOTE y axis value = (16-bit value)*scale/6553.6 + offset
                    scale = ((maxdata - mindata)*6553.6) / double(intmax('int16') - intmin('int16'));
                    offset = maxdata - double(intmax('int16')) *scale/6553.6;
                                   
                case 'int16'
                    % real = (16-bit ADC value)*scale + offset
                                        
                    scale = (maxdata -mindata) / double(intmax('int16') - intmin('int16'));
                    
                    offset = maxdata - double(intmax('int16')) * scale;
                    
            end
            
        end
        
        
        function obj = loadobj(s)
            
            obj = WaveformChan;
            obj.ChanInfo_ = s.ChanInfo_;
            obj.Data_ = s.Data_;
            obj.Offset = s.Offset;
            obj.Scale = s.Scale;
            
        end
    end
    
end

