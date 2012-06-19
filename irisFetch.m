classdef irisFetch
    %IRISFETCH allows seamless access to data stored within the  IRIS-DMC
    %
    %Each method retrieves data from the IRIS-DMC, returned in array of
    %matlab structures.  To see the available methods, type:
    %
    %  methods(irisFetch)
    %
    %
    %OVERVIEW OF METHODS:
    %
    %  IRISFETCH.Traces() retrieves sac-equivelent waveform with channel
    %  metadata from the IRIS-DMC
    %
    %    Example:
    %      tr = IRISFETCH.Traces('IU','ANMO','10','BHZ',...
    %           '2010-05-10 03:00:00','2010-05-10 03:04:00')
    %
    %
    %  IRISFETCH.Stations() retrieves station metadata from the IRIS-DMC.
    %  This data may be retrieved at a variety of detail levels. From
    %  broadest to most specific, these are NETWORK, STATION, CHANNEL, and
    %  RESPONSE.
    %
    %    Example:
    %      s = IRISFETCH.Stations('channel','*','ANMO','10','BH?')
    %
    %
    %  IRISFETCH.Events() retrieves event metadata from the IRIS-DMC
    %    Example:
    %      ev = IRISFETCH.Events('MinimumMagnitude',6.0,...
    %          'minimumLatitude',45,'maximumLatitude', 60,...
    %          'minimumLongitude', -150,'maximumLongitude', -90)
    %
    %
    %
    %REQUIREMENTS:
    %
    %  IRISFETCH requires the latest version of the Web Services Library,
    %  which is a java .jar file available from:
    %
    %  http://www.iris.edu/manuals/javawslibrary/#download
    %
    %  This jar file must be added to your MATLAB path, which may be done
    %  in a variety of ways.  One common way is to include a javaaddpath
    %  statement in the startup.m file.  For more details, consult MATLAB's
    %  documentation for 'Bringing Java Classes and Methods into MATLAB
    %  Workspace'.
    %
    %
    %FOR FURTHER GUIDANCE:
    %
    %  A discussion about using MATLAB to access data from the IRIS-DMC can
    %  be found at:
    %
    %  http://www.iris.edu/manuals/javawslibrary/matlab/
    %
    %
    %
    %see also IRISFETCH.TRACES IRISFETCH.STATIONS IRISFETCH.EVENTS
    %IRISFETCH.FLATTENTOCHANNEL JAVAADDPATH
    
    % Celso Reyes, Rich Karstens
    % IRIS-DMC
    % February 2012
    
    properties
    end %properties
    
    methods(Static)
        function v = version()
            v = '1.0.0';
        end
        
        function ts = Traces(network, station, location, channel, startDateStr, endDateStr, quality, verbosity )
            %irisFetch.Traces accesses waveform with associated channel metadata
            %
            %USAGE
            % tr = irisFetch.Traces(network, station, location,
            % channel, startDate, endDate) will use channel and date
            % criteria to retrieve one or more seismic traces, which are
            % stored as structures containing typical SAC-equivalent
            % metadata.
            %
            % network, station, location, and channel all accept '*' and
            % '?' wildcards, as well as comma separated lists.
            %
            % tr = irisFetch.Traces(net, sta, loc, ch, startt, endt, quality)
            % allows you to specify a quality factor, such as 'B' or 'M'.
            % If not specified, quality defaults to 'B'
            %
            %
            %NOTICE:
            %  startDate and endDate must be formatted thusly:
            %      'YYYY-MM-DD hh:mm:ss'
            %      'YYYY-MM-DD hh:mm:ss.sss'
            %
            %EXAMPLES:
            %
            %    Example 1:
            %      % get 3 components for a 4-minute time period,
            %      % using the '?' wildcard.
            %      ts = irisFetch.Traces('IU','ANMO','10','BH?',...
            %           '2010-05-10 03:00:00','2010-05-10 03:04:00')
            %
            %    Example 2:
            %      % get z-channels for a comma-separated list of locations
            %      % while specifying a quality of 'B'
            %      ts = irisFetch.Traces('IU','ANMO','00,10','BHZ',...
            %           '2010-05-10 03:00:00','2010-05-10 03:04:00', 'B')
            %
            %    Example 3:
            %      % get z-data for all BHZ stations that belong to the IU
            %      % network and have a location code of '00'
            %      ts = irisFetch.Traces('IU','*','00','BHZ',...
            %           '2010-05-10 03:00:00','2010-05-10 03:04:00')
            
            if ~exist('verbosity', 'var')
                verbosity = false;
            end
            
            if ~exist('quality', 'var')
                quality = 'B';
            end
            
            try
                traces = edu.iris.dmc.ws.extensions.fetch.TraceData.fetchTraces(network, station, location, channel, startDateStr, endDateStr, quality, verbosity);
                ts = irisFetch.convertTraces(traces);
                clear traces;
            catch je
                if strcmp(je.identifier, 'MATLAB:undefinedVarOrClass')
                    errtext=['The Web Services library does not appear to be in the javaclasspath.\n',...
                        'Please download the latest version from \n',...
                        'http://www.iris.edu/manuals/javawslibrary/download/IRIS-WS-latest.jar\n ',...
                        'and then add it to your classpath. \n'];
                    warning('IRISFETCH:NoIrisWSJarInstalled',errtext);
                    isSilent = true; %suppress messages within the connectTo...
                    success = irisFetch.connectTo_IRIS_WS_jar(isSilent);
                    if success
                        disp('irisFetch.connectTo_IRIS_WS_jar() has was able to connect you to the appropriate java library. Continuing...');
                        try
                            traces = edu.iris.dmc.ws.extensions.fetch.TraceData.fetchTraces(network, station, location, channel, startDateStr, endDateStr, quality, verbosity);
                            ts = irisFetch.convertTraces(traces);
                            clear traces;
                        catch je2
                            rethrow(je2)
                        end
                    else
                        error('IRISFETCH:Traces:UnableToInstallIrisWSJar',...
                            'irisFetch was unable to recover, Please download and add the latest IRIS-WS-JAR to your javaclasspath');
                    end
                else
                    fprintf('Exception occured in IRIS Web Services Library: %s\n', je.message);
                end
            end
        end
        
        
        function [networkStructure, urlParams] = Stations(detailLevel, network, station, location, channel, varargin)
            %irisFetch.Stations retrieves station metadata from IRIS-DMC
            %
            %USAGE:
            %  s = irisFetch.Stations(detail, network, station, location,
            %  Channel),  where detail is one of "NETWORK", "STATION",
            %  "CHANNEL", or "RESPONSE".  These five parameters are
            %  required for all queries, but may be wildcarded by using []
            %  for their values.
            %
            %  Network, station, location, and channel parameters are
            %  passed directly to the java library, so lists (separated by
            %  commmas) and wildcards (? and *) are accepted.
            %
            %  [s, myParams] = irisFetch.Stations( ... ) also returns the
            %  URL parameters that were used to make the query.
            %
            %  s = irisFetch.Stations( ... , param1, value1 [, ...]])
            %  allows any number of parameter-value pairs to be included in
            %  the selection criteria.  All of the StationCriteria 'set'
            %  methods are supported as parameter pairs, as are a couple
            %  special parameter shortcuts.
            %
            % To determine which parameters can be set,
            %    1) type the following:
            %         methods(edu.iris.dmc.ws.criteria.StationCriteria)
            %    2) look at the list:
            %         All methods starting with "set" are accessible via
            %         the parameter list.
            %
            %   Example:
            %     n = irisFetch.Stations(....,'endbefore',now, 'maxlongitude',-100)
            %
            %     which would invoke the following setters:
            %       stationCriteria.setEndBefore( now )
            %       stationCriteria.setMaxLongitude(-100)
            %
            %
            %
            %  Usable parameters are listed below.  For detailed
            %  descriptions of their effect and use, consult the station
            %  webservice webpage, available at:
            %
            %  http://www.iris.edu/ws/station/
            %
            %PARAMETER LIST (as of 2/1/2012)
            %  'MinLatitude', 'MaxLatitude', 'MinLongitude',
            %  'MaxLongitude', 'Latitude', 'Longitude',
            %  'MinRadius','MaxRadius', 'StartAfter', 'EndAfter',
            %  'StartBefore', 'EndBefore', 'StartTime', 'EndTime',
            %  'UpdatedAfter'
            %
            %CONVENIENCE PARAMETERS
            %   'boxcoordinates'    : [minLat, minLon, maxLat, maxLon]
            %                           % use NaN as a wildcard
            %   'radialcoordinates' : 1x4 double :
            %                           [Lat, Lon, MaxRadius, MinRadius]
            %                           % MinRadius is optional
            %
            %
            %ADDITIONAL DETAILS: USING CELLS
            %  Network, station, location, and channel may be strings or
            %  cell arrays. Each element of the cell array is added to the
            %  search criteria as "or".  That is:
            %
            %      irisFetch.Stations(detail,'AB,BC,CD',...)
            %
            %  is equivalent to
            %
            %      irisFetch.Stations(detail,{'AB','BC','CD'},...)
            %
            %  However, the former makes one call to addNetwork() while
            %  the latter makes three calls to addNetwork(). The net effect
            %  may be the same, but the execution is different.
            %
            %
            %WORKING WITH THE RESULTS
            %  The results are returned in a structure tree, with the same
            %  hierarchy found in the StationXML.  To make this easier to
            %  work with within matlab, you can use the
            %  irisFetch.flattenToChannel routine to shuffle the results
            %  into a 1xN array of Channels.  This only works if the detail
            %  level was "Channel" or "Response".
            %
            %SEE ALSO IRISFETCH.FLATTENTOCHANNEL
            
            %END OF HELP
            
            %-------------------------------------------------------------
            % An alternate base URL can be specified by providing the
            % parameter pair (...,'BASEURL',alternateURL)
            %
            % This is mostly useful in testing, and is likely not relevent
            % for users
            %-------------------------------------------------------------
            
            import edu.iris.dmc.*
            
            if nargin==1 && strcmpi(detailLevel,'help')
                disp('HELP request recognized, but not implemented');
                return
            elseif nargin < 5
                error('not enough arguments.%d',nargin);
            end
            
            try
                outputLevel = ws.criteria.OutputLevel.(upper(detailLevel));
            catch
                error('IRISFETCH:invalidOutputLevel',...
                    'The selected outputLevel [''%s''] was not recognized.',...
                    upper(detailLevel));
            end
            
            
            indexOffsetOfBASEURL=find(strcmpi(varargin(1:2:end),'BASEURL'),1,'first') * 2;
            try
                baseURL = varargin{indexOffsetOfBASEURL}
            catch
            end
            
            
            
            serviceManager = ws.service.ServiceUtil.getInstance();
            serviceManager.setAppName('MATLAB/irisFetch');
            if exist('baseURL','var')
                varargin(indexOffsetOfBASEURL-1:indexOffsetOfBASEURL) = [];
                service = serviceManager.getStationService(baseURL);
            else
                service = serviceManager.getStationService();
            end
            criteria = ws.criteria.StationCriteria;
            
            %----------------------------------------------------------
            % Deal with the Station/Network/Channel/Location parameters
            % These are treated separately, as they're "add" & not "set"
            % Each may handle multiple strings (as a cell array)
            %----------------------------------------------------------
            
            
            criteria = irisFetch.addCriteria(criteria, network, 'addNetwork');
            criteria = irisFetch.addCriteria(criteria, station, 'addStation');
            criteria = irisFetch.addCriteria(criteria, location,'addLocation');
            criteria = irisFetch.addCriteria(criteria, channel, 'addChannel');
            
            criteria = irisFetch.setCriteria(criteria, varargin);
            
            j_networks = service.fetch(criteria, outputLevel);
            networkStructure = irisFetch.parseCollection(j_networks);
            
            if nargout == 2
                urlParams = criteria.toUrlParams;
            end
        end
        
        %%
        function [events, urlParams] = Events(varargin)
            %irisFetch.Events retrieves event data from the IRIS-DMC
            %
            %USAGE:
            %  ev = irisFetch.Events(param1, value1 [, ...]) retrieves
            %  event data from the IRIS-DMC database, returning it as a
            %  matlab structure.  An arbitrary number of parameter-value
            %  pairs may be specified in order to narrow down the search
            %  results.
            %
            %  [ev, myParams] = irisFetch.Events( ... ) also returns the
            %  URL parameters that were used to make the query.
            %
            %  Usable parameters are listed below.  For detailed
            %  descriptions of their effect and use, consult the webservice
            %  webpage for events, available at:
            %
            %  http://www.iris.edu/ws/event/
            %
            %PARAMETER LIST (as of 2/1/2012)
            %  'EventId'
            %  'FetchLimit'
            %  'MinLatitude'
            %  'MaxLatitude'
            %  'MinLongitude'
            %  'MaxLongitude'
            %  'Latitude'
            %  'Longitude'
            %  'MinRadius'
            %  'MaxRadius'
            %  'StartTime'
            %  'EndTime'
            %  'UpdatedAfter'
            %  'MagnitudeType'
            %  'MinimumMagnitude'
            %  'MaximumMagnitude'
            %  'catalogContains'
            %  'contributorContains'
            %  'includeArrivals'
            %  'includeallMagnitudes'
            %  'preferredOnly'
            %
            %
            %CONVENIENCE PARAMETERS
            %   'boxcoordinates'    : [minLat, minLon, maxLat, maxLon]   % use NaN as a wildcard
            %   'radialcoordinates' : [Lat, Lon, MaxRadius, MinRadius]   % MinRadius is optional
            %
            
            %END OF HELP
            
            %-------------------------------------------------------------
            % An alternate base URL can be specified by providing the
            % parameter pair (...,'BASEURL',alternateURL)
            %
            % This is mostly useful in testing, and is likely not relevent
            % for users
            %-------------------------------------------------------------
            
            import edu.iris.dmc.*
            
            serviceManager = ws.service.ServiceUtil.getInstance();
            serviceManager.setAppName('MATLAB/irisFetch');
            
            
            indexOffsetOfBASEURL=find(strcmpi(varargin(1:2:end),'BASEURL'),1,'first') * 2;
            
            try
                baseURL = varargin{indexOffsetOfBASEURL}
            catch
            end
            
            if exist('baseURL','var')
                varargin(indexOffsetOfBASEURL-1:indexOffsetOfBASEURL) = [];
                service = serviceManager.getEventService(baseURL);
            else
                service = serviceManager.getEventService();
            end
            
            criteria = ws.criteria.EventCriteria;
            criteria = irisFetch.setCriteria(criteria, varargin);
            disp('fetching from IRIS-DMC')
            j_events = service.fetch(criteria);
            disp('parsing into MATLAB structures')
            events = irisFetch.parseCollection(j_events);
            if nargout == 2
                urlParams = criteria.toUrlParams;
            end
        end
        
        
        %%
        function success = connectTo_IRIS_WS_jar(isSilent)
            %irisFetch.connectTo_IRIS_WS_jar tries to set up the jar for
            %use within MATLAB for this session
            %
            %USAGE:
            %  success = connectTo_IRIS_WS_jar(isSilent)
            %
            %  This routine searches the javaclasspath for the IRIS-WS jar
            %  file. If it does not exist, then it will try to access the
            %  latest jar over the internet.
            
            if ~exist('isSilent','var')
                isSilent=false;
            end
            success = false;
            %Check for required jar file for winston
            try
                % store the java class path for inspection
                jcp = javaclasspath('-all');
                
            catch er
                disp('Java is not enabled on this machine.  The Web Services Library will not work.');
                return
            end
            
            RequiredFiles = {'IRIS-WS-'};
            
            introuble = false;
            
            for FN = RequiredFiles
                if isempty(strfind([jcp{:}],FN{1}))
                    if ~isSilent
                        disp(['Missing ' FN{1}]);
                    end
                    introuble = true;
                end
            end
            
            if introuble
                if ~isSilent
                    disp('please add the IRIS-WS-latest.jar file to your javaclasspath');
                    disp('ex.  javaaddpath(''/usr/local/somewhere/IRIS-WS-latest.jar'');');
                end
                
                surrogate_jar = 'http://www.iris.edu/manuals/javawslibrary/download/IRIS-WS-latest.jar';
                
                [~,success] = urlread(surrogate_jar);%can we read the usgs.jar? if not don't bother to add it.
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % MODIFY FOLLOWING LINE TO POINT TO YOUR LOCAL
                % IRIS-WS-{something}.jar
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if success
                    javaaddpath(surrogate_jar);
                else
                    warning('irisFetch:noDefaultJar',...
                        'Unable to access the default jar.  Please add the the appropriate IRIS-WS jar file to your javaclasspath.');
                end
            end;
            success = true;
        end
        
        
        function channelList = flattenToChannel(networkTree)
            %irisFetch.flattenToChannel flattens the structure returned by irisFetch.Stations
            %
            %
            %USAGE
            %  flatStruct = irisFetch.flattenToChannel(networkTree)
            %
            %This takes the hierarchy returned by irisFetch.Stations, and
            %returns a 1xN array of channels.
            %
            % networkTree is a nested structure of the following format
            %   network.Station.Epoch.Channel.Epoch.[etc]
            %
            % flatStruct is an array containing ALL channel epochs, along
            % with unique identifying information from the parent levels,
            % such as networkTree.code, networkTree.station.code, etc.
            %
            %WORKING with the channelList
            %
            % Example: grabbing all BHZ channels from network IU
            %   myChannels = channelList({strcmp(channelList.NetworkCode},'IU') & ...
            %                   strcmp(channelList.ChannelCode, 'BHZ')
            %
            % Example: Find out how many IU networked stations were
            %          retrieved
            %   sum(strcmp({channelList.NetworkCode},'IU'))
            
            
            % moving from:  network -> station -> station epoch ->
            % channel -> channel epoch -> etc.
            
            % moving to: flat channels
            
            % hard-coded.
            %first, loop through and get rid of excess fields
            if ~isa(networkTree,'struct')
                error('Cannot Flatten a non-structure');
            end
            for eachNetwork = 1 : numel(networkTree)
                stationCodes = {networkTree(eachNetwork).Stations.Code};
                for eachStation = 1 : numel(networkTree(eachNetwork).Stations)
                    stationSites = {networkTree(eachNetwork).Stations(eachStation).Epochs.Site};
                    for eachStationEpoch = 1 : numel(networkTree(eachNetwork).Stations(eachStation).Epochs)
                        channelCodes = {networkTree(eachNetwork).Stations(eachStation).Epochs(eachStationEpoch).Channels.Code};
                        locationCodes = {networkTree(eachNetwork).Stations(eachStation).Epochs(eachStationEpoch).Channels.Location};
                        for eachChannel = 1 : numel(networkTree(eachNetwork).Stations(eachStation).Epochs(eachStationEpoch).Channels)
                            networkTree(eachNetwork).Stations(eachStation).Epochs(eachStationEpoch).Channels(eachChannel).Epochs = rmfield(networkTree(eachNetwork).Stations(eachStation).Epochs(eachStationEpoch).Channels(eachChannel).Epochs,'Class');
                            theseEpochs = networkTree(eachNetwork).Stations(eachStation).Epochs(eachStationEpoch).Channels(eachChannel).Epochs;
                            [theseEpochs.NetworkCode] = deal(networkTree(eachNetwork).Code);
                            [theseEpochs.NetworkDescription] = deal(networkTree(eachNetwork).Description);
                            [theseEpochs.StationCode] = deal(stationCodes{eachStation});
                            [theseEpochs.ChannelCode] = deal(channelCodes{eachChannel});
                            [theseEpochs.LocationCode]= deal(locationCodes{eachChannel});
                            [theseEpochs.Site] = deal(stationSites{eachStationEpoch});
                            networkTree(eachNetwork).Stations(eachStation).Epochs(eachStationEpoch).Channels(eachChannel).Epochs = theseEpochs;
                        end %eachChannel
                        networkTree(eachNetwork).Stations(eachStation).Epochs(eachStationEpoch).Channels = networkTree(eachNetwork).Stations(eachStation).Epochs(eachStationEpoch).Channels.Epochs;
                        % now, the structure is
                        % net -> sta -> epoch -> chan
                        networkTree(eachNetwork).Stations(eachStation).Epochs(eachStationEpoch).Channels = rmfield(networkTree(eachNetwork).Stations(eachStation).Epochs(eachStationEpoch).Channels,{'SensorChannelEpoch', 'SensorClass', 'SensitivityClass'});
                    end %eachStationEpoch
                    % fold the structure back even further...
                    % net -> sta -> chan
                    networkTree(eachNetwork).Stations(eachStation).Channels = deal([networkTree(eachNetwork).Stations(eachStation).Epochs.Channels]);
                end %eachStation
                networkTree(eachNetwork).Stations = rmfield(networkTree(eachNetwork).Stations,'Epochs');
                networkTree(eachNetwork).Channels = deal([networkTree(eachNetwork).Stations.Channels]);
            end %eachNetwork
            channelList = deal([networkTree.Channels]);
            %
        end
        
        
        
    end % static methods
    
    
    %%
    methods(Static, Access=protected)
        function mts = convertTraces(traces)
            %irisFetch.convertTraces converts traces from java to a matlab structure
            %USAGE:
            %  mts = convertTraces(traces) where TRACES a java trace
            %  class.
            
            
            for i = 1:length(traces)
                mt.network = char(traces(i).getNetwork());
                mt.station = char(traces(i).getStation());
                mt.location = char(traces(i).getLocation());
                mt.channel = char(traces(i).getChannel());
                
                mt.quality = char(traces(i).getQuality());
                
                mt.latitude = traces(i).getLatitude();
                mt.longitude = traces(i).getLongitude();
                mt.elevation = traces(i).getElevation();
                mt.depth = traces(i).getDepth();
                mt.azimuth = traces(i).getAzimuth();
                mt.dip = traces(i).getDip();
                
                mt.sensitivity = traces(i).getSensitivity();
                mt.sensitivityFrequency = traces(i).getSensitivityFrequency();
                
                mt.instrument = char(traces(i).getInstrument());
                mt.sensitivityUnits = char(traces(i).getSensitivityUnits());
                mt.data = traces(i).getData();
                
                mt.sampleCount = traces(i).getSampleCount();
                
                startDateString = char(traces(i).getStartTime().toString());
                endDateString = char(traces(i).getEndTime().toString());
                
                mt.startTime = datenum(startDateString, 'yyyy-mm-dd HH:MM:SS.FFF');
                mt.endTime = datenum(endDateString, 'yyyy-mm-dd HH:MM:SS.FFF');
                
                mts(i) = mt;
            end
        end
        
        
        %----------------------------------------------------------------
        % DATE conversion routines
        %----------------------------------------------------------------
        function javadate = mdate2jdate(matlabdate)
            %mdate2jdate converts a matlab date to a java Date class
            sdf = java.text.SimpleDateFormat('yyyy-MM-dd hh:mm:ss.SSS');
            sdf.setTimeZone(java.util.TimeZone.getTimeZone('GMT'));
            matlabDateString = datestr(matlabdate,'yyyy-mm-dd HH:MM:SS.FFF');
            javadate = sdf.parse(matlabDateString);
        end
        
        function matlabdate = jdate2mdate(javadate)
            %jdate2mdate converts a java Date class to a matlab datenum
            % WARNING: NO MILLISECONDS!!!!
            try
                matlabdate = datenum(...
                    javadate.getYear+1900,javadate.getMonth+1,javadate.getDate,...
                    javadate.getHours,javadate.getMinutes,javadate.getSeconds);
                matlabdate=datestr(matlabdate,31);
            catch
                matlabdate = [];
            end
        end
        
        %----------------------------------------------------------------
        % Look for GET / SET methods for the class.
        %----------------------------------------------------------------
        function [M, argType] = getSettableFields(obj)
            % strip the first 3 letters off the field ('set')
            M = irisFetch.getSetters(obj);
            for n=1:numel(M)
                M(n) = {M{n}(4:end)};
            end
        end
        
        
        function [M, argType] = getSetters(obj)
            [M, argType] = irisFetch.getMethods(obj,'set');
        end
        
        function [M, argType] = getGetters(obj)
            [M, argType] = irisFetch.getMethods(obj,'get');
        end
        
        function [methodList, argType] = getMethods(obj,searchPrefix)
            persistent className methodsAndArguments
            if isempty(className)
                className = {''};
                methodsAndArguments = {{''},{''}};
            end
            
            thisClass = class(obj);
            TF = strcmp(thisClass, className);
            if any(TF) % shortcut if this has been done before
                methodList = methodsAndArguments{TF,1};
                argType = methodsAndArguments{TF,2};
                return
            else
                loc = numel(className)+1;
            end
            
            
            argType = {}; %methodList = [];
            M = methods(obj);
            M2 = methods(obj,'-full');
            idx = false(size(M));
            for n=1:numel(M)
                if length(M{n}) < 3
                    continue
                end
                if strfind(M{n}(1:length(searchPrefix)),searchPrefix)
                    idx(n)=true;
                end
            end
            methodList = M(idx);
            
            %now get the argument class
            argList = M2(idx);
            p1=strfind(argList,'(');
            p2=strfind(argList,')');
            for n=1:numel(argList)
                argType(n) = {argList{n}(p1{n}+1:p2{n}-1)};
            end
            
            className(loc) = {thisClass};
            methodsAndArguments(loc,1) = {methodList};
            methodsAndArguments(loc,2) = {argType};
            
        end
        
        %----------------------------------------------------------------
        %
        %----------------------------------------------------------------
        function myStruct = parseCollection(obj)
            % ensure that obj is a java collection, such as
            % java.util.ArrayList
            itemCount  = obj.size;
            myStruct = [];
            
            prevClass = '';
            
            for n = 1 : itemCount
                thisObj = obj.get(n-1);
                %thisObj = obj.get(0);
                %obj.remove(0);
                
                myClass = class(thisObj);
                
                if ~strcmp(prevClass,myClass)
                    prevClass = myClass;
                    getterList = irisFetch.getGetters(thisObj);
                end
                
                %deal with special (known) collections
                switch myClass
                    case 'edu.iris.dmc.ws.station.model.ComplexNumber'
                        myStruct(n) = complex(...
                            double(thisObj.getReal),...
                            double(thisObj.getImaginary));
                        continue;
                    otherwise
                end
                
                for idx = 1 : numel(getterList)
                    thisMethod = getterList{idx};
                    if strcmpi(thisMethod,'get')
                        %this method is just "get"; skip it.
                        continue;
                    end
                    thisField = thisMethod(4:end); %get rid of 'get'
                    % if strcmpi(thisMethod,'getClass') || strcmpi(thisMethod,'getAny')
                    if strcmpi(thisMethod,'getAny')
                        continue
                    end
                    objectInField = thisObj.(thisMethod);
                    
                    thisClass = class(objectInField);
                    switch thisClass
                        case {'java.util.ArrayList','java.lang.String',...
                                'java.util.Date','java.lang.Double',...
                                'java.lang.Integer','java.math.BigInteger'}
                            myStruct(n).(thisField) = irisFetch.parseJavaClass(objectInField);
                            
                        case {'edu.iris.dmc.ws.station.model.Sensitivity', 'edu.iris.dmc.ws.station.model.Sensor'}
                            %prefix = irisFetch.getClassNameWithoutPackage(objectInField);
                            prefix = thisClass(find(thisClass=='.',1,'last')+1:end);
                            littleGetterList = irisFetch.getGetters(objectInField);
                            for idx2 = 1 : numel(littleGetterList)
                                thisLittleMethod = littleGetterList{idx2};
                                thisSubField = [prefix, thisLittleMethod(4:end)]; %get rid of 'get'
                                % if strcmpi(thisLittleMethod,'getClass')|| strcmpi(thisMethod,'getAny')
                                if strcmpi(thisMethod,'getAny')
                                    continue
                                end
                                myStruct(n).(thisSubField) = ...
                                    irisFetch.parseJavaClass(objectInField.(thisLittleMethod));
                            end
                            
                        otherwise
                            if ~isjava(objectInField)
                                myStruct(n).(thisField) = objectInField;
                            elseif isa(objectInField,'java.lang.Class')
                                myStruct(n).(thisField) = char(objectInField.toString);
                            else
                                % When this point is reached, we're looking
                                % at a pointer to the parent. if we were to
                                % parse it, we'd get a recursion error.
                            end
                    end
                    
                end
            end
            obj.clear;
        end
        
        function criteria = addCriteria(criteria, value, addMethod)
            %used to add Sta, Net, Loc, Chan to criteria
            % For example:
            %   singleAddToCriteria(criteria, {'OKCF','MNO?'},'addStation')
            % will invoke:
            %   criteria.addStation('OKCF').addStation('MNO?')
            if isempty(value)
                return %do nothing
            end
            if ~iscell(value)
                % it's probably a string
                criteria.(addMethod)(value);
            else
                % a cell may have multiple values, use 'em all.
                for n=1:numel(value)
                    criteria.(addMethod)(value{n})
                end
            end
        end
        
        
        function criteria = setBoxCoordinates(criteria, thisValue)
            % setBoxCoordinates (minLat, maxLat, minLon, maxLon)
            % values of 'NAN' are ignored
            if numel(thisValue) ~=4
                error('IRISFETCH:setBoxCoordinates:InvalidParameterCount',...
                    'Expected [minLat, maxLat, minLon, maxLon]');
            end
            setMethods = {'setMinimumLatitude','setMaximumLatitude',...
                'setMinimumLongitude','setMaximumLongitude'};
            for n=1:numel(setMethods)
                if ~isnan(thisValue(n))
                    criteria.(setMethods{n})(java.lang.Double(thisValue(n)));
                end
            end
        end
        
        
        function val = parseJavaClass(obj)
            % returns a matlab representative of some basic Java classes
            switch class(obj)
                case 'java.util.ArrayList'
                    val = irisFetch.parseCollection(obj);
                case 'java.lang.String'
                    val = char(obj);
                case 'java.util.Date'
                    val = irisFetch.jdate2mdate(obj);
                case {'java.lang.Double','java.lang.Integer',...
                        'java.math.BigInteger'}
                    val = double(obj);
                otherwise %add the java object in anyway.
                    val=obj;
            end
        end
        
        function criteria = setCriteria(criteria, paramList)
            
            %----------------------------------------------------------
            % The following code allows for open-ended search criteria
            % allowing it to change whenever the java library changes
            %
            % Instead of hard-coding each Setter, I query the java class to
            % find out its methods, then keep the ones that start with
            % "set".
            %
            % I also find out what the input parameters are for each, and
            % use that to properly create/cast the data.  Without doing
            % this, neither dates nor doubles would work. Boo.
            %----------------------------------------------------------
            
            % Get a list of parameters, their set functions, and input
            % types, and do it outside the loop so they are not needlessly
            % rerun
            
            [allSetMethods, argType] = irisFetch.getSetters(criteria);
            settableFieldnames = irisFetch.getSettableFields(criteria);
            allMethods = methods(criteria);
            
            while ~isempty(paramList) && numel(paramList) >= 2
                
                % Grab the parameter pair, then remove from parameter list
                thisParam = paramList{1};
                thisValue = paramList{2};
                paramList(1:2)=[];
                
                indexOfMethod = strcmpi(thisParam,settableFieldnames);
                if any(indexOfMethod)
                    
                    setMethod = allSetMethods{indexOfMethod};
                    switch argType{indexOfMethod}
                        
                        case 'java.util.Date'
                            criteria.(setMethod)(irisFetch.mdate2jdate(thisValue));
                            
                        case 'java.lang.Double'
                            criteria.(setMethod)(java.lang.Double(thisValue));
                            
                        otherwise
                            disp('Unanticipated argument type... trying');
                            criteria.(setMethod)(thisValue);
                    end
                    continue;
                end
                % we shall only pass this point if existing methods were not used
                
                switch lower(thisParam)
                    %handle special cases
                    case 'boxcoordinates'
                        criteria = irisFetch.setBoxCoordinates(criteria, thisValue);
                        
                    case 'radialcoordinates'
                        criteria.setLatitude(thisValue(1));
                        criteria.setLongitude(thisValue(2));
                        criteria.setMaximumRadius(thisValue(3));
                        if numel(thisValue) ==4 && ~isnan(thisValue(4))
                            criteria.setMinimumRadius(thisValue(4));
                        end
                        
                    otherwise
                        % this will blow up if java doesn't recongize
                        % thisValue
                        criteria.(allMethods{strcmpi(allMethods,thisParam)})(thisValue);
                end
            end
        end
        
        function s = getClassNameWithoutPackage(obj)
            s=char(obj.getClass.toString);
            s(1:find(s=='.',1,'last')) = [];
        end
    end %static protected methods
end

