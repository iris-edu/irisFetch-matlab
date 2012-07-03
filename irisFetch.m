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
   %           '2010-02-27 06:30:00','2010-02-27 10:30:00')
   %
   %
   %  IRISFETCH.Stations() retrieves station metadata from the IRIS-DMC.
   %  This data may be retrieved at a variety of detail levels. From
   %  broadest to most specific, these are NETWORK, STATION, CHANNEL, and
   %  RESPONSE.
   %
   %    Example:
   %      s = IRISFETCH.Stations('channel','*','ANMO','10','BH?')
   %      % the above returns a tree-like structure, now make it more manageable
   %      s = IRISFETCH.flattenToChannel(s); % convert a 1xN array of channel epochs.
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
   %IRISFETCH.FLATTENTOCHANNEL IRISFETCH.FLATTENTOSTATION
   %IRISFETCH.TESTCOMPATIBILITY IRISFETCH.VERSION JAVAADDPATH
   
   % Celso Reyes, Rich Karstens
   % IRIS-DMC
   % February 2012
   
   % 2012 June 18 r1.3.2
   % spelling fix in parse, and initialized  
   %
   % 2012 June 14 r1.3.1
   % Fixed problem where Traces.sacpz.units was not converted from java
   % strings
   %
   % 2012 June 8, r1.3.0
   % Changed Traces routine to match IRIS-WS v1.5's new conventions.
   %  - This change should be transparent to matlab users (except new
   %    version of library must be downloaded)
   %  - verbosity set through routine, not as parameter
   %  - fetch sacpz is required
   % added sacpz subfield: units.
   % fixed bug in Event where radius didn't convert from double to
   % java.lang.Double, resulting in error.
   % added the App name to Traces.
   % added flattenToStation routine
   %
   % Modified the Traces routine: Added authorization ability and ability to get
   % poles and zeros responses.  Reworked the parameter list to be more
   % flexible. Completely overhauled the java->matlab parsing for accuracy and speed.
   %
   % 2012 Mar 8, Fixed problem where empty traces created a "data not
   % assigned" style message. Removed java objects from the returned
   % structures... these have already been parsed into fields.
   %
   % 2012 Feb 29, Fixed date parsing issues, allowing accuracy to
   % millisecond level. Also fixed problem where filter numerators not
   % showing up properly.
   %
   % 2012 Feb 27, 1.1.1 fixed problem where channel epochs were missing from
   % flattened stations, and added ability to translate into both Long and
   % Integer java types.
   %
   % 2012 Feb 13, added SampleRate to trace structure, improved error
   % catching.
   %
   % 2012 Feb 9, minor documentation update
   
   
   properties
   end %properties
   
   methods(Static)
      function v = version()
         % return the version number of irisFetch
         v = '1.3.2';
      end
      
  
      function ts = Traces(network, station, location, channel, startDateStr, endDateStr, varargin )
         %irisFetch.Traces accesses waveform with associated channel metadata
         %
         %USAGE
         % tr = irisFetch.Traces(network, station, location,
         % channel, startDate, endDate) will use channel and date
         % criteria to retrieve one or more seismic traces, which are
         % stored as structures containing typical SAC-equivalent
         % metadata.
         %
         % startDate and endDate must be formatted thusly:
         %      'YYYY-MM-DD hh:mm:ss' or 'YYYY-MM-DD hh:mm:ss.sss'
         %
         % network, station, location, and channel all accept '*' and
         % '?' wildcards, as well as comma separated lists.
         %
         % tr = irisFetch.Traces(..., quality)
         % allows you to specify a quality factor, such as 'B' or 'M'.
         % If not specified, quality defaults to 'B'
         %
         % tr = irisFetch.Traces(..., 'includePZ')
         % will retrieve poles and zeroes as well.  These are the same
         % poles/zeros as found from http://www.iris.edu/ws/sacpz/
         %
         % tr = irisFetch.Traces(..., 'verbose')
         % provides additional debugging information
         %
         % tr = irisFetch.Traces(..., usernameAndPassword )
         % allows authorized users to access restricted data.
         % usernameAndPassword must be a cell containing the username and
         % password.  
         %    Sample:
         %      unamepwd = {'myname@iris.edu', 'mypassword'}
         %
         %
         %
         %EXAMPLES:
         %
         %    Example 1:
         %      % get 3 components for a 4-minute time period,
         %      % using the '?' wildcard and also retrieving response info.
         %      ts = irisFetch.Traces('IU','ANMO','10','BH?',...
         %           '2010-02-27 06:30:00','2010-02-27 10:30:00','includePZ')
         %
         %    Example 2:
         %      % get z-channels for a comma-separated list of stations
         %      % while specifying a quality of 'B'
         %      ts = irisFetch.Traces('IU','ANMO,ANTO,YSS','00','BHZ',...
         %           '2010-02-27 06:30:00','2010-02-27 10:30:00', 'B')
         %
         %    Example 3:
         %      % get z-data for all BHZ stations that belong to the IU
         %      % network and have a location code of '00'
         %      ts = irisFetch.Traces('IU','*','00','BHZ',...
         %           '2010-02-27 06:30:00','2010-02-27 10:30:00')
         %
         %ABOUT THE RETURNED TRACE
         %  The returned trace(s) will be an array of structures, with
         %  each structure containing fields full of information. For
         %  example, if I retrieve N traces, I will have the following
         %  structure:
         %
         %  1xN struct array with fields:
         %     network
         %     station
         %     location
         %     channel
         %     quality
         %     latitude
         %     longitude
         %     elevation
         %     depth
         %     azimuth
         %     dip
         %     sensitivity
         %     sensitivityFrequency
         %     instrument
         %     sensitivityUnits
         %     data
         %     sampleCount
         %     sampleRate
         %     startTime
         %     endTime
         %
         %   COMMON MANIPULATIONS
         %     To access the date as text, use datestr().  For trace(s)
         %     ts the full date with milliseconds can be seen using:
         %     datestr(ts, ,'YYYY-MM-DD hh:mm:ss.FFF')
         %
         %     To scale the data, divide each trace's data by its
         %     sensitivity.  The resulting values are in
         %     sensitivityUnits.
         %
         % SEE ALSO datestr
         
         getsacpz = false;
         verbosity = false;
         authorize = false;
         quality = 'B';
         
         % deal with a variety of possible parameters
         for n=1:numel(varargin)
            thisParameter = varargin{n};
            switch class(thisParameter)
               case 'cell'
                  % parameter is the username/pwd combo
                  assert(numel(thisParameter)==2 &&...
                     ischar(thisParameter{1}) && ...
                     ischar(thisParameter{2}),...
                     ['A cell passed as an optional parameter is assumed',...
                     ' to contain credentials. ',...
                     'eg. {''myname'',''mypassword''}.']);
                  username = thisParameter{1};
                  userpwd = thisParameter{2};
                  authorize = true;
               case 'char'
                  switch upper(thisParameter)
                     case {'D','R','Q','M','B'}
                        quality = thisParameter;
                     case {'INCLUDEPZ'}
                        getsacpz = true;
                     case {'VERBOSE'}
                        verbosity = true;
                     otherwise
                        error('IRISFETCH:Trace:unrecognizedParameter',...
                           'The text you included as an optional parameter did not parse to either a qualitytype (D,R,Q,M,B) or ''INCLUDEPZ'' or ''VERBOSE'''); 
                  end
               case 'logical'
                  verbosity = thisParameter;
                  % old usage, may be deprecated in the future.
               otherwise
                  error('IRISFETCH:Trace:unrecognizedParameter',...
                           'The optional parameter wasn''t recognized. %s', class(thisParameter)); 
            end
         end
                  
         if isnumeric(startDateStr)
            startDateStr = datestr(startDateStr,'yyyy-mm-dd HH:MM:SS' );
         end
         
         if isnumeric(endDateStr)
            endDateStr = datestr(endDateStr,'yyyy-mm-dd HH:MM:SS' );
         end
            
         location = strrep(location,' ','-');
         
         try
           tracedata = edu.iris.dmc.ws.extensions.fetch.TraceData();
         catch er
            switch er.identifier
               case 'MATLAB:undefinedVarOrClass'
                  
                  isSilent = true; %suppress messages within the connectTo...
                  errtext=['The Web Services library does not appear to be in the javaclasspath.\n',...
                     'Please download the latest version from \n',...
                     'http://www.iris.edu/manuals/javawslibrary/#download\n ',...
                     'and then add it to your classpath. \n'];
                  warning('IRISFETCH:NoIrisWSJarInstalled',errtext);
                  success = irisFetch.connectTo_IRIS_WS_jar(isSilent);
                  if ~success
                     error('IRISFETCH:Traces:UnableToInstallIrisWSJar',...
                        'irisFetch was unable to recover, please download and add the latest IRIS-WS-JAR to your javaclasspath');
                  end
                     disp('irisFetch.connectTo_IRIS_WS_jar() has was able to connect you to the appropriate java library. Continuing...');
                  tracedata = edu.iris.dmc.ws.extensions.fetch.TraceData();
               otherwise
                  rethrow(er)
                     
            end
         end
         tracedata.setAppName(['MATLAB:irisFetch/' irisFetch.version()]);
         
         try
            % only library 1.5 and greater will successfully do this
            tracedata.setVerbosity(verbosity);
         catch er
            % if the error is due to a bad library version, then recommend
            % updating the library.
            
            %otherwise
            rethrow(er);
         end
         
         
         try
            if authorize
               traces = tracedata.fetchTraces(network, station, location, channel, startDateStr, endDateStr, quality, getsacpz, username, userpwd);
            else % not authorizing
               traces = tracedata.fetchTraces(network, station, location, channel, startDateStr, endDateStr, quality, getsacpz);
            end
            
            ts = irisFetch.convertTraces(traces);
            clear traces;
         catch je
            switch je.identifier
               case 'MATLAB:undefinedVarOrClass'
                  errtext=['The Web Services library does not appear to be in the javaclasspath.\n',...
                     'Please download the latest version from \n',...
                     'http://www.iris.edu/manuals/javawslibrary/#download\n ',...
                     'and then add it to your classpath. \n'];
                  warning('IRISFETCH:NoIrisWSJarInstalled',errtext);
                  isSilent = true; %suppress messages within the connectTo...
                  success = irisFetch.connectTo_IRIS_WS_jar(isSilent);
                  if success
                     disp('irisFetch.connectTo_IRIS_WS_jar() has was able to connect you to the appropriate java library. Continuing...');
                     try
                        if authorize
                           traces = tracedata.fetchTraces(network, station, location, channel, startDateStr, endDateStr, quality, getsacpz, username, userpwd);
                        else % not authorizing
                           traces = tracedata.fetchTraces(network, station, location, channel, startDateStr, endDateStr, quality, getsacpz);
                        end
                        ts = irisFetch.convertTraces(traces);
                        clear traces;
                        
                     catch je2
                        rethrow(je2)
                     end
                  else
                     error('IRISFETCH:Traces:UnableToInstallIrisWSJar',...
                        'irisFetch was unable to recover, please download and add the latest IRIS-WS-JAR to your javaclasspath');
                  end
               case 'MATLAB:Java:GenericException'
                  if any(strfind(je.message,'URLNotFoundException'))
                     % we got a 404 from somewhere. (based on ice.net)
                     error('IRISFETCH:Trace:URLNotFoundException',...
                        'Trace found no requested data. Instead, it ran in to the following error:\n%s',...
                        je.message);
                  else
                     rethrow(je)
                  end
               otherwise
                  fprintf('Exception occured in IRIS Web Services Library: %s\n', je.message);
                  rethrow(je)
                  
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
         import edu.iris.dmc.ws.station.model.*
         
         if nargin==1 && strcmpi(detailLevel,'help')
            disp('HELP request recognized, but not implemented');
            return
         elseif nargin < 5
            error('not enough arguments.%d',nargin);
         end
         
         try
            outputLevel = ws.criteria.OutputLevel.(upper(detailLevel));
         catch je
            switch je.identifier
               case 'MATLAB:undefinedVarOrClass'
                  error('IRISFETCH:NoIrisWSJarInstalled',...
                     'The necessary IRIS-WS java library was not recognized or found. Please ensure it is on your javaclasspath');
               case 'MATLAB:subscripting:classHasNoPropertyOrMethod'
                  error('IRISFETCH:invalidOutputLevel',...
                     'The selected outputLevel [''%s''] was not recognized.',...
                     upper(detailLevel));
               otherwise
                  rethrow(je);
            end
         end
         
         
         indexOffsetOfBASEURL=find(strcmpi(varargin(1:2:end),'BASEURL'),1,'first') * 2;
         try
            baseURL = varargin{indexOffsetOfBASEURL};
         catch
            % don't do anything
         end
         
         
         
         serviceManager = ws.service.ServiceUtil.getInstance();
         serviceManager.setAppName(['MATLAB:irisFetch/' irisFetch.version()]);
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
         
         try
            j_networks = service.fetch(criteria, outputLevel);
         catch je
            if strfind(je.message,'ServiceNotSupportedException')
               error('IRISFETCH:ServiceNotSupportedByLibrary',...
                  'The IRIS-WS java library version doesn''t support the requested station service version');
            else
               rethrow(je)
            end
         end
         
         networkStructure = irisFetch.parse(j_networks);
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
         %  'MinimumDepth'
         %  'MaximumDepth'
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
         serviceManager.setAppName(['MATLAB:irisFetch/' irisFetch.version()]);
         
         
         indexOffsetOfBASEURL=find(strcmpi(varargin(1:2:end),'BASEURL'),1,'first') * 2;
         
         try
            baseURL = varargin{indexOffsetOfBASEURL};
         catch
            % don't do anything
         end
         
         if exist('baseURL','var')
            varargin(indexOffsetOfBASEURL-1:indexOffsetOfBASEURL) = [];
            service = serviceManager.getEventService(baseURL);
         else
            service = serviceManager.getEventService();
         end
         
         criteria = ws.criteria.EventCriteria;
         criteria = irisFetch.setCriteria(criteria, varargin);
         if nargout == 2
            urlParams = criteria.toUrlParams;
         end
         disp('fetching from IRIS-DMC')
         j_events = service.fetch(criteria);
         disp('parsing into MATLAB structures')
         %tic;events = irisFetch.parseCollection(j_events);toc
         %disp(events)
         events = irisFetch.parse(j_events);
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
         minimumJarVersion = '1.5';
         
         
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
               disp('ex.  javaaddpath(''/usr/local/somewhere/IRIS-WS.jar'');');
            end            
            
         surrogate_jar = 'http://www.iris.edu/manuals/javawslibrary/matlab/IRIS-WS-1.5-matlab.jar';
         
            [~,success] = urlread(surrogate_jar);%can we read the .jar? if not don't bother to add it.
            if success
               javaaddpath(surrogate_jar);
            else
               warning('irisFetch:noDefaultJar',...
                  'Unable to access the default jar.  Please download and add the latest IRIS-WS-JAR to your javaclasspath.');
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
         %returns a 1xN array of channels  (channel epochs, technically).
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
         
         makeAlphabetical = false; %leave fields in alphabetical order, or reorder according to common-sense
         
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
                    % networkTree(eachNetwork).Stations(eachStation).Epochs(eachStationEpoch).Channels(eachChannel).Epochs = rmfield(networkTree(eachNetwork).Stations(eachStation).Epochs(eachStationEpoch).Channels(eachChannel).Epochs,'Class');
                     theseEpochs = networkTree(eachNetwork).Stations(eachStation).Epochs(eachStationEpoch).Channels(eachChannel).Epochs;
                     [theseEpochs.NetworkCode] = deal(networkTree(eachNetwork).Code);
                     [theseEpochs.NetworkDescription] = deal(networkTree(eachNetwork).Description);
                     [theseEpochs.StationCode] = deal(stationCodes{eachStation});
                     [theseEpochs.ChannelCode] = deal(channelCodes{eachChannel});
                     [theseEpochs.LocationCode]= deal(locationCodes{eachChannel});
                     [theseEpochs.Site] = deal(stationSites{eachStationEpoch});
                     networkTree(eachNetwork).Stations(eachStation).Epochs(eachStationEpoch).Channels(eachChannel).Epochs = theseEpochs;
                  end %eachChannel
                  networkTree(eachNetwork).Stations(eachStation).Epochs(eachStationEpoch).Channels = [networkTree(eachNetwork).Stations(eachStation).Epochs(eachStationEpoch).Channels.Epochs];
                  % now, the structure is
                  % net -> sta -> epoch -> chan
               end %eachStationEpoch
               % fold the structure back even further...
               % net -> sta -> chan
               networkTree(eachNetwork).Stations(eachStation).Channels = deal([networkTree(eachNetwork).Stations(eachStation).Epochs.Channels]);
            end %eachStation
            networkTree(eachNetwork).Stations = rmfield(networkTree(eachNetwork).Stations,'Epochs');
            networkTree(eachNetwork).Channels = deal([networkTree(eachNetwork).Stations.Channels]);
         end %eachNetwork
         channelList = deal([networkTree.Channels]);
         
         % one more thing... now bring Sensitivity and Sensor up, since
         % they are 1x1 structs
         for n=1:numel(channelList)
            if isstruct(channelList(n).Sensitivity)
            sens=channelList(n).Sensitivity;
            fn = fieldnames(sens);
            for m=1:numel(fn);
               channelList(n).(fn{m}) = sens.(fn{m});
            end       
            else
               % Sensitivy wasn't included!
            end
            if isstruct(channelList(n).Sensor)
            sensor=channelList(n).Sensor;
            fn = fieldnames(sensor);
            for m=1:numel(fn);
               channelList(n).(fn{m}) = sensor.(fn{m});
            end
            else
               % Sensor wasn't included!
            end
         end
         channelList = rmfield(channelList,{'Sensitivity','Sensor'});
         %         
         
         if ~makeAlphabetical
            % now, reorder to make it visually coherent.
            descriptorstuff={'NetworkCode';'StationCode';'LocationCode';'ChannelCode';'NetworkDescription';'Site'};
            positionalstuff={'Latitude';'Longitude';'Elevation';'Depth';'Azimuth';'Dip'};
            otherstuff={'SampleRate';'StartDate';'EndDate'};
            fieldsattop=[descriptorstuff; positionalstuff; otherstuff];
            
            fn = fieldnames(channelList);
            
            fieldsattop = fieldsattop(ismember(fieldsattop,fn)); %ensure fields exist
            
            for n=1:numel(fieldsattop);
               fn(strcmp(fn,fieldsattop(n))) = [];
            end
            neworder = [fieldsattop; fn];
            channelList = orderfields(channelList, neworder);
         end
      end
      
      function flatStruct = flattenToStation(networkTree)
         %irisFetch.flattenToStation flattens the structure returned by irisFetch.Stations
         %
         %
         %USAGE
         %  flatStruct = irisFetch.flattenToStation(networkTree)
         %
         %This takes the hierarchy returned by irisFetch.Stations, and
         %returns a 1xN array of stations (station epochs, technically).
         %
         % networkTree is a nested structure of the following format
         %   network.Station.Epoch.[etc]
         %
         % flatStruct is an array containing ALL station epochs, along
         % with unique identifying information from the parent levels,
         % such as networkTree.code, networkTree.station.code, etc.
         %
         % If the 
         
         if isempty(networkTree)
            flatStruct = networkTree;
            return
         end
         if isempty([networkTree.Stations])
            flatStruct = networkTree;
            return
         end
         
         makeAlphabetical = false; %leave fields in alphabetical order, or reorder according to common-sense
         
         for netIdx = 1: numel(networkTree)
            % add the Network code to the Stations
            [networkTree(netIdx).Stations.NetworkCode] = deal(networkTree(netIdx).Code);
            [networkTree(netIdx).Stations.NetworkDescription] = deal(networkTree(netIdx).Description);
         end
         flatStruct = [networkTree.Stations];
         for staIdx = 1 : numel(flatStruct);
            [flatStruct(staIdx).Epochs.StationCode] = deal(flatStruct(staIdx).Code);
            [flatStruct(staIdx).Epochs.NetworkCode] = deal(flatStruct(staIdx).NetworkCode);
            [flatStruct(staIdx).Epochs.NetworkDescription] = deal(flatStruct(staIdx).NetworkDescription);
         end
         flatStruct = [flatStruct.Epochs];
         if isempty([flatStruct.Channels]); flatStruct = rmfield(flatStruct,'Channels'); end
         
         if ~makeAlphabetical
            % now, reorder to make it visually coherent.
            descriptorstuff={'NetworkCode';'StationCode';'LocationCode';'ChannelCode';'NetworkDescription';'Site'};
            positionalstuff={'Latitude';'Longitude';'Elevation';'Depth';'Azimuth';'Dip'};
            otherstuff={'SampleRate';'StartDate';'EndDate'};
            fieldsattop=[descriptorstuff; positionalstuff; otherstuff];
            
            fn = fieldnames(flatStruct);
            
            fieldsattop = fieldsattop(ismember(fieldsattop,fn)); %ensure fields exist
            
            for n=1:numel(fieldsattop);
               fn(strcmp(fn,fieldsattop(n))) = [];
            end
            neworder = [fieldsattop; fn];
            flatStruct = orderfields(flatStruct, neworder);
         end
         
      end
      
   end % static methods
   
   
   
   %%
   methods(Static, Access=protected)
      function d = jArrayList2complex(jArrayList)
         % for use on ArrayList objects containing things with getReal() and getImaginary()
         %  edu.iris.dmc.ws.sacpz.model.Pole
         %  edu.iris.dmc.ws.sacpz.model.Zero
         %  edu.iris.dmc.ws.station.model.ComplexNumber
         r= zeros(jArrayList.size(),2);
         for n=1:jArrayList.size()
            r(n,:)=double([jArrayList.get(n-1).getReal(), jArrayList.get(n-1).getImaginary()]);
         end
         
         if any(r(:,2))
            d = complex(r(:,1),r(:,2));
         else
            d = r(:,1);
         end
      end
      
      function mts = convertTraces(traces)
         %irisFetch.convertTraces converts traces from java to a matlab structure
         %USAGE:
         %  mts = convertTraces(traces) where TRACES a java trace
         %  class.
         
         blankSacPZ = struct('units','','constant',[],'poles',[],'zeros',[]);
         
         blankTrace = struct('network','','station','','location',''...
            ,'channel','','quality','',...
            'latitude',0,'longitude',0,'elevation',0,'depth',0,...
            'azimuth',0,'dip',0,...
            'sensitivity',0,'sensitivityFrequency',0,...
            'instrument','','sensitivityUnits','UNK',...
            'data',[],'sampleCount',0,'sampleRate',nan,...
            'startTime',0,'endTime',0,'sacpz',blankSacPZ);
         mts=blankTrace;
         for i = 1:length(traces)
            mt=blankTrace;
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
            mt.sampleRate = traces(i).getSampleRate();
            
            startDateString = char(traces(i).getStartTime().toString());
            endDateString = char(traces(i).getEndTime().toString());
            
            mt.startTime = datenum(startDateString, 'yyyy-mm-dd HH:MM:SS.FFF');
            mt.endTime = datenum(endDateString, 'yyyy-mm-dd HH:MM:SS.FFF');
            
            try
               jsacpz = traces(i).getSacpz();
            catch er
               if strcmp(er.identifier,'MATLAB:noSuchMethodOrField')
                  warning('IRISFETCH:convertTraces:noGetSacPZmethod',...
                     'probably using older verision of the ws-library. please retrieve the latest version');
                  jsacpz = [];
               else
                  rethrow(er)
               end
            end
            if ~isempty(jsacpz)
               sacpz.units = char(traces(i).getSacpz().getInputUnit());
               sacpz.constant = traces(i).getSacpz().getConstant();
               sacpz.poles= irisFetch.jArrayList2complex(traces(i).getSacpz().getPoles());
               sacpz.zeros= irisFetch.jArrayList2complex(traces(i).getSacpz().getZeros());
               mt.sacpz = sacpz;
            end
            mts(i) = mt;
         end
      end
      
      
      %----------------------------------------------------------------
      % DATE conversion routines
      % 
      % 1970-01-01 is datenum 719529; there are 86400000 ms in a day.
      %----------------------------------------------------------------
      function javadate = mdate2jdate(matlabdate)
         %mdate2jdate converts a matlab date to a java Date class
         % TRUNCATES TO Milliseconds
         javadate = java.util.Date;
         javadate.setTime((datenum(matlabdate) - 719529) * 86400000);
      end
      
      function matlabdate = jdate2mdate(javadate)
         % jdate2mdate converts a java Date class to a matlab datenum
         try 
            matlabdate = 719529 + javadate.getTime() / 86400000;
            matlabdate=datestr(matlabdate,'yyyy-mm-dd HH:MM:SS.FFF');
         catch je
            warning(je)
            matlabdate = [];
         end
      end
      
      %----------------------------------------------------------------
      % Look for GET / SET methods for the class.
      %----------------------------------------------------------------
      function M = getSettableFields(obj)
         % strip the first 3 letters off the field ('set')
         M = irisFetch.getSetters(obj);
         for n=1:numel(M)
            M(n) = {M{n}(4:end)};
         end
      end
            
      function [M, argType] = getSetters(obj)
         [M, argType] = irisFetch.getMethods(obj,'set');
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
         idx = strncmp(searchPrefix,M, length(searchPrefix));
         methodList = M(idx);
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
 %%     
      %================================================================
      %----------------------------------------------------------------
      % BEGIN: PARSING ROUTINES
      %----------------------------------------------------------------
 
       function [getterList, fieldList] = getMethodsAndFields(obj)
          
          
         % this function uses a cache to speed up the retrieval of
         % get_methods and fieldnames.
         
         persistent className 
         persistent methodL 
         persistent fieldL
         
         if isempty(className)
            className = {''};
            methodL = {''};
            fieldL = {''};
         end
         
         thisClass = class(obj);
         
         TF = strcmp(thisClass, className);
         
         if any(TF) % shortcut if this has been done before
            getterList = methodL{TF};
            fieldList = fieldL{TF};
            return
         else
            loc = numel(className)+1;
         end

         allMethods = methods(obj);
         getterList = allMethods(strncmp('get',allMethods, 3));
         
         % filter classes need to have class names for them to make sense
         % to the user.
         if isa(obj,'edu.iris.dmc.ws.station.model.Filter')            
            getterList = getterList(...
               ~( strcmp('get',getterList) | ...
               strcmp('getAny',getterList) ));
            
         else
            getterList = getterList(...
               ~( strcmp('get',getterList) | ...
               strcmp('getClass',getterList) | ...
               strcmp('getAny',getterList) ));
         end
         % eliminate recursions
         switch thisClass
            case 'edu.iris.dmc.ws.station.model.Station'
               n = strcmp(getterList,'getNetwork');
            case 'edu.iris.dmc.ws.station.model.StationEpoch'
               n = strcmp(getterList,'getStation');
            case 'edu.iris.dmc.ws.station.model.Channel'
               n = strcmp(getterList,'getStationEpoch');
            case 'edu.iris.dmc.ws.station.model.ChannelEpoch'
               n = strcmp(getterList,'getChannel');
            case 'edu.iris.dmc.ws.station.model.Response'
               n = strcmp(getterList,'getChannelEpoch');
            case 'edu.iris.dmc.ws.station.model.Sensor'
               n = strcmp(getterList,'getChannelEpoch');
            otherwise
               n=[];
         end
         getterList(n) = [];
         
         fieldList = strrep(getterList,'get',''); %get rid of 'get'
         
         className(loc) = {thisClass};
         methodL(loc) = {getterList};
         fieldL(loc) = {fieldList};
      end
     
      function myStruct = parseObjectViaGetMethods(thisObj)
         % This routine should only be called for objects. Not for arrays
         % NOTE: assumes a single/scalar object.
         myStruct=[];
         [getterList, fieldnameList] = irisFetch.getMethodsAndFields(thisObj);
         
         for idx = 1 : numel(getterList)
            value = thisObj.(getterList{idx});
            if isjava(value)
               myStruct.(fieldnameList{idx}) = irisFetch.parse(value);
            else
               myStruct.(fieldnameList{idx}) = value;
            end
         end
         %
         % disp(myStruct);
         %
         
      end %function_parseObjectViaGetMethods
      
      
      function myGuts = parse(obj)
         % parse this object, first dealing with java built-in classes,
         % then proceeding to the iris classes
         
         myClass = class(obj);
         firstbit = myClass(1:find(myClass == '.',1,'first')-1);
         %switch strtok(myClass,'.')
         switch firstbit
            case 'java' % deal with a built-in java class
               switch myClass
                  case {'java.lang.String'}
                     myGuts = char(obj);
                  case 'java.lang.Class'
                     myGuts = char(obj.getCanonicalName);
                     
                  case {'java.lang.Double',...
                        'java.lang.Integer',...
                        'java.math.BigInteger',...
                        'java.lang.Long'}
                     myGuts = obj.doubleValue;
                     
                  case 'java.util.ArrayList' % this was an array of arrays
                     % myGuts = irisFetch.parseCollection(obj);
                     % myGuts = irisFetch.parseArrayList(obj, level);
                     
                     if obj.isEmpty();
                        myGuts = [];
                     else
                        
                        
                        switch class(obj.get(0))
                           
                           case {'edu.iris.dmc.ws.sacpz.model.Pole', ...
                                 'edu.iris.dmc.ws.sacpz.model.Zero', ...
                                 'edu.iris.dmc.ws.station.model.ComplexNumber'}
                              myGuts = irisFetch.jArrayList2complex(obj);
                              
                           case {'java.lang.Double',...
                                 'java.lang.Integer',...
                                 'java.math.BigInteger',...
                                 'java.lang.Long',...
                                 'double'}
                              myGuts = double(obj.toArray(javaArray('java.lang.Double',obj.size())));
                              
                           otherwise
                              for n = obj.size:-1:1
                                 mG = irisFetch.parse(obj.get(n-1));
                                 try
                                    myGuts(n) = mG;
                                 catch er
                                    if strcmp(er.identifier, 'MATLAB:heterogenousStrucAssignment')
                                       f = fieldnames(mG);
                                       for z=1:numel(f)
                                          myGuts(n).(f{z}) = mG.(f{z});
                                       end
                                    else
                                       rethrow (er)
                                    end
                                 end
                                 
                              end
                        end
                     end %endif obj is empty
                     
                  case 'java.util.Date' % 1970-01-01 is datenum 719529; 86400000 ms in a day.
                     myGuts = irisFetch.jdate2mdate(obj);
                     %matlabdate = 719529 + obj.getTime() / 86400000;
                     %myGuts=datestr(matlabdate,'yyyy-mm-dd HH:MM:SS.FFF');
                                          
                  otherwise
                     disp(['not sure how to deal with JAVA class :', myClass]);
               end
               
            case 'edu' % deal with one of our classes
               switch myClass
                  case {'edu.iris.dmc.ws.station.model.Sensitivity','edu.iris.dmc.ws.station.model.Sensor'}
                     
                     % in versions of irisFetch prior to 1.2, these were 
                     % add these particular classes to the existing structure
                     % WITHOUT nesting to another level. Instead, use the
                     % class name as the prefix.
                     
                     % beware of recursions!!!
                     [getterList, fieldnameList] = irisFetch.getMethodsAndFields(obj);
                     prefix = myClass( find( myClass=='.', 1, 'last') + 1:end);         
                     fieldnameList = strcat(prefix, fieldnameList);
                     
                     for idx = 1 : numel(getterList)
                        myGuts.(fieldnameList{idx}) = irisFetch.parse(obj.(getterList{idx}));
                     end
                     
                  otherwise
                     myGuts = irisFetch.parseObjectViaGetMethods(obj);
               end
            otherwise
               % take best guess
               myGuts = irisFetch.parseObjectViaGetMethods(obj);
         end
               
      end % fuction_parse
      
      function myguts = parseArrayList(obj)
         % it is known that obj is of class java.util.ArrayList
         if obj.isEmpty(); myguts = []; return; end;
         switch class(obj.get(0))
            
            case {'edu.iris.dmc.ws.sacpz.model.Pole', ...
                  'edu.iris.dmc.ws.sacpz.model.Zero', ...
                  'edu.iris.dmc.ws.station.model.ComplexNumber'}
               myguts = irisFetch.jArrayList2complex(obj);
               
            case {'java.lang.Double',...
                  'java.lang.Integer',...
                  'java.math.BigInteger',...
                  'java.lang.Long',...
                  'double'}
               myguts = double(obj.toArray(javaArray('java.lang.Double',obj.size())));
               
            otherwise
               for n = obj.size:-1:1
                  myguts(n) = irisFetch.parse(obj.get(n-1));
               end
         end
         
      end
      
      %----------------------------------------------------------------
      % END: PARSING ROUTINES
      %----------------------------------------------------------------
      %================================================================
      %%

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
                     criteria.toUrlParams;
                  case 'java.lang.Double'
                     criteria.(setMethod)(java.lang.Double(thisValue));
                     
                  case 'java.lang.Long'
                     criteria.(setMethod)(java.lang.Long(thisValue));
                     
                  case 'java.lang.Integer'
                     criteria.(setMethod)(java.lang.Integer(thisValue));
                     
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
                  criteria.setLatitude(java.lang.Double(thisValue(1)));
                  criteria.setLongitude(java.lang.Double(thisValue(2)));
                  criteria.setMaximumRadius(java.lang.Double(thisValue(3)));
                  if numel(thisValue) ==4 && ~isnan(thisValue(4))
                     criteria.setMinimumRadius(java.lang.Double(thisValue(4)));
                  end
                  
               otherwise
                  % this will blow up if java doesn't recongize
                  % thisValue
                  criteria.(allMethods{strcmpi(allMethods,thisParam)})(thisValue);
            end
         end
      end
      
      %--------------------------------------------------------------------
   
   end %static protected methods
end

