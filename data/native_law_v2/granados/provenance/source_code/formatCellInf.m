function formInf = formatCellInf(cellInf, varargin)
% FORMATCELLINF Format cellInf for analysis and apply selected filters
%
%   formInf = formatCellInf(cellInf) returns a struct where each field
%   contains a separate channel of cellInf. Channel names default to 
%   'Channel 1', ... 'Channel n'. 'cellInf' can be a struct array or a
%   pre-formatted formInf, in which case no re-structuring is applied. In
%   any case, any sparse fields in the cellInf are converted to full fields
%   with zero values set to NaN.
%   
%   formInf = filterCellInf(cellInf, ..., 'Param1', val1, ...) enables the
%   specification of optional parameter name/value pairs. Parameters are:
%       'ChannelNames' -- a cellstr of channel names to use for each slice
%       of an array-like cellInf (will not be applied if the cellInf is
%       already formatted in a struct with channel-name fields).
%
%       'Fields' -- a cellstr of fields to keep, or a struct of cellstr 
%       specifying the fields to keep for each channel.
%
%       'CellFilter' -- a logical index specifying cells to keep.
%
%       'TimeFilter' -- a logical index specifying time points to keep.
%
%       'Origin' -- the time point considered to be the origin. This will
%       be included as a parameter of the formInf structure, and used to
%       offset the values in the times array if there is one.

% Check input
ip = inputParser;
ip.addRequired('cellInf', @(x) isstruct(x));
ip.addParameter('ChannelNames',...
    strcat({'Channel'}, arrayfun(@num2str,1:length(cellInf),'Uni',0)),...
    @(x) iscellstr(x) && length(x)==length(cellInf));
ip.addParameter('Fields', [], @(x) iscellstr(x) || isstruct(x));
ip.addParameter('CellFilter', [], @(x) islogical(x) && sum(size(x)==1)==1 && ismatrix(x));
ip.addParameter('Origin', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
ip.addParameter('TimeFilter', [], @(x) islogical(x) && sum(size(x)==1)==1 && ismatrix(x));
ip.parse(cellInf, varargin{:})

% Initialise output structure
formInf = struct();

% Convert the cellInf to a struct of cellInf indexed by channel name if it
% is not already in that form:
if length(cellInf)==1 && all(structfun(@(x) isstruct(x),cellInf))
    chNames = fieldnames(cellInf);
else
    chNames = ip.Results.ChannelNames;
    arrayInf = cellInf;
    cellInf = struct();
    for ch=1:length(chNames)
        cellInf.(chNames{ch}) = arrayInf(ch);
    end
end
allfields = struct2cell(structfun(@fieldnames,cellInf,'Uni',0));

if isempty(chNames)
    % cellInf is empty, so return with empty formInf
    return
end

% Create a field selection struct based on the channel names
fieldsel = ip.Results.Fields;
if isstruct(fieldsel)
    assert(all(ismember(fieldnames(fieldsel),chNames)),...
        'All fields in the "Fields" struct should be valid channel names');
else
    if isempty(fieldsel)
        fieldsel = unique(cat(1,allfields{:}));
    end
    fieldsel = structfun(@(x) fieldsel,cellInf,'Uni',0);
end
% 'fieldsel' specifies output channels:
channels = fieldnames(fieldsel);

% Load filtering parameters
cellfilt = ip.Results.CellFilter;
origin = ip.Results.Origin;
timefilt = ip.Results.TimeFilter;

if length(cellfilt)==length(timefilt) && ~isempty(cellfilt)
    warning(['The number of cells and number of time points are identical. ',...
        'Identification of vector types may be ambiguous, and cell filtering ',...
        'will take precedence.']);
end

% Some fields may only be included in one channel. These are considered to
% be general fields that are consistent across all channels, so we make a
% struct that contains an entry for each field that is present only once:
commonfields = allfields{1};
for f=allfields', commonfields = intersect(commonfields,f{1}); end
genInf = struct();
for f=commonfields'
    hasfield = structfun(@(c) isempty(c.(f{1})), cellInf);
    if length(hasfield)>1 && sum(hasfield)==1
        genInf.(f{1}) = cellInf.(chNames{hasfield}).(f{1});
    end
end

% Determine the mean time at the 'Origin' time point
timeOffset = 0;
if ~isempty(origin)
    for ch=1:length(channels)
        chInf = cellInf.(channels{ch});
        if isfield(chInf,'times') && ~isempty(chInf.times)
            timeOffset = mean(chInf.times(:,origin));
            break
        end
    end
end

% If both a time filter and origin are specified, work out the
% time-filtered origin:
if ~isempty(timefilt) && ~isempty(origin)
    originmask = false(size(timefilt));
    originmask(origin) = true;
    origin = find(originmask(timefilt),1);
end

for ch=1:length(channels)
    channel = channels{ch};
    formInf.(channel) = struct();
    if ismember('times', fieldsel.(channel))
        if ~isempty(origin)
            formInf.(channel).origin = origin;
        elseif isfield(cellInf.(channel), 'origin')
            formInf.(channel).origin = cellInf.(channel).origin;
        end
    end
    chfields = fieldsel.(channel);
    for f=1:length(chfields)
        field = chfields{f};
        var = cellInf.(channel).(field);
        
        % Skip fields in the cellInf that are empty
        if isempty(var)
            if isfield(genInf,field), var = genInf.(field);
            else continue; end
        end
        
        % Apply cell and time filters
        if length(size(var))==2 && all(size(var)>1)
            % matrix
            if size(var,1)==length(cellfilt)
                % cell * X matrix
                var = var(cellfilt,:);
            end
            if size(var,2)==length(timefilt)
                % X * time matrix
                var = var(:,timefilt);
            end
        elseif length(var)==length(cellfilt) && sum(size(var)==1)==1
            % cell vector
            var = var(cellfilt);
        elseif length(var)==length(timefilt) && sum(size(var)==1)==1
            % time vector
            var = var(timefilt);
        end
        if issparse(var)
            var = full(var);
            if isnumeric(var), var(var==0) = NaN; end
        end
        if strcmp(field,'times')
            var = var - timeOffset;
        end
        formInf.(channel).(field) = var;
    end
end

end