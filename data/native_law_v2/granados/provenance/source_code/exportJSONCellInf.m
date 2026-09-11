function exportJSONCellInf(cellInfCollection, filename, varargin)
% EXPORTJSONCELLINF Export cellInf data to a JSON file
%
%   exportJSONCellInf(cellInfCollection, filename) exports fields and 
%   channels of the cellInfCollection to a JSON file after formatting with
%   the FORMATCELLINF function. 
%   
%   'cellInfCollection' can be an EXPERIMENTTRACKING object (i.e., a 
%   cExperiment), a cellInf array obtained from an EXPERIMENTTRACKING 
%   object, a nested struct hierarchy of such cellInf arrays, or a nested 
%   struct hierarchy of formatted cellInf structs as obtained from the 
%   FORMATCELLINF function.
%   
%   The formatted structure (nested or otherwise) will be written to the
%   file 'filename' as JSON (in text format).
%   
%   Any additional parameter name/value pairs will be passed to the
%   FORMATCELLINF function.
%
%   See also: FORMATCELLINF, EXPERIMENTTRACKING

% Check input
ip = inputParser;
ip.addRequired('cellInfCollection', @(x) isstruct(x) || isa(x,'experimentTracking'));
ip.addRequired('filename', @(x) ischar(x));
ip.KeepUnmatched = true;
ip.parse(cellInfCollection, filename, varargin{:})

formatCellInfArgs = [fieldnames(ip.Unmatched),struct2cell(ip.Unmatched)]';

if isa(cellInfCollection, 'experimentTracking')
    cellInfCollection = cellInfCollection.cellInf;
end

% Recursively apply 'formatCellInf' to cellInfs in the cellInfCollection:
cellInfCollection = formatCollection(cellInfCollection);

% Export the resulting struct to the specified file as JSON:
savejson('', cellInfCollection, 'FileName', filename,...
    'ParseLogical', true, 'NaN', 'null');

    function out = formatCollection(collection)
        assert(isstruct(collection),...
            'The cellInfCollection is badly formatted.');
        
        % Check if this level of the nested struct is cellInf-like:
        iscontainer = all(structfun(@(x) ...
            isstruct(x) && length(x)==1, collection));
        formattedInf = false;
        if iscontainer
            formattedInf = any(structfun(@(s) ...
                any(structfun(@(x) ~isstruct(x), s)), collection));
        end
        if length(collection)>1 || ~iscontainer || formattedInf
            % Looks like a cellInf, so format:
            out = formatCellInf(collection, formatCellInfArgs{:});
        else
            % Looks like a struct of structs, so recurse:
            out = structfun(@formatCollection,collection,'Uni',0);
        end
    end
end