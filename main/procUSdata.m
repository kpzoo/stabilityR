%% Process data from Texas state for Fig S1
function [tdate, Ideme, nDeme, nFull] = procUSdata(data, percthresh)

% Find dates and ensure none missing
tdate = unique(data.date); 
nday = length(tdate); 
if ~issorted(tdate) || nday ~= days(minus(max(tdate), min(tdate)))+1
    error('Date data incorrect');
end

% Separate by county and get timeseries
countyName = unique(data.county); 
nDeme = length(countyName); 

% Variables to extract from table
Idata = cell(1, nDeme); totcase = zeros(1, nDeme);
tstart = NaT([1, nDeme]); countydate = Idata; ts = Idata;
tend = tstart; idstart = totcase; idend = totcase;

% Search across demes and get epidemic curves
for ii = 1:nDeme
    id = strcmp(data.county, countyName{ii});
    % Entries of cases and date for county
    Idata{ii} = data.cases(id);
    countydate{ii} = data.date(id);

    % Statistics of cumulative cases and start
    tstart(ii) = min(countydate{ii});
    tend(ii) = max(countydate{ii});
    totcase(ii) = max(Idata{ii});

    % Construct time in ids from dates
    idstart(ii) = find(tdate == tstart(ii));
    idend(ii) = find(tdate == tend(ii));

    % Find date ids of the county time series
    ts{ii} = idstart(ii):idend(ii);
end

% Take most important demes by percent of cases
val = totcase/sum(totcase);
[valdeme, iddeme] = sort(val, 'descend');
% Stop as cross total case threshold
idstop = find(cumsum(valdeme) > percthresh, 1, 'first');

% Demes to extract data from
Ideme = cell(1, idstop); tdeme = Ideme; tdatedeme = Ideme;
for ii = 1:idstop
    Ideme{ii} = Idata{iddeme(ii)};
    tdeme{ii} = ts{iddeme(ii)};
    tdatedeme{ii} = countydate{iddeme(ii)};
end

% Extend each selected time series across full timeline
Ideme_full = zeros(idstop, nday);     
for ii = 1:idstop
    % Find indices of this deme's dates within full tdate
    [~, loc] = ismember(tdatedeme{ii}, tdate);    
    % Fill in counts at correct positions, leave others as zero
    Ideme_full(ii, loc) = Ideme{ii};
end

% These are cumulative cases so account for difference
tdate = tdate(2:end); Ideme = zeros(idstop, nday-1);
for ii = 1:idstop
    Ideme(ii, :) = diff(Ideme_full(ii, :));
    % Remove any negative values to 0
    Ideme(ii, Ideme(ii, :) < 0) = 0;
    % Add smoothing (trailing)
    Ideme(ii, :) = round(movmean(Ideme(ii, :), [6 0]));
end
% Rename variables for output
nFull = nDeme; nDeme = idstop; 

