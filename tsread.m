function [ freq, P, t, R0 ] = tsread( fileName )
% [ freq, P, t, R0 ] = tsread( fileName )
%
% Parses multiport network parameters in Touchstone format
% Now only handles up to 4 ports
%

freq = [];
P    = [];
t    = 's';
R0   = 50;

f = fopen( fileName, 'rt');
    
if f < 1
    return;
end

multiplier = 1e6;
format = 'ma';
R0 = 50;

% One row array per frequency point as read from the file
fpdata = { };

while ~feof(f)
    linestr = fgetl(f);

    if isempty(linestr) || '!' == linestr(1)
        % skip empty lines and comments
    elseif '#' == linestr(1)
        prevtok = '';
        linestr = linestr(2:end);
        while !isempty( strtrim( linestr ) )
            [ token, linestr ] = strtok(linestr, ' ');
            lowertok = lower(token);
            if strcmp('ri', lowertok) || strcmp('ma', lowertok)
                format = lowertok;
            elseif strcmp('z', lowertok) || strcmp('y', lowertok) || strcmp('s', lowertok)
                t = lowertok;
            elseif strcmp(prevtok, 'r')
                R0 = str2double( token );
            elseif strcmp('hz', lowertok)
                % do nothing, we assume hz
            elseif strcmp('r', lowertok)
                % do nothing, value will be parsed next
            else
                error( 'Unrecognized format option %s', token )
            end
            prevtok = lowertok;
        end
    else
        linevals = [ ];
        while !isempty( strtrim( linestr ) )
            [ token, linestr ] = strtok(linestr, ' ');
            linevals = [ linevals, str2double( token ) ];
        end
        if mod(length(linevals),2) % odd -- new frequency starts here
            fpdata{end+1} = [ ];
        end
        fpdata{end} = [ fpdata{end}, linevals ];
    end
end

fclose(f);

% Now reformat the data
for fp=fpdata
    freq = [ freq ; fp{1}(1) ];
    rval = fp{1}(2:end);
    if strcmp( format, 'ma' )
        vals = rval(1:2:end-1).*exp(j*deg2rad(rval(2:2:end)));
    elseif strcmp( format, 'ri' )
        vals = rval(1:2:end-1) + j*rval(2:2:end);
    else
        error('Unrecognized format')
    end
    N = round( sqrt( length(vals) ) );
    V = zeros(N, N);
    for r=1:N
        for  c=1:N
            % 1 or 2 ports - column-major, row-major otherwise
            if N <= 2
                vi = r + (c-1)*N;
            else
                vi = (r-1)*N + c;
            end
            V( r, c ) = vals( vi );
        end
    end
    P = cat(3, P, V);
end


