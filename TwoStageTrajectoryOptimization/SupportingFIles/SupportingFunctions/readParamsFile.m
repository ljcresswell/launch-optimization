function params = readParamsFile(filename)
%Reads parameter text file into a struct
%
% params = readParamsFile(filename)

    fid = fopen(filename, 'r');
    if fid == -1
        error('Could not open file: %s', filename);
    end

    params = struct();

    while ~feof(fid)
        line = strtrim(fgetl(fid));

        % Skip empty lines or comments
        if isempty(line) || startsWith(line, '%')
            continue;
        end

        % Remove inline comments (if any)
        line = strtok(line, '%');

        % Split on '='
        tokens = strsplit(line, '=');
        if numel(tokens) ~= 2
            continue; % skip malformed lines
        end

        name = strtrim(tokens{1});
        valueStr = strtrim(tokens{2});

        % Convert numeric value if possible
        value = str2num(valueStr); %#ok<ST2NM>
        if isempty(value)
            value = valueStr; % fallback to string
        end

        % Store in struct
        params.(name) = value;
    end
    fclose(fid);

    % Determine rEndAtmo
    [PresSeaLevel, ~] = atmosphere(0);
    targetP = params.atmoEndPercent * PresSeaLevel;
    h = 0:100:150000;
    [poh, ~] = atmosphere(h);
    idx = find(poh<targetP,1,"first");
    if isempty(idx)
        error('Error: Please Raise atmoEndPercent.');
    end
    params.rEndAtmo = h(idx) + params.radiusEarth;
    if params.rEndAtmo >= params.a*(1-params.e)
        error('Error: Orbit Intersects Atmosphere.')
    end
end

