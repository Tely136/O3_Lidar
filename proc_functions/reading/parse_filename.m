function dt = parse_filename(name)
%PARSE_FILENAME Parse date/time encoded in filename prefixes.
%   dt = PARSE_FILENAME(name) accepts a string scalar or string array.
%   Invalid or missing names return datetime(1900,1,1).

    fallback = datetime(1900, 1, 1);

    name = string(name);
    dt = repmat(fallback, size(name));

    for i = 1:numel(name)
        try
            if ismissing(name(i)) || strlength(name(i)) == 0
                continue
            end

            baseName = extractBefore(name(i), ".");
            if baseName == ""
                baseName = name(i);
            end

            if strlength(baseName) < 9
                continue
            end

            yyStr = extractBetween(baseName, 3, 4);
            mmHex = extractBetween(baseName, 5, 5);
            ddStr = extractBetween(baseName, 6, 7);
            hhStr = extractBetween(baseName, 8, 9);

            yy = str2double(yyStr);
            dd = str2double(ddStr);
            hh = str2double(hhStr);
            mm = hex2dec(char(mmHex));

            if any(isnan([yy, dd, hh])) || mm < 1 || mm > 12
                continue
            end

            dt(i) = datetime(2000 + yy, mm, dd, hh, 0, 0);
        catch
            dt(i) = fallback;
        end
    end
end
