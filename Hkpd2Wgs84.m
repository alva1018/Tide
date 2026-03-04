function H = Hkpd2Wgs84(h, pos)
% Hkpd2Wgs84 Convert heights from HKPD datum to WGS84 datum
%   input: h - heights in HKPD
%          pos - [lat, lon, alt] of the location in degrees and meters
%   output: H - heights in WGS84

    if pos(2) > 114.2
            N = 2.42;
    else
            N = 2.38; 
    end
    H = h + N;
end