function [week, tow] = utc2gps(utc)
    % input: utc - datetime object in UTC time
    % output: week - GPS week number
    %         tow  - Time of Week in seconds
    
    % transform UTC datetime to Julian Date
    Y = year(utc); M = month(utc); D = day(utc);
    if M <= 2, Y = Y - 1; M = M + 12; end
    JD = floor(365.25 * Y) + floor(30.6001 * (M + 1)) + D + 1720981.5 + ...
         hour(utc)/24 + minute(utc)/1440 + second(utc)/86400;
     
    % GPS epoch start date in Julian Date
    JD_GPS_START = 2444244.5;
    
    diff_days = JD - JD_GPS_START;
    week = floor(diff_days / 7);
    tow = (diff_days - week * 7) * 86400;
    
    % add leap seconds (as of 2024, there are 18 leap seconds)
    leap_seconds = 18; 
    tow = tow + leap_seconds;
    
    % adjust tow and week if tow exceeds one week
    if tow >= 604800
        tow = tow - 604800;
        week = week + 1;
    end
end