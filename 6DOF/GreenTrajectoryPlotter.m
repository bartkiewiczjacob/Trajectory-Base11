function GreenTrajectoryPlotter(Lat,Long,Alt) % Units: (°,°,meters)

% If vector lengths are the same, create the .KML file
if (length(Lat) == length(Long)) && (length(Long) == length(Alt))
    
    % Filename to write .KML file to
    filename = 'Temp_TrajectoryPlot.kml';
    
    filenameTaken = exist(filename,'file');
    count = 0;
    while filenameTaken == 2
        count = count + 1;
        filename = strcat('Temp_TrajectoryPlot',num2str(count),'.kml');
        filenameTaken = exist(filename,'file');
    end
    
    % Open file
    fileID = fopen(filename,'w');
    
    % Determine number of points to plot
    numPoints = length(Lat);
    
    % Begin writing .KML file
    fprintf(fileID,'<?xml version="1.0" encoding="UTF-8"?>\n');
    fprintf(fileID,'<kml xmlns="http://www.opengis.net/kml/2.2" ');
    fprintf(fileID,'xmlns:gx="http://www.google.com/kml/ext/2.2" ');
    fprintf(fileID,'xmlns:kml="http://www.opengis.net/kml/2.2" ');
    fprintf(fileID,'xmlns:atom="http://www.w3.org/2005/Atom">\n');
    fprintf(fileID,'<Document>\n');
    fprintf(fileID,'\t<name>%s</name>\n',filename);
    fprintf(fileID,'\t<Placemark>\n');
    fprintf(fileID,'\t\t<Style>\n');
    fprintf(fileID,'\t\t\t<LineStyle>\n');
    fprintf(fileID,'\t\t\t\t<color>6414f000</color>\n');
    fprintf(fileID,'\t\t\t\t<width>2</width>\n');
    fprintf(fileID,'\t\t\t</LineStyle>\n');
    fprintf(fileID,'\t\t</Style>\n');
    fprintf(fileID,'\t\t<LineString>\n');
    fprintf(fileID,'\t\t\t<tessellate>1</tessellate>\n');
    fprintf(fileID,'\t\t\t<altitudeMode>absolute</altitudeMode>\n');
    fprintf(fileID,'\t\t\t<coordinates>\n');
    
    % Add coordinates to file [Longitude,Latitude,Altitude] (°,°,meters)
    for i = 1:numPoints
        fprintf(fileID,'\t\t\t\t%.10f,%.10f,%.10f\n',Long(i),Lat(i),Alt(i));
    end
    
    % Finish writing .KML file
    fprintf(fileID,'\t\t\t</coordinates>\n');
    fprintf(fileID,'\t\t</LineString>\n');
    fprintf(fileID,'\t</Placemark>\n');
    fprintf(fileID,'</Document>\n');
    fprintf(fileID,'</kml>\n');
    
    % Close .KML file
    fclose(fileID);
    
    % Open .KML file in Google Earth
    winopen(filename)
    
% Vector lengths not the same -> print error message and return
else
    fprintf('PlotFlightPath: Vectors passed not all same length.\n')
    return
    
end