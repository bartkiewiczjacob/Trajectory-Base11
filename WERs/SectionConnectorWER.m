function connectorWeight = SectionConnectorWER(diam1,diam2)

%Preset connector values that can be set to change for different stresses
%if need be
wallThickness = 1/4;
connectorLength = 4;
materialDens = 0.1; %Density in lb/in^3

%Calculates Volumes for lower and upper connector to connect two different
%diameter sections
upperSectionVol = (connectorLength/2)*pi()*((diam2/2)^2-((diam2/2)-wallThickness)^2);
lowerSectionVol = (connectorLength/2)*pi()*((diam1/2)^2-((diam1/2)-wallThickness)^2);

%Calculates total volume to find weight using density
totalVol = upperSectionVol+lowerSectionVol;

%Calculates weight of connector
connectorWeight = totalVol*materialDens;
end

