function connectorWeight = SectionConnectorWER(diameter,wallThickness,sectionCount)

%Preset connector values that can be set to change for different stresses
%if need be
connectorThickness = 1/8;   %Thickness of the connecting portion of the section connector
midThickness = .25; %Thickness of the section divider in the middle of the connector
materialDens = 0.1; %Density in lb/in^3 for aluminum
connectionLength = 2;   %Length one side of the connector extends into tubing

rocketID = diameter - wallThickness;   %Calculates inner diameter for section
connectorID = rocketID - connectorThickness;    %calculates the inner diameter of the connector itself

midSectionVol = pi()*((diameter/2)^2)*midThickness; %Volume of the middle section divider

connectionVol = 2*pi()*((rocketID^2)-(connectorID^2));  %Volume of both sides that connect the rocket.

%Calculates total volume to find weight using density
totalVol = midSectionVol+connectionVol;

%Calculates weight of connector for total number of connectors
connectorWeight = (sectionCount - 1)*totalVol*materialDens;
end

