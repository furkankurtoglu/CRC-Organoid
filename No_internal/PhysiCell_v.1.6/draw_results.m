clear
close all
clc

%%
cd C:\Users\Furkan\Documents\GitHub\CRC-Organoid\No_internal\PhysiCell_v.1.6\output

%%
s=what;
MatFiles = s.mat;
OutMatFiles = MatFiles(contains(MatFiles,'micro'));
OutMatFiles(1) = [];
OutMatFiles(1) = [];
for i = 1:length(OutMatFiles)
    OutMatFiles{i}=OutMatFiles{i}(1:14);
end

edge_cell_energy = zeros(1,length(OutMatFiles));
center_cell_energy = zeros(1,length(OutMatFiles));

center_glucose = zeros(1,length(OutMatFiles));
center_oxygen = zeros(1,length(OutMatFiles));

edge_glucose = zeros(1,length(OutMatFiles));
edge_oxygen = zeros(1,length(OutMatFiles));



%%
NoofCells=zeros(1,length(OutMatFiles));

% for i = 1:2
for i = 1:length(OutMatFiles)

    xmlname=strcat(OutMatFiles{i},'.xml');
    MCDS = read_MultiCellDS_xml(xmlname);
    % simple_plot( MCDS, MCDS.discrete_cells.live_cells, 'r' )
    %%
%     
%     % make it easier to work with the cell positions;
%     P = MCDS.discrete_cells.state.position;
%     
%     % find type 1 cells
%     ind1 = find( MCDS.discrete_cells.metadata.type == 1 );
%     ind2 = find( MCDS.discrete_cells.metadata.type == 2 );
%     plot3( P(ind1,1), P(ind1,2), P(ind1,3), 'ro' )
%     hold on
%     plot3( P(ind2,1), P(ind2,2), P(ind2,3), 'bo' )
%     hold off
%     axis image
%     axis(100*[-5 5 -8 8 -5 5] )
%     v = [0 70 -180];
%     [caz,cel]=view(v);
%     xlabel('micron')
%     ylabel('micron')
%     zlabel('micron')
%     
%     NoCells=length(MCDS.discrete_cells.live_cells)-10000;
%     NoofCells(i)=NoCells;
%     dim = [.2 .5 .3 .3];
%     time=num2str(i*60/24/60);  
%     annotation('textbox',dim,'String',NoCells,'FitBoxToText','on');
%     title(['Time = ',time,' days'])
%     set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
%     figure(1)
%     
%     %% Oxygen
    %% Glucose
    k = find( MCDS.mesh.Z_coordinates == 10 );
    figure(2)
    contourf( MCDS.mesh.X(:,:,k), MCDS.mesh.Y(:,:,k), MCDS.continuum_variables(1).data(:,:,k) , 20 ) ;

    axis image
    colorbar
    xlabel( sprintf( 'x (%s)' , MCDS.metadata.spatial_units) );
    ylabel( sprintf( 'y (%s)' , MCDS.metadata.spatial_units) );

    title( sprintf('%s (%s) at t = %3.2f %s, z = %3.2f %s', MCDS.continuum_variables(1).name , ...
        MCDS.continuum_variables(1).units , ...
        MCDS.metadata.current_time , ...
        MCDS.metadata.time_units, ...
        MCDS.mesh.Z_coordinates(k), ...
        MCDS.metadata.spatial_units ) );
    %% Glucose
    k = find( MCDS.mesh.Z_coordinates == 10 );
    figure(3)
    contourf( MCDS.mesh.X(:,:,k), MCDS.mesh.Y(:,:,k), MCDS.continuum_variables(2).data(:,:,k) , 20 ) ;

    axis image
    colorbar
    xlabel( sprintf( 'x (%s)' , MCDS.metadata.spatial_units) );
    ylabel( sprintf( 'y (%s)' , MCDS.metadata.spatial_units) );

    title( sprintf('%s (%s) at t = %3.2f %s, z = %3.2f %s', MCDS.continuum_variables(2).name , ...
        MCDS.continuum_variables(2).units , ...
        MCDS.metadata.current_time , ...
        MCDS.metadata.time_units, ...
        MCDS.mesh.Z_coordinates(k), ...
        MCDS.metadata.spatial_units ) );

    %% Glutamine
%     k = find( MCDS.mesh.Z_coordinates == 10 );
%     figure(4)
%     contourf( MCDS.mesh.X(:,:,k), MCDS.mesh.Y(:,:,k), MCDS.continuum_variables(3).data(:,:,k) , 20 ) ;
%     
%     axis image
%     colorbar
%     xlabel( sprintf( 'x (%s)' , MCDS.metadata.spatial_units) );
%     ylabel( sprintf( 'y (%s)' , MCDS.metadata.spatial_units) );
%     
%     title( sprintf('%s (%s) at t = %3.2f %s, z = %3.2f %s', MCDS.continuum_variables(3).name , ...
%         MCDS.continuum_variables(3).units , ...
%         MCDS.metadata.current_time , ...
%         MCDS.metadata.time_units, ...
%         MCDS.mesh.Z_coordinates(k), ...
%         MCDS.metadata.spatial_units ) );
% %     
%     
     %% Lactate
%     k = find( MCDS.mesh.Z_coordinates == 10 );
%     figure(5)
%     contourf( MCDS.mesh.X(:,:,k), MCDS.mesh.Y(:,:,k), MCDS.continuum_variables(4).data(:,:,k) , 20 ) ;
%     
%     axis image
%     colorbar
%     xlabel( sprintf( 'x (%s)' , MCDS.metadata.spatial_units) );
%     ylabel( sprintf( 'y (%s)' , MCDS.metadata.spatial_units) );
%     
%     title( sprintf('%s (%s) at t = %3.2f %s, z = %3.2f %s', MCDS.continuum_variables(4).name , ...
%         MCDS.continuum_variables(4).units , ...
%         MCDS.metadata.current_time , ...
%         MCDS.metadata.time_units, ...
%         MCDS.mesh.Z_coordinates(k), ...
%         MCDS.metadata.spatial_units ) );
%      
 %   edge_cell_energy(i) = MCDS.discrete_cells.custom.energy(8624);
%    edge_oxygen(i) = MCDS.discrete_cells.custom.sensed_oxygen(8624);
   % edge_glucose(i) = MCDS.discrete_cells.custom.sensed_glucose(8624);
  %  
  %  center_cell_energy(i) = MCDS.discrete_cells.custom.energy(9031);
  %  center_oxygen(i) = MCDS.discrete_cells.custom.sensed_oxygen(9031);
  %  center_glucose(i) = MCDS.discrete_cells.custom.sensed_glucose(9031);

end 
%%
figure()
plot(edge_cell_energy)
title('edge cell energy')
hold off
figure()
plot(center_cell_energy)
title('center cell energy')
%%
figure()
plot(edge_oxygen)
title('edge oxygen')

figure()
plot(edge_glucose)
title('edge glucose')


figure()
plot(center_oxygen)
title('center oxygen')

figure()
plot(center_glucose)
title('center glucose')


%%
cd C:\Users\Furkan\Documents\GitHub\CRC-Organoid\No_internal\PhysiCell_v.1.6