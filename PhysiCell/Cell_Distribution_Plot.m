clear
close all
clc

%%
cd output/


%%
MCDS = read_MultiCellDS_xml( 'output00000000.xml'); 
% simple_plot( MCDS, MCDS.discrete_cells.live_cells, 'r' )
%%

% make it easier to work with the cell positions; 
P = MCDS.discrete_cells.state.position;
 
% find type 1 cells
ind1 = find( MCDS.discrete_cells.metadata.type == 1 ); 
ind2 = find( MCDS.discrete_cells.metadata.type == 2 ); 
plot3( P(ind1,1), P(ind1,2), P(ind1,3), 'ro' )
hold on
plot3( P(ind2,1), P(ind2,2), P(ind2,3), 'bo' )
hold off
axis image
axis(100*[-5 5 -8 8 -5 5] )


%% Oxygen
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
k = find( MCDS.mesh.Z_coordinates == 10 );
figure(4)
contourf( MCDS.mesh.X(:,:,k), MCDS.mesh.Y(:,:,k), MCDS.continuum_variables(3).data(:,:,k) , 20 ) ;

axis image
colorbar 
xlabel( sprintf( 'x (%s)' , MCDS.metadata.spatial_units) ); 
ylabel( sprintf( 'y (%s)' , MCDS.metadata.spatial_units) );

title( sprintf('%s (%s) at t = %3.2f %s, z = %3.2f %s', MCDS.continuum_variables(3).name , ...
     MCDS.continuum_variables(3).units , ...
     MCDS.metadata.current_time , ...
     MCDS.metadata.time_units, ... 
     MCDS.mesh.Z_coordinates(k), ...
     MCDS.metadata.spatial_units ) ); 


%% Lactate
k = find( MCDS.mesh.Z_coordinates == 10 );
figure(5)
contourf( MCDS.mesh.X(:,:,k), MCDS.mesh.Y(:,:,k), MCDS.continuum_variables(4).data(:,:,k) , 20 ) ;

axis image
colorbar 
xlabel( sprintf( 'x (%s)' , MCDS.metadata.spatial_units) ); 
ylabel( sprintf( 'y (%s)' , MCDS.metadata.spatial_units) );

title( sprintf('%s (%s) at t = %3.2f %s, z = %3.2f %s', MCDS.continuum_variables(4).name , ...
     MCDS.continuum_variables(4).units , ...
     MCDS.metadata.current_time , ...
     MCDS.metadata.time_units, ... 
     MCDS.mesh.Z_coordinates(k), ...
     MCDS.metadata.spatial_units ) ); 
cd ..