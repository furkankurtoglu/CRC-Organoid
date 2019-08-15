clear
clc

cd output/


%%
MCDS = read_MultiCellDS_xml( 'output00000000.xml');

k = find( MCDS.mesh.Z_coordinates == 10 );
%%
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
 
%%
cd ..