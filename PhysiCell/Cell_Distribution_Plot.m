clear
clc

cd output/


%%
MCDS = read_MultiCellDS_xml( 'output00000000.xml'); 
% simple_plot( MCDS, MCDS.discrete_cells.live_cells, 'r' )
%%

% make it easier to work with the cell positions; 
P = MCDS.discrete_cells.state.position;
 
% find type 1 cells
ind1 = find( MCDS.discrete_cells.metadata.type == 1 ); 
hold on
plot3( P(ind1,1), P(ind1,2), P(ind1,3), 'ro' )
hold off
axis( 1000*[-1 1 -1 1 -1.6 1.6] )

%%
cd ..