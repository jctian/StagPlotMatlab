clear

start_dir = 'D:\Polybox\ETH\2020Fall\Research\Plotting\PlotSurfaceFields';
%directory = 'D:\WSLFiles\Output\Venus\210424_Default_oldisl\210424_1e20\+op'
directory = 'F:\StagOutput\Thesis\210812_Oldiff0.3Ma\+op'
file_name = 'Venus';
plate_analyse = strcat(directory, '\',file_name, '_plates_analyse.dat');


starting_frame = 1400;
num = 1;

total_frame = size(importdata(plate_analyse).data, 1)-1

%total_frame = 330

cd(directory)
include.topography=true;
include.age       =true;
include.strainrate=true;
include.crustthick=true;
include.geoid     =false;
include.surfacemelt = true;
include.crustmelt = true; % sum of melt within crust
include.crustage = true;

struct2table(include)

for i = starting_frame:num:total_frame
   PlotSurfFields_function(directory, file_name, i, include)  
   disp(strcat('Ploting frame ', num2str(i)))
end

cd(start_dir)

