%This plots the results for Project 3 part_b, the object oriented solar
%system

clear 
[File,Path]=uigetfile('*.txt','MultiSelect','off');
 
 str=sprintf('%s', [Path File]);                            %makes str be the name of file (along with its path)
 format shortG                                              %change formating so doesn't show 0's for e-11 values. 
 data= load (str);                                          %load the .txt file into matrix called "data"
   
    
 planet_index = data(:,2);               %this column are the indices which identify the planets 
 
 %total_values = size(planet_index);     
 
 planet1_indices = find(planet_index==1);   %returns a column vector of the indices (locations) of where the values in the vector planet_index==1
 planet2_indices = find(planet_index==2);
 
 x1 = data(planet1_indices,4);    %import into x1 the values from data which have row #'s = planet1_indices, and are in column 4 (easy!)
 x2 = data(planet2_indices,4);
 y1 = data(planet1_indices,5);
 y2 = data(planet2_indices,5);
 z1 = data(planet1_indices,6);
 z2 = data(planet2_indices,6);

 
 h = plot3(x1,y1,z1);   %3D plot is called by plot3
 set(h,'LineWidth',1.5);                              
 hold on     
 xlabel({'x'});
 ylabel({'y'});
 plot3(x2,y2,z2);
 