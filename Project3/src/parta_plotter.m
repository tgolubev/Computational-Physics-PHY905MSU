%This plots the results for Project 3 part_a, non-object oriented earth-sun
%system.

clear 
[File,Path]=uigetfile('*.txt','MultiSelect','off');
 
 data = importdata([Path File], ' ', 0);     %importdata(FILENAME, DELIM, NHEADERLINES) loads data into a struct from 
                                             %ASCII file FILENAME, reading numeric data starting from line NHEADERLINES+1.
                                             %struct has fields data and
                                             %textdata
    %to extract the data from the struct, use data.data
    
    %HOWEVER, if heave NHEADERLINES=0, then the importdata doesn't return a
    %struct, but returns a regular matrix, so extract data as in below.

 x = data(:,2);             
 y = data(:,3);
 
 h = plot(x,y);
 set(h,'LineWidth',1.5);                              
 hold on     
 xlabel({'x'});
 ylabel({'y'});
 