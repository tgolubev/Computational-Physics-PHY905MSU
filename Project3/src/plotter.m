%This plots the wavefunction for two electrons (project2 partd result) as
%functions of the relative coordinate r and different values of w_r.

clear 
[File,Path]=uigetfile('*.txt','MultiSelect','off');
 
 data = importdata([Path File], ' ', 0);     %importdata(FILENAME, DELIM, NHEADERLINES) loads data into a struct from 
                                             %ASCII file FILENAME, reading numeric data starting from line NHEADERLINES+1.
                                             %struct has fields data and
                                             %textdata
    

 x = data(:,2);             %to extract the data from the struct, use data.data
 y = data(:,3);
 
 h = plot(x,y);
 set(h,'LineWidth',1.5);                              
 hold on     
 xlabel({'x'});
 ylabel({'y'});
 