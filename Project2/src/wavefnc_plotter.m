%This plots the wavefunction for two electrons (project2 partd result) as
%functions of the relative coordinate r and different values of w_r.

clear 
[File,Path]=uigetfile('*.txt','MultiSelect','off');
 
 data = importdata([Path File], ' ', 2);     %importdata(FILENAME, DELIM, NHEADERLINES) loads data into a struct from 
                                             %ASCII file FILENAME, reading numeric data starting from line NHEADERLINES+1.
                                             %struct has fields data and
                                             %textdata
    

 r = data.data(:,1);             %to extract the data from the struct, use data.data
 wavefnc1 = data.data(:,2);
 wavefnc2 = data.data(:,3);
 wavefnc3 = data.data(:,4);
 wavefnc4 = data.data(:,5);
 
 h = plot(r,wavefnc1);
 set(h,'LineWidth',1.5);                              
 hold on     
 xlabel({'r'});
 ylabel({'Wavefnc(r)'});
 plot(r, wavefnc2);
 plot(r, wavefnc3);
 plot(r, wavefnc4);
