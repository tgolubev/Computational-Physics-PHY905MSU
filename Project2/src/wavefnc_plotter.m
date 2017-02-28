%This plots the wavefunction for two electrons (project2 partd result) as
%functions of the relative coordinate r and different values of w_r.

clear 
[File,Path]=uigetfile('*.txt','MultiSelect','on');

N = numel(File);                                              %Counts the number of files that are in the cell array "File" that was outputted by uigetfile. 

for  num=1:N                                                  %Repeat loop for each file (each single curve analyzed seperately)

   name= File(1,num);
   str=sprintf('%s', [Path name{1}]);                         %makes str be the name of file (along with its path)
   format shortG                                              %change formating so doesn't show 0's for e-11 values. 

 
   data = importdata(str, ' ', 2);                            %importdata(FILENAME, DELIM, NHEADERLINES) loads data into a struct from 
                                                              %ASCII file FILENAME, reading numeric data starting from line NHEADERLINES+1.
                                                              %struct has fields data and text data
                                          
   r_values(:,num) = data.data(:,1);                          %to extract the data from the struct, use data.data
   wavefnc_values(:,num) = data.data(:,2);
 
   h = plot(r_values(:,num),wavefnc_values(:,num));           %plot the r_values and wavefnc_value that correspond to the current file
   set(h,'LineWidth',1.5);                              
   hold on     
   xlabel({'Relative Coordinate r'});
   ylabel({'Wavefnc(r)'});
   
end

