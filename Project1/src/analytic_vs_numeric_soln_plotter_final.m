%Numerical vs. analytic result plotter for project1
%Plots all the results from the selected output files on the same figure
%along with the analytic result

clear                                                         %clear variables
[File,Path]=uigetfile('*.txt','MultiSelect','on');

N = numel(File);                                              %This counts the number of files that are in the cell array "File" that was outputted by uigetfile. 


for  num=1:N                                                  %Repeat loop for each file (each single curve analyzed seperately)

   name= File(1,num);
   str=sprintf('%s', [Path name{1}]);                         %makes str be the name of file (along with its path)
   format shortG                                              %change formating so doesn't show 0's for e-11 values. 
   data= load (str);                                          %load the .txt file into matrix called "data"
    
   x_data = data(:,1); 
   num_pts = size(x_data,1);
   
   numeric_soln = data(:,2);
   analytic_soln = data(:,3); 

x_data_shifted = zeros(num_pts+2,1);
x_data_shifted(2:end-1) = x_data(1:end);                     %reshift the x_data by down 1 to be able to add a x=o as 1st value
x_data_shifted(1,1) = 0;                                     %add endpoints
x_data_shifted(end,1) = 1;

 if num==N                                                   %Only need to plot analytic solution once, and do it for the last file (this will be file w/ greatest # of pts since files are labeled w/ numbers
                                                             %of increasing iteration/n value that was used.                          
    analytic_soln_shifted = zeros(num_pts+2,1);              %do the reshifting only within this loop because only need analytic soln. once
    analytic_soln_shifted(2:end-1) = analytic_soln(1:end);   %reshift the data by down 1 to be able to add a x=o as 1st value
    analytic_soln_shifted(1,1) = 0;                          %add endpoints
    analytic_soln_shifted(end,1) = 0;   
     h = plot(x_data, analytic_soln); 
     set(h,'Color','b')
     hold on
     set(h,'LineWidth',1.5);
   else
 end

numeric_soln_shifted = zeros(num_pts+2,1);
numeric_soln_shifted(2:end-1) = numeric_soln(1:end);         %reshift the x_data by down 1 to be able to add a x=o as 1st value
numeric_soln_shifted(1,1) = 0;                               %boundary conditions
numeric_soln_shifted(end,1) = 0; 

 h = plot(x_data_shifted,numeric_soln_shifted);
 set(h,'LineWidth',1.5);                              
 hold on             
 xlabel({'x'});
 ylabel({'u(x)'});
 
end

hold off                                                    %so doesn't add curves from different program runs
   