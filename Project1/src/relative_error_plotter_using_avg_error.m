%Plotter for PHY905 program1 results.
%Can select multiple result files at once. Finds the maximum relative error
%value for each result file (each iteration) and uses this to plot the
%relative error as fnc. of step size.
%For each plot line desired, need to run the code again. It has 'hold on'
%enabled, so each new generated line will be added to the currently opened
%figure.

%Coded by Timofey Golubev

clear                                                         %clear variables
[File,Path]=uigetfile('*.txt','MultiSelect','on');

N = numel(File);                                              %Counts the number of files that are in the cell array "File" that was outputted by uigetfile. 

for  num=1:N                                                  %Repeat loop for each file (each single curve analyzed seperately)

   name= File(1,num);
   str=sprintf('%s', [Path name{1}]);                         %makes str be the name of file (along with its path)
   format shortG                                              %change formating so doesn't show 0's for e-11 values. 
   data= load (str);                                          %load the .txt file into matrix called "data"
    
   rel_error = data(:,4);                                     %log was already taken in cpp program
   avg_error(num,1) = mean(rel_error);                        %finds the maximum error in the column for that specific iteration and puts it into row = num(corresponding to which iteration/file its from)
   n = size(rel_error,1)+1;                                   %From the number of rows (1st dimension of rel_error: synthax is size(variable, dimension) where for 
                                                              %dimension enter 1 or 2 (row or column) in the results file, add 1 to get matrix size that was used.
                 
   step_size(num,1) = 1/n;                                    %can't use h() because that's reserved for handles in Matlab
   log_h = log10(step_size);
   
end

 set(0,'defaulttextinterpreter','tex')                       %Sets .tex format as the default text interpreter. This displays i.e. 'log_{10}' with 10 being a subscript on labels. 
                                                             %Use i.e. 'words^{superscript}' to make a superscript.
                                        
 g = plot(log_h, avg_error);
 set(g,'Color','b')
 set(g,'LineWidth',1.5);                              
 hold on                                                     %here we purposely leave hold on, so can run this code again to generate another line and add it to existing figure                                                    
 
 xlabel({'log_{10}(h)'});
 ylabel({'log_{10}(error)'});

 
