%Plotter for PHY905 program2 eigenvector results.
%This is for input files which contain only  1 eigenvector each.
%This can be easily modified for plotting any data from .txt file which has
%2 columns (x and y values) by simply changing the NHEADERLINES in
%importdata() function and the plot properties.
%Can select multiple result files at once. 

%Coded by Timofey Golubev

clear                                                         %clear variables
[File,Path]=uigetfile('*.txt','MultiSelect','off');

N = numel(File);                                              %Counts the number of files that are in the cell array "File" that was outputted by uigetfile. 

for  num=1:N   
    
                                                
   name= File(1,num);
   
   str=sprintf('%s', [Path name{1}]);                         %makes str be the name of file (along with its path)
   format shortG                                              %change formating so doesn't show 0's for e-11 values.                     
   
   data = importdata(str, ' ', 2);              %importdata(FILENAME, DELIM, NHEADERLINES) loads data into a struct from 
                                                %ASCII file FILENAME, reading numeric data starting from line NHEADERLINES+1.
                                                %struct has fields data and textdata
                                                
    x_values =  data.data(:,1);  
    y_values =  data.data(:,2);
   
   h = plot(x_values, y_values, 'LineWidth',1.5);  
   hold on;
end
    
   set(gcf, 'PaperPositionMode', 'manual');              %Makes sure that when resize figure box while viewing, the actual figure size doesn't change
                                                         %Ensures that all saved figures have consistent size
                                                        
   set(gca,'fontsize',20, 'fontname', 'Times');
   axes1 = gca;   %create an axes object
   axes1.Position = [0.13 0.11 0.775 0.815];                     %axes the Position parameter of the object "axes1"
   title('\omega_{r} = ', 'FontSize', 24, 'FontName', 'Times');  %'FontWeight', 'normal' would make font non-bold
   xlabel({'Relative Coordinate r (Bohr Radii)'},'FontSize', 22, 'FontName','Times');
   ylabel({'Ground State Wavefunction (r)'}, 'FontSize', 22, 'FontName','Times');
   %Set Legend names to file names
   name1 = char(File(1,1));
   name2 = char(File(1,2));
   Legend = legend(name1, name2);                         %define Legend as an object
   legend boxoff                                         %remove the box around legend
   set(Legend, 'FontSize', 20, 'FontName', 'Times');     %set properties of legend
   hold off;