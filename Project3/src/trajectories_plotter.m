%This plots the results for Project 3 for any number of planets that are in
%the output file. It creates both a 3D and 2D plot.
clear 
[File,Path]=uigetfile('*.txt','MultiSelect','off');
 
 str=sprintf('%s', [Path File]);                            %makes str be the name of file (along with its path)
 format shortG                                              %change formating so doesn't show 0's for e-11 values. 
 data= load (str);                                          %load the .txt file into matrix called "data"
   
    
 planet_index = data(:,2);               %this column are the indices which identify the planets 
 
 %total_values = size(planet_index);     
 total_planets = max(planet_index);
 for i=1:(total_planets)
     planet_indices(:,i)=find(planet_index==i);  %put into ith column of planet_indices matrix indices for ith planet
     x(:,i)=data(planet_indices(:,i),4);
     y(:,i)=data(planet_indices(:,i),5);
     z(:,i)=data(planet_indices(:,i),6);
 end
     
for i=2:(total_planets)
 h = plot3(x(:,i),y(:,i),z(:,i));   %3D plot is called by plot3                          
 hold on     
end
set(h,'LineWidth',1.5);   
title('Solar System Orbits', 'FontSize', 24, 'FontName', 'Times');
xlabel({'x (AU)'},'FontSize', 22, 'FontName','Times');
ylabel({'y (AU)'},'FontSize', 22, 'FontName','Times');
zlabel({'z (AU)'},'FontSize', 22, 'FontName','Times');
set(gcf, 'PaperPositionMode', 'manual');              %Makes sure that when resize figure box while viewing, the actual figure size doesn't change
                                                         %Ensures that all saved figures have consistent size
set(gca,'fontsize',20, 'fontname', 'Times');

 figure                    %create a seperate figure window for the 2D plot
for i=2fu:(total_planets)
 g = plot(x(:,i),y(:,i));                               
 hold on     
end
set(g,'LineWidth',1.5); xlabel({'x (AU)'},'FontSize', 22, 'FontName','Times');
title('Solar System Orbits', 'FontSize', 24, 'FontName', 'Times');
ylabel({'y (AU)'},'FontSize', 22, 'FontName','Times');
set(gcf, 'PaperPositionMode', 'manual');              %Makes sure that when resize figure box while viewing, the actual figure size doesn't change
                                                         %Ensures that all saved figures have consistent size
set(gca,'fontsize',20, 'fontname', 'Times');
%hold off;
     

