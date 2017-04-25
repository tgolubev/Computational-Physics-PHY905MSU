%This plots the results for Project 3 part e.

clear 
[File,Path]=uigetfile('*.txt','MultiSelect','off');
 
 str=sprintf('%s', [Path File]);                            %makes str be the name of file (along with its path)
 format shortG                                              %change formating so doesn't show 0's for e-11 values. 
 data= load (str);                                          %load the .txt file into matrix called "data"
   
    
 time = data(:,1);               %this column are the indices which identify the planets 
 %THE UNITS FOR TIME ARE WEIRD: CONVERT EITHER HERE OR IN CPP CODE!
 kinetic_E = data(:,2);
 potential_E = data(:,3);
 total_E = kinetic_E + potential_E;
 temperature = data(:,4);
 diffusion_coeff = data(:,5);
 
 %total_values = size(planet_index);     
 
 %energies analysis
 average_KE = mean(kinetic_E)
 average_PE = mean(potential_E)
 average_total_energy = mean(total_E)
 KE_amplitude = 0.5*(max(kinetic_E)-min(kinetic_E))
 PE_amplitude = 0.5*(max(potential_E)-min(potential_E))
 total_energy_amplitude = 0.5*(max(total_E)-min(total_E))
 
 rel_KE_change = KE_amplitude/average_KE
 rel_PE_change = PE_amplitude/abs(average_PE)
 rel_total_change = total_energy_amplitude/abs(average_total_energy)
 
 set(gcf, 'PaperPositionMode', 'manual');              %Makes sure that when resize figure box while viewing, the actual figure size doesn't change
                                                         %Ensures that all saved figures have consistent size
                                                        
 
 %axes1 = gca;   %create an axes object
 %axes1.Position = [ 0.11 0.775 0.815]; 

 
 h = plot(time,diffusion_coeff);   %3D plot is called by plot3
 set(h,'LineWidth',1.5);                              
 hold on     
 set(gca,'fontsize',20, 'fontname', 'Times');   %sets the size of tick mark numbers on axes
 xlabel({'time (fs)'});
 ylabel({'Diffusion Coefficient'});
 hold off          %to not add more plot data to this figure window
 
 figure;     %to create new figure window
 g = plot(time,kinetic_E);
 hold on
 set(gca,'fontsize',20, 'fontname', 'Times');   %sets the size of tick mark numbers on axes
 plot(time,potential_E);
 title('System Energies vs. Time', 'FontSize', 24, 'FontName', 'Times');
 xlabel({'Time (fs)'},'FontSize', 22, 'FontName','Times');
 ylabel({'Energy (arb. units)'},'FontSize', 22, 'FontName','Times');
 hold off
 
%NOTE: KE and PE have sharp change near time=0