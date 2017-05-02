%This plots results statistics.txt file and analyses the energy
%oscillations.

clear 
[File,Path]=uigetfile('*.txt','MultiSelect','off');
 
 str=sprintf('%s', [Path File]);                            %makes str be the name of file (along with its path)
 format shortG                                              %change formating so doesn't show 0's for e-11 values. 
 data= load (str);                                          %load the .txt file into matrix called "data"
   
    
 time = data(:,1);               %this column are the indices which identify the planets 
 kinetic_E = data(:,2);
 potential_E = data(:,3);
 total_E = kinetic_E + potential_E;
 temperature = data(:,4);
 diffusion_coeff = data(:,5);   
 
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
                                                       
 h = plot(time,diffusion_coeff);  
 hold on     
 set(gca,'fontsize',24, 'fontname', 'Times');   %sets the size of tick mark numbers on axes
 xlabel({'Time (s)'});
 ylabel({'{Diffusion Coefficient} (m^2/s)'});
 title('Diffusion Coefficient vs. Time', 'FontSize', 26, 'FontName', 'Times');
 %ylabel({'{Diffusion Coefficient} ($m^2/s$)'},'Interpreter','tex');   %allows to use latex synthax in labels
 hold off          %to not add more plot data to this figure window

 figure;     %to create new figure window
 g = plot(time,kinetic_E);
 hold on
 set(gca,'fontsize',24, 'fontname', 'Times');   
 plot(time,potential_E);
 plot(time,total_E);
 xlim([-0.1e-9 inf]);   %set left axis limit so can see the inital rise in energy at system initialization
 title('System Energies vs. Time', 'FontSize', 26, 'FontName', 'Times');
 xlabel({'Time (s)'},'FontSize', 24, 'FontName','Times');
 ylabel({'Energy (J)'},'FontSize', 24, 'FontName','Times');
 Legend = legend('Kinetic Energy', 'Potential Energy', 'Total Energy');  %define Legend as an object
 legend boxoff                                         %remove the box around legend
 set(Legend, 'FontSize', 20, 'FontName', 'Times');     %set properties of legend
 hold off
 
 figure; 
 k = plot(time,temperature);
 hold on 
 left_axis_limit = -0.1e-10; 
 set(k,'LineWidth',1);      
 xlim([left_axis_limit inf]);
 ylim([-inf inf]);  %tell it to auto reset the axes based on what's on the plot
 set(gca,'fontsize',24, 'fontname', 'Times');  
 xlabel({'Time (s)'});
 ylabel({'Temperature (K)'});
 title('System Temperature vs. Time', 'FontSize', 26, 'FontName', 'Times');
 hold off  