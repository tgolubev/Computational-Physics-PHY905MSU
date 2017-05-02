%This attempts to find the transition temperature from the jump in energies

clear 
[File,Path]=uigetfile('*.txt','MultiSelect','off');
 
 str=sprintf('%s', [Path File]);                            %makes str be the name of file (along with its path)
 format shortG                                              %change formating so doesn't show 0's for e-11 values. 
 data= load (str);                                          %load the .txt file into matrix called "data"
   
    
 time = data(:,1);              
 kinetic_E = data(:,2);
 potential_E = data(:,3);
 total_E = kinetic_E + potential_E;
 temperature = data(:,4);
 diffusion_coeff = data(:,5);
 size(planet_index);     
 
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
                                                       
 %------------------------------------------------------------
 %fit the Potential energy with a line up to transition point: where E
 %jumps: Use this to find the temperature at phase transition
 
 %find the transition point and temperature by linear fitting from left graph
 %boundary and gradually including more points to fit over, until find  the
 %best fit. Best fit should corrspond to when fit over just the 1st segment
 %of the graph before have phase transition.
 
 
  clearvars potentialE_Rsquared                                                         
  clearvars phase_transition_temp
  
  minpts = 500;                                                                     %minimum # of points over which to do a fit (i.e. don't want to fit a line btw. just 2 pts)
  left_start_index = 2;   %skip the 1st data point since could be too low
      
  for j = left_start_index+minpts:0.3*size(time)      %we know the range before transition will be on  left 1/2 of the graph         
      %NOTE: FOR NOW NEED TO MODIFY THIS FACTOR MULTIPLYING SIZE(TIME) SO
      %FIT OVER RIGHT REGION,  OTHERWISE IT DOESN'T WORK PROPERTLY!
             
      time_data = time(left_start_index:j,1);
      potentialE_data = potential_E(left_start_index:j,1);                                                   
      [~, potentialE_stat] = polyfit(time_data, potentialE_data, 1);               %find linear fit parameters. Format is [temo_fit, temp_stat] = polyfit(x_value, y_values, degree_of_polynomial) where temp_fit gives the m and b values of linear fit
                                                                                  %since I don't use this temp_fit output (this one is just to calculate R^2) insert ~ to not calculate temp_fit here
      potentialE_Rsquared(j,1) =  1 - potentialE_stat.normr^2 / norm(potentialE_data-mean(potentialE_data))^2;    %calculate Rsquared and save values in column vector
 
  end
 
  [value,location] = max(potentialE_Rsquared);                                            %returns value and location (index) of maximum VB_Rsquared
                                                                 
  %redo the fit which gave maximum R^2
  time_data = time(1:location,1);
  potentialE_data = potential_E(1:location,1);                                                          
  [potentialE_fit, potentialE_stat] = polyfit(time_data, potentialE_data, 1);                   
  potentialE_Rsquared_final = 1 - potentialE_stat.normr^2 / norm(potentialE_data-mean(potentialE_data))^2 
  phase_transition_temp = temperature(location)
 
 syms x;
 potentialE_fitline = potentialE_fit(1)*x+potentialE_fit(2);
 
 h = plot(time,diffusion_coeff);   
 set(h,'LineWidth',1.5);                              
 hold on     
 set(gca,'fontsize',20, 'fontname', 'Times');   %sets the size of tick mark numbers on axes
 xlabel({'Time (s)'});
 ylabel({'Diffusion Coefficient'});
 hold off          %to not add more plot data to this figure window
 
 
 figure; 
 k = plot(time,temperature);
 hold on 
 left_axis_limit = -0.1e-9;
 %right_fit_limit = time(k);
 set(k,'LineWidth',1);      
 k=refline(left_axis_limit,phase_transition_temp);         %add reference line 
 set(k,'Color','r')
 set(k,'LineWidth',2); 
 set(k,'LineStyle',':')
 xlim([left_axis_limit inf]);
 ylim([-inf inf]);  %tell it to auto reset the axes based on what's on the plot
 set(gca,'fontsize',20, 'fontname', 'Times');   
 xlabel({'Time (s)'});
 ylabel({'Temperature (K)'});
 title('System Temperature vs. Time', 'FontSize', 24, 'FontName', 'Times');
 hold off  
 
 figure;     %create new figure window
 g = plot(time,kinetic_E);
 hold on
 g = plot(time,potential_E);
 g = plot(time,total_E);
 hold on
 g = plot(time, potentialE_fit(1)*time + potentialE_fit(2))
 %g = ezplot(potentialE_fitline, [0, inf]); 
 set(g,'Color','b');
 set(g,'LineWidth',1.5); 
 xlim([-0.1e-9 inf]);   %set left axis limit so can see the inital rise in energy at system initialization
 ylim([-inf inf]);
 title('System Energies vs. Time', 'FontSize', 24, 'FontName', 'Times');
 xlabel({'Time (s)'},'FontSize', 22, 'FontName','Times');
 ylabel({'Energy (J)'},'FontSize', 22, 'FontName','Times');
 set(gca,'fontsize',20, 'fontname', 'Times');   %sets the size of tick mark numbers on axes
 %Legend = legend('Kinetic Energy', 'Potential Energy', 'Total Energy');                         %define Legend as an object
 %legend boxoff                                         %remove the box around legend
 %set(Legend, 'FontSize', 20, 'FontName', 'Times');     %set properties of legend
 hold off
 