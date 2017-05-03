%This plots the results for gradual heating case, analyzes the energy
%oscillations and finds the melting point by a fit of the temperature vs.
%time graph.

%Timofey Golubev

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
 %fit the temperature with a line up to transition point
 
 %find the transition point and temperature by linear fitting from left graph
 %boundary and gradually including more points to fit over, until find  the
 %best fit. Best fit should corrspond to when fit over just the 1st segment
 %of the graph before have phase transition.
 
 
  clearvars temp_Rsquared                                                          
  clearvars m
  minpts = 500;           %minimum # of points over which to do a fit (i.e. don't want to fit a line btw. just 2 pts)
  left_start_index = 2;   %skip the 1st data point since could be too low
      
  for j = left_start_index+minpts:0.2*size(time)      %we know the range before transition will be on  left 1/2 of the graph    
      %NOTE: USE 0.8*size(time) for if jump is further to the right
             
      time_data = time(left_start_index:j,1);
      temp_data = temperature(left_start_index:j,1);                         
      [~, temp_stat] = polyfit(time_data, temp_data, 1);                    %find linear fit parameters. Format is [temo_fit, temp_stat] = polyfit(x_value, y_values, degree_of_polynomial) where temp_fit gives the m and b values of linear fit
                                                                           %since I don't use this temp_fit output (this one is just to calculate R^2) insert ~ to not calculate temp_fit here
      temp_Rsquared(j,1) =  1 - temp_stat.normr^2 / norm(temp_data-mean(temp_data))^2;    %calculate Rsquared and save values in column vector
 
  end
 
  [value,location] = max(temp_Rsquared);         %returns value and location (index) of maximum VB_Rsquared
                                                               
  %redo the fit which gave maximum R^2
  time_data = time(1:location,1);
  temp_data = temperature(1:location,1);             
  [temp_fit, temp_stat] = polyfit(time_data, temp_data, 1);                   
  temp_Rsquared_final =  1 - temp_stat.normr^2 / norm(temp_data-mean(temp_data))^2 
  phase_transition_temp = mean(temperature(location-6:location))       %take average of temp's surrounding phase transition area to get accurate estimate b/c temp. oscillates
 
 %------------------------------------------------------------------------------------------------------------------------
 %Different attempt
 %temperature_avg = mean(temperature(1:5));  %start off averging temperature
%  transition_index = 0;   %initialize
%  approx_transition_temp = 0;
%  for i=1:size(time)
%      %mean(temperature(i:i+5)-temperature_avg)
%      if temperature(i+1)-temperature(i) < -30
%      %if mean(temperature(i:i+5)-temperature_avg)<-15  %looking for significant temp. drop at phase transition pt.
%          transition_index = i;   %transition is ~halfway in the range which  average over
%          %approx_transition_temp = mean(temperature(i:i+2)); %take the mean of the 1st half of this range: ie shoulud be before transition
%          break;         %don't need to look further
%      end
%      %temperature_avg = mean(temperature(i:i+5)); %update temperature_avg for next 10 points
%      
%  end     
 %times_before_transition = time(1:transition_index);
 %temps_before_transition = temperature(1:transition_index); 
 %[temp_fit, temp_stat] = polyfit(times_before_transition, temps_before_transition,1)
 %---------------------------------------------------------------------------------------------------------------
 
 syms x;
 temp_fitline = temp_fit(1)*x+temp_fit(2);
 
 h = plot(time,diffusion_coeff);   %3D plot is called by plot3
 set(h,'LineWidth',1.5);                              
 hold on     
 set(gca,'fontsize',26, 'fontname', 'Times');   %sets the size of tick mark numbers on axes
 xlabel({'Time (s)'});
 ylabel({'Diffusion Coefficient (m^2/s)'});
 hold off          %to not add more plot data to this figure window
 
 figure;     %create new figure window
 g = plot(time,kinetic_E);
 hold on
 set(gca,'fontsize',26, 'fontname', 'Times');   
 plot(time,potential_E);
 plot(time,total_E);
 xlim([-0.1e-9 inf]);   %set left axis limit so can see the inital rise in energy at system initialization
 %title('System Energies vs. Time', 'FontSize', 30, 'FontName', 'Times');
 xlabel({'Time (s)'},'FontSize', 30, 'FontName','Times');
 ylabel({'Energy (J)'},'FontSize', 30, 'FontName','Times');
 Legend = legend('Kinetic Energy', 'Potential Energy', 'Total Energy');   %define Legend as an object
 legend boxoff                                         %remove the box around legend
 set(Legend, 'FontSize', 20, 'FontName', 'Times');     %set properties of legend
 hold off
 
 figure; 
 k = plot(time,temperature);
 hold on 
 left_axis_limit = -0.1e-9;
 k=refline(left_axis_limit,phase_transition_temp);     %add reference line 
 set(k,'Color','r')
 set(k,'LineWidth',2); 
 set(k,'LineStyle',':')
 %right_fit_limit = time(k);
 k = ezplot(temp_fitline, [0, time(location)]);  
 set(k,'LineWidth',1.5);      
 xlim([left_axis_limit inf]);
 ylim([-inf inf]);  %tell it to auto reset the axes based on what's on the plot
 set(gca,'fontsize',26, 'fontname', 'Times');  
 xlabel({'Time (s)'});
 ylabel({'Temperature (K)'});
 %title('System Temperature vs. Time', 'FontSize', 30, 'FontName', 'Times');
 hold off  