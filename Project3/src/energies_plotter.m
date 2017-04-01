%plot energies

clear 
[File,Path]=uigetfile('*.txt','MultiSelect','off');
 
 str=sprintf('%s', [Path File]);                            %makes str be the name of file (along with its path)
 format shortG                                              %change formating so doesn't show 0's for e-11 values. 
 data= load (str);                                          %load the .txt file into matrix called "data"
 
 time = data(:,1);
 total_kinetic = data(:,2);
 total_potential = data(:,3);
 total_energy = data(:,4);
 average_KE = mean(total_kinetic)
 average_PE = mean(total_potential)
 average_total_energy = mean(total_energy)
 KE_amplitude = 0.5*(max(total_kinetic)-min(total_kinetic))
 PE_amplitude = 0.5*(max(total_potential)-min(total_potential))
 total_energy_amplitude = 0.5*(max(total_energy)-min(total_energy))
 
 rel_KE_change = KE_amplitude/average_KE
 rel_PE_change = PE_amplitude/abs(average_PE)
 rel_total_change = total_energy_amplitude/abs(average_total_energy)
 
 
 %initial_KE = data(1,2)
 %initial_PE = data(1,3)
 %initial_total_energy = data(1,4)
 
 h = plot(time, total_kinetic);   
 set(h,'LineWidth',1.5);                              
 hold on     
 xlabel({'t (years)'});
 ylabel({'Energy (arb. units'});
 plot(time, total_potential);
 plot(time, total_energy);
 hold off