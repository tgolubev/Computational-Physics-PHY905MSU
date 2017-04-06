%Mercury precession
%This script is specific for analyzing the output from project 3 part g
%(the mercury-sun relativistic system w/ no other planets present)
%It Finds the initial perihelion (position when closest to sun)
%and then compares this with perihelion of the last simulated orbit of
%mercury (last 88days of the simulation). Then the shift in perihelions is
%used to calculate the angle related to precession of mercury.

%Timofey Golubev

clear 
[File,Path]=uigetfile('*.txt','MultiSelect','off');
str=sprintf('%s', [Path File]);                            %makes str be the name of file (along with its path)
format shortG                                              %change formating so doesn't show 0's for e-11 values. 
data= load (str);                                          %load the .txt file into matrix called "data"  
    
planet_index = data(:,2);                                               %this column are the indices which identify the planets 
mercury_indices=find(planet_index==2);                                  %find the values corresponding to mercury (index 1 = sun, index 2 = mercury)
time_step = data(mercury_indices(2),1)-data(mercury_indices(1),1)       %subtract 2 adjacent times corresponding to mercury from data file to find the time step
%RECALL that the time step values are repeated for each planet, which is
%why are using mercury_indices to subtract 2 adjacent times.
                                                                    
steps_per_orbit = round((87.969/365.25)/time_step)   %1 orbit of mercury is 88 days, use round to round to  the nearest integer of steps
 
%  Test that perihelion calculation works for the 1st orbit(since we set the intial
%  perihelion in C++ code
%NOTE:this test is only possible if output the initial orbit's positions to
%the output file!!
%  for i=1:steps_per_orbit
%      x=data(mercury_indices(i),4);
%      y=data(mercury_indices(i),5);
%      z=data(mercury_indices(i),6);
%      initial_radii(i) =  sqrt(x^2+y^2+z^2);
%  end
%  [initial_perihelion, index] = min(initial_radii)    
% %This is being calculated correctly, can verify that these equals the value put into C++ code for initial mercury position
%  x_initial_perihelion = data(mercury_indices(index),4)
%  y_initial_perihelion = data(mercury_indices(index),5)
%  z_initial_perihelion = data(mercury_indices(index),6)
%  initial_perihelion_angle = atand(abs(y_initial_perihelion/x_initial_perihelion))   %verify that initial perihelion angle is 0 
 
total_steps = size(mercury_indices,1);  %size along the 1st dimension (# of rows)
last_orbit_start = total_steps-steps_per_orbit;
j=1;
for i=last_orbit_start:total_steps
    x=data(mercury_indices(i),3);
    y=data(mercury_indices(i),4);
    %z=data(mercury_indices(i),5);        %don't use z because they are 0 %anyway
     
    final_radii(j) = sqrt(x^2+y^2);       %USE j as the index, because need to fill the array from 1:steps_per_orbit. Otherwise if use i, it makes array size=total_steps and 
                                          % most of array is just empty (it starts filling it with element # = last_orbit_start (auto initialized elements to 0!!( so when do
                                          % min(final_radii) it yeilds 0!.
                                                                                    
    j=j+1;
end
[final_perihelion, index]=min(final_radii) %this index it returns is j from final_radii(j) 
perihelion_index = index-1+last_orbit_start;   %j=1 when i=last_orbit_start. 
%I have tested that the indices indeed give the minimum radius.
 
x_final_perihelion = data(mercury_indices(perihelion_index),3);
y_final_perihelion = data(mercury_indices(perihelion_index),4);
%z_final_perihelion = data(mercury_indices(perihelion_index),5);
 
final_peri_test = sqrt(x_final_perihelion^2+y_final_perihelion^2)  %ensure the final perihelion radius is same as initial perihelion we input
 
perihelion_angle = atand(abs(y_final_perihelion/x_final_perihelion));  %atand returns arctan in degrees 
angles_in_arcsec = perihelion_angle*3600
 