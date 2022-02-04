function field_analysis(FOLDER, prefix, res)
%function field_analysis(folder, prefix, res)
%
%Script to analyze fields of meep output.  We assume meep length
%constant of 1 um used always. The user needs to provide the length
%resolution (samples per micron), as well as the folder with meep
%data.
%
%Inputs
%--------
% folder --> folder with meep data
% prefix --> prefix for data file
% res --> spatial resolution (points/um)
%
%Output
%--------
% All data saved in mat file format into FOLDER, with title 'field_data.mat'.

%--------------------
% Analysis Code
%--------------------
tip_file = [FOLDER, '/', prefix, '-tEy_xy.h5'];
bg_file = [FOLDER, '/', prefix, '-bEy_xy.h5'];
output_file = [FOLDER, '/', prefix, '-field_data.mat'];

c = 299792458*1e6*1e-15; %speed of light in microns/femtosecond
S = 0.5; %Courant field factor for the simulation
fcen = 2.6923;
df = 4.6154;
res=300 #Resolution

sample = 15;
y_tip_bottom = 148;%found through analysis of the data...
pad_factor = 5; %Factor for padding when taking the fft.  How many
                %times length of the triangle response to take.

%Load the Ey data for the triangle:
%Data is structured as:
% dim1 --> time
% dim2 --> y
% dim3 --> x
load(tip_file);
ey_t = ey;

%Now we find the maximum field...which should correlate to tip
%apex:
[val loc_x] = max(max(max(ey_t)));
[val loc_y] = max(max(ey_t(:, :, loc_x)));
[val loc_t] = max(ey_t(:, loc_y, loc_x));

%Load the Ey data for the background:
load(bg_file);
ey_b = ey;

%Remove the ey data.
clear ey;

%Get the sizes of the data
%Arranged as time x y_data x x_data
background_size = size(ey_b);
triangle_size = size(ey_t);

%X and Y should be the same:
x = (1:background_size(3))/res; %x in microns
y = (1:background_size(2))/res; %y in microns

%Times are different:
del_T = 1/(fcen + 0.5*df)/sample/c;%Convert sample spacing into
                                   %actual time in fs...
t_background = del_T*(1:background_size(1));
t_triangle = del_T*(1:triangle_size(1));

ey_t_tip = ey_t(:, loc_y, loc_x);
ey_t_tip_norm = ey_t_tip./max(abs(ey_t_tip));

ey_b_tip = ey_b(:,loc_y, loc_x);
ey_b_tip_norm = ey_b_tip./max(abs(ey_b_tip));


%Image the tip field max (just to double-check location
figure(1);
imagesc(x, y, squeeze(abs(ey_t(loc_t, :, :))));
size_ey_t = size(ey_t);
x_marg = ceil(size_ey_t(3)*0.1);
y_marg = ceil(size_ey_t(2)*0.1);
xlim([x(loc_x - x_marg), x(loc_x + x_marg)]);
ylim([y(loc_y - y_marg), y(loc_y + y_marg)]);

%Now plot the field as a function of time:
figure(2);
plot(t_triangle, ey_t_tip_norm, 'linewidth', 3);
hold all;
plot(t_background, ey_b_tip_norm, 'linewidth', 3);
field_enhancement = max(abs(ey_t_tip))./max(abs(ey_b_tip))
xlabel('Time (fs)', 'FontSize', 15);
ylabel('Normalized Field (a.u.)', 'FontSize', 15);


%Prepare variables for frequency analysis:
N = pad_factor*length(ey_t_tip); %Number of terms
T = pad_factor*(t_triangle(end) - t_triangle(1)); %Total timespan
f = (1/T)*[0:N/2]; %Frequency axis
lambda = c./f; %wavelength axis


%Now for frequency analysis with the triangle:
ey_t_tip_pad = zeros(N, 1);
ey_t_tip_pad(1:length(ey_t_tip)) = ey_t_tip;
ey_t_tip_f = fft(ey_t_tip_pad);
ey_t_tip_f = ey_t_tip_f(1:length(f));


%Now for frequency analysis of the background:
ey_b_tip_pad = zeros(N, 1);
ey_b_tip_pad(1:length(ey_b_tip)) = ey_b_tip;
ey_b_tip_f = fft(ey_b_tip_pad);
ey_b_tip_f = ey_b_tip_f(1:length(f));

%Now plot the transfer function:
range = find((lambda > 0.2)&(lambda < 3.0));
TF = ey_t_tip_f(range)./ey_b_tip_f(range);
figure(3);
plot(lambda(range), abs(TF).^2, 'linewidth', 3);
xlabel('Wavelength (micron)');
ylabel('|TF|^2');

figure(4);
plot(lambda(range), unwrap(angle(TF)), 'linewidth', 3);
xlabel('Wavelength (micron)');
ylabel('angle(TF)');


%Save the data:
wavelength = lambda(range);
transfer_function = TF;
tip_field = ey_t_tip;
incident_field = ey_b_tip;
t_tip = t_triangle;
t_incident = t_background;

%Save as MATLAB 7 binary:
save('-7', output_file, 'wavelength', 'transfer_function', ...
	'field_enhancement', 'tip_field', 'incident_field', ...
	't_tip', 't_incident');
