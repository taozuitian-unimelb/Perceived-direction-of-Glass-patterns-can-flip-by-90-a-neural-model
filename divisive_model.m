% Wavelength_pixel = 20;
% Orientation = 0;
% %g = gabor(Wavelength,Orientation,'SpatialAspectRatio',0.6,'SpatialFrequencyBandwidth',7);
%clear
imageSize_pixel = [300,300]*4; %
imageSize_degee = 30;
imageSize_degree_pixel_raw = imageSize_degee/300;
imageSize_degree_pixel = imageSize_degee/imageSize_pixel(1);
%Wavelength = 30;
fs = 2.16; % cycle/degree
Ts_degree = 1/fs;
SD_degree = 0.16; 
Wavelength_pixel = Ts_degree/imageSize_degree_pixel;
theta = 0;
SD_pixel = SD_degree / imageSize_degree_pixel;
aspect_ratio = 0.6;
phase = pi/2;
SpatialKernel = gabor2D(imageSize_pixel, Wavelength_pixel, theta, SD_pixel, aspect_ratio, phase);
X1 = 525:675;
Y = 525:675;
X2 = X1 - floor(Wavelength_pixel*2.2);  %%%
X3 = X1 + floor(Wavelength_pixel*2.2);
excigarbor1 = SpatialKernel(X1, Y);
excigarbor2 = SpatialKernel(X2, Y);
excigarbor3 = SpatialKernel(X3, Y);
% figure;
% imshow(excigarbor1,[])
% figure;
% imshow(excigarbor2,[])
% figure;
% imshow(excigarbor3,[])
% hold on 
%%
pattern_size = 1200;
center_x = 597;
center_y = 593;
X = 600; Y = 600; %r = [3,6,9];
dot_radius = 2;
canvas = zeros(pattern_size);
r_degrees = [0.5 * Ts_degree, 1*Ts_degree, 2 *Ts_degree];
n = 3;
r_pixels_range =  r_degrees / imageSize_degree_pixel;
theta_gap = pi/6;
for x = 1:pattern_size
    for y = 1:pattern_size
        if (x - center_x)^2 + (y - center_y)^2 <= dot_radius^2
            canvas(y, x) = 1; 
        end
    end
end
endstopping_R_final = [];

for translation_angle = 0:theta_gap:pi/2 
    ratio = translation_angle/(pi/2);      %(r_degrees(n)/Ts_degree)
    canvas_dotpair = canvas;
    r_pixels = r_pixels_range(n);
    X_translated = round(center_x + r_pixels * sin(translation_angle));
    Y_translated = round(center_y + r_pixels * cos(translation_angle));
    % Draw the black dot on the canvas at the translated position
    for x = 1:pattern_size
        for y = 1:pattern_size
            if (x - X_translated)^2 + (y - Y_translated)^2 <= dot_radius^2
                canvas_dotpair(y, x) = -1; % Set the pixel to 1 (black) 
            end
        end
    end
    Canvas_cropped = canvas_dotpair(525:675,525:675);
    if translation_angle > -1
        figure 
        imshow(excigarbor1 + excigarbor2 + excigarbor3, [-1,1])
        %print(gcf,'garbor.png','-dpng')
    end
    %%calresponse
    exci_R1 = max(sum(Canvas_cropped.*excigarbor1,'all')+1,0);
    exci_R2 = max(sum(Canvas_cropped.*excigarbor2,'all')+1,0);
    exci_R3 = max(sum(Canvas_cropped.*excigarbor3,'all')+1,0);
    %endstopping_R = exci_R1 / (3+(exci_R1 + 1*(exci_R2 + exci_R3)));
    endstopping_R = exci_R1 / (1+(exci_R1 + 5*(exci_R2 + exci_R3)));
    %endstopping_R = exci_R - 0.5*(1/(1.5*1.5))*inhi_R; %random dot
    %endstopping_R_final = [endstopping_R_final, max(endstopping_R, 0)];
    endstopping_R_final = [endstopping_R_final, endstopping_R];

end

if n == 1
    rm = max(endstopping_R_final);
end

circulardata = [endstopping_R_final flip(endstopping_R_final) endstopping_R_final flip(endstopping_R_final)];
theta = pi/2-([0:theta_gap:pi/2, pi/2:theta_gap:pi, pi:theta_gap:1.5*pi, 1.5*pi:theta_gap:2*pi]);
figure
polarplot(theta,circulardata,'LineWidth',2);
rlim([0 1])
ax = gca;
ax.ThetaTick = [];
%ax.ThetaTick = ['0°' '90°' '180°' '270°'];
print(gcf,['polar_dotsep_' num2str(round(r_degrees(n),2)) '_contrast.png' ],'-dpng')