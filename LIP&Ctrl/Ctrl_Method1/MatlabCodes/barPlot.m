clear all;clc;close all;
x = [0.4 0.5 0.6 0.7 0.8 0.9 1 1.1 1.2 1.3 1.4 1.5];
y = [435 390;
     430 390;
     420 370;
     410 360;
     390 350;
     400 350;
     390 330;
     380 330;
     380 320;
     370 310;
     370 310;
     360 290];
bar(x, y)
 
ylabel("maximum Impulse (N.s)")
xlabel("velocity (m/s)")
legend("3mass model","1mass model")