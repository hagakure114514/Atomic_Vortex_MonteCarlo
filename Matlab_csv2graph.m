data = readmatrix('traj_p50um_v1cm.csv');              % ファイルが現在のフォルダーまたは MATLAB パス上のフォルダーにない場合は、filename に絶対パス名または相対パス名を指定してね。

x=data(:,1);
y=data(:,2);
z=data(:,3);
E=data(:,4);
vz=data(:,5);

f1 = figure;                                %traj plot
f2 = figure;                                %Energy plot
f3 = figure;                                %traj plot
f4 = figure;                                %zv plot

M=14134;

figure(f1);
plot(x*1e6,y*1e6);
hold on
xlabel(' x[um] '); ylabel(' y[um] ');
box on
plot(x(1)*1e6,y(1)*1e6, 'ko')
plot(x(M)*1e6,y(M)*1e6, 'ro')
ax = gca;
ax.FontSize = 13;
axis equal;
% xlim([-110 110])
xlim([-400 400])
ylim([-400 400])
xticks([-400 -200 0 200 400])
yticks([-400 -200 0 200 400])


figure(f2);
t=0:50:1000000;
t = t(1 : M);
hold on
xlabel(' Time[ms] '); ylabel(' Total Energy [mK] ');
box on
plot(t*1e-3,abs(E)*1e3);
ax = gca;
ax.FontSize = 13;

figure(f3);
plot(y*1e6,z*1e2);
hold on
xlabel(' y[um] '); ylabel(' z[cm] ');
plot(y(1)*1e6,z(1)*1e2, 'ko')
plot(y(M)*1e6,z(M)*1e2, 'ro')
box on
ax = gca;
ax.FontSize = 13;
% axis equal;
% xlim([-110 110])
xlim([-400 400])
ylim([-30 25])
xticks([-400 -200 0 200 400])
yticks([-26 -20 -10 0 10 20])


figure(f4);
hold on
xlabel(' Time[ms] '); ylabel(' v_z [cm/s] ');
box on
plot(t*1e-3,vz*1e2);
ax = gca;
ax.FontSize = 13;

