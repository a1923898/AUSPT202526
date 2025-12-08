function plotRoom(roomDimensions,rx,tx,figHandle)
figure(figHandle)
X = [0;roomDimensions(1);roomDimensions(1);0;0];
Y = [0;0;roomDimensions(2);roomDimensions(2);0];
Z = [0;0;0;0;0];
figure;
hold on;
plot3(X,Y,Z,"k",LineWidth=1.5);
plot3(X,Y,Z+roomDimensions(3),"k",LineWidth=1.5);
set(gca,"View",[-28,35]);
for k=1:length(X)-1
    plot3([X(k);X(k)],[Y(k);Y(k)],[0;roomDimensions(3)],"k",LineWidth=1.5);
end
grid on
xlabel("X (m)")
ylabel("Y (m)")
zlabel("Z (m)")

plot3(tx(1,1),tx(1,2),tx(1,3),"go",LineWidth=2)
plot3(tx(2,1),tx(2,2),tx(2,3),"ro",LineWidth=2)
num_rx = size(rx,1);
for rx_id = 1:num_rx
    plot3(rx(rx_id,1),rx(rx_id,2),rx(rx_id,3),"bx",LineWidth=1)
end
end