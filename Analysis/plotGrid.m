A =imread('MMchambercrop.png');
success = [2 2 3 1 3 0 -2 -3;1 0 0 -1 -3 -4 -2 -2];
bone = [2 3 4 2 3; 2 1 0 -2 -5];
pixels= size(A,1);

figure
subplot(1,2,1)
image(A)
axis square
subplot(1,2,2)
hold on
plot(success(1,:), success(2,:), 'o')
plot(bone(1,:),bone(2,:),'x')
grid on
ylim([-6 6])
xlim([-6 6])
axis square
h = circle(0,0,7)
ax.YTick = -6:6;
ax.XTick = -6:6;


