function preprocess(s,t)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Resample first as histogram does not account for weights



histograms = [];


subplot(4,2,1)
% hold off
histograms = [histograms histogram(t.x)];
hax = axis;
hold on
plot([s.x,s.x],[hax(3),hax(4)], 'r')
xlabel('x')
subplot(4,2,2)
% hold off
histograms = [histograms histogram(t.y)];
hax = axis;
hold on
plot([s.y,s.y],[hax(3),hax(4)], 'r')
xlabel('y')
subplot(4,2,3)
% hold off
histograms = [histograms histogram(t.z)];
hax = axis;
hold on
plot([s.z,s.z],[hax(3),hax(4)], 'r')
xlabel('z')
subplot(4,2,4)
% hold off
histograms = [histograms histogram(t.Q)];
hax = axis;
hold on
plot([s.Q,s.Q],[hax(3),hax(4)], 'r')
xlabel('Q')
subplot(4,2,5)
% hold off
histograms = [histograms histogram(t.u)];
hax = axis;
hold on
plot([s.u,s.u],[hax(3),hax(4)], 'r')
xlabel('u')
subplot(4,2,6)
% hold off
histograms = [histograms histogram(t.phi*180/pi)];
hax = axis;
hold on
plot([s.phi*180/pi,s.phi*180/pi],[hax(3),hax(4)], 'r')
xlabel('\phi')
subplot(4,2,7)
% hold off
histograms = [histograms histogram(t.ci)];
hax = axis;
hold on
plot([s.ci,s.ci],[hax(3),hax(4)], 'r')
xlabel('ci')
subplot(4,2,8)
% hold off
histograms = [histograms histogram(t.cii)];
hax = axis;
hold on
plot([s.cii,s.cii],[hax(3),hax(4)], 'r')
xlabel('cii')
end

