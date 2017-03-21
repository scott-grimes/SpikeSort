function [output, flow] = nn(input,n)
gap1 = input.*n.l1; %row is output from neuron

for i = 1: length(n.l1)
    gap2(:,i) = gap1(i)+gap1(i)*n.l2; %row is output from neuron. colum is input from neuron number from gap1
end

flow1 = gap2+gap2.*n.l3(1);
flow2 = gap2+gap2.*n.l3(2);
flow = [flow1;flow2];
output1 = max(max(flow1));
output2 = max(max(flow2));
if output1>output2
    output = 1;
else
    output=2;
end