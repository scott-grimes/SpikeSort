close all; clear all; clc
if 1==1
    clear all
    n.l1 = rand(3,1)*2;
    n.l2 = rand(7,1)*2;
    n.l3 = rand(2,1)*2;
    save('neural_matrix.mat','n')
end


runs = 1000000
num_correct =zeros(1,100);
num_wrong = zeros(1,100);
num_actual_correct = zeros(1,100);
sucess_coef = 1;
punishment_coef_bad = .9; %.9
punishment_coef_good = 1.3;  %1.3
load('neural_matrix.mat')
for i = 1:runs
input = randi(2);
[output flow]=nn(input,n);
actual = input;
teh_error = min(min(abs((actual-flow)./actual)));
if teh_error<(.3 - i*(.001-.2)/(runs/100-1))
      n.l1 = n.l1.*sucess_coef;
      n.l2 = n.l2.*sucess_coef;
      n.l3 = n.l3.*(1/sucess_coef);
num_wrong(i) = 0;
    num_correct(i) = 1;
    if min(min(abs((actual-flow)./actual)))<.01
    num_actual_correct(i) = 1;
    end
else
   [second,first]= find(min(min(abs((actual-flow)./actual)))==abs((actual-flow)./actual));
   %third of first
   %second of second
   n.l1 = n.l1.*punishment_coef_bad;
   n.l1(first) = n.l1(first).*punishment_coef_good;
   n.l2 = n.l2.*punishment_coef_bad;
   n.l2(second) = n.l2(second).*punishment_coef_good;
    num_wrong(i) = 1;
    num_correct(i) = 0;
    num_actual_correct(i) = 0;
end
if randi(1000)==50 && i > 102
   clc
   percent_correct = 100*sum(num_actual_correct)/(length(num_actual_correct));
   fprintf('%4.2i Percent Correct\n%4.2i Percent Complete',round(percent_correct),round(i*100/runs));
end


end
    save('neural_matrix.mat','n')
    plot(1:runs,cumsum(num_correct)./(cumsum(num_correct+num_wrong)))
    figure
    plot(1:runs,cumsum(num_correct),1:runs,cumsum(num_wrong),'r')
    %plot(1:100000,cumsum(num_correct),'blue',1:100000,cumsum(num_wrong),'r')
