function [] = qcheck(Nstate,Lstate,Tstate,order,n)

%%Compares original states vs reverse lifted states and post truncated
%%states

% Nstate: Vanilla States
% Lstate: Lifted and returned states (for sanity checks)
% Tstate: Truncated and returned state
% order: order of FFT - May not be necessary

[N,m] = size(Nstate);
Nstate(3,:) = wrapToPi(Nstate(3,:));
Nstate(4,:) = wrapToPi(Nstate(4,:));
if n == 0
    n = N;
    start = 1;
    datit = sprintf('State Comparison:2 groups of %d states showing the normal and truncated return-lifted states',n);
else
    start = n;
    datit = sprintf('State Comparison: Original state and truncated return of state %d',n);
end
cute = real(ifft(Lstate,order,1));
cute(3,:) = wrapToPi(cute(3,:));
cute(4,:) = wrapToPi(cute(4,:));
papa = real(ifft(Tstate,order,1));
papa(3,:) = wrapToPi(papa(3,:));
papa(4,:) = wrapToPi(papa(4,:));


%plotting original data


figure
set(gca,'ColorOrderIndex',1)
t = 1:1:m-1;
for z = start:n
    plot(t,Nstate(z,t),'HandleVisibility','off')
    hold on

end
set(gca,'ColorOrderIndex',1)
%t = 1:25:m-1;
for q =start:n
    plot(t,cute(q,t),'-.','HandleVisibility','off');
    hold on
end
set(gca,'ColorOrderIndex',1)
%t = 1:45:m-1;
 for v = start:n
    plot(t,papa(v,t),':','HandleVisibility','off');
    hold on 

 end
 
%datit = sprintf('State Comparison:3 groups of %d states showing the normal, return-lifted and truncated return-lifted states',n);
title(datit)
legend
set(gca,'ColorOrderIndex',1)

%lets add some sparkle to this plot
t = 1:33:m-1;
for z = start:n
    plot(t,Nstate(z,t),'d')
    hold on

end
set(gca,'ColorOrderIndex',1)
t = 1:56:m-1;
for q =start:n
    plot(t,cute(q,t),'v');
    hold on
end
set(gca,'ColorOrderIndex',1)
t = 1:25:m-1;
 for v = start:n
    plot(t,papa(v,t),'*');
    hold on 

 end
%  pause
 
end

