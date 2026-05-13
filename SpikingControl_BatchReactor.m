clear all; close all;

set(0,'defaultFigureWindowStyle','docked');
set(0,'defaultTextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(0,'DefaultAxesFontSize', 14);
set(0,'DefaultTextFontSize', 14);

%% System variables
A = [
    1.38     -0.2077  6.715   -5.676
    -0.5814  -4.29    0       0.675  
    1.067    4.273    -6.654  5.893
    0.048    4.273    1.343   -2.104
];

B = [
    0      0
    5.679  0
    1.136  -3.146
    1.136  0
];

C = [
    1 0 1 -1
    0 1 0 0
];

K = [-0.5 -2;5 0.5];

spike_amplitude = [1 4; 3 0.3]./25;
threshold = abs(spike_amplitude./K);
spike_sign = sign(K);

system=ss(A,B,C,0);

N=length(A);
M=2*numel(K);

ix=1:N;
ieta=N+1:N+M;

rng(3);

% Simulation parameters
x0       = rand(4,1)*10;
xi00     = [x0; zeros(M,1)];
Tend     = 10;
Jend     = 5000;

tt=0:1e-4:Tend;
%%
colors=orderedcolors('gem');

[T1, X1, TRIGGERS1] = simulateSpikySystem(system,Tend,Jend,xi00,threshold,spike_amplitude,spike_sign);
disp('Simulated system 1...');
[T2, X2, TRIGGERS2] = simulateSpikySystem(system,Tend,Jend,xi00,threshold./4,spike_amplitude./4,spike_sign);
disp('Simulated system 2...');
[T3, X3, TRIGGERS3] = simulateSpikySystem(system,Tend,Jend,xi00,threshold./15,spike_amplitude./15,spike_sign);
disp('Simulated system 3...');

total_spikes_1=sum(cellfun(@length,TRIGGERS1),'all');
total_spikes_2=sum(cellfun(@length,TRIGGERS2),'all');
total_spikes_3=sum(cellfun(@length,TRIGGERS3),'all');

ssCL=ss(A+B*K*C,B,C,0);
[Y_CL,T_CL,X_CL] = initial(ssCL,x0,tt);


%%
Y1=C*X1(:,ix)';
Y2=C*X2(:,ix)';
Y3=C*X3(:,ix)';

U_CL=K*Y_CL';

X1_norm=vecnorm(X1(:,ix)')';
X2_norm=vecnorm(X2(:,ix)')';
X3_norm=vecnorm(X3(:,ix)')';
X_CL_norm=vecnorm(X_CL')';

gamma = norm(B,2) + integral(@(t)compute_gamma_int(t,A+B*K*C,B),0,Inf);

E1 = max(vecnorm(interp1(T1,X1(:,ix),T_CL)'-X_CL'));
E2 = max(vecnorm(interp1(T2,X2(:,ix),T_CL)'-X_CL'));
E3 = max(vecnorm(interp1(T3,X3(:,ix),T_CL)'-X_CL'));

UB1 = gamma*vecnorm(2*spike_amplitude*ones(size(spike_amplitude,2),1));
UB2 = gamma*vecnorm(2*spike_amplitude./4*ones(size(spike_amplitude,2),1));
UB3 = gamma*vecnorm(2*spike_amplitude./15*ones(size(spike_amplitude,2),1));

fprintf('\n\n');
disp(table([total_spikes_1;total_spikes_2;total_spikes_3],[UB1;UB2;UB3],[E1;E2;E3],'VariableNames',{'Total # spikes','Guaranteed bound','Simulated bound'},'RowNames',{'I','II','III'}));

%%
clf;close all;
figure(1);hold on;axis normal;grid on;
plot(Y1(1,:),Y1(2,:),'LineWidth',1.75,'Color',colors(1,:));
plot(Y2(1,:),Y2(2,:),'LineWidth',1.75,'Color',colors(2,:));
plot(Y3(1,:),Y3(2,:),'LineWidth',1.75,'Color',colors(4,:));
plot(Y_CL(:,1),Y_CL(:,2),'k:','LineWidth',2);
xlabel('$y_1$');
ylabel('$y_2$');
legend('Spiky controller 1','Spiky controller 2','Spiky controller 3','Continuous SOF');
% cleanfigure('targetResolution',1200);


figure(2);
subplot(4,2,1:2);hold on;grid on;
plot(T1,vecnorm(X1(:,ix),2,2),'Color',colors(1,:),'LineWidth',1.75);
plot(T2,vecnorm(X2(:,ix),2,2),'Color',colors(2,:),'LineWidth',1.75);
plot(T3,vecnorm(X3(:,ix),2,2),'Color',colors(4,:),'LineWidth',1.75);
plot(T_CL,vecnorm(X_CL(:,ix),2,2),'k:','LineWidth',2);
xlabel('Time [s]');
ylabel('$|x(t)|$');
legend('Spiky controller 1','Spiky controller 2','Spiky controller 3', 'Continuous SOF');
ax=gca();
ax.YScale='log';

subplot(4,2,3);hold on;grid on;
plot(T1,Y1(1,:),'Color',colors(1,:),'LineWidth',1.75);
plot(T2,Y2(1,:),'Color',colors(2,:),'LineWidth',1.75);
plot(T3,Y3(1,:),'Color',colors(4,:),'LineWidth',1.75);
plot(T_CL,Y_CL(:,1),'k:','LineWidth',2);
xlabel('Time [s]');
ylabel('$y_1(t)$');

subplot(4,2,4);hold on;grid on;
plot(T1,Y1(2,:),'Color',colors(1,:),'LineWidth',1.75);
plot(T2,Y2(2,:),'Color',colors(2,:),'LineWidth',1.75);
plot(T3,Y3(2,:),'Color',colors(4,:),'LineWidth',1.75);
plot(T_CL,Y_CL(:,2),'k:','LineWidth',2);
xlabel('Time [s]');
ylabel('$y_1(t)$');

for ii=1:size(K,1)
    for ij=1:size(K,2)
        subplot(4,2,4+ij);hold on;grid on;xlim([0 10]);
        for ik=1:2
            stem(cell2mat(TRIGGERS1(ii,ij,ik)),ones(length(cell2mat(TRIGGERS1(ii,ij,ik))),1).*spike_amplitude(ii,ij).*(-1)^ik,'Color',colors(1,:),'LineWidth',1.75);
            stem(cell2mat(TRIGGERS2(ii,ij,ik)),ones(length(cell2mat(TRIGGERS2(ii,ij,ik))),1).*spike_amplitude(ii,ij)./4.*(-1)^ik,'Color',colors(2,:),'LineWidth',1.75);
            stem(cell2mat(TRIGGERS3(ii,ij,ik)),ones(length(cell2mat(TRIGGERS3(ii,ij,ik))),1).*spike_amplitude(ii,ij)./15.*(-1)^ik,'Color',colors(4,:),'LineWidth',1.75);
        end
    end
end
subplot(4,2,5);ylabel('$u_1$ (spiky)');xlabel('Time [s]');
subplot(4,2,6);ylabel('$u_2$ (spiky)');xlabel('Time [s]');

subplot(4,2,7);hold on;grid on;
plot(T_CL,U_CL(1,:),'k:','LineWidth',2);
ylabel('$u_1$ (SOF)');xlabel('Time [s]');

subplot(4,2,8);hold on;grid on;
plot(T_CL,U_CL(2,:),'k:','LineWidth',2);
ylabel('$u_2$ (SOF)');xlabel('Time [s]');
% cleanfigure('targetResolution',1200);


function [T, X, TRIGGER] = simulateSpikySystem(sys,Tend,Jend,xi00,threshold,spike_amplitude,spike_sign)
    if ~issystem(sys)
        error("First argument should be a system.");
    end
    sys=ss(sys);
    A=sys.A;
    B=sys.B;
    C=sys.C;

    N=length(A);
    
    dimK = [size(B,2), size(C,1)];
    G=cell([dimK 2]);

    M=2*(dimK(1)*dimK(2));
    
    ix=1:N;
    ieta=N+1:N+M;
    
    Cbar=[];
    for ii = 1:dimK(2)
        Cbar=[Cbar; repmat(C(ii,:),dimK(1),1)];
    end

    F=@(t,xi)[
        A*xi(ix);
        max(0,[Cbar;-Cbar]*xi(ix));
    ];

    for ii=1:dimK(1)
        for ij=1:dimK(2)
            for ik=1:2
                G{ii,ij,ik}=@(t,xi)[
                    xi(ix)-(-1)^ik.*B(:,ii)*spike_amplitude(ii,ij)*spike_sign(ii,ij);
                    reshape(oneij(size(G),sub2ind([dimK 2],ii,ij,ik)),[],1).*xi(ieta);
                ];
            end
        end
    end
    %% Solver options
    options = odeset('RelTol',1e-8,'AbsTol',1e-9,'MaxStep',1e-3,'Events',@(t,xi)event_spike(t,xi,ieta,threshold));

    %% Simulation loop
    xi0=xi00;
    te      = 0;
    T       = zeros(20000000,1);
    X       = zeros(20000000,length(xi00));
    
    TRIGGER=cell([dimK 2]);
    
    jj=0;
    ti=1;
    while te<Tend & jj<Jend
        [t,xi,te,xie,~]     =   ode45(F,[te Tend],xi0,options);
        if ~isempty(te)
            te = te(1);
            xie = xie(1,:);
            vl=sum(t<te);
            T(ti:ti+vl-1)     =   t(t<te);%
            X(ti:ti+vl-1,:)   =   xi(t<te,:);
            ti=ti+vl;
        else
            vl=length(t);
            T(ti:ti+vl-1)     =   t;
            X(ti:ti+vl-1,:)   =   xi;
            ti=ti+vl;
        end
        
        if ~isempty(te)
            xi0=xie';
            while any(xi0(ieta)>=[threshold(:);threshold(:)])
                for ii = 1:dimK(1)
                    for ij = 1:dimK(2)
                        for ik=1:2
                            idx=sub2ind([dimK 2],ii,ij,ik);
                            if xi0(ieta(idx))>=threshold(ii,ij)
                                TRIGGER{ii,ij,ik}  =   [TRIGGER{ii,ij,ik} te];
                                jj=jj+1;
                                xi0         =   G{ii,ij,ik}(te,xi0);
                            end
                        end
                    end
                end
            end
        end
    end
    T=T(1:ti-1);
    X=X(1:ti-1,:);
end

function [value, isterminal, directions] = event_spike(~, xi, ieta, threshold)
    value       = xi(ieta)-[threshold(:); threshold(:)];
    isterminal  = ones(8,1);
    directions   = ones(8,1);
end

function M=oneij(dim,idx)
    M=ones(dim);
    M(idx)=0;
end

function int = compute_gamma_int(t, A_CL, B)
    fun=@(t)norm(A_CL*expm(A_CL*t)*B,2); %Accepts scalar arguments *only*, which is why we need this wrapper function
    int=zeros(size(t));
    parfor ti = 1:length(t)
        int(ti) = fun(t(ti));
    end
end