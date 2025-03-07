% This script implements the line attractor with noise and generate samples
% that are used to compute SSI and TSI values, and fit GLMs. 
% Employed to construct Figures S2 and S3

clearvars
close all

v_min=5;    %minimum velocity of the range cm/seg
v_max=28;   %maximum velocity of the range cm/seg
v_pasos=100;%number of steps in the velocity range
l_total=150; %cm
t_total=10;  %seg

sigma_min=20; %minimum sigma in firing probability function
sigma_max=60; %maximum sigma in firing probability function
eta_min=0.5;%minimum eta in the range of eta's
eta_max=0.5;%maximum eta in the range of eta's

eta_mean=eta_min+(eta_max-eta_min)/2;

n_target=50;    %choosen neuron

n_bines_sigma=30; %steps from min sigma to max sigma
n_bines_eta = 10;%steps from min eta to max eta (number of repetitions if limits are equal)
ventana_t=100;  %window length in which mean activation is computed
FR_max=50/1000; %maximum activation rate per ms

sigma_list=linspace(sigma_min,sigma_max,n_bines_sigma);
eta_list=linspace(eta_min,eta_max,n_bines_eta);

%% data generation
z_save=[];
z_H0_save=[];
v_save=[];
p_save=[];
t_save=[];
for c1=1:length(sigma_list)
    k2=1;
    for c2=1:length(eta_list)
        [t_disparo, p_disparo, z_t_disparo, z_p_disparo, v, v_n]=compute_field_centres(n_target, eta_list(c2), v_min, v_max, v_pasos);
        p_disparo_media=mean(p_disparo);
        t_disparo_media=mean(t_disparo);

        z=[];
        zAux=[];
        zAux_H0=[];
        vAux=[];
        pAux=[];
        tAux=[];
        for c3=1:length(v)
            for t=1:t_total*1000
                pos=v(c3)*t/1000;
                n=v_n(c3)*t/1000;
                n_H0=v(c3)*t/1000;

                prob_disp = exp(-((n_target-n)^2)/(sigma_list(c1)^2))*FR_max;
                prob_disp_H0 = exp(-((n_target-n_H0)^2)/(sigma_list(c1)^2))*FR_max;
                zAux(t,c3)=(rand<prob_disp)+0;
                zAux_H0(t,c3)=(rand<prob_disp_H0)+0;

                
                vAux(t,c3)=v(c3);
                pAux(t,c3)=pos;
                tAux(t,c3)=t/1000;
                
            end
        end
        zAux(zAux<0)=0;
        zAux_H0(zAux_H0<0)=0;
        k1=1;
        for t=1:ventana_t:size(zAux,1)-ventana_t
            z_save(k1,:,c1,c2)=mean(zAux(t:t+ventana_t-1,:),1);
            z_H0_save(k1,:,c1,c2)=mean(zAux_H0(t:t+ventana_t-1,:),1);
            v_save(k1,:,c1,c2)=v;
            p_save(k1,:,c1,c2)=pAux(t+ventana_t/2,:);
            t_save(k1,:,c1,c2)=tAux(t+ventana_t/2,:);
            k1=k1+1;
        end
    end
end

%% Computation of indices SSI, TSI and TSI_H0 (TSI under pure place cell model)
% Employed to construct Figure S2

n_bines=50; %number of bins to parse space, time and velocity
p_max=l_total;
t_max=t_total;
p_edges=linspace(0,p_max,n_bines+1);
t_edges=linspace(0,t_max,n_bines+1);
v_edges=linspace(0,v_max,n_bines+1);
[N,edges,p_bines]=histcounts(p_save,p_edges);
[N,edges,t_bines]=histcounts(t_save,t_edges);
[N,edges,v_bines]=histcounts(v_save,v_edges);

media_p_bin=[];
media_t_bin=[];
media_p_total=[];
media_t_total=[];
media_p_H0_bin=[];
media_t_H0_bin=[];
media_p_H0_total=[];
media_t_H0_total=[];
z_p=[];
z_t=[];
z_p_H0=[];
z_t_H0=[];
z_v=[];
z_v_H0=[];

p_bin_max=max(p_bines(:));
t_bin_max=max(t_bines(:));
v_bin_max=max(v_bines(:));
for c1=1:length(sigma_list)
    k2=1;
    for c2=1:length(eta_list)
        z_save_aux=z_save(:,:,c1,c2);
        z_H0_save_aux=z_H0_save(:,:,c1,c2);
        for i=1:p_bin_max
            ind=p_bines(:,:,c1,c2)==i;
            media_p_bin(c1,c2,i)=mean(z_save_aux(ind));
            media_p_H0_bin(c1,c2,i)=mean(z_H0_save_aux(ind));
            
            z_save_aux2=z_save_aux;
            z_save_aux2(~ind)=NaN;
            z_p(c1,c2,i,:)=nanmean(z_save_aux2,1);
            if any(isnan(z_p(:)))
                1;
            end
            
            z_save_aux2=z_H0_save_aux;
            z_save_aux2(~ind)=NaN;
            z_p_H0(c1,c2,i,:)=nanmean(z_save_aux2,1);
            

        end
        media_p_total(c1,c2)=mean(media_p_bin(c1,c2,:));
        media_p_H0_total(c1,c2)=mean(media_p_H0_bin(c1,c2,:));
        
        for i=1:t_bin_max
            ind=t_bines(:,:,c1,c2)==i;
            media_t_bin(c1,c2,i)=mean(z_save_aux(ind));
            media_t_H0_bin(c1,c2,i)=mean(z_H0_save_aux(ind));
            
            
            z_save_aux2=z_save_aux;
            z_save_aux2(~ind)=NaN;
            z_t(c1,c2,i,:)=nanmean(z_save_aux2,1);
            
            
            z_save_aux2=z_H0_save_aux;
            z_save_aux2(~ind)=NaN;
            z_t_H0(c1,c2,i,:)=nanmean(z_save_aux2,1);
        end
        media_t_total(c1,c2)=mean(media_t_bin(c1,c2,:));
        media_t_H0_total(c1,c2)=mean(media_t_H0_bin(c1,c2,:));
        
        
        for i=1:v_bin_max
            ind=v_bines(:,:,c1,c2)==i;
            media_v_bin(c1,c2,i)=mean(z_save_aux(ind));
            media_v_H0_bin(c1,c2,i)=mean(z_H0_save_aux(ind));
            

            z_save_aux2=z_save_aux;
            z_save_aux2(~ind)=NaN;
            z_v(c1,c2,i,:)=nanmean(z_save_aux2,1);
            
            
            z_save_aux2=z_H0_save_aux;
            z_save_aux2(~ind)=NaN;
            z_v_H0(c1,c2,i,:)=nanmean(z_save_aux2,1);
        end
        media_t_total(c1,c2)=mean(media_t_bin(c1,c2,:));
        media_t_H0_total(c1,c2)=mean(media_t_H0_bin(c1,c2,:));
        
        
    end
end

lambda=media_p_bin./media_p_total;
SSI=nanmean(lambda.*log2(lambda),3);

lambda=media_t_bin./media_t_total;
TSI=nanmean(lambda.*log2(lambda),3);

lambda=media_t_H0_bin./media_t_H0_total;
TSI_H0=nanmean(lambda.*log2(lambda),3);


figure;hold on
plot(SSI(:),TSI(:),'.')
plot(SSI(:),TSI_H0(:),'.r')
xlabel('SSI')
ylabel('TSI')
legend('linear attractor','pure place cell')



%% Testing if the space vs speed and time vs speed functions are linear or 
% nonlinear. Employed to construct Figure S3

ind_var=1;
z_p_aux=squeeze(z_p(ind_var,:,:,:));
z_t_aux=squeeze(z_t(ind_var,:,:,:));

[X,Y,Z]=meshgrid(1:size(z_p_aux,1),1:size(z_p_aux,2),1:size(z_p_aux,3));

X=permute(X,[2,1,3]);
Y=permute(Y,[2,1,3]);
Z=permute(Z,[2,1,3]);
aux=[X(:),Y(:),Z(:)];

v_aux=aux(:,3);
t_aux=aux(:,2);
p_aux=aux(:,2);

z_t_aux2=z_t_aux(:);
z_p_aux2=z_p_aux(:);

z_v_t=zeros(v_pasos,n_bines);
z_v_p=zeros(v_pasos,n_bines);
for i=1:size(v_aux,1)
    z_v_t(v_aux(i),t_aux(i))=z_v_t(v_aux(i),t_aux(i))+z_t_aux2(i);
    z_v_p(v_aux(i),p_aux(i))=z_v_p(v_aux(i),p_aux(i))+z_p_aux2(i);
    
end

figure;imagesc(z_v_t)
xlabel('tiempo')
ylabel('velocidad')

figure;imagesc(z_v_p)
xlabel('posicion')
ylabel('velocidad')


COM_t=sum(z_v_t.*repmat(1:size(z_v_t,2),size(z_v_t,1),1),2)./sum(z_v_t,2);
COM_t=flipud(COM_t*t_max/n_bines);
bin_2_3=round(length(COM_t)*2/3);
v_para_reg=linspace(v_min,v_max,length(COM_t));

ind_low=1:bin_2_3;
Xreg = v_para_reg(ind_low)';
Yreg = COM_t(ind_low);

mdl_t=fitlm(Xreg,Yreg);
coefs_t = mdl_t.Coefficients.Estimate;
COM_t_pred = coefs_t(1) + coefs_t(2)*v_max;
COM_t_v_max=COM_t(end);
error_t = COM_t_v_max - COM_t_pred;
COM_t_pred_plot = coefs_t(1) + coefs_t(2)*v_para_reg;

color_marker=ones(1,3)*0;

figure
plot(COM_t,v_para_reg,'Color',color_marker)
hold on
plot(COM_t_pred_plot ,v_para_reg, '--k')
xlim([0, t_max])
ylim([v_min, v_max])
xlabel('time (s)')
ylabel('velocity (cm/s)')
title(['observed-predicted = ',num2str(error_t)])


COM_p=nansum(z_v_p.*repmat(1:size(z_v_p,2),size(z_v_p,1),1),2)./nansum(z_v_p,2);
COM_p=flipud(COM_p*p_max/n_bines);
bin_2_3=round(length(COM_p)*2/3);
v_para_reg=linspace(v_min,v_max,length(COM_p));

ind_low=1:bin_2_3;
Xreg = v_para_reg(ind_low)';
Yreg = COM_p(ind_low);


mdl_p=fitlm(Xreg,Yreg);
coefs_p = mdl_p.Coefficients.Estimate;
COM_p_pred = coefs_p(1) + coefs_p(2)*v_max;
COM_p_v_max=COM_p(end);
error_p = COM_p_v_max - COM_p_pred;
COM_p_pred_plot = coefs_p(1) + coefs_p(2)*v_para_reg;

figure

plot(COM_p,v_para_reg,'Color',color_marker)
hold on
plot(COM_p_pred_plot ,v_para_reg, '--k')
xlim([0, p_max])
ylim([v_min, v_max])
xlabel('space (cm)')
ylabel('velocity (cm/s)')
title(['observed-predicted = ',num2str(error_p)])


%% Fitting GLMs (position, time) model (predictor variable X_pt) and 
% (position,velocity) model (predictor variable X_pv)
% Returns DeltaBIC


X_pt = [p_save(:),t_save(:)];
X_pv = [p_save(:),v_save(:)];
y = z_save(:);

X_pt=(X_pt-min(X_pt,[],1));
X_pt=X_pt./max(X_pt,[],1);

X_pv=X_pv-min(X_pv,[],1);
X_pv=X_pv./max(X_pv,[],1);

y=y*1000;
y=y/max(y)*15;

nIterMax=20000;
statset('MaxIter',nIterMax)
statset('TolX',1e-6)
mdl_pt=fitglm(X_pt,y,'poly55','Distribution','poisson');

BIC_pt=mdl_pt.ModelCriterion.BIC;

statset('MaxIter',nIterMax)
statset('TolX',1e-6)
mdl_pv=fitglm(X_pv,y,'poly55','Distribution','poisson');

BIC_pv=mdl_pv.ModelCriterion.BIC;

mdl_pt.Rsquared
mdl_pv.Rsquared

%delta BIC between the two models. 
DeltaBIC = (BIC_pt-BIC_pv)/(BIC_pt+BIC_pv); 
display(['DeltaBIC = ',num2str(DeltaBIC)])

