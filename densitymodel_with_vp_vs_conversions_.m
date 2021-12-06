function [rho,zieq,t]=densitymodel(Tavg,bdot,rhos,z,model)
% Tavg: 10m temperature in celcius
% bdot: accumulation rate in mwe/yr or (kg/m2/yr)
% rhos: surface density in kg/m3
% z: depth in true_metres 
rhoi=917; rhoc=550; rhow=2000; R=8.314;

if nargin==0
   Tavg=-1.5;
    bdot=(0.1*rhoi)  %around 176
    rhos=400;
    z=(0:200)';
end

if nargin<5
    model=2; %default model is HL-johnsen
elseif ischar(model)
    tmodel=find(strcmpi(model,{'hlj' 'hl' 'lz' 'h' 'nh'}));
    if isempty(tmodel)
        tmodel=find(strcmpi(model,{'hljohnsen' 'herronlangway' 'lizwally' 'helsen' 'nabarroherring'}));
    end
    model=tmodel;
end

Tavg=Tavg+273.15;
T=Tavg;%-2;

switch model
    case 1  %herron-langway with Johnsen et al 2000 corrections.
        %Small corrections to HL model which are not in arthern et al. 2010
        c0=0.85*11*(bdot/rhow)*exp(-10160./(R*Tavg));
        c1=1.15*575*sqrt(bdot/rhow)*exp(-21400./(R*Tavg));
    case 2 %herron-langway 1980 (arthern eq. 2)
        c0=11*(bdot/rhow)*exp(-10160./(R*Tavg));
        c1=575*sqrt(bdot/rhow)*exp(-21400./(R*Tavg));
    case 3 %%li zwally 2004 (arthern eq. 2)
        c0=(bdot/rhoi)*(139.21-0.542*Tavg)*8.36*(273.15-T).^-2.061;
        c1=c0;
    case 4 %Helsen et al. 2008 (arthern eq. 2)
        c0=(bdot/rhoi)*(76.138-0.28965*Tavg)*8.36*(273.15-T).^-2.061;
        c1=c0;
    case 5 %Nabarro-Herring (Arthern et al 2010 eq. 4)
        g=9.82;Ec=60e3;Eg=42.2e3;
        c0=0.07*bdot*g*exp(-Ec./(R*T)+Eg/(R*Tavg));
        c1=0.03*bdot*g*exp(-Ec./(R*T)+Eg/(R*Tavg));
    otherwise
        error('unknown density model')
end


k0=c0./bdot; %~g4
k1=c1./bdot;

%critical depth at which rho=rhoc
zc=(log(rhoc/(rhoi-rhoc))-log(rhos/(rhoi-rhos)))/(k0*rhoi); %g6

ix=z<=zc; %find the z's above and below zc
upix=find(ix); %indices above zc
dnix=find(~ix); %indices below zc

q=zeros(size(z)); %pre-allocate some space for q
q(dnix)=exp(k1*rhoi*(z(dnix)-zc)+log(rhoc/(rhoi-rhoc))); %g7
q(upix)=exp(k0*rhoi*z(upix)+log(rhos/(rhoi-rhos))); %g7
rho=q.*rhoi./(1+q); %[g8]


if nargout>1
    %only calculate this if you want zieq
    tc=(log(rhoi-rhos)-log(rhoi-rhoc))/c0; %age at rho=rhoc [g17]
    t=zeros(size(z)); %pre allocate a vector for age as a function of z
    t(upix)=(log(rhoi-rhos)-log(rhoi-rho(upix)))/c0; % [g16] above zc
    t(dnix)=(log(rhoi-rhoc)-log(rhoi-rho(dnix)))/c1+tc; % [g16] below zc
    zieq=t*bdot/rhoi; %[g15]
end



%want to calc vp now using kohnen

Vice = 3578;
rho_ice = 917;
Vs_ice =  1800; %(Dewart, 1970; Polom et al., 2014)
for i = 1:numel(rho)
    Vp(i) = Vice - 2250 * ((rho_ice/rho(i))-1)^(50/61);

end

%now calc Vs using DIEZ et al 2014
for i = 1:numel(rho)
    Vs(i) = Vs_ice - (950*((rho_ice/rho(i))-1)^(100/117));
end


     fid = fopen('HIC_HL_2_200.txt','w');
     for i = 1:numel(rho)
         fprintf(fid,'%f %f %f %f \n', z(i), Vp(i),rho(i), Vs(i));
     end
      fclose(fid) ;


    if nargout==0
    %if you dont want any of the output arguments then you probably want a plot.
    close all
    plot(rho,z,rhoc,zc,'b.')
    axis ij %flip y-axis
    xlabel('Density (kg/m^3)')
    ylabel('Depth (m)')

   
   clear rho
end

