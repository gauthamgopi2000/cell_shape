


name = 'coords.csv';
tbl = readtable(name);


%%
m = table2array(tbl);
logfile = readtable('log.csv');
logarray = table2array(logfile);
params = logarray;
param_name = sprintf('KAA_%.2f_KAA_a_%.2f_a_cut_c%.2f_B_%.2f_KA_%.2f_Kp_%.2f_Kb_%.2f',...
    params(1), params(2), params(3), ...
    params(4), params(5), params(6), params(7));

dumpname = strcat('dump',param_name, '.dat') ;
%L=1000.;
box_ls = m(:, 15);
c = m(:,2:12);
Nparticles = m(:,8);
c(:, 16) = m(:, 13);
c(:, 17) = m(:, 14);



index =  1;
start =  1;
Np = Nparticles(start);
last = start + Np-1;
% At = zeros([
% 
% Nt-1,2]);
ts = m(:,1);
total_r = size(ts,1);
Nt = floor(total_r/Np);

energies = zeros([1, Nt]);
Ts = zeros([1,Nt]);
Ntypes = max(m(:,3));
nper = Np/Ntypes;
mean_time_tension = zeros([1,Nt]);

minimum_dist = 100000;
for t =1:Nt


    Np = Nparticles(start);
    box = box_ls(start);
    %box = 70;


    %%%%%calculating overlapped area for packing fraction calculation
    % xs = m(start:last, 4);
    % ys = m(start:last, 5);
    % D0 = 2*m(1,7);
    % Aoverlap = 0;
    % for ii = 1:(Np-1)
    %     for jj = (ii+1):Np
    %         dx = xs(ii) - xs(jj);
    %         dy = ys(ii) - ys(jj);
    %         d = sqrt(dx.^2 + dy.^2);
    %         if(d<D0)
    %             l = (D0-d)/2.;
    %             rad = D0/2;
    %             theta = 2*acos(d/(2*rad));
    %             Ae = 0.5*rad*rad*(theta-sin(theta));
    %             Aoverlap = Aoverlap + Ae;
    %         end
    %     end
    % end

    %%another method for overlap calculation

    xs = m(start:last, 4);
    ys = m(start:last, 5);
    types = m(start:last,3);

    for dd = 1:(length(xs)-1)
        for ee = (dd+1):length(xs)
            if(types(ee)==types(dd))
                continue;
            end
            dxi = pbc_separation(xs(dd)-xs(ee), box_ls(start));
            dyi = pbc_separation(ys(dd)-ys(ee), box_ls(start));
            dst = sqrt(dxi^2+dyi^2);
            if dst<minimum_dist
                minimum_dist = dst;
            end

        end
    end

   
    D0 = 2*m(1,7);
    Aoverlap = 0;
    L_curr = box_ls(start);
    regions = 0;
    dpstart = 1;
    dpend = dpstart+nper-1;
    
    
    for n = 1:Ntypes

        for ii = dpstart:dpend
            jj = ii+1;
            if(ii==dpend)
                jj = dpstart;
            end

            
            dx = xs(ii) - xs(jj);
            %dx = pbc_separation(dx, L_curr);
            dy = ys(ii) - ys(jj);
            % dy = pbc_separation(dy, L_curr);
            d = sqrt(dx.^2 + dy.^2);
            
            if(d<D0)
                l = (D0-d)/2.;
                rad = D0/2;
                theta = 2*acos(d/(2*rad));
                Ae = rad.^2*acos(d/(2*rad)) -0.25*sqrt(d^2*(4*rad.^2-d.^2));
                Aoverlap = Aoverlap + Ae;
                regions = regions+1;
            end
            
        end
        dpstart = dpstart+nper;
        dpend = dpend+nper;
    end






        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%packing fraction calculation
    idstart = start;
    current_m = m(start:last, :);
    apolys = 0;
    atotal = 0;
    for ntype = 1:Ntypes

    
        idend = idstart + nper-1;

        r0s = m(idstart:idend, 7);
        r0 = r0s(1);
        apolys =  apolys + m(idstart,9);
        atotal = atotal + m(idstart,9) + 0.5*pi*r0*r0*(nper+2);
        


        ls_dp = m(idstart:idend,11);
        l0s_dp = m(idstart:idend,12);
        c_curr = c(idstart:idend, 14);

        mean_ls_dp = mean(ls_dp);
        mean_l0s_dp = mean(l0s_dp);
        mean_c_curr = mean(c_curr);

        c(idstart:idend, 12) = mean_ls_dp;
        c(idstart:idend, 13) = mean_l0s_dp;
        c(idstart:idend, 15) = mean_c_curr;
        idstart = idstart+nper;

        
    end
    atotal  =  atotal - Aoverlap;
    L_curr = box_ls(start);
    phi = atotal/L_curr.^2;
    c(start:last, 14) = phi;


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




    xyz = c(start:last,:);
    mean_l_time = mean(c(start:last, 12));
    mean_l0_time = mean(c(start:last, 13));
    mean_time_tension(t) = mean_l_time-mean_l0_time;
    time  = ts(start);
    %xy = xyz(:, 1:2);
    % T = alphaShape(xy, inf);
    % A = area(T);
    % At(t,1) = t;
    % At(t, 2)=A;
    lammps_write(xyz,box,Np,time, dumpname);
    start =  start + Np;
    index =  index + Np;
    
    last = start+Np-1;
end
% plot(1:Nt, mean_time_tension);
% xlabel('\it Time', 'FontSize',20);
% ylabel('\it \langle l - l_0 \rangle', 'FontSize',20);
% xlim([50, 350])


%varname = strcat('k',num2str(k));
%R.(varname)  = [k, At(end,2)/At(end,1)];
%save('rate.mat', "R");
%pgon=polyshape(xyz(:,1),xyz(:,2));