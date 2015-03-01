function MCtest()

    % 1 kT/nm =  4.114 pN
    pNperkT = 4.114;
    kTperpN = 1/pNperkT;
    eVperkT = 0.027;
    
    % total number of links
    N = 100;
    % link length
    param_a = 0.5;
    % paper used h = 100kT/L0^2, close to what I have
    %params['ks'] = 100.0/a^2
    % this one is equivalent to 1000pN elastic modulus
    %param_ks = 1000.0*kTperpN/param_a;
    param_ks = 120;
    % persistence length
    param_Lp = 1.6; % nm
    % bending modulus
    param_kb = param_Lp/param_a;
    % total length
    param_L = N*param_a;
    
    % charge per base & link
    eperbase = 0.12;
    eperlink = eperbase*2*param_a; % ~2 bases/nm

    
    % assumes fmt (r,z,V)
    data = csvread('C:\Users\tamas\Dropbox\research\calculations\comsol\biopore.csv');
    
    rs = unique(data(:,1));
    zs = unique(data(:,2));
    Vs = reshape(data(:,3),[numel(rs),numel(zs)]);
    %Vs = fliplr(Vs);
    
    Uwork = griddedInterpolant({rs,zs},Vs);
    
    % now parse into pore opening outline
    nan0 = find(any(isnan(Vs),1),1,'first');
    nan1 = find(any(isnan(Vs),1),1,'last');
    nanmid = find(any(isnan(Vs),2),1,'first');
    r0 = rs(find(isnan(Vs(:,nan0)),1,'first'));
    r1 = rs(find(isnan(Vs(:,nan1)),1,'first'));
    rmid = rs(nanmid);
    z0 = zs(nan0);
    z1 = zs(nan1);
    zmid = 0;
    
    plot([r0+10, r0, rmid, r1, r1+10],[z0, z0, zmid, z1, z1])
    hold on
    plot(-[r0+10, r0, rmid, r1, r1+10],[z0, z0, zmid, z1, z1])
    pp = [];
    
    function u=Utotal(x)
        % vector of displacements
        dx = diff(x);
        % and their lengths
        ds = sqrt(sum(dx.^2,2));
        % stretching contribution to energy
        u = 0.5*param_ks*sum((ds-param_a).^2);
        % unit vector displacements
        ts = dx./repmat(ds,[1,3]);
        % and bending contribution
        u = u - param_kb*sum(sum(ts(1:end-1,:).*ts(2:end,:)));
        % and finally get the work
        r = sqrt(sum(x(:,1:2).^2,2));
        uw = sum(Uwork([r,x(:,3)]));
        % this is in volts
        % now multiply it by effective charge per link
        uw = uw * eperlink;
        % now it is in eV, convert to kT
        uw = uw / eVperkT;
        u = u + uw;
        if (isnan(u))
            u = 1e100;
        end
    end

    nsamp = 1e5;
    burnin = 1e2;
    thin = 4;
    
    % make big 3d array to hold configurations
    Xs = zeros([nsamp, N, 3]);
    Bases = zeros(nsamp,1);
    % and initialize current config
    X = zeros(N,3);
    % to a line at r = 0
    X(:,3) = -(1:N)*param_a + 8;
    % get the energy
    U = Utotal(X);
    % now loop through and sample configuration space
    curind = 1;
    % number since last one stored
    nsince = -burnin*thin*N;
    % accepted and total
    nacc = [0,0,0];
    ntot = [0,0,0];
    
    % random x proposal distribution, based on stretchiness and stuff
    delta = 1.0*sqrt(2/param_ks);
    % random rotation angle scaling
    dtheta = 0.01*sqrt(2/param_kb);
    % crankshaft angle step size, go all out
    dcrank = 2*pi;
    
    % bookkeeping
    t0 = tic();
    tlast = t0;
    fprintf('Starting %d iterations...\n',nsamp);
    
    while (1)
        % propose new configuration
        Xnew = X;
        
        p = rand();
        if (p < 0.33)
            % randomly perturb (displace) a handful of beads, just not the first one
            for k=1:5
                i = randi(N-1);
                xt = Xnew((i+1):end,:);
                xt = xt + 2*delta*repmat(rand(1,3)-0.5,[size(xt,1),1]);
                Xnew((i+1):end,:) = xt;
            end
            mtype = 1;
        elseif (p < 0.66)
            % pick random bead to perturb
            i = randi(N);
            % make global move
            R = rot_rand(dtheta);
            % rotate the rest of the chain. doing it this way to make it easier to
            % keep end fixed
            xt = Xnew((i+1):end,:);
            x0 = repmat(Xnew(i,:),[size(xt,1), 1]);
            Xnew((i+1):end,:) = x0 + (xt-x0)*R;
            mtype = 2;
        else
            % pick two random beads, not the same
            i = randi(N);
            j = i;
            while (i == j)
                j = randi(N);
            end
            % put them in order
            if (i>j)
                t = i;
                i = j;
                j = t;
            end
            % get the axis between them
            axis = Xnew(j,:) - Xnew(i,:);
            % and normalize it
            axis = axis/sqrt(sum(axis.^2));
            % create random rotation matrix along that axis
            R = rot_aa(axis,2*(rand()-0.5)*dcrank);
            % rotate the rest of the chain. doing it this way to make it easier to
            % keep end fixed
            xt = Xnew((i+1):j,:);
            x0 = repmat(Xnew(i,:),[size(xt,1), 1]);
            Xnew((i+1):j,:) = x0 + (xt-x0)*R;
            mtype = 3;
        end
            
        % now calculate and compare energies
        Unew = Utotal(Xnew);
        
        if rand() < exp(U-Unew)
            % accept proposal
            X = Xnew;
            U = Unew;
            nacc(mtype) = nacc(mtype) + 1;
        end

        if nsince == 0
            tlast = tic;
        end
            
        % do we savae current conformation or not?
        if nsince >= thin*N
            Xs(curind,:,:) = X;
            nsince = 0;
            
            [~,Bases(curind)] = min(abs(X(:,3)));
            fprintf('Mean: %0.1f, Std: %0.2f\n',mean(Bases(1:curind)),std(Bases(1:curind)));

            curind = curind + 1
            
            if mod(curind,1) == 0
                %{
                plot3(X(:,1),X(:,2),X(:,3))
                xlim([-10,10])
                ylim(xlim)
                zlim([-10,20]);
                %}
                delete(pp);
                pp = plot(X(:,1),X(:,3));
                xlim([-10,10])
                ylim([-20,10]);
                drawnow
            end
        end
            
        ntot(mtype) = ntot(mtype) + 1;
        nsince = nsince + 1;
        if (curind == nsamp)
            break
        end
    end

    fprintf(1,'Accept ratios: %0.2f\n',100*nacc/ntot);
    fprintf(1,'Done in %0.2f sec\n',(toc-t0));

end

