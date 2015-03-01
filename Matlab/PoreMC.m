classdef PoreMC < handle

    properties
        % loaded values from csv
        Rs
        Zs
        Vs
        % gridded interpolant of above
        Vfun
        % pore outline points
        PoreRs
        PoreZs
        % Sim parameters
        Params
        % output arrays
        Xs
        Bases
        % and current values
        X
        U
        Index
        % statistics tracking
        StepsAcc
        StepsTot
    end
    
    properties (Hidden,Constant)
        % 1 kT/nm =  4.114 pN
        pNperkT = 4.114;
        kTperpN = 1/PoreMC.pNperkT;
        eVperkT = 0.027;
    end
    
    methods
        function obj = PoreMC(varargin)
            
            % load precalculated voltages from COMSOL-generated csv
            % assumes fmt (r,z,V)
            data = csvread('C:\Users\szalay\Dropbox\research\calculations\comsol\biopore.csv');

            obj.Rs = unique(data(:,1));
            obj.Zs = unique(data(:,2));
            obj.Vs = reshape(data(:,3),[numel(obj.Rs),numel(obj.Zs)]);
            
            obj.Vfun = griddedInterpolant({obj.Rs,obj.Zs},obj.Vs);
            
            % now parse into pore opening outline
            vn = isnan(obj.Vs);
            nan0 = find(any(vn,1),1,'first');
            nan1 = find(any(vn,1),1,'last');
            nanmid = find(any(vn,2),1,'first');
            r0 = obj.Rs(find(vn(:,nan0),1,'first'));
            r1 = obj.Rs(find(vn(:,nan1),1,'first'));
            rmid = obj.Rs(nanmid);
            z0 = obj.Zs(nan0);
            z1 = obj.Zs(nan1);
            zmid = 0;

            obj.PoreRs = [r0+10, r0, rmid, r1, r1+10];
            obj.PoreZs = [z0, z0, zmid, z1, z1];
            
            % now load/validate all of the parameters
            p = inputParser;
            
            % no required arguments
            % optional sequence argument goes first, no name
            addOptional(p,'sequence',repmat('A',[1 100]),@ischar);
            
            % then all of our named parameter arguments
            % these are the MCMC parameters
            addParameter(p,'samples',1000,@isnumeric);
            addParameter(p,'burnin',100,@isnumeric);
            addParameter(p,'thin',4,@isnumeric);
            % and these are the DNA parameters
            addParameter(p,'linklength',0.5,@isnumeric);
            addParameter(p,'persistence',1.6,@isnumeric);
            addParameter(p,'eperbase',0.12,@isnumeric);
            % this is in kT/nm^2 or such, from Dessinges et al.
            addParameter(p,'kstretch',120,@isnumeric);
            % DNA-pore interaction energy
            addParameter(p,'uinter',[0 -1.5 1 3],@isnumeric);
            % and interaction distance falloff
            addParameter(p,'dinter',0.2,@isnumeric);
            
            parse(p,varargin{:});
            
            obj.Params = p.Results;
            obj.Params.sequence = nt2int(obj.Params.sequence);
            obj.Params.N = numel(obj.Params.sequence);
            obj.Params.eperlink = obj.Params.eperbase*2 ...
                            *obj.Params.linklength; % ~2 bases/nm
            obj.Params.kbend = obj.Params.persistence/obj.Params.linklength;
            
            % make big 3d array to hold configurations
            obj.Xs = zeros([obj.Params.samples, obj.Params.N, 3]);
            obj.Bases = zeros(obj.Params.samples,1);
            % and initialize current config
            obj.X = zeros(obj.Params.N,3);
            % to a line at r = 0
            obj.X(:,3) = -(1:obj.Params.N)*obj.Params.linklength + 7;
            % get the energy
            obj.U = obj.Utotal(obj.X);
            
            obj.Index = 0;
            
            
            % stats
            obj.StepsAcc = [0 0 0];
            obj.StepsTot = [0 0 0];
            
        end
        
        
        function Plot(obj, hax)
            
            if nargin < 2
                figure(1);
                clf
                hax = axes();
            end
            
            fill(obj.PoreRs,obj.PoreZs,[0.5 0.8 0.2],'Parent',hax);
            hold on
            fill(-obj.PoreRs,obj.PoreZs,[0.5 0.8 0.2],'Parent',hax);
            plot(hax,obj.X(:,1),obj.X(:,3),'k');
            xlim(10*[-1 1]);
            ylim([-10 10]);
            daspect([1 1 1])
            hold off
        end
        
        
        function u = Utotal(obj,x)
            
            % vector of displacements
            dx = diff(x);
            % and their lengths
            ds = sqrt(sum(dx.^2,2));
            % stretching contribution to energy
            u = 0.5*obj.Params.kstretch*sum((ds-obj.Params.linklength).^2);
            % unit vector displacements
            ts = dx./repmat(ds,[1,3]);
            % and bending contribution
            u = u - obj.Params.kbend*sum(sum(ts(1:end-1,:).*ts(2:end,:)));
            % and finally get the work
            r = sqrt(sum(x(:,1:2).^2,2));
            uw = sum(obj.Vfun([r,x(:,3)]));
            % this is in volts
            % now multiply it by effective charge per link
            uw = uw * obj.Params.eperlink;
            % now it is in eV, convert to kT
            uw = uw / PoreMC.eVperkT;
            u = u + uw;
            if (isnan(u))
                u = 1e100;
            end
            % now interactions
            dinters = exp(-0.5*(x(:,3)/obj.Params.dinter).^2);
            uinters = obj.Params.uinter(obj.Params.sequence)*dinters;
            u = u + uinters;
            
        end
        
        function SingleStep(obj)
            
            % random x proposal distribution, based on stretchiness and stuff
            delta = 0.2*sqrt(2/obj.Params.kstretch);
            dscale = [1 1 4];
            % random rotation angle scaling
            dtheta = 0.1*sqrt(2/obj.Params.kbend);
            % crankshaft angle step size, go all out
            dcrank = 2*pi;

            % propose new configuration
            Xnew = obj.X;

            % shorthand, cause Matlab
            N = obj.Params.N;

            % pick from various steps
            p = rand();
            if (p < 0.33)
                % randomly perturb (displace) a handful of beads, just not the first one
                for k=1:5
                    i = randi(N-1);
                    xt = Xnew((i+1):end,:);
                    drand = (rand(1,3)-0.5).*dscale;
                    xt = xt + 2*delta*repmat(drand,[size(xt,1),1]);
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
            Unew = obj.Utotal(Xnew);

            if rand() < exp(obj.U-Unew)
                % accept proposal
                obj.X = Xnew;
                obj.U = Unew;
                obj.StepsAcc(mtype) = obj.StepsAcc(mtype) + 1;
            end
            obj.StepsTot(mtype) = obj.StepsTot(mtype) + 1;

        end
        
        function isDone = Next(obj)
            
            % Takes a single output step, consists of N*thin baby steps
            
            isDone = false;
            
            % (except if no burnin yet, does a burnin first)
            if obj.Index == 0
                for i=1:obj.Params.burnin*obj.Params.thin*obj.Params.N
                    obj.SingleStep();
                end
                obj.Index = 1;
            elseif obj.Index > obj.Params.samples
                isDone = true;
                return
            end
            
            for i=1:obj.Params.thin*obj.Params.N
                obj.SingleStep();
            end
            obj.Xs(obj.Index,:,:) = obj.X;
            % weighted avg.
            dinters = exp(-0.5*(obj.X(:,3)/0.5).^2);
            dinters = dinters / sum(dinters);
            obj.Bases(obj.Index) = sum(dinters.*(1:numel(dinters))');
            %[~,obj.Bases(obj.Index)] = min(abs(obj.X(:,3)));
            fprintf('Mean: %0.1f, Std: %0.2f, Accs: %0.2f,%0.2f,%0.2f\n', ...
                mean(obj.Bases(1:obj.Index)),std(obj.Bases(1:obj.Index)), ...
                obj.StepsAcc ./ obj.StepsTot);
            obj.Index = obj.Index + 1;
            
        end

    end
    
end

