function GC_olshanii_v2(varargin)
%Copyright Dmitry Yampolsky
%October 2016

%Creating entanglement using integrals of motion
%Maxim Olshanii, Thibault Scoquart, Dmitry Yampolsky, Vanja Dunjko, and Steven Glenn Jackson
%Phys. Rev. A 97, 013630 – Published 26 January 2018

%Version 2 uses epsilon to define mass distribution and doesnt need wrapper
%to exhaust this variable

%MS stands for mass spectrum
%MS_r - MS_realization
%CofV - center of velocity

digits(32);

if nargin == 1
    N_balls = varargin{1};
else
    N_balls =4;
end

GC_init(); %positions, velocities,time

MS = rand(1,N_balls);% to be specified;
MS_size = 3000;
MS = ones(1,N_balls);

N_MS_runs = 50;
rightmost_mass = 1;

MS_base = (N_balls * (N_balls + 1))./([1:N_balls].*([1:N_balls] + 1));

% moved epsilon init
visualize_switch = true;
vis_timescale=.5;
y_positions = zeros(1,N_balls);%for visualization
radii = ones(N_balls,1)*.02;%for visualization
circ_colors = rand(N_balls,3);

n_steps=10000000;

presision_vpa=false;

if (presision_vpa)
    precisionF = @vpa;
else
    precisionF =@(x) x;
end

if false
    my_mass_flunc_F =@(varargin) dpoissrnd;
else
    my_mass_flunc_F =@(varargin) unifrnd(varargin{1},varargin{2});
end

if visualize_switch
    theFigure=figure('Position',[1 500 1600 500]);
    theAxes = axes;
    axis(theAxes,[-5 x_position(end)+N_balls*3 0 1]);
    axis equal
    axis off

    xlabel({sprintf('Rightmost mass = %d',rightmost_mass),...% num2str(rightmost_mass),...
        sprintf('N balls  = %d',N_balls)},'FontSize',16);
    hold on
end

current_drawingH=[];
collision_mark_handle=[];
collision_at = -1;

%main loop
loop_time=0;
start_time = cputime;

if false
    profile ON
end

eps_range_steps = 50;
eps_range = linspace(0.001, 0.2,eps_range_steps);

GC_init_handle = @()GC_init();
energy_handle = @()energy();
qS_handle = @()qS();
vizualization_handle = @()vizualization();

for Neps_ctr = 1:length(eps_range)%different epsilons

    start_time = cputime;
    next_collision_dT_index=0;


    t_stops=[];
    q=[];
    MS_history=zeros(N_MS_runs,N_balls);
    v_stops=[];
    x_stops=[];
    CofVarray=[];
    Nsteps=[];

    eps_ctr =  eps_range(Neps_ctr);
    std_Eps =  eps_ctr / sqrt(3);

    fprintf(['eps = ' num2str(std_Eps) '\n' ]);

    dmn = std_Eps*sqrt(MS_base);
    distribution_radius=sqrt(3)*dmn;

    distribution_interval_min = MS_base - distribution_radius;%used by my_mass_flunc_F
    distribution_interval_max = MS_base + distribution_radius;

    for MS_ctr = 1:N_MS_runs

        x_position = (1:N_balls);
        x_position = (x_position .* (x_position+1) * .5).^.5;
        v_velo = zeros(1,N_balls);
        v_velo(end) = -1;% last ball is headed towards the rest at speed of 1. This is the energy source

        MS_r = my_mass_flunc_F(distribution_interval_min,distribution_interval_max);

        MS_history(MS_ctr,:) = MS_r;

        offCoM=mean(MS_r(1:end-1));

        %propogation loop
        stuck_iteration_indicator = false;

        for step_ctr = 1:n_steps% for explicit calculation steps, only

            v = (([0 v_velo(1:end-1)]) - (v_velo));
            dT = ((x_position - [0 x_position(1:end-1)])./ v);

            if prod(logical(dT<=0.0 | isinf(dT)))
                %all balls have been "launched"
                break
            end

            dT(dT<=0) = Inf;

            [next_collision_dT, next_collision_dT_index]  = min(dT);
            collision_at=next_collision_dT_index;%for vis

            x_position_tmp = precisionF((x_position) + (v_velo) .* (next_collision_dT));
            x_position = x_position_tmp;%positions are replaced with velocity changes below, instead of being derived here, to handle precision at zero

            if false
                fprintf([' ' num2str(next_collision_dT_index) ' '  num2str(x_position(next_collision_dT_index)) '\n']);
            end

            %Velocity changes for collided balls.
            %They are currently in contact.
            if next_collision_dT_index == 1

                x_position(next_collision_dT_index) =  0.0;
                v_velo(next_collision_dT_index) = -v_velo(next_collision_dT_index);

            else

                x_position(next_collision_dT_index) =   x_position(next_collision_dT_index-1);%****
                v_velo_tmp=v_velo;

                if true

                    v_velo_tmp(next_collision_dT_index-1) =((MS_r(next_collision_dT_index-1) .* v_velo(next_collision_dT_index-1) -...
                        MS_r(next_collision_dT_index) .* v_velo(next_collision_dT_index-1) +...
                        2 * MS_r(next_collision_dT_index) .* v_velo(next_collision_dT_index)) ./ ...
                        (MS_r(next_collision_dT_index-1) + MS_r(next_collision_dT_index)));

                    v_velo_tmp(next_collision_dT_index) = ((MS_r(next_collision_dT_index) .* v_velo(next_collision_dT_index) -...
                        MS_r(next_collision_dT_index-1) .* v_velo(next_collision_dT_index) +...
                        2 * MS_r(next_collision_dT_index-1) .* v_velo(next_collision_dT_index-1)) ./ ...
                        (MS_r(next_collision_dT_index-1) + MS_r(next_collision_dT_index)));

                end

                v_velo=v_velo_tmp;

            end

            CofV = sum(v_velo(1:end-1) .* offCoM)./sum(offCoM); %ot be averaged of mass spcetrum

            if visualize_switch
                vizualization_handle();
            end

        end%end propagation

        CofVarray(end+1) = (sum(v_velo(1:end-1) .* offCoM)./sum(offCoM)); %to be averaged of mass spcetrum
        v_stops(end+1,:) = (v_velo(:));
        x_stops(end+1,:)=  (x_position(:));
        Nsteps(end+1)=step_ctr;

    end

    %add new epsilon value realization results to data
    energy_handle();
    qS_handle();

    output_struct(Neps_ctr).N_balls = N_balls;
    output_struct(Neps_ctr).Nsteps = Nsteps;
    output_struct(Neps_ctr).MS_history = MS_history;
    output_struct(Neps_ctr).v_stops = v_stops;
    output_struct(Neps_ctr).CofVarray = CofVarray;
    output_struct(Neps_ctr).x_stops = x_stops;
    output_struct(Neps_ctr).t_stops = t_stops;
    output_struct(Neps_ctr).q = q;
    output_struct(Neps_ctr).std_Eps = std_Eps;

end%all runs all episilon

%if file doesnt extis - save
filename = ['V2_' num2str(N_balls) 'balls' '.mat'];
if ~(exist(filename, 'file') == 2)
    save(filename,'output_struct');
end

%save to ws
assignin('base', ['output_struct_'   num2str(N_balls)  'balls_'  num2str(eps_range_steps-1) 'eps' ], output_struct);

%illustration
if true
    illustration();
end

if false
    profile VIEWER
end

    function GC_init

        t_time_global=0;
        x_position = ([1:N_balls]);
        x_position = (x_position .* (x_position+1) * .5).^.5;
        v_velo = zeros(1,N_balls);
        v_velo(end) = -1;% last ball is headed towards the rest at speed of 1. This is the energy source

    end

    function vizualization

        delete(current_drawingH);

        axis equal
        for drawctr=1:N_balls
            if drawctr==collision_at

                delete(collision_mark_handle);
                collision_mark_handle = plot([x_position(drawctr) x_position(drawctr)],[1,3],'r');
                current_drawingH(drawctr) = viscircles([x_position(drawctr)+.05 0],.2,'EdgeColor',circ_colors(drawctr,:),'DrawBackgroundCircle',false);
                collision_at = -1;
            else
                current_drawingH(drawctr) = viscircles([x_position(drawctr) 0],.2,'EdgeColor',circ_colors(drawctr,:),'DrawBackgroundCircle',false);
            end
        end
        drawnow
        pause(.005);
    end

    function illustration
        figure
        hist_res=25;
        ex=zeros(N_balls,hist_res);%10 - default hist resoluition
        exhist=hist([v_stops(:,1) v_stops(:,2)]);
        for ctr = [1,N_balls] % 1:N_balls
            ex(ctr,:)=hist(v_stops(:,ctr),hist_res);
        end

        [a,b] = hist([v_stops(:,1) ; v_stops(:,end)] ,50);
        [a,~]=hist(v_stops(:,1),b);
        h = bar(b,a,'BarWidth',1,'FaceColor',[0.6000 0.2000 0]);
        hold
        [a,~]=hist(v_stops(:,end),b);
        h = bar(b,a,'BarWidth',1);

        title('Final velocity distribution');

        xlabel({sprintf('Red is wallside ball, blue is rigtmost, rightmost mass = %d',rightmost_mass),...% num2str(rightmost_mass),...
            sprintf('N balls  = %d, N samples = %d',N_balls,N_MS_runs)});


    end

    function energy

        qn_exp=((1+N_balls))./(N_balls.*(1:N_balls).*(1+(1:N_balls)));
        ens=(MS_history.*(v_stops.^2))./2;
        et=((MS_history(:,end).*(1)./2));
        etrep=repmat(et,[1,N_balls]);
        q=ens./etrep;

    end

    function qS

        qsmean=mean(q);
        S = sum(qsmean.*log(1./qsmean));

    end

end
