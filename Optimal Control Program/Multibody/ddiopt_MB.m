function output = ddiopt_MB(setup,lambda)

setup.assist.nP = numel(setup.bound.upper.phase);
setup.mesh.flag = 1;
for q=1:setup.assist.nP
    try
        M = legslb(setup.phase(q).mesh_number+1);
        setup.mesh.phase(q).colpoints = setup.phase(q).mesh_points*ones(setup.phase(q).mesh_number,1);
        M = 0.5*(M+1);
        setup.mesh.phase(q).fraction = M(2:end)-M(1:end-1);
    catch
        M = legslb(setup.mesh_number+1);
        setup.mesh.phase(q).colpoints = setup.mesh_points*ones(setup.mesh_number,1);
        M = 0.5*(M+1);
        setup.mesh.phase(q).fraction = M(2:end)-M(1:end-1);
    end
end

if numel(setup.bound.lower.phase(1).velocity)~=numel(setup.bound.lower.phase(1).position)
    error('Number of position variables must equal to the number of velocities')
end

for i=1:setup.assist.nP
    if min(size(setup.bound.lower.phase(i).velocity))~=1
    setup.phase(i).assist.n = 2*min(size(setup.bound.lower.phase(i).velocity));
    setup.phase(i).assist.vv = 1;
    else
    setup.phase(i).assist.n = 2*numel(setup.bound.lower.phase(i).velocity);
    setup.phase(i).assist.vv = 0;
    end
    if min(size(setup.bound.lower.phase(i).position))~=1
        setup.phase(i).assist.pv = 1;
    else
        setup.phase(i).assist.pv = 0;
    end
    if min(size(setup.bound.lower.phase(i).control))~=1
        setup.phase(i).assist.p = min(size(setup.bound.lower.phase(i).control));
        setup.phase(i).assist.cv = 1;
    else
        setup.phase(i).assist.p = numel(setup.bound.lower.phase(i).control);
        setup.phase(i).assist.cv = 0;
    end
try
    setup.phase(i).assist.cp = numel(setup.bound.lower.phase(i).path); 
catch
    setup.phase(i).assist.cp = 0;
end
if isfield(setup.bound.lower.phase(i),'state')
setup.phase(i).assist.w = numel(setup.bound.lower.phase(i).state);
else
    setup.phase(i).assist.w = 0;
end
end
    try
        setup.assist.cn = numel(setup.bound.lower.pconstraints);
    catch
        setup.assist.cn = 0;
    end
try
    setup.assist.v = numel(setup.bound.lower.parameter);
catch
    setup.assist.v = 0;
end

% setup.assist.stP = zeros(1,setup.assist.nP);
% if isfield(setup,'stabilization_constraints')
%     for i=1:numel(setup.stabilization_phases)
%         try
%             setup.assist.stP(setup.stabilization_phases(i)) = 1;
%         catch
%             error(['There are no such phase as ',num2str(setup.stabilization_phases(i)),'.'])
%         end
%     end
% end

global D
[D,setup,I,P] = Differentiation_matrix(setup);
[ Tm,setup ] = sparse_timeM(setup);
xg = initial_guessM(setup);
number = 0;
for q=1:setup.assist.nP
istatepoints = sparse((1:sum(setup.mesh.phase(q).colpoints)+1)'*ones(1,setup.phase(q).assist.n),ones(sum(setup.mesh.phase(q).colpoints)+1,1)*(1:setup.phase(q).assist.n),1+number:setup.phase(q).assist.n*(sum(setup.mesh.phase(q).colpoints)+1)+number);
setup.mesh.phase(q).ipositionpoints = istatepoints(:,1:0.5*size(istatepoints,2));
setup.mesh.phase(q).ivelocitypoints = istatepoints(:,0.5*size(istatepoints,2)+1:size(istatepoints,2));
setup.mesh.phase(q).positionpoints = setup.mesh.phase(q).ipositionpoints(1:end-1,:);
setup.mesh.phase(q).velocitypoints = setup.mesh.phase(q).ivelocitypoints(1:end-1,:);
if not(setup.phase(q).assist.w == 0)
setup.mesh.phase(q).istatepoints = sparse((1:sum(setup.mesh.phase(q).colpoints)+1)'*ones(1,setup.phase(q).assist.w),ones(sum(setup.mesh.phase(q).colpoints)+1,1)*(1:setup.phase(q).assist.w),(1:setup.phase(q).assist.w*(sum(setup.mesh.phase(q).colpoints)+1))+istatepoints(end));
setup.mesh.phase(q).statepoints = setup.mesh.phase(q).istatepoints(1:end-1,:);
setup.mesh.phase(q).controlpoints = sparse((1:sum(setup.mesh.phase(q).colpoints))'*ones(1,setup.phase(q).assist.p),ones(sum(setup.mesh.phase(q).colpoints),1)*(1:setup.phase(q).assist.p),(1:setup.phase(q).assist.p*(sum(setup.mesh.phase(q).colpoints)))+setup.mesh.phase(q).istatepoints(end));
else
setup.mesh.phase(q).controlpoints = sparse((1:sum(setup.mesh.phase(q).colpoints))'*ones(1,setup.phase(q).assist.p),ones(sum(setup.mesh.phase(q).colpoints),1)*(1:setup.phase(q).assist.p),(1:setup.phase(q).assist.p*(sum(setup.mesh.phase(q).colpoints)))+istatepoints(end));
end
setup.mesh.phase(q).initialtimepoint = setup.mesh.phase(q).controlpoints(end)+1;
setup.mesh.phase(q).finaltimepoint = setup.mesh.phase(q).controlpoints(end)+2;
number = setup.mesh.phase(q).finaltimepoint;
if isfield(setup.phase,'stabilization_constraints')
    if isfield(setup.phase,'stabilization_method')
    else
        setup.phase.stabilization_method = 'Differentiation';
    end
end
end
setup.mesh.parameterpoint = (number+1:number+setup.assist.v)';
global x_old dx_old
x_old = rand(numel(xg),1);
dx_old = rand(numel(xg),1);

try
    constraints = Constraints_MRPM(xg,setup);
catch ME
    link = ['<a href="matlab:opentoline(',strrep(ME.stack(1).file,'\','\\'),',',num2str(ME.stack(1).line),',0)">line ',num2str(ME.stack(1).line),'</a>'];
    error(ME.identifier,['Constraints function could not be evaluated because:\n',ME.message,'\nIn ',ME.stack(1).name,' and ',link,'.'])
end
try
    Jfinal = Jacobian_MRPM(xg,setup);
catch ME
    link = ['<a href="matlab:opentoline(',strrep(ME.stack(1).file,'\','\\'),',',num2str(ME.stack(1).line),',0)">line ',num2str(ME.stack(1).line),'</a>'];
    error(ME.identifier,['Jacobian function could not be evaluated because:\n',ME.message,'\nIn ',ME.stack(1).name,' and ',link,'.'])
end
try
    gradient = Gradient_MRPM(xg,setup);
catch ME
    link = ['<a href="matlab:opentoline(',strrep(ME.stack(1).file,'\','\\'),',',num2str(ME.stack(1).line),',0)">line ',num2str(ME.stack(1).line),'</a>'];
    error(ME.identifier,['Gradient function could not be evaluated because:\n',ME.message,'\nIn ',ME.stack(1).name,' and ',link,'.'])
end
try
    objective = Objective_RPM(xg,setup);
catch ME
    link = ['<a href="matlab:opentoline(',strrep(ME.stack(1).file,'\','\\'),',',num2str(ME.stack(1).line),',0)">line ',num2str(ME.stack(1).line),'</a>'];
    error(ME.identifier,['Objective function could not be evaluated because:\n',ME.message,'\nIn ',ME.stack(1).name,' and ',link,'.'])
end

options.lb = [];
options.ub = [];
options.cl = [];
options.cu = [];
for q=1:setup.assist.nP
    if setup.phase(q).assist.pv==0
    for i=1:0.5*setup.phase(q).assist.n
    options.lb = [options.lb;setup.bound.lower.phase(q).initial.position(i)];
    options.ub = [options.ub;setup.bound.upper.phase(q).initial.position(i)];
    options.lb = [options.lb;setup.bound.lower.phase(q).position(i)*ones(sum(setup.mesh.phase(q).colpoints)-1,1)];
    options.ub = [options.ub;setup.bound.upper.phase(q).position(i)*ones(sum(setup.mesh.phase(q).colpoints)-1,1)];
    options.lb = [options.lb;setup.bound.lower.phase(q).final.position(i)];
    options.ub = [options.ub;setup.bound.upper.phase(q).final.position(i)];
    options.cl = [options.cl;zeros(sum(setup.mesh.phase(q).colpoints),1)];
    options.cu = [options.cu;zeros(sum(setup.mesh.phase(q).colpoints),1)];
    end
    else
    for i=1:0.5*setup.phase(q).assist.n
%     options.lb = [options.lb;setup.bound.lower.phase(q).initial.position(i)];
%     options.ub = [options.ub;setup.bound.upper.phase(q).initial.position(i)];
    options.lb = [options.lb;setup.bound.lower.phase(q).position(:,i)];
    options.ub = [options.ub;setup.bound.upper.phase(q).position(:,i)];
%     options.lb = [options.lb;setup.bound.lower.phase(q).final.position(i)];
%     options.ub = [options.ub;setup.bound.upper.phase(q).final.position(i)];
    options.cl = [options.cl;zeros(sum(setup.mesh.phase(q).colpoints),1)];
    options.cu = [options.cu;zeros(sum(setup.mesh.phase(q).colpoints),1)];
    end
    end
    if setup.phase(q).assist.vv==0
    for i=1:0.5*setup.phase(q).assist.n
    options.lb = [options.lb;setup.bound.lower.phase(q).initial.velocity(i)];
    options.ub = [options.ub;setup.bound.upper.phase(q).initial.velocity(i)];
    options.lb = [options.lb;setup.bound.lower.phase(q).velocity(i)*ones(sum(setup.mesh.phase(q).colpoints)-1,1)];
    options.ub = [options.ub;setup.bound.upper.phase(q).velocity(i)*ones(sum(setup.mesh.phase(q).colpoints)-1,1)];
    options.lb = [options.lb;setup.bound.lower.phase(q).final.velocity(i)];
    options.ub = [options.ub;setup.bound.upper.phase(q).final.velocity(i)];
    options.cl = [options.cl;zeros(sum(setup.mesh.phase(q).colpoints),1)];
    options.cu = [options.cu;zeros(sum(setup.mesh.phase(q).colpoints),1)];
    end
    else
    for i=1:0.5*setup.phase(q).assist.n
%     options.lb = [options.lb;setup.bound.lower.phase(q).initial.position(i)];
%     options.ub = [options.ub;setup.bound.upper.phase(q).initial.position(i)];
    options.lb = [options.lb;setup.bound.lower.phase(q).velocity(:,i)];
    options.ub = [options.ub;setup.bound.upper.phase(q).velocity(:,i)];
%     options.lb = [options.lb;setup.bound.lower.phase(q).final.position(i)];
%     options.ub = [options.ub;setup.bound.upper.phase(q).final.position(i)];
    options.cl = [options.cl;zeros(sum(setup.mesh.phase(q).colpoints),1)];
    options.cu = [options.cu;zeros(sum(setup.mesh.phase(q).colpoints),1)];
    end
    end
    for i=1:setup.phase(q).assist.w
    options.lb = [options.lb;setup.bound.lower.phase(q).initial.state(i)];
    options.ub = [options.ub;setup.bound.upper.phase(q).initial.state(i)];
    options.lb = [options.lb;setup.bound.lower.phase(q).state(i)*ones(sum(setup.mesh.phase(q).colpoints)-1,1)];
    options.ub = [options.ub;setup.bound.upper.phase(q).state(i)*ones(sum(setup.mesh.phase(q).colpoints)-1,1)];
    options.lb = [options.lb;setup.bound.lower.phase(q).final.state(i)];
    options.ub = [options.ub;setup.bound.upper.phase(q).final.state(i)];
    options.cl = [options.cl;zeros(sum(setup.mesh.phase(q).colpoints),1)];
    options.cu = [options.cu;zeros(sum(setup.mesh.phase(q).colpoints),1)];
    end
    if setup.phase(q).assist.cv==0
    for i=1:setup.phase(q).assist.p
        options.lb = [options.lb;setup.bound.lower.phase(q).control(i)*ones(sum(setup.mesh.phase(q).colpoints),1)];
        options.ub = [options.ub;setup.bound.upper.phase(q).control(i)*ones(sum(setup.mesh.phase(q).colpoints),1)];
    end
    else
    for i=1:setup.phase(q).assist.p
        options.lb = [options.lb;setup.bound.lower.phase(q).control(:,i)];
        options.ub = [options.ub;setup.bound.upper.phase(q).control(:,i)];
    end
    end
    options.lb = [options.lb;setup.bound.lower.phase(q).initial.time];
    options.ub = [options.ub;setup.bound.upper.phase(q).initial.time];
    options.lb = [options.lb;setup.bound.lower.phase(q).final.time];
    options.ub = [options.ub;setup.bound.upper.phase(q).final.time];
    if isfield(setup.bound.lower.phase(q),'path')
    options.cl = [options.cl;reshape(ones(sum(setup.mesh.phase(q).colpoints),1)*reshape(setup.bound.lower.phase(q).path,1,[]),[],1)];
    options.cu = [options.cu;reshape(ones(sum(setup.mesh.phase(q).colpoints),1)*reshape(setup.bound.upper.phase(q).path,1,[]),[],1)];
    end
    if isfield(setup.phase(q),'stabilization_constraints')
        if not(isempty(setup.phase(q).stabilization_constraints))
        if strcmp(setup.phase(q).stabilization_method,'Integration')
    options.cl = [options.cl;zeros(setup.phase(q).stabilization_constraints*(sum(setup.mesh.phase(q).colpoints)+1),1)];
    options.cu = [options.cu;zeros(setup.phase(q).stabilization_constraints*(sum(setup.mesh.phase(q).colpoints)+1),1)];
        else
    options.cl = [options.cl;zeros(setup.phase(q).stabilization_constraints*(sum(setup.mesh.phase(q).colpoints)),1)];
    options.cu = [options.cu;zeros(setup.phase(q).stabilization_constraints*(sum(setup.mesh.phase(q).colpoints)),1)];
        end
        end
    end
end
try
options.cl = [options.cl;reshape(setup.bound.lower.pconstraints,[],1)];
options.cu = [options.cu;reshape(setup.bound.upper.pconstraints,[],1)];
catch
end
try
options.lb = [options.lb;reshape(setup.bound.lower.parameter,[],1)];
options.ub = [options.ub;reshape(setup.bound.upper.parameter,[],1)];
catch
end

try
    hessian = Hessian_MRPM(xg,12,rand(numel(options.cl),1),setup);
catch ME
    link = ['<a href="matlab:opentoline(',strrep(ME.stack(1).file,'\','\\'),',',num2str(ME.stack(1).line),',0)">line ',num2str(ME.stack(1).line),'</a>'];
    error(ME.identifier,['Hessian function could not be evaluated because:\n',ME.message,'\nIn ',ME.stack(1).name,' and ',link,'.'])
end

%% IPOPT

%   options.ipopt.hessian_approximation = 'limited-memory';
%   options.ipopt.print_level            = 0;
  options.ipopt.max_iter               = 5000;
  options.ipopt.mu_strategy            = 'adaptive';
  options.ipopt.tol                    = 1e-7;
  if nargin >1
      options = initial_multipliers(setup,lambda,options);
  end
%   options.ipopt.derivative_test       = 'first-order';
  Jacstr = Jacobian_MRPM(rand(numel(xg),1),setup)~=0;
  Hesstr = Hessian_MRPM(rand(numel(xg),1),12,rand(numel(options.cl),1),setup)~=0;

  
%   global cont
%   cont = 0;
global flag
if isempty(flag)
    flag.fig = [];
end
if ishandle(flag.fig)
    delete(flag.pl1)
    flag.h1 = newplot(flag.fig);
    flag.x = 0;
    flag.y11 = objective;
    flag.pl1 = line('XData',flag.x,'YData',flag.y11,'Color',[0,0,1],'Parent',flag.h1,'linewidth',2);
%     close(flag.fig)
else
flag.fig = figure('Name','optimization');
% flag.h1 = subplot(2,1,1,'hold','on');
% hold on
flag.h1 = newplot([]);
% flag.pl1 = line('parent',flag.h1,'color',[0,0,0],'marker','o','erase','none', ...
%               'xdata',[],'ydata',[]);
% flag.p11 = plot(x(1),y11(1));
flag.pl1 = plot(0,objective);
% set(flag.h1,'XData',x(1),'YData',y11(1))
% ax = get(flag.pl1,'Parent');
%      hold(ax)
% flag.h2 = subplot(2,1,2,'hold','on');
% flag.p21 = plot(x(1),y21(1));
% flag.p22 = plot(x(1),y22(1));
% set(flag.h1,'XData',x(1),'YData',y11(1))
flag.x = 0;
flag.y11 = objective;
end
  % The callback functions.
  funcs.objective         = @(x) Objective_MRPM(x,setup);
  funcs.constraints       = @(x) Constraints_MRPM(x,setup);
  funcs.gradient          = @(x) Gradient_MRPM(x,setup);
  funcs.jacobian          = @(x) Jacobian_MRPM(x,setup);
  funcs.jacobianstructure = @() Jacstr;
  funcs.hessian           = @(x,sigma,lambda) Hessian_MRPM(x,sigma,lambda,setup);
  funcs.hessianstructure  = @() Hesstr;
  funcs.iterfunc          = @stop_run;
  
  % Run IPOPT.
%   global final_iter
%   final_iter = 0;
%   cleanupObj = onCleanup(@() cleanMeUp(setup));
  [x,info] = ipopt(xg,funcs,options);
  
for q=1:setup.assist.nP
setup.solution.phase(q).time = [double(dtime(x(setup.mesh.phase(q).initialtimepoint),x(setup.mesh.phase(q).finaltimepoint),2,[1,2],setup,q));x(setup.mesh.phase(q).finaltimepoint)];
setup.solution.phase(q).position = x(setup.mesh.phase(q).ipositionpoints);
setup.solution.phase(q).velocity = x(setup.mesh.phase(q).ivelocitypoints);
if setup.phase(q).assist.w ~= 0
    setup.solution.phase(q).state = x(setup.mesh.phase(q).istatepoints);
end
setup.solution.phase(q).control = x(setup.mesh.phase(q).controlpoints);
setup.solution.phase(q).control = [setup.solution.phase(q).control;setup.solution.phase(q).control(end,:)];
end
setup.solution.parameter = x(setup.mesh.parameterpoint);

setup.solution.info = info;
setup.solution.vector = x;
setup.solution.objective = Objective_MRPM(x,setup);
setup.solution.constraints.vector = Constraints_MRPM(x,setup);
cont = 0;
for q=1:setup.assist.nP
    for i=1:0.5*setup.phase(q).assist.n
        setup.solution.constraints.phase(q).position(:,i) = setup.solution.constraints.vector(cont+1+sum(setup.mesh.phase(q).colpoints)*(i-1):cont+sum(setup.mesh.phase(q).colpoints)*i);
        setup.solution.multipliers.constraints.phase(q).position(:,i) = setup.solution.info.lambda(cont+1+sum(setup.mesh.phase(q).colpoints)*(i-1):cont+sum(setup.mesh.phase(q).colpoints)*i);
    end
    setup.solution.multipliers.lower.phase(q).position = setup.solution.info.zl(setup.mesh.phase(q).ipositionpoints);
    setup.solution.multipliers.upper.phase(q).position = setup.solution.info.zu(setup.mesh.phase(q).ipositionpoints);
    cont = cont + sum(setup.mesh.phase(q).colpoints)*i;
    for i=1:0.5*setup.phase(q).assist.n
        setup.solution.constraints.phase(q).velocity(:,i) = setup.solution.constraints.vector(cont+1+sum(setup.mesh.phase(q).colpoints)*(i-1):cont+sum(setup.mesh.phase(q).colpoints)*i);
        setup.solution.multipliers.constraints.phase(q).velocity(:,i) = setup.solution.info.lambda(cont+1+sum(setup.mesh.phase(q).colpoints)*(i-1):cont+sum(setup.mesh.phase(q).colpoints)*i);
    end
    setup.solution.multipliers.lower.phase(q).velocity = setup.solution.info.zl(setup.mesh.phase(q).ivelocitypoints);
    setup.solution.multipliers.upper.phase(q).velocity = setup.solution.info.zu(setup.mesh.phase(q).ivelocitypoints);
    cont = cont + sum(setup.mesh.phase(q).colpoints)*i;
    if setup.phase(q).assist.w~=0
    for i=1:setup.phase(q).assist.w
        setup.solution.constraints.phase(q).state(:,i) = setup.solution.constraints.vector(cont+1+sum(setup.mesh.phase(q).colpoints)*(i-1):cont+sum(setup.mesh.phase(q).colpoints)*i);
        setup.solution.multipliers.constraints.phase(q).state(:,i) = setup.solution.info.lambda(cont+1+sum(setup.mesh.phase(q).colpoints)*(i-1):cont+sum(setup.mesh.phase(q).colpoints)*i);
    end
    setup.solution.multipliers.lower.phase(q).state = setup.solution.info.zl(setup.mesh.phase(q).istatepoints);
    setup.solution.multipliers.upper.phase(q).state = setup.solution.info.zu(setup.mesh.phase(q).istatepoints);
    cont = cont + sum(setup.mesh.phase(q).colpoints)*i;
    end
    setup.solution.multipliers.lower.phase(q).control = setup.solution.info.zl(setup.mesh.phase(q).controlpoints);
    setup.solution.multipliers.upper.phase(q).control = setup.solution.info.zu(setup.mesh.phase(q).controlpoints);
    setup.solution.multipliers.lower.phase(q).time = setup.solution.info.zl([setup.mesh.phase(q).initialtimepoint,setup.mesh.phase(q).finaltimepoint]);
    setup.solution.multipliers.upper.phase(q).time = setup.solution.info.zu([setup.mesh.phase(q).initialtimepoint,setup.mesh.phase(q).finaltimepoint]);
    if isfield(setup.bound.lower.phase(q),'path')
        for i=1:setup.phase(q).assist.cp
        setup.solution.constraints.phase(q).path(:,i) = setup.solution.constraints.vector(cont+1+sum(setup.mesh.phase(q).colpoints)*(i-1):cont+sum(setup.mesh.phase(q).colpoints)*i);
        setup.solution.multipliers.constraints.phase(q).path(:,i) = setup.solution.info.lambda(cont+1+sum(setup.mesh.phase(q).colpoints)*(i-1):cont+sum(setup.mesh.phase(q).colpoints)*i);
        end
        cont = cont + sum(setup.mesh.phase(q).colpoints)*i;
    end
    if isfield(setup.phase(q),'stabilization_constraints')
        if not(isempty(setup.phase(q).stabilization_constraints))
            for i=1:setup.phase(q).stabilization_constraints
        if strcmp(setup.phase(q).stabilization_method,'Integration')
        setup.solution.constraints.phase(q).stabilization(:,i) = setup.solution.constraints.vector(cont+1+(sum(setup.mesh.phase(q).colpoints)+1)*(i-1):cont+(sum(setup.mesh.phase(q).colpoints)+1)*i);
        setup.solution.multipliers.constraints.phase(q).stabilization(:,i) = setup.solution.info.lambda(cont+1+(sum(setup.mesh.phase(q).colpoints)+1)*(i-1):cont+(sum(setup.mesh.phase(q).colpoints)+1)*i);
        else
        setup.solution.constraints.phase(q).stabilization(:,i) = setup.solution.constraints.vector(cont+1+sum(setup.mesh.phase(q).colpoints)*(i-1):cont+sum(setup.mesh.phase(q).colpoints)*i);
        setup.solution.multipliers.constraints.phase(q).stabilization(:,i) = setup.solution.info.lambda(cont+1+sum(setup.mesh.phase(q).colpoints)*(i-1):cont+sum(setup.mesh.phase(q).colpoints)*i);
        end
            end
            cont = cont + sum(setup.mesh.phase(q).colpoints)*i;
        end
    end
end
if setup.assist.cn~=0
    setup.solution.constraints.pconstraints = setup.solution.constraints.vector(cont+1:end);
    setup.solution.multipliers.constraints.pconstraints = setup.solution.info.lambda(cont+1:end);
end
if setup.assist.v~=0
    setup.solution.multipliers.lower.parameters = setup.solution.info.zl(setup.mesh.parameterpoint);
    setup.solution.multipliers.upper.parameters = setup.solution.info.zu(setup.mesh.parameterpoint);
end

output = setup;

function cleanMeUp(setup)

global final_iter
x = final_iter;
for q=1:setup.assist.nP
setup.solution.phase(q).time = [double(dtime(x(setup.mesh.phase(q).initialtimepoint),x(setup.mesh.phase(q).finaltimepoint),2,[1,2],setup,q));x(setup.mesh.phase(q).finaltimepoint)];
setup.solution.phase(q).position = x(setup.mesh.phase(q).ipositionpoints);
setup.solution.phase(q).velocity = x(setup.mesh.phase(q).ivelocitypoints);
if setup.phase(q).assist.w ~= 0
    setup.solution.phase(q).state = x(setup.mesh.phase(q).istatepoints);
end
setup.solution.phase(q).control = x(setup.mesh.phase(q).controlpoints);
setup.solution.phase(q).control = [setup.solution.phase(q).control;setup.solution.phase(q).control(end,:)];
end
setup.solution.parameter = x(setup.mesh.parameterpoint);

assignin('base','final_iter',setup)