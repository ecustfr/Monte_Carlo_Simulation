
rng();
n           = 0;
nc          = 4;
temperature = 1;
density     = 0.75;
inertia     = 1;
Length      = 0;
r_cut       = 2.5;
dr_max      = 0.15;

nblock      = 10;
nstep       = 100;

n=4*nc^3; %nc：一个方向的晶格数量， 4：晶胞中的原子数(FCC)

box=(n/density)^(1/3);%约化盒子长度

[r,~]=Initialize(n,box,Length,'fcc');%初始化位型

r = r/box;
r = r - round(r); %周期性边界条件

total = Potential_Total_LJ(box, r_cut, r);
assert(~total.ovr,'Overlap in initial configuration')

m_ratio = 0;

%Read_Cnf_Atoms();%读取原子数据
[varibles] = calculate_varibles(box,r_cut,r,n,temperature,total,m_ratio);
%run_begin(); %run_begin就是打印一坨变量

all_varibles.m_r   = 0;
all_varibles.e_c   = 0;
all_varibles.p_c   = 0;
all_varibles.e_f   = 0;
all_varibles.p_f   = 0;
all_varibles.t_c   = 0;
all_varibles.times = 0;

zeros_varibles =   all_varibles;
%all_varibles=over_varibles;
%block 可以理解为切片
disp("Monte Carlo simulation starts:")
disp("-------------------------------")
for blk = 1:nblock
    over_varibles=zeros_varibles; %清空感兴趣的变量
   
    for stp = 1:nstep
        moves = 0;
        for i = 1:n
            rj          = r;
            rj(i,:)     = [];
            partial_old = potential_lj_1( r(i,:), box, r_cut, rj );
            ri_new      = random_translate_vector( dr_max/box, r(i,:));
            ri_new      = ri_new - round(ri_new);
            partial_new = potential_lj_1( ri_new, box, r_cut, rj );
            
            if ~partial_new.ovr
                delta = partial_new.pot - partial_old.pot;
                delta = delta / temperature; 
           
                if metropolis(delta)
                    total  = potential_add(total,partial_new)  ;
                    total  = potential_sub(total,partial_old);
                    r(i,:) = ri_new;
                    moves  = moves+1;
                end
                
            end
        move_ratio    = moves / n;
        varibles      = calculate_varibles(box,r_cut,r,n,temperature,total,move_ratio);
        over_varibles = blk_add( over_varibles, varibles);  
        
        %disp(delta)
        end
    end
     disp(["block:",num2str(blk)])
     disp("-------------------------------")
     all_varibles=blk_end(all_varibles,over_varibles); %把这一坨block 加入到总结果了,同时要展示出来
     %save_configuration
end

run_over(all_varibles);%展示结果并保存数据
total_2 = potential ( box, r_cut, r ) ;

function new = random_translate_vector(unit, old)
    new = (2*rand(1,3)-1)*unit+old;
end

function option = metropolis(delta)
    thresold = 10;
    if delta > thresold
        option = false ; %reject
    elseif delta < 0 
        option = true  ; %accept
    else
        option = (exp(-delta)>rand(1));
    end
end

function [varibles]   = calculate_varibles(box,r_cut,r,n,temperature,total,m_ratio)

vol    = box^3;
rho    = n/vol;
fsq    = force_sq( box, r_cut, r);
%注意区分感兴趣的变量是瞬时变量or not

varibles.m_r    = m_ratio;
%Ideal gas contribution plus cut (but not shifted) PE divided by N
%Internal energy per atom for simulated, cut, potential
varibles.e_c    = 1.5*temperature + total.pot/n ;


%Internal energy per atom for full potential with LRC
%LRC plus ideal gas contribution plus cut (but not shifted) PE divided by N
varibles.e_f    = potential_lrc(rho, r_cut) + 1.5*temperature+total.pot/n;

%Pressure for simulated, cut, potential
%Delta correction plus ideal gas contribution plus total virial divided by V
varibles.p_c    = pressure_delta(rho, r_cut) + rho*temperature + total.vir/vol;

%Pressure for full potential with LRC
%LRC plus ideal gas contribution plus total virial divided by V
varibles.p_f    = pressure_lrc(rho, r_cut)   + rho*temperature + total.vir/vol;

%Configurational temperature
%Total squared force divided by total Laplacian
varibles.t_c    = fsq/total.lap;

varibles.times  = 1;

end

function varibles     = blk_add(varibles_1,varibles_2)
    fields      = fieldnames(varibles_1);
    for i = 1:length(fields)
        k              = fields(i);
        s              = k{1};
        varibles_1.(s) = varibles_1.(s) + varibles_2.(s);
    end
    varibles = varibles_1;
end

function varibles_end = blk_end(all_varibles, over_varibles)
   
    varibles_end = blk_add(all_varibles,over_varibles);
    fields       = fieldnames(all_varibles);
    for i = 1:length(fields)
        k = fields(i);
        s = k{1};    
        if i ~= length(fields)
            disp([s,':' ,num2str(over_varibles.(s)/over_varibles.times)] );
        else
            disp([s,':',num2str(over_varibles.(s))]);
        end
    end
    disp("-------------------------------")
end

function run_over(all_varibles)
disp("across all the simulation")
fields       = fieldnames(all_varibles);
for i = 1:length(fields)
    k = fields(i);
    s = k{1};
    if i~=length(fields)
        disp([s,':' ,num2str(all_varibles.(s)/all_varibles.times),'\t'] );
    else
        disp([s,':',num2str(all_varibles.(s))]);
    end
end
disp("-------------------------------")
disp("Simulation is over")
end

function [total] = Potential_Total_LJ(box ,r_cut ,r)
    [n,d] = size(r);
    total.pot = 0;
    total.vir = 0;
    total.lap = 0;
    total.ovr = false;

    for i=1:n-1
        partial = potential_lj_1(  r(i,:), box, r_cut, r( (i+1:1:end),: ) );
        if partial.ovr
            total.ovr=True;
            break;
        end
    total  = potential_add(total,partial);
    end

end

function partial = potential_lj_1(ri, box, r_cut, r)
    [nj,d]=size(r);
    
    sr2_over     = 1.77; %Overlap threshold (pot>100) ?
    r_cut_box    = r_cut / box;
    r_cut_box_sq = r_cut_box^2;
    box_sq       = box^2;

    rij           = ri-r;
    rij           = rij - round(rij);
    rij_sq        = sum(rij.^2,2);
    in_range      = rij_sq < r_cut_box_sq;
    rij_sq        = rij_sq*box_sq; %Now in sigma=1 units
    sr2           = 1./rij_sq; %(sigma/rij)**2, only if in range
    sr2(~in_range) = 0;
    ovr           = sr2>sr2_over;
    
    if any(ovr) %如果有重叠

        %弹出一个警告
        partial.pot = 0;
        partial.vir = 0;
        partial.lap = 0;
        partial.ovr = any(ovr);
        return;
    end

    sr6  = sr2.^3;
    sr12 = sr6.^2;
    pot  = sr12 - sr6; 
    vir  = pot + sr12;  % LJ pair virials 
    lap  = (22*sr12 - 5*sr6).*sr2;           % LJ pair Laplacians 
    
    partial.pot = sum(pot)*4;
    partial.vir = sum(vir)*24/3;
    partial.lap = sum(lap)*24*2;
    partial.ovr = any(ovr);

    
end

function potential_sum = potential_add(potential_1,potential_2)

    potential_sum.pot  = potential_1.pot + potential_2.pot;
    potential_sum.vir  = potential_1.vir + potential_2.vir;
    potential_sum.lap  = potential_1.lap + potential_2.lap;
    potential_sum.ovr  = potential_1.ovr || potential_2.ovr;
end

function potential_sub = potential_sub(potential_1,potential_2)

    potential_sub.pot  = potential_1.pot - potential_2.pot;
    potential_sub.vir  = potential_1.vir - potential_2.vir;
    potential_sub.lap  = potential_1.lap - potential_2.lap;
    potential_sub.ovr  = potential_1.ovr || potential_2.ovr;
end

function f = force_sq( box, r_cut, r)
    [n,d]=size(r);
    r_cut_box    = r_cut / box;
    r_cut_box_sq = r_cut_box.^2 ;
    box_sq       = box^2;
    f            = zeros(n,d);
    
    for i = 1:n-1
        rij           = r(i,:) - r((i+1:1:n),:); %j>1
        rij           = rij - round(rij);
        rij_sq        = sum(rij.^2,2); % 不用约化单位
        in_range      = rij_sq < r_cut_box_sq;
        rij_sq        = rij_sq * box_sq;
        rij           = rij * box;     % 不用约化单位
        sr2           = 1./rij_sq;
        sr2(~in_range) = 0;
        sr6           = sr2.^3;
        sr12          = sr6.^2;
        fij           = sr2.*(2*sr12 - sr6);
        fij           = rij.*fij;
        f(i,:)        = f(i,:)+sum(fij,1);
        f((i+1:1:n) ,:)      = f((i+1:1:n) ,:)-fij;
    end
    
    f = f*24;
    f = sum(sum(f.^2));
end
    
function p_lrc = potential_lrc( density, r_cut)
    sr3   = 1./r_cut.^3;
    p_lrc = pi*(8/9*sr3.^3 - 8/3*sr3)*density;
end

function pre_lrc = pressure_lrc( density, r_cut)
    sr3     = 1/r_cut.^3;
    pre_lrc = pi*( (32/9)*sr3^3 - 8/3 * sr3) *density^2;
end

function pre_delta = pressure_delta( density, r_cut)
    sr3       = 1/r_cut^3;
    pre_delta = pi*(8/3)*(sr3.^3 - sr3)*density^2; 
end