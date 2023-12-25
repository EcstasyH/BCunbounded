using Pkg
Pkg.activate(".")

using SumOfSquares
using DynamicPolynomials
using MosekTools
using LinearAlgebra

# Parameters 
ϵ = 10^(-5)
λ = -1
sostol = 4

@polyvar z

# Using Mosek as the SDP solver
solver = optimizer_with_attributes(Mosek.Optimizer, MOI.Silent() => true)

function bc_suf(deg, num_tech)
    # synthesize BC by using the sufficient condition
    # deg: degree of BC template
    # num_tech: whether to use intermediate enhancement techniques (=1) or not (=0)

    # compute init and unsafe region
    mi = length(gi)
    mu = length(gu)
    init = @set(gi[1]>=0)
    for i = 2:mi
        init = intersect(init, @set(gi[i]>=0))
    end
    unsafe = @set(gu[1]>=0)
    for i = 2:mu
        unsafe = intersect(unsafe, @set(gu[i]>=0))
    end
        
    model = SOSModel(solver)
    monos = monomials(vars, 0:deg)
    
    if tech == 1
        # use scaledmonomial basis
        monos = ScaledMonomialBasis(monos)
    end
    
    @variable(model, B, Poly(monos))
    @variable(model, γ)
    dBdt = dot(differentiate(B, vars), f)
    
    @constraint(model, - B - γ >= 0, domain = init, maxdegree =  maxdegree(B)+sostol)
    @constraint(model, B - γ - ϵ >= 0, domain = unsafe, maxdegree =  maxdegree(B)+sostol)
    @constraint(model, λ*B - dBdt - γ >= 0, maxdegree =  maxdegree(dBdt)+sostol)
    @objective(model,Max,γ)
    
    JuMP.optimize!(model)
    if (JuMP.has_values(model)&&SumOfSquares.value(γ)>= 0)
        B_val = SumOfSquares.value(B)
        if tech==1
            coef_list = coefficients(B_val);
            for i in 1:length(coef_list)
                if coef_list[i] <= ϵ && coef_list[i] >= -ϵ
                    coef_list[i] = 0
                end
            end
            B_val = dot(coef_list, monomials(B_val))
        end
        return B_val
    else
        return 0
    end
end

function homo(f)
    # homogenize f w.r.t. variable z (denoted x_0 in paper)
    
    f_homo = 0
    d = maxdegree(f)
    for t in terms(f)
        f_homo += t*z^(d-degree(t))
    end
    return f_homo
end

function bc_nec(deg, tech)
    # synthesize BC by using the necessary condition (homogenization formulation)
    # deg: degree of BC template
    # num_tech: whether to use intermediate enhancement techniques (=1) or not (=0)

    mi = length(gi)
    mu = length(gu)
    init = @set(gi[1]>=0)
    for i = 2:mi
        init = intersect(init, @set(gi[i]>=0))
    end
    unsafe = @set(gu[1]>=0)
    for i = 2:mu
        unsafe = intersect(unsafe, @set(gu[i]>=0))
    end
    
    init_h = @set(homo(gi[1])>=0)
    for i = 2:length(gi)
        init_h = intersect(init_h, @set(homo(gi[i])>=0))
    end
    unsafe_h = @set(homo(gu[1])>=0)
    for i = 2:length(gu)
        unsafe_h = intersect(unsafe_h, @set(homo(gu[i])>=0))
    end

    θ = z^2
    for i = 1:length(vars)
        θ = θ + vars[i]^2
    end

    init_h = intersect(init_h, @set(θ==1));
    unsafe_h = intersect(unsafe_h, @set(θ==1));

    model = SOSModel(solver)
    monos = monomials(vars, 0:deg)
    
    if tech == 1
        monos = ScaledMonomialBasis(monos)    
    end
        
    @variable(model, B, Poly(monos))
    @variable(model, γ)
    dBdt = dot(differentiate(B, vars), f)

    @constraint(model, - homo(B+γ) + ϵ >= 0 , domain = intersect(init_h,@set(z>=0)), maxdegree = maxdegree(B)+sostol)
    @constraint(model,   homo(B-γ) + ϵ >= 0 , domain = intersect(unsafe_h,@set(z>=0)), maxdegree =  maxdegree(B)+sostol)    
    @constraint(model,   homo(λ*B - dBdt - γ) + ϵ >= 0, domain = intersect(@set(θ==1),@set(z>=0)), maxdegree =  maxdegree(dBdt)+sostol)        
    @objective(model, Max, γ)
        
    JuMP.optimize!(model)
    if (JuMP.has_values(model)&&SumOfSquares.value(γ)>=0)
        #println("A feasible solution is found! Optimal Value: ",SumOfSquares.value(γ))
        B_val_h = SumOfSquares.value(B)
        if tech == 1
            coef_list = coefficients(B_val_h);
            for i in 1:length(coef_list)
                if coef_list[i] <= ϵ && coef_list[i] >= -ϵ
                    coef_list[i] = 0
                end
            end
            B_val_h = dot(coef_list, monomials(B_val_h))
        end
        return B_val_h
    else 
        return 0
    end
end

# run all benchmarks with specified BC degrees

tech = 1

benchmarks = ["vector-1", "vector-2", "barrier-1", "barrier-2", "lie-der-1", "lie-der-2", 
    "arch1-1", "arch1-2", "arch2-1", "arch2-2", "arch3-1", "arch3-2", "arch4-1", "arch4-2",
    "nagumo-1", "nagumo-2", "lotka-1", "lotka-2", "lorenz-1", "lorenz-2", "lyapunov-1", "lyapunov-2"]

for name in benchmarks
    include("./Benchmarks/"*name*".jl");
    
    # record system information (used for verification)
    file = open("./Results/systems/"*name*".txt", "w");
    for k = [vars,f,gi,gu]
        write(file, "{")
        for i = 1:length(k)-1
            write(file, string(k[i])*",")
        end
        write(file, string(last(k))*"}\n")
    end
    close(file)

    println(name)
    # sufficient condition
    file = open("./Results/sufficient/"*name*".txt", "w");
    stats = @timed B = bc_suf(bc_deg,tech)
    write(file, Base.replace(string(B),"e"=>"*10^")*"\n")
    write(file, string(stats.time)*"\n")
    println("suf: ", stats.time)
    close(file)
   
    # necessary condition
    file = open("./Results/necessary/"*name*".txt", "w");
    stats = @timed B = bc_nec(bc_deg,tech)
    write(file, Base.replace(string(B),"e"=>"*10^")*"\n")
    write(file, string(stats.time)*"\n")
    println("nec: ", stats.time)
    close(file)
end

