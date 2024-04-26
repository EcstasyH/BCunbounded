using Pkg
Pkg.activate(".")

using JuMP
using MosekTools
using DynamicPolynomials
using MultivariatePolynomials
using TSSOS
using LinearAlgebra

ϵ = 10^(-5) # ϵ_e in paper
λ = -1      # parameter in exponential type barrier certificate
@polyvar x0 # homogenization variable
@polyvar u
@polyvar v
sos_tol = 1
error = 5

function bc_bound(deg)
    # synthesize BC by using the sound characterization for bounded domains
    # deg: degree of BC template
    
    model = Model(optimizer_with_attributes(Mosek.Optimizer))
    set_optimizer_attribute(model, MOI.Silent(), true)
    
    B, Bc, Bb = add_poly!(model, vars, deg)    
    dBdt = dot(differentiate(B, vars), f)
    d_relax = div(deg+1,2)

    model,info1 = add_psatz!(model, -B, vars, gi, [], div(deg+sos_tol,2), QUIET=true, CS=false, TS=false, Groebnerbasis=true)
    model,info2 = add_psatz!(model, B-ϵ , vars, gu, [], div(deg+sos_tol,2), QUIET=true, CS=false, TS=false, Groebnerbasis=true)
    model,info3 = add_psatz!(model, λ*B-dBdt , vars, [], [], div(maxdegree(λ*B-dBdt)+sos_tol,2), QUIET=true, CS=false, TS=false, Groebnerbasis=true)

    optimize!(model)
    status = termination_status(model)
    Bc = value.(Bc)  
    for i in 1:length(Bc)
        Bc[i] = round(Bc[i]; digits = error)
        # if Bc[i] <= ϵ && Bc[i] >= -ϵ
        #     Bc[i] = 0
        # end
    end
    return Bc'*Bb
end

function bc_complete(deg)
    # synthesize BC by using the complete characterization (polynomial case)
    # deg: degree of BC template

    model = Model(optimizer_with_attributes(Mosek.Optimizer))
    set_optimizer_attribute(model, MOI.Silent(), true)
    
    B, Bc, Bb = add_poly!(model, vars, deg)    
    dBdt = dot(differentiate(B, vars), f)

    model,info1 = add_psatz!(model, -homogenize(B, x0), [x0; vars], [homogenize.(gi, x0); x0], [1-sum([x0;vars].^2)], div(deg+sos_tol,2), QUIET=true, CS=false, TS=false, Groebnerbasis=true)
    model,info2 = add_psatz!(model, homogenize(B-ϵ, x0), [x0; vars], [homogenize.(gu, x0); x0], [1-sum([x0;vars].^2)], div(deg+sos_tol,2), QUIET=true, CS=false, TS=false, Groebnerbasis=true)
    model,info3 = add_psatz!(model, homogenize(λ*B-dBdt, x0), [x0; vars], [x0], [1-sum([x0;vars].^2)], div(maxdegree(λ*B-dBdt)+sos_tol,2), QUIET=true, CS=false, TS=false, Groebnerbasis=true)

    optimize!(model)
    status = termination_status(model)
    Bc = value.(Bc)  
    for i in 1:length(Bc)
        Bc[i] = round(Bc[i]; digits = error)
        # if Bc[i] <= ϵ && Bc[i] >= -ϵ
        #     Bc[i] = 0
        # end
    end
    return Bc'*Bb
end

function bc_completesemi(deg)
    # synthesize BC by using the complete characterization (non-polynomial case)
    # deg: degree of BC template

    model = Model(optimizer_with_attributes(Mosek.Optimizer))
    set_optimizer_attribute(model, MOI.Silent(), true)
    
   
    B1, Bc1, Bb1 = add_poly!(model, vars, deg)    
    B2, Bc2, Bb2 = add_poly!(model, vars, deg)    
    dB1dt = dot(differentiate(B1, vars), f)
    dB2dt = dot(differentiate(B2, vars), f)
    G = λ*(B1+u*B2) - dB1dt - u*dB2dt - v*B2*dot(vars, f)

    model,info1 = add_psatz!(model, -homogenize(B1+u*B2, x0), [x0; vars; u], [homogenize.(gi, x0); x0; u], [u^2-sum([x0;vars].^2), 1-sum([x0;vars;u].^2)], div(deg+1+sos_tol,2), QUIET=true, CS=false, TS=false, Groebnerbasis=true)
    model,info2 = add_psatz!(model, homogenize(B1+u*B2-ϵ, x0), [x0; vars; u], [homogenize.(gu, x0); x0; u], [u^2-sum([x0;vars].^2), 1-sum([x0;vars;u].^2)], div(deg+1+sos_tol,2), QUIET=true, CS=false, TS=false, Groebnerbasis=true)
    model,info3 = add_psatz!(model, homogenize(G, x0), [x0; vars; u; v], [x0;u], [u^2-sum([x0;vars].^2), u*v-x0^2, 1-sum([x0;vars;u;v].^2)], div(maxdegree(G)+sos_tol,2), QUIET=true, CS=false, TS=false, Groebnerbasis=true)

    optimize!(model)
    status = termination_status(model)
    Bc1 = value.(Bc1)  
    for i in 1:length(Bc1)
        Bc1[i] = round(Bc1[i]; digits = error)
        # if Bc1[i] <= ϵ && Bc1[i] >= -ϵ
        #     Bc1[i] = 0
        # end
    end
    Bc2 = value.(Bc2)  
    for i in 1:length(Bc2)
        Bc2[i] = round(Bc2[i]; digits = error)
        # if Bc2[i] <= ϵ && Bc2[i] >= -ϵ
        #     Bc2[i] = 0
        # end
    end
    return Bc1'*Bb1+u*Bc2'*Bb2
end


benchmarks = ["arch1-2", "arch2-1", "arch2-2", "arch3-1", "arch3-2", "arch4-1", "arch4-2",
    "nagumo-1", "nagumo-2", "lotka-1", "lotka-2", "lorenz-1", "lorenz-2"]

for name in benchmarks
    println(name) 
    include("./Benchmarks/"*name*".jl");
    
    # print system
    file = open("./Results/systems/"*name*".txt", "w");
    for k = [vars,f,gi,gu]
        write(file, "{")
        for i = 1:length(k)-1
            write(file, string(k[i])*",")
        end
        write(file, string(last(k))*"}\n")
    end
    close(file)
    
    # print sufficient condition results
    file = open("./Results/sound/"*name*".txt", "w");
    for deg = 1:6
        stats = @timed B = bc_bound(deg)
        write(file, Base.replace(string(B),"e"=>"*10^")*"\n")
        write(file, string(stats.time)*"\n") 
        # println("1:",stats.time)
    end
    close(file)
    
    # print homogenization approach results
    file = open("./Results/complete/"*name*".txt", "w");
    for deg = 1:6
        stats = @timed B = bc_complete(deg)
        write(file, Base.replace(string(B),"e"=>"*10^")*"\n")
        write(file, string(stats.time)*"\n") 
        # println("2:",stats.time)
    end
    close(file)
    
    # print homogenization approach results
    file = open("./Results/completesemi/"*name*".txt", "w");
    for deg = 1:4
        stats = @timed B = bc_completesemi(deg)
        write(file, Base.replace(string(B),"e"=>"*10^")*"\n")
        write(file, string(stats.time)*"\n") 
        # println("3:",stats.time)
    end
    close(file)
end

