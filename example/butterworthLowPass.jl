# Imports
using Reduce
@force using Reduce.Algebra
using Calculus
using Revise

using Plots
gr();

import mna

# Netlist
Z₀ = 50;
Vₛ = 1;

testbench = Dict([ "Vs" => (:V, 1, 0, Vₛ)
                 , "Rs" => (:R, 1, 2, Z₀)
                 , "Rl" => (:R, 3, 0, Z₀)
                 ]);

lowPass = Dict([ "C1" => (:C, 2, 0, 15e-12)
               , "C2" => (:C, 3, 0, 15e-12)
               , "L1" => (:L, 2, 3, 82e-9)
               ]);

netlist = merge(testbench, lowPass);

# Modified Nodal Analysis
x = mna.analyze(netlist);
H = mna.transferFunction(x);
ϕ = (ω) -> angle.(H(ω));
τg = (ω) -> derivative(ϕ, ω);

# Sample over Frequency
samples = 201;
f = range(1, stop=2e9, length=samples);
ω = 2π .* f;
tf = hcat(H.(ω)...);

# Get S-Parameters
S₁₁ = ((tf[2,:] ./ -tf[4,:]) .- Z₀) ./ ((tf[2,:] ./ -tf[4,:]) .+ Z₀);
S₂₁ = (tf[3,:] ./ tf[2,:]) .* (1 .+ S₁₁);

S₁₁dB = 20 .* log10.(abs.(S₁₁));
S₂₁dB = 20 .* log10.(abs.(S₂₁));

S₂₁ps = rad2deg.(angle.(S₂₁));

plot( f, S₂₁dB, lab="S₂₁ [dB]", w=1, minorgrid=true
    , xaxis=("f", (1e7, 3e9), :log10)
    , yaxis=("gain [dB]", (-60, 10)) );
plot!( f, S₁₁dB, lab="S₁₁ [dB]", w=1)

plot( f, S₂₁ps
     , lab="Phase S₂₁ (deg)", w=1
     , yaxis=("Phase S₂₁ [deg]")
     , xaxis=("f", (1e7, 3e9), :log10))

# Get Group Delay
groupDelay = -hcat(τg.(ω)...)[3,:];

plot( f, groupDelay
    , lab="Group Delay [s]", w=1
    , yaxis=("Group Delay [s]")
    , xaxis=("f", (1e7, 3e9), :log10))
