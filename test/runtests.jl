using mna
using Test

@testset "Butterworth Low Pass Filter (cuttoff = 200MHz)" begin
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

    x = mna.analyze(netlist);
    H = mna.transferFunction(x);

    samples = 201;
    f = range(1, stop=2e9, length=samples);
    ω = 2π .* f;
    tf = hcat(H.(ω)...);

    S₁₁ = ((tf[2,:] ./ -tf[4,:]) .- Z₀) ./ ((tf[2,:] ./ -tf[4,:]) .+ Z₀);
    S₂₁ = (tf[3,:] ./ tf[2,:]) .* (1 .+ S₁₁);

    S₁₁dB = 20 .* log10.(abs.(S₁₁));
    S₂₁dB = 20 .* log10.(abs.(S₂₁));

    fc = 200e6;     # Target: 200MHz
    ftol = 1;       # Tolerance: can be off by 1Hz (depends on sample rate)
    fc_ = first(f[S₂₁dB .< -3]);
    @test (fc_ - 1) <= fc <= (fc_ + 1)
end
