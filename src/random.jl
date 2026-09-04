
using Pkg
# Pkg.activate(".")s
using GLMakie
using OrdinaryDiffEq
using Interpolations
using ImageFiltering
using LaTeXStrings
using Optim
using Statistics
using Random

Random.seed!(1234);
include("utils.jl")
include("sim.jl")
include("settings.jl")


algs = [algPulsed0, algPulsed, algChirpPulsed, algSin0, algSin,algChirp]
colorAlgs = [colors[5], colors[7], colors[4], colors[5], colors[7], colors[4]]
###################################################################

# get maximum allowed excitation energy based on the maximum current and coil resistance for ChirpSin algorithm
algChirp.f_TX_Array[1] = ω_nat/ (2*π) 							# initial frequency of excitation (in Hz)
sequence.ψ_TX_Array[1] = 0.0 												# initial phase of excitation (in rad)
sequence.excitationAlg = algChirp
setIandBAmpsForSequence!(sequence, BAmp_Max, 0.0, R, factorTeslaToCurrent, factorCurrentToTesla, algChirp)
deflectionAngleAllFrames, deflectionAngleAllFramesDot = simulate_meas(frameTime, sensor, sequence);
E_max = getEmaxForAlg(BAmp_Max, factorTeslaToCurrent, sequence, algChirp, R)

noSamples = 100
maxDeflectionAnglesDistribution = zeros(length(algs), noSamples)
for i in 1:noSamples
	@info "Running simulation $i of $noSamples"
	for (algIndex, alg) in enumerate(algs) 
		# initial settings for first frame
		alg.f_TX_Array[1] = ω_nat/ (2*π) 							# initial frequency of excitation (in Hz)
		sequence.ψ_TX_Array[1] = 0.0 												# initial phase of excitation (in rad)
		param = "$(shortName(alg))"
		sequence.excitationAlg = alg
		# if isa(alg, ChirpSin_Excitation) # TODO später löschen nur jetzt so schneller
		# 	setIandBAmpsForSequence!(sequence, BAmp_Max, E_max, R, factorTeslaToCurrent, factorCurrentToTesla, algSin)
		# elseif  isa(alg, ChirpRect_Excitation)
		# 	setIandBAmpsForSequence!(sequence, BAmp_Max, E_max, R, factorTeslaToCurrent, factorCurrentToTesla, algPulsed)
		# else
			setIandBAmpsForSequence!(sequence, BAmp_Max, E_max, R, factorTeslaToCurrent, factorCurrentToTesla, sequence.excitationAlg)
		# end


		# run simulation
		deflectionAngleAllFrames, deflectionAngleAllFramesDot = simulate_meas(frameTime, sensor, sequence, useRandom = true);


		maxDeflectionAnglesDistribution[algIndex,i] = rad2deg.(maximum(deflectionAngleAllFrames, dims=1))[end]

	end

end

@info "Max deflection angles distribution over $noSamples runs:"
for (algIndex, alg) in enumerate(algs)
	@info "$(shortName(alg)): mean = $(round(mean(maxDeflectionAnglesDistribution[algIndex,:]), digits=2)) deg, std = $(round(std(maxDeflectionAnglesDistribution[algIndex,:]), digits=2)) deg, max = $(round(maximum(maxDeflectionAnglesDistribution[algIndex,:]), digits=2)) deg, min = $(round(minimum(maxDeflectionAnglesDistribution[algIndex,:]), digits=2)) deg"
end
# make violin plot of the max deflection angles distribution
figViolin = Figure(size = (800, 800))

axViolin = Axis(
    figViolin[1, 1],
    xlabel = "Excitation algorithm",
    ylabel = L"\text{max. def. angle}\,/\,°",
    title = "Distribution of max. deflection angles over $noSamples runs"
)

for (i, alg) in enumerate(algs)
    violin!(
        axViolin,
        fill(i, size(maxDeflectionAnglesDistribution, 2)),
        maxDeflectionAnglesDistribution[i, :],
        color = colorAlgs[i]
    )
end

axViolin.xticks = (1:length(algs), [shortName(alg) for alg in algs])

figViolin

# save figures
save("results/$(numPeriodsForTX)T_MMR_ViolinPlot.png", 	figViolin)
