
using Pkg
Pkg.activate(".")
Pkg.instantiate()

using GLMakie
using OrdinaryDiffEq
using Interpolations
using ImageFiltering
using LaTeXStrings
using Optim
using Statistics


include("utils.jl")
include("sim.jl")
include("settings.jl")


algs = [algPulsed0, algPulsed, algChirpPulsed, algSin0, algSin,algChirp]
colorAlgs = [colors[5], colors[7], colors[4], colors[5], colors[7], colors[4]]
###################################################################
plotFrame = sequence.numFrames

# get maximum allowed excitation energy based on the maximum current and coil resistance for ChirpSin algorithm
algChirp.f_TX_Array[1] = ω_nat/ (2*π) 							# initial frequency of excitation (in Hz)
sequence.ψ_TX_Array[1] = 0.0 												# initial phase of excitation (in rad)
sequence.excitationAlg = algChirp
setIandBAmpsForSequence!(sequence, BAmp_Max, 0.0, R, factorTeslaToCurrent, factorCurrentToTesla, algChirp)
deflectionAngleAllFrames, deflectionAngleAllFramesDot = simulate_meas(frameTime, sensor, sequence);
E_max = getEmaxForAlg(BAmp_Max, factorTeslaToCurrent, sequence, algChirp, R)


figFreqMaxAngle = Figure(size = (800,800))
axFreq = Axis(figFreqMaxAngle[1,1], xlabel = "time / ms", ylabel = "f / Hz", title = "frequency trace")

axMaxAngle = Axis(figFreqMaxAngle[3,1], xlabel = "frame number", ylabel = L"\text{max. def. angle}\,/\,°", title = L"%$(numPeriodsForTX)T_\mathrm{MMR}")
xlims!(axMaxAngle, 1, sequence.numFrames)
ylims!(axMaxAngle,20, 40)

rowIndex = 1
figSignal = Figure(size = (1600,800))
axtop = Axis(figSignal[rowIndex,1:3], xlabel = "time / s", ylabel = L" \varphi\,/\, °", xlabelpadding = -2)
xmin, xmax = 0, frameTime    
xlims!(axtop, xmin, xmax)

axsL = [Axis(figSignal[i+rowIndex, 1],
			xlabel = (i == length(algs) ? "time / ms" : ""),
			ylabel =  L" \dot{\varphi}\,/\,\mathrm{Rad/s}",
			ylabelcolor = colors[1]) for i in 1:length(algs)]
axsR = [Axis(figSignal[i+rowIndex, 1],
				yaxisposition = :right,
			ylabel = L"\mathrm{B} \,/\, \mathrm{µT}",
			ylabelcolor = colors[3]) for i in 1:length(algs)]

axsLMid = [Axis(figSignal[i+rowIndex, 2],
			xlabel = (i == length(algs) ? "time / ms" : ""),
			ylabel =  L" \dot{\varphi}\,/\,\mathrm{Rad/s}",
			ylabelcolor = colors[1]) for i in 1:length(algs)]
axsRMid = [Axis(figSignal[i+rowIndex, 2],
				yaxisposition = :right,
			ylabel = L"\mathrm{B} \,/\, \mathrm{µT}",
			ylabelcolor = colors[3]) for i in 1:length(algs)]

axsLend = [Axis(figSignal[i+rowIndex, 3],
			xlabel = (i == 2 ? "time (ms)" : ""),
			ylabel =  "φ (°)",
			ylabelcolor = colors[1]) for i in 1:length(algs)]
axsRend = [Axis(figSignal[i+rowIndex, 3],
				yaxisposition = :right,
			ylabel = "B / µT",
			ylabelcolor = colors[3]) for i in 1:length(algs)]

window_area = round(Int, 2 * pi * 2 / sensor.ω_nat * sequence.fs) # samples for two periods of the natural frequency

for (algIndex, alg) in enumerate(algs) 
	# initial settings for first frame
	alg.f_TX_Array[1] = ω_nat/ (2*π) 							# initial frequency of excitation (in Hz)
	sequence.ψ_TX_Array[1] = 0.0 												# initial phase of excitation (in rad)
	param = "$(shortName(alg))"
	sequence.excitationAlg = alg
	setIandBAmpsForSequence!(sequence, BAmp_Max, E_max, R, factorTeslaToCurrent, factorCurrentToTesla, sequence.excitationAlg)
	# run simulation
	deflectionAngleAllFrames, deflectionAngleAllFramesDot = simulate_meas(frameTime, sensor, sequence);

	# Calculate excitation-response phase offset at each zero crossing.
	phaseDiff = Float64[]
	txTimes = range(0, sequence.t_TX, step=1/sequence.fs)
	data = deflectionAngleAllFramesDot[1:length(txTimes),plotFrame]
	for i=1:length(data)-1
		if data[i]*data[i+1] < 0
			# Linear interpolation to improve zero crossing estimation
			t1 = txTimes[i]
			t2 = txTimes[i+1]
			u1 = data[i]
			u2 = data[i+1]
			t_crossing = t1 + (0 - u1) * (t2 - t1) / (u2 - u1)

			# calculate phase difference at this time point
			if isa(alg, ChirpSin_Excitation) || isa(alg, ChirpRect_Excitation)
				# for chirp-based algorithms we need to take into account the frequency change over time for the phase calculation
				phaseTX = sequence.ψ_TX_Array[plotFrame] + 2pi * (alg.m[plotFrame]/2*t_crossing^2+alg.f_TX_Array[plotFrame] * t_crossing ) # phase at crossing time point for chirp-based algorithms
			else
				phaseTX = alg.f_TX_Array[plotFrame] * 2pi * t_crossing + sequence.ψ_TX_Array[plotFrame]
			end
			phaseAtCrossing = mod2pi(phaseTX)
			responsePhaseAtCrossing = (u2 - u1) > 0 ? 0.0 : pi # since we are at a zero crossing, the response phase is either 0.0 or pi depending on the direction of crossing
			push!(phaseDiff, atan(sin(phaseAtCrossing - responsePhaseAtCrossing), cos(phaseAtCrossing - responsePhaseAtCrossing)))
		end
	end
	@info "Phase difference for $param: mean = $(rad2deg(round(mean(phaseDiff), digits=2))) deg, std = $(rad2deg(round(std(phaseDiff), digits=2))) deg, max = $(rad2deg(round(maximum(abs.(phaseDiff)), digits=2))) deg"

	##########################################################
  # Plotting 
	##########################################################
	# plot the complete time signals for the chosen frame
  if isa(alg, FAvgRect_Excitation)
		t_Frame = range(0, sequence.numSamplesPerFrame/fs , step=1/sequence.fs)[1:end-1]
		lines!(axtop, t_Frame, rad2deg.(deflectionAngleAllFrames[:,plotFrame]), color = colors[1])
    vspan!(axtop, 0, sequence.t_TX, color = (colors[3], 0.5))
		ylims!(axtop, -40, 40)
		axtop.yticks = [-40, 0, 40]
	end

	# plot zoomed-in views of time signals for the chosen frame
	t_FrameTX = range(0, sequence.t_TX, step=1/sequence.fs)
	phiDotFrame = deflectionAngleAllFramesDot[1:length(t_FrameTX),plotFrame]
	fTX = alg.f_TX_Array[plotFrame]
	ψ_TX = sequence.ψ_TX_Array[plotFrame]
	B = alg.B.(t_FrameTX, sequence.BAmp[plotFrame],fTX, ψ_TX, alg.m[plotFrame], alg.D, sequence.t_TX)

	# start of TX window
	axZoomStart = axsL[algIndex]
	axZoomStart_R = axsR[algIndex]
	configure_zoom_axes!(axZoomStart, axZoomStart_R, algIndex, length(algs), BAmp_Max, "$(shortName(alg))")
	plot_zoom_window!(axZoomStart, axZoomStart_R, t_FrameTX, phiDotFrame, B, 1, 1 + window_area)
	
	# middle of TX window
	axZoomMid = axsLMid[algIndex]
	axZoomMid_R = axsRMid[algIndex]
	configure_zoom_axes!(axZoomMid, axZoomMid_R, algIndex, length(algs), BAmp_Max, "$(shortName(alg))")
	midIndex = round(Int, length(t_FrameTX)/2)
	plot_zoom_window!(axZoomMid, axZoomMid_R, t_FrameTX, phiDotFrame, B, midIndex, midIndex + window_area)

	# end of TX window
	axZoomEnd = axsLend[algIndex]
	axZoomEnd_R = axsRend[algIndex]
	configure_zoom_axes!(axZoomEnd, axZoomEnd_R, algIndex, length(algs), BAmp_Max, "$(shortName(alg))")
	endIndex = length(t_FrameTX)
	plot_zoom_window!(axZoomEnd, axZoomEnd_R, t_FrameTX, phiDotFrame, B, endIndex - window_area + 1, endIndex)

  if isa(alg, ChirpSin_Excitation)
		txTimes = range(0, sequence.t_TX, step=1/sequence.fs)
		zeroCrossings = Float64[]
		data = deflectionAngleAllFrames[1:length(txTimes),plotFrame]
		# Estimate the instantaneous response frequency from successive zero crossings.
		for i=1:length(data)-1
			if data[i]*data[i+1] < 0
				# Linear interpolation to improve zero crossing estimation
				t1 = txTimes[i]
				t2 = txTimes[i+1]
				u1 = data[i]
				u2 = data[i+1]
				t_crossing = t1 + (0 - u1) * (t2 - t1) / (u2 - u1)
				push!(zeroCrossings, t_crossing)
			end
		end
		times_ = zeroCrossings[1:end-1] .+ diff(zeroCrossings) / 2
		freqs_ = 1 ./ diff(zeroCrossings) / 2
		freqs_ = imfilter(freqs_, Kernel.gaussian((4,)));
		freqs = extrapolate(interpolate((times_,), freqs_, Gridded(Linear())), Linear())
		times_inter = range(times_[1], times_[end], length=sequence.numTxSamples)
		freqsInterp = freqs(times_inter)
    lines!(axFreq, times_inter .* 1e3 , freqsInterp, color=colors[1], label = "MMR")
    freqEstimChirp(m,t) = m*t + algChirp.f_TX_Array[plotFrame]
    lines!(axFreq, times_inter .* 1e3, freqEstimChirp.(algChirp.m[plotFrame], times_inter), color=colors[4], label="ChirpSin", linestyle=:dash)
    startFreq = algChirp.f_TX_Array[plotFrame]
    endFreq = startFreq + algChirp.m[plotFrame]*sequence.t_TX
    freqEstimSin0(t) = startFreq
    lines!(axFreq, times_inter .* 1e3, freqEstimSin0.(times_inter), color=colors[6], label="F0Sin", linestyle=:dash)
    freqEstimSinMean(t) = (startFreq + endFreq)/2
    lines!(axFreq, times_inter .* 1e3, freqEstimSinMean.(times_inter), color=colors[7], label="FAvgSin", linestyle=:dash)
    xlims!(axFreq, 0, sequence.t_TX * 1e3)
    ylims!(axFreq, 785, 815)
    Legend(figFreqMaxAngle[2,1], axFreq, framevisible = false,        
    backgroundcolor = :transparent,orientation = :horizontal, nbanks = 1,)
  end
  maxDeflectionAngles = rad2deg.(maximum(deflectionAngleAllFrames, dims=1))
	lines!(axMaxAngle, 1:sequence.numFrames, maxDeflectionAngles[:], label = "$(param)", color=colorAlgs[algIndex],linestyle = (algIndex <4 ? :dash : :solid))
end
Legend(figFreqMaxAngle[4,1], axMaxAngle, framevisible = false,        
		backgroundcolor = :transparent,orientation = :horizontal, nbanks = 1,)

alg = algPulsed
	
dutyCycles = 0.4:0.1:1.0

maxDeflectionAnglesDutyCycle = zeros(length(dutyCycles))
@info "Duty cycle results: "
for (Dindex, D) in enumerate(dutyCycles)	
	algPulsed.D = D
	# initial settings for first frame
	algPulsed.f_TX_Array[1] = ω_nat/ (2*π) 							# initial frequency of excitation (in Hz)
	sequence.ψ_TX_Array[1] = 0.0 												# initial phase of excitation (in rad)
	param = "D = $(round(D, digits=1))"#
	sequence.excitationAlg = algPulsed
	setIandBAmpsForSequence!(sequence, BAmp_Max, E_max, R, factorTeslaToCurrent, factorCurrentToTesla, sequence.excitationAlg)

	# run simulation
	deflectionAngleAllFrames, ~ = simulate_meas(frameTime, sensor, sequence)
  maxDeflectionAnglesDutyCycle[Dindex] = maximum(rad2deg.(maximum(deflectionAngleAllFrames, dims=1)[:]))
	@info "D = $(round(D, digits=1)): max deflection angle = $(round(maxDeflectionAnglesDutyCycle[Dindex], digits=2)) deg"
end

# save figures
save("results/$(numPeriodsForTX)T_MMR_FrequencyAndMaxAngle.png", figFreqMaxAngle)
save("results/$(numPeriodsForTX)T_MMR_Signals.png", figSignal)