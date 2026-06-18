

abstract type AbstractExcitationAlg end

mutable struct FZeroSin_Excitation <: AbstractExcitationAlg
  B::Function                       # waveform function
  m::Vector{Float64}                # (not used for fixed frequency excitation, included for consistency with ChirpSin_Excitation)
  f_TX_Array::Vector{Float64}       # frequency of excitation for each frame (in Hz)
  D::Float16                        # duty cycle for pulsed excitation
end
shortName(alg::FZeroSin_Excitation) =  "F0Sin"
mutable struct FAvgSin_Excitation <: AbstractExcitationAlg
  B::Function                 # waveform function    
  m::Vector{Float64}          # slope (not used for fixed frequency excitation, included for consistency with ChirpSin_Excitation)
  f_TX_Array::Vector{Float64} # frequency of excitation for each frame (in Hz)
  D::Float16                  # duty cycle for pulsed excitation
end
shortName(alg::FAvgSin_Excitation) = "FAvgSin"
mutable struct FAvgRect_Excitation <: AbstractExcitationAlg
  B::Function                 # waveform function
  m::Vector{Float64}          # slope (not used for pulsed excitation, included for consistency with ChirpSin_Excitation)
  f_TX_Array::Vector{Float64} # frequency of excitation for each frame (in Hz)
  D::Float16                  # duty cycle for pulsed excitation
end
shortName(alg::FAvgRect_Excitation) = "FAvgRect"

mutable struct FZero_Excitation <: AbstractExcitationAlg
  B::Function                 # waveform function  
  m::Vector{Float64}          # slope (not used for fixed frequency excitation, included for consistency with ChirpSin_Excitation)
  f_TX_Array::Vector{Float64} # frequency of excitation for each frame (in Hz)
  D::Float16                  # duty cycle for pulsed excitation
end
shortName(alg::FZero_Excitation) = "F0Rect" 
@kwdef mutable struct ChirpSin_Excitation <: AbstractExcitationAlg
  B::Function                 # waveform function
  m::Vector{Float64}          # slope of the chirp
  f_TX_Array::Vector{Float64} # frequency of excitation for each frame (in Hz)
  D::Float16                  # duty cycle for pulsed excitation
end
shortName(alg::ChirpSin_Excitation) = "ChirpSin"
@kwdef mutable struct ChirpRect_Excitation <: AbstractExcitationAlg
  B::Function                 # waveform function
  m::Vector{Float64}          # slope of the chirp
  f_TX_Array::Vector{Float64} # frequency of excitation for each frame (in Hz)
  D::Float16                  # duty cycle for pulsed excitation
end
shortName(alg::ChirpRect_Excitation) = "ChirpRect"

mutable struct MMRSequence
  numFrames::Int                            # number of frames in the sequence      
  numTxSamples::Int                         # number of samples for transmit window  
  numRxSamples::Int                         # number of samples for receive window
  numSamplesPerFrame::Int                   # total number of samples per frame (numTxSamples + numRxSamples)
  ψ_TX_Array::Vector{Float64}               # phase of excitation for each frame
  frequencyRange::Tuple{Float64, Float64}   # frequency range for spectrogram analysis (in Hz)
  fs::Float64                               # sampling frequency
  BAmp::Vector{Float64}                     # amplitude of excitation field (in Tesla)
  excitationAlg::Union{AbstractExcitationAlg, Nothing}      # excitation algorithm
  IAmp_MagV::Vector{Float64}                # amplitude of excitation current (in A)
  t_TX::Float64                             # duration of transmit window (in seconds)
end

mutable struct MMRSensor
  ω_nat::Float64                      # resonance angular frequency (in rad/s)
  r_rotor::Float64                    # radius of rotor (in m)
  v_rotor::Float64                    # volume of rotor (in m^3)
  m_rotor::Float64                    # mass of rotor (in kg)
  ρ_rotor::Float64                    # density of rotor (in kg/m^3)
  inertia_rotor::Float64              # inertia of rotor (in kg m^2)
  remanentMagnetization::Float64      # remanent magnetization of rotor (in uT)
  Q::Float64                          # quality factor
  ms_rotor::Float64                   # magnetic moment of dipole - rotator
end

"""
    I_RMS(Iamp, alg) -> Float64
    I_RMS(Iamp, sequence, currFrame, alg) -> Float64

Computes the root mean square (RMS) excitation current amplitude.
For pulsed waveforms the duty cycle scales the amplitude; for sinusoidal waveforms
a factor of 1/√2 is applied.
"""
function I_RMS(Iamp, alg::FAvgRect_Excitation)
  return Iamp * sqrt(alg.D[1])
end

function I_RMS(Iamp, alg::FZero_Excitation)
  return Iamp * sqrt(alg.D[1])
end

function I_RMS(Iamp, alg::ChirpRect_Excitation)
  return Iamp * sqrt(alg.D[1])
end

function I_RMS(Iamp, sequence, currFrame, alg::ChirpRect_Excitation)
  times_ = range(0, sequence.t_TX, length=10000)
  signal = alg.B.(times_, Iamp, alg.f_TX_Array[currFrame], sequence.ψ_TX_Array[currFrame], alg.m[currFrame], alg.D, sequence.t_TX)
  return sqrt(sum(signal.^2) / length(signal))
end

function I_RMS(Iamp, alg::FZeroSin_Excitation)
  return Iamp / sqrt(2)
end

function I_RMS(Iamp, alg::FAvgSin_Excitation)
  return Iamp / sqrt(2)
end

function I_RMS(Iamp, sequence, currFrame, alg::ChirpSin_Excitation)
  times_ = range(0, sequence.t_TX, length=10000)
  signal = alg.B.(times_, Iamp, alg.f_TX_Array[currFrame], sequence.ψ_TX_Array[currFrame], alg.m[currFrame], alg.D, sequence.t_TX)
  return sqrt(sum(signal.^2) / length(signal))
end

"""
    IFromE_max(E_max, R, t_TX, numFrames, alg) -> Float64

Computes the peak excitation current amplitude required to stay within the total
energy budget `E_max` given resistance `R`, transmit window
duration `t_TX`, and the number of frames. The conversion from RMS to peak
depends on the excitation waveform type.
"""
function IFromE_max(E_max, R, t_TX, numFrames, alg::FAvgRect_Excitation)
  I_RMS = sqrt(E_max / (R * t_TX * numFrames))
  return I_RMS / sqrt(alg.D)
end

function IFromE_max(E_max, R, t_TX, numFrames, alg::FZero_Excitation)
  I_RMS = sqrt(E_max / (R * t_TX * numFrames))
  return I_RMS / sqrt(alg.D)
end

function IFromE_max(E_max, R, t_TX, numFrames, alg::ChirpRect_Excitation)
  I_RMS = sqrt(E_max / (R * t_TX * numFrames))
  return I_RMS / sqrt(alg.D)
end

function IFromE_max(E_max, R, t_TX, numFrames, ::FZeroSin_Excitation)
  I_RMS = sqrt(E_max / (R * t_TX * numFrames))
  return I_RMS * sqrt(2)
end

function IFromE_max(E_max, R, t_TX, numFrames, ::FAvgSin_Excitation)
  I_RMS = sqrt(E_max / (R * t_TX * numFrames))
  return I_RMS * sqrt(2)
end


"""
    setIandBAmpsForSequence!(sequence, BAmp_Max, E_max, R,
                             factorTeslaToCurrent, factorCurrentToTesla, alg)

Sets the current amplitude and corresponding magnetic field amplitude for all frames so that the
total energy dissipated in resistance `R` does not exceed `E_max`.
For chirp-based excitation the amplitude is constrained to `BAmp_Max`;
for pulsed/sinusoidal excitation it is derived analytically from `E_max`.
For `ChirpRect_Excitation` the amplitude is found numerically.
"""
function setIandBAmpsForSequence!(sequence, BAmp_Max, E_max, R, factorTeslaToCurrent, factorCurrentToTesla, alg::FZeroSin_Excitation)
  Iv = IFromE_max(E_max / numFrames, R, sequence.t_TX, 1, alg)
  sequence.IAmp_MagV .= Iv	                                                # amplitude of excitation current
  sequence.BAmp .= sequence.IAmp_MagV .* factorCurrentToTesla[1] 			# amplitude of excitation field 
end					

function setIandBAmpsForSequence!(sequence, BAmp_Max, E_max, R, factorTeslaToCurrent, factorCurrentToTesla, alg::FAvgSin_Excitation)
  Iv = IFromE_max(E_max/numFrames, R, sequence.t_TX, 1, alg)
  sequence.IAmp_MagV .= Iv 	                                                # amplitude of excitation current
  sequence.BAmp .= sequence.IAmp_MagV .* factorCurrentToTesla[1] 			# amplitude of excitation field
end							


function setIandBAmpsForSequence!(sequence, BAmp_Max, E_max, R, factorTeslaToCurrent, factorCurrentToTesla, alg::FAvgRect_Excitation)
  Iv = IFromE_max(E_max/numFrames, R, sequence.t_TX, 1, alg)
  sequence.IAmp_MagV .= Iv 	                                                    # amplitude of excitation current
  sequence.BAmp .= sequence.IAmp_MagV .* factorCurrentToTesla[1] 					# amplitude of excitation field
  E_total = R * I_RMS(Iv, alg)^2 * sequence.t_TX * numFrames
  @info "Difference total energy consumption and E_max: $(E_total - E_max)"			
end	

function setIandBAmpsForSequence!(sequence, BAmp_Max, E_max, R, factorTeslaToCurrent, factorCurrentToTesla, alg::FZero_Excitation)
  Iv = IFromE_max(E_max/numFrames, R, sequence.t_TX, 1, alg)
  sequence.IAmp_MagV .= Iv 	                                                    # amplitude of excitation current

  sequence.BAmp .= sequence.IAmp_MagV .* factorCurrentToTesla[1] 					# amplitude of excitation field
  E_total = R * I_RMS(Iv, alg)^2 * sequence.t_TX * numFrames
  @info "Difference total energy consumption and E_max: $(E_total - E_max)"				
end	


function setIandBAmpsForSequence!(sequence, BAmp_Max, E_max_Chirp, R, factorTeslaToCurrent, factorCurrentToTesla, alg::ChirpRect_Excitation)
  function f(p)
    IAmp_MagV = p
    sequence.IAmp_MagV .= IAmp_MagV 	                                     
    sequence.BAmp .= sequence.IAmp_MagV .* factorCurrentToTesla[1] 			
    # simulate measurement to compute the actual energy consumption for this amplitude setting
    ~, ~ = simulate_meas(frameTime, sensor,sequence );
    E_max = getEmaxForAlg(sequence.BAmp, factorTeslaToCurrent, sequence, alg, R)
    return (E_max - E_max_Chirp)^2
  end
  # as it will be in a similar range as for the averaged frequency pulsed excitation, we set the boundaries close to the analytical solution  
  lower = sequence.IAmp_MagV[1]*0.8
  upper = sequence.IAmp_MagV[1]*1.2
  result = optimize(x -> f([x]), lower, upper)
  Iv = Optim.minimizer(result)
  E_max = getEmaxForAlg(sequence.BAmp, factorTeslaToCurrent, sequence, alg, R)
  sequence.IAmp_MagV .= Iv                                                        # amplitude of excitation current
  sequence.BAmp .= sequence.IAmp_MagV .* factorCurrentToTesla[1] 						    # amplitude of excitation field
  @info "Difference total energy consumption and E_max: $(E_max - E_max_Chirp)"							
end	

function setIandBAmpsForSequence!(sequence, BAmp_Max, E_max, R, factorTeslaToCurrent, factorCurrentToTesla, alg ::ChirpSin_Excitation)
  sequence.IAmp_MagV .= BAmp_Max * factorTeslaToCurrent[1] 	                      # amplitude of excitation current
  sequence.BAmp .= sequence.IAmp_MagV .* factorCurrentToTesla[1] 						# amplitude of excitation field
end			

"""
    getEmaxForAlg(BAmp, factorTeslaToCurrent, sequence, alg, R) -> Float64

Computes the total energy dissipated in resistance `R` over all frames
of `sequence` by numerically integrating the waveform for each frame.
`BAmp` is the per-frame field amplitude array (in Tesla) and `factorTeslaToCurrent`
converts it to current.
"""
function getEmaxForAlg(BAmp, factorTeslaToCurrent, sequence, alg::ChirpSin_Excitation, R)
  E_max = 0.0
  for frame_ in 1:sequence.numFrames
    BAmp = sequence.BAmp[frame_]
    E_max += R * I_RMS(BAmp * factorTeslaToCurrent[1], sequence, frame_, alg)^2 * sequence.t_TX
  end
  return E_max
end

function getEmaxForAlg(BAmp, factorTeslaToCurrent, sequence, alg::ChirpRect_Excitation, R)
  E_max = 0.0
  for frame_ in 1:sequence.numFrames
    E_max += R * I_RMS(BAmp[frame_] * factorTeslaToCurrent[1], sequence, frame_, alg)^2 * sequence.t_TX
  end
  return E_max
end

"""
    getNewPhaseAndFreqs(sequence, data, rxTimes) -> (f_start, f_end, phase)

Estimates the instantaneous frequency and phase of the received signal based on
zero-crossing detection. Linear interpolation is used to refine the crossing times
beyond the sample grid. The phase is computed from the time elapsed between the last zero crossing and the
end of the receive window, corrected for a potential negative final half-cycle.
"""
function getNewPhaseAndFreqs(sequence, data, rxTimes)
  zeroCrossings = Float64[]
  for i in 1:length(data)-1
    if data[i] * data[i+1] < 0
      # Linear interpolation to refine the zero-crossing time
      t1, t2 = rxTimes[i], rxTimes[i+1]
      u1, u2 = data[i], data[i+1]
      t_crossing = t1 + (0 - u1) * (t2 - t1) / (u2 - u1)
      push!(zeroCrossings, t_crossing)
    end
  end
  freqs_ = 1 ./ diff(zeroCrossings) / 2
  freqs_ = imfilter(freqs_, Kernel.gaussian((4,)))
  deltaT = sequence.numSamplesPerFrame / sequence.fs - zeroCrossings[end]
  if data[end] < 0
    deltaT += 1 / (2 * freqs_[end])
  end
  phase = mod2pi(2 * π * freqs_[end] * deltaT)
  return freqs_[1], freqs_[end], phase
end



# Configure one left/right axis pair used by the zoom subplots.
function configure_zoom_axes!(axL, axR, algIndex, numAlgs, BAmp_Max, ylabelTickLabel)
	axL.yticks = [-4000, 0, 4000]
	ylims!(axL, -4000, 4000)
	axL.title = ylabelTickLabel

	axR.yticks = [-BAmp_Max * 1e6, 0, BAmp_Max * 1e6]
	ylims!(axR, -BAmp_Max * 1e6, BAmp_Max * 1e6)
	hidexdecorations!(axR)
	axR.xgridvisible = false
	axR.ygridvisible = false

	if algIndex != numAlgs
		axL.xticklabelsvisible = false
	end
end

# Plot one zoom window for phi-dot and B while keeping both y-axes aligned.
function plot_zoom_window!(axL, axR, t, phiDot, B, startIdx, stopIdx)
	t_ms = t[startIdx:stopIdx] .* 1e3

	lines!(axL, t_ms, phiDot[startIdx:stopIdx], color = colors[1])
	lines!(axR, t_ms, B[startIdx:stopIdx] .* 1e6, color = colors[3])
	xlims!(axL, t_ms[1], t_ms[end])
	xlims!(axR, t_ms[1], t_ms[end])
end
	