#########################################################################################

colors = [
    RGBf(0/255, 73/255, 146/255),     # blue
    RGBf(239/255, 123/255, 5/255),    # orange (dark)
    RGBf(138/255, 189/255, 36/255),   # green
    RGBf(178/255, 34/255, 41/255),    # red
    RGBf(170/255, 156/255, 143/255),  # mocca
    RGBf(87/255, 87/255, 86/255),     # black (text)
    RGBf(255/255, 223/255, 0/255),    # yellow
    RGBf(104/255, 195/255, 205/255),  # TUHH
    RGBf(45/255, 198/255, 214/255),   # TUHH
    RGBf(193/255, 216/255, 237/255)
]

#########################################################################################
# Coil and sensor constants.
R = 4.0                                  # coil resistance in Ohm
factorCurrentToTesla = 0.0008794831466724841  # T/A
factorTeslaToCurrent = 1 / factorCurrentToTesla

# Sensor parameters for the MMRS model.
ω_nat = 2 * π * 814                     # resonance angular frequency in rad/s
r_rotor = 0.5e-3                        # rotor radius in m
v_rotor = 4 / 3 * π * r_rotor^3         # rotor volume in m^3
m_rotor = 4.14e-6                       # rotor mass in kg
ρ_rotor = m_rotor / v_rotor             # rotor density in kg/m^3
inertia_rotor = 2 / 5 * ρ_rotor * r_rotor^2 * v_rotor
remanentMagnetization = 1.44e6          # remanent magnetization in uT
Q = 852                                 # quality factor
ms_rotor = v_rotor * remanentMagnetization
sensor = MMRSensor(
    ω_nat,
    r_rotor,
    v_rotor,
    m_rotor,
    ρ_rotor,
    inertia_rotor,
    remanentMagnetization,
    Q,
    ms_rotor,
)

# Sequence parameters.
fs = 100000.0                           # sampling frequency in Hz
numPeriodsForTX = 60
numFrames = 6                           # number of frames to simulate
t_TX = numPeriodsForTX * (2 * π / sensor.ω_nat)
t_RX = 0.3                              # receive window duration in s
numTxSamples = trunc(Int64, t_TX * fs)
numRxSamples = trunc(Int64, t_RX * fs)
numSamplesPerFrame = numTxSamples + numRxSamples
frameTime = numSamplesPerFrame / fs     # frame duration in s
frequencyRange = (0.0, 2000.0)          # spectrum analysis range in Hz
BAmp_Max = 40e-6                        # peak excitation field amplitude in Tesla
sequence = MMRSequence(
    numFrames,
    numTxSamples,
    numRxSamples,
    numSamplesPerFrame,
    zeros(Float64, numFrames),
    frequencyRange,
    fs,
    zeros(Float64, numFrames),
    nothing,
    zeros(Float64, numFrames),
    t_TX,
)
#########################################################################################
# Excitation algorithms.

D = 0.7 # duty cycle for pulsed excitation algorithms.

fixed_frequency_sin_waveform(t, BAmp, fTX, ψ_TX, m, D, endTimes) =
    BAmp * sin(2 * π * fTX * t + ψ_TX) * (t < endTimes)

chirp_sin_waveform(t, BAmp, fTX, ψ_TX, m, D, endTimes) =
    BAmp * sin(2 * π * (m / 2 * t^2 + fTX * t) + ψ_TX) * (t < endTimes)

function fixed_frequency_rect_waveform(t, BAmp, fTX, ψ_TX, m, D, endTimes)
    ϕ = mod(t * fTX + ψ_TX / (2π), 1.0)
    signal = (0.25 - D / 4 ≤ ϕ < 0.25 + D / 4) ? BAmp :
             (0.75 - D / 4 ≤ ϕ < 0.75 + D / 4) ? -BAmp :
             0.0
    return signal * (t < endTimes)
end

function chirp_rect_waveform(t, BAmp, fTX, ψ_TX, m, D, endTimes)
    ϕ = mod(m / 2 * t^2 + fTX * t + ψ_TX / (2π), 1.0)
    signal = (0.25 - D / 4 ≤ ϕ < 0.25 + D / 4) ? BAmp :
             (0.75 - D / 4 ≤ ϕ < 0.75 + D / 4) ? -BAmp :
             0.0
    return signal * (t < endTimes)
end

algSin = FAvgSin_Excitation(fixed_frequency_sin_waveform, zeros(Float64, numFrames), zeros(Float64, numFrames), 0.0)
algSin0 = FZeroSin_Excitation(fixed_frequency_sin_waveform, zeros(Float64, numFrames), zeros(Float64, numFrames), 0.0)
algChirp = ChirpSin_Excitation(chirp_sin_waveform, zeros(Float64, numFrames), zeros(Float64, numFrames), 0.0)
algChirp.m[1] = -100

algPulsed = FAvgRect_Excitation(fixed_frequency_rect_waveform, zeros(Float64, numFrames), zeros(Float64, numFrames), D)
algPulsed0 = FZero_Excitation(fixed_frequency_rect_waveform, zeros(Float64, numFrames), zeros(Float64, numFrames), D)
algChirpPulsed = ChirpRect_Excitation(
chirp_rect_waveform,
zeros(Float64, numFrames),
zeros(Float64, numFrames),
D,
)
algChirpPulsed.m[1] = -100