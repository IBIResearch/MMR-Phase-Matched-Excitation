function simulate_meas(frameTime, sensor, sequence)
	deflectionAngleAllFrames = zeros(Float32, sequence.numSamplesPerFrame,sequence.numFrames)
	deflectionAngleAllFramesDot = zeros(Float32, sequence.numSamplesPerFrame,sequence.numFrames)

	u₀ = [0.0, 0.0] # initial conditions for ODE solver
	alg = sequence.excitationAlg

	for frame_ in 1:sequence.numFrames
		# set up excitation field 
		ψ_TX = sequence.ψ_TX_Array[frame_] 
		f_Tx = alg.f_TX_Array[frame_] 

		BAmp = sequence.BAmp[frame_]
		tspan = (0.0, frameTime)

		function ode(du, u, p, t)
				θ = u[1]
				dθ = u[2]
				du[1] = dθ
				du[2] =  -sensor.ω_nat^2 * sin(θ) - sensor.ω_nat/sensor.Q * dθ + sensor.ms_rotor/sensor.inertia_rotor*(alg.B(t, BAmp,f_Tx, ψ_TX, alg.m[frame_],alg.D, sequence.t_TX)*cos(θ))	
		end

		# Pass to solvers
		prob = ODEProblem(ode, u₀, tspan)
		t = range(0.0, step=1/sequence.fs, length=sequence.numSamplesPerFrame)
		sol = solve(prob, Tsit5(), tstops=t)

    # get deflection angle
		deflectionAngleAllFrames[:,frame_] = [sol(t_)[1] for t_ in t]
		deflectionAngleAllFramesDot[:,frame_] = [sol(t_)[2] for t_ in t]

		# update excitation parameters for next frame based on the deflection angle in the current frame
		if frame_ < sequence.numFrames			
			rxTimes = range(sequence.t_TX+1/sequence.fs, length = sequence.numRxSamples, step=1/sequence.fs)
			f_startRx,f_endRX, ψ_endRX = getNewPhaseAndFreqs(sequence,deflectionAngleAllFrames[end-length(rxTimes)+1:end,frame_],rxTimes)
			if isa(alg, ChirpSin_Excitation) || isa(alg, ChirpRect_Excitation)
				alg.f_TX_Array[frame_+1] = f_endRX
				alg.m[frame_+1] = (f_startRx - f_endRX)/sequence.t_TX # update slope for next frame
			elseif isa(alg,FZeroSin_Excitation) || isa(alg, FZero_Excitation)
				alg.f_TX_Array[frame_+1] = f_endRX
			elseif  isa(alg, FAvgSin_Excitation) || isa(alg, FAvgRect_Excitation) 
				alg.f_TX_Array[frame_+1] = (f_endRX + f_startRx)/2 # use the average of the initial part frequency estimation and the estimation from the last part 
			else
				error("Unknown excitation algorithm type")
			end

			sequence.ψ_TX_Array[frame_+1] = ψ_endRX + pi/2 

		end
		u₀ = sol(t[end]) # use the last state as initial condition for the next frame
	
	end

	return deflectionAngleAllFrames, deflectionAngleAllFramesDot
end