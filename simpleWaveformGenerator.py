import numpy as np
import matplotlib.pyplot as plt
# from sympy import rf
from External import adiabatic
from External import slr 
# print(slr.__file__)

## Adapted and extended from the Pulpy library by J.B. Martin (https://github.com/jonbmartin/pulpy/tree/master)

## Gradient waveform designers

def trap_grad(ramp_up, ramp_down, plateau, dt):
    r"""Trapezoidal gradient designer. Design for specific ramp up, ramp down and 
    plateau times.

    Args:
        ramp_up (float): ramp up time in sec
        ramp_down (float): ramp down time in sec
        plateau (float): plateau time in sec
        dt (float): sample time in sec

    Returns:
        2-element tuple containing

        - **trap** (*array*): gradient waveform normalized between 0 an 1.
        - **amplitude** (*int*): gradient amplitude in mT/m. Fixed to 1 for this function.
        - **ramp_up_time** (*float*): ramp up time in us.
        - **ramp_down_time** (*float*): ramp down time in us.
        - **plateau_time** (*float*): plateau time in us.

    """

    # make a flat portion of magnitude 1 and enough area for the swath
    pts = np.floor(plateau / dt)
    flat = np.ones((1, int(pts)))

    # make attack and decay ramps
    ramp_up_pts = int(np.ceil(ramp_up / dt))
    # print(f"Ramp up points: {ramp_up_pts}, Ramp down points: {int(np.ceil(ramp_down / dt))}")
    ramp_up = np.linspace(0, ramp_up_pts, num=ramp_up_pts + 1) / ramp_up_pts * np.max(flat)
    ramp_down_pts = int(np.ceil(ramp_down / dt))
    ramp_down = np.linspace(ramp_down_pts, 0, num=ramp_down_pts + 1) / ramp_down_pts * np.max(flat)

    trap = np.concatenate((ramp_up, np.squeeze(flat), ramp_down))

    ramp_up_time = ramp_up_pts * dt * 1e6  # convert to us
    ramp_down_time = ramp_down_pts * dt * 1e6  # convert to us
    plateau_time = (len(trap)) * dt * 1e6 - (ramp_up_time + ramp_down_time)
    amplitude = 1  # fixed amplitude for this function

    return np.expand_dims(trap, axis=0), amplitude, int(ramp_up_time), int(ramp_down_time), int(plateau_time)

def min_trap_grad(area, gmax, slew, dt): 
    r"""Minimal duration trapezoidal gradient designer. Design for target area
    under the flat portion (for non-ramp-sampled pulses)

    Args:
        area (float): pulse area in mT/m/s
        gmax (float): maximum gradient in mT/m
        slew (float): max slew rate in T/m/s
        dt (float): sample time in sec

    Returns:
        2-element tuple containing

        - **normalized_trap** (*array*): normalize gradient waveform between 0 an 1.
        - **amplitude** (*int*): gradient amplitude in mT/m.
        - **ramp_up_time** (*float*): ramp up time in us.
        - **ramp_down_time** (*float*): ramp down time in us.
        - **plateau_time** (*float*): plateau time in us.

    """
    slew = slew * 1e-3  # convert slew rate to mT/m/s

    if np.abs(area) > 0:
        # we get the solution for plateau amp by setting derivative of
        # duration as a function of amplitude to zero and solving
        a = np.sqrt(slew * area / 2)

        # finish design with discretization
        # make a flat portion of magnitude a and enough area for the swath
        pts = np.floor(area / a / dt)
        flat = np.ones((1, int(pts)))
        flat = flat / np.sum(flat) * area / dt
        if np.max(flat) > gmax:
            flat = np.ones((1, int(np.ceil(area / gmax / dt))))
            flat = flat / np.sum(flat) * area / dt

        # make attack and decay ramps
        ramppts = int(np.ceil(np.max(flat) / slew / dt))
        ramp_up = np.linspace(0, ramppts, num=ramppts + 1) / ramppts * np.max(flat)
        ramp_dn = np.linspace(ramppts, 0, num=ramppts + 1) / ramppts * np.max(flat)

        trap = np.concatenate((ramp_up, np.squeeze(flat), ramp_dn))

    else:
        # negative-area trap requested?
        trap, ramppts = 0, 0

    amplitude = np.max(np.abs(trap))
    normalized_trap = trap / amplitude

    print("ramp points", ramppts)
    print("plateau points", (len(normalized_trap) - 2 * ramppts))
    print("amplitude", amplitude)

    ramp_up_time = ramppts * dt * 1e6  # convert to 1e-6 us
    ramp_down_time = ramppts * dt * 1e6  # convert to 1e-6 us
    plateau_time = (len(normalized_trap) - 2 * ramppts) * dt * 1e6  # convert to 1e-6 us

    return np.expand_dims(normalized_trap, axis=0), amplitude, int(ramp_up_time), int(ramp_down_time), int(plateau_time)


def ramp_sampled_trap_grad(area, gmax, slew, dt, *args):
    r"""General trapezoidal gradient designer for total target area
    (for rewinders)

    Args:
        area (float): pulse area in mT/m/s
        gmax (float): maximum gradient in mT/m
        slew (float): max slew rate in T/m/s
        dt (float): sample time in sec

    Returns:
        2-element tuple containing

        - **normalized_trap** (*array*): normalize gradient waveform between 0 an 1.
        - **amplitude** (*int*): gradient amplitude in mT/m.
        - **ramp_up_time** (*float*): ramp up time in us.
        - **ramp_down_time** (*float*): ramp down time in us.
        - **plateau_time** (*float*): plateau time in us.

    """

    slew = slew * 1e-3  # convert slew rate to mT/m/s

    if len(args) < 5:
        # in case we are making a rewinder
        rampsamp = 1

    if np.abs(area) > 0:
        if rampsamp:
            ramppts = int(np.ceil(gmax / slew / dt))
            triareamax = ramppts * dt * gmax
            sign = 1
            if area<0:
                area = -area
                sign = -1
            if triareamax > np.abs(area):
                # triangle pulse
                newgmax = np.sqrt(np.abs(area) * slew)
                ramppts = int(np.ceil(newgmax / slew / dt))
                ramp_up = np.linspace(0, ramppts, num=ramppts + 1) / ramppts
                ramp_dn = np.linspace(ramppts, 0, num=ramppts + 1) / ramppts
                pulse = sign * np.concatenate((ramp_up, ramp_dn))
            else:
                # trapezoid pulse
                nflat = int(np.ceil((area - triareamax) / gmax / dt / 2) * 2)
                # nflat = np.abs(nflat) # this is a test
                ramp_up = np.linspace(0, ramppts, num=ramppts + 1) / ramppts
                ramp_dn = np.linspace(ramppts, 0, num=ramppts + 1) / ramppts
                pulse = sign * np.concatenate((ramp_up, np.ones(nflat), ramp_dn))

            trap = pulse * (area / (sum(pulse) * dt))

        else:
            # make a flat portion of magnitude gmax
            # and enough area for the entire swath
            flat = np.ones(1, np.ceil(area / gmax / dt))
            flat = flat / sum(flat) * area / dt
            flat_top = np.max(flat)

            # make attack and decay ramps
            ramppts = int(np.ceil(np.max(flat) / slew / dt))
            ramp_up = np.linspace(0, ramppts, num=ramppts + 1) / ramppts * flat_top
            ramp_dn = np.linspace(ramppts, 0, num=ramppts + 1) / ramppts * flat_top
            trap = np.concatenate((ramp_up, flat, ramp_dn))

    else:
        trap, ramppts = 0, 0

    amplitude = np.max(np.abs(trap))
    normalized_trap = trap / amplitude

    ramp_up_time = ramppts * dt * 1e6  # convert to us
    ramp_down_time = ramppts * dt * 1e6  # convert to us
    plateau_time = (len(normalized_trap) - 2 * ramppts) * dt * 1e6  # convert to us

    return np.expand_dims(normalized_trap, axis=0), amplitude, int(ramp_up_time), int(ramp_down_time), int(plateau_time)

def test_trap_grad():
    trap, ampl, ramp_up_time, ramp_down_time, plateau_duration = trap_grad(200, 200, 2560, 1e-5)
    # print(f"Ramp Up Time: {ramp_up_time} us, Ramp Down Time: {ramp_down_time} us, Plateau Time: {plateau_duration} us")
    plt.plot(trap[0]*ampl)
    plt.title("Trapezoidal Gradient Waveform")
    plt.xlabel("Time (10us)")
    plt.ylabel("Gradient (mT/m)")
    plt.show()

def test_min_trap_grad():   
    morm_trap, ampl, ramp_up_time, ramp_down_time, plateau_duration = min_trap_grad(200e-5, 20, 200, 1e-5)
    # print(f"Ramp Up Time: {ramp_up_time} us, Ramp Down Time: {ramp_down_time} us, Plateau Time: {plateau_duration} us")
    plt.plot(morm_trap[0]*ampl)
    plt.title("Minimal Duration Trapezoidal Gradient Waveform")
    plt.xlabel("Time (10us)")
    plt.ylabel("Gradient (mT/m)")
    plt.show()

def test_ramp_sampled_trap_grad():      
    norm_trap, ampl, ramp_up_time, ramp_down_time, plateau_duration = ramp_sampled_trap_grad(200e-5, 20, 200, 1e-5)
    # print(f"Ramp Up Time: {ramp_up_time} us, Ramp Down Time: {ramp_down_time} us, Plateau Time: {plateau_duration} us")
    plt.plot(norm_trap[0]*ampl)
    plt.title("Ramp Sampled Trapezoidal Gradient Waveform")
    plt.xlabel("Time (10us)")
    plt.ylabel("Gradient (mT/m)")
    plt.show()

# Run tests
# test_trap_grad()
# test_min_trap_grad()
# test_ramp_sampled_trap_grad()



## RF pulse designers

def pulse_designer(pulse_type, args):
    r"""General RF pulse designer. This function is a placeholder for future
    implementations of various RF pulse designs.

    Args:
        pulse_type (str): type of RF pulse to design: available options include 'slr', 'sinc', 'adiabatic', etc.
        args (list): additional arguments for the specific pulse design.

    Returns:
        None: This function currently does not return any value.
    """
    print(f"Designing {pulse_type} RF pulse with arguments: {args}")
    match pulse_type:
        case 'slr':
            tb = args[0]        # RF pulse time-bandwidth product
            N =  args[1]        # number of samples
            d1 =  args[2]       # magnetization passband ripple level
            d2 =  args[3]       # magnetization stopband ripple level
            p_type =  args[4]   # RF pulse type: 
                                # 'st' (small-tip excitation),
                                # 'ex' (pi/2excitation pulse), 
                                # 'se' (spin-echo pulse), 
                                # 'inv' (inversion), 
                                # or 'sat' (pi/2 saturation pulse).
            f_type =  args[5]   # filter type:
                                # 'ms' (sinc), 
                                # 'pm'(Parks-McClellan equal-ripple), 
                                # 'min' (minphase using factored pm),
                                # 'max' (maxphase using factored pm), 
                                # 'ls' (least squares).
            pulse = slr.dzrf(N, tb, p_type, f_type, d1, d2, True)
            envelope = pulse / np.max(np.abs(pulse))
            phase_value = 0
            rf = envelope * np.exp(1j * phase_value)
            magnitude = np.abs(rf)
            phase = np.angle(rf)
        case 'sinc':
            n = args[0]  # number of samples  
            m = args[1]  # number of side lobes
            pulse = slr.msinc(n, m)
            envelope = pulse / np.max(np.abs(pulse))
            phase_value = 0
            rf = envelope * np.exp(1j * phase_value)
            magnitude = np.abs(rf)
            phase = np.angle(rf)
        case 'adiabatic':
            type = args[0]  # adiabatic pulse type: 'bir4', 'wurst', 'hyperbolic'
            n = args[1]  # number of samples  (should be a multiple of 4)
            match type:
                case 'bir4':
                    beta = args[2]  # AM waveform parameter.
                    kappa = args[3]  # FM waveform parameter.
                    theta = args[4]  # flip angle in radians.
                    dw0 = args[5]  # FM waveform scaling (radians/s).
                    pulse, freq_mod = adiabatic.bir4(n, beta, kappa, theta, dw0)
                case 'wurst':
                    n_fac = args[2]  # power to exponentiate to within AM term. 
                                     # ~20 or greater istypical.
                    bw = args[3]  # bandwidth
                    dur = args[4]  # duration in seconds
                    pulse, freq_mod = adiabatic.wurst(n, n_fac, bw, dur)
                case 'hyperbolic':
                    beta = args[2]  # AM waveform parameter.
                    mu = args[3]  # a constant, determines amplitude of frequency sweep.
                    dur = args[4]  # duration in seconds
                    pulse, freq_mod = adiabatic.hypsec(n, beta, mu, dur)
            
            dt = 1e-5
            cumulative_phase = np.cumsum(freq_mod) * dt 
            wrapped_phase = (cumulative_phase + np.pi) % (2 * np.pi) - np.pi
            magnitude = np.abs(pulse) / np.max(np.abs(pulse))  # Normalize the pulse
            phase = wrapped_phase  # Use the wrapped phase as the phase of the pulse
        
    return magnitude, phase

## Running the tests
# rfType = "SLR"
# magnitude, phase = pulse_designer("slr", [8, 512, 0.01, 0.01, 'ex', 'ls'])
# rfType = "Sinc"
# magnitude, phase = pulse_designer("sinc", [64, 2])
# rfType = "adiabatic wurst"
# magnitude, phase = pulse_designer("adiabatic", ['wurst', 512, 40, 40e3, 2e-3])
# rfType = "adiabatic bir4"
# # magnitude, phase = pulse_designer("adiabatic", ['bir4', 512, 10, np.arctan(20), np.pi/4, 100*np.pi/1e-5/512])
# rfType = "adiabatic hyperbolic"
# magnitude, phase = pulse_designer("adiabatic", ['hyperbolic', 512, 800, 4.9, 0.012])
# fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
# ax1.plot(magnitude, label='Magnitude')
# ax1.legend()
# ax1.set_title("RF Waveform " + rfType)
# ax2.plot(phase, color='orange', label='Phase')
# ax2.legend()
# plt.tight_layout()
# plt.show()
